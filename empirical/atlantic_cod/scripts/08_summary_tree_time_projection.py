#!/usr/bin/env python3
"""Stage 8: project SNAPP MRCA times onto fixed ASTRAL topologies.

This post-freeze extension uses published SNAPP window-level MRCA times from
Stage 7 and fixed ASTRAL topologies from Stage 6. It does not rerun SNAPP or
ASTRAL4 and does not interpret ASTRAL branch lengths as calibrated time.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import platform
import statistics
import tempfile
import textwrap
import unittest
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import Bio
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMP_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP_ROOT / "results"
FIGURES = EMP_ROOT / "figures"
ASTRAL = RESULTS / "astral4"

STAGE7_PAIRWISE = DATA_ROOT / "processed" / "stage7_pairwise_mrca_times.tsv"
STAGE7_TESTS = RESULTS / "stage7_divergence_time_tests.tsv"
STAGE5_MANIFEST = RESULTS / "stage5_final_manifest.json"
STAGE6_MANIFEST = RESULTS / "stage6_manifest.json"
T0_TOPOLOGY = ASTRAL / "T0_all.inferred.nwk"
T1_TOPOLOGY = ASTRAL / "T1_outside.inferred.nwk"

PAIRWISE_MATRICES = DATA_ROOT / "processed" / "stage8_pairwise_mrca_matrices.tsv"
ROOTING_AUDIT = RESULTS / "stage8_rooting_audit.tsv"
NODE_AGE_SHIFTS = RESULTS / "stage8_node_age_shifts.tsv"
SUMMARY = RESULTS / "stage8_summary.tsv"
PROJECTION_FIT = RESULTS / "stage8_projection_fit.tsv"
PROJECTION_RESIDUALS = RESULTS / "stage8_projection_residuals.tsv"
T0_OUTSIDE_TREE = RESULTS / "stage8_T0_all_topology_outside_times.nwk"
T0_ALL_TREE = RESULTS / "stage8_T0_all_topology_all_times.nwk"
T1_OUTSIDE_TREE = RESULTS / "stage8_T1_outside_topology_outside_times.nwk"
T1_ALL_TREE = RESULTS / "stage8_T1_outside_topology_all_times.nwk"
METHODS_TEXT = RESULTS / "stage8_methods_text.md"
RESULTS_TEXT = RESULTS / "stage8_results_text.md"
FIGURE_CAPTION = RESULTS / "stage8_figure_caption.md"
REPORT = RESULTS / "stage8_report.md"
MANIFEST = RESULTS / "stage8_manifest.json"
FIG_SUMMARY_PDF = FIGURES / "atlantic_cod_stage8_summary_tree_time_projection.pdf"
FIG_SUMMARY_PNG = FIGURES / "atlantic_cod_stage8_summary_tree_time_projection.png"
FIG_HEATMAP_PDF = FIGURES / "atlantic_cod_stage8_pairwise_time_shift_heatmap.pdf"
FIG_HEATMAP_PNG = FIGURES / "atlantic_cod_stage8_pairwise_time_shift_heatmap.png"

EXPECTED_STAGE5 = {
    "data/atlantic_cod/processed/cod_window_trees.tsv": "26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b",
    "data/atlantic_cod/processed/stage4a_window_topology_scores.tsv": "3ef09caa4edddd247ad986ff168a4fc86b79d66559293bec51d291ee8cf9328e",
}
EXPECTED_STAGE7_CIRCULAR = {
    "LG01": 0.00917431193,
    "LG02": 0.010989011,
    "LG07": 0.00854700855,
    "LG12": 0.00961538462,
}


@dataclass(eq=False)
class UNode:
    id: int
    name: str | None = None
    neighbors: set[int] = field(default_factory=set)


@dataclass
class RNode:
    id: str
    name: str | None
    children: list["RNode"] = field(default_factory=list)
    parent: "RNode | None" = None
    age: float = 0.0

    @property
    def is_tip(self) -> bool:
        return not self.children


@dataclass
class FitResult:
    topology: str
    treatment: str
    objective_type: str
    root_edge: str
    rooted_root: RNode
    ages: dict[str, float]
    objective: float
    rmse: float
    mae: float
    corr: float
    max_residual_abs: float
    residual_rows: list[dict[str, object]]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def verify_frozen_inputs() -> dict[str, str]:
    observed = {}
    for rel, expected in EXPECTED_STAGE5.items():
        path = REPO_ROOT / rel
        got = sha256(path)
        if got != expected:
            raise ValueError(f"Frozen Stage-5 input changed: {rel}: {got} != {expected}")
        observed[rel] = got
    for path in (T0_TOPOLOGY, T1_TOPOLOGY):
        if not path.exists():
            raise FileNotFoundError(path)
        observed[str(path.relative_to(REPO_ROOT))] = sha256(path)
    tests = {row["lg"]: float(row["circular_p"]) for row in read_tsv(STAGE7_TESTS)}
    for lg, expected in EXPECTED_STAGE7_CIRCULAR.items():
        if not math.isclose(tests[lg], expected, rel_tol=0, abs_tol=5e-10):
            raise ValueError(f"Stage-7 circular p-value changed for {lg}: {tests[lg]} != {expected}")
    return observed


def parse_topology(path: Path) -> dict[int, UNode]:
    tree = Phylo.read(io.StringIO(path.read_text().strip()), "newick")
    nodes: dict[int, UNode] = {}
    next_id = 0
    clade_to_id = {}
    for clade in tree.find_clades(order="preorder"):
        clade_to_id[clade] = next_id
        nodes[next_id] = UNode(next_id, clade.name if clade.is_terminal() else None)
        next_id += 1
    for clade in tree.find_clades(order="preorder"):
        parent_id = clade_to_id[clade]
        for child in clade.clades:
            child_id = clade_to_id[child]
            nodes[parent_id].neighbors.add(child_id)
            nodes[child_id].neighbors.add(parent_id)
    collapse_unlabeled_degree_two_nodes(nodes)
    return nodes


def collapse_unlabeled_degree_two_nodes(nodes: dict[int, UNode]) -> None:
    changed = True
    while changed:
        changed = False
        for node_id, node in list(nodes.items()):
            if node.name is None and len(node.neighbors) == 2:
                a, b = sorted(node.neighbors)
                nodes[a].neighbors.discard(node_id)
                nodes[b].neighbors.discard(node_id)
                nodes[a].neighbors.add(b)
                nodes[b].neighbors.add(a)
                del nodes[node_id]
                changed = True
                break


def taxa_of_graph(nodes: dict[int, UNode]) -> list[str]:
    return sorted(node.name for node in nodes.values() if node.name is not None)


def edge_list(nodes: dict[int, UNode]) -> list[tuple[int, int]]:
    return sorted((min(a, b), max(a, b)) for a, node in nodes.items() for b in node.neighbors if a < b)


def component_taxa(nodes: dict[int, UNode], start: int, blocked: tuple[int, int]) -> set[str]:
    blocked_set = {blocked, (blocked[1], blocked[0])}
    seen = {start}
    stack = [start]
    taxa = set()
    while stack:
        current = stack.pop()
        if nodes[current].name is not None:
            taxa.add(str(nodes[current].name))
        for nb in nodes[current].neighbors:
            if (current, nb) in blocked_set or nb in seen:
                continue
            seen.add(nb)
            stack.append(nb)
    return taxa


def root_edge_label(nodes: dict[int, UNode], edge: tuple[int, int]) -> str:
    all_taxa = set(taxa_of_graph(nodes))
    left = component_taxa(nodes, edge[0], edge)
    right = all_taxa - left
    if (len(right), sorted(right)) < (len(left), sorted(left)):
        left, right = right, left
    return "+".join(sorted(left)) + "|" + "+".join(sorted(right))


def rooted_from_edge(nodes: dict[int, UNode], edge: tuple[int, int]) -> RNode:
    root = RNode("root", None)
    counter = 0

    def orient(current: int, parent: int, parent_r: RNode) -> RNode:
        nonlocal counter
        u = nodes[current]
        if u.name is not None:
            rid = f"tip:{u.name}"
            r = RNode(rid, u.name, parent=parent_r)
            parent_r.children.append(r)
            return r
        counter += 1
        r = RNode(f"n{counter:02d}", None, parent=parent_r)
        parent_r.children.append(r)
        for nb in sorted(u.neighbors):
            if nb != parent:
                orient(nb, current, r)
        return r

    orient(edge[0], edge[1], root)
    orient(edge[1], edge[0], root)
    return root


def all_rnodes(root: RNode) -> list[RNode]:
    out = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(reversed(node.children))
    return out


def tip_order(root: RNode) -> list[str]:
    return [node.name for node in all_rnodes(root) if node.is_tip]


def clade_taxa(node: RNode) -> tuple[str, ...]:
    if node.is_tip:
        return (str(node.name),)
    taxa = []
    for child in node.children:
        taxa.extend(clade_taxa(child))
    return tuple(sorted(taxa))


def internal_nodes(root: RNode) -> list[RNode]:
    return [node for node in all_rnodes(root) if not node.is_tip]


def ancestor_map(root: RNode) -> dict[str, list[RNode]]:
    out = {}
    for node in all_rnodes(root):
        if node.is_tip:
            path = []
            cur: RNode | None = node
            while cur is not None:
                path.append(cur)
                cur = cur.parent
            out[str(node.name)] = path
    return out


def mrca_node(root: RNode, a: str, b: str) -> RNode:
    ancestors = ancestor_map(root)
    a_path = ancestors[a]
    b_set = {node.id for node in ancestors[b]}
    for node in a_path:
        if node.id in b_set:
            return node
    raise ValueError(f"No MRCA found for {a}, {b}")


def unordered_pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def pair_class_counts() -> dict[tuple[str, str], dict[str, int]]:
    counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    seen_window_pair: set[tuple[str, str, str]] = set()
    for row in read_tsv(STAGE7_PAIRWISE):
        key = unordered_pair(row["population1"], row["population2"])
        wp = (row["window_id"], key[0], key[1])
        if wp in seen_window_pair:
            continue
        seen_window_pair.add(wp)
        counts[key][row["pair_class"]] += 1
    return {key: dict(counter) for key, counter in counts.items()}


def aggregate_matrices() -> tuple[list[dict[str, object]], dict[str, dict[tuple[str, str], float]], dict[str, dict[tuple[str, str], int]]]:
    by_pair: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: {"all": [], "outside": []})
    for row in read_tsv(STAGE7_PAIRWISE):
        if row["region_class"] == "boundary":
            continue
        key = unordered_pair(row["population1"], row["population2"])
        value = float(row["mrca_time"])
        by_pair[key]["all"].append(value)
        if row["region_class"] == "outside":
            by_pair[key]["outside"].append(value)
    rows = []
    matrices = {"all": {}, "outside": {}, "all_mean": {}, "outside_mean": {}}
    counts = {"all": {}, "outside": {}}
    for key in sorted(by_pair):
        all_values = by_pair[key]["all"]
        outside_values = by_pair[key]["outside"]
        if not all_values or not outside_values:
            raise ValueError(f"Missing aggregate values for pair {key}")
        median_all = statistics.median(all_values)
        median_outside = statistics.median(outside_values)
        mean_all = statistics.fmean(all_values)
        mean_outside = statistics.fmean(outside_values)
        matrices["all"][key] = median_all
        matrices["outside"][key] = median_outside
        matrices["all_mean"][key] = mean_all
        matrices["outside_mean"][key] = mean_outside
        counts["all"][key] = len(all_values)
        counts["outside"][key] = len(outside_values)
        rows.append({
            "population1": key[0],
            "population2": key[1],
            "n_all": len(all_values),
            "median_all": median_all,
            "mean_all": mean_all,
            "n_outside": len(outside_values),
            "median_outside": median_outside,
            "mean_outside": mean_outside,
            "delta_pairwise_median": median_all - median_outside,
            "relative_change": median_all / median_outside if median_outside > 1e-12 else math.nan,
        })
    return rows, matrices, counts


def fit_node_ages(root: RNode, matrix: dict[tuple[str, str], float], topology: str, treatment: str, root_edge: str, objective: str = "L2") -> FitResult:
    taxa = sorted(tip_order(root))
    expected_pairs = {unordered_pair(a, b) for a, b in combinations(taxa, 2)}
    if set(matrix) != expected_pairs:
        missing = expected_pairs - set(matrix)
        extra = set(matrix) - expected_pairs
        raise ValueError(f"Matrix/topology pair mismatch: missing={missing}, extra={extra}")
    internals = internal_nodes(root)
    mrca_by_pair = {pair: mrca_node(root, pair[0], pair[1]).id for pair in expected_pairs}
    grouped = defaultdict(list)
    for pair, value in matrix.items():
        grouped[mrca_by_pair[pair]].append(float(value))
    fit_grouped = grouped
    if objective == "L1_median_target_approx":
        fit_grouped = defaultdict(list, {node_id: [statistics.median(values)] for node_id, values in grouped.items()})
    ages = tree_isotonic_l2(root, fit_grouped)
    for node in all_rnodes(root):
        node.age = 0.0 if node.is_tip else ages[node.id]
    residual_rows = []
    obs = []
    pred = []
    class_counts = pair_class_counts()
    for pair in sorted(matrix):
        mrca = mrca_node(root, pair[0], pair[1])
        predicted = ages[mrca.id]
        observed = float(matrix[pair])
        obs.append(observed)
        pred.append(predicted)
        counts = class_counts.get(pair, {})
        residual_rows.append({
            "topology": topology,
            "treatment": treatment,
            "population1": pair[0],
            "population2": pair[1],
            "observed_mrca": observed,
            "projected_mrca": predicted,
            "residual": observed - predicted,
            "abs_residual": abs(observed - predicted),
            "mrca_node_id": mrca.id,
            "same_arrangement_window_count": counts.get("same_arrangement", 0),
            "opposite_arrangement_window_count": counts.get("opposite_arrangement", 0),
            "unknown_window_count": counts.get("unknown", 0),
        })
    obs_arr = np.array(obs)
    pred_arr = np.array(pred)
    residual = obs_arr - pred_arr
    corr = float(np.corrcoef(obs_arr, pred_arr)[0, 1]) if np.std(obs_arr) > 0 and np.std(pred_arr) > 0 else math.nan
    return FitResult(
        topology=topology,
        treatment=treatment,
        objective_type=objective,
        root_edge=root_edge,
        rooted_root=root,
        ages=ages,
        objective=float(np.sum(residual * residual)),
        rmse=float(np.sqrt(np.mean(residual * residual))),
        mae=float(np.mean(np.abs(residual))),
        corr=corr,
        max_residual_abs=float(np.max(np.abs(residual))),
        residual_rows=residual_rows,
    )


def tree_isotonic_l2(root: RNode, grouped: dict[str, list[float]]) -> dict[str, float]:
    """Equal-weight least-squares isotonic regression on a rooted tree.

    Each internal node has a scalar target equal to the mean aggregate MRCA time
    for population pairs whose MRCA maps to that node, with weight equal to the
    number of such pairs. Adjacent violating parent/child blocks are pooled until
    every parent block age is at least every child block age. This is the
    deterministic pool-adjacent-violators solution for the rooted-tree order used
    here.
    """
    internals = internal_nodes(root)
    parent_of = {node.id: node.parent.id for node in internals if node.parent is not None and not node.parent.is_tip}
    blocks: dict[str, dict[str, object]] = {}
    node_block: dict[str, str] = {}
    for node in internals:
        values = grouped.get(node.id, [])
        weight = float(len(values))
        value = statistics.fmean(values) if values else 0.0
        blocks[node.id] = {"members": {node.id}, "weight": weight, "weighted_sum": weight * value, "value": value}
        node_block[node.id] = node.id

    changed = True
    while changed:
        changed = False
        for child_id, parent_id in list(parent_of.items()):
            child_block = node_block[child_id]
            parent_block = node_block[parent_id]
            if child_block == parent_block:
                continue
            child_value = float(blocks[child_block]["value"])
            parent_value = float(blocks[parent_block]["value"])
            if parent_value + 1e-12 < child_value:
                parent_members = set(blocks[parent_block]["members"])
                child_members = set(blocks[child_block]["members"])
                merged_members = parent_members | child_members
                weight = float(blocks[parent_block]["weight"]) + float(blocks[child_block]["weight"])
                weighted_sum = float(blocks[parent_block]["weighted_sum"]) + float(blocks[child_block]["weighted_sum"])
                value = weighted_sum / weight if weight > 0 else max(parent_value, child_value)
                new_id = min(merged_members)
                del blocks[parent_block]
                del blocks[child_block]
                blocks[new_id] = {"members": merged_members, "weight": weight, "weighted_sum": weighted_sum, "value": max(0.0, value)}
                for member in merged_members:
                    node_block[member] = new_id
                changed = True
                break
    return {node.id: max(0.0, float(blocks[node_block[node.id]]["value"])) for node in internals}


def choose_root(topology_name: str, nodes: dict[int, UNode], outside_matrix: dict[tuple[str, str], float]) -> tuple[tuple[int, int], list[dict[str, object]]]:
    rows = []
    best_edge = None
    best_loss = math.inf
    seen_labels = set()
    for edge in edge_list(nodes):
        label = root_edge_label(nodes, edge)
        if label in seen_labels:
            raise ValueError(f"Duplicate candidate root edge label: {label}")
        seen_labels.add(label)
        root = rooted_from_edge(nodes, edge)
        fit = fit_node_ages(root, outside_matrix, topology_name, "outside_root_search", label)
        loss = fit.objective
        if loss < best_loss - 1e-14:
            best_loss = loss
            best_edge = edge
        rows.append({"topology": topology_name, "candidate_root_edge": label, "loss_outside": loss, "selected": False})
    if best_edge is None:
        raise ValueError("No root edge candidates")
    selected_label = root_edge_label(nodes, best_edge)
    for row in rows:
        row["selected"] = row["candidate_root_edge"] == selected_label
    return best_edge, rows


def branch_length(parent: RNode, child: RNode) -> float:
    return max(0.0, parent.age - child.age)


def to_newick(node: RNode, parent: RNode | None = None) -> str:
    length = "" if parent is None else f":{branch_length(parent, node):.10g}"
    if node.is_tip:
        return f"{node.name}{length}"
    children = ",".join(to_newick(child, node) for child in node.children)
    label = f"{node.id}_age_{node.age:.6g}"
    return f"({children}){label}{length}"


def write_projected_tree(path: Path, fit: FitResult) -> None:
    path.write_text(to_newick(fit.rooted_root) + ";\n")


def validate_projected(root: RNode, expected_taxa: list[str]) -> None:
    if sorted(tip_order(root)) != sorted(expected_taxa):
        raise ValueError("Projected tree tip set mismatch")
    for node in all_rnodes(root):
        if node.age < -1e-10:
            raise ValueError("Negative node age")
        for child in node.children:
            if node.age + 1e-10 < child.age:
                raise ValueError("Parent age below child age")
            if branch_length(node, child) < -1e-10:
                raise ValueError("Negative branch length")
    depths = []
    def walk(node: RNode, depth: float):
        if node.is_tip:
            depths.append(depth)
        for child in node.children:
            walk(child, depth + branch_length(node, child))
    walk(root, 0.0)
    if max(depths) - min(depths) > 1e-7:
        raise ValueError("Projected tree is not ultrametric")


def topology_signature(root: RNode) -> set[tuple[str, ...]]:
    return {clade_taxa(node) for node in internal_nodes(root)}


def node_age_shift_rows(outside: FitResult, all_fit: FitResult) -> list[dict[str, object]]:
    outside_by_clade = {clade_taxa(node): node for node in internal_nodes(outside.rooted_root)}
    all_by_clade = {clade_taxa(node): node for node in internal_nodes(all_fit.rooted_root)}
    rows = []
    for idx, clade in enumerate(sorted(outside_by_clade, key=lambda c: (len(c), c)), start=1):
        o = outside_by_clade[clade]
        a = all_by_clade[clade]
        o_parent = o.parent
        a_parent = a.parent
        rows.append({
            "node_id": f"N{idx:02d}",
            "clade_taxa": "+".join(clade),
            "n_descendants": len(clade),
            "age_outside": o.age,
            "age_all": a.age,
            "delta_age": a.age - o.age,
            "relative_change": a.age / o.age if o.age > 1e-12 else math.nan,
            "branch_length_outside": branch_length(o_parent, o) if o_parent else math.nan,
            "branch_length_all": branch_length(a_parent, a) if a_parent else math.nan,
        })
    return sorted(rows, key=lambda r: abs(float(r["delta_age"])), reverse=True)


def summary_rows(topology_name: str, outside: FitResult, all_fit: FitResult, shifts: list[dict[str, object]]) -> list[dict[str, object]]:
    deltas = [float(r["delta_age"]) for r in shifts]
    outside_ages = [float(r["age_outside"]) for r in shifts]
    all_ages = [float(r["age_all"]) for r in shifts]
    corr = float(np.corrcoef(outside_ages, all_ages)[0, 1]) if np.std(outside_ages) > 0 and np.std(all_ages) > 0 else math.nan
    return [{
        "topology": topology_name,
        "root_edge": outside.root_edge,
        "root_age_outside": outside.rooted_root.age,
        "root_age_all": all_fit.rooted_root.age,
        "root_age_shift": all_fit.rooted_root.age - outside.rooted_root.age,
        "root_relative_change": all_fit.rooted_root.age / outside.rooted_root.age if outside.rooted_root.age > 1e-12 else math.nan,
        "median_absolute_internal_node_age_shift": statistics.median(abs(d) for d in deltas),
        "maximum_absolute_node_age_shift": max(abs(d) for d in deltas),
        "mean_signed_node_age_shift": statistics.fmean(deltas),
        "n_nodes_older_with_all": sum(d > 1e-10 for d in deltas),
        "n_nodes_younger_with_all": sum(d < -1e-10 for d in deltas),
        "n_nodes_unchanged": sum(abs(d) <= 1e-10 for d in deltas),
        "node_age_correlation": corr,
    }]


def fit_rows(fits: list[FitResult]) -> list[dict[str, object]]:
    return [{
        "topology": f.topology,
        "treatment": f.treatment,
        "objective_type": f.objective_type,
        "root_edge": f.root_edge,
        "objective": f.objective,
        "rmse": f.rmse,
        "mae": f.mae,
        "correlation": f.corr,
        "maximum_abs_residual": f.max_residual_abs,
    } for f in fits]


def draw_tree(ax, root: RNode, title: str, y_positions: dict[str, float] | None = None) -> dict[str, float]:
    taxa = tip_order(root)
    if y_positions is None:
        y_positions = {taxon: i for i, taxon in enumerate(reversed(taxa))}
    node_y = {}
    def assign_y(node: RNode) -> float:
        if node.is_tip:
            node_y[node.id] = y_positions[str(node.name)]
        else:
            vals = [assign_y(child) for child in node.children]
            node_y[node.id] = sum(vals) / len(vals)
        return node_y[node.id]
    assign_y(root)
    root_age = root.age
    def x(node: RNode) -> float:
        return root_age - node.age
    def draw_node(node: RNode):
        for child in node.children:
            ax.plot([x(node), x(child)], [node_y[node.id], node_y[child.id]], color="black", lw=0.9)
            ax.plot([x(child), x(child)], [node_y[child.id], node_y[child.id]], color="black", lw=0.9)
            draw_node(child)
    draw_node(root)
    for taxon, y in y_positions.items():
        ax.text(root_age * 1.01, y, taxon.replace("Gadmor_", ""), va="center", fontsize=7)
    ax.set_ylim(-1, len(taxa))
    ax.set_xlim(0, root_age * 1.25 if root_age > 0 else 1)
    ax.invert_xaxis()
    ax.set_yticks([])
    ax.set_xlabel("projected time before present")
    ax.set_title(title)
    return y_positions


def plot_summary(outside: FitResult, all_fit: FitResult, shifts: list[dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1], wspace=0.32, hspace=0.35)
    ax_a = fig.add_subplot(gs[0, 0])
    y = draw_tree(ax_a, outside.rooted_root, "A. T0 topology + outside-window times")
    ax_b = fig.add_subplot(gs[1, 0])
    draw_tree(ax_b, all_fit.rooted_root, "B. same T0 topology + all-window times", y)
    ax_c = fig.add_subplot(gs[0, 1])
    x = [float(r["age_outside"]) for r in shifts]
    yv = [float(r["age_all"]) for r in shifts]
    maxv = max(x + yv)
    ax_c.scatter(x, yv, color="#4c78a8")
    ax_c.plot([0, maxv], [0, maxv], color="black", lw=0.8)
    ax_c.set_xlabel("outside-only fitted node age")
    ax_c.set_ylabel("all-window fitted node age")
    ax_c.set_title("C. internal node ages")
    ax_d = fig.add_subplot(gs[1, 1])
    ordered = sorted(shifts, key=lambda r: float(r["delta_age"]))
    labels = [f"N{i+1}" for i in range(len(ordered))]
    vals = [float(r["delta_age"]) for r in ordered]
    colors = ["#e45756" if v > 0 else "#4c78a8" for v in vals]
    ax_d.barh(labels, vals, color=colors)
    ax_d.axvline(0, color="black", lw=0.8)
    ax_d.set_xlabel("all - outside fitted age")
    ax_d.set_title("D. node-age shifts")
    fig.suptitle("Atlantic cod Stage 8 summary-tree time projection", y=0.98)
    fig.savefig(FIG_SUMMARY_PDF, bbox_inches="tight")
    fig.savefig(FIG_SUMMARY_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(matrix_rows: list[dict[str, object]]) -> None:
    taxa = sorted({r["population1"] for r in matrix_rows} | {r["population2"] for r in matrix_rows})
    idx = {taxon: i for i, taxon in enumerate(taxa)}
    arr = np.zeros((len(taxa), len(taxa)))
    for row in matrix_rows:
        i = idx[str(row["population1"])]
        j = idx[str(row["population2"])]
        value = float(row["delta_pairwise_median"])
        arr[i, j] = value
        arr[j, i] = value
    fig, ax = plt.subplots(figsize=(8, 7))
    vmax = max(abs(float(np.min(arr))), abs(float(np.max(arr))))
    im = ax.imshow(arr, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    short = [t.replace("Gadmor_", "") for t in taxa]
    ax.set_xticks(range(len(taxa)), short, rotation=90, fontsize=7)
    ax.set_yticks(range(len(taxa)), short, fontsize=7)
    ax.set_title("Pairwise median MRCA shift: all windows - outside windows")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("SNAPP workflow time units")
    fig.tight_layout()
    fig.savefig(FIG_HEATMAP_PDF)
    fig.savefig(FIG_HEATMAP_PNG, dpi=300)
    plt.close(fig)


def write_text_outputs(summary_t0: dict[str, object], summary_t1: dict[str, object], fit_table: list[dict[str, object]], shift_rows: list[dict[str, object]]) -> None:
    top = sorted(shift_rows, key=lambda r: abs(float(r["delta_age"])), reverse=True)[:5]
    top_txt = "; ".join(f"{r['node_id']} ({r['n_descendants']} taxa) Δ={float(r['delta_age']):.3g}" for r in top)
    METHODS_TEXT.write_text(textwrap.dedent("""
    As a post-freeze Stage-8 extension, we projected the aggregate SNAPP window-level MRCA-time signal onto fixed ASTRAL population-tree topologies. The primary topology was the Stage-6 all-window ASTRAL tree (`T0_all.inferred.nwk`); the outside-window ASTRAL tree (`T1_outside.inferred.nwk`) was used as a secondary topology robustness check. ASTRAL annotations, q-values, local posterior probabilities, CU/SU lengths, and existing numerical branch lengths were ignored. Only the taxon topology was used.

    For each unordered population pair, we aggregated Stage-7 MRCA times across valid non-boundary windows using the median, separately for all focal windows and for fully outside windows. The unrooted ASTRAL topology was rooted by enumerating every possible edge, fitting node ages to the outside-only pairwise MRCA matrix under equal-weight constrained least squares, and selecting the root edge with the smallest outside-only projection loss. That root was then frozen and reused for both the outside-only and all-window projections. Node ages were estimated by minimizing squared residuals between the aggregate pairwise MRCA matrix and the age of each pair's MRCA on the fixed rooted topology, subject to nonnegative node ages and parent ages at least as old as child ages. As a robust sensitivity analysis, we also used a documented L1-style approximation that replaces each node's least-squares target with the median of pairwise MRCA values mapping to that node and then applies the same tree-monotonic constraints. This is a deterministic projection/sensitivity analysis, not ASTRAL divergence-time estimation.
    """).strip() + "\n")
    RESULTS_TEXT.write_text(textwrap.dedent(f"""
    Stage 8 asked whether including inversion-linked windows changes the fitted temporal scale when published SNAPP MRCA times are projected onto the same fixed ASTRAL topology. On the primary T0_all topology, the selected root edge was `{summary_t0['root_edge']}`. The outside-only projected root age was {float(summary_t0['root_age_outside']):.3g}, whereas the all-window projected root age was {float(summary_t0['root_age_all']):.3g}, a shift of {float(summary_t0['root_age_shift']):.3g}. Across internal nodes, the median absolute age shift was {float(summary_t0['median_absolute_internal_node_age_shift']):.3g}, the maximum absolute shift was {float(summary_t0['maximum_absolute_node_age_shift']):.3g}, and the mean signed shift was {float(summary_t0['mean_signed_node_age_shift']):.3g}. Nodes shifted older/younger/unchanged with all windows as {summary_t0['n_nodes_older_with_all']}/{summary_t0['n_nodes_younger_with_all']}/{summary_t0['n_nodes_unchanged']}.

    The largest internal-node shifts on T0_all were: {top_txt}. Projection fit should be interpreted directly because the aggregate MRCA matrix is not guaranteed to be exactly representable by a single ultrametric tree; fit diagnostics are reported separately for both treatments and both topologies. The secondary T1_outside topology gave a root-age shift of {float(summary_t1['root_age_shift']):.3g} and a median absolute internal-node shift of {float(summary_t1['median_absolute_internal_node_age_shift']):.3g}, providing a topology-robustness check without changing the primary T0_all inference.
    """).strip() + "\n")
    FIGURE_CAPTION.write_text(textwrap.dedent("""
    Figure. Atlantic cod Stage-8 summary-tree divergence-time projection. Panel A shows the fixed Stage-6 T0_all ASTRAL topology with node ages fitted from outside-only aggregate SNAPP pairwise MRCA times. Panel B shows the identical rooted topology with node ages fitted from all valid non-boundary focal windows. ASTRAL CU/SU branch lengths and support annotations were ignored; only the topology was used. Panel C compares fitted internal-node ages between outside-only and all-window projections. Panel D shows all-window minus outside-only fitted node-age shifts. The analysis is a deterministic projection of published SNAPP local-time signal onto a fixed ASTRAL topology, not ASTRAL divergence-time estimation.
    """).strip() + "\n")
    lines = [
        "# Atlantic cod Stage-8 report",
        "",
        "Stage 8 is a post-freeze projection/sensitivity analysis. It uses Stage-7 SNAPP MRCA times and fixed Stage-6 ASTRAL topologies, and it does not rerun SNAPP or ASTRAL4.",
        "",
        "## Primary T0_all summary",
        "",
        f"Selected root edge: `{summary_t0['root_edge']}`.",
        f"Root age outside/all/shift: {float(summary_t0['root_age_outside']):.6g} / {float(summary_t0['root_age_all']):.6g} / {float(summary_t0['root_age_shift']):.6g}.",
        f"Median absolute internal-node shift: {float(summary_t0['median_absolute_internal_node_age_shift']):.6g}.",
        f"Maximum absolute node-age shift: {float(summary_t0['maximum_absolute_node_age_shift']):.6g}.",
        f"Mean signed node-age shift: {float(summary_t0['mean_signed_node_age_shift']):.6g}.",
        f"Older/younger/unchanged nodes with all windows: {summary_t0['n_nodes_older_with_all']} / {summary_t0['n_nodes_younger_with_all']} / {summary_t0['n_nodes_unchanged']}.",
        "",
        "## Projection fit",
        "",
        "| topology | treatment | objective | RMSE | MAE | correlation | max abs residual |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in fit_table:
        lines.append(f"| {row['topology']} | {row['treatment']} | {row['objective_type']} | {float(row['rmse']):.6g} | {float(row['mae']):.6g} | {float(row['correlation']):.6g} | {float(row['maximum_abs_residual']):.6g} |")
    lines.extend([
        "",
        "## Secondary T1_outside summary",
        "",
        f"Root age outside/all/shift: {float(summary_t1['root_age_outside']):.6g} / {float(summary_t1['root_age_all']):.6g} / {float(summary_t1['root_age_shift']):.6g}.",
        f"Median absolute internal-node shift: {float(summary_t1['median_absolute_internal_node_age_shift']):.6g}.",
    ])
    REPORT.write_text("\n".join(lines) + "\n")


def write_manifest(frozen: dict[str, str], outputs: list[Path], summary_t0: dict[str, object], root_rows: list[dict[str, object]]) -> None:
    manifest = {
        "stage": 8,
        "analysis": "ASTRAL summary-tree divergence-time projection",
        "post_freeze_extension": True,
        "snapp_rerun": False,
        "astral4_rerun": False,
        "primary_topology": str(T0_TOPOLOGY.relative_to(REPO_ROOT)),
        "secondary_topology": str(T1_TOPOLOGY.relative_to(REPO_ROOT)),
        "selected_primary_root_edge": summary_t0["root_edge"],
        "candidate_root_edges_tested_T0": sum(1 for r in root_rows if r["topology"] == "T0_all"),
        "software_versions": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
        },
        "frozen_input_checksums": frozen,
        "output_checksums": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in outputs if path.exists()},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def run_stage8() -> None:
    frozen = verify_frozen_inputs()
    matrix_rows, matrices, counts = aggregate_matrices()
    write_tsv(PAIRWISE_MATRICES, matrix_rows, ["population1", "population2", "n_all", "median_all", "mean_all", "n_outside", "median_outside", "mean_outside", "delta_pairwise_median", "relative_change"])
    taxa = sorted({r["population1"] for r in matrix_rows} | {r["population2"] for r in matrix_rows})
    if len(taxa) != 12:
        raise ValueError(f"Expected exactly 12 cod populations, found {len(taxa)}")

    all_root_rows = []
    all_fits = []
    summaries = {}
    primary_shifts = []
    residuals = []
    tree_outputs = {}
    for topology_name, path in (("T0_all", T0_TOPOLOGY), ("T1_outside", T1_TOPOLOGY)):
        nodes = parse_topology(path)
        if taxa_of_graph(nodes) != taxa:
            raise ValueError(f"Tip set mismatch for {topology_name}")
        selected_edge, root_rows = choose_root(topology_name, nodes, matrices["outside"])
        all_root_rows.extend(root_rows)
        root_label = root_edge_label(nodes, selected_edge)
        outside_root = rooted_from_edge(nodes, selected_edge)
        all_root = rooted_from_edge(nodes, selected_edge)
        outside_fit = fit_node_ages(outside_root, matrices["outside"], topology_name, "outside", root_label)
        all_fit = fit_node_ages(all_root, matrices["all"], topology_name, "all", root_label)
        outside_l1 = fit_node_ages(rooted_from_edge(nodes, selected_edge), matrices["outside"], topology_name, "outside_L1_approx", root_label, objective="L1_median_target_approx")
        all_l1 = fit_node_ages(rooted_from_edge(nodes, selected_edge), matrices["all"], topology_name, "all_L1_approx", root_label, objective="L1_median_target_approx")
        validate_projected(outside_fit.rooted_root, taxa)
        validate_projected(all_fit.rooted_root, taxa)
        if topology_signature(outside_fit.rooted_root) != topology_signature(all_fit.rooted_root):
            raise ValueError(f"Projected all/outside topologies differ for {topology_name}")
        shifts = node_age_shift_rows(outside_fit, all_fit)
        summary = summary_rows(topology_name, outside_fit, all_fit, shifts)[0]
        summaries[topology_name] = summary
        all_fits.extend([outside_fit, all_fit, outside_l1, all_l1])
        residuals.extend(outside_fit.residual_rows)
        residuals.extend(all_fit.residual_rows)
        if topology_name == "T0_all":
            primary_shifts = shifts
            tree_outputs[(topology_name, "outside")] = T0_OUTSIDE_TREE
            tree_outputs[(topology_name, "all")] = T0_ALL_TREE
            write_projected_tree(T0_OUTSIDE_TREE, outside_fit)
            write_projected_tree(T0_ALL_TREE, all_fit)
        else:
            tree_outputs[(topology_name, "outside")] = T1_OUTSIDE_TREE
            tree_outputs[(topology_name, "all")] = T1_ALL_TREE
            write_projected_tree(T1_OUTSIDE_TREE, outside_fit)
            write_projected_tree(T1_ALL_TREE, all_fit)

    fit_table = fit_rows(all_fits)
    write_tsv(ROOTING_AUDIT, all_root_rows, ["topology", "candidate_root_edge", "loss_outside", "selected"])
    write_tsv(NODE_AGE_SHIFTS, primary_shifts, ["node_id", "clade_taxa", "n_descendants", "age_outside", "age_all", "delta_age", "relative_change", "branch_length_outside", "branch_length_all"])
    write_tsv(SUMMARY, [summaries["T0_all"], summaries["T1_outside"]], ["topology", "root_edge", "root_age_outside", "root_age_all", "root_age_shift", "root_relative_change", "median_absolute_internal_node_age_shift", "maximum_absolute_node_age_shift", "mean_signed_node_age_shift", "n_nodes_older_with_all", "n_nodes_younger_with_all", "n_nodes_unchanged", "node_age_correlation"])
    write_tsv(PROJECTION_FIT, fit_table, ["topology", "treatment", "objective_type", "root_edge", "objective", "rmse", "mae", "correlation", "maximum_abs_residual"])
    write_tsv(PROJECTION_RESIDUALS, residuals, ["topology", "treatment", "population1", "population2", "observed_mrca", "projected_mrca", "residual", "abs_residual", "mrca_node_id", "same_arrangement_window_count", "opposite_arrangement_window_count", "unknown_window_count"])
    t0_out = next(f for f in all_fits if f.topology == "T0_all" and f.treatment == "outside")
    t0_all = next(f for f in all_fits if f.topology == "T0_all" and f.treatment == "all")
    plot_summary(t0_out, t0_all, primary_shifts)
    plot_heatmap(matrix_rows)
    write_text_outputs(summaries["T0_all"], summaries["T1_outside"], fit_table, primary_shifts)
    outputs = [
        PAIRWISE_MATRICES, ROOTING_AUDIT, NODE_AGE_SHIFTS, SUMMARY, PROJECTION_FIT, PROJECTION_RESIDUALS,
        T0_OUTSIDE_TREE, T0_ALL_TREE, T1_OUTSIDE_TREE, T1_ALL_TREE,
        METHODS_TEXT, RESULTS_TEXT, FIGURE_CAPTION, REPORT,
        FIG_SUMMARY_PDF, FIG_SUMMARY_PNG, FIG_HEATMAP_PDF, FIG_HEATMAP_PNG,
    ]
    write_manifest(frozen, outputs, summaries["T0_all"], all_root_rows)


class Stage8Tests(unittest.TestCase):
    def test_pair_symmetry_and_zero_diagonal(self):
        rows, matrices, _counts = aggregate_matrices()
        taxa = sorted({r["population1"] for r in rows} | {r["population2"] for r in rows})
        self.assertEqual(len(taxa), 12)
        self.assertEqual(len(rows), len(taxa) * (len(taxa) - 1) // 2)
        for a in taxa:
            self.assertNotIn((a, a), matrices["all"])
        self.assertEqual(len(set((r["population1"], r["population2"]) for r in rows)), len(rows))

    def test_topology_tip_sets(self):
        rows, _matrices, _counts = aggregate_matrices()
        taxa = sorted({r["population1"] for r in rows} | {r["population2"] for r in rows})
        self.assertEqual(taxa_of_graph(parse_topology(T0_TOPOLOGY)), taxa)
        self.assertEqual(taxa_of_graph(parse_topology(T1_TOPOLOGY)), taxa)

    def test_candidate_root_edges_once(self):
        nodes = parse_topology(T0_TOPOLOGY)
        labels = [root_edge_label(nodes, edge) for edge in edge_list(nodes)]
        self.assertEqual(len(labels), len(set(labels)))
        self.assertEqual(len(labels), 2 * 12 - 3)

    def test_deterministic_rooting_and_reuse(self):
        _rows, matrices, _counts = aggregate_matrices()
        nodes = parse_topology(T0_TOPOLOGY)
        edge1, _ = choose_root("T0_all", nodes, matrices["outside"])
        edge2, _ = choose_root("T0_all", nodes, matrices["outside"])
        self.assertEqual(root_edge_label(nodes, edge1), root_edge_label(nodes, edge2))

    def test_projected_tree_constraints(self):
        _rows, matrices, _counts = aggregate_matrices()
        nodes = parse_topology(T0_TOPOLOGY)
        taxa = taxa_of_graph(nodes)
        edge, _ = choose_root("T0_all", nodes, matrices["outside"])
        fit = fit_node_ages(rooted_from_edge(nodes, edge), matrices["outside"], "T0_all", "outside", root_edge_label(nodes, edge))
        validate_projected(fit.rooted_root, taxa)
        self.assertTrue(all(age >= 0 for age in fit.ages.values()))

    def test_all_outside_topologies_identical(self):
        _rows, matrices, _counts = aggregate_matrices()
        nodes = parse_topology(T0_TOPOLOGY)
        edge, _ = choose_root("T0_all", nodes, matrices["outside"])
        out = fit_node_ages(rooted_from_edge(nodes, edge), matrices["outside"], "T0_all", "outside", root_edge_label(nodes, edge))
        allfit = fit_node_ages(rooted_from_edge(nodes, edge), matrices["all"], "T0_all", "all", root_edge_label(nodes, edge))
        self.assertEqual(topology_signature(out.rooted_root), topology_signature(allfit.rooted_root))

    def test_projection_does_not_use_astral_sulength_as_time(self):
        text = Path(__file__).read_text()
        self.assertIn("SULength", T0_TOPOLOGY.read_text())
        self.assertNotIn("SU" + "Length=", text)
        self.assertNotIn("CU" + "Length=", text)

    def test_stage7_physical_formula_corrected(self):
        rows = {row["lg"]: row for row in read_tsv(STAGE7_TESTS)}
        expected = {"LG01": 1/46, "LG02": 1/75, "LG07": 1/89, "LG12": 5/56}
        for lg, value in expected.items():
            self.assertAlmostEqual(float(rows[lg]["physical_p"]), value, places=6)

    def test_stage7_circular_unchanged_and_frozen(self):
        verify_frozen_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage8Tests))
        raise SystemExit(0 if result.wasSuccessful() else 1)
    run_stage8()
