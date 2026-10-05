#!/usr/bin/env python3
"""Shared Stage 2 helpers for house-mouse t-complex controls."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import random
import re
import shutil
import subprocess
import zipfile
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio import Phylo

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "house_mouse_t_complex"
DATA_ROOT = REPO_ROOT / "data" / "house_mouse_t_complex"
PROCESSED = DATA_ROOT / "processed"
METADATA = DATA_ROOT / "metadata"
RESULTS = EMPIRICAL_ROOT / "results"
FIGURES = EMPIRICAL_ROOT / "figures"
STAGE2_DATA = PROCESSED / "stage2"
STAGE2_RESULTS = RESULTS / "stage2"

Q_SPECIES = "Q_SPECIES"
Q_T_ALT = "Q_T_ALT"
Q_OTHER = "Q_OTHER"
Q_UNRESOLVED = "Q_UNRESOLVED"

Q_SPLITS = {
    Q_SPECIES: frozenset(
        [
            frozenset(["Mus_musculus_musculus", "Mus_musculus_castaneus"]),
            frozenset(["Mus_musculus_domesticus", "Mus_spretus"]),
        ]
    ),
    Q_T_ALT: frozenset(
        [
            frozenset(["Mus_musculus_domesticus", "Mus_musculus_musculus"]),
            frozenset(["Mus_musculus_castaneus", "Mus_spretus"]),
        ]
    ),
    Q_OTHER: frozenset(
        [
            frozenset(["Mus_musculus_domesticus", "Mus_musculus_castaneus"]),
            frozenset(["Mus_musculus_musculus", "Mus_spretus"]),
        ]
    ),
}

SPECIES_LABELS = {
    "Mus musculus domesticus": "Mus_musculus_domesticus",
    "Mus musculus musculus": "Mus_musculus_musculus",
    "Mus musculus castaneus": "Mus_musculus_castaneus",
    "Mus spretus": "Mus_spretus",
}


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
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
            writer.writerow({field: format_value(row.get(field, "")) for field in fields})


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.10g}"
    if isinstance(value, (list, tuple, set, frozenset)):
        return ",".join(str(x) for x in value)
    return str(value)


def load_stage1() -> tuple[list[dict[str, str]], list[str], list[dict[str, str]]]:
    meta = read_tsv(PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv")
    trees = [line.strip() for line in (PROCESSED / "house_mouse_t_complex_ml_5kb.tre").read_text().splitlines() if line.strip()]
    mapping = read_tsv(METADATA / "tip_mapping.tsv")
    if len(meta) != len(trees):
        raise RuntimeError(f"Metadata/tree count mismatch: {len(meta)} != {len(trees)}")
    return meta, trees, mapping


def mapping_by_tip(mapping: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["tree_tip"]: row for row in mapping}


def tips_by_species(mapping: list[dict[str, str]], treatment: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for row in mapping:
        confidence = row.get("mapping_confidence", "")
        if confidence not in {"strong", "exact"}:
            continue
        subspecies = row["subspecies"]
        t_status = row["t_status"]
        include = treatment == "ALL_TIPS"
        if treatment == "STANDARD_ONLY":
            include = t_status in {"standard_noncarrier", "outgroup_not_t_haplotype"}
        elif treatment == "T_ONLY":
            include = t_status in {"pseudo-t_haplotype", "outgroup_not_t_haplotype"}
        if include:
            out[subspecies].append(row["tree_tip"])
    for key in out:
        out[key].sort()
    return dict(out)


def all_required_groups_present(groups: dict[str, list[str]]) -> bool:
    required = ["Mus musculus domesticus", "Mus musculus musculus", "Mus musculus castaneus", "Mus spretus"]
    return all(groups.get(group) for group in required)


def tree_adjacency_and_tips(newick: str) -> tuple[list[list[int]], dict[str, int]]:
    tree = Phylo.read(io.StringIO(newick), "newick")
    adjacency: list[list[int]] = []
    tip_index: dict[str, int] = {}
    node_id: dict[int, int] = {}

    def get_id(clade) -> int:
        key = id(clade)
        if key not in node_id:
            node_id[key] = len(adjacency)
            adjacency.append([])
        return node_id[key]

    for parent in tree.find_clades(order="level"):
        pid = get_id(parent)
        if parent.is_terminal():
            tip_index[parent.name] = pid
        for child in parent.clades:
            cid = get_id(child)
            adjacency[pid].append(cid)
            adjacency[cid].append(pid)
    return adjacency, tip_index


def topological_distance_matrix(newick: str, tips: list[str]) -> np.ndarray:
    adjacency, tip_to_node = tree_adjacency_and_tips(newick)
    missing = [tip for tip in tips if tip not in tip_to_node]
    if missing:
        raise RuntimeError(f"Tree missing expected tips: {missing[:5]}")
    n = len(tips)
    dist = np.zeros((n, n), dtype=np.int16)
    for i, tip in enumerate(tips):
        start = tip_to_node[tip]
        seen = {start: 0}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for nb in adjacency[node]:
                if nb not in seen:
                    seen[nb] = seen[node] + 1
                    queue.append(nb)
        for j, other in enumerate(tips):
            dist[i, j] = seen[tip_to_node[other]]
    return dist


def quartet_counts_for_tree(newick: str, groups: dict[str, list[str]]) -> dict[str, int]:
    required = ["Mus musculus domesticus", "Mus musculus musculus", "Mus musculus castaneus", "Mus spretus"]
    if not all_required_groups_present(groups):
        return {Q_SPECIES: 0, Q_T_ALT: 0, Q_OTHER: 0, Q_UNRESOLVED: 0}
    tips = sorted({tip for group in required for tip in groups[group]})
    index = {tip: i for i, tip in enumerate(tips)}
    dist = topological_distance_matrix(newick, tips)
    d = np.array([index[t] for t in groups["Mus musculus domesticus"]], dtype=np.int32)
    m = np.array([index[t] for t in groups["Mus musculus musculus"]], dtype=np.int32)
    c = np.array([index[t] for t in groups["Mus musculus castaneus"]], dtype=np.int32)
    s = np.array([index[t] for t in groups["Mus spretus"]], dtype=np.int32)
    dd, mm, cc, ss = np.meshgrid(d, m, c, s, indexing="ij")
    species = dist[mm, cc] + dist[dd, ss]
    t_alt = dist[dd, mm] + dist[cc, ss]
    other = dist[dd, cc] + dist[mm, ss]
    stacked = np.stack([species, t_alt, other], axis=0)
    minv = stacked.min(axis=0)
    ties = (stacked == minv).sum(axis=0)
    return {
        Q_SPECIES: int(((species == minv) & (ties == 1)).sum()),
        Q_T_ALT: int(((t_alt == minv) & (ties == 1)).sum()),
        Q_OTHER: int(((other == minv) & (ties == 1)).sum()),
        Q_UNRESOLVED: int((ties != 1).sum()),
    }


def score_to_row(counts: dict[str, int]) -> dict[str, object]:
    total = sum(counts.values())
    resolved = counts[Q_SPECIES] + counts[Q_T_ALT] + counts[Q_OTHER]
    if total == 0:
        return {
            "n_resolved_quartets": 0,
            "n_unresolved_quartets": 0,
            "q_species": float("nan"),
            "q_t_alt": float("nan"),
            "q_other": float("nan"),
            "q_unresolved": float("nan"),
            "delta_species_alt": float("nan"),
        }
    q_species = counts[Q_SPECIES] / total
    q_t_alt = counts[Q_T_ALT] / total
    return {
        "n_resolved_quartets": resolved,
        "n_unresolved_quartets": counts[Q_UNRESOLVED],
        "q_species": q_species,
        "q_t_alt": q_t_alt,
        "q_other": counts[Q_OTHER] / total,
        "q_unresolved": counts[Q_UNRESOLVED] / total,
        "delta_species_alt": q_species - q_t_alt,
    }


def canonical_split(parts: Iterable[Iterable[str]]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    a, b = [tuple(sorted(p)) for p in parts]
    return tuple(sorted([a, b]))  # type: ignore[return-value]


def four_taxon_topology_from_newick(newick: str) -> str:
    tree = Phylo.read(io.StringIO(newick), "newick")
    taxa = sorted(t.name for t in tree.get_terminals())
    if set(taxa) != set(SPECIES_LABELS.values()):
        raise ValueError(f"Expected four species labels, observed {taxa}")
    dist = {}
    for a in taxa:
        for b in taxa:
            if a < b:
                dist[(a, b)] = len(tree.trace(a, b))
    pairings = {
        Q_SPECIES: dist[tuple(sorted(["Mus_musculus_musculus", "Mus_musculus_castaneus"]))]
        + dist[tuple(sorted(["Mus_musculus_domesticus", "Mus_spretus"]))],
        Q_T_ALT: dist[tuple(sorted(["Mus_musculus_domesticus", "Mus_musculus_musculus"]))]
        + dist[tuple(sorted(["Mus_musculus_castaneus", "Mus_spretus"]))],
        Q_OTHER: dist[tuple(sorted(["Mus_musculus_domesticus", "Mus_musculus_castaneus"]))]
        + dist[tuple(sorted(["Mus_musculus_musculus", "Mus_spretus"]))],
    }
    minv = min(pairings.values())
    winners = [k for k, v in pairings.items() if v == minv]
    return winners[0] if len(winners) == 1 else Q_UNRESOLVED


def tree_bipartitions(newick: str) -> set[tuple[tuple[str, ...], tuple[str, ...]]]:
    tree = Phylo.read(io.StringIO(newick), "newick")
    all_tips = frozenset(t.name for t in tree.get_terminals())
    splits = set()
    for clade in tree.find_clades(order="postorder"):
        tips = frozenset(t.name for t in clade.get_terminals())
        if 1 < len(tips) < len(all_tips) - 1:
            other = all_tips - tips
            splits.add(canonical_split([tips, other]))
    return splits


def rf_distance(newick_a: str, newick_b: str) -> int:
    a = tree_bipartitions(newick_a)
    b = tree_bipartitions(newick_b)
    return len(a - b) + len(b - a)


def run_astral(astral: Path, tree_file: Path, map_file: Path, output: Path, threads: int = 2, seed: int = 20261005) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    log = output.with_suffix(".log")
    cmd = [
        str(astral),
        "-i",
        str(tree_file),
        "-a",
        str(map_file),
        "-o",
        str(output),
        "-u",
        "2",
        "-t",
        str(threads),
        "--length",
        "CULength",
        "--seed",
        str(seed),
    ]
    with log.open("w") as handle:
        subprocess.run(cmd, cwd=REPO_ROOT, stdout=handle, stderr=subprocess.STDOUT, check=True)
    return parse_astral_result(output, log) | {"command": " ".join(cmd), "tree_file": rel(tree_file), "map_file": rel(map_file)}


def find_astral() -> Path:
    env = os.environ.get("ASTRAL4_BIN")
    candidates = []
    if env:
        candidates.append(Path(env))
    found = shutil.which("astral4")
    if found:
        candidates.append(Path(found))
    candidates.append(REPO_ROOT / "external" / "ASTER" / "bin" / "astral4")
    for path in candidates:
        if path.exists() and os.access(path, os.X_OK):
            return path
    raise RuntimeError("No astral4 executable found in ASTRAL4_BIN, PATH, or external/ASTER/bin/astral4")


def parse_astral_result(tree_path: Path, log_path: Path) -> dict[str, object]:
    newick = tree_path.read_text().strip()
    log = log_path.read_text()
    result: dict[str, object] = {"newick": newick, "tree_output": rel(tree_path), "log": rel(log_path)}
    for line in log.splitlines():
        if line.startswith("#Genetrees:"):
            result["n_gene_trees"] = int(line.split(":", 1)[1].strip())
        elif line.startswith("#Species:"):
            result["n_species"] = int(line.split(":", 1)[1].strip())
        elif line.startswith("Final Tree:"):
            result["final_tree_unannotated"] = line.split(":", 1)[1].strip()
        elif line.startswith("Score:"):
            result["astral_score"] = line.split(":", 1)[1].strip()
    annotations = parse_first_annotation(newick)
    result.update(annotations)
    try:
        result["topology"] = four_taxon_topology_from_newick(str(result.get("final_tree_unannotated", newick)))
    except Exception:
        result["topology"] = "NA"
    return result


def parse_first_annotation(newick: str) -> dict[str, object]:
    match = re.search(r"\'\[([^\]]+)\]\'", newick)
    out: dict[str, object] = {"CU_length": "NA", "SU_length": "NA", "localPP": "NA", "q1": "NA", "q2": "NA", "q3": "NA"}
    if not match:
        return out
    for part in match.group(1).split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        if key == "CULength":
            out["CU_length"] = value
        elif key == "SULength":
            out["SU_length"] = value
        elif key == "localPP":
            out["localPP"] = value
        elif key in {"q1", "q2", "q3"}:
            out[key] = value
    return out


def prune_newick(newick: str, keep: set[str]) -> str:
    tree = Phylo.read(io.StringIO(newick), "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep:
            tree.prune(terminal)
    out = io.StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip()


def write_tree_file(path: Path, newicks: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(n.strip() for n in newicks) + "\n")


def write_map(path: Path, mapping_rows: list[dict[str, str]], tips: Iterable[str], field: str) -> None:
    by_tip = mapping_by_tip(mapping_rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for tip in sorted(tips):
            value = by_tip[tip][field].replace(" ", "_")
            handle.write(f"{tip}\t{value}\n")


def deterministic_balanced_samples(mapping_rows: list[dict[str, str]], n_replicates: int, seed: int) -> list[dict[str, object]]:
    rng = random.Random(seed)
    by_species_status: dict[tuple[str, str], list[str]] = defaultdict(list)
    spretus = []
    for row in mapping_rows:
        tip = row["tree_tip"]
        species = row["subspecies"]
        status = row["t_status"]
        if species == "Mus spretus":
            spretus.append(tip)
        elif status == "standard_noncarrier":
            by_species_status[(species, "standard")].append(tip)
        elif status == "pseudo-t_haplotype":
            by_species_status[(species, "t")].append(tip)
    for key in by_species_status:
        by_species_status[key].sort()
    spretus.sort()
    species_order = ["Mus musculus castaneus", "Mus musculus domesticus", "Mus musculus musculus"]
    rows = []
    for rep in range(n_replicates):
        selected_spretus = rng.sample(spretus, 3)
        standard = list(selected_spretus)
        t_only = list(selected_spretus)
        mixed = list(selected_spretus)
        for species in species_order:
            st = rng.sample(by_species_status[(species, "standard")], 3)
            tt = rng.sample(by_species_status[(species, "t")], 3)
            standard.extend(st)
            t_only.extend(tt)
            mixed.extend(st + tt)
        rows.extend(
            [
                {"replicate": rep, "treatment": "B0_STANDARD_MATCHED", "selected_tips": sorted(standard)},
                {"replicate": rep, "treatment": "B1_T_MATCHED", "selected_tips": sorted(t_only)},
                {"replicate": rep, "treatment": "B2_MIXED_BALANCED", "selected_tips": sorted(mixed)},
            ]
        )
    return rows


def spatial_thin_records(meta: list[dict[str, str]], spacing: int, phase: int) -> list[int]:
    bins: dict[int, tuple[int, int]] = {}
    for i, row in enumerate(meta):
        start = int(row["start_bp"])
        shifted = start - 5_000_000 - phase
        if shifted < 0:
            continue
        b = shifted // spacing
        if b not in bins or start < bins[b][0]:
            bins[b] = (start, i)
    return [idx for _start, idx in sorted(bins.values())]


def aggregate_quartet_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    by_treatment: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_treatment[row["treatment"]].append(row)
    out = []
    for treatment, subset in sorted(by_treatment.items()):
        n = len(subset)
        out.append(
            {
                "treatment": treatment,
                "n_loci": n,
                "mean_q_species": sum(float(r["q_species"]) for r in subset) / n,
                "mean_q_t_alt": sum(float(r["q_t_alt"]) for r in subset) / n,
                "mean_q_other": sum(float(r["q_other"]) for r in subset) / n,
                "mean_q_unresolved": sum(float(r["q_unresolved"]) for r in subset) / n,
                "mean_delta_species_alt": sum(float(r["delta_species_alt"]) for r in subset) / n,
            }
        )
    return out
