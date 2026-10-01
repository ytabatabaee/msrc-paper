#!/usr/bin/env python3
"""Retrieve and normalize the published fire-ant TWISST/local-tree dataset.

Stage 1 creates topology-neutral canonical tables for the 213 published windows.
It verifies Stage-0 frozen provenance, retrieves the three commit-pinned
topology-bearing upstream files when needed, validates Newick/weight structure,
and writes normalized handoff tables. It does not classify topologies into the
Stage-0 focal quartet and does not calculate focal support statistics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import tempfile
import unittest
import urllib.request
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"
RAW_TWISST_DIR = DATA_ROOT / "raw" / "upstream" / "Topology weighting" / "results" / "2021-08-07-twisst"

STAGE0_FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
STAGE0_COORDINATES = DATA_ROOT / "raw" / "upstream" / "Topology weighting" / "results" / "window_coordinates"
STAGE0_PROVENANCE = DATA_ROOT / "metadata" / "upstream_repository_provenance.json"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"

INPUT_TREE = RAW_TWISST_DIR / "input.tree"
TOPOLOGIES_OUTPUT = RAW_TWISST_DIR / "topologies_output.trees"
WEIGHTS_OUTPUT = RAW_TWISST_DIR / "weights_output.csv"

WINDOW_TREES = DATA_ROOT / "processed" / "stage1_window_trees.tsv"
TWISST_TOPOLOGIES = DATA_ROOT / "processed" / "stage1_twisst_topologies.tsv"
TWISST_WEIGHTS_SPARSE = DATA_ROOT / "processed" / "stage1_twisst_weights_sparse.tsv"
TWISST_WEIGHT_SUMMARY = DATA_ROOT / "processed" / "stage1_twisst_weight_summary.tsv"
WINDOW_INDEX = DATA_ROOT / "processed" / "stage1_window_index.tsv"

TOPOLOGY_CROSSCHECK = EMPIRICAL_ROOT / "results" / "stage1_topology_crosscheck.tsv"
SOURCE_PROVENANCE = EMPIRICAL_ROOT / "results" / "stage1_source_provenance.md"
VALIDATION_SUMMARY = EMPIRICAL_ROOT / "results" / "stage1_validation_summary.tsv"
STAGE1_MANIFEST = EMPIRICAL_ROOT / "results" / "stage1_manifest.json"

UPSTREAM_REPO = "https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression"
UPSTREAM_COMMIT = "bdc7823941680fa560e61b6467af9f77950f64b0"
UPSTREAM_BRANCH = "master"
DOWNLOAD_DATE = "2026-09-30"
EXPECTED_STAGE0_FOCAL_SHA = "5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6"
EXPECTED_STAGE0_COORDINATE_SHA = "c46139b1a3a33c8e22f6b1ee85cbe08d2d7dfd23b47f0263d1fab579e157d791"
EXPECTED_REGION_ORDER = ["chr1", "chr16A", "chr16B", "chr16_supergene"]
EXPECTED_REGION_COUNTS = {"chr1": 117, "chr16A": 42, "chr16B": 2, "chr16_supergene": 52}
EXPECTED_GROUPS = [
    "geminata",
    "saevissima",
    "pusillignis",
    "invicta/macdonaghi_Sb",
    "invicta/macdonaghi_SB",
    "richteri_Sb",
    "richteri_SB",
]
EXPECTED_TOPOLOGY_IDS = [f"topo{i}" for i in range(1, 946)]
RAW_TOPOLOGY_FILES = {
    INPUT_TREE: {
        "original_filename": "Topology weighting/results/2021-08-07-twisst/input.tree",
        "url": f"https://raw.githubusercontent.com/wurmlab/2021-fire-ant-social-supergene-introgression/{UPSTREAM_COMMIT}/Topology%20weighting/results/2021-08-07-twisst/input.tree",
        "purpose": "Published 213 full individual-level RAxML local trees, one non-empty Newick tree per upstream window.",
    },
    TOPOLOGIES_OUTPUT: {
        "original_filename": "Topology weighting/results/2021-08-07-twisst/topologies_output.trees",
        "url": f"https://raw.githubusercontent.com/wurmlab/2021-fire-ant-social-supergene-introgression/{UPSTREAM_COMMIT}/Topology%20weighting/results/2021-08-07-twisst/topologies_output.trees",
        "purpose": "Published inventory of 945 possible seven-group TWISST topologies.",
    },
    WEIGHTS_OUTPUT: {
        "original_filename": "Topology weighting/results/2021-08-07-twisst/weights_output.csv",
        "url": f"https://raw.githubusercontent.com/wurmlab/2021-fire-ant-social-supergene-introgression/{UPSTREAM_COMMIT}/Topology%20weighting/results/2021-08-07-twisst/weights_output.csv",
        "purpose": "Published TWISST weights for each of 945 grouped topologies across 213 local-tree windows.",
    },
}
OUTPUT_PATHS = [
    WINDOW_TREES,
    TWISST_TOPOLOGIES,
    TWISST_WEIGHTS_SPARSE,
    TWISST_WEIGHT_SUMMARY,
    WINDOW_INDEX,
    TOPOLOGY_CROSSCHECK,
    SOURCE_PROVENANCE,
    VALIDATION_SUMMARY,
    STAGE1_MANIFEST,
]
FORBIDDEN_STAGE1_OUTPUT_FIELDS = {
    "q_s",
    "q_h",
    "q_3",
    "d",
    "delta_d",
    "species_split_support",
    "haplotype_split_support",
    "third_split_support",
    "quartet_class",
    "species_support",
    "haplotype_support",
    "third_support",
}
TSV_LINETERMINATOR = "\n"


class Stage1Error(RuntimeError):
    """Raised when Stage-1 validation fails."""


@dataclass(frozen=True)
class ParsedTree:
    tips: tuple[str, ...]
    splits: tuple[str, ...]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def dec_to_text(value: Decimal) -> str:
    if value == value.to_integral_value():
        return str(value.quantize(Decimal(1)))
    return format(value.normalize(), "f")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else str(row.get(field, "")) for field in fields})


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def verify_stage0() -> None:
    observed_focal = sha256(STAGE0_FOCAL_QUARTET)
    if observed_focal != EXPECTED_STAGE0_FOCAL_SHA:
        raise Stage1Error(f"Stage-0 focal-quartet checksum changed: {observed_focal}")
    observed_coord = sha256(STAGE0_COORDINATES)
    if observed_coord != EXPECTED_STAGE0_COORDINATE_SHA:
        raise Stage1Error(f"Stage-0 coordinate checksum changed: {observed_coord}")
    data = json.loads(STAGE0_PROVENANCE.read_text())
    expected = {
        "repository_url": UPSTREAM_REPO,
        "default_branch": UPSTREAM_BRANCH,
        "head_commit_sha": UPSTREAM_COMMIT,
    }
    for key, value in expected.items():
        if data.get(key) != value:
            raise Stage1Error(f"Stage-0 upstream provenance {key} changed: {data.get(key)!r}")


def download_if_needed(skip_download: bool = False) -> None:
    RAW_TWISST_DIR.mkdir(parents=True, exist_ok=True)
    for path, meta in RAW_TOPOLOGY_FILES.items():
        if path.exists() and path.stat().st_size > 0:
            continue
        if skip_download:
            raise Stage1Error(f"missing raw file in --skip-download mode: {path.relative_to(REPO_ROOT)}")
        urllib.request.urlretrieve(meta["url"], path)


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.i = 0

    def parse(self) -> object:
        if not self.text:
            raise Stage1Error("empty Newick tree")
        node = self.parse_subtree()
        self.skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ";":
            raise Stage1Error("Newick tree missing terminal semicolon")
        self.i += 1
        self.skip_ws()
        if self.i != len(self.text):
            raise Stage1Error("unexpected content after Newick semicolon")
        return node

    def skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def parse_subtree(self) -> object:
        self.skip_ws()
        if self.i >= len(self.text):
            raise Stage1Error("unexpected end of Newick tree")
        if self.text[self.i] == "(":
            self.i += 1
            children = [self.parse_subtree()]
            while True:
                self.skip_ws()
                if self.i >= len(self.text):
                    raise Stage1Error("unterminated internal node")
                if self.text[self.i] == ",":
                    self.i += 1
                    children.append(self.parse_subtree())
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise Stage1Error(f"unexpected Newick character {self.text[self.i]!r}")
            self.consume_optional_node_label()
            self.consume_optional_branch_length()
            if len(children) < 2:
                raise Stage1Error("internal node has fewer than two children")
            return tuple(children)
        label = self.consume_label()
        if not label:
            raise Stage1Error("empty tip label")
        self.consume_optional_branch_length()
        return label

    def consume_label(self) -> str:
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in "(),:;":
            self.i += 1
        return self.text[start:self.i].strip()

    def consume_optional_node_label(self) -> None:
        self.skip_ws()
        while self.i < len(self.text) and self.text[self.i] not in ",:();":
            self.i += 1

    def consume_optional_branch_length(self) -> None:
        self.skip_ws()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            start = self.i
            while self.i < len(self.text) and self.text[self.i] not in ",();":
                self.i += 1
            if self.i == start:
                raise Stage1Error("empty branch length")


def collect_tips(node: object, tips: list[str]) -> None:
    if isinstance(node, str):
        tips.append(node)
        return
    for child in node:
        collect_tips(child, tips)


def tree_edges_and_leaves(node: object) -> tuple[dict[int, set[int]], dict[int, str]]:
    graph: dict[int, set[int]] = {}
    leaves: dict[int, str] = {}
    next_id = 0

    def add(current: object, parent: int | None = None) -> int:
        nonlocal next_id
        node_id = next_id
        next_id += 1
        graph.setdefault(node_id, set())
        if parent is not None:
            graph[node_id].add(parent)
            graph[parent].add(node_id)
        if isinstance(current, str):
            leaves[node_id] = current
        else:
            for child in current:
                add(child, node_id)
        return node_id

    add(node)
    return graph, leaves


def unrooted_splits(node: object) -> tuple[str, ...]:
    graph, leaves_by_id = tree_edges_and_leaves(node)
    all_labels = frozenset(leaves_by_id.values())
    leaf_ids = set(leaves_by_id)
    splits: set[str] = set()

    for a, neighbors in graph.items():
        for b in neighbors:
            if a > b:
                continue
            stack = [a]
            seen = {b}
            component_ids: set[int] = set()
            while stack:
                current = stack.pop()
                if current in seen:
                    continue
                seen.add(current)
                component_ids.add(current)
                stack.extend(graph[current] - seen)
            left = frozenset(leaves_by_id[i] for i in component_ids & leaf_ids)
            if len(left) <= 1 or len(left) >= len(all_labels) - 1:
                continue
            right = all_labels - left
            left_text = ",".join(sorted(left))
            right_text = ",".join(sorted(right))
            splits.add("|".join(sorted([left_text, right_text])))
    return tuple(sorted(splits))


def parse_newick(text: str) -> ParsedTree:
    root = NewickParser(text).parse()
    tips: list[str] = []
    collect_tips(root, tips)
    if not tips:
        raise Stage1Error("Newick tree has no tips")
    if len(tips) != len(set(tips)):
        duplicates = sorted(label for label, count in Counter(tips).items() if count > 1)
        raise Stage1Error(f"Newick tree has duplicate tip labels: {', '.join(duplicates)}")
    return ParsedTree(tips=tuple(tips), splits=unrooted_splits(root))


def parse_coordinates() -> list[dict[str, str]]:
    fields, rows = read_tsv(STAGE0_COORDINATES)
    required = ["chrom", "start", "end", "mid", "region", "window", "gene1", "gene2", "gene3", "gene4"]
    missing = [field for field in required if field not in fields]
    if missing:
        raise Stage1Error(f"coordinate table missing fields: {', '.join(missing)}")
    if len(rows) != 213:
        raise Stage1Error(f"expected 213 coordinate rows, observed {len(rows)}")
    expected_sequence = []
    for region in EXPECTED_REGION_ORDER:
        expected_sequence.extend([region] * EXPECTED_REGION_COUNTS[region])
    observed_sequence = [row["region"] for row in rows]
    if observed_sequence != expected_sequence:
        raise Stage1Error("coordinate region sequence does not match documented upstream construction order")
    for region in EXPECTED_REGION_ORDER:
        region_windows = [int(row["window"]) for row in rows if row["region"] == region]
        expected_windows = list(range(1, 1 + 4 * len(region_windows), 4))
        if region_windows != expected_windows:
            raise Stage1Error(f"{region} window IDs do not follow documented four-gene increments")
    return rows


def parse_local_trees() -> tuple[list[dict[str, object]], dict[str, object]]:
    lines = nonempty_lines(INPUT_TREE)
    if len(lines) != 213:
        raise Stage1Error(f"expected 213 local trees, observed {len(lines)}")
    rows: list[dict[str, object]] = []
    n_tips_values: list[int] = []
    tip_set_hashes: set[str] = set()
    for i, line in enumerate(lines, 1):
        parsed = parse_newick(line)
        if len(parsed.tips) < 4:
            raise Stage1Error(f"local tree {i} has fewer than four tips")
        tip_hash = sha256_text("\n".join(sorted(parsed.tips)) + "\n")
        tree_hash = sha256_text(line + "\n")
        n_tips_values.append(len(parsed.tips))
        tip_set_hashes.add(tip_hash)
        rows.append(
            {
                "window_index": i,
                "tree_newick": line,
                "n_tips": len(parsed.tips),
                "tip_set_sha256": tip_hash,
                "tree_sha256": tree_hash,
            }
        )
    return rows, {
        "min_tree_tips": min(n_tips_values),
        "max_tree_tips": max(n_tips_values),
        "n_distinct_tip_sets": len(tip_set_hashes),
    }


def parse_twisst_topologies() -> tuple[list[dict[str, object]], dict[str, tuple[str, ...]]]:
    lines = nonempty_lines(TOPOLOGIES_OUTPUT)
    if len(lines) != 945:
        raise Stage1Error(f"expected 945 TWISST topologies, observed {len(lines)}")
    rows: list[dict[str, object]] = []
    signatures: dict[str, tuple[str, ...]] = {}
    seen_signatures: set[tuple[str, ...]] = set()
    expected_group_set = set(EXPECTED_GROUPS)
    for i, line in enumerate(lines, 1):
        topology_id = f"topo{i}"
        parsed = parse_newick(line)
        if len(parsed.tips) != 7:
            raise Stage1Error(f"{topology_id} has {len(parsed.tips)} groups, expected 7")
        if set(parsed.tips) != expected_group_set:
            raise Stage1Error(f"{topology_id} group set does not match expected seven TWISST groups")
        if parsed.splits in seen_signatures:
            raise Stage1Error(f"{topology_id} duplicates a previous labeled unrooted topology")
        seen_signatures.add(parsed.splits)
        signatures[topology_id] = parsed.splits
        rows.append(
            {
                "topology_id": topology_id,
                "topology_number": i,
                "tree_newick": line,
                "n_groups": 7,
                "group_set_sha256": sha256_text("\n".join(sorted(parsed.tips)) + "\n"),
                "tree_sha256": sha256_text(line + "\n"),
            }
        )
    return rows, signatures


def parse_weights() -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, tuple[str, ...]], dict[str, object]]:
    comment_topologies: dict[str, tuple[str, ...]] = {}
    header: list[str] | None = None
    data_rows: list[list[str]] = []
    comment_re = re.compile(r"^#(topo[0-9]+)\s+(.+)$")
    with WEIGHTS_OUTPUT.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                match = comment_re.match(line)
                if not match:
                    raise Stage1Error(f"malformed weight comment line: {line[:80]}")
                topology_id, tree = match.groups()
                comment_topologies[topology_id] = parse_newick(tree).splits
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
            else:
                data_rows.append(parts)
    if header is None:
        raise Stage1Error("weights_output.csv has no header")
    if header != EXPECTED_TOPOLOGY_IDS:
        raise Stage1Error("weights_output.csv topology columns are not exactly topo1-topo945")
    if sorted(comment_topologies, key=lambda item: int(item.removeprefix("topo"))) != EXPECTED_TOPOLOGY_IDS:
        raise Stage1Error("weights_output.csv comment topology IDs are not exactly topo1-topo945")
    if len(data_rows) != 213:
        raise Stage1Error(f"expected 213 weight rows, observed {len(data_rows)}")

    sparse_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    totals: list[Decimal] = []
    for window_index, values in enumerate(data_rows, 1):
        if len(values) != 945:
            raise Stage1Error(f"weight row {window_index} has {len(values)} columns, expected 945")
        total = Decimal(0)
        n_nonzero = 0
        max_weight = Decimal(0)
        for topology_id, raw_value in zip(header, values, strict=True):
            try:
                weight = Decimal(raw_value)
            except InvalidOperation as exc:
                raise Stage1Error(f"non-numeric weight in row {window_index}, {topology_id}: {raw_value!r}") from exc
            if weight < 0:
                raise Stage1Error(f"negative weight in row {window_index}, {topology_id}")
            total += weight
            if weight > 0:
                n_nonzero += 1
                max_weight = max(max_weight, weight)
                sparse_rows.append(
                    {
                        "window_index": window_index,
                        "topology_id": topology_id,
                        "raw_weight": dec_to_text(weight),
                    }
                )
        if total <= 0:
            raise Stage1Error(f"weight row {window_index} has non-positive total weight")
        totals.append(total)
        summary_rows.append(
            {
                "window_index": window_index,
                "total_raw_weight": dec_to_text(total),
                "n_nonzero_topologies": n_nonzero,
                "max_raw_weight": dec_to_text(max_weight),
            }
        )
    return sparse_rows, summary_rows, comment_topologies, {
        "min_total_weight": min(totals),
        "max_total_weight": max(totals),
    }


def crosscheck_topology_definitions(
    inventory_signatures: dict[str, tuple[str, ...]], comment_signatures: dict[str, tuple[str, ...]]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for topology_id in EXPECTED_TOPOLOGY_IDS:
        inventory = inventory_signatures.get(topology_id)
        comment = comment_signatures.get(topology_id)
        status = "match" if inventory == comment else "mismatch"
        rows.append(
            {
                "topology_id": topology_id,
                "inventory_signature_sha256": sha256_text("\n".join(inventory or ()) + "\n"),
                "weight_comment_signature_sha256": sha256_text("\n".join(comment or ()) + "\n"),
                "status": status,
                "notes": "" if status == "match" else "unrooted labeled split signatures differ",
            }
        )
    return rows


def build_window_tables(
    coordinates: list[dict[str, str]],
    local_tree_rows: list[dict[str, object]],
    weight_summary_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    tree_by_index = {int(row["window_index"]): row for row in local_tree_rows}
    summary_by_index = {int(row["window_index"]): row for row in weight_summary_rows}
    window_tree_rows: list[dict[str, object]] = []
    index_rows: list[dict[str, object]] = []
    for idx, coord in enumerate(coordinates, 1):
        tree = tree_by_index[idx]
        summary = summary_by_index[idx]
        base = {
            "window_index": idx,
            "chrom": coord["chrom"],
            "start": coord["start"],
            "end": coord["end"],
            "mid": coord["mid"],
            "region": coord["region"],
            "upstream_window_id": coord["window"],
            "gene1": coord["gene1"],
            "gene2": coord["gene2"],
            "gene3": coord["gene3"],
            "gene4": coord["gene4"],
        }
        window_tree_rows.append(
            {
                **base,
                "tree_newick": tree["tree_newick"],
                "n_tips": tree["n_tips"],
                "tip_set_sha256": tree["tip_set_sha256"],
                "tree_sha256": tree["tree_sha256"],
            }
        )
        index_rows.append(
            {
                **base,
                "tree_sha256": tree["tree_sha256"],
                "total_raw_weight": summary["total_raw_weight"],
                "n_nonzero_topologies": summary["n_nonzero_topologies"],
            }
        )
    return window_tree_rows, index_rows


def update_download_manifest() -> None:
    fields, rows = read_tsv(DOWNLOAD_MANIFEST)
    existing = {row["local_path"]: row for row in rows}
    for path, meta in RAW_TOPOLOGY_FILES.items():
        rel = path.relative_to(REPO_ROOT).as_posix()
        existing[rel] = {
            "source": "GitHub wurmlab/2021-fire-ant-social-supergene-introgression",
            "original_filename": meta["original_filename"],
            "url_or_doi": meta["url"],
            "download_date": DOWNLOAD_DATE,
            "size_bytes": str(path.stat().st_size),
            "sha256": sha256(path),
            "purpose": f"{meta['purpose']} Upstream commit: {UPSTREAM_COMMIT}.",
            "local_path": rel,
        }
    ordered_rows = []
    seen = set()
    for row in rows:
        key = row["local_path"]
        ordered_rows.append(existing[key])
        seen.add(key)
    for path in RAW_TOPOLOGY_FILES:
        rel = path.relative_to(REPO_ROOT).as_posix()
        if rel not in seen:
            ordered_rows.append(existing[rel])
    write_tsv(DOWNLOAD_MANIFEST, ordered_rows, fields)


def write_source_provenance(metrics: dict[str, object], raw_hashes: dict[str, str]) -> None:
    raw_lines = [f"- `{path}`: `{digest}`" for path, digest in raw_hashes.items()]
    text = [
        "# Fire ants chromosome 16 Stage 1 source provenance",
        "",
        "Stage 1 retrieves and normalizes the published local-tree/TWISST dataset. It does not classify topologies into the Stage-0 focal quartet and does not calculate focal support statistics.",
        "",
        "## Frozen upstream source",
        "",
        f"- Repository: `{UPSTREAM_REPO}`",
        f"- Branch recorded in Stage 0: `{UPSTREAM_BRANCH}`",
        f"- Frozen commit: `{UPSTREAM_COMMIT}`",
        "",
        "## Raw downloaded files",
        "",
        *raw_lines,
        "",
        "## Published construction order",
        "",
        "The upstream `Topology weighting/README.md` constructs `results/2021-08-07-twisst/input.tree` by appending windows from `chr1`, then `chr16A`, then `chr16B`, then the chromosome-16 supergene windows. Within each region, trees are appended in increasing upstream window order from `results/window_coordinates`.",
        "",
        "Stage 1 maps coordinate rows to `input.tree` lines by this documented construction order. It does not infer or optimize the correspondence from topology.",
        "",
        "## Object distinction",
        "",
        "- `input.tree` contains 213 full individual-level RAxML local trees.",
        "- `topologies_output.trees` contains the 945 possible seven-group topologies used by TWISST.",
        "- `weights_output.csv` gives the TWISST weight of those grouped topologies for each local tree/window.",
        "",
        "The 945 topology inventory is not a set of 945 loci or 945 gene trees.",
        "",
        "## Validation counts",
        "",
        f"- Coordinate rows: `{metrics['coordinate_rows']}`",
        f"- Local trees: `{metrics['local_tree_rows']}`",
        f"- TWISST grouped topologies: `{metrics['twisst_topologies']}`",
        f"- Weight rows: `{metrics['weight_rows']}`",
        f"- Weight columns: `{metrics['weight_columns']}`",
        f"- Seven-group vocabulary: `{', '.join(EXPECTED_GROUPS)}`",
        f"- Minimum local-tree tips: `{metrics['min_tree_tips']}`",
        f"- Maximum local-tree tips: `{metrics['max_tree_tips']}`",
        f"- Distinct local-tree tip sets: `{metrics['n_distinct_tip_sets']}`",
        "",
        "If tip sets differ among windows, the canonical table records per-window `n_tips` and `tip_set_sha256`; no pruning or imputation is performed in Stage 1.",
    ]
    SOURCE_PROVENANCE.write_text("\n".join(text) + "\n")


def write_validation_summary(metrics: dict[str, object]) -> None:
    rows = [{"metric": key, "value": value, "notes": ""} for key, value in metrics.items()]
    write_tsv(VALIDATION_SUMMARY, rows, ["metric", "value", "notes"])


def write_manifest(metrics: dict[str, object], raw_hashes: dict[str, str]) -> None:
    processed = {path.relative_to(REPO_ROOT).as_posix(): sha256(path) for path in OUTPUT_PATHS if path != STAGE1_MANIFEST}
    data = {
        "stage": "Stage 1",
        "description": "Retrieve, validate, and normalize published local-tree/TWISST dataset without focal-quartet support calculations.",
        "stage0_focal_quartet_checksum": EXPECTED_STAGE0_FOCAL_SHA,
        "stage0_coordinate_checksum": EXPECTED_STAGE0_COORDINATE_SHA,
        "upstream_repository": UPSTREAM_REPO,
        "upstream_branch": UPSTREAM_BRANCH,
        "upstream_commit": UPSTREAM_COMMIT,
        "raw_topology_file_checksums": raw_hashes,
        "processing_script": {
            "path": "empirical/fire_ants_chr16/scripts/01_retrieve_normalize_twisst.py",
            "sha256": sha256(Path(__file__)),
        },
        "processed_output_checksums": processed,
        "expected_row_counts": {
            "coordinate_rows": 213,
            "local_tree_rows": 213,
            "twisst_topologies": 945,
            "weight_rows": 213,
            "weight_columns": 945,
            "region_counts": EXPECTED_REGION_COUNTS,
        },
        "observed_metrics": {key: str(value) for key, value in metrics.items()},
        "inputs": {
            "stage0_focal_quartet": STAGE0_FOCAL_QUARTET.relative_to(REPO_ROOT).as_posix(),
            "stage0_coordinates": STAGE0_COORDINATES.relative_to(REPO_ROOT).as_posix(),
            "input_tree": INPUT_TREE.relative_to(REPO_ROOT).as_posix(),
            "topologies_output": TOPOLOGIES_OUTPUT.relative_to(REPO_ROOT).as_posix(),
            "weights_output": WEIGHTS_OUTPUT.relative_to(REPO_ROOT).as_posix(),
        },
        "outputs": [path.relative_to(REPO_ROOT).as_posix() for path in OUTPUT_PATHS],
        "prohibited_stage1_products": sorted(FORBIDDEN_STAGE1_OUTPUT_FIELDS),
    }
    STAGE1_MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate_output_schemas() -> None:
    for path in [WINDOW_TREES, TWISST_TOPOLOGIES, TWISST_WEIGHTS_SPARSE, TWISST_WEIGHT_SUMMARY, WINDOW_INDEX, TOPOLOGY_CROSSCHECK, VALIDATION_SUMMARY]:
        fields, _ = read_tsv(path)
        lowered = {field.lower() for field in fields}
        forbidden = lowered & FORBIDDEN_STAGE1_OUTPUT_FIELDS
        if forbidden:
            raise Stage1Error(f"{path.relative_to(REPO_ROOT)} contains prohibited Stage-1 fields: {', '.join(sorted(forbidden))}")


def run_stage1(skip_download: bool = False) -> dict[str, object]:
    verify_stage0()
    download_if_needed(skip_download=skip_download)
    update_download_manifest()

    coordinates = parse_coordinates()
    local_tree_rows, tree_metrics = parse_local_trees()
    topology_rows, inventory_signatures = parse_twisst_topologies()
    sparse_rows, weight_summary_rows, comment_signatures, weight_metrics = parse_weights()
    crosscheck_rows = crosscheck_topology_definitions(inventory_signatures, comment_signatures)
    mismatches = sum(1 for row in crosscheck_rows if row["status"] != "match")
    if mismatches:
        raise Stage1Error(f"topology definition cross-check found {mismatches} mismatches")

    window_tree_rows, window_index_rows = build_window_tables(coordinates, local_tree_rows, weight_summary_rows)
    region_counts = Counter(row["region"] for row in coordinates)
    sparse_totals: dict[int, Decimal] = Counter()
    for row in sparse_rows:
        sparse_totals[int(row["window_index"])] += Decimal(str(row["raw_weight"]))
    summary_totals = {int(row["window_index"]): Decimal(str(row["total_raw_weight"])) for row in weight_summary_rows}
    if sparse_totals != summary_totals:
        raise Stage1Error("sparse weights do not reconstruct each original row total exactly")

    write_tsv(
        WINDOW_TREES,
        window_tree_rows,
        [
            "window_index",
            "chrom",
            "start",
            "end",
            "mid",
            "region",
            "upstream_window_id",
            "gene1",
            "gene2",
            "gene3",
            "gene4",
            "tree_newick",
            "n_tips",
            "tip_set_sha256",
            "tree_sha256",
        ],
    )
    write_tsv(TWISST_TOPOLOGIES, topology_rows, ["topology_id", "topology_number", "tree_newick", "n_groups", "group_set_sha256", "tree_sha256"])
    write_tsv(TWISST_WEIGHTS_SPARSE, sparse_rows, ["window_index", "topology_id", "raw_weight"])
    write_tsv(TWISST_WEIGHT_SUMMARY, weight_summary_rows, ["window_index", "total_raw_weight", "n_nonzero_topologies", "max_raw_weight"])
    write_tsv(
        WINDOW_INDEX,
        window_index_rows,
        [
            "window_index",
            "chrom",
            "start",
            "end",
            "mid",
            "region",
            "upstream_window_id",
            "gene1",
            "gene2",
            "gene3",
            "gene4",
            "tree_sha256",
            "total_raw_weight",
            "n_nonzero_topologies",
        ],
    )
    write_tsv(TOPOLOGY_CROSSCHECK, crosscheck_rows, ["topology_id", "inventory_signature_sha256", "weight_comment_signature_sha256", "status", "notes"])

    raw_hashes = {path.relative_to(REPO_ROOT).as_posix(): sha256(path) for path in RAW_TOPOLOGY_FILES}
    metrics: dict[str, object] = {
        "coordinate_rows": len(coordinates),
        "local_tree_rows": len(local_tree_rows),
        "weight_rows": len(weight_summary_rows),
        "weight_columns": 945,
        "twisst_topologies": len(topology_rows),
        "unique_twisst_topologies": len(set(inventory_signatures.values())),
        "chr1_windows": region_counts["chr1"],
        "chr16A_windows": region_counts["chr16A"],
        "chr16B_windows": region_counts["chr16B"],
        "chr16_supergene_windows": region_counts["chr16_supergene"],
        "min_tree_tips": tree_metrics["min_tree_tips"],
        "max_tree_tips": tree_metrics["max_tree_tips"],
        "n_distinct_tip_sets": tree_metrics["n_distinct_tip_sets"],
        "min_total_weight": dec_to_text(weight_metrics["min_total_weight"]),
        "max_total_weight": dec_to_text(weight_metrics["max_total_weight"]),
        "topology_crosscheck_mismatches": mismatches,
    }
    write_validation_summary(metrics)
    write_source_provenance(metrics, raw_hashes)
    write_manifest(metrics, raw_hashes)
    validate_output_schemas()
    return metrics


class Stage1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.metrics = run_stage1(skip_download=True)

    def test_stage0_focal_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE0_FOCAL_QUARTET), EXPECTED_STAGE0_FOCAL_SHA)

    def test_upstream_commit_unchanged(self) -> None:
        data = json.loads(STAGE0_PROVENANCE.read_text())
        self.assertEqual(data["head_commit_sha"], UPSTREAM_COMMIT)

    def test_downloaded_hashes_deterministic(self) -> None:
        self.assertEqual(sha256(INPUT_TREE), "0b0d64217170610b2bd3e8dd1481138b032ef033589887b18a883099c1780e66")
        self.assertEqual(sha256(TOPOLOGIES_OUTPUT), "36a58a44eb730972bc5e4c2003fc9cdd9d1f8681c9ba2c14e053daa2706bd57c")
        self.assertEqual(sha256(WEIGHTS_OUTPUT), "617edfa9a1d0ae35ed7b0203c73e8a4eb78a6bb40bf32cdaae80ddb5b4182c7c")

    def test_coordinate_rows(self) -> None:
        self.assertEqual(len(parse_coordinates()), 213)

    def test_local_tree_count_and_parse(self) -> None:
        rows, _ = parse_local_trees()
        self.assertEqual(len(rows), 213)

    def test_no_duplicate_tips_per_local_tree(self) -> None:
        for line in nonempty_lines(INPUT_TREE):
            parsed = parse_newick(line)
            self.assertEqual(len(parsed.tips), len(set(parsed.tips)))

    def test_twisst_topology_definitions(self) -> None:
        rows, _ = parse_twisst_topologies()
        self.assertEqual(len(rows), 945)
        self.assertTrue(all(row["n_groups"] == 7 for row in rows))

    def test_grouped_topology_group_set(self) -> None:
        for line in nonempty_lines(TOPOLOGIES_OUTPUT):
            parsed = parse_newick(line)
            self.assertEqual(set(parsed.tips), set(EXPECTED_GROUPS))

    def test_topology_ids(self) -> None:
        rows, _ = parse_twisst_topologies()
        self.assertEqual([row["topology_id"] for row in rows], EXPECTED_TOPOLOGY_IDS)

    def test_grouped_topologies_unique(self) -> None:
        _, signatures = parse_twisst_topologies()
        self.assertEqual(len(signatures), len(set(signatures.values())))

    def test_weight_columns_and_rows(self) -> None:
        _, summaries, _, _ = parse_weights()
        fields = WEIGHTS_OUTPUT.read_text().splitlines()[945].split("\t")
        self.assertEqual(fields, EXPECTED_TOPOLOGY_IDS)
        self.assertEqual(len(summaries), 213)

    def test_weights_numeric_nonnegative_positive_totals(self) -> None:
        _, summaries, _, _ = parse_weights()
        for row in summaries:
            self.assertGreater(Decimal(str(row["total_raw_weight"])), 0)

    def test_weight_comment_crosscheck(self) -> None:
        _, inventory = parse_twisst_topologies()
        _, _, comments, _ = parse_weights()
        rows = crosscheck_topology_definitions(inventory, comments)
        self.assertTrue(all(row["status"] == "match" for row in rows))

    def test_coordinate_tree_weight_correspondence(self) -> None:
        self.assertEqual(self.metrics["coordinate_rows"], 213)
        self.assertEqual(self.metrics["local_tree_rows"], 213)
        self.assertEqual(self.metrics["weight_rows"], 213)

    def test_expected_region_counts(self) -> None:
        self.assertEqual(self.metrics["chr1_windows"], 117)
        self.assertEqual(self.metrics["chr16A_windows"], 42)
        self.assertEqual(self.metrics["chr16B_windows"], 2)
        self.assertEqual(self.metrics["chr16_supergene_windows"], 52)

    def test_sparse_weights_reconstruct_totals(self) -> None:
        _, sparse = read_tsv(TWISST_WEIGHTS_SPARSE)
        _, summaries = read_tsv(TWISST_WEIGHT_SUMMARY)
        totals: dict[int, Decimal] = Counter()
        for row in sparse:
            totals[int(row["window_index"])] += Decimal(row["raw_weight"])
        expected = {int(row["window_index"]): Decimal(row["total_raw_weight"]) for row in summaries}
        self.assertEqual(totals, expected)

    def test_repeated_runs_byte_stable(self) -> None:
        before = {path: sha256(path) for path in OUTPUT_PATHS}
        run_stage1(skip_download=True)
        after = {path: sha256(path) for path in OUTPUT_PATHS}
        self.assertEqual(before, after)

    def test_no_focal_quartet_support_statistic_produced(self) -> None:
        validate_output_schemas()


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage1Tests)
    with tempfile.TemporaryDirectory(prefix="fire_ants_stage1_tests_"):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true", help="regenerate outputs from existing immutable raw files")
    parser.add_argument("--run-tests", action="store_true", help="run embedded Stage-1 tests after normalization")
    args = parser.parse_args(argv)
    try:
        metrics = run_stage1(skip_download=args.skip_download)
    except Stage1Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.run_tests:
        return 0 if run_tests() else 1
    print(f"Wrote {WINDOW_INDEX.relative_to(REPO_ROOT)}")
    print(f"Wrote {TWISST_WEIGHTS_SPARSE.relative_to(REPO_ROOT)}")
    print(f"Validated {metrics['local_tree_rows']} local trees and {metrics['twisst_topologies']} TWISST topologies")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
