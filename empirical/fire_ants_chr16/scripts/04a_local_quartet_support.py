#!/usr/bin/env python3
"""Stage 4A fire-ant local quartet-support analysis.

This is the first formal unblinded local-topology stage. It classifies the
published seven-group TWISST topologies by the independently frozen focal
quartet, aggregates published TWISST weights by local window, and writes
descriptive tables and raw support figures. It performs no inferential test,
P-value, circular shift, smoothing, boundary optimization, or ASTRAL/ASTER run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import sys
import tempfile
import unittest
from collections import Counter, deque
from decimal import Decimal, getcontext
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


getcontext().prec = 50

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"

STAGE0_FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
STAGE1_MANIFEST = EMPIRICAL_ROOT / "results" / "stage1_manifest.json"
STAGE2_MANIFEST = EMPIRICAL_ROOT / "results" / "stage2_manifest.json"
STAGE2_FOCAL_PARTITION = DATA_ROOT / "processed" / "stage2_focal_partition.tsv"
STAGE3_MANIFEST = EMPIRICAL_ROOT / "results" / "stage3_manifest.json"
STAGE3_BACKGROUND_PARTITION = DATA_ROOT / "processed" / "stage3_background_partition.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"

STAGE1_TOPOLOGIES = DATA_ROOT / "processed" / "stage1_twisst_topologies.tsv"
STAGE1_WEIGHTS_SPARSE = DATA_ROOT / "processed" / "stage1_twisst_weights_sparse.tsv"
STAGE1_WEIGHT_SUMMARY = DATA_ROOT / "processed" / "stage1_twisst_weight_summary.tsv"
STAGE1_WINDOW_INDEX = DATA_ROOT / "processed" / "stage1_window_index.tsv"

TOPOLOGY_CLASSES = DATA_ROOT / "processed" / "stage4a_topology_classes.tsv"
WINDOW_SUPPORT = DATA_ROOT / "processed" / "stage4a_window_quartet_support.tsv"
REGION_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4a_region_summary.tsv"
PRIMARY_CONTRAST = EMPIRICAL_ROOT / "results" / "stage4a_primary_contrast.tsv"
REPORT = EMPIRICAL_ROOT / "results" / "stage4a_report.md"
FIGURE_CAPTION = EMPIRICAL_ROOT / "results" / "stage4a_figure_caption.md"
MANIFEST = EMPIRICAL_ROOT / "results" / "stage4a_manifest.json"

CHR16_SUPPORT_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_quartet_support.pdf"
CHR16_SUPPORT_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_quartet_support.png"
CHR1_SUPPORT_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr1_quartet_support_control.pdf"
CHR1_SUPPORT_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr1_quartet_support_control.png"
CHR16_D_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_D_track.pdf"
CHR16_D_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_D_track.png"

EXPECTED_STAGE0_FOCAL_SHA = "5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6"
EXPECTED_STAGE1_MANIFEST_SHA = "fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b"
EXPECTED_STAGE2_MANIFEST_SHA = "2b56ac7ea4a4d0a31588b803b30896522af4e90bc1e91bf23db6d0b9db80e427"
EXPECTED_STAGE2_FOCAL_PARTITION_SHA = "52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002"
EXPECTED_STAGE3_MANIFEST_SHA = "4ed41d1558dcc055efbb6522efe11b273062a187a7335a2975597fdbe0216fd7"
EXPECTED_STAGE3_BACKGROUND_PARTITION_SHA = "fa4cc04afcc40deb5e6b3d5ca76dc55b3720dfe070f7d90c3939e253b8006286"
EXPECTED_STAGE1_TOPOLOGIES_SHA = "030f6fb41d708fab2856b01b9f460f27e66a5171f7ba85084a6a7bf5e3a97ec0"
EXPECTED_STAGE1_WEIGHTS_SHA = "b9686750586419ae48e445c9c021991fa25eee6dc4957e55f4a12f0128ce9e58"
EXPECTED_STAGE1_WINDOW_INDEX_SHA = "96ed3fcca3d351a96b2434f17aa7510349d94138d11d60ccb56d45d1b57e237d"
EXPECTED_STAGE1_WEIGHT_SUMMARY_SHA = "7df2198149a456db8b80cb4d1d3c5d71c1aa3743821c1cdf408bf5c4432c7273"
EXPECTED_REGION_COUNTS = {"chr1": 117, "chr16A": 42, "chr16B": 2, "chr16_supergene": 52}
EXPECTED_GROUPS = {
    "geminata",
    "saevissima",
    "pusillignis",
    "invicta/macdonaghi_Sb",
    "invicta/macdonaghi_SB",
    "richteri_Sb",
    "richteri_SB",
}
ROLE_LABELS = {
    "invicta/macdonaghi_SB": "A",
    "invicta/macdonaghi_Sb": "B",
    "richteri_SB": "C",
    "richteri_Sb": "D",
}
SPLIT_TO_CLASS = {
    "AB|CD": "species",
    "AC|BD": "haplotype",
    "AD|BC": "third",
}
CLASS_TO_SPLIT = {value: key for key, value in SPLIT_TO_CLASS.items()}
CLASS_FIELDS = {"species": "W_species", "haplotype": "W_haplotype", "third": "W_third"}
Q_TOLERANCE = 1e-12
TIE_TOLERANCE = 1e-12
TSV_LINETERMINATOR = "\n"
DECIMAL_PLACES = Decimal("0.000000000000000")
OUTPUTS = [
    TOPOLOGY_CLASSES,
    WINDOW_SUPPORT,
    REGION_SUMMARY,
    PRIMARY_CONTRAST,
    REPORT,
    FIGURE_CAPTION,
    CHR16_SUPPORT_PDF,
    CHR16_SUPPORT_PNG,
    CHR1_SUPPORT_PDF,
    CHR1_SUPPORT_PNG,
    CHR16_D_PDF,
    CHR16_D_PNG,
    MANIFEST,
]


class Stage4AError(RuntimeError):
    """Raised when Stage-4A validation fails."""


class ParsedTree:
    def __init__(self, root: object, tips: tuple[str, ...]) -> None:
        self.root = root
        self.tips = tips


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def fmt_decimal(value: Decimal) -> str:
    return str(value.quantize(DECIMAL_PLACES))


def fmt_float(value: float) -> str:
    return f"{value:.12f}"


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


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.i = 0

    def parse(self) -> object:
        if not self.text:
            raise Stage4AError("empty Newick tree")
        node = self.parse_subtree()
        self.skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ";":
            raise Stage4AError("Newick tree missing terminal semicolon")
        self.i += 1
        self.skip_ws()
        if self.i != len(self.text):
            raise Stage4AError("unexpected content after Newick semicolon")
        return node

    def skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def parse_subtree(self) -> object:
        self.skip_ws()
        if self.i >= len(self.text):
            raise Stage4AError("unexpected end of Newick tree")
        if self.text[self.i] == "(":
            self.i += 1
            children = [self.parse_subtree()]
            while True:
                self.skip_ws()
                if self.i >= len(self.text):
                    raise Stage4AError("unterminated internal node")
                if self.text[self.i] == ",":
                    self.i += 1
                    children.append(self.parse_subtree())
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise Stage4AError(f"unexpected Newick character {self.text[self.i]!r}")
            self.consume_optional_label()
            self.consume_optional_branch_length()
            if len(children) < 2:
                raise Stage4AError("internal node has fewer than two children")
            return tuple(children)
        label = self.consume_tip_label()
        if not label:
            raise Stage4AError("empty tip label")
        self.consume_optional_branch_length()
        return label

    def consume_tip_label(self) -> str:
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in "(),:;":
            self.i += 1
        return self.text[start:self.i].strip()

    def consume_optional_label(self) -> None:
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
                raise Stage4AError("empty branch length")


def collect_tips(node: object, tips: list[str]) -> None:
    if isinstance(node, str):
        tips.append(node)
        return
    for child in node:
        collect_tips(child, tips)


def parse_newick(text: str) -> ParsedTree:
    root = NewickParser(text).parse()
    tips: list[str] = []
    collect_tips(root, tips)
    if len(tips) != len(set(tips)):
        duplicates = sorted(label for label, count in Counter(tips).items() if count > 1)
        raise Stage4AError(f"duplicate tip labels: {', '.join(duplicates)}")
    return ParsedTree(root=root, tips=tuple(tips))


def graph_and_leaves(node: object) -> tuple[dict[int, set[int]], dict[int, str]]:
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


def focal_distance_sum(graph: dict[int, set[int]], leaf_ids: dict[str, int], left: tuple[str, str], right: tuple[str, str]) -> int:
    return graph_distance(graph, leaf_ids[left[0]], leaf_ids[left[1]]) + graph_distance(graph, leaf_ids[right[0]], leaf_ids[right[1]])


def graph_distance(graph: dict[int, set[int]], start: int, end: int) -> int:
    queue = deque([(start, 0)])
    seen = {start}
    while queue:
        node, distance = queue.popleft()
        if node == end:
            return distance
        for neighbor in graph[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, distance + 1))
    raise Stage4AError("tree graph is disconnected")


def classify_focal_quartet(newick: str) -> tuple[str, str]:
    parsed = parse_newick(newick)
    if set(parsed.tips) != EXPECTED_GROUPS and not set(ROLE_LABELS).issubset(set(parsed.tips)):
        raise Stage4AError("tree lacks required focal groups")
    graph, leaves_by_id = graph_and_leaves(parsed.root)
    id_by_label = {label: node_id for node_id, label in leaves_by_id.items()}
    missing = [label for label in ROLE_LABELS if label not in id_by_label]
    if missing:
        raise Stage4AError(f"tree missing focal groups: {', '.join(missing)}")
    leaf_ids = {ROLE_LABELS[label]: id_by_label[label] for label in ROLE_LABELS}
    sums = {
        "AB|CD": focal_distance_sum(graph, leaf_ids, ("A", "B"), ("C", "D")),
        "AC|BD": focal_distance_sum(graph, leaf_ids, ("A", "C"), ("B", "D")),
        "AD|BC": focal_distance_sum(graph, leaf_ids, ("A", "D"), ("B", "C")),
    }
    min_value = min(sums.values())
    winners = [split for split, value in sums.items() if value == min_value]
    if len(winners) != 1:
        raise Stage4AError(f"ambiguous focal quartet classification: {sums}")
    split = winners[0]
    return split, SPLIT_TO_CLASS[split]


def verify_prior_stages() -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    if sha256(STAGE0_FOCAL_QUARTET) != EXPECTED_STAGE0_FOCAL_SHA:
        raise Stage4AError("Stage-0 focal quartet checksum changed")
    if sha256(STAGE1_MANIFEST) != EXPECTED_STAGE1_MANIFEST_SHA:
        raise Stage4AError("Stage-1 manifest checksum changed")
    if sha256(STAGE2_MANIFEST) != EXPECTED_STAGE2_MANIFEST_SHA:
        raise Stage4AError("Stage-2 manifest checksum changed")
    if sha256(STAGE2_FOCAL_PARTITION) != EXPECTED_STAGE2_FOCAL_PARTITION_SHA:
        raise Stage4AError("Stage-2 focal partition checksum changed")
    if sha256(STAGE3_MANIFEST) != EXPECTED_STAGE3_MANIFEST_SHA:
        raise Stage4AError("Stage-3 manifest checksum changed")
    if sha256(STAGE3_BACKGROUND_PARTITION) != EXPECTED_STAGE3_BACKGROUND_PARTITION_SHA:
        raise Stage4AError("Stage-3 background partition checksum changed")
    if sha256(STAGE1_TOPOLOGIES) != EXPECTED_STAGE1_TOPOLOGIES_SHA:
        raise Stage4AError("Stage-1 topology table checksum changed")
    if sha256(STAGE1_WEIGHTS_SPARSE) != EXPECTED_STAGE1_WEIGHTS_SHA:
        raise Stage4AError("Stage-1 sparse weights checksum changed")
    if sha256(STAGE1_WINDOW_INDEX) != EXPECTED_STAGE1_WINDOW_INDEX_SHA:
        raise Stage4AError("Stage-1 window index checksum changed")
    if sha256(STAGE1_WEIGHT_SUMMARY) != EXPECTED_STAGE1_WEIGHT_SUMMARY_SHA:
        raise Stage4AError("Stage-1 weight summary checksum changed")

    stage1 = json.loads(STAGE1_MANIFEST.read_text())
    expected_metrics = {
        "local_tree_rows": "213",
        "min_tree_tips": "267",
        "max_tree_tips": "267",
        "n_distinct_tip_sets": "1",
        "twisst_topologies": "945",
        "weight_rows": "213",
        "weight_columns": "945",
        "topology_crosscheck_mismatches": "0",
    }
    for key, expected in expected_metrics.items():
        if str(stage1.get("observed_metrics", {}).get(key)) != expected:
            raise Stage4AError(f"Stage-1 metric {key} expected {expected}, observed {stage1.get('observed_metrics', {}).get(key)}")
    for rel_path, expected_sha in stage1.get("raw_topology_file_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage4AError(f"Stage-1 raw checksum changed: {rel_path}")
    for rel_path, expected_sha in stage1.get("processed_output_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage4AError(f"Stage-1 processed checksum changed: {rel_path}")

    stage2 = json.loads(STAGE2_MANIFEST.read_text())
    for rel_path, expected_sha in stage2.get("stage2_output_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage4AError(f"Stage-2 output checksum changed: {rel_path}")

    stage3 = json.loads(STAGE3_MANIFEST.read_text())
    if stage3.get("stage3_background_partition_matches_stage0_species_split") is not True:
        raise Stage4AError("Stage-3 background partition no longer matches Stage-0 species split")
    for rel_path in stage3.get("outputs", []):
        if rel_path == "empirical/fire_ants_chr16/results/stage3_manifest.json":
            continue
        recorded = None
        if rel_path == "data/fire_ants_chr16/processed/stage3_background_partition.tsv":
            recorded = stage3["stage3_background_partition_checksum"]
        elif rel_path == "data/fire_ants_chr16/metadata/background_tree_provenance.tsv":
            recorded = stage3["background_tree_provenance_checksum"]
        elif rel_path == "data/fire_ants_chr16/processed/stage3_background_tree.nwk":
            recorded = stage3["normalized_background_tree_checksum"]
        if recorded and sha256(REPO_ROOT / rel_path) != recorded:
            raise Stage4AError(f"Stage-3 output checksum changed: {rel_path}")
    return stage1, stage2, stage3


def load_topology_classes() -> tuple[list[dict[str, object]], dict[str, str]]:
    fields, rows = read_tsv(STAGE1_TOPOLOGIES)
    required = {"topology_id", "topology_number", "tree_newick", "tree_sha256"}
    if not required.issubset(fields):
        raise Stage4AError("stage1_twisst_topologies.tsv missing required fields")
    if len(rows) != 945:
        raise Stage4AError(f"expected 945 grouped topologies, observed {len(rows)}")
    out_rows: list[dict[str, object]] = []
    class_by_topology: dict[str, str] = {}
    for row in rows:
        parsed = parse_newick(row["tree_newick"])
        if set(parsed.tips) != EXPECTED_GROUPS:
            raise Stage4AError(f"{row['topology_id']} does not contain the expected seven groups")
        split, klass = classify_focal_quartet(row["tree_newick"])
        out_rows.append(
            {
                "topology_id": row["topology_id"],
                "topology_number": row["topology_number"],
                "focal_quartet_split": split,
                "quartet_class": klass,
                "tree_sha256": row["tree_sha256"],
                "classification_method": "unrooted_graph_distances_four_point",
            }
        )
        class_by_topology[row["topology_id"]] = klass
    counts = Counter(row["quartet_class"] for row in out_rows)
    if counts != Counter({"species": 315, "haplotype": 315, "third": 315}):
        raise Stage4AError(f"unexpected topology class counts: {dict(counts)}")
    return out_rows, class_by_topology


def load_weight_summaries() -> dict[int, dict[str, int]]:
    fields, rows = read_tsv(STAGE1_WEIGHT_SUMMARY)
    required = {"window_index", "total_raw_weight", "n_nonzero_topologies"}
    if not required.issubset(fields):
        raise Stage4AError("stage1_twisst_weight_summary.tsv missing required fields")
    summaries = {
        int(row["window_index"]): {
            "total_raw_weight": int(row["total_raw_weight"]),
            "n_nonzero_topologies": int(row["n_nonzero_topologies"]),
        }
        for row in rows
    }
    if len(summaries) != 213:
        raise Stage4AError(f"expected 213 weight summary rows, observed {len(summaries)}")
    totals = {row["total_raw_weight"] for row in summaries.values()}
    if totals != {155843072}:
        raise Stage4AError(f"unexpected total raw weight values: {sorted(totals)}")
    return summaries


def aggregate_window_support(class_by_topology: dict[str, str]) -> list[dict[str, object]]:
    _, windows = read_tsv(STAGE1_WINDOW_INDEX)
    _, weights = read_tsv(STAGE1_WEIGHTS_SPARSE)
    summaries = load_weight_summaries()
    if len(windows) != 213:
        raise Stage4AError(f"expected 213 window index rows, observed {len(windows)}")
    region_counts = Counter(row["region"] for row in windows)
    if region_counts != Counter(EXPECTED_REGION_COUNTS):
        raise Stage4AError(f"unexpected region counts: {dict(region_counts)}")
    sums: dict[int, dict[str, int]] = {
        int(row["window_index"]): {"species": 0, "haplotype": 0, "third": 0} for row in windows
    }
    for row in weights:
        idx = int(row["window_index"])
        topology_id = row["topology_id"]
        weight = int(row["raw_weight"])
        sums[idx][class_by_topology[topology_id]] += weight
    out_rows: list[dict[str, object]] = []
    for row in windows:
        idx = int(row["window_index"])
        w_species = sums[idx]["species"]
        w_haplotype = sums[idx]["haplotype"]
        w_third = sums[idx]["third"]
        w_total = summaries[idx]["total_raw_weight"]
        if w_species + w_haplotype + w_third != w_total:
            raise Stage4AError(f"window {idx} class weights do not sum to total")
        q_s = Decimal(w_species) / Decimal(w_total)
        q_h = Decimal(w_haplotype) / Decimal(w_total)
        q_3 = Decimal(w_third) / Decimal(w_total)
        q_sum = q_s + q_h + q_3
        if abs(float(q_sum - Decimal(1))) > Q_TOLERANCE:
            raise Stage4AError(f"window {idx} q values do not sum to one")
        d_value = q_h - q_s
        supports = {"species": q_s, "haplotype": q_h, "third": q_3}
        sorted_supports = sorted(supports.items(), key=lambda item: item[1], reverse=True)
        if sorted_supports[0][1] - sorted_supports[1][1] <= Decimal(str(TIE_TOLERANCE)):
            dominant_class = "tie"
        else:
            dominant_class = sorted_supports[0][0]
        dominant_support = sorted_supports[0][1]
        second_support = sorted_supports[1][1]
        margin = dominant_support - second_support
        out_rows.append(
            {
                "window_index": idx,
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "mid": row["mid"],
                "region": row["region"],
                "upstream_window_id": row["upstream_window_id"],
                "gene1": row["gene1"],
                "gene2": row["gene2"],
                "gene3": row["gene3"],
                "gene4": row["gene4"],
                "W_species": w_species,
                "W_haplotype": w_haplotype,
                "W_third": w_third,
                "W_total": w_total,
                "q_S": fmt_decimal(q_s),
                "q_H": fmt_decimal(q_h),
                "q_3": fmt_decimal(q_3),
                "D": fmt_decimal(d_value),
                "dominant_class": dominant_class,
                "dominant_support": fmt_decimal(dominant_support),
                "second_support": fmt_decimal(second_support),
                "support_margin": fmt_decimal(margin),
            }
        )
    return out_rows


def decimal_values(rows: list[dict[str, object]], field: str) -> list[Decimal]:
    return [Decimal(str(row[field])) for row in rows]


def mean_decimal(values: list[Decimal]) -> Decimal:
    return sum(values, Decimal(0)) / Decimal(len(values))


def median_decimal(values: list[Decimal]) -> Decimal:
    ordered = sorted(values)
    n = len(ordered)
    if n % 2 == 1:
        return ordered[n // 2]
    return (ordered[n // 2 - 1] + ordered[n // 2]) / Decimal(2)


def region_summary_rows(window_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    regions: list[tuple[str, list[dict[str, object]]]] = []
    for region in ["chr1", "chr16A", "chr16B", "chr16_supergene"]:
        regions.append((region, [row for row in window_rows if row["region"] == region]))
    regions.insert(3, ("chr16_outside", [row for row in window_rows if row["region"] in {"chr16A", "chr16B"}]))
    rows: list[dict[str, object]] = []
    for region, subset in regions:
        if not subset:
            raise Stage4AError(f"region summary has no rows for {region}")
        n = len(subset)
        dominant_counts = Counter(row["dominant_class"] for row in subset)
        summary = {
            "region": region,
            "n_windows": n,
            "mean_q_S": fmt_decimal(mean_decimal(decimal_values(subset, "q_S"))),
            "median_q_S": fmt_decimal(median_decimal(decimal_values(subset, "q_S"))),
            "mean_q_H": fmt_decimal(mean_decimal(decimal_values(subset, "q_H"))),
            "median_q_H": fmt_decimal(median_decimal(decimal_values(subset, "q_H"))),
            "mean_q_3": fmt_decimal(mean_decimal(decimal_values(subset, "q_3"))),
            "median_q_3": fmt_decimal(median_decimal(decimal_values(subset, "q_3"))),
            "mean_D": fmt_decimal(mean_decimal(decimal_values(subset, "D"))),
            "median_D": fmt_decimal(median_decimal(decimal_values(subset, "D"))),
            "n_species_dominant": dominant_counts["species"],
            "n_haplotype_dominant": dominant_counts["haplotype"],
            "n_third_dominant": dominant_counts["third"],
            "n_ties": dominant_counts["tie"],
            "fraction_species_dominant": fmt_decimal(Decimal(dominant_counts["species"]) / Decimal(n)),
            "fraction_haplotype_dominant": fmt_decimal(Decimal(dominant_counts["haplotype"]) / Decimal(n)),
            "fraction_third_dominant": fmt_decimal(Decimal(dominant_counts["third"]) / Decimal(n)),
        }
        rows.append(summary)
    outside = next(row for row in rows if row["region"] == "chr16_outside")
    if int(outside["n_windows"]) != 44:
        raise Stage4AError("chr16_outside must contain exactly 44 pooled windows")
    return rows


def contrast_rows(window_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    def subset(region_set: set[str]) -> list[dict[str, object]]:
        return [row for row in window_rows if row["region"] in region_set]

    def contrast_row(name: str, inside: list[dict[str, object]], outside: list[dict[str, object]]) -> dict[str, object]:
        means_inside = {field: mean_decimal(decimal_values(inside, field)) for field in ["D", "q_S", "q_H", "q_3"]}
        means_outside = {field: mean_decimal(decimal_values(outside, field)) for field in ["D", "q_S", "q_H", "q_3"]}
        return {
            "comparison": name,
            "n_inside": len(inside),
            "n_outside": len(outside),
            "mean_D_inside": fmt_decimal(means_inside["D"]),
            "mean_D_outside": fmt_decimal(means_outside["D"]),
            "delta_D": fmt_decimal(means_inside["D"] - means_outside["D"]),
            "mean_q_S_inside": fmt_decimal(means_inside["q_S"]),
            "mean_q_S_outside": fmt_decimal(means_outside["q_S"]),
            "delta_q_S": fmt_decimal(means_inside["q_S"] - means_outside["q_S"]),
            "mean_q_H_inside": fmt_decimal(means_inside["q_H"]),
            "mean_q_H_outside": fmt_decimal(means_outside["q_H"]),
            "delta_q_H": fmt_decimal(means_inside["q_H"] - means_outside["q_H"]),
            "mean_q_3_inside": fmt_decimal(means_inside["q_3"]),
            "mean_q_3_outside": fmt_decimal(means_outside["q_3"]),
            "delta_q_3": fmt_decimal(means_inside["q_3"] - means_outside["q_3"]),
        }

    supergene = subset({"chr16_supergene"})
    chr16_outside = subset({"chr16A", "chr16B"})
    chr1 = subset({"chr1"})
    if len(supergene) != 52 or len(chr16_outside) != 44 or len(chr1) != 117:
        raise Stage4AError("unexpected contrast group sizes")
    return [
        contrast_row("chr16_supergene_vs_chr16_outside", supergene, chr16_outside),
        contrast_row("chr16_supergene_vs_chr1", supergene, chr1),
    ]


def region_span(region_id: str) -> tuple[float, float]:
    _, rows = read_tsv(REGION_MANIFEST)
    for row in rows:
        if row["region_id"] == region_id:
            return float(row["coordinate_start"]) / 1_000_000, float(row["coordinate_end"]) / 1_000_000
    raise Stage4AError(f"region not found in manifest: {region_id}")


def rows_for_chrom(window_rows: list[dict[str, object]], chrom: str) -> list[dict[str, object]]:
    return [row for row in window_rows if row["chrom"] == chrom]


def plot_support_tracks(window_rows: list[dict[str, object]], chrom: str, pdf_path: Path, png_path: Path, shade_supergene: bool) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    colors = {"q_S": "#1b9e77", "q_H": "#d95f02", "q_3": "#7570b3"}
    labels = {"q_S": "q_S species", "q_H": "q_H haplotype", "q_3": "q_3 third"}
    region_order = ["chr1"] if chrom == "chr1" else ["chr16A", "chr16B", "chr16_supergene"]
    for field in ["q_S", "q_H", "q_3"]:
        first = True
        for region in region_order:
            subset = [row for row in window_rows if row["region"] == region]
            if not subset:
                continue
            subset = sorted(subset, key=lambda row: float(row["mid"]))
            x = [float(row["mid"]) / 1_000_000 for row in subset]
            y = [float(row[field]) for row in subset]
            ax.plot(
                x,
                y,
                marker="o",
                markersize=3,
                linewidth=0.8,
                alpha=0.82,
                color=colors[field],
                label=labels[field] if first else None,
            )
            first = False
    if shade_supergene:
        start, end = region_span("chr16_supergene")
        ax.axvspan(start, end, color="#d9d9d9", alpha=0.35, label="author-designated supergene analysis region / observed BUSCO-window span")
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"{chrom} physical midpoint (Mb)")
    ax.set_ylabel("quartet support")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.set_title(f"Fire ants {chrom} TWISST-derived focal quartet support")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(pdf_path, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def plot_d_track(window_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.3))
    for region in ["chr16A", "chr16B", "chr16_supergene"]:
        subset = [row for row in window_rows if row["region"] == region]
        subset = sorted(subset, key=lambda row: float(row["mid"]))
        x = [float(row["mid"]) / 1_000_000 for row in subset]
        y = [float(row["D"]) for row in subset]
        ax.plot(x, y, marker="o", markersize=3, linewidth=0.8, color="#4c78a8")
    start, end = region_span("chr16_supergene")
    ax.axvspan(start, end, color="#d9d9d9", alpha=0.35, label="author-designated supergene analysis region / observed BUSCO-window span")
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("chr16 physical midpoint (Mb)")
    ax.set_ylabel("D = q_H - q_S")
    ax.set_title("Fire ants chr16 focal quartet contrast")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    fig.tight_layout()
    with PdfPages(CHR16_D_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(CHR16_D_PNG, dpi=200)
    plt.close(fig)


def make_figures(window_rows: list[dict[str, object]]) -> None:
    plot_support_tracks(window_rows, "chr16", CHR16_SUPPORT_PDF, CHR16_SUPPORT_PNG, shade_supergene=True)
    plot_support_tracks(window_rows, "chr1", CHR1_SUPPORT_PDF, CHR1_SUPPORT_PNG, shade_supergene=False)
    plot_d_track(window_rows)


def write_caption() -> None:
    FIGURE_CAPTION.write_text(
        "# Stage 4A chromosome-16 quartet-support figure caption\n\n"
        "Local quartet support across fire-ant chromosome 16, derived by summing the authors' published TWISST topology weights according to the independently frozen focal quartet. "
        "`q_S` is support for the frozen species-history quartet split `AB|CD`, where `A=invicta/macdonaghi_SB`, `B=invicta/macdonaghi_Sb`, `C=richteri_SB`, and `D=richteri_Sb`. "
        "`q_H` is support for the frozen SB/Sb haplotype-partition split `AC|BD`, and `q_3` is support for the third quartet resolution `AD|BC`. "
        "Each point is one published four-BUSCO-gene window; points are connected only within the independently defined source regions and no smoothing, fitted boundary, or inferential test is applied. "
        "The shaded interval marks the author-designated chromosome-16 supergene analysis region / observed BUSCO-window span, not exact inversion breakpoints.\n"
    )


def write_report(region_rows: list[dict[str, object]], contrast: list[dict[str, object]], class_counts: Counter) -> None:
    region_table = "\n".join(
        f"- `{row['region']}`: n={row['n_windows']}, mean q_S={row['mean_q_S']}, mean q_H={row['mean_q_H']}, mean q_3={row['mean_q_3']}, mean D={row['mean_D']}"
        for row in region_rows
    )
    dominant_table = "\n".join(
        f"- `{row['region']}`: species={row['n_species_dominant']}, haplotype={row['n_haplotype_dominant']}, third={row['n_third_dominant']}, ties={row['n_ties']}"
        for row in region_rows
    )
    primary = next(row for row in contrast if row["comparison"] == "chr16_supergene_vs_chr16_outside")
    secondary = next(row for row in contrast if row["comparison"] == "chr16_supergene_vs_chr1")
    REPORT.write_text(
        "# Fire ants chromosome 16 Stage 4A report\n\n"
        "## Purpose\n\n"
        "Stage 4A is the first formal unblinded local-topology stage. It classifies the 945 published TWISST group topologies according to the frozen focal quartet and aggregates published TWISST weights into per-window `q_S`, `q_H`, `q_3`, and `D=q_H-q_S` values.\n\n"
        "Stage 4A is descriptive only. No permutation test, circular shift, block-placement test, bootstrap hypothesis test, P-value, smoothing, fitted boundary, or ASTRAL/ASTER run was performed.\n\n"
        "## Frozen inputs\n\n"
        f"- Stage 0 focal quartet checksum: `{EXPECTED_STAGE0_FOCAL_SHA}`\n"
        f"- Stage 1 manifest checksum: `{EXPECTED_STAGE1_MANIFEST_SHA}`\n"
        f"- Stage 2 focal partition checksum: `{EXPECTED_STAGE2_FOCAL_PARTITION_SHA}`\n"
        f"- Stage 3 background partition checksum: `{EXPECTED_STAGE3_BACKGROUND_PARTITION_SHA}`\n\n"
        "## Topology classification\n\n"
        "Each seven-group TWISST topology was treated as unrooted. The three non-focal groups were ignored for the focal quartet, and the induced split among A/B/C/D was classified by unweighted graph-distance four-point sums.\n\n"
        f"Class counts: species={class_counts['species']}, haplotype={class_counts['haplotype']}, third={class_counts['third']}.\n\n"
        "## Window-level support calculation\n\n"
        "For each of 213 windows, raw TWISST weights were summed within the three topology classes to give `W_species`, `W_haplotype`, and `W_third`. Normalized supports were calculated as class weight divided by total class weight; every window satisfied `q_S + q_H + q_3 = 1` within `1e-12`.\n\n"
        "Dominant classes use a tie tolerance of `1e-12`.\n\n"
        "## Descriptive region summaries\n\n"
        f"{region_table}\n\n"
        "Dominant-class counts:\n\n"
        f"{dominant_table}\n\n"
        "## Primary descriptive contrast\n\n"
        f"Primary chr16 contrast `chr16_supergene_vs_chr16_outside`: `delta_D = {primary['delta_D']}` (`n_inside={primary['n_inside']}`, `n_outside={primary['n_outside']}`).\n\n"
        f"Secondary descriptive control `chr16_supergene_vs_chr1`: `delta_D = {secondary['delta_D']}` (`n_inside={secondary['n_inside']}`, `n_outside={secondary['n_outside']}`).\n\n"
        "## Figures\n\n"
        f"- `{CHR16_SUPPORT_PDF.relative_to(REPO_ROOT)}` and `{CHR16_SUPPORT_PNG.relative_to(REPO_ROOT)}`\n"
        f"- `{CHR1_SUPPORT_PDF.relative_to(REPO_ROOT)}` and `{CHR1_SUPPORT_PNG.relative_to(REPO_ROOT)}`\n"
        f"- `{CHR16_D_PDF.relative_to(REPO_ROOT)}` and `{CHR16_D_PNG.relative_to(REPO_ROOT)}`\n\n"
        "## Interpretation\n\n"
        "The formal Stage-4A question is whether local quartet support shifts from the independently frozen background species partition toward the independently frozen SB/Sb haplotype partition within the author-defined supergene region. This descriptive analysis may demonstrate an association between a recombination-suppressed structural haplotype and a coherent alternative genealogy regime, but it does not identify the historical mechanism responsible for sharing that haplotype across species.\n\n"
        "The source paper's introgression caveat remains in force: Stage 4A does not claim that introgression is absent, that MSRC without gene flow explains the fire-ant history, or that topology proves MSRC.\n\n"
        "## Limitations\n\n"
        "`chr16B` contains only two windows and is reported separately for transparency. The primary outside comparator pools the actual `chr16A + chr16B` windows. Formal spatial inference is reserved for Stage 4B.\n\n"
        "## Stage boundary\n\n"
        "Stop after Stage 4A. Stage 4B will perform the predeclared spatial-null inference.\n"
    )


def write_manifest(stage1: dict[str, object], stage2: dict[str, object], stage3: dict[str, object], class_counts: Counter, window_rows: list[dict[str, object]]) -> None:
    figure_checksums = {
        path.relative_to(REPO_ROOT).as_posix(): sha256(path)
        for path in [CHR16_SUPPORT_PDF, CHR16_SUPPORT_PNG, CHR1_SUPPORT_PDF, CHR1_SUPPORT_PNG, CHR16_D_PDF, CHR16_D_PNG]
    }
    data = {
        "stage": "Stage 4A",
        "description": "First formal unblinded local quartet-support analysis; descriptive only.",
        "stage0_focal_quartet_checksum": sha256(STAGE0_FOCAL_QUARTET),
        "stage1_manifest_checksum": sha256(STAGE1_MANIFEST),
        "stage2_manifest_checksum": sha256(STAGE2_MANIFEST),
        "stage2_focal_partition_checksum": sha256(STAGE2_FOCAL_PARTITION),
        "stage3_manifest_checksum": sha256(STAGE3_MANIFEST),
        "stage3_background_partition_checksum": sha256(STAGE3_BACKGROUND_PARTITION),
        "stage1_topology_table_checksum": sha256(STAGE1_TOPOLOGIES),
        "stage1_sparse_weight_checksum": sha256(STAGE1_WEIGHTS_SPARSE),
        "stage1_window_index_checksum": sha256(STAGE1_WINDOW_INDEX),
        "stage4a_output_checksums": {
            TOPOLOGY_CLASSES.relative_to(REPO_ROOT).as_posix(): sha256(TOPOLOGY_CLASSES),
            WINDOW_SUPPORT.relative_to(REPO_ROOT).as_posix(): sha256(WINDOW_SUPPORT),
            REGION_SUMMARY.relative_to(REPO_ROOT).as_posix(): sha256(REGION_SUMMARY),
            PRIMARY_CONTRAST.relative_to(REPO_ROOT).as_posix(): sha256(PRIMARY_CONTRAST),
            REPORT.relative_to(REPO_ROOT).as_posix(): sha256(REPORT),
            FIGURE_CAPTION.relative_to(REPO_ROOT).as_posix(): sha256(FIGURE_CAPTION),
        },
        "figure_checksums": figure_checksums,
        "processing_script": {
            "path": "empirical/fire_ants_chr16/scripts/04a_local_quartet_support.py",
            "sha256": sha256(Path(__file__)),
        },
        "class_counts": dict(class_counts),
        "n_windows": len(window_rows),
        "region_counts": dict(Counter(row["region"] for row in window_rows)),
        "q_sum_tolerance": Q_TOLERANCE,
        "tie_tolerance": TIE_TOLERANCE,
        "inferential_tests_performed": False,
        "p_values_calculated": False,
        "smoothing_applied_to_figures": False,
        "stage1_manifest_reference": stage1.get("observed_metrics", {}),
        "stage2_manifest_reference": {
            "focal_species_partition": stage2.get("focal_species_partition"),
            "focal_haplotype_partition": stage2.get("focal_haplotype_partition"),
        },
        "stage3_manifest_reference": {
            "stage3_background_partition_matches_stage0_species_split": stage3.get("stage3_background_partition_matches_stage0_species_split"),
        },
    }
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate_no_pvalue_columns() -> None:
    bad_patterns = ("p_value", "pvalue", "p_", "qvalue", "q_value")
    for path in [TOPOLOGY_CLASSES, WINDOW_SUPPORT, REGION_SUMMARY, PRIMARY_CONTRAST]:
        fields, _ = read_tsv(path)
        lowered = [field.lower() for field in fields]
        for field in lowered:
            if any(pattern == field or field.startswith(pattern) for pattern in bad_patterns):
                raise Stage4AError(f"{path.relative_to(REPO_ROOT)} contains inferential column {field}")


def run_stage4a() -> dict[str, object]:
    stage1, stage2, stage3 = verify_prior_stages()
    topology_class_rows, class_by_topology = load_topology_classes()
    class_counts = Counter(row["quartet_class"] for row in topology_class_rows)
    write_tsv(TOPOLOGY_CLASSES, topology_class_rows, ["topology_id", "topology_number", "focal_quartet_split", "quartet_class", "tree_sha256", "classification_method"])
    window_rows = aggregate_window_support(class_by_topology)
    write_tsv(
        WINDOW_SUPPORT,
        window_rows,
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
            "W_species",
            "W_haplotype",
            "W_third",
            "W_total",
            "q_S",
            "q_H",
            "q_3",
            "D",
            "dominant_class",
            "dominant_support",
            "second_support",
            "support_margin",
        ],
    )
    region_rows = region_summary_rows(window_rows)
    write_tsv(
        REGION_SUMMARY,
        region_rows,
        [
            "region",
            "n_windows",
            "mean_q_S",
            "median_q_S",
            "mean_q_H",
            "median_q_H",
            "mean_q_3",
            "median_q_3",
            "mean_D",
            "median_D",
            "n_species_dominant",
            "n_haplotype_dominant",
            "n_third_dominant",
            "n_ties",
            "fraction_species_dominant",
            "fraction_haplotype_dominant",
            "fraction_third_dominant",
        ],
    )
    contrast = contrast_rows(window_rows)
    write_tsv(
        PRIMARY_CONTRAST,
        contrast,
        [
            "comparison",
            "n_inside",
            "n_outside",
            "mean_D_inside",
            "mean_D_outside",
            "delta_D",
            "mean_q_S_inside",
            "mean_q_S_outside",
            "delta_q_S",
            "mean_q_H_inside",
            "mean_q_H_outside",
            "delta_q_H",
            "mean_q_3_inside",
            "mean_q_3_outside",
            "delta_q_3",
        ],
    )
    make_figures(window_rows)
    write_caption()
    validate_no_pvalue_columns()
    write_report(region_rows, contrast, class_counts)
    write_manifest(stage1, stage2, stage3, class_counts, window_rows)
    return {
        "class_counts": class_counts,
        "n_windows": len(window_rows),
        "region_counts": Counter(row["region"] for row in window_rows),
        "manifest_sha": sha256(MANIFEST),
    }


class Stage4ATests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.metrics = run_stage4a()

    def test_prior_stage_checksums_validate(self) -> None:
        verify_prior_stages()

    def test_exactly_945_grouped_topologies_loaded(self) -> None:
        _, rows = read_tsv(STAGE1_TOPOLOGIES)
        self.assertEqual(len(rows), 945)

    def test_every_topology_has_required_seven_groups(self) -> None:
        _, rows = read_tsv(STAGE1_TOPOLOGIES)
        for row in rows:
            self.assertEqual(set(parse_newick(row["tree_newick"]).tips), EXPECTED_GROUPS)

    def test_all_945_focal_quartets_classify(self) -> None:
        _, rows = read_tsv(TOPOLOGY_CLASSES)
        self.assertEqual(len(rows), 945)
        self.assertTrue(all(row["quartet_class"] in {"species", "haplotype", "third"} for row in rows))

    def test_no_ambiguous_topology_classification(self) -> None:
        _, rows = read_tsv(TOPOLOGY_CLASSES)
        self.assertEqual(sum(row["quartet_class"] == "ambiguous" for row in rows), 0)

    def test_class_counts_balanced(self) -> None:
        _, rows = read_tsv(TOPOLOGY_CLASSES)
        self.assertEqual(Counter(row["quartet_class"] for row in rows), Counter({"species": 315, "haplotype": 315, "third": 315}))

    def test_classification_invariant_to_rooting(self) -> None:
        cases = [
            "((A,B),(C,D));",
            "(C,D,(A,B));",
            "((C,D),(A,B));",
        ]
        mapped = [case.translate(str.maketrans({"A": "A", "B": "B", "C": "C", "D": "D"})) for case in cases]
        for case in mapped:
            self.assertEqual(classify_synthetic(case), "species")

    def test_classification_invariant_to_child_order(self) -> None:
        self.assertEqual(classify_synthetic("((B,A),(D,C));"), "species")
        self.assertEqual(classify_synthetic("((C,A),(D,B));"), "haplotype")

    def test_synthetic_species_quartet(self) -> None:
        self.assertEqual(classify_synthetic("((A,B),(C,D));"), "species")

    def test_synthetic_haplotype_quartet(self) -> None:
        self.assertEqual(classify_synthetic("((A,C),(B,D));"), "haplotype")

    def test_synthetic_third_quartet(self) -> None:
        self.assertEqual(classify_synthetic("((A,D),(B,C));"), "third")

    def test_synthetic_extra_taxa_branch_lengths_labels(self) -> None:
        self.assertEqual(classify_synthetic("((out1:0.1,(A:1,B:1)99:0.2),(out2:0.3,(C:1,D:1)88:0.4));"), "species")
        self.assertEqual(classify_synthetic("(out1,((A,C),(out2,(B,D))));"), "haplotype")

    def test_exactly_213_windows_processed(self) -> None:
        _, rows = read_tsv(WINDOW_SUPPORT)
        self.assertEqual(len(rows), 213)

    def test_region_counts(self) -> None:
        _, rows = read_tsv(WINDOW_SUPPORT)
        counts = Counter(row["region"] for row in rows)
        self.assertEqual(counts["chr1"], 117)
        self.assertEqual(counts["chr16A"], 42)
        self.assertEqual(counts["chr16B"], 2)
        self.assertEqual(counts["chr16_supergene"], 52)

    def test_weight_sums_equal_total(self) -> None:
        _, rows = read_tsv(WINDOW_SUPPORT)
        for row in rows:
            self.assertEqual(int(row["W_species"]) + int(row["W_haplotype"]) + int(row["W_third"]), int(row["W_total"]))

    def test_q_values_sum_to_one(self) -> None:
        _, rows = read_tsv(WINDOW_SUPPORT)
        for row in rows:
            total = float(row["q_S"]) + float(row["q_H"]) + float(row["q_3"])
            self.assertLessEqual(abs(total - 1.0), Q_TOLERANCE)

    def test_q_values_and_d_ranges(self) -> None:
        _, rows = read_tsv(WINDOW_SUPPORT)
        for row in rows:
            for field in ["q_S", "q_H", "q_3"]:
                self.assertGreaterEqual(float(row[field]), 0.0)
                self.assertLessEqual(float(row[field]), 1.0)
            self.assertGreaterEqual(float(row["D"]), -1.0)
            self.assertLessEqual(float(row["D"]), 1.0)

    def test_chr16_outside_contains_44_windows(self) -> None:
        _, rows = read_tsv(REGION_SUMMARY)
        outside = next(row for row in rows if row["region"] == "chr16_outside")
        self.assertEqual(int(outside["n_windows"]), 44)

    def test_no_inferential_pvalue_column(self) -> None:
        validate_no_pvalue_columns()

    def test_no_smoothing_applied_to_figures(self) -> None:
        manifest = json.loads(MANIFEST.read_text())
        self.assertFalse(manifest["smoothing_applied_to_figures"])

    def test_repeated_runs_byte_stable(self) -> None:
        before = {path: sha256(path) for path in OUTPUTS}
        run_stage4a()
        after = {path: sha256(path) for path in OUTPUTS}
        self.assertEqual(before, after)


def classify_synthetic(newick: str) -> str:
    label_map = {
        "A": "invicta/macdonaghi_SB",
        "B": "invicta/macdonaghi_Sb",
        "C": "richteri_SB",
        "D": "richteri_Sb",
    }
    text = re.sub(
        r"(?<=[(,])([^():,;]+)(?=[:),;])",
        lambda match: label_map.get(match.group(1), match.group(1)),
        newick,
    )
    if not text.strip().endswith(";"):
        text = text.strip() + ";"
    _, klass = classify_focal_quartet(text)
    return klass


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage4ATests)
    with tempfile.TemporaryDirectory(prefix="fire_ants_stage4a_tests_"):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run embedded Stage-4A tests after writing outputs")
    args = parser.parse_args(argv)
    try:
        metrics = run_stage4a()
    except Stage4AError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.run_tests:
        return 0 if run_tests() else 1
    print(f"Wrote {WINDOW_SUPPORT.relative_to(REPO_ROOT)}")
    print(f"Classified {sum(metrics['class_counts'].values())} topologies")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
