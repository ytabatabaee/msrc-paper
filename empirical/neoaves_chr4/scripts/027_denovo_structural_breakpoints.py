#!/usr/bin/env python3
"""Infer de novo chr4 structural breakpoints from complete maf2synteny blocks.

This stage is intentionally structural-only. During inference it reads only
``rearrangements/totab_chr4.csv.xz`` and never reads prior outlier windows,
Stage-2 structural intervals, or Stage-3 genealogy change points.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
import shutil
import statistics
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"
INPUT_PATH = REPO_ROOT / "rearrangements" / "totab_chr4.csv.xz"
REFERENCE_DESCRIPTION = "GCF_000002315.5_GalGal6.chr4"

ACCESSION_TO_SPECIES = {
    "GCA_017639555": "Stork",
    "GCA_009819775": "Flamingo",
    "GCA_901699155": "Dove",
    "GCA_009769525": "Sandgrouse",
    "GCA_009769465": "Turaco",
    "GCA_017976375": "Cuckoo",
}
FOCAL_SPECIES = ["Stork", "Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo"]
PLOT_ORDER = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]

FORBIDDEN_INFERENCE_INPUTS = {
    "rearrangements/outlier-regions-by-maf2synteny.tsv",
    "rearrangements/outlier-regions-by-maf2synteny-summary.tsv",
    "rearrangements/outlier_plot.R",
    "empirical/neoaves_chr4/results/structural_intervals.tsv",
    "empirical/neoaves_chr4/results/canonical_structural_regions.tsv",
    "empirical/neoaves_chr4/results/canonical_structural_boundaries.tsv",
    "empirical/neoaves_chr4/results/locus_table_chr4.tsv",
    "empirical/neoaves_chr4/results/genealogy_change_points.tsv",
    "empirical/neoaves_chr4/results/genealogy_change_points_raw.tsv",
    "empirical/neoaves_chr4/results/chr4_change_points_vs_structure.tsv",
    "empirical/neoaves_chr4/figures/chr4_change_points_vs_structure.pdf",
    "empirical/neoaves_chr4/figures/chr4_change_points_vs_structure.png",
}

BLOCKS_OUT = RESULTS_DIR / "denovo_chr4_blocks.tsv"
ADJ_OUT = RESULTS_DIR / "denovo_adjacent_block_features.tsv"
RUNS_OUT = RESULTS_DIR / "denovo_synteny_runs.tsv"
RAW_BP_OUT = RESULTS_DIR / "denovo_species_breakpoints_raw.tsv"
RAW_DISRUPTIONS_OUT = RESULTS_DIR / "denovo_raw_structural_disruptions.tsv"
SPECIES_BP_OUT = RESULTS_DIR / "denovo_species_breakpoints.tsv"
CONSENSUS_OUT = RESULTS_DIR / "denovo_consensus_structural_breakpoints.tsv"
BREAKPOINT_SETS_OUT = RESULTS_DIR / "denovo_breakpoint_sets.tsv"
SEGMENTS_OUT = RESULTS_DIR / "denovo_structural_segments.tsv"
SENSITIVITY_OUT = RESULTS_DIR / "denovo_breakpoint_sensitivity.tsv"
MANIFEST_OUT = RESULTS_DIR / "denovo_breakpoint_manifest.json"
STAGE27B_MANIFEST_OUT = RESULTS_DIR / "denovo_breakpoint_manifest_stage27b.json"
STAGE27B_CONSENSUS_OUT = RESULTS_DIR / "denovo_consensus_structural_breakpoints_stage27b.tsv"
STAGE27_AUDIT_OUT = RESULTS_DIR / "stage27b_vs_stage27c_consensus.tsv"
MAP_PDF = FIGURES_DIR / "denovo_chr4_synteny_map.pdf"
MAP_PNG = FIGURES_DIR / "denovo_chr4_synteny_map.png"
CONS_PDF = FIGURES_DIR / "denovo_consensus_breakpoints.pdf"
CONS_PNG = FIGURES_DIR / "denovo_consensus_breakpoints.png"
INVALID_MANIFEST_OUT = RESULTS_DIR / "denovo_breakpoint_manifest_invalid_orientation.json"

DEFAULT_SETTINGS = {
    "lenient": {"k_mad": 2.5, "min_score": 2.5},
    "default": {"k_mad": 3.5, "min_score": 3.5},
    "stringent": {"k_mad": 5.0, "min_score": 5.0},
}


@dataclass
class Block:
    species: str
    accession: str
    block_id: str
    reference_chrom: str
    reference_start: int
    reference_end: int
    reference_strand: str
    reference_seq_id: str
    reference_original_start: int
    reference_original_end: int
    query_chrom: str
    query_start: int
    query_end: int
    query_strand: str
    query_seq_id: str
    query_original_start: int
    query_original_end: int
    query_description: str
    relative_orientation: str
    block_length: int


@dataclass
class Run:
    species: str
    run_id: str
    reference_start: int
    reference_end: int
    query_chrom: str
    query_start: int
    query_end: int
    relative_orientation: str
    n_raw_blocks: int
    source_block_ids: list[str]


def relative_orientation(reference_strand: str, query_strand: str) -> str:
    return "+" if reference_strand == query_strand else "-"


def fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    if isinstance(value, (list, tuple)):
        return ";".join(map(str, value))
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: fmt(row.get(col, "")) for col in columns})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assert_inference_inputs_allowed(paths: list[Path]) -> None:
    repo = REPO_ROOT.resolve()
    forbidden = {str((repo / p).resolve()) for p in FORBIDDEN_INFERENCE_INPUTS}
    for path in paths:
        resolved = str(path.resolve())
        if resolved in forbidden:
            raise RuntimeError(f"Forbidden inference input requested: {path}")
    if len(paths) != 1 or paths[0].resolve() != INPUT_PATH.resolve():
        raise RuntimeError("De novo inference must use only rearrangements/totab_chr4.csv.xz")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def iqr(values: list[int]) -> tuple[float, float, float]:
    if not values:
        return 0.0, 0.0, 0.0
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[0], q[1], q[2]


def mad_threshold(values: list[float], k: float, floor: float) -> float:
    if not values:
        return floor
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    return max(floor, med + k * 1.4826 * mad)


def parse_complete_blocks(input_path: Path) -> tuple[list[Block], int]:
    assert_inference_inputs_allowed([input_path])
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with lzma.open(input_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Seq_id", "Strand", "Start", "End", "Length", "Block", "sp", "Size", "Description"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {input_path}: {sorted(missing)}")
        for row in reader:
            if row["sp"] in ACCESSION_TO_SPECIES:
                grouped[(row["sp"], row["Block"])].append(row)

    blocks: list[Block] = []
    for (accession, block_id), rows in grouped.items():
        ref = [r for r in rows if r["Description"] == REFERENCE_DESCRIPTION]
        query = [r for r in rows if r["Description"] != REFERENCE_DESCRIPTION]
        if len(ref) != 1 or len(query) != 1:
            continue
        rr, qr = ref[0], query[0]
        rs, re = sorted((int(rr["Start"]), int(rr["End"])))
        qs, qe = sorted((int(qr["Start"]), int(qr["End"])))
        blocks.append(
            Block(
                species=ACCESSION_TO_SPECIES[accession],
                accession=accession,
                block_id=block_id,
                reference_chrom="chr4",
                reference_start=rs,
                reference_end=re,
                reference_strand=rr["Strand"],
                reference_seq_id=rr["Seq_id"],
                reference_original_start=int(rr["Start"]),
                reference_original_end=int(rr["End"]),
                query_chrom=query_chrom_label(qr["Description"]),
                query_start=qs,
                query_end=qe,
                query_strand=qr["Strand"],
                query_seq_id=qr["Seq_id"],
                query_original_start=int(qr["Start"]),
                query_original_end=int(qr["End"]),
                query_description=qr["Description"],
                relative_orientation=relative_orientation(rr["Strand"], qr["Strand"]),
                block_length=int(qr["Length"]),
            )
        )
    blocks.sort(key=lambda b: (FOCAL_SPECIES.index(b.species), b.reference_start, b.reference_end, int(b.block_id)))
    return blocks, len(grouped)


def query_chrom_label(description: str) -> str:
    """Keep scaffold/chromosome accessions intact while dropping assembly text."""
    if ".pri." in description:
        return description.split(".pri.", 1)[1]
    if "bStrTur1.1." in description:
        return description.split("bStrTur1.1.", 1)[1]
    parts = description.split(".")
    if len(parts) >= 2 and parts[-2] in {"1", "2"} and not parts[-1].startswith("chr"):
        return ".".join(parts[-2:])
    return parts[-1]


def block_rows(blocks: list[Block]) -> list[dict[str, object]]:
    return [b.__dict__.copy() for b in blocks]


def expected_query_gap(left, right) -> int:
    orientation = left.relative_orientation
    if orientation == "-":
        return left.query_start - right.query_end
    return right.query_start - left.query_end


def adjacency_features(left, right) -> dict[str, object]:
    same_chrom = left.query_chrom == right.query_chrom
    strand_left = getattr(left, "query_strand", "")
    strand_right = getattr(right, "query_strand", "")
    ref_strand_left = getattr(left, "reference_strand", "")
    ref_strand_right = getattr(right, "reference_strand", "")
    rel_left = left.relative_orientation
    rel_right = right.relative_orientation
    ref_gap = right.reference_start - left.reference_end
    if same_chrom and rel_left == rel_right:
        qgap = expected_query_gap(left, right)
        consistent = qgap >= -10_000
        reversal = qgap < -10_000
        overlap = max(0, -qgap)
        jump = abs(qgap)
    elif same_chrom:
        qgap = min(abs(right.query_start - left.query_end), abs(left.query_start - right.query_end))
        consistent = False
        reversal = True
        overlap = 0
        jump = qgap
    else:
        qgap = ""
        consistent = False
        reversal = False
        overlap = 0
        jump = ""
    discrepancy = abs(jump - ref_gap) if isinstance(jump, int) else ""
    return {
        "species": left.species,
        "left_id": getattr(left, "block_id", getattr(left, "run_id", "")),
        "right_id": getattr(right, "block_id", getattr(right, "run_id", "")),
        "reference_left_end": left.reference_end,
        "reference_right_start": right.reference_start,
        "reference_gap_bp": ref_gap,
        "query_gap_bp": qgap,
        "same_query_chromosome": int(same_chrom),
        "reference_strand_left": ref_strand_left,
        "reference_strand_right": ref_strand_right,
        "strand_left": strand_left,
        "strand_right": strand_right,
        "query_strand_left": strand_left,
        "query_strand_right": strand_right,
        "relative_orientation_left": rel_left,
        "relative_orientation_right": rel_right,
        "strand_switch": int(rel_left != rel_right),
        "relative_orientation_switch": int(rel_left != rel_right),
        "query_order_consistent": int(consistent),
        "query_order_reversal": int(reversal),
        "query_overlap": overlap,
        "query_jump_bp": jump,
        "reference_query_gap_discrepancy": discrepancy,
        "scaffold_or_chromosome_switch": int(not same_chrom),
    }


def adjacent_block_features(blocks: list[Block]) -> list[dict[str, object]]:
    rows = []
    by_species: dict[str, list[Block]] = defaultdict(list)
    for block in blocks:
        by_species[block.species].append(block)
    for species in FOCAL_SPECIES:
        bs = by_species.get(species, [])
        for left, right in zip(bs, bs[1:]):
            rows.append(adjacency_features(left, right))
    return rows


def empirical_gap_report(adj_rows: list[dict[str, object]]) -> dict[str, object]:
    ref_gaps = [max(0, int(r["reference_gap_bp"])) for r in adj_rows]
    query_gaps = [int(r["query_jump_bp"]) for r in adj_rows if str(r["query_jump_bp"]) != ""]
    discrep = [int(r["reference_query_gap_discrepancy"]) for r in adj_rows if str(r["reference_query_gap_discrepancy"]) != ""]
    q1, med, q3 = iqr(ref_gaps)
    qq1, qmed, qq3 = iqr(query_gaps)
    d1, dmed, d3 = iqr(discrep)
    return {
        "reference_gap_q1": q1,
        "reference_gap_median": med,
        "reference_gap_q3": q3,
        "reference_gap_p95": percentile(ref_gaps, 0.95),
        "query_jump_q1": qq1,
        "query_jump_median": qmed,
        "query_jump_q3": qq3,
        "query_jump_p95": percentile(query_gaps, 0.95),
        "gap_discrepancy_q1": d1,
        "gap_discrepancy_median": dmed,
        "gap_discrepancy_q3": d3,
        "gap_discrepancy_p95": percentile(discrep, 0.95),
    }


def percentile(values: list[int], p: float) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    idx = min(len(values) - 1, max(0, int(round((len(values) - 1) * p))))
    return float(values[idx])


def choose_default_tolerances(adj_rows: list[dict[str, object]], max_reference_gap: int | None, max_query_gap: int | None) -> tuple[int, int]:
    ref_gaps = [max(0, int(r["reference_gap_bp"])) for r in adj_rows]
    query_jumps = [int(r["query_jump_bp"]) for r in adj_rows if str(r["query_jump_bp"]) != "" and int(r["query_jump_bp"]) >= 0]
    ref = max_reference_gap if max_reference_gap is not None else int(max(100_000, percentile(ref_gaps, 0.95)))
    query = max_query_gap if max_query_gap is not None else int(max(250_000, percentile(query_jumps, 0.95)))
    return ref, query


def can_merge(left: Block, right: Block, max_reference_gap: int, max_query_gap: int) -> bool:
    feat = adjacency_features(left, right)
    if not feat["same_query_chromosome"] or feat["strand_switch"]:
        return False
    if not feat["query_order_consistent"]:
        return False
    if int(feat["reference_gap_bp"]) > max_reference_gap:
        return False
    if int(feat["query_jump_bp"]) > max_query_gap:
        return False
    if int(feat["reference_query_gap_discrepancy"]) > max(max_reference_gap, max_query_gap):
        return False
    return True


def build_synteny_runs(blocks: list[Block], max_reference_gap: int, max_query_gap: int) -> list[Run]:
    runs: list[Run] = []
    by_species: dict[str, list[Block]] = defaultdict(list)
    for block in blocks:
        by_species[block.species].append(block)
    for species in FOCAL_SPECIES:
        current: list[Block] = []
        run_idx = 1
        for block in by_species.get(species, []):
            if not current or can_merge(current[-1], block, max_reference_gap, max_query_gap):
                current.append(block)
                continue
            runs.append(make_run(species, run_idx, current))
            run_idx += 1
            current = [block]
        if current:
            runs.append(make_run(species, run_idx, current))
    return runs


def make_run(species: str, idx: int, blocks: list[Block]) -> Run:
    return Run(
        species=species,
        run_id=f"{species[:3].upper()}_R{idx:04d}",
        reference_start=min(b.reference_start for b in blocks),
        reference_end=max(b.reference_end for b in blocks),
        query_chrom=blocks[0].query_chrom,
        query_start=min(b.query_start for b in blocks),
        query_end=max(b.query_end for b in blocks),
        relative_orientation=blocks[0].relative_orientation,
        n_raw_blocks=len(blocks),
        source_block_ids=[b.block_id for b in blocks],
    )


def run_rows(runs: list[Run]) -> list[dict[str, object]]:
    return [r.__dict__.copy() for r in runs]


def robust_scales(run_adj: list[dict[str, object]]) -> tuple[float, float]:
    jumps = [float(r["query_jump_bp"]) for r in run_adj if str(r["query_jump_bp"]) != ""]
    discrep = [float(r["reference_query_gap_discrepancy"]) for r in run_adj if str(r["reference_query_gap_discrepancy"]) != ""]
    return max(1.0, mad_threshold(jumps, 3.5, 100_000.0)), max(1.0, mad_threshold(discrep, 3.5, 100_000.0))


def run_length(run: Run) -> int:
    return max(0, run.reference_end - run.reference_start)


def same_regime(left: Run, right: Run) -> bool:
    return left.query_chrom == right.query_chrom and left.relative_orientation == right.relative_orientation


def flank_support(runs: list[Run], transition_idx: int) -> tuple[int, int, int, int]:
    left = runs[transition_idx]
    right = runs[transition_idx + 1]
    left_runs = left.n_raw_blocks
    left_bp = run_length(left)
    i = transition_idx - 1
    while i >= 0 and same_regime(runs[i], left):
        left_runs += runs[i].n_raw_blocks
        left_bp += run_length(runs[i])
        i -= 1
    right_runs = right.n_raw_blocks
    right_bp = run_length(right)
    i = transition_idx + 2
    while i < len(runs) and same_regime(right, runs[i]):
        right_runs += runs[i].n_raw_blocks
        right_bp += run_length(runs[i])
        i += 1
    return left_bp, right_bp, left_runs, right_runs


def score_run_transition(row: dict[str, object], jump_scale: float, discrepancy_scale: float) -> dict[str, object]:
    chrom = 5.0 * int(row["scaffold_or_chromosome_switch"])
    orient = 3.0 * int(row["strand_switch"])
    reversal = 3.0 * int(row["query_order_reversal"])
    jump_raw = float(row["query_jump_bp"]) if str(row["query_jump_bp"]) != "" else 0.0
    discrep_raw = float(row["reference_query_gap_discrepancy"]) if str(row["reference_query_gap_discrepancy"]) != "" else 0.0
    jump = min(3.0, jump_raw / jump_scale)
    discrep = min(2.0, discrep_raw / discrepancy_scale)
    score = chrom + orient + reversal + jump + discrep
    evidence = []
    if chrom:
        evidence.append("chromosome_switch")
    if orient:
        evidence.append("orientation_switch")
    if reversal:
        evidence.append("order_reversal")
    if jump >= 1:
        evidence.append("large_query_jump")
    if discrep >= 1:
        evidence.append("gap_discrepancy")
    if not evidence:
        evidence.append("weak_gap_signal")
    out = row.copy()
    out.update(
        {
            "score_chromosome_switch": chrom,
            "score_orientation_switch": orient,
            "score_order_reversal": reversal,
            "score_query_jump": jump,
            "score_gap_discrepancy": discrep,
            "breakpoint_score": score,
            "evidence": ",".join(evidence),
        }
    )
    return out


def build_raw_disruptions(runs: list[Run], min_flank_bp: int, min_flank_blocks: int) -> list[dict[str, object]]:
    run_adj = []
    by_species: dict[str, list[Run]] = defaultdict(list)
    for run in runs:
        by_species[run.species].append(run)
    for species in FOCAL_SPECIES:
        rs = by_species.get(species, [])
        for idx, (left, right) in enumerate(zip(rs, rs[1:])):
            row = adjacency_features(left, right)
            row["left_run_id"] = left.run_id
            row["right_run_id"] = right.run_id
            row["left_run_index"] = idx + 1
            row["right_run_index"] = idx + 2
            row["left_query_chrom"] = left.query_chrom
            row["right_query_chrom"] = right.query_chrom
            row["left_relative_orientation"] = left.relative_orientation
            row["right_relative_orientation"] = right.relative_orientation
            row["reference_position"] = int(round((left.reference_end + right.reference_start) / 2))
            left_bp, right_bp, left_blocks, right_blocks = flank_support(rs, idx)
            jump_bp = int(row["query_jump_bp"]) if str(row["query_jump_bp"]) != "" else 0
            gap_discrepancy = int(row["reference_query_gap_discrepancy"]) if str(row["reference_query_gap_discrepancy"]) != "" else 0
            regime_change_signal = (
                not same_regime(left, right)
                or jump_bp >= min_flank_bp
                or gap_discrepancy >= min_flank_bp
                or int(row["query_order_reversal"])
            )
            row["left_support_bp"] = left_bp
            row["right_support_bp"] = right_bp
            row["left_support_blocks"] = left_blocks
            row["right_support_blocks"] = right_blocks
            row["persistent_regime_change"] = int(
                left_bp >= min_flank_bp
                and right_bp >= min_flank_bp
                and left_blocks >= min_flank_blocks
                and right_blocks >= min_flank_blocks
                and regime_change_signal
            )
            run_adj.append(row)
    jump_scale, discrepancy_scale = robust_scales(run_adj)
    return [score_run_transition(row, jump_scale, discrepancy_scale) for row in run_adj]


def call_breakpoints(disruptions: list[dict[str, object]], threshold_setting: str, k_mad: float, min_score: float) -> list[dict[str, object]]:
    threshold = mad_threshold([float(r["breakpoint_score"]) for r in disruptions], k_mad, min_score)
    calls = []
    counts = Counter()
    for row in disruptions:
        persistent = int(row["persistent_regime_change"])
        strong_categorical = int(row["scaffold_or_chromosome_switch"]) or int(row["strand_switch"]) or (
            int(row["query_order_reversal"]) and float(row["breakpoint_score"]) >= min_score
        )
        quantitative = float(row["breakpoint_score"]) >= threshold
        if persistent and (strong_categorical or quantitative):
            counts[row["species"]] += 1
            calls.append(
                {
                    "species": row["species"],
                    "breakpoint_id": f"{row['species']}_{threshold_setting}_raw_{counts[row['species']]:04d}",
                    "reference_position": row["reference_position"],
                    "left_run_id": row["left_run_id"],
                    "right_run_id": row["right_run_id"],
                    "left_run_index": row["left_run_index"],
                    "right_run_index": row["right_run_index"],
                    "left_query_chrom": row["left_query_chrom"],
                    "right_query_chrom": row["right_query_chrom"],
                    "left_relative_orientation": row["left_relative_orientation"],
                    "right_relative_orientation": row["right_relative_orientation"],
                    "chromosome_switch": row["scaffold_or_chromosome_switch"],
                    "orientation_switch": row["strand_switch"],
                    "order_reversal": row["query_order_reversal"],
                    "query_jump_bp": row["query_jump_bp"],
                    "gap_discrepancy": row["reference_query_gap_discrepancy"],
                    "persistent_regime_change": row["persistent_regime_change"],
                    "left_support_bp": row["left_support_bp"],
                    "right_support_bp": row["right_support_bp"],
                    "left_support_blocks": row["left_support_blocks"],
                    "right_support_blocks": row["right_support_blocks"],
                    "score_chromosome_switch": row["score_chromosome_switch"],
                    "score_orientation_switch": row["score_orientation_switch"],
                    "score_order_reversal": row["score_order_reversal"],
                    "score_query_jump": row["score_query_jump"],
                    "score_gap_discrepancy": row["score_gap_discrepancy"],
                    "breakpoint_score": row["breakpoint_score"],
                    "threshold": threshold,
                    "threshold_setting": threshold_setting,
                    "robustness": "",
                    "evidence": row["evidence"],
                }
            )
    return calls


def merge_species_calls(raw_calls: list[dict[str, object]], within_bp: int) -> list[dict[str, object]]:
    calls = [r for r in raw_calls if r["threshold_setting"] == "default"]
    all_threshold_pos = defaultdict(list)
    for r in raw_calls:
        all_threshold_pos[(r["species"], r["left_run_id"], r["right_run_id"])].append(r["threshold_setting"])
    merged = []
    for species in FOCAL_SPECIES:
        rows = sorted([r for r in calls if r["species"] == species], key=lambda r: int(r["reference_position"]))
        clusters: list[list[dict[str, object]]] = []
        for row in rows:
            close = bool(clusters) and int(row["reference_position"]) - int(clusters[-1][-1]["reference_position"]) <= within_bp
            adjacent_in_run_space = bool(clusters) and int(row["left_run_index"]) <= int(clusters[-1][-1]["right_run_index"]) + 1
            if not close or not adjacent_in_run_space:
                clusters.append([row])
            else:
                clusters[-1].append(row)
        for idx, cluster in enumerate(clusters, 1):
            weights = [max(0.001, float(r["breakpoint_score"])) for r in cluster]
            pos = int(round(sum(int(r["reference_position"]) * w for r, w in zip(cluster, weights)) / sum(weights)))
            best = max(cluster, key=lambda r: float(r["breakpoint_score"])).copy()
            key_settings = set()
            for r in cluster:
                key_settings.update(all_threshold_pos[(r["species"], r["left_run_id"], r["right_run_id"])])
            robustness = "high" if {"lenient", "default", "stringent"} <= key_settings else ("medium" if len(key_settings) >= 2 else "low")
            best.update(
                {
                    "breakpoint_id": f"{species}_BP{idx:03d}",
                    "reference_position": pos,
                    "threshold_setting": "default",
                    "robustness": robustness,
                    "contributing_raw_calls": ";".join(r["breakpoint_id"] for r in cluster),
                    "n_contributing_raw_calls": len(cluster),
                }
            )
            merged.append(best)
    return merged


def robustness_rank(value: object) -> int:
    return {"high": 2, "medium": 1, "low": 0}.get(str(value), -1)


def representative_sort_key(row: dict[str, object]) -> tuple[float, float, float, int, int]:
    support_bp = min(int(row.get("left_support_bp", 0) or 0), int(row.get("right_support_bp", 0) or 0))
    support_blocks = min(int(row.get("left_support_blocks", 0) or 0), int(row.get("right_support_blocks", 0) or 0))
    return (
        float(row.get("breakpoint_score", 0.0) or 0.0),
        robustness_rank(row.get("robustness", "")),
        float(row.get("n_contributing_raw_calls", 0) or 0),
        support_bp,
        support_blocks,
    )


def median_position(positions: list[int]) -> int:
    return int(round(statistics.median(positions)))


def cluster_consensus(calls: list[dict[str, object]], cross_species_bp: int) -> list[dict[str, object]]:
    rows = sorted(calls, key=lambda r: int(r["reference_position"]))
    clusters: list[list[dict[str, object]]] = []
    for row in rows:
        pos = int(row["reference_position"])
        if not clusters or pos - min(int(r["reference_position"]) for r in clusters[-1]) > cross_species_bp:
            clusters.append([row])
        else:
            clusters[-1].append(row)
    consensus = []
    for idx, cluster in enumerate(clusters, 1):
        by_species: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in cluster:
            by_species[str(row["species"])].append(row)
        retained = []
        removed = []
        for species_name in sorted(by_species, key=FOCAL_SPECIES.index):
            species_rows = by_species[species_name]
            best = max(species_rows, key=representative_sort_key)
            retained.append(best)
            removed.extend(r for r in species_rows if r is not best)
        positions = [int(r["reference_position"]) for r in retained]
        scores = [float(r["breakpoint_score"]) for r in retained]
        species = [str(r["species"]) for r in retained]
        evidence = sorted(set(",".join(str(r["evidence"]) for r in retained).split(",")))
        robust_count = len({r["species"] for r in retained if r.get("robustness") == "high"})
        if len(species) >= 4:
            confidence = "strong_consensus"
        elif 2 <= len(species) <= 3:
            confidence = "moderate_consensus"
        else:
            confidence = "species_specific"
        span = max(positions) - min(positions)
        assert span <= cross_species_bp, f"Consensus cluster exceeds span limit: {span} > {cross_species_bp}"
        consensus.append(
            {
                "consensus_breakpoint_id": f"DSB_{idx:03d}",
                "reference_position": median_position(positions),
                "min_position": min(positions),
                "max_position": max(positions),
                "cluster_span_bp": span,
                "n_species_support": len(species),
                "species_support": ",".join(species),
                "n_total_candidate_calls": len(cluster),
                "n_same_species_calls_removed": len(removed),
                "removed_same_species_call_ids": ";".join(str(r["breakpoint_id"]) for r in removed),
                "mean_score": statistics.mean(scores),
                "median_score": statistics.median(scores),
                "evidence_types": ",".join(e for e in evidence if e),
                "robust_species_count": robust_count,
                "confidence_class": confidence,
                "source_breakpoint_ids": ";".join(str(r["breakpoint_id"]) for r in retained),
                "candidate_breakpoint_ids": ";".join(str(r["breakpoint_id"]) for r in cluster),
            }
        )
    return consensus


def breakpoint_set_rows(consensus: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in consensus:
        n_species = int(row["n_species_support"])
        confidence = str(row["confidence_class"])
        rows.append(
            {
                "consensus_breakpoint_id": row["consensus_breakpoint_id"],
                "reference_position": row["reference_position"],
                "n_species_support": n_species,
                "confidence_class": confidence,
                "in_primary": int(n_species >= 2),
                "in_stringent": int(confidence == "strong_consensus"),
                "in_inclusive": 1,
            }
        )
    return rows


def preserve_stage27b_outputs() -> tuple[list[dict[str, str]], dict[str, object] | None]:
    current_consensus = read_tsv(CONSENSUS_OUT) if CONSENSUS_OUT.exists() else []
    current_manifest = json.loads(MANIFEST_OUT.read_text()) if MANIFEST_OUT.exists() else None
    if current_consensus and not STAGE27B_CONSENSUS_OUT.exists():
        shutil.copyfile(CONSENSUS_OUT, STAGE27B_CONSENSUS_OUT)
    if current_manifest is not None and not STAGE27B_MANIFEST_OUT.exists():
        shutil.copyfile(MANIFEST_OUT, STAGE27B_MANIFEST_OUT)
    old_consensus = read_tsv(STAGE27B_CONSENSUS_OUT) if STAGE27B_CONSENSUS_OUT.exists() else current_consensus
    old_manifest = json.loads(STAGE27B_MANIFEST_OUT.read_text()) if STAGE27B_MANIFEST_OUT.exists() else current_manifest
    return old_consensus, old_manifest


def audit_stage27b_vs_stage27c(old_consensus: list[dict[str, str]], new_consensus: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    matched_new: set[str] = set()
    old_sources = {
        row["consensus_breakpoint_id"]: set(filter(None, str(row.get("source_breakpoint_ids", "")).split(";")))
        for row in old_consensus
    }
    new_sources = {
        str(row["consensus_breakpoint_id"]): set(filter(None, str(row.get("candidate_breakpoint_ids", row.get("source_breakpoint_ids", ""))).split(";")))
        for row in new_consensus
    }
    for old in old_consensus:
        old_id = old["consensus_breakpoint_id"]
        old_pos = int(float(old["reference_position"]))
        overlaps = [
            new
            for new in new_consensus
            if old_sources[old_id] and old_sources[old_id] & new_sources[str(new["consensus_breakpoint_id"])]
        ]
        if not overlaps:
            rows.append(
                {
                    "old_breakpoint_id": old_id,
                    "old_position": old_pos,
                    "new_breakpoint_id": "",
                    "new_position": "",
                    "distance_bp": "",
                    "old_species_support": old.get("species_support", ""),
                    "new_species_support": "",
                    "status": "removed",
                }
            )
            continue
        status = "split" if len(overlaps) > 1 else "retained"
        for new in overlaps:
            new_id = str(new["consensus_breakpoint_id"])
            matched_new.add(new_id)
            new_pos = int(new["reference_position"])
            if status != "split" and new_pos != old_pos:
                status = "shifted"
            rows.append(
                {
                    "old_breakpoint_id": old_id,
                    "old_position": old_pos,
                    "new_breakpoint_id": new_id,
                    "new_position": new_pos,
                    "distance_bp": abs(new_pos - old_pos),
                    "old_species_support": old.get("species_support", ""),
                    "new_species_support": new.get("species_support", ""),
                    "status": status,
                }
            )
    for new in new_consensus:
        new_id = str(new["consensus_breakpoint_id"])
        if new_id in matched_new:
            continue
        rows.append(
            {
                "old_breakpoint_id": "",
                "old_position": "",
                "new_breakpoint_id": new_id,
                "new_position": new["reference_position"],
                "distance_bp": "",
                "old_species_support": "",
                "new_species_support": new.get("species_support", ""),
                "status": "newly_resolved",
            }
        )
    return rows


def make_segments(consensus: list[dict[str, object]], chrom_length: int) -> list[dict[str, object]]:
    sorted_bp = sorted(consensus, key=lambda r: int(r["reference_position"]))
    starts = [1] + [int(r["reference_position"]) + 1 for r in sorted_bp]
    ends = [int(r["reference_position"]) for r in sorted_bp] + [chrom_length]
    rows = []
    for idx, (start, end) in enumerate(zip(starts, ends), 1):
        rows.append(
            {
                "segment_id": f"DSEG_{idx:03d}",
                "start": start,
                "end": end,
                "length_bp": max(0, end - start + 1),
                "left_breakpoint_id": sorted_bp[idx - 2]["consensus_breakpoint_id"] if idx > 1 else "",
                "right_breakpoint_id": sorted_bp[idx - 1]["consensus_breakpoint_id"] if idx <= len(sorted_bp) else "",
            }
        )
    return rows


def sensitivity_analysis(blocks: list[Block], base_args: argparse.Namespace, frozen_consensus: list[dict[str, object]]) -> list[dict[str, object]]:
    ref_gaps = sorted({max(50_000, base_args.max_reference_gap // 2), base_args.max_reference_gap, base_args.max_reference_gap * 2})
    query_gaps = sorted({max(100_000, base_args.max_query_gap // 2), base_args.max_query_gap, base_args.max_query_gap * 2})
    flank_bps = sorted({max(50_000, base_args.min_flank_bp // 2), base_args.min_flank_bp, base_args.min_flank_bp * 2})
    flank_blocks = sorted({max(1, base_args.min_flank_blocks - 1), base_args.min_flank_blocks, base_args.min_flank_blocks + 1})
    within_vals = sorted({max(50_000, base_args.within_species_merge_bp // 2), base_args.within_species_merge_bp})
    cross_vals = sorted({max(100_000, base_args.cross_species_cluster_bp // 2), base_args.cross_species_cluster_bp, base_args.cross_species_cluster_bp * 2})
    setting_vals = ["lenient", "default", "stringent"]
    recovered: list[dict[str, object]] = []
    run_id = 0
    for rg in ref_gaps:
        for qg in query_gaps:
            runs = build_synteny_runs(blocks, rg, qg)
            for flank_bp in flank_bps:
                for flank_block in flank_blocks:
                    disruptions = build_raw_disruptions(runs, flank_bp, flank_block)
                    all_events = []
                    for setting in setting_vals:
                        cfg = DEFAULT_SETTINGS[setting]
                        all_events.extend(call_breakpoints(disruptions, setting, cfg["k_mad"], cfg["min_score"]))
                    for within in within_vals:
                        merged = merge_species_calls(all_events, within)
                        for cross in cross_vals:
                            run_id += 1
                            consensus = cluster_consensus(merged, cross)
                            for row in consensus:
                                out = row.copy()
                                out.update(
                                    {
                                        "sensitivity_run": run_id,
                                        "max_reference_gap": rg,
                                        "max_query_gap": qg,
                                        "min_flank_bp": flank_bp,
                                        "min_flank_blocks": flank_block,
                                        "within_species_merge_bp": within,
                                        "cross_species_cluster_bp": cross,
                                    }
                                )
                                recovered.append(out)
    return summarize_sensitivity(recovered, base_args.cross_species_cluster_bp, frozen_consensus)


def summarize_sensitivity(recovered: list[dict[str, object]], cluster_bp: int, frozen_consensus: list[dict[str, object]] | None = None) -> list[dict[str, object]]:
    n_runs = len({r["sensitivity_run"] for r in recovered})
    if frozen_consensus is not None:
        rows = []
        for bp in frozen_consensus:
            pos = int(bp["reference_position"])
            hits = [r for r in recovered if abs(int(r["reference_position"]) - pos) <= cluster_bp]
            runs = {r["sensitivity_run"] for r in hits}
            conf = Counter(r["confidence_class"] for r in hits)
            positions = [int(r["reference_position"]) for r in hits]
            species_counts = [int(r["n_species_support"]) for r in hits]
            stable_conf = sum(1 for r in hits if r["confidence_class"] == bp["confidence_class"])
            position_sd = statistics.pstdev(positions) if len(positions) > 1 else 0.0
            rows.append(
                {
                    "consensus_breakpoint_id": bp["consensus_breakpoint_id"],
                    "reference_position": pos,
                    "sensitivity_breakpoint_id": bp["consensus_breakpoint_id"],
                    "median_reference_position": median_position(positions) if hits else "",
                    "min_position": min(positions) if hits else "",
                    "max_position": max(positions) if hits else "",
                    "recovered_runs": len(runs),
                    "total_parameter_runs": n_runs,
                    "recovery_fraction": len(runs) / n_runs if n_runs else 0.0,
                    "position_sd_or_range": f"sd={position_sd:.6g};range={max(positions) - min(positions) if positions else ''}",
                    "median_species_support": statistics.median(species_counts) if species_counts else "",
                    "confidence_stability": stable_conf / len(hits) if hits else 0.0,
                    "modal_confidence_class": conf.most_common(1)[0][0] if conf else "",
                }
            )
        return rows
    clusters: list[list[dict[str, object]]] = []
    for row in sorted(recovered, key=lambda r: int(r["reference_position"])):
        if not clusters or int(row["reference_position"]) - int(clusters[-1][-1]["reference_position"]) > cluster_bp:
            clusters.append([row])
        else:
            clusters[-1].append(row)
    rows = []
    for idx, cluster in enumerate(clusters, 1):
        positions = [int(r["reference_position"]) for r in cluster]
        runs = {r["sensitivity_run"] for r in cluster}
        species_sets = Counter(r["species_support"] for r in cluster)
        conf = Counter(r["confidence_class"] for r in cluster).most_common(1)[0][0]
        rows.append(
            {
                "sensitivity_breakpoint_id": f"SENS_{idx:03d}",
                "median_reference_position": median_position(positions),
                "min_position": min(positions),
                "max_position": max(positions),
                "recovered_runs": len(runs),
                "total_parameter_runs": n_runs,
                "recovery_fraction": len(runs) / n_runs if n_runs else 0.0,
                "modal_species_support": species_sets.most_common(1)[0][0],
                "modal_confidence_class": conf,
            }
        )
    return rows


def plot_outputs(runs: list[Run], species_calls: list[dict[str, object]], consensus: list[dict[str, object]], chrom_length: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping figures because matplotlib is unavailable: {exc}", file=sys.stderr)
        return

    colors = defaultdict(lambda: "#737373")
    palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22"]
    for chrom in sorted({r.query_chrom for r in runs}):
        colors[chrom] = palette[len(colors) % len(palette)]
    call_by_species = defaultdict(list)
    for row in species_calls:
        call_by_species[row["species"]].append(row)

    fig, ax = plt.subplots(figsize=(13, 4.8))
    for y, species in enumerate(PLOT_ORDER):
        rs = [r for r in runs if r.species == species]
        for run in rs:
            height = 0.32
            ybase = y - height / 2
            ax.broken_barh([(run.reference_start, run.reference_end - run.reference_start)], (ybase, height), facecolors=colors[run.query_chrom], edgecolors="none")
            if run.relative_orientation == "-":
                ax.plot([run.reference_start, run.reference_end], [y, y], color="black", linewidth=0.35)
        for bp in call_by_species.get(species, []):
            ax.plot([int(bp["reference_position"]), int(bp["reference_position"])], [y - 0.42, y + 0.42], color="black", linewidth=0.7)
    ax.set_xlim(0, chrom_length)
    ax.set_ylim(-0.7, len(PLOT_ORDER) - 0.3)
    ax.set_yticks(range(len(PLOT_ORDER)))
    ax.set_yticklabels(PLOT_ORDER)
    ax.set_xlabel("GalGal6 chr4 coordinate (bp)")
    ax.set_title("De novo chr4 synteny runs and structural breakpoints")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(MAP_PDF)
    fig.savefig(MAP_PNG, dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 2.8))
    class_styles = {
        "strong_consensus": {"color": "#1b7837", "linewidth": 2.4, "linestyle": "-", "label": "strong"},
        "moderate_consensus": {"color": "#2166ac", "linewidth": 1.7, "linestyle": "-", "label": "moderate"},
        "species_specific": {"color": "#7f7f7f", "linewidth": 0.9, "linestyle": "--", "label": "species-specific"},
    }
    seen_labels = set()
    for row in consensus:
        pos = int(row["reference_position"])
        style = class_styles[str(row["confidence_class"])]
        label = style["label"] if style["label"] not in seen_labels else None
        seen_labels.add(style["label"])
        ax.axvline(pos, color=style["color"], linewidth=style["linewidth"], linestyle=style["linestyle"], alpha=0.9, label=label)
        if int(row["n_species_support"]) >= 2:
            ax.plot(pos, 0.45, marker="o", markersize=2.8, color=style["color"])
        ax.text(pos, 0.55, row["consensus_breakpoint_id"], rotation=90, va="bottom", ha="center", fontsize=7)
    ax.set_xlim(0, chrom_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("GalGal6 chr4 coordinate (bp)")
    ax.set_title("Consensus structural breakpoints inferred from chr4 synteny only")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(loc="upper right", frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(CONS_PDF)
    fig.savefig(CONS_PNG, dpi=220)
    plt.close(fig)


def ref_minus_fraction(blocks: list[Block]) -> dict[str, float]:
    out = {}
    for species in FOCAL_SPECIES:
        rows = [b for b in blocks if b.species == species]
        out[species] = (sum(1 for b in rows if b.reference_strand == "-") / len(rows)) if rows else 0.0
    return out


def invalid_orientation_run_counts(blocks: list[Block], max_reference_gap: int, max_query_gap: int) -> Counter:
    fixed = []
    for block in blocks:
        clone = Block(**block.__dict__)
        clone.relative_orientation = clone.query_strand
        fixed.append(clone)
    return Counter(r.species for r in build_synteny_runs(fixed, max_reference_gap, max_query_gap))


def write_manifest(
    args: argparse.Namespace,
    blocks: list[Block],
    runs: list[Run],
    disruptions: list[dict[str, object]],
    species_calls: list[dict[str, object]],
    consensus: list[dict[str, object]],
    sensitivity: list[dict[str, object]],
    gap_report: dict[str, object],
    old_run_counts: Counter,
) -> None:
    sensitivity_by_id = {r["consensus_breakpoint_id"]: r["recovery_fraction"] for r in sensitivity}
    classes = Counter(str(r["confidence_class"]) for r in consensus)
    sets = breakpoint_set_rows(consensus)
    manifest = {
        "stage": "027_denovo_structural_breakpoints",
        "version": "2.7c_final_frozen_structural_breakpoints",
        "freeze_status": "final frozen structural breakpoint manifest for Stage 3A",
        "git_commit": git_commit(),
        "code_version": git_commit(),
        "command": " ".join(sys.argv),
        "parameters": {
            "max_reference_gap": args.max_reference_gap,
            "max_query_gap": args.max_query_gap,
            "min_flank_bp": args.min_flank_bp,
            "min_flank_blocks": args.min_flank_blocks,
            "within_species_merge_bp": args.within_species_merge_bp,
            "cross_species_cluster_bp": args.cross_species_cluster_bp,
            "settings": DEFAULT_SETTINGS,
        },
        "flank_semantics": "raw_blocks",
        "cross_species_clustering": "complete-span constrained",
        "max_one_call_per_species": True,
        "consensus_coordinate": "median",
        "random_seed": None,
        "input_file": str(INPUT_PATH.relative_to(REPO_ROOT)),
        "input_sha256": file_sha256(INPUT_PATH),
        "forbidden_inference_inputs": sorted(FORBIDDEN_INFERENCE_INPUTS),
        "independence_statement": "All de novo structural breakpoint calls were inferred from complete totab_chr4.csv.xz blocks before any published outlier or genealogy table was read.",
        "empirical_gap_distributions": gap_report,
        "number_raw_blocks": len(blocks),
        "number_synteny_runs": len(runs),
        "raw_block_count_per_species": Counter(b.species for b in blocks),
        "reference_minus_fraction_per_species": ref_minus_fraction(blocks),
        "invalid_query_strand_orientation_run_count_per_species": old_run_counts,
        "synteny_run_count_per_species": Counter(r.species for r in runs),
        "raw_disruption_count_per_species": Counter(r["species"] for r in disruptions),
        "merged_breakpoint_count_per_species": Counter(r["species"] for r in species_calls),
        "total_consensus_breakpoint_count": len(consensus),
        "strong_consensus_count": classes["strong_consensus"],
        "moderate_consensus_count": classes["moderate_consensus"],
        "species_specific_count": classes["species_specific"],
        "primary_breakpoint_count": sum(int(r["in_primary"]) for r in sets),
        "stringent_breakpoint_count": sum(int(r["in_stringent"]) for r in sets),
        "inclusive_breakpoint_count": sum(int(r["in_inclusive"]) for r in sets),
        "consensus_breakpoint_coordinates": {r["consensus_breakpoint_id"]: r["reference_position"] for r in consensus},
        "final_breakpoint_coordinates": [r["reference_position"] for r in consensus],
        "consensus_sensitivity_recovery_fraction": sensitivity_by_id,
    }
    MANIFEST_OUT.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=dict) + "\n")


def git_commit() -> str:
    import subprocess

    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def run_pipeline(args: argparse.Namespace) -> dict[str, object]:
    old_consensus, _old_manifest = preserve_stage27b_outputs()
    if MANIFEST_OUT.exists() and not INVALID_MANIFEST_OUT.exists():
        old_manifest = json.loads(MANIFEST_OUT.read_text())
        old_manifest["invalid_reason"] = "Previous Stage 2.7 manifest used query strand alone as orientation; superseded by reference-relative orientation."
        INVALID_MANIFEST_OUT.write_text(json.dumps(old_manifest, indent=2, sort_keys=True, default=dict) + "\n")
    blocks, _ = parse_complete_blocks(INPUT_PATH)
    adj = adjacent_block_features(blocks)
    if args.max_reference_gap is None or args.max_query_gap is None:
        args.max_reference_gap, args.max_query_gap = choose_default_tolerances(adj, args.max_reference_gap, args.max_query_gap)
    gap_report = empirical_gap_report(adj)
    old_run_counts = invalid_orientation_run_counts(blocks, args.max_reference_gap, args.max_query_gap)
    runs = build_synteny_runs(blocks, args.max_reference_gap, args.max_query_gap)
    disruptions = build_raw_disruptions(runs, args.min_flank_bp, args.min_flank_blocks)
    event_calls = []
    for setting, cfg in DEFAULT_SETTINGS.items():
        event_calls.extend(call_breakpoints(disruptions, setting, cfg["k_mad"], cfg["min_score"]))
    species_calls = merge_species_calls(event_calls, args.within_species_merge_bp)
    consensus = cluster_consensus(species_calls, args.cross_species_cluster_bp)
    chrom_length = max(b.reference_end for b in blocks)
    segments = make_segments(consensus, chrom_length)
    sensitivity = sensitivity_analysis(blocks, args, consensus)
    sensitivity_by_id = {r["consensus_breakpoint_id"]: r["recovery_fraction"] for r in sensitivity}
    for row in consensus:
        row["parameter_recovery_fraction"] = sensitivity_by_id.get(row["consensus_breakpoint_id"], 0.0)
    sets = breakpoint_set_rows(consensus)
    audit = audit_stage27b_vs_stage27c(old_consensus, consensus)

    write_tsv(BLOCKS_OUT, block_rows(blocks), list(Block.__annotations__.keys()))
    write_tsv(ADJ_OUT, adj, list(adj[0].keys()) if adj else [])
    write_tsv(RUNS_OUT, run_rows(runs), list(Run.__annotations__.keys()))
    disruption_columns = [
        "species", "reference_position", "left_run_id", "right_run_id", "left_query_chrom",
        "right_query_chrom", "left_relative_orientation", "right_relative_orientation",
        "reference_gap_bp", "query_gap_bp", "same_query_chromosome", "relative_orientation_switch",
        "query_order_consistent", "query_order_reversal", "query_overlap", "query_jump_bp",
        "reference_query_gap_discrepancy", "scaffold_or_chromosome_switch",
        "score_chromosome_switch", "score_orientation_switch", "score_order_reversal",
        "score_query_jump", "score_gap_discrepancy", "breakpoint_score", "evidence",
        "persistent_regime_change", "left_support_bp", "right_support_bp",
        "left_support_blocks", "right_support_blocks",
    ]
    write_tsv(RAW_DISRUPTIONS_OUT, disruptions, disruption_columns)
    bp_columns = [
        "species", "breakpoint_id", "reference_position", "left_run_id", "right_run_id",
        "left_query_chrom", "right_query_chrom", "left_relative_orientation", "right_relative_orientation",
        "chromosome_switch", "orientation_switch", "order_reversal", "query_jump_bp",
        "gap_discrepancy", "persistent_regime_change", "left_support_bp", "right_support_bp",
        "left_support_blocks", "right_support_blocks", "score_chromosome_switch", "score_orientation_switch",
        "score_order_reversal", "score_query_jump", "score_gap_discrepancy",
        "breakpoint_score", "threshold", "threshold_setting", "robustness", "evidence",
        "contributing_raw_calls", "n_contributing_raw_calls",
    ]
    write_tsv(RAW_BP_OUT, event_calls, bp_columns)
    write_tsv(SPECIES_BP_OUT, species_calls, bp_columns)
    write_tsv(
        CONSENSUS_OUT,
        consensus,
        [
            "consensus_breakpoint_id", "reference_position", "min_position", "max_position",
            "cluster_span_bp", "n_species_support", "species_support", "n_total_candidate_calls",
            "n_same_species_calls_removed", "removed_same_species_call_ids", "mean_score", "median_score",
            "evidence_types", "robust_species_count", "parameter_recovery_fraction",
            "confidence_class", "source_breakpoint_ids", "candidate_breakpoint_ids",
        ],
    )
    write_tsv(
        BREAKPOINT_SETS_OUT,
        sets,
        [
            "consensus_breakpoint_id", "reference_position", "n_species_support",
            "confidence_class", "in_primary", "in_stringent", "in_inclusive",
        ],
    )
    write_tsv(SEGMENTS_OUT, segments, ["segment_id", "start", "end", "length_bp", "left_breakpoint_id", "right_breakpoint_id"])
    write_tsv(
        SENSITIVITY_OUT,
        sensitivity,
        [
            "consensus_breakpoint_id", "reference_position", "sensitivity_breakpoint_id", "median_reference_position", "min_position", "max_position",
            "recovered_runs", "total_parameter_runs", "recovery_fraction", "position_sd_or_range",
            "median_species_support", "confidence_stability", "modal_confidence_class",
        ],
    )
    write_tsv(
        STAGE27_AUDIT_OUT,
        audit,
        [
            "old_breakpoint_id", "old_position", "new_breakpoint_id", "new_position",
            "distance_bp", "old_species_support", "new_species_support", "status",
        ],
    )
    plot_outputs(runs, species_calls, consensus, chrom_length)
    write_manifest(args, blocks, runs, disruptions, species_calls, consensus, sensitivity, gap_report, old_run_counts)
    return {
        "blocks": blocks,
        "runs": runs,
        "disruptions": disruptions,
        "raw_calls": event_calls,
        "species_calls": species_calls,
        "consensus": consensus,
        "sensitivity": sensitivity,
    }


class SyntheticTests(unittest.TestCase):
    def make_blocks(self, specs: list[tuple[int, int, str, int, int, str, str, str]]) -> list[Block]:
        rows = []
        for i, (rs, re, qc, qs, qe, ref_strand, query_strand, _label) in enumerate(specs, 1):
            rows.append(
                Block(
                    "Stork", "x", str(i), "chr4", rs, re, ref_strand, "ref", rs, re,
                    qc, qs, qe, query_strand, "q", qs, qe, qc,
                    relative_orientation(ref_strand, query_strand), re - rs
                )
            )
        return rows

    def test_relative_orientation_all_strand_combinations(self):
        self.assertEqual(relative_orientation("+", "+"), "+")
        self.assertEqual(relative_orientation("+", "-"), "-")
        self.assertEqual(relative_orientation("-", "+"), "-")
        self.assertEqual(relative_orientation("-", "-"), "+")

    def test_ref_minus_query_plus_relative_minus(self):
        block = self.make_blocks([(1, 100, "q1", 900, 1000, "-", "+", "")])[0]
        self.assertEqual(block.relative_orientation, "-")

    def test_ref_minus_query_minus_relative_plus(self):
        block = self.make_blocks([(1, 100, "q1", 1, 100, "-", "-", "")])[0]
        self.assertEqual(block.relative_orientation, "+")

    def default_calls(self, blocks: list[Block], min_flank_bp: int = 200, min_flank_blocks: int = 2) -> tuple[list[Run], list[dict[str, object]], list[dict[str, object]]]:
        runs = build_synteny_runs(blocks, 100_000, 250_000)
        disruptions = build_raw_disruptions(runs, min_flank_bp, min_flank_blocks)
        calls = call_breakpoints(disruptions, "default", 3.5, 3.5)
        return runs, disruptions, calls

    def test_collinear_no_breakpoint(self):
        blocks = self.make_blocks([(1, 100, "q1", 1, 100, "+", "+", ""), (150, 250, "q1", 150, 250, "+", "+", "")])
        runs, disruptions, calls = self.default_calls(blocks)
        self.assertEqual(len(runs), 1)
        self.assertEqual(disruptions, [])
        self.assertEqual(calls, [])

    def test_fragmented_reverse_collinearity_collapses(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 901, 1000, "+", "-", ""),
            (150, 250, "q1", 751, 850, "+", "-", ""),
            (300, 400, "q1", 601, 700, "+", "-", ""),
        ])
        runs, _, calls = self.default_calls(blocks)
        self.assertEqual(len(runs), 1)
        self.assertEqual(runs[0].relative_orientation, "-")
        self.assertEqual(calls, [])

    def test_isolated_short_orientation_flip_raw_not_event(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 330, "q1", 900, 930, "+", "-", ""),
            (400, 500, "q1", 400, 500, "+", "+", ""),
            (550, 650, "q1", 550, 650, "+", "+", ""),
        ])
        runs = build_synteny_runs(blocks, 100_000, 250_000)
        disruptions = build_raw_disruptions(runs, 200, 2)
        calls = call_breakpoints(disruptions, "default", 3.5, 3.5)
        self.assertGreaterEqual(len(disruptions), 2)
        self.assertEqual(calls, [])

    def test_long_persistent_inversion_two_event_boundaries(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 400, "q1", 1000, 1100, "+", "-", ""),
            (450, 550, "q1", 850, 950, "+", "-", ""),
            (600, 700, "q1", 600, 700, "+", "+", ""),
            (750, 850, "q1", 750, 850, "+", "+", ""),
        ])
        _, _, calls = self.default_calls(blocks, min_flank_bp=150, min_flank_blocks=2)
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(int(c["persistent_regime_change"]) for c in calls))

    def test_persistent_chromosome_switch_retained(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 400, "q2", 300, 400, "+", "+", ""),
            (450, 550, "q2", 450, 550, "+", "+", ""),
        ])
        _, _, calls = self.default_calls(blocks, min_flank_bp=150, min_flank_blocks=2)
        self.assertEqual(len(calls), 1)
        self.assertEqual(int(calls[0]["chromosome_switch"]), 1)

    def test_local_query_coordinate_noise_no_breakpoint(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 170, 270, "+", "+", ""),
            (300, 400, "q1", 300, 400, "+", "+", ""),
        ])
        runs, _, calls = self.default_calls(blocks)
        self.assertEqual(len(runs), 1)
        self.assertEqual(calls, [])

    def test_large_query_jump_breakpoint_with_persistence(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 400, "q1", 2_000_000, 2_000_100, "+", "+", ""),
            (450, 550, "q1", 2_000_150, 2_000_250, "+", "+", ""),
        ])
        runs = build_synteny_runs(blocks, 100_000, 250_000)
        disruptions = build_raw_disruptions(runs, 150, 2)
        calls = call_breakpoints(disruptions, "default", 1.0, 1.0)
        self.assertEqual(len(calls), 1)
        self.assertIn("large_query_jump", calls[0]["evidence"])

    def test_flank_support_counts_raw_blocks_not_synteny_runs(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 400, "q2", 300, 400, "+", "+", ""),
            (450, 550, "q2", 450, 550, "+", "+", ""),
        ])
        runs = build_synteny_runs(blocks, 100_000, 250_000)
        disruptions = build_raw_disruptions(runs, 150, 2)
        self.assertEqual(len(runs), 2)
        self.assertEqual(disruptions[0]["left_support_blocks"], 2)
        self.assertEqual(disruptions[0]["right_support_blocks"], 2)
        self.assertEqual(disruptions[0]["persistent_regime_change"], 1)

    def test_one_long_run_with_multiple_blocks_satisfies_flank_support(self):
        blocks = self.make_blocks([
            (1, 100, "q1", 1, 100, "+", "+", ""),
            (150, 250, "q1", 150, 250, "+", "+", ""),
            (300, 400, "q1", 300, 400, "+", "+", ""),
            (450, 550, "q2", 450, 550, "+", "+", ""),
            (600, 700, "q2", 600, 700, "+", "+", ""),
            (750, 850, "q2", 750, 850, "+", "+", ""),
        ])
        runs = build_synteny_runs(blocks, 100_000, 250_000)
        disruptions = build_raw_disruptions(runs, 250, 3)
        self.assertEqual(len(runs), 2)
        self.assertEqual(disruptions[0]["left_support_blocks"], 3)
        self.assertEqual(disruptions[0]["right_support_blocks"], 3)
        self.assertEqual(disruptions[0]["persistent_regime_change"], 1)

    def test_fragmented_collinearity_collapses(self):
        blocks = self.make_blocks([(i * 100, i * 100 + 50, "q1", i * 100, i * 100 + 50, "+", "+", "") for i in range(1, 20)])
        runs, _, calls = self.default_calls(blocks)
        self.assertEqual(len(runs), 1)
        self.assertEqual(calls, [])

    def test_nearby_species_calls_merge(self):
        raw = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 1000, "left_run_id": "l1", "right_run_id": "r1", "left_run_index": 1, "right_run_index": 2, "breakpoint_score": 5, "threshold_setting": "default", "evidence": "chromosome_switch"},
            {"species": "Stork", "breakpoint_id": "b", "reference_position": 1050, "left_run_id": "l2", "right_run_id": "r2", "left_run_index": 2, "right_run_index": 3, "breakpoint_score": 4, "threshold_setting": "default", "evidence": "chromosome_switch"},
        ]
        merged = merge_species_calls(raw, 100)
        self.assertEqual(len(merged), 1)

    def test_consensus_clustering_nearby_species_calls(self):
        calls = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 1000, "breakpoint_score": 5, "evidence": "chromosome_switch", "robustness": "high"},
            {"species": "Flamingo", "breakpoint_id": "b", "reference_position": 1100, "breakpoint_score": 6, "evidence": "orientation_switch", "robustness": "high"},
        ]
        consensus = cluster_consensus(calls, 200)
        self.assertEqual(len(consensus), 1)
        self.assertEqual(consensus[0]["n_species_support"], 2)

    def test_consensus_clustering_is_complete_span_constrained(self):
        calls = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 0, "breakpoint_score": 5, "evidence": "chromosome_switch", "robustness": "high"},
            {"species": "Flamingo", "breakpoint_id": "b", "reference_position": 150_000, "breakpoint_score": 6, "evidence": "orientation_switch", "robustness": "high"},
            {"species": "Dove", "breakpoint_id": "c", "reference_position": 300_000, "breakpoint_score": 7, "evidence": "orientation_switch", "robustness": "high"},
        ]
        consensus = cluster_consensus(calls, 200_000)
        self.assertEqual(len(consensus), 2)
        self.assertTrue(all(int(c["cluster_span_bp"]) <= 200_000 for c in consensus))

    def test_same_species_dedup_and_unique_species_support(self):
        calls = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 1000, "breakpoint_score": 5, "evidence": "chromosome_switch", "robustness": "low", "left_support_bp": 1000, "right_support_bp": 1000, "left_support_blocks": 2, "right_support_blocks": 2},
            {"species": "Stork", "breakpoint_id": "b", "reference_position": 1050, "breakpoint_score": 8, "evidence": "chromosome_switch", "robustness": "high", "left_support_bp": 1000, "right_support_bp": 1000, "left_support_blocks": 2, "right_support_blocks": 2},
            {"species": "Flamingo", "breakpoint_id": "c", "reference_position": 1100, "breakpoint_score": 6, "evidence": "orientation_switch", "robustness": "high", "left_support_bp": 1000, "right_support_bp": 1000, "left_support_blocks": 2, "right_support_blocks": 2},
        ]
        consensus = cluster_consensus(calls, 200)
        self.assertEqual(consensus[0]["n_species_support"], 2)
        self.assertEqual(consensus[0]["n_total_candidate_calls"], 3)
        self.assertEqual(consensus[0]["n_same_species_calls_removed"], 1)
        self.assertEqual(consensus[0]["removed_same_species_call_ids"], "a")
        self.assertEqual(consensus[0]["source_breakpoint_ids"], "b;c")

    def test_median_consensus_position(self):
        calls = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 1000, "breakpoint_score": 5, "evidence": "chromosome_switch", "robustness": "high"},
            {"species": "Flamingo", "breakpoint_id": "b", "reference_position": 1200, "breakpoint_score": 50, "evidence": "orientation_switch", "robustness": "high"},
            {"species": "Dove", "breakpoint_id": "c", "reference_position": 1600, "breakpoint_score": 6, "evidence": "orientation_switch", "robustness": "high"},
        ]
        consensus = cluster_consensus(calls, 1000)
        self.assertEqual(consensus[0]["reference_position"], 1200)

    def test_breakpoint_set_membership_is_deterministic(self):
        consensus = [
            {"consensus_breakpoint_id": "a", "reference_position": 1, "n_species_support": 4, "confidence_class": "strong_consensus"},
            {"consensus_breakpoint_id": "b", "reference_position": 2, "n_species_support": 2, "confidence_class": "moderate_consensus"},
            {"consensus_breakpoint_id": "c", "reference_position": 3, "n_species_support": 1, "confidence_class": "species_specific"},
        ]
        rows = breakpoint_set_rows(consensus)
        self.assertEqual([(r["in_primary"], r["in_stringent"], r["in_inclusive"]) for r in rows], [(1, 1, 1), (1, 0, 1), (0, 0, 1)])

    def test_species_specific_events_separate(self):
        calls = [
            {"species": "Stork", "breakpoint_id": "a", "reference_position": 1000, "breakpoint_score": 5, "evidence": "chromosome_switch", "robustness": "high"},
            {"species": "Flamingo", "breakpoint_id": "b", "reference_position": 5000, "breakpoint_score": 6, "evidence": "orientation_switch", "robustness": "high"},
        ]
        consensus = cluster_consensus(calls, 200)
        self.assertEqual(len(consensus), 2)
        self.assertTrue(all(c["confidence_class"] == "species_specific" for c in consensus))

    def test_forbidden_input_guard(self):
        with self.assertRaises(RuntimeError):
            assert_inference_inputs_allowed([REPO_ROOT / "empirical/neoaves_chr4/results/genealogy_change_points.tsv"])
        with self.assertRaises(RuntimeError):
            assert_inference_inputs_allowed([REPO_ROOT / "empirical/neoaves_chr4/results/locus_table_chr4.tsv"])
        with self.assertRaises(RuntimeError):
            assert_inference_inputs_allowed([REPO_ROOT / "rearrangements/outlier-regions-by-maf2synteny.tsv"])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-reference-gap", type=int, default=None, help="Maximum gap for merging collinear raw blocks; default is max(100 kb, empirical p95).")
    parser.add_argument("--max-query-gap", type=int, default=None, help="Maximum query jump for merging collinear raw blocks; default is max(250 kb, empirical p95).")
    parser.add_argument("--min-flank-bp", type=int, default=200_000)
    parser.add_argument("--min-flank-blocks", type=int, default=2, help="Minimum raw maf2synteny blocks supporting each flank of an event-scale breakpoint.")
    parser.add_argument("--within-species-merge-bp", type=int, default=100_000)
    parser.add_argument("--cross-species-cluster-bp", type=int, default=200_000)
    parser.add_argument("--run-tests", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SyntheticTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            return 1
    result = run_pipeline(args)
    print("De novo structural breakpoint inference complete.")
    print(f"Raw blocks: {len(result['blocks'])}")
    print(f"Synteny runs: {len(result['runs'])}")
    print(f"Raw structural disruptions: {len(result['disruptions'])}")
    print(f"Merged species breakpoints: {len(result['species_calls'])}")
    print(f"Consensus breakpoints: {len(result['consensus'])}")
    print(f"Manifest: {MANIFEST_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
