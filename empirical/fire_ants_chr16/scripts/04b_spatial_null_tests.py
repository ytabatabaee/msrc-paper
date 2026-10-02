#!/usr/bin/env python3
"""Stage 4B fire-ant chromosome-16 spatial-null inference.

This stage tests whether the Stage-4A shift toward SB/Sb haplotype-partition
quartet support is unusually aligned with the independently frozen chromosome-16
supergene region. It performs the predeclared exact circular-shift null and a
deterministic physical-coordinate-aware same-width interval sensitivity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import unittest
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


getcontext().prec = 60

DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"

STAGE0_FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
STAGE1_MANIFEST = EMPIRICAL_ROOT / "results" / "stage1_manifest.json"
STAGE2_MANIFEST = EMPIRICAL_ROOT / "results" / "stage2_manifest.json"
STAGE2_FOCAL_PARTITION = DATA_ROOT / "processed" / "stage2_focal_partition.tsv"
STAGE3_MANIFEST = EMPIRICAL_ROOT / "results" / "stage3_manifest.json"
STAGE3_BACKGROUND_PARTITION = DATA_ROOT / "processed" / "stage3_background_partition.tsv"
STAGE4A_MANIFEST = EMPIRICAL_ROOT / "results" / "stage4a_manifest.json"
STAGE4A_WINDOW_SUPPORT = DATA_ROOT / "processed" / "stage4a_window_quartet_support.tsv"
STAGE4A_PRIMARY_CONTRAST = EMPIRICAL_ROOT / "results" / "stage4a_primary_contrast.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"

STAGE1_TOPOLOGIES = DATA_ROOT / "processed" / "stage1_twisst_topologies.tsv"
STAGE1_WEIGHTS_SPARSE = DATA_ROOT / "processed" / "stage1_twisst_weights_sparse.tsv"
STAGE1_WEIGHT_SUMMARY = DATA_ROOT / "processed" / "stage1_twisst_weight_summary.tsv"
STAGE1_WINDOW_INDEX = DATA_ROOT / "processed" / "stage1_window_index.tsv"

CHR16_PHYSICAL_ORDER = EMPIRICAL_ROOT / "results" / "stage4b_chr16_physical_order.tsv"
CIRCULAR_SHIFTS = EMPIRICAL_ROOT / "results" / "stage4b_circular_shifts.tsv"
PRIMARY_TEST = EMPIRICAL_ROOT / "results" / "stage4b_primary_test.tsv"
COORDINATE_PLACEMENTS = EMPIRICAL_ROOT / "results" / "stage4b_coordinate_placements.tsv"
COORDINATE_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4b_coordinate_summary.tsv"
COORDINATE_LENGTH_WEIGHTED = EMPIRICAL_ROOT / "results" / "stage4b_coordinate_length_weighted.tsv"
COORDINATE_LENGTH_WEIGHTED_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4b_coordinate_length_weighted_summary.tsv"
BOUNDARY_CONTEXT = EMPIRICAL_ROOT / "results" / "stage4b_boundary_context.tsv"
MAIN_TABLE = EMPIRICAL_ROOT / "results" / "stage4b_main_table.tsv"
REPORT = EMPIRICAL_ROOT / "results" / "stage4b_report.md"
MANIFEST = EMPIRICAL_ROOT / "results" / "stage4b_manifest.json"
README = EMPIRICAL_ROOT / "README.md"

CIRCULAR_NULL_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_circular_null.pdf"
CIRCULAR_NULL_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_circular_null.png"
COORDINATE_NULL_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_coordinate_null.pdf"
COORDINATE_NULL_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_coordinate_null.png"
COORDINATE_LENGTH_WEIGHTED_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_coordinate_length_weighted.pdf"
COORDINATE_LENGTH_WEIGHTED_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_coordinate_length_weighted.png"
SPATIAL_TEST_PDF = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_spatial_test.pdf"
SPATIAL_TEST_PNG = EMPIRICAL_ROOT / "figures" / "fire_ants_chr16_spatial_test.png"

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
EXPECTED_STAGE4A_WINDOW_SUPPORT_SHA = "b8ff3ea95ba31ae5553f56f0aeb45a393bdb0ba2b19081fc1aca822dfdd70758"
EXPECTED_STAGE4A_MANIFEST_SHA = "d1fa79815bbe4d1f45a2af7f5a1d5487d39fe3041ccac4f79f961a998ab99bf1"
EXPECTED_REGION_COUNTS = {"chr1": 117, "chr16A": 42, "chr16B": 2, "chr16_supergene": 52}
EXPECTED_CLASS_COUNTS = {"species": 315, "haplotype": 315, "third": 315}
EXPECTED_DELTA_D = Decimal("1.610109863495346")
EXPECTED_LENGTH_WEIGHTED_P = Decimal("0.07832901837")

DECIMAL_PLACES = Decimal("0.000000000000000")
P_DECIMAL_PLACES = Decimal("0.0000000000")
TSV_LINETERMINATOR = "\n"
TIE_TOLERANCE = Decimal("0.000000000001")


class Stage4BError(RuntimeError):
    """Raised when Stage-4B validation or inference fails."""


@dataclass(frozen=True)
class AnalysisResults:
    window_rows: list[dict[str, object]]
    chr16_rows: list[dict[str, object]]
    region_spans: dict[str, dict[str, object]]
    observed: dict[str, Decimal]
    circular_rows: list[dict[str, object]]
    primary_row: dict[str, object]
    coordinate_rows: list[dict[str, object]]
    coordinate_summary: dict[str, object]
    coordinate_length_weighted_rows: list[dict[str, object]]
    coordinate_length_weighted_summary: dict[str, object]
    boundary_rows: list[dict[str, object]]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def fmt_decimal(value: Decimal) -> str:
    return str(value.quantize(DECIMAL_PLACES))


def fmt_p(value: Decimal) -> str:
    return str(value.quantize(P_DECIMAL_PLACES))


def fmt_bool(value: bool) -> str:
    return "true" if value else "false"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            out = {}
            for field in fields:
                value = row.get(field, "")
                if isinstance(value, bool):
                    out[field] = fmt_bool(value)
                elif value is None:
                    out[field] = ""
                else:
                    out[field] = str(value)
            writer.writerow(out)


def mean_decimal(values: list[Decimal]) -> Decimal:
    if not values:
        raise Stage4BError("cannot take mean of an empty set")
    return sum(values, Decimal("0")) / Decimal(len(values))


def delta_for_rows(inside: list[dict[str, object]], outside: list[dict[str, object]], field: str = "D") -> tuple[Decimal, Decimal, Decimal]:
    mean_inside = mean_decimal([row[field] for row in inside])
    mean_outside = mean_decimal([row[field] for row in outside])
    return mean_inside - mean_outside, mean_inside, mean_outside


def require_sha(path: Path, expected: str, label: str) -> None:
    observed = sha256(path)
    if observed != expected:
        raise Stage4BError(f"{label} checksum changed: expected {expected}, observed {observed}")


def verify_recorded_checksums(entries: dict[str, str], label: str) -> None:
    for path_text, expected in sorted(entries.items()):
        path = REPO_ROOT / path_text
        if not path.exists():
            raise Stage4BError(f"{label} recorded path is missing: {path_text}")
        require_sha(path, expected, f"{label}: {path_text}")


def verify_prior_stages() -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    require_sha(STAGE0_FOCAL_QUARTET, EXPECTED_STAGE0_FOCAL_SHA, "Stage-0 focal quartet")
    require_sha(STAGE1_MANIFEST, EXPECTED_STAGE1_MANIFEST_SHA, "Stage-1 manifest")
    require_sha(STAGE2_MANIFEST, EXPECTED_STAGE2_MANIFEST_SHA, "Stage-2 manifest")
    require_sha(STAGE2_FOCAL_PARTITION, EXPECTED_STAGE2_FOCAL_PARTITION_SHA, "Stage-2 focal partition")
    require_sha(STAGE3_MANIFEST, EXPECTED_STAGE3_MANIFEST_SHA, "Stage-3 manifest")
    require_sha(STAGE3_BACKGROUND_PARTITION, EXPECTED_STAGE3_BACKGROUND_PARTITION_SHA, "Stage-3 background partition")
    require_sha(STAGE1_TOPOLOGIES, EXPECTED_STAGE1_TOPOLOGIES_SHA, "Stage-1 topology table")
    require_sha(STAGE1_WEIGHTS_SPARSE, EXPECTED_STAGE1_WEIGHTS_SHA, "Stage-1 sparse weights")
    require_sha(STAGE1_WINDOW_INDEX, EXPECTED_STAGE1_WINDOW_INDEX_SHA, "Stage-1 window index")
    require_sha(STAGE1_WEIGHT_SUMMARY, EXPECTED_STAGE1_WEIGHT_SUMMARY_SHA, "Stage-1 weight summary")
    require_sha(STAGE4A_WINDOW_SUPPORT, EXPECTED_STAGE4A_WINDOW_SUPPORT_SHA, "Stage-4A window support")
    require_sha(STAGE4A_MANIFEST, EXPECTED_STAGE4A_MANIFEST_SHA, "Stage-4A manifest")

    stage1 = json.loads(STAGE1_MANIFEST.read_text())
    expected_metrics = {
        "coordinate_rows": "213",
        "local_tree_rows": "213",
        "min_tree_tips": "267",
        "max_tree_tips": "267",
        "n_distinct_tip_sets": "1",
        "twisst_topologies": "945",
        "unique_twisst_topologies": "945",
        "weight_rows": "213",
        "weight_columns": "945",
        "topology_crosscheck_mismatches": "0",
        "chr1_windows": "117",
        "chr16A_windows": "42",
        "chr16B_windows": "2",
        "chr16_supergene_windows": "52",
    }
    for key, expected in expected_metrics.items():
        observed = str(stage1.get("observed_metrics", {}).get(key))
        if observed != expected:
            raise Stage4BError(f"Stage-1 metric {key} expected {expected}, observed {observed}")
    verify_recorded_checksums(stage1.get("raw_topology_file_checksums", {}), "Stage-1 raw")
    verify_recorded_checksums(stage1.get("processed_output_checksums", {}), "Stage-1 processed")

    stage2 = json.loads(STAGE2_MANIFEST.read_text())
    verify_recorded_checksums(stage2.get("stage2_output_checksums", {}), "Stage-2 output")

    stage3 = json.loads(STAGE3_MANIFEST.read_text())
    if stage3.get("stage3_background_partition_matches_stage0_species_split") is not True:
        raise Stage4BError("Stage-3 background partition no longer matches Stage-0 species split")
    stage3_recorded = {
        "data/fire_ants_chr16/processed/stage3_background_partition.tsv": stage3["stage3_background_partition_checksum"],
        "data/fire_ants_chr16/metadata/background_tree_provenance.tsv": stage3["background_tree_provenance_checksum"],
        "data/fire_ants_chr16/processed/stage3_background_tree.nwk": stage3["normalized_background_tree_checksum"],
    }
    verify_recorded_checksums(stage3_recorded, "Stage-3 output")

    stage4a = json.loads(STAGE4A_MANIFEST.read_text())
    if stage4a.get("inferential_tests_performed") is not False:
        raise Stage4BError("Stage-4A manifest no longer records inferential_tests_performed=false")
    if stage4a.get("p_values_calculated") is not False:
        raise Stage4BError("Stage-4A manifest no longer records p_values_calculated=false")
    if stage4a.get("class_counts") != EXPECTED_CLASS_COUNTS:
        raise Stage4BError(f"Stage-4A class counts changed: {stage4a.get('class_counts')}")
    if stage4a.get("n_windows") != 213:
        raise Stage4BError("Stage-4A window count changed")
    if stage4a.get("region_counts") != EXPECTED_REGION_COUNTS:
        raise Stage4BError(f"Stage-4A region counts changed: {stage4a.get('region_counts')}")
    verify_recorded_checksums(stage4a.get("stage4a_output_checksums", {}), "Stage-4A output")
    verify_recorded_checksums(stage4a.get("figure_checksums", {}), "Stage-4A figure")

    return stage1, stage2, stage3, stage4a


def load_window_support() -> list[dict[str, object]]:
    fields, rows = read_tsv(STAGE4A_WINDOW_SUPPORT)
    required = {"window_index", "chrom", "start", "end", "mid", "region", "q_S", "q_H", "q_3", "D", "dominant_class"}
    missing = required.difference(fields)
    if missing:
        raise Stage4BError(f"Stage-4A window support missing fields: {sorted(missing)}")
    parsed: list[dict[str, object]] = []
    for row in rows:
        parsed.append(
            {
                **row,
                "window_index": int(row["window_index"]),
                "start": Decimal(row["start"]),
                "end": Decimal(row["end"]),
                "mid": Decimal(row["mid"]),
                "q_S": Decimal(row["q_S"]),
                "q_H": Decimal(row["q_H"]),
                "q_3": Decimal(row["q_3"]),
                "D": Decimal(row["D"]),
            }
        )
    if len(parsed) != 213:
        raise Stage4BError(f"expected 213 Stage-4A windows, observed {len(parsed)}")
    if Counter(row["region"] for row in parsed) != Counter(EXPECTED_REGION_COUNTS):
        raise Stage4BError("Stage-4A window region counts do not match frozen counts")
    return parsed


def load_region_spans() -> dict[str, dict[str, object]]:
    fields, rows = read_tsv(REGION_MANIFEST)
    required = {"region_id", "chromosome", "coordinate_start", "coordinate_end", "coordinate_definition"}
    missing = required.difference(fields)
    if missing:
        raise Stage4BError(f"region manifest missing fields: {sorted(missing)}")
    spans: dict[str, dict[str, object]] = {}
    for row in rows:
        spans[row["region_id"]] = {
            **row,
            "coordinate_start": Decimal(row["coordinate_start"]),
            "coordinate_end": Decimal(row["coordinate_end"]),
        }
    for region in ["chr16A", "chr16B", "chr16_supergene"]:
        if region not in spans:
            raise Stage4BError(f"region manifest missing {region}")
    return spans


def verify_stage4a_contrast(window_rows: list[dict[str, object]]) -> dict[str, Decimal]:
    supergene = [row for row in window_rows if row["region"] == "chr16_supergene"]
    outside = [row for row in window_rows if row["region"] in {"chr16A", "chr16B"}]
    if len(supergene) != 52 or len(outside) != 44:
        raise Stage4BError("Stage-4A primary chr16 contrast has wrong supergene/outside counts")
    delta_d, mean_supergene, mean_outside = delta_for_rows(supergene, outside, "D")
    delta_q_s, mean_q_s_supergene, mean_q_s_outside = delta_for_rows(supergene, outside, "q_S")
    delta_q_h, mean_q_h_supergene, mean_q_h_outside = delta_for_rows(supergene, outside, "q_H")
    delta_q_3, mean_q_3_supergene, mean_q_3_outside = delta_for_rows(supergene, outside, "q_3")
    if delta_d.quantize(DECIMAL_PLACES) != EXPECTED_DELTA_D:
        raise Stage4BError(f"Stage-4A delta_D expected {EXPECTED_DELTA_D}, observed {fmt_decimal(delta_d)}")

    _, contrast_rows = read_tsv(STAGE4A_PRIMARY_CONTRAST)
    primary = next((row for row in contrast_rows if row["comparison"] == "chr16_supergene_vs_chr16_outside"), None)
    if primary is None:
        raise Stage4BError("Stage-4A primary contrast row is missing")
    if Decimal(primary["delta_D"]) != delta_d.quantize(DECIMAL_PLACES):
        raise Stage4BError("Stage-4A primary contrast does not reproduce from frozen window table")

    return {
        "delta_D": delta_d,
        "mean_D_supergene": mean_supergene,
        "mean_D_outside": mean_outside,
        "delta_q_S": delta_q_s,
        "mean_q_S_supergene": mean_q_s_supergene,
        "mean_q_S_outside": mean_q_s_outside,
        "delta_q_H": delta_q_h,
        "mean_q_H_supergene": mean_q_h_supergene,
        "mean_q_H_outside": mean_q_h_outside,
        "delta_q_3": delta_q_3,
        "mean_q_3_supergene": mean_q_3_supergene,
        "mean_q_3_outside": mean_q_3_outside,
    }


def physical_chr16_rows(window_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    chr16 = sorted([row for row in window_rows if row["chrom"] == "chr16"], key=lambda row: (row["mid"], row["window_index"]))
    if len(chr16) != 96:
        raise Stage4BError(f"expected 96 chr16 windows, observed {len(chr16)}")
    if any(chr16[i]["mid"] > chr16[i + 1]["mid"] for i in range(len(chr16) - 1)):
        raise Stage4BError("physical chr16 rows are not sorted by increasing midpoint")
    region_counts = Counter(row["region"] for row in chr16)
    expected_counts = Counter({"chr16A": 42, "chr16_supergene": 52, "chr16B": 2})
    if region_counts != expected_counts:
        raise Stage4BError(f"chr16 physical region counts changed: {region_counts}")
    compressed = []
    for row in chr16:
        if not compressed or compressed[-1] != row["region"]:
            compressed.append(row["region"])
    if compressed != ["chr16A", "chr16_supergene", "chr16B"]:
        raise Stage4BError(f"physical chr16 region sequence is not A -> supergene -> B: {compressed}")
    original_chr16 = [row for row in window_rows if row["chrom"] == "chr16"]
    if [row["window_index"] for row in original_chr16] == [row["window_index"] for row in chr16]:
        raise Stage4BError("Stage 4B would use Stage-1 concatenation order; physical midpoint sort did not change order")
    return chr16


def write_physical_order(chr16_rows: list[dict[str, object]]) -> None:
    rows = []
    for rank, row in enumerate(chr16_rows, start=1):
        rows.append(
            {
                "physical_rank": rank,
                "window_index": row["window_index"],
                "mid": str(row["mid"]),
                "region": row["region"],
                "D": fmt_decimal(row["D"]),
            }
        )
    write_tsv(CHR16_PHYSICAL_ORDER, rows, ["physical_rank", "window_index", "mid", "region", "D"])


def calculate_circular_null(chr16_rows: list[dict[str, object]], observed_delta: Decimal) -> tuple[list[dict[str, object]], dict[str, object]]:
    d_values = [row["D"] for row in chr16_rows]
    q_s_values = [row["q_S"] for row in chr16_rows]
    q_h_values = [row["q_H"] for row in chr16_rows]
    q_3_values = [row["q_3"] for row in chr16_rows]
    mask = [row["region"] == "chr16_supergene" for row in chr16_rows]
    if sum(mask) != 52 or len(mask) - sum(mask) != 44:
        raise Stage4BError("circular mask does not have 52 supergene and 44 outside positions")

    rows: list[dict[str, object]] = []
    n = len(chr16_rows)
    observed_delta_q_s = observed_delta_q_h = observed_delta_q_3 = None
    for shift in range(n):
        shifted_d = [d_values[(i - shift) % n] for i in range(n)]
        shifted_q_s = [q_s_values[(i - shift) % n] for i in range(n)]
        shifted_q_h = [q_h_values[(i - shift) % n] for i in range(n)]
        shifted_q_3 = [q_3_values[(i - shift) % n] for i in range(n)]
        inside_d = [value for value, is_inside in zip(shifted_d, mask, strict=True) if is_inside]
        outside_d = [value for value, is_inside in zip(shifted_d, mask, strict=True) if not is_inside]
        mean_inside = mean_decimal(inside_d)
        mean_outside = mean_decimal(outside_d)
        delta_d = mean_inside - mean_outside
        delta_q_s = mean_decimal([value for value, is_inside in zip(shifted_q_s, mask, strict=True) if is_inside]) - mean_decimal(
            [value for value, is_inside in zip(shifted_q_s, mask, strict=True) if not is_inside]
        )
        delta_q_h = mean_decimal([value for value, is_inside in zip(shifted_q_h, mask, strict=True) if is_inside]) - mean_decimal(
            [value for value, is_inside in zip(shifted_q_h, mask, strict=True) if not is_inside]
        )
        delta_q_3 = mean_decimal([value for value, is_inside in zip(shifted_q_3, mask, strict=True) if is_inside]) - mean_decimal(
            [value for value, is_inside in zip(shifted_q_3, mask, strict=True) if not is_inside]
        )
        exceeds = delta_d >= observed_delta - TIE_TOLERANCE
        if shift == 0:
            observed_delta_q_s = delta_q_s
            observed_delta_q_h = delta_q_h
            observed_delta_q_3 = delta_q_3
        rows.append(
            {
                "shift": shift,
                "is_observed_alignment": shift == 0,
                "delta_D": fmt_decimal(delta_d),
                "mean_D_supergene_mask": fmt_decimal(mean_inside),
                "mean_D_outside_mask": fmt_decimal(mean_outside),
                "delta_q_S": fmt_decimal(delta_q_s),
                "delta_q_H": fmt_decimal(delta_q_h),
                "delta_q_3": fmt_decimal(delta_q_3),
                "exceeds_or_equals_observed": exceeds,
            }
        )
    if len(rows) != 96:
        raise Stage4BError("circular null did not evaluate exactly 96 shifts")
    if Decimal(rows[0]["delta_D"]) != observed_delta.quantize(DECIMAL_PLACES):
        raise Stage4BError("circular shift 0 does not reproduce observed delta_D")
    sorted_deltas = sorted((Decimal(row["delta_D"]) for row in rows), reverse=True)
    observed_rank = 1 + sum(value > observed_delta + TIE_TOLERANCE for value in sorted_deltas)
    n_ge = sum(1 for row in rows if row["exceeds_or_equals_observed"])
    p_exact = Decimal(n_ge) / Decimal(96)
    primary = {
        "test": "exact_circular_shift_null",
        "observed_delta_D": fmt_decimal(observed_delta),
        "n_exact_alignments": 96,
        "observed_rank": observed_rank,
        "n_ge_observed": n_ge,
        "p_one_sided": fmt_p(p_exact),
        "alternative": "Delta_D > 0; D higher inside the independently frozen supergene region",
        "null_description": "exact 96-alignment circular-shift null rotating the physical chr16 D track while keeping the frozen 52-window supergene mask fixed",
        "comparison_tolerance": str(TIE_TOLERANCE),
        "observed_delta_q_S": fmt_decimal(observed_delta_q_s),
        "observed_delta_q_H": fmt_decimal(observed_delta_q_h),
        "observed_delta_q_3": fmt_decimal(observed_delta_q_3),
    }
    return rows, primary


def membership_key(rows: list[dict[str, object]], start: Decimal, width: Decimal) -> tuple[int, ...]:
    end = start + width
    return tuple(row["window_index"] for row in rows if start <= row["mid"] <= end)


def coordinate_domain_events(
    chr16_rows: list[dict[str, object]],
    region_spans: dict[str, dict[str, object]],
) -> dict[str, object]:
    supergene_span = region_spans["chr16_supergene"]
    observed_start = supergene_span["coordinate_start"]
    observed_end = supergene_span["coordinate_end"]
    width = observed_end - observed_start
    if width != Decimal("16237060"):
        raise Stage4BError(f"unexpected observed supergene width {width}")

    domain_start = min(row["start"] for row in chr16_rows)
    domain_end = max(row["end"] for row in chr16_rows)
    min_start = domain_start
    max_start = domain_end - width
    if not (min_start <= observed_start <= max_start):
        raise Stage4BError("observed supergene interval is outside the chr16 analysis domain")

    mids = [row["mid"] for row in chr16_rows]
    raw_events = sorted(set(mids + [mid - width for mid in mids]))
    valid_events = [value for value in raw_events if min_start <= value <= max_start]
    breakpoints = sorted(set([min_start, max_start, *valid_events]))
    if breakpoints[0] != min_start or breakpoints[-1] != max_start:
        raise Stage4BError("coordinate event partition does not include both start-domain endpoints")
    return {
        "observed_start": observed_start,
        "observed_end": observed_end,
        "width": width,
        "domain_start": domain_start,
        "domain_end": domain_end,
        "min_start": min_start,
        "max_start": max_start,
        "mids": mids,
        "raw_events": raw_events,
        "valid_events": valid_events,
        "breakpoints": breakpoints,
    }


def coordinate_null(
    chr16_rows: list[dict[str, object]],
    region_spans: dict[str, dict[str, object]],
    observed_delta: Decimal,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    domain = coordinate_domain_events(chr16_rows, region_spans)
    observed_start = domain["observed_start"]
    observed_end = domain["observed_end"]
    width = domain["width"]
    domain_start = domain["domain_start"]
    domain_end = domain["domain_end"]
    min_start = domain["min_start"]
    max_start = domain["max_start"]
    event_starts = domain["raw_events"]
    valid_events = domain["valid_events"]
    breakpoints = domain["breakpoints"]
    representative_starts: dict[Decimal, str] = {min_start: "domain_start", max_start: "domain_end_minus_width", observed_start: "observed_interval"}
    for value in valid_events:
        representative_starts[value] = "boundary_event"
    for left, right in zip(breakpoints, breakpoints[1:], strict=False):
        if left < right:
            representative_starts[(left + right) / Decimal("2")] = "between_boundary_events"

    valid_starts = sorted(value for value in representative_starts if min_start <= value <= max_start)
    dedup: dict[tuple[int, ...], dict[str, object]] = {}
    for start in valid_starts:
        key = membership_key(chr16_rows, start, width)
        if not key or len(key) == len(chr16_rows):
            continue
        interval_end = start + width
        existing = dedup.get(key)
        is_observed_key = key == membership_key(chr16_rows, observed_start, width)
        row = {
            "representative_start": start,
            "representative_start_kind": representative_starts[start],
            "interval_start": observed_start if is_observed_key else start,
            "interval_end": observed_end if is_observed_key else interval_end,
            "candidate_start_values": [start],
            "candidate_start_kinds": [representative_starts[start]],
            "is_observed_interval": is_observed_key,
        }
        if existing is None:
            dedup[key] = row
        else:
            existing["candidate_start_values"].append(start)
            existing["candidate_start_kinds"].append(representative_starts[start])
            if is_observed_key:
                existing["interval_start"] = observed_start
                existing["interval_end"] = observed_end
                existing["representative_start"] = observed_start
                existing["representative_start_kind"] = "observed_interval"
                existing["is_observed_interval"] = True

    placement_rows: list[dict[str, object]] = []
    rows_by_index = {row["window_index"]: row for row in chr16_rows}
    for placement_id, (key, audit) in enumerate(
        sorted(dedup.items(), key=lambda item: (item[1]["interval_start"], item[0])), start=1
    ):
        inside = [rows_by_index[index] for index in key]
        outside = [row for row in chr16_rows if row["window_index"] not in set(key)]
        delta, mean_inside, mean_outside = delta_for_rows(inside, outside, "D")
        placement_rows.append(
            {
                "placement_id": placement_id,
                "interval_start": str(audit["interval_start"]),
                "interval_end": str(audit["interval_end"]),
                "interval_width": str(width),
                "n_inside_windows": len(inside),
                "n_outside_windows": len(outside),
                "inside_window_indices": ";".join(str(index) for index in key),
                "delta_D": fmt_decimal(delta),
                "mean_D_inside": fmt_decimal(mean_inside),
                "mean_D_outside": fmt_decimal(mean_outside),
                "is_observed_interval": audit["is_observed_interval"],
                "exceeds_or_equals_observed": delta >= observed_delta - TIE_TOLERANCE,
                "representative_start_kind": audit["representative_start_kind"],
                "n_representative_starts_collapsed": len(audit["candidate_start_values"]),
                "collapsed_start_min": str(min(audit["candidate_start_values"])),
                "collapsed_start_max": str(max(audit["candidate_start_values"])),
                "collapsed_start_kinds": ";".join(sorted(set(audit["candidate_start_kinds"]))),
            }
        )

    observed_rows = [row for row in placement_rows if row["is_observed_interval"]]
    if len(observed_rows) != 1:
        raise Stage4BError(f"expected exactly one observed coordinate placement, observed {len(observed_rows)}")
    observed_key = set(int(value) for value in observed_rows[0]["inside_window_indices"].split(";"))
    observed_inside = [row for row in chr16_rows if row["window_index"] in observed_key]
    if Counter(row["region"] for row in observed_inside) != Counter({"chr16_supergene": 52}):
        raise Stage4BError("observed coordinate interval does not select exactly 52 chr16_supergene windows")
    if Decimal(observed_rows[0]["delta_D"]) != observed_delta.quantize(DECIMAL_PLACES):
        raise Stage4BError("observed coordinate interval does not reproduce observed delta_D")

    sorted_deltas = sorted((Decimal(row["delta_D"]) for row in placement_rows), reverse=True)
    observed_rank = 1 + sum(value > observed_delta + TIE_TOLERANCE for value in sorted_deltas)
    n_ge = sum(1 for row in placement_rows if row["exceeds_or_equals_observed"])
    p_coordinate = Decimal(n_ge) / Decimal(len(placement_rows))
    summary = {
        "test": "unique_membership_state_coordinate_sensitivity",
        "observed_delta_D": fmt_decimal(observed_delta),
        "coordinate_width_source": rel(REGION_MANIFEST),
        "observed_interval_start": str(observed_start),
        "observed_interval_end": str(observed_end),
        "interval_width": str(width),
        "analysis_domain_start": str(domain_start),
        "analysis_domain_end": str(domain_end),
        "analysis_domain_definition": "minimum source window start and maximum source window end among frozen chr16 Stage-4A windows",
        "candidate_boundary_events": len(event_starts),
        "valid_boundary_events": len(valid_events),
        "valid_starts_considered": len(valid_starts),
        "n_unique_membership_states": len(placement_rows),
        "rank_unique_membership": observed_rank,
        "n_ge_observed": n_ge,
        "p_unique_membership": fmt_p(p_coordinate),
        "alternative": "Delta_D > 0",
        "enumeration_rule": "unique candidate starts induced by s=x_i and s=x_i-L boundary events, chr16 domain boundaries, observed interval, and deterministic between-event representatives; exact inside-window sets deduplicated",
        "weighting": "equal weight per distinct sampled-window membership state, regardless of the physical start-coordinate length producing that state",
        "arbitrary_grid_spacing_used": False,
    }
    return placement_rows, summary


def coordinate_length_weighted_null(
    chr16_rows: list[dict[str, object]],
    region_spans: dict[str, dict[str, object]],
    observed_delta: Decimal,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    domain = coordinate_domain_events(chr16_rows, region_spans)
    width = domain["width"]
    min_start = domain["min_start"]
    max_start = domain["max_start"]
    breakpoints = domain["breakpoints"]
    total_length = max_start - min_start
    if total_length <= 0:
        raise Stage4BError("coordinate start domain has non-positive length")

    rows_by_index = {row["window_index"]: row for row in chr16_rows}
    segment_rows: list[dict[str, object]] = []
    extreme_length = Decimal("0")
    for segment_id, (left, right) in enumerate(zip(breakpoints, breakpoints[1:]), start=1):
        length = right - left
        if length <= 0:
            raise Stage4BError("coordinate event partition contains a non-positive segment")
        representative_start = (left + right) / Decimal("2")
        key = membership_key(chr16_rows, representative_start, width)
        if not key or len(key) == len(chr16_rows):
            raise Stage4BError("continuous coordinate segment has an invalid inside/outside split")
        inside = [rows_by_index[index] for index in key]
        outside = [row for row in chr16_rows if row["window_index"] not in set(key)]
        delta, mean_inside, mean_outside = delta_for_rows(inside, outside, "D")
        exceeds = delta >= observed_delta - TIE_TOLERANCE
        if exceeds:
            extreme_length += length
        segment_rows.append(
            {
                "segment_id": segment_id,
                "start_min": str(left),
                "start_max": str(right),
                "segment_length": str(length),
                "representative_start": str(representative_start),
                "n_inside_windows": len(inside),
                "n_outside_windows": len(outside),
                "inside_window_indices": ";".join(str(index) for index in key),
                "delta_D": fmt_decimal(delta),
                "mean_D_inside": fmt_decimal(mean_inside),
                "mean_D_outside": fmt_decimal(mean_outside),
                "exceeds_or_equals_observed": exceeds,
            }
        )

    summed_length = sum((Decimal(row["segment_length"]) for row in segment_rows), Decimal("0"))
    if summed_length != total_length:
        raise Stage4BError(f"coordinate segments sum to {summed_length}, expected {total_length}")
    p_length_weighted = extreme_length / total_length
    if abs(p_length_weighted - EXPECTED_LENGTH_WEIGHTED_P) > Decimal("0.00000000001"):
        raise Stage4BError(f"length-weighted p expected near {EXPECTED_LENGTH_WEIGHTED_P}, observed {p_length_weighted}")
    summary = {
        "observed_delta_D": fmt_decimal(observed_delta),
        "interval_width": str(width),
        "start_domain_min": str(min_start),
        "start_domain_max": str(max_start),
        "total_start_domain_length": str(total_length),
        "n_constant_membership_segments": len(segment_rows),
        "extreme_start_domain_length": str(extreme_length),
        "p_length_weighted": fmt_p(p_length_weighted),
        "alternative": "Delta_D > 0",
        "null_description": "continuous coordinate sensitivity with interval start uniformly distributed on the valid physical start-coordinate domain; constant-membership segments weighted by physical start-coordinate length",
        "event_rule": "events are generated only from chr16 window midpoints x_i, x_i-L, and the start-domain endpoints",
        "arbitrary_grid_spacing_used": False,
    }
    return segment_rows, summary


def boundary_context(chr16_rows: list[dict[str, object]], region_spans: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    boundaries = {
        "left": region_spans["chr16_supergene"]["coordinate_start"],
        "right": region_spans["chr16_supergene"]["coordinate_end"],
    }
    out: list[dict[str, object]] = []
    for name, coordinate in boundaries.items():
        left = [row for row in chr16_rows if row["mid"] < coordinate]
        right = [row for row in chr16_rows if row["mid"] > coordinate]
        selected = [("left", row) for row in left[-5:]] + [("right", row) for row in right[:5]]
        rank_by_window = {row["window_index"]: rank for rank, row in enumerate(chr16_rows, start=1)}
        for side, row in selected:
            out.append(
                {
                    "boundary": name,
                    "side": side,
                    "physical_rank": rank_by_window[row["window_index"]],
                    "window_index": row["window_index"],
                    "mid": str(row["mid"]),
                    "distance_to_boundary": str(abs(row["mid"] - coordinate)),
                    "region": row["region"],
                    "q_S": fmt_decimal(row["q_S"]),
                    "q_H": fmt_decimal(row["q_H"]),
                    "q_3": fmt_decimal(row["q_3"]),
                    "D": fmt_decimal(row["D"]),
                    "dominant_class": row["dominant_class"],
                }
            )
    return out


def write_primary_outputs(results: AnalysisResults) -> None:
    write_physical_order(results.chr16_rows)
    write_tsv(
        CIRCULAR_SHIFTS,
        results.circular_rows,
        [
            "shift",
            "is_observed_alignment",
            "delta_D",
            "mean_D_supergene_mask",
            "mean_D_outside_mask",
            "delta_q_S",
            "delta_q_H",
            "delta_q_3",
            "exceeds_or_equals_observed",
        ],
    )
    write_tsv(
        PRIMARY_TEST,
        [results.primary_row],
        [
            "test",
            "observed_delta_D",
            "n_exact_alignments",
            "observed_rank",
            "n_ge_observed",
            "p_one_sided",
            "alternative",
            "null_description",
            "comparison_tolerance",
            "observed_delta_q_S",
            "observed_delta_q_H",
            "observed_delta_q_3",
        ],
    )
    write_tsv(
        COORDINATE_PLACEMENTS,
        results.coordinate_rows,
        [
            "placement_id",
            "interval_start",
            "interval_end",
            "interval_width",
            "n_inside_windows",
            "n_outside_windows",
            "inside_window_indices",
            "delta_D",
            "mean_D_inside",
            "mean_D_outside",
            "is_observed_interval",
            "exceeds_or_equals_observed",
            "representative_start_kind",
            "n_representative_starts_collapsed",
            "collapsed_start_min",
            "collapsed_start_max",
            "collapsed_start_kinds",
        ],
    )
    write_tsv(
        COORDINATE_SUMMARY,
        [results.coordinate_summary],
        [
            "test",
            "observed_delta_D",
            "coordinate_width_source",
            "observed_interval_start",
            "observed_interval_end",
            "interval_width",
            "analysis_domain_start",
            "analysis_domain_end",
            "analysis_domain_definition",
            "candidate_boundary_events",
            "valid_boundary_events",
            "valid_starts_considered",
            "n_unique_membership_states",
            "rank_unique_membership",
            "n_ge_observed",
            "p_unique_membership",
            "alternative",
            "enumeration_rule",
            "weighting",
            "arbitrary_grid_spacing_used",
        ],
    )
    write_tsv(
        COORDINATE_LENGTH_WEIGHTED,
        results.coordinate_length_weighted_rows,
        [
            "segment_id",
            "start_min",
            "start_max",
            "segment_length",
            "representative_start",
            "n_inside_windows",
            "n_outside_windows",
            "inside_window_indices",
            "delta_D",
            "mean_D_inside",
            "mean_D_outside",
            "exceeds_or_equals_observed",
        ],
    )
    write_tsv(
        COORDINATE_LENGTH_WEIGHTED_SUMMARY,
        [results.coordinate_length_weighted_summary],
        [
            "observed_delta_D",
            "interval_width",
            "start_domain_min",
            "start_domain_max",
            "total_start_domain_length",
            "n_constant_membership_segments",
            "extreme_start_domain_length",
            "p_length_weighted",
            "alternative",
            "null_description",
        ],
    )
    write_tsv(
        BOUNDARY_CONTEXT,
        results.boundary_rows,
        [
            "boundary",
            "side",
            "physical_rank",
            "window_index",
            "mid",
            "distance_to_boundary",
            "region",
            "q_S",
            "q_H",
            "q_3",
            "D",
            "dominant_class",
        ],
    )
    main_row = {
        "comparison": "chr16_supergene_vs_chr16_outside",
        "n_supergene": 52,
        "n_outside": 44,
        "mean_D_supergene": fmt_decimal(results.observed["mean_D_supergene"]),
        "mean_D_outside": fmt_decimal(results.observed["mean_D_outside"]),
        "delta_D": fmt_decimal(results.observed["delta_D"]),
        "circular_rank": results.primary_row["observed_rank"],
        "circular_n": results.primary_row["n_exact_alignments"],
        "circular_p": results.primary_row["p_one_sided"],
        "coordinate_unique_rank": results.coordinate_summary["rank_unique_membership"],
        "coordinate_unique_n": results.coordinate_summary["n_unique_membership_states"],
        "coordinate_unique_p": results.coordinate_summary["p_unique_membership"],
        "coordinate_length_weighted_p": results.coordinate_length_weighted_summary["p_length_weighted"],
        "interpretation": "positive_spatial_alignment",
    }
    write_tsv(
        MAIN_TABLE,
        [main_row],
        [
            "comparison",
            "n_supergene",
            "n_outside",
            "mean_D_supergene",
            "mean_D_outside",
            "delta_D",
            "circular_rank",
            "circular_n",
            "circular_p",
            "coordinate_unique_rank",
            "coordinate_unique_n",
            "coordinate_unique_p",
            "coordinate_length_weighted_p",
            "interpretation",
        ],
    )


def plot_circular_null(circular_rows: list[dict[str, object]], observed_delta: Decimal) -> None:
    rows = sorted(circular_rows, key=lambda row: Decimal(row["delta_D"]))
    x = list(range(1, len(rows) + 1))
    y = [float(Decimal(row["delta_D"])) for row in rows]
    observed_y = float(observed_delta)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    colors = ["#b23a48" if row["is_observed_alignment"] else "#3d5a80" for row in rows]
    ax.scatter(x, y, s=26, color=colors, linewidth=0, alpha=0.9)
    ax.axhline(observed_y, color="#b23a48", linewidth=1.2, label="observed Delta_D")
    ax.set_xlabel("ranked circular alignment")
    ax.set_ylabel("Delta_D")
    ax.set_title("Exact 96-alignment circular-shift null")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    CIRCULAR_NULL_PDF.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(CIRCULAR_NULL_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(CIRCULAR_NULL_PNG, dpi=220)
    plt.close(fig)


def plot_coordinate_null(coordinate_rows: list[dict[str, object]], observed_delta: Decimal) -> None:
    rows = sorted(coordinate_rows, key=lambda row: Decimal(row["interval_start"]))
    x = [(float(Decimal(row["interval_start"])) + float(Decimal(row["interval_end"]))) / 2_000_000 for row in rows]
    y = [float(Decimal(row["delta_D"])) for row in rows]
    colors = ["#b23a48" if row["is_observed_interval"] else "#4c78a8" for row in rows]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.scatter(x, y, s=24, color=colors, linewidth=0, alpha=0.9)
    ax.axhline(float(observed_delta), color="#b23a48", linewidth=1.2, label="observed supergene interval")
    ax.set_xlabel("candidate interval midpoint on chr16 (Mb)")
    ax.set_ylabel("Delta_D")
    ax.set_title("Physical-coordinate same-width interval sensitivity")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    with PdfPages(COORDINATE_NULL_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(COORDINATE_NULL_PNG, dpi=220)
    plt.close(fig)


def plot_coordinate_length_weighted(
    segment_rows: list[dict[str, object]],
    observed_delta: Decimal,
    region_spans: dict[str, dict[str, object]],
) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.3))
    for row in segment_rows:
        start = float(Decimal(row["start_min"])) / 1_000_000
        end = float(Decimal(row["start_max"])) / 1_000_000
        y = float(Decimal(row["delta_D"]))
        color = "#b23a48" if row["exceeds_or_equals_observed"] else "#4c78a8"
        ax.hlines(y, start, end, color=color, linewidth=2.0, alpha=0.88)
    observed_start = float(region_spans["chr16_supergene"]["coordinate_start"]) / 1_000_000
    ax.axhline(float(observed_delta), color="#b23a48", linewidth=1.2, label="observed Delta_D")
    ax.axvline(observed_start, color="#222222", linewidth=1.0, linestyle=":", label="observed supergene start")
    ax.set_xlabel("interval start coordinate on chr16 (Mb)")
    ax.set_ylabel("Delta_D")
    ax.set_title("Uniform physical-start coordinate sensitivity")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    with PdfPages(COORDINATE_LENGTH_WEIGHTED_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(COORDINATE_LENGTH_WEIGHTED_PNG, dpi=220)
    plt.close(fig)


def plot_spatial_summary(chr16_rows: list[dict[str, object]], circular_rows: list[dict[str, object]], region_spans: dict[str, dict[str, object]], observed_delta: Decimal) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.0))
    ax = axes[0]
    x = [float(row["mid"]) / 1_000_000 for row in chr16_rows]
    y = [float(row["D"]) for row in chr16_rows]
    ax.plot(x, y, color="#3d5a80", marker="o", markersize=3.2, linewidth=0.9)
    start = float(region_spans["chr16_supergene"]["coordinate_start"]) / 1_000_000
    end = float(region_spans["chr16_supergene"]["coordinate_end"]) / 1_000_000
    ax.axvspan(start, end, color="#d9d9d9", alpha=0.45, label="frozen supergene span")
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_ylim(-1.08, 1.08)
    ax.set_xlabel("chr16 physical midpoint (Mb)")
    ax.set_ylabel("D = q_H - q_S")
    ax.set_title("A. Observed chr16 D track")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1]
    rows = sorted(circular_rows, key=lambda row: Decimal(row["delta_D"]))
    ax.scatter(range(1, len(rows) + 1), [float(Decimal(row["delta_D"])) for row in rows], s=24, color="#3d5a80", linewidth=0)
    ax.axhline(float(observed_delta), color="#b23a48", linewidth=1.2, label="observed Delta_D")
    ax.set_xlabel("ranked circular alignment")
    ax.set_ylabel("Delta_D")
    ax.set_title("B. Exact 96-alignment circular-shift null")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    with PdfPages(SPATIAL_TEST_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(SPATIAL_TEST_PNG, dpi=220)
    plt.close(fig)


def make_figures(results: AnalysisResults) -> None:
    plot_circular_null(results.circular_rows, results.observed["delta_D"])
    plot_coordinate_null(results.coordinate_rows, results.observed["delta_D"])
    plot_coordinate_length_weighted(results.coordinate_length_weighted_rows, results.observed["delta_D"], results.region_spans)
    plot_spatial_summary(results.chr16_rows, results.circular_rows, results.region_spans, results.observed["delta_D"])


def write_report(results: AnalysisResults) -> None:
    dominance_counts = Counter(row["dominant_class"] for row in results.chr16_rows if row["region"] == "chr16_supergene")
    outside_dominance_counts = Counter(row["dominant_class"] for row in results.chr16_rows if row["region"] in {"chr16A", "chr16B"})
    boundary_preview = "\n".join(
        f"- {row['boundary']} {row['side']}: rank {row['physical_rank']}, window {row['window_index']}, region {row['region']}, mid {row['mid']}, D {row['D']}, dominant {row['dominant_class']}"
        for row in results.boundary_rows
    )
    REPORT.write_text(
        "# Fire ants chromosome 16 Stage 4B report\n\n"
        "## Purpose\n\n"
        "Stage 4B tests whether the Stage-4A shift from frozen species-history quartet support toward SB/Sb haplotype-partition support is unusually aligned with the independently frozen chromosome-16 supergene region.\n\n"
        "## Frozen inputs\n\n"
        f"- Stage-0 focal quartet checksum: `{sha256(STAGE0_FOCAL_QUARTET)}`\n"
        f"- Stage-1 manifest checksum: `{sha256(STAGE1_MANIFEST)}`\n"
        f"- Stage-2 focal partition checksum: `{sha256(STAGE2_FOCAL_PARTITION)}`\n"
        f"- Stage-3 background partition checksum: `{sha256(STAGE3_BACKGROUND_PARTITION)}`\n"
        f"- Stage-4A manifest checksum: `{sha256(STAGE4A_MANIFEST)}`\n"
        f"- Stage-4A window-support checksum: `{sha256(STAGE4A_WINDOW_SUPPORT)}`\n\n"
        "All recorded Stage-1 raw and processed checksums, Stage-2 outputs, Stage-3 outputs, and Stage-4A outputs were rechecked before inference.\n\n"
        "## Physical chromosome-16 ordering\n\n"
        "The Stage-1/TWISST source files concatenate chr16 regions as `chr16A`, `chr16B`, then `chr16_supergene`. Stage 4B does not use that row order for spatial inference. It reconstructs the chromosome-16 track by `chrom == chr16` and increasing physical midpoint.\n\n"
        "The verified physical order is `chr16A -> chr16_supergene -> chr16B`, with 96 total chr16 windows: 42 `chr16A`, 52 `chr16_supergene`, and 2 `chr16B`.\n\n"
        "## Observed Stage-4A contrast\n\n"
        f"The predeclared response is `D(w)=q_H(w)-q_S(w)`. The observed contrast is `Delta_D = mean(D_supergene) - mean(D_chr16 outside) = {fmt_decimal(results.observed['delta_D'])}`.\n\n"
        f"Component changes are descriptive: `Delta q_H = {fmt_decimal(results.observed['delta_q_H'])}`, `Delta q_S = {fmt_decimal(results.observed['delta_q_S'])}`, and `Delta q_3 = {fmt_decimal(results.observed['delta_q_3'])}`.\n\n"
        "## Primary exact circular-shift null\n\n"
        "The primary null rotates the complete ordered chr16 `D` vector through all 96 circular alignments while keeping the frozen 52-window supergene mask fixed. This preserves 52 inside and 44 outside windows, all observed `D` values, and the ordered local genealogy track. It destroys only the alignment between that track and the frozen supergene annotation. Ties to the observed statistic are counted with tolerance `1e-12`.\n\n"
        f"Observed rank among 96 shifts: `{results.primary_row['observed_rank']}`. `n_ge_observed = {results.primary_row['n_ge_observed']}`. Exact one-sided `p = {results.primary_row['p_one_sided']}`.\n\n"
        "## Physical-coordinate sensitivity\n\n"
        f"The same-width interval uses the frozen Stage-0 supergene span from `{rel(REGION_MANIFEST)}`: start `{results.coordinate_summary['observed_interval_start']}`, end `{results.coordinate_summary['observed_interval_end']}`, width `{results.coordinate_summary['interval_width']}` bp. The chr16 analysis domain is defined conservatively as the minimum source window start and maximum source window end among frozen chr16 Stage-4A windows: `{results.coordinate_summary['analysis_domain_start']}` to `{results.coordinate_summary['analysis_domain_end']}`.\n\n"
        "Secondary sensitivity 1 is the original unique-membership-state coordinate sensitivity. It enumerates starts induced by fixed-width interval boundary events `s=x_i` and `s=x_i-L`, plus domain boundaries, the observed interval, and deterministic between-event representatives, then deduplicates exact inside-window membership sets. Each distinct sampled-window membership pattern receives equal weight regardless of how much physical start-coordinate range produces it.\n\n"
        f"Unique membership states: `{results.coordinate_summary['n_unique_membership_states']}`. Observed rank: `{results.coordinate_summary['rank_unique_membership']}`. `n_ge_observed = {results.coordinate_summary['n_ge_observed']}`. Equal-weight unique-membership `p = {results.coordinate_summary['p_unique_membership']}`.\n\n"
        "Secondary sensitivity 2 is the continuous physical-coordinate null. It partitions the allowed interval-start domain exactly at events generated only by `x_i`, `x_i-L`, and the domain endpoints. Open intervals between adjacent events are weighted by their physical start-coordinate length, which corresponds to drawing the interval start uniformly over the valid physical domain.\n\n"
        f"Continuous start domain: `{results.coordinate_length_weighted_summary['start_domain_min']}` to `{results.coordinate_length_weighted_summary['start_domain_max']}` bp, total length `{results.coordinate_length_weighted_summary['total_start_domain_length']}` bp. Extreme start-coordinate length: `{results.coordinate_length_weighted_summary['extreme_start_domain_length']}` bp. Length-weighted physical-coordinate `p = {results.coordinate_length_weighted_summary['p_length_weighted']}`.\n\n"
        "No arbitrary coordinate grid is used in either coordinate sensitivity. These coordinate analyses are secondary sensitivities, not the primary P-value.\n\n"
        "## Boundary context\n\n"
        "Nearest-window context around the independently frozen supergene boundaries is reported without fitting a change point or moving boundaries:\n\n"
        f"{boundary_preview}\n\n"
        "## Results\n\n"
        f"Primary: exact 96-alignment circular shift, observed rank `{results.primary_row['observed_rank']}`, `p = {results.primary_row['p_one_sided']}`. Secondary sensitivity 1: equal-weight unique membership states, `p = {results.coordinate_summary['p_unique_membership']}`. Secondary sensitivity 2: continuous uniform physical interval start, `p = {results.coordinate_length_weighted_summary['p_length_weighted']}`.\n\n"
        "The primary window-space circular null places the observed supergene alignment second among 96 possible alignments. A coordinate-aware analysis remains supportive but is more conservative when candidate placements are weighted by the amount of physical start-coordinate space producing each sampled-window membership pattern.\n\n"
        "Dominance-pattern summary is descriptive only. Chr16 outside windows are "
        f"{outside_dominance_counts['species']}/44 species-dominant. Chr16 supergene windows are species={dominance_counts['species']}, haplotype={dominance_counts['haplotype']}, third={dominance_counts['third']}.\n\n"
        "## Interpretation\n\n"
        "The shift from species-history quartet support toward SB/Sb haplotype-partition support is unusually aligned with the independently defined chromosome-16 supergene region. The supergene region is spatially associated with a coherent alternative genealogy regime relative to collinear chromosome-16 windows.\n\n"
        "## Introgression caveat\n\n"
        "The source study's recurrent-introgression interpretation remains an explicit caveat. These spatial nulls do not prove MSRC, do not show that inversions caused the genealogy, and do not rule out introgression.\n\n"
        "## Limitations\n\n"
        "The primary null is in window space. The coordinate-aware sensitivity preserves physical width but changes the number of sampled windows in candidate intervals. Both analyses use the Stage-4A TWISST-derived local support table and do not re-infer trees.\n\n"
        "## Stage boundary\n\n"
        "Stage 4B stops after spatial-null inference. No ASTRAL/ASTER run was performed, no MSRC model was fit, and no supergene boundary was optimized.\n"
    )


def update_readme() -> None:
    text = README.read_text()
    old = (
        "Stage 0 design freeze is complete.\n\n"
        "Stage 1 public local-tree/TWISST dataset retrieval and normalization is complete. It creates canonical topology-neutral tables for the 213 published windows. It does not classify topologies into the focal quartet, calculate support statistics, compare regions, run ASTRAL/ASTER, fit MSRC parameters, or draw biological conclusions.\n\n"
        "Stage 2 biological group/state provenance freeze is complete. It records the upstream `Species` + `Supergene.Variant` rule used to construct the seven TWISST groups, freezes the focal species and haplotype partitions, and audits sample labels without using topology support.\n\n"
        "Stage 3 independent background species-history freeze is complete. It records the published chromosome 1-15 ASTRAL species-history baseline, maps it to the frozen Stage-2 focal groups, and confirms the resulting background split matches the Stage-0 `species_split`.\n\n"
        "Stage 4A first formal unblinded local quartet-support analysis is complete. It classifies the 945 published TWISST group topologies by the frozen focal quartet, aggregates TWISST weights for the 213 published windows, calculates `q_S`, `q_H`, `q_3`, and `D`, and writes descriptive summaries and raw support figures. It does not perform spatial-null inference, calculate P-values, smooth support tracks, optimize boundaries, run ASTRAL/ASTER, or draw causal conclusions.\n"
    )
    new = (
        "Stage 0 complete.\n\n"
        "Stage 1 complete.\n\n"
        "Stage 2 complete.\n\n"
        "Stage 3 complete.\n\n"
        "Stage 4A complete.\n\n"
        "Stage 4B complete — spatial-null inference. It implements the frozen exact circular-shift null and the frozen physical-coordinate-aware same-width interval sensitivity. It does not optimize supergene boundaries, run ASTRAL/ASTER, fit MSRC parameters, or identify a historical mechanism.\n"
    )
    if old in text:
        text = text.replace(old, new)
    elif new not in text:
        raise Stage4BError("README stage-status block was not found in the expected frozen form")
    old_plan = (
        "Stage 5: Robustness/final manuscript figure and analysis freeze.\n\n"
        "Stage 6: ASTER/ASTRAL4 inference-sensitivity analysis: background versus supergene versus all local windows and progressive supergene downweighting.\n"
    )
    new_plan = (
        "Stage 5: Robustness/manuscript freeze.\n\n"
        "Stage 6: ASTER/ASTRAL4 inference-sensitivity analysis, run separately after Stage 5: background versus supergene versus all local windows and progressive supergene downweighting.\n"
    )
    if old_plan in text:
        text = text.replace(old_plan, new_plan)
    elif new_plan not in text:
        raise Stage4BError("README planned-stage block was not found in the expected frozen form")
    README.write_text(text)


def write_manifest(stage1: dict[str, object], stage2: dict[str, object], stage3: dict[str, object], stage4a: dict[str, object], results: AnalysisResults) -> None:
    result_outputs = [
        CHR16_PHYSICAL_ORDER,
        CIRCULAR_SHIFTS,
        PRIMARY_TEST,
        COORDINATE_PLACEMENTS,
        COORDINATE_SUMMARY,
        COORDINATE_LENGTH_WEIGHTED,
        COORDINATE_LENGTH_WEIGHTED_SUMMARY,
        BOUNDARY_CONTEXT,
        MAIN_TABLE,
        REPORT,
        README,
    ]
    figure_outputs = [
        CIRCULAR_NULL_PDF,
        CIRCULAR_NULL_PNG,
        COORDINATE_NULL_PDF,
        COORDINATE_NULL_PNG,
        COORDINATE_LENGTH_WEIGHTED_PDF,
        COORDINATE_LENGTH_WEIGHTED_PNG,
        SPATIAL_TEST_PDF,
        SPATIAL_TEST_PNG,
    ]
    data = {
        "stage": "Stage 4B",
        "description": "Predeclared spatial-null inference for frozen fire-ant chromosome-16 Stage-4A support track.",
        "parent_checksums": {
            "stage0_focal_quartet": sha256(STAGE0_FOCAL_QUARTET),
            "stage1_manifest": sha256(STAGE1_MANIFEST),
            "stage2_manifest": sha256(STAGE2_MANIFEST),
            "stage2_focal_partition": sha256(STAGE2_FOCAL_PARTITION),
            "stage3_manifest": sha256(STAGE3_MANIFEST),
            "stage3_background_partition": sha256(STAGE3_BACKGROUND_PARTITION),
            "stage4a_manifest": sha256(STAGE4A_MANIFEST),
            "stage4a_window_support": sha256(STAGE4A_WINDOW_SUPPORT),
        },
        "stage1_recorded_checksums": {
            "raw_topology_file_checksums": stage1.get("raw_topology_file_checksums", {}),
            "processed_output_checksums": stage1.get("processed_output_checksums", {}),
        },
        "stage2_recorded_checksums": stage2.get("stage2_output_checksums", {}),
        "stage3_recorded_checksums": {
            "stage3_background_partition": stage3.get("stage3_background_partition_checksum"),
            "background_tree_provenance": stage3.get("background_tree_provenance_checksum"),
            "normalized_background_tree": stage3.get("normalized_background_tree_checksum"),
        },
        "stage4a_output_checksums": stage4a.get("stage4a_output_checksums", {}),
        "stage4a_figure_checksums": stage4a.get("figure_checksums", {}),
        "observed_delta_D": fmt_decimal(results.observed["delta_D"]),
        "physical_chr16_ordering_rule": "chrom == chr16; sort ascending by physical midpoint; do not use Stage-1 concatenation order",
        "physical_chr16_region_sequence": ["chr16A", "chr16_supergene", "chr16B"],
        "n_chr16_windows": len(results.chr16_rows),
        "n_circular_shifts": results.primary_row["n_exact_alignments"],
        "exact_circular_p_value": results.primary_row["p_one_sided"],
        "observed_circular_rank": results.primary_row["observed_rank"],
        "coordinate_width_source": rel(REGION_MANIFEST),
        "physical_interval_width": results.coordinate_summary["interval_width"],
        "coordinate_unique_membership_null": True,
        "coordinate_uniform_start_null": True,
        "coordinate_null_enumeration_rule": results.coordinate_summary["enumeration_rule"],
        "number_unique_coordinate_placements": results.coordinate_summary["n_unique_membership_states"],
        "unique_membership_coordinate_p_value": results.coordinate_summary["p_unique_membership"],
        "length_weighted_coordinate_p_value": results.coordinate_length_weighted_summary["p_length_weighted"],
        "coordinate_start_domain_min": results.coordinate_length_weighted_summary["start_domain_min"],
        "coordinate_start_domain_max": results.coordinate_length_weighted_summary["start_domain_max"],
        "coordinate_total_start_domain_length": results.coordinate_length_weighted_summary["total_start_domain_length"],
        "coordinate_extreme_start_domain_length": results.coordinate_length_weighted_summary["extreme_start_domain_length"],
        "output_checksums": {rel(path): sha256(path) for path in result_outputs},
        "figure_checksums": {rel(path): sha256(path) for path in figure_outputs},
        "processing_script": {
            "path": rel(Path(__file__)),
            "sha256": sha256(Path(__file__)),
        },
        "comparison_tolerance": str(TIE_TOLERANCE),
        "astER_or_astral_run": False,
        "msrc_model_fit": False,
        "boundaries_optimized": False,
        "arbitrary_coordinate_grid_used": False,
        "inferential_tests_performed": [
            "exact_circular_shift_null_primary",
            "unique_membership_state_coordinate_sensitivity",
            "continuous_uniform_start_coordinate_sensitivity",
        ],
    }
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def build_results() -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object], AnalysisResults]:
    stage1, stage2, stage3, stage4a = verify_prior_stages()
    window_rows = load_window_support()
    region_spans = load_region_spans()
    observed = verify_stage4a_contrast(window_rows)
    chr16_rows = physical_chr16_rows(window_rows)
    circular_rows, primary_row = calculate_circular_null(chr16_rows, observed["delta_D"])
    coordinate_rows, coordinate_summary = coordinate_null(chr16_rows, region_spans, observed["delta_D"])
    coordinate_length_weighted_rows, coordinate_length_weighted_summary = coordinate_length_weighted_null(chr16_rows, region_spans, observed["delta_D"])
    boundary_rows = boundary_context(chr16_rows, region_spans)
    results = AnalysisResults(
        window_rows=window_rows,
        chr16_rows=chr16_rows,
        region_spans=region_spans,
        observed=observed,
        circular_rows=circular_rows,
        primary_row=primary_row,
        coordinate_rows=coordinate_rows,
        coordinate_summary=coordinate_summary,
        coordinate_length_weighted_rows=coordinate_length_weighted_rows,
        coordinate_length_weighted_summary=coordinate_length_weighted_summary,
        boundary_rows=boundary_rows,
    )
    return stage1, stage2, stage3, stage4a, results


def run_stage4b(update_readme_file: bool = True) -> dict[str, object]:
    stage1, stage2, stage3, stage4a, results = build_results()
    write_primary_outputs(results)
    make_figures(results)
    write_report(results)
    if update_readme_file:
        update_readme()
    write_manifest(stage1, stage2, stage3, stage4a, results)
    return {
        "observed_delta_D": fmt_decimal(results.observed["delta_D"]),
        "circular_rank": results.primary_row["observed_rank"],
        "circular_p": results.primary_row["p_one_sided"],
        "coordinate_unique_rank": results.coordinate_summary["rank_unique_membership"],
        "coordinate_unique_p": results.coordinate_summary["p_unique_membership"],
        "unique_coordinate_placements": results.coordinate_summary["n_unique_membership_states"],
        "coordinate_length_weighted_p": results.coordinate_length_weighted_summary["p_length_weighted"],
        "manifest_sha": sha256(MANIFEST),
    }


class Stage4BTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.stage1, cls.stage2, cls.stage3, cls.stage4a, cls.results = build_results()

    def test_prior_stage_checksums_pass(self) -> None:
        self.assertEqual(sha256(STAGE0_FOCAL_QUARTET), EXPECTED_STAGE0_FOCAL_SHA)
        self.assertEqual(sha256(STAGE2_FOCAL_PARTITION), EXPECTED_STAGE2_FOCAL_PARTITION_SHA)
        self.assertEqual(sha256(STAGE3_BACKGROUND_PARTITION), EXPECTED_STAGE3_BACKGROUND_PARTITION_SHA)
        self.assertTrue(self.stage3["stage3_background_partition_matches_stage0_species_split"])

    def test_observed_stage4a_delta_reproduces(self) -> None:
        self.assertEqual(self.results.observed["delta_D"].quantize(DECIMAL_PLACES), EXPECTED_DELTA_D)

    def test_chr16_physical_order(self) -> None:
        self.assertEqual(len(self.results.chr16_rows), 96)
        self.assertEqual(Counter(row["region"] for row in self.results.chr16_rows), Counter({"chr16A": 42, "chr16_supergene": 52, "chr16B": 2}))
        self.assertTrue(all(self.results.chr16_rows[i]["mid"] <= self.results.chr16_rows[i + 1]["mid"] for i in range(95)))
        compressed = []
        for row in self.results.chr16_rows:
            if not compressed or compressed[-1] != row["region"]:
                compressed.append(row["region"])
        self.assertEqual(compressed, ["chr16A", "chr16_supergene", "chr16B"])

    def test_original_stage1_row_order_not_used(self) -> None:
        original_chr16 = [row for row in self.results.window_rows if row["chrom"] == "chr16"]
        self.assertNotEqual([row["window_index"] for row in original_chr16], [row["window_index"] for row in self.results.chr16_rows])

    def test_circular_null_properties(self) -> None:
        self.assertEqual(len(self.results.circular_rows), 96)
        self.assertEqual(Decimal(self.results.circular_rows[0]["delta_D"]), EXPECTED_DELTA_D)
        self.assertEqual(self.results.primary_row["observed_rank"], 2)
        self.assertEqual(self.results.primary_row["p_one_sided"], "0.0208333333")
        mask = [row["region"] == "chr16_supergene" for row in self.results.chr16_rows]
        self.assertEqual(sum(mask), 52)
        self.assertEqual(len(mask) - sum(mask), 44)
        original = sorted(row["D"] for row in self.results.chr16_rows)
        shift_96 = [self.results.chr16_rows[(i - 96) % 96]["D"] for i in range(96)]
        self.assertEqual(sorted(shift_96), original)
        shift_0 = [self.results.chr16_rows[i]["D"] for i in range(96)]
        self.assertEqual(shift_96, shift_0)
        self.assertEqual(self.results.primary_row["n_exact_alignments"], 96)

    def test_coordinate_sensitivity_properties(self) -> None:
        self.assertEqual(self.results.coordinate_summary["interval_width"], "16237060")
        self.assertEqual(self.results.coordinate_summary["n_unique_membership_states"], 97)
        self.assertEqual(self.results.coordinate_summary["rank_unique_membership"], 3)
        self.assertEqual(self.results.coordinate_summary["p_unique_membership"], "0.0309278351")
        observed = [row for row in self.results.coordinate_rows if row["is_observed_interval"]]
        self.assertEqual(len(observed), 1)
        indices = set(int(value) for value in observed[0]["inside_window_indices"].split(";"))
        inside = [row for row in self.results.chr16_rows if row["window_index"] in indices]
        self.assertEqual(Counter(row["region"] for row in inside), Counter({"chr16_supergene": 52}))
        rows2, summary2 = coordinate_null(self.results.chr16_rows, self.results.region_spans, self.results.observed["delta_D"])
        self.assertEqual(rows2, self.results.coordinate_rows)
        self.assertEqual(summary2, self.results.coordinate_summary)
        membership_sets = [row["inside_window_indices"] for row in self.results.coordinate_rows]
        self.assertEqual(len(membership_sets), len(set(membership_sets)))
        self.assertIn("s=x_i", self.results.coordinate_summary["enumeration_rule"])
        self.assertFalse(self.results.coordinate_summary["arbitrary_grid_spacing_used"])

    def test_coordinate_length_weighted_properties(self) -> None:
        domain = coordinate_domain_events(self.results.chr16_rows, self.results.region_spans)
        self.assertEqual(domain["width"], Decimal("16237060"))
        self.assertEqual(domain["min_start"], Decimal("23916"))
        self.assertEqual(domain["max_start"], Decimal("12706521"))
        expected_events = sorted(
            set(
                [domain["min_start"], domain["max_start"]]
                + [
                    event
                    for row in self.results.chr16_rows
                    for event in (row["mid"], row["mid"] - domain["width"])
                    if domain["min_start"] <= event <= domain["max_start"]
                ]
            )
        )
        self.assertEqual(domain["breakpoints"], expected_events)
        segment_lengths = [Decimal(row["segment_length"]) for row in self.results.coordinate_length_weighted_rows]
        self.assertEqual(sum(segment_lengths, Decimal("0")), domain["max_start"] - domain["min_start"])
        self.assertEqual(self.results.coordinate_length_weighted_summary["total_start_domain_length"], "12682605")
        self.assertEqual(self.results.coordinate_length_weighted_summary["extreme_start_domain_length"], "993416")
        self.assertEqual(self.results.coordinate_length_weighted_summary["p_length_weighted"], "0.0783290184")
        self.assertAlmostEqual(
            float(Decimal(self.results.coordinate_length_weighted_summary["p_length_weighted"])),
            float(EXPECTED_LENGTH_WEIGHTED_P),
            places=10,
        )
        extreme_segments = [row for row in self.results.coordinate_length_weighted_rows if row["exceeds_or_equals_observed"]]
        self.assertEqual([(row["start_min"], row["start_max"]) for row in extreme_segments], [("11667901", "11700079"), ("11700079", "12432739"), ("12432739", "12661317")])
        for row in self.results.coordinate_length_weighted_rows:
            left = Decimal(row["start_min"])
            right = Decimal(row["start_max"])
            midpoint_key = membership_key(self.results.chr16_rows, (left + right) / Decimal("2"), domain["width"])
            first_quarter_key = membership_key(self.results.chr16_rows, left + (right - left) / Decimal("4"), domain["width"])
            third_quarter_key = membership_key(self.results.chr16_rows, left + (right - left) * Decimal("3") / Decimal("4"), domain["width"])
            self.assertEqual(midpoint_key, first_quarter_key)
            self.assertEqual(midpoint_key, third_quarter_key)
        length_weighted_numerator = sum(Decimal(row["segment_length"]) for row in self.results.coordinate_length_weighted_rows if row["exceeds_or_equals_observed"])
        count_weighted_numerator = sum(1 for row in self.results.coordinate_length_weighted_rows if row["exceeds_or_equals_observed"])
        self.assertEqual(length_weighted_numerator, Decimal("993416"))
        self.assertNotEqual(
            Decimal(count_weighted_numerator) / Decimal(len(self.results.coordinate_length_weighted_rows)),
            Decimal(self.results.coordinate_length_weighted_summary["p_length_weighted"]),
        )
        rows2, summary2 = coordinate_length_weighted_null(self.results.chr16_rows, self.results.region_spans, self.results.observed["delta_D"])
        self.assertEqual(rows2, self.results.coordinate_length_weighted_rows)
        self.assertEqual(summary2, self.results.coordinate_length_weighted_summary)

    def test_no_forbidden_stage4b_actions(self) -> None:
        script = Path(__file__).read_text()
        lowered = script.lower()
        import_lines = [line.strip().lower() for line in script.splitlines() if line.strip().startswith(("import ", "from "))]
        self.assertNotIn("import " + "random", import_lines)
        self.assertNotIn("from " + "random", import_lines)
        self.assertNotIn("numpy" + "." + "random", lowered)
        self.assertNotIn("import " + "subprocess", import_lines)
        self.assertNotIn("subprocess" + ".", lowered)
        if MANIFEST.exists():
            manifest = json.loads(MANIFEST.read_text())
            self.assertFalse(manifest["msrc_model_fit"])
            self.assertFalse(manifest["astER_or_astral_run"])
            self.assertFalse(manifest["boundaries_optimized"])
        self.assertIn("no astral/aster run", lowered)
        function_lines = [line.strip().lower() for line in script.splitlines() if line.strip().startswith("def ")]
        self.assertFalse(any(line.startswith("def fit") for line in function_lines))
        self.assertFalse(any(line.startswith("def optimize") for line in function_lines))

    def test_repeated_outputs_byte_stable_where_practical(self) -> None:
        run_stage4b(update_readme_file=False)
        outputs = [
            CHR16_PHYSICAL_ORDER,
            CIRCULAR_SHIFTS,
            PRIMARY_TEST,
            COORDINATE_PLACEMENTS,
            COORDINATE_SUMMARY,
            COORDINATE_LENGTH_WEIGHTED,
            COORDINATE_LENGTH_WEIGHTED_SUMMARY,
            BOUNDARY_CONTEXT,
            MAIN_TABLE,
            REPORT,
            CIRCULAR_NULL_PDF,
            CIRCULAR_NULL_PNG,
            COORDINATE_NULL_PDF,
            COORDINATE_NULL_PNG,
            COORDINATE_LENGTH_WEIGHTED_PDF,
            COORDINATE_LENGTH_WEIGHTED_PNG,
            SPATIAL_TEST_PDF,
            SPATIAL_TEST_PNG,
            MANIFEST,
        ]
        first = {path: sha256(path) for path in outputs}
        run_stage4b(update_readme_file=False)
        second = {path: sha256(path) for path in outputs}
        self.assertEqual(first, second)


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage4BTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run Stage-4B tests after regenerating outputs")
    args = parser.parse_args(argv)

    summary = run_stage4b(update_readme_file=True)
    if args.run_tests and not run_tests():
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
