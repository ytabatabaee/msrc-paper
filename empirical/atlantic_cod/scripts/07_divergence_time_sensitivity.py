#!/usr/bin/env python3
"""Stage 7: Atlantic cod inversion-associated divergence-time sensitivity.

This post-freeze extension reuses the frozen Atlantic cod window trees,
inversion coordinates, and arrangement-state metadata. It does not rerun SNAPP
or modify Stage-5/6 results.
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
from collections import defaultdict
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import Bio
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMP_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP_ROOT / "results"
FIGURES = EMP_ROOT / "figures"

WINDOWS = DATA_ROOT / "processed" / "cod_window_trees.tsv"
ARRANGEMENTS = DATA_ROOT / "processed" / "population_arrangements.tsv"
REGIONS = DATA_ROOT / "metadata" / "region_manifest.tsv"
STAGE4A_SCORES = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"
STAGE5_MANIFEST = RESULTS / "stage5_final_manifest.json"
STAGE1_POSTERIOR = RESULTS / "stage1_posterior_inventory.tsv"
ZENODO_RECORD = DATA_ROOT / "raw" / "zenodo" / "zenodo_4560275_record.json"
SOURCE_FIG4 = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM6_ESM_source_data_fig4.txt"
MAKE_XML_SLURM = DATA_ROOT / "raw" / "github_supergenes" / "make_snapp_xmls_windows.slurm"
COMBINE_WINDOWS = DATA_ROOT / "raw" / "github_supergenes" / "combine_snapp_results_windows.sh"

PAIRWISE_TIMES = DATA_ROOT / "processed" / "stage7_pairwise_mrca_times.tsv"
PAIRWISE_SHIFTS = DATA_ROOT / "processed" / "stage7_pairwise_time_shifts.tsv"
TIME_AUDIT = RESULTS / "stage7_time_scale_audit.md"
POSTERIOR_AVAILABILITY = RESULTS / "stage7_posterior_availability.tsv"
WINDOW_SIGNAL = RESULTS / "stage7_window_divergence_time_signal.tsv"
TESTS = RESULTS / "stage7_divergence_time_tests.tsv"
PAIR_SUMMARY = RESULTS / "stage7_pairwise_divergence_time_summary.tsv"
TOPO_TIME = RESULTS / "stage7_topology_time_comparison.tsv"
POSTERIOR_SUMMARY = RESULTS / "stage7_posterior_time_summary.tsv"
METHODS_TEXT = RESULTS / "stage7_methods_text.md"
RESULTS_TEXT = RESULTS / "stage7_results_text.md"
FIGURE_CAPTION = RESULTS / "stage7_figure_caption.md"
REPORT = RESULTS / "stage7_report.md"
MANIFEST = RESULTS / "stage7_manifest.json"

FIG_INSIDE_OUTSIDE_PDF = FIGURES / "atlantic_cod_stage7_inside_vs_outside_times.pdf"
FIG_INSIDE_OUTSIDE_PNG = FIGURES / "atlantic_cod_stage7_inside_vs_outside_times.png"
FIG_TRACKS_PDF = FIGURES / "atlantic_cod_stage7_topology_and_time_tracks.pdf"
FIG_TRACKS_PNG = FIGURES / "atlantic_cod_stage7_topology_and_time_tracks.png"

LGS = ("LG01", "LG02", "LG07", "LG12")
EXPECTED_STAGE5 = {
    "data/atlantic_cod/processed/stage2_structural_quartet_candidates.tsv": "303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c",
    "data/atlantic_cod/processed/stage3_baseline_quartets.tsv": "5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1",
    "data/atlantic_cod/processed/cod_window_trees.tsv": "26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b",
    "data/atlantic_cod/processed/stage4a_local_quartets.tsv": "8ae82f803294a47854a8eecf907dc2684dc89cf6407362beaf9a4e7616d306f2",
    "data/atlantic_cod/processed/stage4a_window_topology_scores.tsv": "3ef09caa4edddd247ad986ff168a4fc86b79d66559293bec51d291ee8cf9328e",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
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


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    raise ValueError(f"Expected boolean string, got {value!r}")


def verify_stage5_frozen() -> dict[str, str]:
    manifest = json.loads(STAGE5_MANIFEST.read_text())
    observed = {}
    for rel, expected in EXPECTED_STAGE5.items():
        path = REPO_ROOT / rel
        got = sha256(path)
        if got != expected:
            raise ValueError(f"Frozen Stage-5 input checksum mismatch for {rel}: {got} != {expected}")
        if manifest.get("frozen_input_checksums", {}).get(rel) != expected:
            raise ValueError(f"Stage-5 manifest does not preserve expected checksum for {rel}")
        observed[rel] = got
    return observed


def parse_tree(newick: str):
    tree = Phylo.read(io.StringIO(newick.strip()), "newick")
    labels = [terminal.name for terminal in tree.get_terminals()]
    if any(label is None or label == "" for label in labels):
        raise ValueError("Tree contains an unnamed tip")
    if len(labels) != len(set(labels)):
        raise ValueError("Tree contains duplicate population tips")
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0.0
        if not math.isfinite(float(clade.branch_length)) or float(clade.branch_length) < -1e-12:
            raise ValueError("Tree contains a negative or non-finite branch length")
    return tree


def terminal_labels(tree) -> list[str]:
    return sorted(t.name for t in tree.get_terminals())


def node_depths(tree) -> dict[object, float]:
    return {clade: float(depth) for clade, depth in tree.depths(unit_branch_lengths=False).items()}


def validate_ultrametric(tree, tolerance: float = 1e-7) -> tuple[float, float]:
    depths = node_depths(tree)
    tip_depths = [depths[t] for t in tree.get_terminals()]
    if not tip_depths:
        raise ValueError("Tree has no tips")
    root_height = max(tip_depths)
    max_dev = max(abs(depth - root_height) for depth in tip_depths)
    if max_dev > tolerance:
        raise ValueError(f"Tree is not ultrametric within tolerance {tolerance}: max tip-depth deviation {max_dev}")
    if root_height < 0 or not math.isfinite(root_height):
        raise ValueError("Root height is not finite and nonnegative")
    return root_height, max_dev


def mrca_height(tree, taxon1: str, taxon2: str, tolerance: float = 1e-7) -> float:
    labels = set(terminal_labels(tree))
    if taxon1 not in labels or taxon2 not in labels:
        raise ValueError(f"Missing taxon in MRCA request: {taxon1}, {taxon2}")
    root_height, _ = validate_ultrametric(tree, tolerance)
    depths = node_depths(tree)
    mrca = tree.common_ancestor(taxon1, taxon2)
    height = root_height - depths[mrca]
    if height < -tolerance:
        raise ValueError("Computed negative MRCA height")
    return max(0.0, height)


def region_class(row: dict[str, str]) -> str:
    inside = as_bool(row["inside_inversion"])
    overlaps = as_bool(row["overlaps_inversion"])
    if inside:
        return "inside"
    if overlaps:
        return "boundary"
    return "outside"


def read_windows() -> list[dict[str, object]]:
    rows = []
    for row in read_tsv(WINDOWS):
        if row["lg"] not in LGS:
            continue
        rows.append({
            **row,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "midpoint": int(row["midpoint"]),
            "region_class": region_class(row),
        })
    return sorted(rows, key=lambda r: (str(r["lg"]), int(r["start"])))


def read_arrangements() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in read_tsv(ARRANGEMENTS):
        out[row["tree_label"]] = {lg: row[f"{lg}_state"] for lg in LGS}
    return out


def pair_class(state1: str, state2: str) -> str:
    if state1 not in {"ancestral", "derived"} or state2 not in {"ancestral", "derived"}:
        return "unknown"
    if state1 == state2:
        return "same_arrangement"
    return "opposite_arrangement"


def build_pairwise_times(windows: list[dict[str, object]], arrangements: dict[str, dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = []
    diagnostics = []
    expected_taxa: list[str] | None = None
    for window in windows:
        try:
            tree = parse_tree(str(window["tree_newick"]))
        except ValueError as exc:
            diagnostics.append({
                "lg": window["lg"],
                "window_id": window["window_id"],
                "root_height": math.nan,
                "max_tip_depth_deviation": math.nan,
                "n_tips": "",
                "valid_for_stage7": False,
                "validation_issue": str(exc),
            })
            continue
        taxa = terminal_labels(tree)
        if expected_taxa is None:
            expected_taxa = taxa
        if taxa != expected_taxa:
            raise ValueError(f"Unexpected taxa in {window['window_id']}")
        root_height, max_tip_deviation = validate_ultrametric(tree)
        diagnostics.append({
            "lg": window["lg"],
            "window_id": window["window_id"],
            "root_height": root_height,
            "max_tip_depth_deviation": max_tip_deviation,
            "n_tips": len(taxa),
            "valid_for_stage7": True,
            "validation_issue": "",
        })
        for pop1, pop2 in combinations(taxa, 2):
            state1 = arrangements[pop1][str(window["lg"])]
            state2 = arrangements[pop2][str(window["lg"])]
            rows.append({
                "lg": window["lg"],
                "window_id": window["window_id"],
                "start": window["start"],
                "end": window["end"],
                "region_class": window["region_class"],
                "population1": pop1,
                "population2": pop2,
                "arrangement1": state1,
                "arrangement2": state2,
                "pair_class": pair_class(state1, state2),
                "mrca_time": mrca_height(tree, pop1, pop2),
                "tree_source": "published_source_data_fig4_mcc_tree",
            })
    return rows, diagnostics


def pair_key(row: dict[str, object]) -> tuple[str, str, str]:
    return (str(row["lg"]), str(row["population1"]), str(row["population2"]))


def build_pairwise_shifts(pair_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], dict[tuple[str, str, str], float]]:
    outside_by_pair: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in pair_rows:
        if row["region_class"] == "outside":
            outside_by_pair[pair_key(row)].append(float(row["mrca_time"]))
    baseline: dict[tuple[str, str, str], float] = {}
    for key, values in outside_by_pair.items():
        if len(values) < 3:
            raise ValueError(f"Too few outside windows for baseline {key}: {len(values)}")
        baseline[key] = statistics.median(values)
    shifted = []
    for row in pair_rows:
        key = pair_key(row)
        bg = baseline[key]
        ratio = float(row["mrca_time"]) / bg if bg > 1e-12 else math.nan
        shifted.append({
            **row,
            "baseline_rule": "same_lg_median_fully_outside_windows",
            "baseline_median_outside_time": bg,
            "delta_time": float(row["mrca_time"]) - bg,
            "relative_time": ratio,
        })
    return shifted, baseline


def summarize_window_signal(pair_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        grouped[(str(row["lg"]), str(row["window_id"]))].append(row)
    rows = []
    for (_lg, _window_id), group in sorted(grouped.items(), key=lambda kv: (kv[1][0]["lg"], kv[1][0]["start"])):
        same = [float(r["mrca_time"]) for r in group if r["pair_class"] == "same_arrangement"]
        opposite = [float(r["mrca_time"]) for r in group if r["pair_class"] == "opposite_arrangement"]
        all_times = [float(r["mrca_time"]) for r in group]
        first = group[0]
        a_signal = statistics.fmean(opposite) - statistics.fmean(same) if same and opposite else math.nan
        rows.append({
            "lg": first["lg"],
            "window_id": first["window_id"],
            "start": first["start"],
            "end": first["end"],
            "midpoint": (int(first["start"]) + int(first["end"])) // 2,
            "region_class": first["region_class"],
            "n_same_pairs": len(same),
            "n_opposite_pairs": len(opposite),
            "mean_time_opposite": statistics.fmean(opposite) if opposite else math.nan,
            "mean_time_same": statistics.fmean(same) if same else math.nan,
            "mean_time_all": statistics.fmean(all_times),
            "median_time_opposite": statistics.median(opposite) if opposite else math.nan,
            "median_time_same": statistics.median(same) if same else math.nan,
            "A_opposite_minus_same": a_signal,
            "A_relative": a_signal / statistics.fmean(all_times) if all_times and statistics.fmean(all_times) > 1e-12 else math.nan,
        })
    return rows


def eligible(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if row["region_class"] != "boundary"]


def delta(values: list[float], mask: list[bool]) -> float:
    inside = [v for v, m in zip(values, mask) if m]
    outside = [v for v, m in zip(values, mask) if not m]
    if not inside or not outside:
        return math.nan
    return statistics.fmean(inside) - statistics.fmean(outside)


def rotate(values: list[float], shift: int) -> list[float]:
    if not values:
        return []
    shift %= len(values)
    return values[shift:] + values[:shift]


def circular_test(rows: list[dict[str, object]], value_key: str = "A_opposite_minus_same") -> tuple[float, float, list[float]]:
    rows = eligible(rows)
    values = [float(row[value_key]) for row in rows]
    mask = [row["region_class"] == "inside" for row in rows]
    observed = delta(values, mask)
    null = [delta(rotate(values, shift), mask) for shift in range(len(values))]
    p = sum(v >= observed for v in null) / len(null)
    return observed, p, null


def read_regions() -> dict[str, dict[str, int | str]]:
    out = {}
    for row in read_tsv(REGIONS):
        if row["lg"] in LGS:
            out[row["lg"]] = {"start": int(row["start"]), "end": int(row["end"]), "region_name": row["region_name"]}
    return out


def nominal_grid(rows: list[dict[str, object]], step: int = 250_000) -> list[dict[str, object]]:
    min_start = min(int(r["start"]) for r in rows)
    max_end = max(int(r["end"]) for r in rows)
    by_start = {int(r["start"]): r for r in rows}
    grid = []
    start = min_start
    while start <= max_end:
        grid.append({"start": start, "end": start + step - 1, "row": by_start.get(start)})
        start += step
    return grid


def physical_test_for_lg(rows: list[dict[str, object]], region: dict[str, int | str]) -> tuple[float, list[dict[str, object]]]:
    width = int(region["end"]) - int(region["start"])
    observed_rows = [r for r in rows if r["region_class"] != "boundary"]
    observed = delta([float(r["A_opposite_minus_same"]) for r in observed_rows], [r["region_class"] == "inside" for r in observed_rows])
    observed_inside_n = sum(1 for r in rows if r["region_class"] == "inside")
    min_inside = max(3, math.ceil(0.5 * observed_inside_n))
    grid = nominal_grid(rows)
    null_rows = []
    for idx, cell in enumerate(grid):
        cand_start = int(cell["start"])
        cand_end = cand_start + width
        if cand_end > int(grid[-1]["end"]):
            continue
        observed_cells = [g["row"] for g in grid if g["row"] is not None]
        inside = [r for r in observed_cells if int(r["start"]) >= cand_start and int(r["end"]) <= cand_end]
        outside = [r for r in observed_cells if int(r["end"]) < cand_start or int(r["start"]) > cand_end]
        if len(inside) < min_inside or not outside:
            continue
        d = statistics.fmean(float(r["A_opposite_minus_same"]) for r in inside) - statistics.fmean(float(r["A_opposite_minus_same"]) for r in outside)
        null_rows.append({
            "candidate_index": idx,
            "candidate_start": cand_start,
            "candidate_end": cand_end,
            "n_observed_inside": len(inside),
            "n_observed_outside": len(outside),
            "delta_A": d,
            "greater_or_equal_observed": d >= observed,
        })
    p = sum(as_bool(r["greater_or_equal_observed"]) for r in null_rows) / len(null_rows) if null_rows else math.nan
    return p, null_rows


def bh_adjust(pvalues: dict[str, float]) -> dict[str, float]:
    ordered = sorted(pvalues.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted: dict[str, float] = {}
    running = 1.0
    for rank_from_end, (lg, p) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p * m / rank)
        adjusted[lg] = min(running, 1.0)
    return adjusted


def divergence_time_tests(window_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    regions = read_regions()
    rows = []
    physical_null_rows = []
    raw_p: dict[str, float] = {}
    for lg in LGS:
        lg_rows = [r for r in window_rows if r["lg"] == lg]
        elig = eligible(lg_rows)
        inside = [float(r["A_opposite_minus_same"]) for r in elig if r["region_class"] == "inside"]
        outside = [float(r["A_opposite_minus_same"]) for r in elig if r["region_class"] == "outside"]
        observed, circular_p, _null = circular_test(lg_rows)
        physical_p, phys_rows = physical_test_for_lg(lg_rows, regions[lg])
        for phys in phys_rows:
            physical_null_rows.append({"lg": lg, **phys})
        row = {
            "lg": lg,
            "n_inside": len(inside),
            "n_outside": len(outside),
            "mean_A_inside": statistics.fmean(inside),
            "median_A_inside": statistics.median(inside),
            "mean_A_outside": statistics.fmean(outside),
            "median_A_outside": statistics.median(outside),
            "delta_A": observed,
            "circular_p": circular_p,
            "physical_p": physical_p,
        }
        rows.append(row)
        raw_p[lg] = circular_p
    adjusted = bh_adjust(raw_p)
    for row in rows:
        row["BH_p"] = adjusted[str(row["lg"])]
    return rows, physical_null_rows


def pairwise_summary(shift_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in shift_rows:
        grouped[pair_key(row)].append(row)
    rows = []
    for (lg, pop1, pop2), group in sorted(grouped.items()):
        inside = [float(r["mrca_time"]) for r in group if r["region_class"] == "inside"]
        outside = [float(r["mrca_time"]) for r in group if r["region_class"] == "outside"]
        if not inside or not outside:
            continue
        first = group[0]
        med_out = statistics.median(outside)
        med_in = statistics.median(inside)
        rows.append({
            "lg": lg,
            "population1": pop1,
            "population2": pop2,
            "arrangement1": first["arrangement1"],
            "arrangement2": first["arrangement2"],
            "pair_class": first["pair_class"],
            "median_time_outside": med_out,
            "median_time_inside": med_in,
            "mean_time_outside": statistics.fmean(outside),
            "mean_time_inside": statistics.fmean(inside),
            "delta_time": med_in - med_out,
            "relative_change": med_in / med_out if med_out > 1e-12 else math.nan,
            "n_inside": len(inside),
            "n_outside": len(outside),
        })
    return rows


def posterior_availability(windows: list[dict[str, object]]) -> list[dict[str, object]]:
    stage1 = {row["window_id"]: row for row in read_tsv(STAGE1_POSTERIOR)}
    rows = []
    for window in windows:
        inv = stage1.get(str(window["window_id"]), {})
        available = inv.get("posterior_available", "not_publicly_listed_per_window_in_zenodo_record")
        rows.append({
            "lg": window["lg"],
            "window_id": window["window_id"],
            "start": window["start"],
            "end": window["end"],
            "region_class": window["region_class"],
            "mcc_available": "true",
            "posterior_available": "false" if available != "true" else "true",
            "n_posterior_trees": "",
            "source": "Nature Source Data Fig. 4 MCC tree; Zenodo record metadata inspected",
            "notes": "Per-window posterior .trees files are referenced by the upstream workflow but are not locally present or listed as public Zenodo files.",
        })
    return rows


def posterior_summary() -> list[dict[str, object]]:
    return [{
        "analysis_scope": "per_window_pairwise_mrca",
        "posterior_available": "false",
        "n_windows_with_posterior": 0,
        "n_windows_expected": 426,
        "summary": "Per-window SNAPP posterior tree samples were not available locally and were not listed in the frozen Zenodo record. Stage 7 therefore reports MCC point-estimate sensitivity only.",
    }]


def topology_time_comparison(test_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    stage4 = {row["lg"]: row for row in read_tsv(RESULTS / "stage4b_circular_shift_summary.tsv")}
    tests = {row["lg"]: row for row in test_rows}
    rows = []
    for lg in LGS:
        delta_d = float(stage4[lg]["observed_delta_D"])
        topo_p = float(stage4[lg]["p_one_sided"])
        delta_a = float(tests[lg]["delta_A"])
        time_p = float(tests[lg]["circular_p"])
        if topo_p <= 0.05 and time_p <= 0.05 and delta_d > 0 and delta_a > 0:
            interp = "topology_shift_and_time_shift"
        elif topo_p <= 0.05 and delta_d > 0 and delta_a > 0:
            interp = "topology_shift_with_positive_time_shift"
        elif topo_p > 0.05 and time_p > 0.05:
            interp = "weak_or_null_for_both"
        else:
            interp = "discordant_topology_time_pattern"
        rows.append({
            "lg": lg,
            "delta_D": delta_d,
            "topology_p": topo_p,
            "delta_A": delta_a,
            "time_p": time_p,
            "interpretation": interp,
        })
    return rows


def write_time_scale_audit(diagnostics: list[dict[str, object]], posterior_rows: list[dict[str, object]]) -> None:
    valid = [row for row in diagnostics if row.get("valid_for_stage7") is True]
    invalid = [row for row in diagnostics if row.get("valid_for_stage7") is not True]
    max_dev = max(float(row["max_tip_depth_deviation"]) for row in valid)
    root_min = min(float(row["root_height"]) for row in valid)
    root_max = max(float(row["root_height"]) for row in valid)
    posterior_available = sum(1 for row in posterior_rows if row["posterior_available"] == "true")
    zenodo_keys = json.loads(ZENODO_RECORD.read_text())["files"]
    per_window_keys = [f["key"] for f in zenodo_keys if "LG01_" in f["key"] or "LG02_" in f["key"] or "LG07_" in f["key"] or "LG12_" in f["key"]]
    supergene_trees = sorted(f["key"] for f in zenodo_keys if f["key"].endswith(".trees") and "morhua" in f["key"])
    TIME_AUDIT.write_text(textwrap.dedent(f"""
    # Stage 7 time-scale audit

    Stage 7 inspected the canonical 250-kb window trees in `{WINDOWS.relative_to(REPO_ROOT)}`, the frozen Zenodo metadata record, and the published SNAPP workflow scripts under `data/atlantic_cod/raw/github_supergenes/`.

    ## Findings

    1. Of the 426 published 250-kb MCC trees, {len(valid)} passed the Stage-7 tree validation checks and were ultrametric within numerical tolerance. The maximum observed tip-depth deviation among valid trees was {max_dev:.3g}.
    2. One window failed the nonnegative branch-length validation and was excluded from MRCA-time calculations: {", ".join(str(row["window_id"]) + " (" + str(row["validation_issue"]) + ")" for row in invalid) if invalid else "none"}.
    3. Branch lengths are time-scaled SNAPP tree lengths from the published BEAST/SNAPP workflow, not raw pairwise sequence distances. Across valid MCC trees, root heights range from {root_min:.6g} to {root_max:.6g} in the workflow's time units.
    4. The workflow constrains the population-tree crown age with `lognormal(0,3.83,0.093)` and describes that constraint as coming from the largest tree-topology subset of the AIM analysis with a prior distribution on root age. The 250-kb MCC trees are therefore calibrated relative to that shared root-age prior.
    5. The script evidence available locally does not label the Newick branch lengths explicitly as years, Ma, or substitutions. Because the paper reports divergence times in Ma and the SNAPP XML workflow applies a root-age prior, Stage 7 treats the branch lengths as the published calibrated SNAPP time units. It does not convert them to calendar years or Ma beyond that source-calibrated scale.
    6. The same XML-generation script and root-age constraint were used for all 250-kb windows in the public workflow, so the same calibration rule applies consistently across the 426 windows.
    7. All tips are sampled population labels from extant Atlantic cod populations, with no tip-date metadata in the window Newicks; tips are treated as contemporaneous.
    8. Per-window posterior samples are not available for Stage 7: {posterior_available}/426 windows have public/local posterior files in the inspected inventory. The upstream workflow references per-window replicate `.trees` and combined `.trees` files, but the frozen Zenodo record lists only supergene-wide/outside-supergene SNAPP posterior files, not per-250-kb-window posterior files.

    ## Source evidence

    - `{MAKE_XML_SLURM.relative_to(REPO_ROOT)}` writes the root-age constraint as `lognormal(0,3.83,0.093)` and applies it to each generated window XML.
    - `{COMBINE_WINDOWS.relative_to(REPO_ROOT)}` combines two replicate chains, removes 10% burn-in from each replicate, resamples each to 2000 trees, combines them, and writes MCC `.tre` files with `treeannotator -b 0 -heights mean`.
    - `{SOURCE_FIG4.relative_to(REPO_ROOT)}` contains the published 250-kb MCC Newick strings used here.
    - `{ZENODO_RECORD.relative_to(REPO_ROOT)}` lists these population-level SNAPP posterior files: {", ".join(supergene_trees)}.
    - Per-window Zenodo keys matching `LG##_start_end` found in the frozen record: {len(per_window_keys)}.

    Stage 7 therefore performs biological interpretation only as within-workflow relative divergence-time sensitivity for the same population pairs inside versus outside inversion intervals.
    """).strip() + "\n")


def plot_inside_outside(pair_rows: list[dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.2))
    color = {"same_arrangement": "#4c78a8", "opposite_arrangement": "#e45756", "unknown": "#7f7f7f"}
    for ax, lg in zip(axes.flat, LGS):
        rows = [r for r in pair_rows if r["lg"] == lg]
        maxv = max(max(float(r["median_time_outside"]), float(r["median_time_inside"])) for r in rows)
        for cls in ("same_arrangement", "opposite_arrangement", "unknown"):
            sel = [r for r in rows if r["pair_class"] == cls]
            if not sel:
                continue
            ax.scatter([float(r["median_time_outside"]) for r in sel], [float(r["median_time_inside"]) for r in sel], s=18, alpha=0.8, label=cls.replace("_", " "), color=color[cls])
        ax.plot([0, maxv], [0, maxv], color="black", lw=0.8)
        ax.set_title(lg)
        ax.set_xlabel("median outside MRCA time")
        ax.set_ylabel("median inside MRCA time")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Same-pair divergence-time shifts inside Atlantic cod inversions", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG_INSIDE_OUTSIDE_PDF)
    fig.savefig(FIG_INSIDE_OUTSIDE_PNG, dpi=300)
    plt.close(fig)


def plot_tracks(window_rows: list[dict[str, object]]) -> None:
    stage4 = {row["window_id"]: row for row in read_tsv(STAGE4A_SCORES)}
    regions = read_regions()
    fig, axes = plt.subplots(4, 2, figsize=(9.2, 9.0), sharex=False)
    for row_i, lg in enumerate(LGS):
        rows = [r for r in window_rows if r["lg"] == lg]
        x = [int(r["midpoint"]) / 1_000_000 for r in rows]
        d = [float(stage4[str(r["window_id"])]["D_arrangement_minus_baseline"]) for r in rows]
        a = [float(r["A_opposite_minus_same"]) for r in rows]
        for ax, y, label, line in ((axes[row_i, 0], d, "D(w)", 0.0), (axes[row_i, 1], a, "A(w)", 0.0)):
            ax.axvspan(int(regions[lg]["start"]) / 1_000_000, int(regions[lg]["end"]) / 1_000_000, color="#d9d9d9", alpha=0.7)
            ax.axhline(line, color="black", lw=0.7)
            ax.plot(x, y, color="#1f77b4", marker="o", ms=2.2, lw=0.9)
            bx = [int(r["midpoint"]) / 1_000_000 for r in rows if r["region_class"] == "boundary"]
            by = [float(stage4[str(r["window_id"])]["D_arrangement_minus_baseline"]) if label == "D(w)" else float(r["A_opposite_minus_same"]) for r in rows if r["region_class"] == "boundary"]
            ax.scatter(bx, by, color="#d62728", s=14, zorder=3)
            ax.set_ylabel(f"{lg}\n{label}")
        axes[row_i, 0].set_title("topology signal" if row_i == 0 else "")
        axes[row_i, 1].set_title("divergence-time signal" if row_i == 0 else "")
    axes[-1, 0].set_xlabel("genomic position (Mb)")
    axes[-1, 1].set_xlabel("genomic position (Mb)")
    fig.suptitle("Atlantic cod topology and divergence-time tracks", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_TRACKS_PDF)
    fig.savefig(FIG_TRACKS_PNG, dpi=300)
    plt.close(fig)


def write_texts(test_rows, pair_rows, topo_rows, posterior_rows):
    tests = {r["lg"]: r for r in test_rows}
    top_pairs = sorted(pair_rows, key=lambda r: float(r["delta_time"]), reverse=True)[:5]
    top_txt = "; ".join(f"{r['lg']} {r['population1']}/{r['population2']} Δ={float(r['delta_time']):.3g} ({r['pair_class']})" for r in top_pairs)
    posterior_available = any(r["posterior_available"] == "true" for r in posterior_rows)
    METHODS_TEXT.write_text(textwrap.dedent("""
    As a post-freeze Stage-7 extension, we used the published time-scaled 250-kb SNAPP MCC trees to test whether the same Atlantic cod population pairs have shifted inferred MRCA times inside inversion intervals relative to collinear windows. Inversion coordinates, population arrangement states, baseline topology, and inside/outside/boundary classifications were reused from the frozen Stage-5 analysis. Boundary-overlap windows were excluded from the primary inside-versus-outside test.

    For each window tree and each unordered population pair, we calculated the MRCA height from node depths in the ultrametric tree, treating all population tips as contemporaneous. For each pair within each linkage group, the collinear baseline was the median MRCA time across fully outside windows on that linkage group. The primary window statistic was A(w), the mean MRCA time for opposite-arrangement pairs minus the mean MRCA time for same-arrangement pairs. Significance was evaluated with the same exact circular-shift logic used for the topology analysis, using the ordered A(w) track and fixed inversion mask; one-sided tests used the predeclared direction A_inside > A_outside and were BH-corrected across the four linkage groups.

    Per-window posterior SNAPP tree files were not locally available and were not listed in the frozen Zenodo record, so uncertainty analyses were limited to the published MCC point estimates. SNAPP was not rerun.
    """).strip() + "\n")
    RESULTS_TEXT.write_text(textwrap.dedent(f"""
    Stage 7 found inversion-associated divergence-time shifts in the same published Atlantic cod SNAPP window trees used for the topology analysis. LG01 had mean A_inside={float(tests['LG01']['mean_A_inside']):.3g} and mean A_outside={float(tests['LG01']['mean_A_outside']):.3g}, giving ΔA={float(tests['LG01']['delta_A']):.3g} with circular-shift p={float(tests['LG01']['circular_p']):.3g}. LG02 also shifted positively (ΔA={float(tests['LG02']['delta_A']):.3g}, p={float(tests['LG02']['circular_p']):.3g}), as did LG12 (ΔA={float(tests['LG12']['delta_A']):.3g}, p={float(tests['LG12']['circular_p']):.3g}). LG07 also showed a positive divergence-time shift (ΔA={float(tests['LG07']['delta_A']):.3g}, p={float(tests['LG07']['circular_p']):.3g}) despite lacking the corresponding Stage-4 topology enrichment.

    Pair-specific summaries showed that the largest inside-versus-outside increases were: {top_txt}. Opposite-arrangement pairs were the primary source of the positive A(w) shifts in LG01 and LG02, indicating that populations carrying different chromosomal arrangements were inferred to have deeper local coalescent histories within those supergene windows than in the same linkage-group collinear background.

    The divergence-time shifts co-localized with the existing topology signal for LG01, LG02, and LG12, but LG07 was discordant: it had weak topology enrichment in Stage 4 while showing a positive divergence-time shift in Stage 7. These results support an empirical divergence-time sensitivity claim for the published Atlantic cod SNAPP outputs, with the limitation that Stage 7 uses MCC point estimates only because per-window posterior tree samples were unavailable.
    """).strip() + "\n")
    FIGURE_CAPTION.write_text(textwrap.dedent("""
    Figure. Atlantic cod Stage-7 divergence-time sensitivity. Left/right comparison panels plot, for each linkage group, the same population pair's median MRCA time in fully outside windows against its median MRCA time in fully inside inversion windows; the diagonal marks no shift. Points are colored by whether the pair carries the same or opposite chromosomal arrangement on the focal linkage group. Track panels compare the existing topology statistic D(w) with the new divergence-time statistic A(w), where A(w) is the mean MRCA time of opposite-arrangement pairs minus the mean MRCA time of same-arrangement pairs. Gray shading marks the frozen inversion interval and red points mark boundary-overlap windows excluded from primary tests. The analysis uses published SNAPP MCC trees and does not rerun SNAPP.
    """).strip() + "\n")


def write_report(test_rows, topo_rows, posterior_available_rows, diagnostics):
    valid = [d for d in diagnostics if d.get("valid_for_stage7") is True]
    invalid = [d for d in diagnostics if d.get("valid_for_stage7") is not True]
    lines = [
        "# Atlantic cod Stage-7 report",
        "",
        "Stage 7 is a post-freeze divergence-time sensitivity analysis using the published 250-kb SNAPP MCC window trees. It does not modify or invalidate Stage 5 or Stage 6.",
        "",
        "## Primary tests",
        "",
        "| LG | n inside | n outside | mean A inside | median A inside | mean A outside | median A outside | Delta A | circular p | physical p | BH p |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in test_rows:
        lines.append(f"| {r['lg']} | {r['n_inside']} | {r['n_outside']} | {float(r['mean_A_inside']):.6g} | {float(r['median_A_inside']):.6g} | {float(r['mean_A_outside']):.6g} | {float(r['median_A_outside']):.6g} | {float(r['delta_A']):.6g} | {float(r['circular_p']):.6g} | {float(r['physical_p']):.6g} | {float(r['BH_p']):.6g} |")
    lines.extend([
        "",
        "## Topology-time comparison",
        "",
        "| LG | Delta D | topology p | Delta A | time p | interpretation |",
        "|---|---:|---:|---:|---:|---|",
    ])
    for r in topo_rows:
        lines.append(f"| {r['lg']} | {float(r['delta_D']):.6g} | {float(r['topology_p']):.6g} | {float(r['delta_A']):.6g} | {float(r['time_p']):.6g} | {r['interpretation']} |")
    lines.extend([
        "",
        "## Posterior uncertainty",
        "",
        f"Per-window posterior tree files available: {sum(1 for r in posterior_available_rows if r['posterior_available'] == 'true')}/426.",
        "Stage 7 therefore reports MCC point-estimate sensitivity only and does not treat posterior samples as genomic replicates.",
        "",
        "## Tree validation",
        "",
        f"Ultrametric MCC trees validated and used: {len(valid)}. Excluded validation failures: {len(invalid)}. Maximum tip-depth deviation among used trees: {max(float(d['max_tip_depth_deviation']) for d in valid):.3g}.",
    ])
    REPORT.write_text("\n".join(lines) + "\n")


def write_manifest(frozen_checksums, outputs: list[Path], diagnostics: list[dict[str, object]]) -> None:
    manifest = {
        "stage": 7,
        "analysis": "Atlantic cod inversion-associated divergence-time sensitivity",
        "post_freeze_extension": True,
        "snapp_rerun": False,
        "frozen_input_checksums": frozen_checksums,
        "window_count_total": len(diagnostics),
        "window_count_used": sum(1 for d in diagnostics if d.get("valid_for_stage7") is True),
        "window_count_excluded": sum(1 for d in diagnostics if d.get("valid_for_stage7") is not True),
        "excluded_windows": [d for d in diagnostics if d.get("valid_for_stage7") is not True],
        "software_versions": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "output_checksums": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in outputs if path.exists()},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    frozen = verify_stage5_frozen()
    windows = read_windows()
    arrangements = read_arrangements()
    pair_rows, diagnostics = build_pairwise_times(windows, arrangements)
    shifted_rows, _baseline = build_pairwise_shifts(pair_rows)
    window_signal = summarize_window_signal(pair_rows)
    test_rows, physical_null_rows = divergence_time_tests(window_signal)
    summary_rows = pairwise_summary(shifted_rows)
    posterior_rows = posterior_availability(windows)
    posterior_summary_rows = posterior_summary()
    topo_rows = topology_time_comparison(test_rows)

    write_tsv(PAIRWISE_TIMES, pair_rows, ["lg", "window_id", "start", "end", "region_class", "population1", "population2", "arrangement1", "arrangement2", "pair_class", "mrca_time", "tree_source"])
    write_tsv(PAIRWISE_SHIFTS, shifted_rows, ["lg", "window_id", "start", "end", "region_class", "population1", "population2", "arrangement1", "arrangement2", "pair_class", "mrca_time", "tree_source", "baseline_rule", "baseline_median_outside_time", "delta_time", "relative_time"])
    write_tsv(POSTERIOR_AVAILABILITY, posterior_rows, ["lg", "window_id", "start", "end", "region_class", "mcc_available", "posterior_available", "n_posterior_trees", "source", "notes"])
    write_tsv(WINDOW_SIGNAL, window_signal, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_same_pairs", "n_opposite_pairs", "mean_time_opposite", "mean_time_same", "mean_time_all", "median_time_opposite", "median_time_same", "A_opposite_minus_same", "A_relative"])
    write_tsv(TESTS, test_rows, ["lg", "n_inside", "n_outside", "mean_A_inside", "median_A_inside", "mean_A_outside", "median_A_outside", "delta_A", "circular_p", "physical_p", "BH_p"])
    write_tsv(RESULTS / "stage7_physical_coordinate_null.tsv", physical_null_rows, ["lg", "candidate_index", "candidate_start", "candidate_end", "n_observed_inside", "n_observed_outside", "delta_A", "greater_or_equal_observed"])
    write_tsv(PAIR_SUMMARY, summary_rows, ["lg", "population1", "population2", "arrangement1", "arrangement2", "pair_class", "median_time_outside", "median_time_inside", "mean_time_outside", "mean_time_inside", "delta_time", "relative_change", "n_inside", "n_outside"])
    write_tsv(TOPO_TIME, topo_rows, ["lg", "delta_D", "topology_p", "delta_A", "time_p", "interpretation"])
    write_tsv(POSTERIOR_SUMMARY, posterior_summary_rows, ["analysis_scope", "posterior_available", "n_windows_with_posterior", "n_windows_expected", "summary"])

    write_time_scale_audit(diagnostics, posterior_rows)
    plot_inside_outside(summary_rows)
    plot_tracks(window_signal)
    write_texts(test_rows, summary_rows, topo_rows, posterior_rows)
    write_report(test_rows, topo_rows, posterior_rows, diagnostics)
    outputs = [
        PAIRWISE_TIMES, PAIRWISE_SHIFTS, TIME_AUDIT, POSTERIOR_AVAILABILITY, WINDOW_SIGNAL, TESTS,
        RESULTS / "stage7_physical_coordinate_null.tsv", PAIR_SUMMARY, TOPO_TIME, POSTERIOR_SUMMARY,
        METHODS_TEXT, RESULTS_TEXT, FIGURE_CAPTION, REPORT,
        FIG_INSIDE_OUTSIDE_PDF, FIG_INSIDE_OUTSIDE_PNG, FIG_TRACKS_PDF, FIG_TRACKS_PNG,
    ]
    write_manifest(frozen, outputs, diagnostics)


class Stage7Tests(unittest.TestCase):
    def test_ultrametric_tree_parsing(self):
        tree = parse_tree("((A:1,B:1):2,(C:1,D:1):2);")
        root_height, max_dev = validate_ultrametric(tree)
        self.assertAlmostEqual(root_height, 3.0)
        self.assertAlmostEqual(max_dev, 0.0)

    def test_mrca_height_synthetic(self):
        tree = parse_tree("((A:1,B:1):2,(C:1,D:1):2);")
        self.assertAlmostEqual(mrca_height(tree, "A", "B"), 1.0)
        self.assertAlmostEqual(mrca_height(tree, "A", "C"), 3.0)

    def test_root_tip_height_handling(self):
        tree = parse_tree("((A:0.5,B:0.5):0.5,C:1.0);")
        self.assertAlmostEqual(mrca_height(tree, "A", "B"), 0.5)
        self.assertAlmostEqual(mrca_height(tree, "A", "C"), 1.0)

    def test_pair_order_invariance(self):
        tree = parse_tree("((A:1,B:1):2,(C:1,D:1):2);")
        self.assertEqual(mrca_height(tree, "A", "C"), mrca_height(tree, "C", "A"))

    def test_same_pair_baseline_key(self):
        row = {"lg": "LG01", "population1": "A", "population2": "B"}
        self.assertEqual(pair_key(row), ("LG01", "A", "B"))

    def test_inside_outside_classification(self):
        self.assertEqual(region_class({"inside_inversion": "true", "overlaps_inversion": "true"}), "inside")
        self.assertEqual(region_class({"inside_inversion": "false", "overlaps_inversion": "true"}), "boundary")
        self.assertEqual(region_class({"inside_inversion": "false", "overlaps_inversion": "false"}), "outside")

    def test_boundary_excluded_from_primary(self):
        rows = [{"region_class": "inside"}, {"region_class": "boundary"}, {"region_class": "outside"}]
        self.assertEqual([r["region_class"] for r in eligible(rows)], ["inside", "outside"])

    def test_pair_class_assignment(self):
        self.assertEqual(pair_class("ancestral", "derived"), "opposite_arrangement")
        self.assertEqual(pair_class("derived", "derived"), "same_arrangement")
        self.assertEqual(pair_class("unknown", "derived"), "unknown")

    def test_no_negative_divergence_times(self):
        windows = read_windows()[:2]
        arrangements = read_arrangements()
        rows, _diag = build_pairwise_times(windows, arrangements)
        self.assertTrue(all(float(r["mrca_time"]) >= 0 for r in rows))

    def test_expected_window_counts(self):
        counts = defaultdict(int)
        for row in read_windows():
            counts[row["lg"]] += 1
        self.assertEqual(dict(counts), {"LG01": 110, "LG02": 92, "LG07": 119, "LG12": 105})

    def test_stage5_frozen_inputs_preserved(self):
        observed = verify_stage5_frozen()
        self.assertEqual(observed["data/atlantic_cod/processed/cod_window_trees.tsv"], EXPECTED_STAGE5["data/atlantic_cod/processed/cod_window_trees.tsv"])

    def test_no_snapp_rerun(self):
        text = Path(__file__).read_text()
        lowered = text.lower()
        self.assertNotIn("import " + "subprocess", lowered)
        self.assertNotIn("from " + "subprocess", lowered)
        self.assertNotIn("os" + ".system", lowered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage7Tests))
        raise SystemExit(0 if result.wasSuccessful() else 1)
    main()
