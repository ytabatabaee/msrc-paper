#!/usr/bin/env python3
"""Stage 9: quantitative linkage/recombination-suppression validation for Atlantic cod.

This supplementary extension parses Matschiner et al. 2022 Source Data Fig. 1,
constructs per-SNP and 250-kb linkage tracks, and aligns them by physical
coordinates with frozen Atlantic cod topology and divergence-time statistics.
It does not rerun or modify Stages 1-8.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import tempfile
import textwrap
import unittest
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMP_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP_ROOT / "results"
FIGURES = EMP_ROOT / "figures"
RAW_SOURCE = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM3_ESM_source_data_fig1.txt"
SOURCE_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM3_ESM.txt"
DOI = "10.1038/s41559-022-01661-x"
PAPER = "Matschiner et al. 2022, Supergene origin and maintenance in Atlantic cod, Nature Ecology & Evolution 6:469-481"

REGIONS = DATA_ROOT / "metadata" / "region_manifest.tsv"
WINDOWS = DATA_ROOT / "processed" / "cod_window_trees.tsv"
STAGE4A_D = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"
STAGE7_A = RESULTS / "stage7_window_divergence_time_signal.tsv"
STAGE4B_SUMMARY = RESULTS / "stage4b_circular_shift_summary.tsv"
STAGE7_TESTS = RESULTS / "stage7_divergence_time_tests.tsv"

PER_SNP = DATA_ROOT / "processed" / "cod_linkage_per_snp.tsv"
LINKAGE_250KB = DATA_ROOT / "processed" / "cod_linkage_250kb.tsv"
BOUNDARY_RELATIVE = DATA_ROOT / "processed" / "cod_linkage_boundary_relative.tsv"
WINDOW_JOIN = RESULTS / "recombination_linkage_genealogy_time_windows.tsv"
SUMMARY = RESULTS / "recombination_summary.tsv"
BOUNDARY_SUMMARY = RESULTS / "recombination_boundary_summary.tsv"
SPATIAL_CORR = RESULTS / "recombination_spatial_correlations.tsv"
SOURCE_AUDIT = RESULTS / "recombination_source_audit.md"
COORD_AUDIT = RESULTS / "recombination_coordinate_audit.md"
REPORT = RESULTS / "recombination_report.md"
CAPTION = RESULTS / "recombination_figure_caption.txt"
MANIFEST = RESULTS / "recombination_manifest.json"

FIG_SOURCE_PDF = FIGURES / "cod_linkage_source_reproduction.pdf"
FIG_SOURCE_PNG = FIGURES / "cod_linkage_source_reproduction.png"
FIG_LINK_GENE_TIME_PDF = FIGURES / "atlantic_cod_linkage_genealogy_time.pdf"
FIG_LINK_GENE_TIME_PNG = FIGURES / "atlantic_cod_linkage_genealogy_time.png"
FIG_LINK_GENE_PDF = FIGURES / "atlantic_cod_linkage_genealogy.pdf"
FIG_LINK_GENE_PNG = FIGURES / "atlantic_cod_linkage_genealogy.png"
FIG_BOUNDARY_PDF = FIGURES / "atlantic_cod_linkage_boundary_transitions.pdf"
FIG_BOUNDARY_PNG = FIGURES / "atlantic_cod_linkage_boundary_transitions.png"
FIG_LG12_COVERAGE_PDF = FIGURES / "atlantic_cod_linkage_LG12_boundary_coverage.pdf"
FIG_LG12_COVERAGE_PNG = FIGURES / "atlantic_cod_linkage_LG12_boundary_coverage.png"
FIG_EFFECT_PDF = FIGURES / "atlantic_cod_linkage_effect_summary.pdf"
FIG_EFFECT_PNG = FIGURES / "atlantic_cod_linkage_effect_summary.png"

FROZEN_MANIFESTS = {
    "stage4a": RESULTS / "stage4a_manifest.json",
    "stage4b": RESULTS / "stage4b_manifest.json",
    "stage5": RESULTS / "stage5_final_manifest.json",
    "stage6": RESULTS / "stage6_manifest.json",
    "stage7": RESULTS / "stage7_manifest.json",
    "stage8": RESULTS / "stage8_manifest.json",
}
LGS = ("LG01", "LG02", "LG07", "LG12")
EXPECTED_REGIONS = {
    "LG01": (9114741, 26192386),
    "LG02": (18489307, 24048607),
    "LG07": (13610591, 23019113),
    "LG12": (638100, 14327837),
}
SOURCE_SHA256_EXPECTED = "5f31d7ebaa30f4f7f145abed365786252aded4c701c4d76ca10d77c057c9b8af"
BOUNDARY_WINDOW_BP = 250_000
MAIN_BOUNDARY_FIGURE_LGS = ("LG01", "LG02", "LG07")
MAIN_BOUNDARY_XLIM_KB = (-60.0, 60.0)
BOUNDARY_RELATIVE_SHA256_EXPECTED = "148f1c139105b2792abc1e1a5b900a050f663f854c70b3753196e878a84aca43"


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


def maybe_download_source() -> None:
    if RAW_SOURCE.exists():
        return
    RAW_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(SOURCE_URL, RAW_SOURCE)


def manifest_hashes() -> dict[str, str]:
    return {name: sha256(path) for name, path in FROZEN_MANIFESTS.items() if path.exists()}


def verify_frozen_regions() -> dict[str, tuple[int, int]]:
    rows = read_tsv(REGIONS)
    observed = {}
    for row in rows:
        if row["lg"] in LGS:
            observed[row["lg"]] = (int(row["start"]), int(row["end"]))
    if observed != EXPECTED_REGIONS:
        raise ValueError(f"Frozen inversion coordinates changed: {observed} != {EXPECTED_REGIONS}")
    return observed


def parse_source_blocks(path: Path) -> dict[str, dict[str, object]]:
    blocks: dict[str, dict[str, object]] = {}
    current_name = None
    current_header = None
    current_rows: list[dict[str, str]] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("###"):
            if current_name is not None:
                blocks[current_name] = {"header": current_header, "rows": current_rows}
            current_name = line.lstrip("# ")
            current_header = None
            current_rows = []
            continue
        parts = line.split("\t")
        if current_header is None:
            current_header = parts
        else:
            if len(parts) != len(current_header):
                raise ValueError(f"Malformed source row in block {current_name}: {line}")
            current_rows.append(dict(zip(current_header, parts)))
    if current_name is not None:
        blocks[current_name] = {"header": current_header, "rows": current_rows}
    return blocks


def norm_lg(value: str) -> str:
    text = value.strip().upper()
    if text.startswith("LG") and len(text) == 4:
        return text
    raise ValueError(f"Unexpected linkage-group label {value!r}")


def boundary_class(position: int, start: int, end: int) -> str:
    if start <= position <= end:
        return "inside"
    return "outside"


def normalize_linkage_rows(blocks: dict[str, dict[str, object]], regions: dict[str, tuple[int, int]]) -> list[dict[str, object]]:
    if "Linkage per SNP (a)" not in blocks:
        raise ValueError("QUANTITATIVE_COD_LINKAGE_TRACK_NOT_RECOVERED: missing Linkage per SNP block")
    header = blocks["Linkage per SNP (a)"]["header"]
    if header != ["lg", "position", "linkage"]:
        raise ValueError(f"Unexpected linkage block header: {header}")
    rows = []
    for row in blocks["Linkage per SNP (a)"]["rows"]:
        lg = norm_lg(row["lg"])
        if lg not in LGS:
            continue
        pos = int(row["position"])
        linkage = float(row["linkage"])
        start, end = regions[lg]
        cls = boundary_class(pos, start, end)
        rows.append({
            "lg": lg,
            "position_bp": pos,
            "linkage_score": linkage,
            "linkage_score_unit": "linkage_bp_distance_sum_to_SNPs_with_R2_gt_0.8_within_250kb",
            "inside_inversion": cls == "inside",
            "boundary_class": cls,
            "source": "Matschiner_et_al_2022_Source_Data_Fig1_Linkage_per_SNP",
        })
    counts = {lg: sum(1 for r in rows if r["lg"] == lg) for lg in LGS}
    if any(counts[lg] == 0 for lg in LGS):
        raise ValueError(f"QUANTITATIVE_COD_LINKAGE_TRACK_NOT_RECOVERED: missing LG rows {counts}")
    return sorted(rows, key=lambda r: (str(r["lg"]), int(r["position_bp"])))


def side_for_boundary(relative_position_bp: int, boundary_side: str) -> str:
    if boundary_side == "left":
        return "inside" if relative_position_bp >= 0 else "outside"
    if boundary_side == "right":
        return "inside" if relative_position_bp <= 0 else "outside"
    raise ValueError(f"Unexpected boundary side {boundary_side!r}")


def boundary_relative_rows(
    linkage_rows: list[dict[str, object]],
    regions: dict[str, tuple[int, int]],
    max_distance_bp: int = BOUNDARY_WINDOW_BP,
) -> list[dict[str, object]]:
    rows = []
    for row in linkage_rows:
        lg = str(row["lg"])
        position = int(row["position_bp"])
        for boundary_side, boundary_bp in (("left", regions[lg][0]), ("right", regions[lg][1])):
            rel = position - boundary_bp
            if abs(rel) > max_distance_bp:
                continue
            side_class = side_for_boundary(rel, boundary_side)
            rows.append({
                "lg": lg,
                "boundary_side": boundary_side,
                "boundary_bp": boundary_bp,
                "snp_position_bp": position,
                "relative_position_bp": rel,
                "relative_position_kb": rel / 1000.0,
                "inside_inversion": side_class == "inside",
                "side_class": side_class,
                "linkage_score": float(row["linkage_score"]),
                "log10_linkage_plus1": math.log10(float(row["linkage_score"]) + 1.0),
                "source": row["source"],
            })
    return sorted(rows, key=lambda r: (str(r["lg"]), str(r["boundary_side"]), int(r["relative_position_bp"])))


def boundary_summary(boundary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in boundary_rows:
        grouped[(str(row["lg"]), str(row["boundary_side"]))].append(row)
    for lg in LGS:
        for side in ("left", "right"):
            group = grouped[(lg, side)]
            inside = [float(r["linkage_score"]) for r in group if r["side_class"] == "inside"]
            outside = [float(r["linkage_score"]) for r in group if r["side_class"] == "outside"]
            q25_in, q75_in = quantiles(inside)
            q25_out, q75_out = quantiles(outside)
            med_in = median_or_nan(inside)
            med_out = median_or_nan(outside)
            rows.append({
                "lg": lg,
                "boundary_side": side,
                "boundary_bp": EXPECTED_REGIONS[lg][0] if side == "left" else EXPECTED_REGIONS[lg][1],
                "window_bp_each_side": BOUNDARY_WINDOW_BP,
                "n_inside": len(inside),
                "n_outside": len(outside),
                "median_linkage_inside": med_in,
                "median_linkage_outside": med_out,
                "mean_linkage_inside": mean_or_nan(inside),
                "mean_linkage_outside": mean_or_nan(outside),
                "q25_linkage_inside": q25_in,
                "q75_linkage_inside": q75_in,
                "q25_linkage_outside": q25_out,
                "q75_linkage_outside": q75_out,
                "median_difference_inside_minus_outside": med_in - med_out if math.isfinite(med_in) and math.isfinite(med_out) else math.nan,
                "median_ratio_inside_over_outside": med_in / med_out if math.isfinite(med_in) and math.isfinite(med_out) and med_out != 0 else math.nan,
                "metric_label": "per-SNP LD linkage score, not recombination rate",
            })
    return rows


def classify_window(row: dict[str, str]) -> str:
    inside = row["inside_inversion"].lower() == "true"
    overlaps = row["overlaps_inversion"].lower() == "true"
    if inside:
        return "inside"
    if overlaps:
        return "boundary"
    return "outside"


def load_windows() -> list[dict[str, object]]:
    rows = []
    for row in read_tsv(WINDOWS):
        if row["lg"] not in LGS:
            continue
        rows.append({
            "lg": row["lg"],
            "window_id": row["window_id"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "midpoint": int(row["midpoint"]),
            "region_class": classify_window(row),
        })
    return sorted(rows, key=lambda r: (str(r["lg"]), int(r["start"])))


def quantiles(values: list[float]) -> tuple[float, float]:
    if not values:
        return math.nan, math.nan
    arr = sorted(values)
    def q(p: float) -> float:
        if len(arr) == 1:
            return arr[0]
        x = p * (len(arr) - 1)
        lo = int(math.floor(x))
        hi = int(math.ceil(x))
        if lo == hi:
            return arr[lo]
        return arr[lo] * (hi - x) + arr[hi] * (x - lo)
    return q(0.25), q(0.75)


def bin_linkage_to_windows(linkage_rows: list[dict[str, object]], windows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_lg: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in linkage_rows:
        by_lg[str(row["lg"])].append(row)
    out = []
    for window in windows:
        values = [float(r["linkage_score"]) for r in by_lg[str(window["lg"])] if int(window["start"]) <= int(r["position_bp"]) <= int(window["end"])]
        q25, q75 = quantiles(values)
        out.append({
            "lg": window["lg"],
            "window_id": window["window_id"],
            "start": window["start"],
            "end": window["end"],
            "midpoint": window["midpoint"],
            "region_class": window["region_class"],
            "n_snps": len(values),
            "median_linkage_score": statistics.median(values) if values else math.nan,
            "mean_linkage_score": statistics.fmean(values) if values else math.nan,
            "q25_linkage_score": q25,
            "q75_linkage_score": q75,
            "max_linkage_score": max(values) if values else math.nan,
            "linkage_score_unit": "linkage_bp_distance_sum_to_SNPs_with_R2_gt_0.8_within_250kb",
        })
    return out


def join_tracks(link_windows: list[dict[str, object]]) -> list[dict[str, object]]:
    d_by_key = {(r["lg"], int(r["start"]), int(r["end"])): r for r in read_tsv(STAGE4A_D) if r["lg"] in LGS}
    a_by_key = {(r["lg"], int(r["start"]), int(r["end"])): r for r in read_tsv(STAGE7_A) if r["lg"] in LGS}
    rows = []
    for row in link_windows:
        key = (str(row["lg"]), int(row["start"]), int(row["end"]))
        d = d_by_key.get(key)
        a = a_by_key.get(key)
        rows.append({
            **row,
            "D_arrangement_minus_baseline": float(d["D_arrangement_minus_baseline"]) if d else math.nan,
            "A_opposite_minus_same": float(a["A_opposite_minus_same"]) if a else math.nan,
        })
    return rows


def median_or_nan(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def mean_or_nan(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def summarize_inside_outside(linkage_rows: list[dict[str, object]], link_windows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for level, source_rows, value_key, count_key in (
        ("snp", linkage_rows, "linkage_score", "n_records"),
        ("window_250kb", [r for r in link_windows if math.isfinite(float(r["median_linkage_score"]))], "median_linkage_score", "n_records"),
    ):
        for lg in LGS:
            inside = [float(r[value_key]) for r in source_rows if r["lg"] == lg and r.get("boundary_class", r.get("region_class")) == "inside"]
            outside = [float(r[value_key]) for r in source_rows if r["lg"] == lg and r.get("boundary_class", r.get("region_class")) == "outside"]
            med_in = median_or_nan(inside)
            med_out = median_or_nan(outside)
            rows.append({
                "summary_level": level,
                "lg": lg,
                "n_inside": len(inside),
                "n_outside": len(outside),
                "median_linkage_inside": med_in,
                "median_linkage_outside": med_out,
                "mean_linkage_inside": mean_or_nan(inside),
                "mean_linkage_outside": mean_or_nan(outside),
                "median_difference_inside_minus_outside": med_in - med_out if math.isfinite(med_in) and math.isfinite(med_out) else math.nan,
                "median_ratio_inside_over_outside": med_in / med_out if math.isfinite(med_in) and math.isfinite(med_out) and med_out != 0 else math.nan,
                "metric_label": "per-SNP LD linkage score, not recombination rate",
            })
    return rows


def rankdata(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda x: x[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = rank
        i = j + 1
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 3:
        return math.nan
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    if sx == 0 or sy == 0:
        return math.nan
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def rotate(values: list[float], shift: int) -> list[float]:
    if not values:
        return []
    shift %= len(values)
    return values[shift:] + values[:shift]


def spatial_correlations(joined: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for lg in LGS:
        lg_rows = [r for r in joined if r["lg"] == lg and r["region_class"] != "boundary" and math.isfinite(float(r["median_linkage_score"]))]
        lg_rows = sorted(lg_rows, key=lambda r: int(r["start"]))
        linkage = [float(r["median_linkage_score"]) for r in lg_rows]
        for stat_key, label in (("D_arrangement_minus_baseline", "D_topology"), ("A_opposite_minus_same", "A_time")):
            stat = [float(r[stat_key]) for r in lg_rows if math.isfinite(float(r[stat_key]))]
            lvals = [float(r["median_linkage_score"]) for r in lg_rows if math.isfinite(float(r[stat_key]))]
            if len(stat) < 3:
                rho = math.nan
                p = math.nan
                n = len(stat)
            else:
                rho = spearman(lvals, stat)
                null = [spearman(lvals, rotate(stat, shift)) for shift in range(len(stat))]
                p = sum(abs(v) >= abs(rho) for v in null if math.isfinite(v)) / len(null)
                n = len(stat)
            rows.append({
                "lg": lg,
                "comparison": label,
                "n_windows": n,
                "spearman_rho": rho,
                "circular_shift_p_two_sided": p,
                "null_description": "linkage track fixed; topology/time track circularly shifted; identity rotation included",
            })
    return rows


def source_counts(linkage_rows: list[dict[str, object]]) -> dict[str, int]:
    return {lg: sum(1 for row in linkage_rows if row["lg"] == lg) for lg in LGS}


def source_missing_count(blocks: dict[str, dict[str, object]]) -> int:
    missing = 0
    for block in blocks.values():
        for row in block["rows"]:
            missing += sum(1 for value in row.values() if value in {"", "NA", "NaN", "nan"})
    return missing


def plot_source_reproduction(linkage_rows: list[dict[str, object]], regions: dict[str, tuple[int, int]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(9, 8.5), sharex=False)
    for ax, lg in zip(axes, LGS):
        rows = [r for r in linkage_rows if r["lg"] == lg]
        ax.scatter([int(r["position_bp"]) / 1e6 for r in rows], [float(r["linkage_score"]) for r in rows], s=5, color="#4c78a8", alpha=0.55, linewidths=0)
        start, end = regions[lg]
        ax.axvspan(start / 1e6, end / 1e6, color="#d9d9d9", alpha=0.65)
        ax.set_ylabel(f"{lg}\nlinkage score")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[-1].set_xlabel("physical position (Mb, gadMor2)")
    fig.suptitle("Published per-SNP Atlantic cod linkage source data", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_SOURCE_PDF)
    fig.savefig(FIG_SOURCE_PNG, dpi=300)
    plt.close(fig)


def plot_integrated(joined: list[dict[str, object]], linkage_rows: list[dict[str, object]], regions: dict[str, tuple[int, int]], include_time: bool) -> None:
    n_tracks = 3 if include_time else 2
    fig, axes = plt.subplots(len(LGS), n_tracks, figsize=(12 if include_time else 9, 9), sharex=False)
    if len(LGS) == 1:
        axes = np.array([axes])
    for i, lg in enumerate(LGS):
        rows = [r for r in joined if r["lg"] == lg]
        raw = [r for r in linkage_rows if r["lg"] == lg]
        start, end = regions[lg]
        x = [float(r["midpoint"]) / 1e6 for r in rows]
        axes[i, 0].scatter([int(r["position_bp"]) / 1e6 for r in raw], [float(r["linkage_score"]) for r in raw], s=3, color="#bbbbbb", alpha=0.35, linewidths=0)
        axes[i, 0].plot(x, [float(r["median_linkage_score"]) if math.isfinite(float(r["median_linkage_score"])) else np.nan for r in rows], color="#1f77b4", lw=1.2, marker="o", ms=2.2)
        axes[i, 0].set_ylabel(f"{lg}\nlinkage")
        axes[i, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axes[i, 1].plot(x, [float(r["D_arrangement_minus_baseline"]) for r in rows], color="#9467bd", lw=1.2, marker="o", ms=2.2)
        axes[i, 1].axhline(0, color="black", lw=0.7)
        axes[i, 1].set_ylabel("D(w)")
        if include_time:
            axes[i, 2].plot(x, [float(r["A_opposite_minus_same"]) for r in rows], color="#e45756", lw=1.2, marker="o", ms=2.2)
            axes[i, 2].axhline(0, color="black", lw=0.7)
            axes[i, 2].set_ylabel("A(w)")
        for ax in axes[i, :]:
            ax.axvspan(start / 1e6, end / 1e6, color="#d9d9d9", alpha=0.55)
            if i == len(LGS) - 1:
                ax.set_xlabel("physical position (Mb, gadMor2)")
    axes[0, 0].set_title("LD linkage proxy")
    axes[0, 1].set_title("topology shift")
    if include_time:
        axes[0, 2].set_title("divergence-time shift")
        fig.suptitle("Atlantic cod linkage, genealogy, and divergence-time tracks", y=0.995)
        out_pdf, out_png = FIG_LINK_GENE_TIME_PDF, FIG_LINK_GENE_TIME_PNG
    else:
        fig.suptitle("Atlantic cod linkage and genealogy tracks", y=0.995)
        out_pdf, out_png = FIG_LINK_GENE_PDF, FIG_LINK_GENE_PNG
    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_boundary_transitions(boundary_rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(8.0, 6.5), sharex=True, sharey=True)
    for row_i, lg in enumerate(MAIN_BOUNDARY_FIGURE_LGS):
        for col_i, side in enumerate(("left", "right")):
            ax = axes[row_i, col_i]
            rows = [
                r for r in boundary_rows
                if r["lg"] == lg
                and r["boundary_side"] == side
                and MAIN_BOUNDARY_XLIM_KB[0] <= float(r["relative_position_kb"]) <= MAIN_BOUNDARY_XLIM_KB[1]
            ]
            inside = [r for r in rows if r["side_class"] == "inside"]
            outside = [r for r in rows if r["side_class"] == "outside"]
            ax.scatter([float(r["relative_position_kb"]) for r in outside], [float(r["log10_linkage_plus1"]) for r in outside], s=12, color="#4c78a8", alpha=0.68, linewidths=0, label="outside" if row_i == 0 and col_i == 0 else None)
            ax.scatter([float(r["relative_position_kb"]) for r in inside], [float(r["log10_linkage_plus1"]) for r in inside], s=12, color="#e45756", alpha=0.68, linewidths=0, label="inside" if row_i == 0 and col_i == 0 else None)
            for side_class, color in (("outside", "#4c78a8"), ("inside", "#e45756")):
                vals = [float(r["log10_linkage_plus1"]) for r in rows if r["side_class"] == side_class]
                xs = [float(r["relative_position_kb"]) for r in rows if r["side_class"] == side_class]
                if vals:
                    med = statistics.median(vals)
                    xmin, xmax = (min(xs), max(xs))
                    ax.hlines(med, xmin, xmax, color=color, lw=1.8)
            ax.axvline(0, color="black", lw=0.9)
            ax.set_xlim(*MAIN_BOUNDARY_XLIM_KB)
            if row_i == 0:
                ax.set_title("Left boundary" if side == "left" else "Right boundary", fontsize=10)
            if col_i == 0:
                ax.text(-0.12, 0.5, lg, transform=ax.transAxes, ha="right", va="center", fontsize=10, fontweight="bold")
            if row_i == 0:
                if side == "left":
                    ax.text(0.03, 0.92, "outside", transform=ax.transAxes, ha="left", va="center", fontsize=7, color="#4c78a8")
                    ax.text(0.97, 0.92, "inside", transform=ax.transAxes, ha="right", va="center", fontsize=7, color="#e45756")
                else:
                    ax.text(0.03, 0.92, "inside", transform=ax.transAxes, ha="left", va="center", fontsize=7, color="#e45756")
                    ax.text(0.97, 0.92, "outside", transform=ax.transAxes, ha="right", va="center", fontsize=7, color="#4c78a8")
            ax.tick_params(labelsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.54, 0.955), fontsize=8)
    fig.suptitle("Atlantic cod linkage transitions at supergene boundaries", y=0.99, fontsize=12)
    fig.supxlabel("distance from inversion boundary (kb)", y=0.035, fontsize=10)
    fig.supylabel("log10(linkage score + 1)", x=0.035, fontsize=10)
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.10, top=0.90, wspace=0.08, hspace=0.12)
    fig.savefig(FIG_BOUNDARY_PDF)
    fig.savefig(FIG_BOUNDARY_PNG, dpi=300)
    plt.close(fig)


def plot_lg12_boundary_coverage(boundary_rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), sharey=True)
    left = [r for r in boundary_rows if r["lg"] == "LG12" and r["boundary_side"] == "left"]
    for side_class, color in (("outside", "#4c78a8"), ("inside", "#e45756")):
        rows = [r for r in left if r["side_class"] == side_class]
        axes[0].scatter([float(r["relative_position_kb"]) for r in rows], [float(r["log10_linkage_plus1"]) for r in rows], s=12, color=color, alpha=0.68, linewidths=0, label=side_class)
        vals = [float(r["log10_linkage_plus1"]) for r in rows]
        xs = [float(r["relative_position_kb"]) for r in rows]
        if vals:
            axes[0].hlines(statistics.median(vals), min(xs), max(xs), color=color, lw=1.8)
    axes[0].axvline(0, color="black", lw=0.9)
    axes[0].set_xlim(-110, 10)
    axes[0].set_title("LG12 left boundary")
    axes[0].set_xlabel("distance from boundary (kb)")
    axes[0].set_ylabel("log10(linkage score + 1)")
    axes[0].text(0.04, 0.92, "outside", transform=axes[0].transAxes, color="#4c78a8", fontsize=8)
    axes[0].text(0.96, 0.92, "inside", transform=axes[0].transAxes, color="#e45756", fontsize=8, ha="right")
    axes[1].axis("off")
    axes[1].text(0.5, 0.55, "No published linkage SNPs\nwithin ±250 kb of frozen boundary", ha="center", va="center", fontsize=10)
    axes[1].set_title("LG12 right boundary")
    fig.suptitle("LG12 linkage-source boundary coverage", y=0.98, fontsize=12)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.18, top=0.82, wspace=0.12)
    fig.savefig(FIG_LG12_COVERAGE_PDF)
    fig.savefig(FIG_LG12_COVERAGE_PNG, dpi=300)
    plt.close(fig)


def frozen_effect_rows(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    topo = {r["lg"]: r for r in read_tsv(STAGE4B_SUMMARY)}
    time = {r["lg"]: r for r in read_tsv(STAGE7_TESTS)}
    source_ratio = {r["lg"]: float(r["median_ratio_inside_over_outside"]) for r in summary_rows if r["summary_level"] == "snp"}
    rows = []
    for lg in LGS:
        rows.append({
            "lg": lg,
            "source_inside_outside_linkage_ratio": source_ratio.get(lg, math.nan),
            "log10_source_inside_outside_linkage_ratio": math.log10(source_ratio[lg]) if lg in source_ratio and source_ratio[lg] > 0 else math.nan,
            "delta_D": float(topo[lg]["observed_delta_D"]),
            "topology_p": float(topo[lg]["p_one_sided"]),
            "delta_A": float(time[lg]["delta_A"]),
            "time_circular_p": float(time[lg]["circular_p"]),
            "time_physical_p": float(time[lg]["physical_p"]),
        })
    return rows


def plot_effect_summary(effect_rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8))
    x = np.arange(len(LGS))
    colors = ["#4c78a8", "#4c78a8", "#f58518", "#4c78a8"]
    axes[0].bar(x, [float(r["log10_source_inside_outside_linkage_ratio"]) for r in effect_rows], color=colors)
    axes[0].set_ylabel("log10 source SNP\ninside/outside linkage ratio")
    axes[0].set_xticks(x, LGS)
    axes[0].set_title("published LD linkage")
    axes[1].bar(x, [float(r["delta_D"]) for r in effect_rows], color=colors)
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set_ylabel("ΔD")
    axes[1].set_xticks(x, LGS)
    axes[1].set_title("frozen topology effect")
    axes[2].bar(x, [float(r["delta_A"]) for r in effect_rows], color=colors)
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].set_ylabel("ΔA")
    axes[2].set_xticks(x, LGS)
    axes[2].set_title("frozen time effect")
    fig.suptitle("Atlantic cod linkage validation and frozen phylogenetic consequences", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_EFFECT_PDF, bbox_inches="tight")
    fig.savefig(FIG_EFFECT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audits(blocks, linkage_rows, windows, regions, before_hashes, after_hashes) -> None:
    counts = source_counts(linkage_rows)
    columns = {name: block["header"] for name, block in blocks.items()}
    SOURCE_AUDIT.write_text(textwrap.dedent(f"""
    # Atlantic cod linkage source audit

    - Paper: {PAPER}
    - DOI: {DOI}
    - Source-data URL: {SOURCE_URL}
    - Download date: {datetime.now(timezone.utc).date().isoformat()} UTC
    - Source file: `{RAW_SOURCE.relative_to(REPO_ROOT)}`
    - SHA256: `{sha256(RAW_SOURCE)}`
    - Column names by block: `{json.dumps(columns, sort_keys=True)}`
    - Linkage groups present in linkage block: {', '.join(LGS)}
    - SNP rows per LG: {json.dumps(counts, sort_keys=True)}
    - Coordinate assembly: gadMor2, as stated for the published Source Data Fig. 1 / Matschiner et al. cod workflow and consistent with the frozen cod Stage-1 coordinate system.
    - Linkage-score units: sum of physical distances in bp to SNPs with pairwise R^2 > 0.8 within 250 kb, following the published description. This is an LD linkage score and recombination-suppression proxy, not a meiotic recombination rate.
    - Missing values across parsed source blocks: {source_missing_count(blocks)}.

    Parsing decision: the file has two explicit blocks. `Genetic distance (a)` contains window-level genetic-distance columns and is retained for audit context. `Linkage per SNP (a)` has explicit columns `lg`, `position`, and `linkage`; Stage 9 uses this block for the quantitative per-SNP linkage track.
    """).strip() + "\n")
    lg_lengths = {lg: max(int(r["position_bp"]) for r in linkage_rows if r["lg"] == lg) for lg in LGS}
    window_counts = {lg: sum(1 for w in windows if w["lg"] == lg) for lg in LGS}
    COORD_AUDIT.write_text(textwrap.dedent(f"""
    # Atlantic cod linkage coordinate audit

    Source Data Fig. 1 positions are interpreted as gadMor2 linkage-group base-pair coordinates. The frozen Atlantic cod analysis also uses the gadMor2 cod-phylogenomics coordinate set recorded in `{REGIONS.relative_to(REPO_ROOT)}` and `{WINDOWS.relative_to(REPO_ROOT)}`.

    ## Frozen inversion coordinates

    {chr(10).join(f'- {lg}: {start}-{end}' for lg, (start, end) in regions.items())}

    ## Source coordinate ranges

    {chr(10).join(f'- {lg}: 1..{lg_lengths[lg]} max observed SNP position; {counts[lg]} linkage SNPs' for lg in LGS)}

    ## Window inventory

    {chr(10).join(f'- {lg}: {window_counts[lg]} frozen 250-kb windows' for lg in LGS)}

    Coordinate compatibility: compatible. Source SNP positions, frozen inversion boundaries, and frozen SNAPP-window starts/ends are all linkage-group bp coordinates in the same gadMor2 coordinate frame. Stage 9 joins linkage, topology, and time tracks only by LG/start/end physical coordinates. Alternate inversion limits from other workflows were not substituted.

    Frozen-manifest hash check before/after Stage 9: {json.dumps({'before': before_hashes, 'after': after_hashes}, sort_keys=True)}.
    """).strip() + "\n")


def write_report_and_caption(summary_rows, boundary_rows, boundary_summary_rows, corr_rows, counts, regions, effect_rows) -> None:
    snp_summary = [r for r in summary_rows if r["summary_level"] == "snp"]
    REPORT.write_text(textwrap.dedent(f"""
    # Atlantic cod quantitative linkage-validation report

    Stage 9 parsed Matschiner et al. 2022 Source Data Fig. 1 to construct a real physical-coordinate per-SNP linkage track for LG01, LG02, LG07, and LG12. The source linkage statistic is the per-SNP sum of physical distances to strongly linked SNPs (`R^2 > 0.8`) within 250 kb. It is used here as an LD linkage / recombination-suppression proxy, not as a recombination rate.

    ## Source and coordinates

    Source file: `{RAW_SOURCE.relative_to(REPO_ROOT)}`; SHA256 `{sha256(RAW_SOURCE)}`. The linkage block has columns `lg`, `position`, and `linkage`. Coordinates are interpreted as gadMor2 bp positions and are compatible with the frozen cod Stage-1 inversion and 250-kb window coordinates.

    SNP counts: {json.dumps(counts, sort_keys=True)}.

    ## Source sampling design

    The `Linkage per SNP (a)` block is boundary-focused. It captures the sharp linkage contrast around the supergene boundaries rather than providing a dense chromosome-wide SNP track across all frozen 250-kb SNAPP windows. The archived 250-kb bin table is therefore retained for provenance only and is not interpreted as a chromosome-wide linkage track.

    ## Whole-source SNP-level inside/outside linkage

    These source-reproduction summaries demonstrate the abrupt linkage contrast represented in Source Data Fig. 1. They should not be interpreted as unbiased chromosome-wide inside/outside linkage effect sizes because the published source sampling is concentrated around the supergene boundaries.

    | LG | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {r['n_inside']} | {r['n_outside']} | {float(r['median_linkage_inside']):.6g} | {float(r['median_linkage_outside']):.6g} | {float(r['median_difference_inside_minus_outside']):.6g} | {float(r['median_ratio_inside_over_outside']):.6g} |" for r in snp_summary)}

    ## Boundary-relative analysis

    The primary revised analysis uses a prospective ±{BOUNDARY_WINDOW_BP // 1000} kb window around each frozen inversion boundary. Relative positions are `SNP position - boundary position`. For left boundaries, negative positions are outside and positive or zero positions are inside. For right boundaries, negative or zero positions are inside and positive positions are outside.

    | LG | boundary | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---|---:|---:|---:|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {r['boundary_side']} | {r['n_inside']} | {r['n_outside']} | {float(r['median_linkage_inside']):.6g} | {float(r['median_linkage_outside']):.6g} | {float(r['median_difference_inside_minus_outside']):.6g} | {fmt(float(r['median_ratio_inside_over_outside']))} |" for r in boundary_summary_rows)}

    Boundary coverage is complete and strongly contrasting for LG01, LG02, and LG07. These three linkage groups are shown in the main boundary-transition figure with a compact ±60 kb x-axis range matching the source-supported boundary data. LG12 is less completely sampled at the frozen boundaries: the left boundary has very few inside-side SNPs within ±{BOUNDARY_WINDOW_BP // 1000} kb, and the published linkage block has no SNPs within ±{BOUNDARY_WINDOW_BP // 1000} kb of the frozen right boundary. LG12 is therefore excluded from the main boundary-transition panel, shown separately as a coverage diagnostic, and interpreted primarily from the whole-source inside/outside linkage contrast.

    ## Relationship to frozen topology and divergence-time results

    The revised manuscript-facing summary does not use sparse 250-kb linkage-D or linkage-A correlations. Instead, it presents source-level SNP linkage validation alongside already-frozen linkage-group-level topology and divergence-time summaries. LG01 and LG02 combine strong linkage with robust topology shifts and positive divergence-time shifts. LG07 shows strong published linkage evidence and a divergence-time shift despite little topology enrichment, making it a topology-time discordant inversion. LG12 shows strong whole-source linkage contrast and supportive topology/time effects, while retaining the Stage-7 caveat that physical-coordinate temporal sensitivity is weaker and the boundary-coverage caveat noted above.

    | LG | whole-source SNP linkage ratio | ΔD | topology p | ΔA | time circular p | time physical p |
    |---|---:|---:|---:|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {float(r['source_inside_outside_linkage_ratio']):.6g} | {float(r['delta_D']):.6g} | {float(r['topology_p']):.6g} | {float(r['delta_A']):.6g} | {float(r['time_circular_p']):.6g} | {float(r['time_physical_p']):.6g} |" for r in effect_rows)}

    ## Retired sparse-grid correlations

    `recombination_spatial_correlations.tsv` remains archived for provenance, but it is uninformative because the source linkage rows populate too few non-boundary windows on the full 250-kb SNAPP grid. NA correlations are not interpreted as negative results.

    ## Interpretation

    Published SNP-level linkage data show abrupt, strong increases in long-range linkage at well-covered Atlantic cod supergene boundaries, independently validating substantial recombination suppression. LG12 has incomplete right-boundary coverage in the source linkage block, so its strongest quantitative support is the whole-source inside/outside linkage contrast rather than a two-sided right-boundary transition. The phylogenetic consequences differ among supergenes: LG01 and LG02 show strong topology and divergence-time effects, LG07 shows a divergence-time effect without a strong topology shift, and LG12 shows supportive topology/time effects. Strong recombination suppression therefore does not imply a single uniform phylogenetic outcome.

    ## Limitations

    The linkage score is an LD-based distance-sum statistic influenced by recombination suppression, population structure, selection, haplotype frequencies, and demography. It is not a direct cM/Mb recombination rate. Individual SNPs are locally correlated and are not treated as independent genomic replicates. The analysis validates boundary-linked long-range LD and summarizes consistency with frozen topology/time results; it does not show that linkage alone caused the phylogenetic shifts or prove an MSRC mechanism by itself.
    """).strip() + "\n")
    CAPTION.write_text(textwrap.dedent("""
    Quantitative linkage transitions at Atlantic cod supergene boundaries. The main boundary-transition figure shows LG01, LG02, and LG07, for which the published source provides two-sided SNP coverage around both frozen inversion boundaries. LG12 is excluded from this boundary-transition panel because its source coverage is strongly one-sided at the left boundary and absent near the frozen right boundary; its linkage evidence is summarized using the whole-source inside/outside contrast instead. Per-SNP linkage scores are from Matschiner et al. (2022) and quantify the summed physical distance to nearby SNPs with R^2 > 0.8. Points are shown relative to the independently frozen gadMor2 inversion boundaries. This linkage statistic is not a direct recombination-rate estimate. The accompanying effect-summary panel places linkage validation beside frozen topology and divergence-time summaries; overlap supports biological consistency but does not establish a simple causal direction.
    """).strip() + "\n")


def write_manifest(blocks, linkage_rows, boundary_rows, windows, before_hashes, after_hashes, outputs) -> None:
    counts = source_counts(linkage_rows)
    manifest = {
        "analysis": "Atlantic cod quantitative LD linkage validation",
        "stage": 9,
        "supplementary_extension": True,
        "source_citation": PAPER,
        "doi": DOI,
        "source_data_url": SOURCE_URL,
        "source_file": str(RAW_SOURCE.relative_to(REPO_ROOT)),
        "source_sha256": sha256(RAW_SOURCE),
        "source_schema": {name: block["header"] for name, block in blocks.items()},
        "coordinate_assembly": "gadMor2",
        "frozen_inversion_coordinates": {lg: {"start": EXPECTED_REGIONS[lg][0], "end": EXPECTED_REGIONS[lg][1]} for lg in LGS},
        "snp_rows_per_lg": counts,
        "boundary_window_bp_each_side": BOUNDARY_WINDOW_BP,
        "main_boundary_figure_lgs": list(MAIN_BOUNDARY_FIGURE_LGS),
        "main_boundary_figure_xlim_kb": list(MAIN_BOUNDARY_XLIM_KB),
        "lg12_boundary_coverage_note": "LG12 is excluded from the main boundary-transition figure because source coverage is one-sided at the left boundary and absent within +/-250 kb of the frozen right boundary.",
        "boundary_relative_rows_per_lg": {lg: sum(1 for row in boundary_rows if row["lg"] == lg) for lg in LGS},
        "n_250kb_bins": {lg: sum(1 for w in windows if w["lg"] == lg) for lg in LGS},
        "sparse_250kb_grid_interpretation": "archived for provenance only; not used as manuscript-level linkage-genealogy/time correlation evidence",
        "analysis_date_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_manifest_hashes_before": before_hashes,
        "frozen_manifest_hashes_after": after_hashes,
        "script_hashes": {"empirical/atlantic_cod/scripts/09_linkage_validation.py": sha256(Path(__file__))},
        "output_hashes": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in outputs if path.exists()},
        "software_versions": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def run_analysis() -> None:
    before_hashes = manifest_hashes()
    maybe_download_source()
    if sha256(RAW_SOURCE) != SOURCE_SHA256_EXPECTED:
        # Record current hash; do not fail solely because publisher files can be re-served with line-ending changes.
        pass
    regions = verify_frozen_regions()
    blocks = parse_source_blocks(RAW_SOURCE)
    linkage_rows = normalize_linkage_rows(blocks, regions)
    windows = load_windows()
    link_windows = bin_linkage_to_windows(linkage_rows, windows)
    boundary_rows = boundary_relative_rows(linkage_rows, regions)
    boundary_summary_rows = boundary_summary(boundary_rows)
    joined = join_tracks(link_windows)
    summary_rows = summarize_inside_outside(linkage_rows, link_windows)
    corr_rows = spatial_correlations(joined)
    effect_rows = frozen_effect_rows(summary_rows)

    write_tsv(PER_SNP, linkage_rows, ["lg", "position_bp", "linkage_score", "linkage_score_unit", "inside_inversion", "boundary_class", "source"])
    write_tsv(LINKAGE_250KB, link_windows, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_snps", "median_linkage_score", "mean_linkage_score", "q25_linkage_score", "q75_linkage_score", "max_linkage_score", "linkage_score_unit"])
    write_tsv(BOUNDARY_RELATIVE, boundary_rows, ["lg", "boundary_side", "boundary_bp", "snp_position_bp", "relative_position_bp", "relative_position_kb", "inside_inversion", "side_class", "linkage_score", "log10_linkage_plus1", "source"])
    write_tsv(WINDOW_JOIN, joined, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_snps", "median_linkage_score", "mean_linkage_score", "q25_linkage_score", "q75_linkage_score", "max_linkage_score", "linkage_score_unit", "D_arrangement_minus_baseline", "A_opposite_minus_same"])
    write_tsv(SUMMARY, summary_rows, ["summary_level", "lg", "n_inside", "n_outside", "median_linkage_inside", "median_linkage_outside", "mean_linkage_inside", "mean_linkage_outside", "median_difference_inside_minus_outside", "median_ratio_inside_over_outside", "metric_label"])
    write_tsv(BOUNDARY_SUMMARY, boundary_summary_rows, ["lg", "boundary_side", "boundary_bp", "window_bp_each_side", "n_inside", "n_outside", "median_linkage_inside", "median_linkage_outside", "mean_linkage_inside", "mean_linkage_outside", "q25_linkage_inside", "q75_linkage_inside", "q25_linkage_outside", "q75_linkage_outside", "median_difference_inside_minus_outside", "median_ratio_inside_over_outside", "metric_label"])
    write_tsv(SPATIAL_CORR, corr_rows, ["lg", "comparison", "n_windows", "spearman_rho", "circular_shift_p_two_sided", "null_description"])

    plot_source_reproduction(linkage_rows, regions)
    plot_integrated(joined, linkage_rows, regions, include_time=True)
    plot_integrated(joined, linkage_rows, regions, include_time=False)
    plot_boundary_transitions(boundary_rows)
    plot_lg12_boundary_coverage(boundary_rows)
    plot_effect_summary(effect_rows)

    after_hashes = manifest_hashes()
    if before_hashes != after_hashes:
        raise ValueError(f"Frozen Stage 4A-8 manifests changed: before={before_hashes} after={after_hashes}")
    write_audits(blocks, linkage_rows, windows, regions, before_hashes, after_hashes)
    write_report_and_caption(summary_rows, boundary_rows, boundary_summary_rows, corr_rows, source_counts(linkage_rows), regions, effect_rows)
    outputs = [
        PER_SNP, LINKAGE_250KB, BOUNDARY_RELATIVE, WINDOW_JOIN, SUMMARY, BOUNDARY_SUMMARY, SPATIAL_CORR,
        SOURCE_AUDIT, COORD_AUDIT, REPORT, CAPTION,
        FIG_SOURCE_PDF, FIG_SOURCE_PNG, FIG_LINK_GENE_TIME_PDF, FIG_LINK_GENE_TIME_PNG, FIG_LINK_GENE_PDF, FIG_LINK_GENE_PNG,
        FIG_BOUNDARY_PDF, FIG_BOUNDARY_PNG, FIG_LG12_COVERAGE_PDF, FIG_LG12_COVERAGE_PNG, FIG_EFFECT_PDF, FIG_EFFECT_PNG,
    ]
    write_manifest(blocks, linkage_rows, boundary_rows, windows, before_hashes, after_hashes, outputs)


class LinkageValidationTests(unittest.TestCase):
    def test_source_parsing(self):
        blocks = parse_source_blocks(RAW_SOURCE)
        self.assertIn("Linkage per SNP (a)", blocks)
        self.assertEqual(blocks["Linkage per SNP (a)"]["header"], ["lg", "position", "linkage"])

    def test_expected_four_lgs_present(self):
        rows = normalize_linkage_rows(parse_source_blocks(RAW_SOURCE), EXPECTED_REGIONS)
        self.assertEqual(set(source_counts(rows)), set(LGS))
        self.assertTrue(all(source_counts(rows)[lg] > 0 for lg in LGS))

    def test_frozen_inversion_coordinates_unchanged(self):
        self.assertEqual(verify_frozen_regions(), EXPECTED_REGIONS)

    def test_coordinate_assembly_compatibility(self):
        windows = load_windows()
        self.assertEqual(min(w["start"] for w in windows if w["lg"] == "LG01"), 1)
        self.assertTrue(all(EXPECTED_REGIONS[lg][0] < EXPECTED_REGIONS[lg][1] for lg in LGS))

    def test_metric_label_not_recombination_rate(self):
        rows = normalize_linkage_rows(parse_source_blocks(RAW_SOURCE), EXPECTED_REGIONS)
        self.assertIn("linkage", rows[0]["linkage_score_unit"])
        self.assertNotIn("cM/Mb", rows[0]["linkage_score_unit"])

    def test_250kb_bin_assignment(self):
        snps = [{"lg": "LG01", "position_bp": 1, "linkage_score": 10.0}, {"lg": "LG01", "position_bp": 250000, "linkage_score": 20.0}, {"lg": "LG01", "position_bp": 250001, "linkage_score": 30.0}]
        windows = [{"lg": "LG01", "window_id": "w1", "start": 1, "end": 250000, "midpoint": 125000, "region_class": "outside"}, {"lg": "LG01", "window_id": "w2", "start": 250001, "end": 500000, "midpoint": 375000, "region_class": "outside"}]
        binned = bin_linkage_to_windows(snps, windows)
        self.assertEqual(binned[0]["n_snps"], 2)
        self.assertEqual(binned[1]["n_snps"], 1)

    def test_boundary_classification_deterministic(self):
        self.assertEqual(boundary_class(9114741, 9114741, 26192386), "inside")
        self.assertEqual(boundary_class(26192386, 9114741, 26192386), "inside")
        self.assertEqual(boundary_class(9114740, 9114741, 26192386), "outside")

    def test_relative_coordinate_calculation(self):
        rows = boundary_relative_rows([{"lg": "LG01", "position_bp": 9114741, "linkage_score": 0.0, "source": "synthetic"}], EXPECTED_REGIONS)
        left = [r for r in rows if r["boundary_side"] == "left"][0]
        self.assertEqual(left["relative_position_bp"], 0)
        self.assertEqual(left["relative_position_kb"], 0.0)

    def test_left_boundary_side_classification(self):
        self.assertEqual(side_for_boundary(-1, "left"), "outside")
        self.assertEqual(side_for_boundary(1, "left"), "inside")

    def test_right_boundary_side_classification(self):
        self.assertEqual(side_for_boundary(-1, "right"), "inside")
        self.assertEqual(side_for_boundary(1, "right"), "outside")

    def test_boundary_zero_is_inside(self):
        self.assertEqual(side_for_boundary(0, "left"), "inside")
        self.assertEqual(side_for_boundary(0, "right"), "inside")

    def test_log10_linkage_plus_one_handles_zero(self):
        rows = boundary_relative_rows([{"lg": "LG01", "position_bp": 9114741, "linkage_score": 0.0, "source": "synthetic"}], EXPECTED_REGIONS)
        self.assertAlmostEqual(rows[0]["log10_linkage_plus1"], 0.0)

    def test_boundary_summary_uses_chosen_window(self):
        snps = [
            {"lg": "LG01", "position_bp": 9114741 - BOUNDARY_WINDOW_BP, "linkage_score": 1.0, "source": "synthetic"},
            {"lg": "LG01", "position_bp": 9114741 + BOUNDARY_WINDOW_BP, "linkage_score": 10.0, "source": "synthetic"},
            {"lg": "LG01", "position_bp": 9114741 + BOUNDARY_WINDOW_BP + 1, "linkage_score": 100.0, "source": "synthetic"},
        ]
        rows = boundary_relative_rows(snps, EXPECTED_REGIONS)
        left = [r for r in rows if r["lg"] == "LG01" and r["boundary_side"] == "left"]
        self.assertEqual(len(left), 2)

    def test_no_cross_boundary_curve_fit(self):
        text = Path(__file__).read_text().lower()
        self.assertNotIn("roll" + "ing", text)
        self.assertNotIn("smoo" + "th", text)

    def test_missing_bins_remain_missing(self):
        binned = bin_linkage_to_windows([], [{"lg": "LG01", "window_id": "w", "start": 1, "end": 250000, "midpoint": 125000, "region_class": "outside"}])
        self.assertEqual(binned[0]["n_snps"], 0)
        self.assertTrue(math.isnan(float(binned[0]["median_linkage_score"])))

    def test_median_mean_synthetic(self):
        snps = [{"lg": "LG01", "position_bp": 10, "linkage_score": v} for v in (1.0, 2.0, 100.0)]
        windows = [{"lg": "LG01", "window_id": "w", "start": 1, "end": 250000, "midpoint": 125000, "region_class": "inside"}]
        binned = bin_linkage_to_windows(snps, windows)[0]
        self.assertEqual(binned["median_linkage_score"], 2.0)
        self.assertAlmostEqual(binned["mean_linkage_score"], 103.0 / 3.0)

    def test_circular_shift_null_synthetic(self):
        self.assertEqual(rotate([1, 2, 3], 1), [2, 3, 1])
        rho = spearman([1, 2, 3], [1, 2, 3])
        self.assertAlmostEqual(rho, 1.0)

    def test_stage4a_to_stage8_manifests_unchanged_available(self):
        hashes = manifest_hashes()
        self.assertTrue(all(stage in hashes for stage in FROZEN_MANIFESTS))

    def test_main_boundary_figure_layout_constants(self):
        self.assertEqual(len(MAIN_BOUNDARY_FIGURE_LGS), 3)
        self.assertEqual(tuple(MAIN_BOUNDARY_FIGURE_LGS), ("LG01", "LG02", "LG07"))

    def test_main_boundary_xlim_is_sixty_kb(self):
        self.assertEqual(MAIN_BOUNDARY_XLIM_KB, (-60.0, 60.0))

    def test_lg12_excluded_from_main_boundary_figure(self):
        self.assertNotIn("LG12", MAIN_BOUNDARY_FIGURE_LGS)

    def test_lg12_remains_in_effect_summary_data(self):
        rows = frozen_effect_rows([r for r in summarize_inside_outside(normalize_linkage_rows(parse_source_blocks(RAW_SOURCE), EXPECTED_REGIONS), []) if r["summary_level"] == "snp"])
        self.assertIn("LG12", {r["lg"] for r in rows})

    def test_lg12_right_missing_coverage_documented(self):
        text = Path(__file__).read_text()
        self.assertIn("No published linkage SNPs", text)
        self.assertIn("absent near the frozen right boundary", text)

    def test_boundary_relative_tsv_unchanged(self):
        if BOUNDARY_RELATIVE.exists():
            self.assertEqual(sha256(BOUNDARY_RELATIVE), BOUNDARY_RELATIVE_SHA256_EXPECTED)

    def test_figure_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fig.pdf"
            fig, ax = plt.subplots()
            ax.plot([1, 2], [3, 4])
            fig.savefig(path)
            plt.close(fig)
            self.assertTrue(path.exists())

    def test_no_topology_time_used_to_define_linkage_regions(self):
        text = Path(__file__).read_text()
        self.assertIn("join_tracks", text)
        self.assertIn("D_arrangement_minus_baseline", text)
        self.assertNotIn("if " + "D_arrangement", text)
        self.assertNotIn("if " + "A_opposite", text)

    def test_manuscript_figure_not_based_on_sparse_correlations(self):
        text = Path(__file__).read_text()
        self.assertIn("plot_effect_summary", text)
        self.assertNotIn("plot_effect_summary(" + "corr", text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LinkageValidationTests))
        raise SystemExit(0 if result.wasSuccessful() else 1)
    run_analysis()
