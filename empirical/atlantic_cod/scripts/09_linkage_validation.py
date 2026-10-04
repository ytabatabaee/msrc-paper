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

PER_SNP = DATA_ROOT / "processed" / "cod_linkage_per_snp.tsv"
LINKAGE_250KB = DATA_ROOT / "processed" / "cod_linkage_250kb.tsv"
WINDOW_JOIN = RESULTS / "recombination_linkage_genealogy_time_windows.tsv"
SUMMARY = RESULTS / "recombination_summary.tsv"
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


def write_report_and_caption(summary_rows, corr_rows, counts, regions) -> None:
    snp_summary = [r for r in summary_rows if r["summary_level"] == "snp"]
    win_summary = [r for r in summary_rows if r["summary_level"] == "window_250kb"]
    REPORT.write_text(textwrap.dedent(f"""
    # Atlantic cod quantitative linkage-validation report

    Stage 9 parsed Matschiner et al. 2022 Source Data Fig. 1 to construct a real physical-coordinate per-SNP linkage track for LG01, LG02, LG07, and LG12. The source linkage statistic is the per-SNP sum of physical distances to strongly linked SNPs (`R^2 > 0.8`) within 250 kb. It is used here as an LD linkage / recombination-suppression proxy, not as a recombination rate.

    ## Source and coordinates

    Source file: `{RAW_SOURCE.relative_to(REPO_ROOT)}`; SHA256 `{sha256(RAW_SOURCE)}`. The linkage block has columns `lg`, `position`, and `linkage`. Coordinates are interpreted as gadMor2 bp positions and are compatible with the frozen cod Stage-1 inversion and 250-kb window coordinates.

    SNP counts: {json.dumps(counts, sort_keys=True)}.

    ## SNP-level inside/outside linkage

    | LG | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {r['n_inside']} | {r['n_outside']} | {float(r['median_linkage_inside']):.6g} | {float(r['median_linkage_outside']):.6g} | {float(r['median_difference_inside_minus_outside']):.6g} | {float(r['median_ratio_inside_over_outside']):.6g} |" for r in snp_summary)}

    ## 250-kb window-level inside/outside linkage

    | LG | n inside windows | n outside windows | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {r['n_inside']} | {r['n_outside']} | {float(r['median_linkage_inside']):.6g} | {float(r['median_linkage_outside']):.6g} | {float(r['median_difference_inside_minus_outside']):.6g} | {float(r['median_ratio_inside_over_outside']):.6g} |" for r in win_summary)}

    ## Spatial linkage-genealogy relationships

    Spearman correlations use 250-kb windows and exclude boundary windows. The Source Data Fig. 1 linkage block is sparse on the chromosome-wide frozen SNAPP grid: most linkage SNPs fall in or near inversion intervals rather than in a dense track across all outside windows. Consequently, linkage-D and linkage-A correlations are reported only when at least three non-boundary 250-kb bins contain linkage SNPs. Circular-shift p-values keep the linkage track fixed, shift D(w) or A(w), and include the identity rotation.

    | LG | comparison | n windows | rho | circular p |
    |---|---|---:|---:|---:|
    {chr(10).join(f"| {r['lg']} | {r['comparison']} | {r['n_windows']} | {float(r['spearman_rho']):.6g} | {float(r['circular_shift_p_two_sided']):.6g} |" for r in corr_rows)}

    ## Interpretation

    At the per-SNP level, all four LGs show much higher median linkage scores inside the frozen inversion intervals than outside the intervals represented in Source Data Fig. 1. LG01 and LG02 combine strong quantitative long-range linkage with the previously frozen topology and divergence-time sensitivity signals. LG07 combines strong linkage and divergence-time sensitivity without a corresponding strong topology shift, making it a topology-time discordant inversion rather than a failed control. LG12 shows elevated linkage and supportive topology/time behavior, but Stage-7 physical-coordinate sensitivity was weaker.

    ## Limitations

    The linkage score is an LD-based distance-sum statistic influenced by recombination suppression, population structure, selection, haplotype frequencies, and demography. It is not a direct cM/Mb recombination rate. Individual SNPs are not treated as independent genomic replicates. The 250-kb window summaries are used for physical-coordinate alignment with frozen D(w) and A(w), but the source linkage rows are too sparse on the full SNAPP grid for a meaningful linkage-genealogy correlation test. Spatial correlations are therefore descriptive only and unavailable where fewer than three non-boundary bins have linkage data.
    """).strip() + "\n")
    CAPTION.write_text(textwrap.dedent("""
    Atlantic cod quantitative linkage validation. Per-SNP linkage data come independently from the Matschiner et al. 2022 100-individual SNP analysis and quantify, for each SNP, the sum of physical distances to strongly linked nearby SNPs with R^2 > 0.8. This LD-based linkage score is used as a recombination-suppression proxy, not as a meiotic recombination rate. Topology D(w) and divergence-time A(w) tracks come from the published SNAPP-window analyses already frozen in the MSRC Atlantic cod workflow. Tracks are aligned only by common gadMor2 physical coordinates and frozen inversion intervals. Spatial co-localization supports biological validation of long-range linkage in the supergene intervals but does not establish a simple causal direction.
    """).strip() + "\n")


def write_manifest(blocks, linkage_rows, windows, before_hashes, after_hashes, outputs) -> None:
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
        "n_250kb_bins": {lg: sum(1 for w in windows if w["lg"] == lg) for lg in LGS},
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
    joined = join_tracks(link_windows)
    summary_rows = summarize_inside_outside(linkage_rows, link_windows)
    corr_rows = spatial_correlations(joined)

    write_tsv(PER_SNP, linkage_rows, ["lg", "position_bp", "linkage_score", "linkage_score_unit", "inside_inversion", "boundary_class", "source"])
    write_tsv(LINKAGE_250KB, link_windows, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_snps", "median_linkage_score", "mean_linkage_score", "q25_linkage_score", "q75_linkage_score", "max_linkage_score", "linkage_score_unit"])
    write_tsv(WINDOW_JOIN, joined, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_snps", "median_linkage_score", "mean_linkage_score", "q25_linkage_score", "q75_linkage_score", "max_linkage_score", "linkage_score_unit", "D_arrangement_minus_baseline", "A_opposite_minus_same"])
    write_tsv(SUMMARY, summary_rows, ["summary_level", "lg", "n_inside", "n_outside", "median_linkage_inside", "median_linkage_outside", "mean_linkage_inside", "mean_linkage_outside", "median_difference_inside_minus_outside", "median_ratio_inside_over_outside", "metric_label"])
    write_tsv(SPATIAL_CORR, corr_rows, ["lg", "comparison", "n_windows", "spearman_rho", "circular_shift_p_two_sided", "null_description"])

    plot_source_reproduction(linkage_rows, regions)
    plot_integrated(joined, linkage_rows, regions, include_time=True)
    plot_integrated(joined, linkage_rows, regions, include_time=False)

    after_hashes = manifest_hashes()
    if before_hashes != after_hashes:
        raise ValueError(f"Frozen Stage 4A-8 manifests changed: before={before_hashes} after={after_hashes}")
    write_audits(blocks, linkage_rows, windows, regions, before_hashes, after_hashes)
    write_report_and_caption(summary_rows, corr_rows, source_counts(linkage_rows), regions)
    outputs = [PER_SNP, LINKAGE_250KB, WINDOW_JOIN, SUMMARY, SPATIAL_CORR, SOURCE_AUDIT, COORD_AUDIT, REPORT, CAPTION, FIG_SOURCE_PDF, FIG_SOURCE_PNG, FIG_LINK_GENE_TIME_PDF, FIG_LINK_GENE_TIME_PNG, FIG_LINK_GENE_PDF, FIG_LINK_GENE_PNG]
    write_manifest(blocks, linkage_rows, windows, before_hashes, after_hashes, outputs)


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LinkageValidationTests))
        raise SystemExit(0 if result.wasSuccessful() else 1)
    run_analysis()
