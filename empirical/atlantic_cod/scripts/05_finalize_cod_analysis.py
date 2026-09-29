#!/usr/bin/env python3
"""Finalize the Atlantic cod empirical analysis for manuscript use.

Stage 5 verifies all frozen inputs, independently reproduces Stage-4B headline
statistics from frozen Stage-4A scores, adds a physical-coordinate-aware
sensitivity null, and writes manuscript-ready figures/tables/text. It does not
alter Stage 2-4 frozen definitions or run new exploratory tests.
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
import sys
import tempfile
import textwrap
import unittest
from collections import defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import Bio


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"
RESULTS = EMPIRICAL_ROOT / "results"
FIGURES = EMPIRICAL_ROOT / "figures"

STAGE2 = DATA_ROOT / "processed" / "stage2_structural_quartet_candidates.tsv"
STAGE3 = DATA_ROOT / "processed" / "stage3_baseline_quartets.tsv"
LOCAL_WINDOWS = DATA_ROOT / "processed" / "cod_window_trees.tsv"
STAGE4A_LOCAL = DATA_ROOT / "processed" / "stage4a_local_quartets.tsv"
STAGE4A_SCORES = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"
REGIONS = DATA_ROOT / "metadata" / "region_manifest.tsv"

EXPECTED = {
    STAGE2: "303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c",
    STAGE3: "5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1",
    LOCAL_WINDOWS: "26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b",
    STAGE4A_LOCAL: "8ae82f803294a47854a8eecf907dc2684dc89cf6407362beaf9a4e7616d306f2",
    STAGE4A_SCORES: "3ef09caa4edddd247ad986ff168a4fc86b79d66559293bec51d291ee8cf9328e",
}

STAGE4B_OUTPUTS = [
    RESULTS / "stage4b_circular_shift_summary.tsv",
    RESULTS / "stage4b_circular_shift_null.tsv",
    RESULTS / "stage4b_primary_pvalues.tsv",
    RESULTS / "stage4b_component_tests.tsv",
    RESULTS / "stage4b_block_placement_summary.tsv",
    RESULTS / "stage4b_boundary_sensitivity.tsv",
    RESULTS / "stage4b_leave_one_population_out.tsv",
    RESULTS / "stage4b_topology_run_summary.tsv",
    RESULTS / "stage4b_global_test.tsv",
]

PHYS_NULL = RESULTS / "stage5_physical_coordinate_null.tsv"
PHYS_SUMMARY = RESULTS / "stage5_physical_coordinate_summary.tsv"
ROBUSTNESS = RESULTS / "stage5_robustness_summary.tsv"
MANIFEST = RESULTS / "stage5_final_manifest.json"
FINAL_ANALYSIS = RESULTS / "FINAL_ANALYSIS.md"
MAIN_FIG_PDF = FIGURES / "atlantic_cod_main.pdf"
MAIN_FIG_PNG = FIGURES / "atlantic_cod_main.png"
SUPP_FIG_PDF = FIGURES / "atlantic_cod_supplement.pdf"
SUPP_FIG_PNG = FIGURES / "atlantic_cod_supplement.png"
MAIN_CAPTION = RESULTS / "atlantic_cod_main_figure_caption.md"
MAIN_TABLE_TSV = RESULTS / "atlantic_cod_main_table.tsv"
MAIN_TABLE_MD = RESULTS / "atlantic_cod_main_table.md"
RESULTS_TEXT = RESULTS / "atlantic_cod_results_text.md"
METHODS_TEXT = RESULTS / "atlantic_cod_methods_text.md"
BORNHOLM_NOTE = RESULTS / "atlantic_cod_bornholm_note.md"

LGS = ("LG01", "LG02", "LG07", "LG12")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_frozen() -> dict[str, str]:
    observed = {}
    for path, expected in EXPECTED.items():
        value = sha256(path)
        if value != expected:
            raise ValueError(f"Frozen checksum mismatch for {path}: {value} != {expected}")
        observed[str(path.relative_to(REPO_ROOT))] = value
    return observed


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fields})


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6f}"
    return str(value)


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    raise ValueError(value)


def load_scores() -> list[dict[str, object]]:
    rows = []
    for row in read_tsv(STAGE4A_SCORES):
        converted: dict[str, object] = dict(row)
        for key in ("start", "end", "midpoint", "n_informative_quartets"):
            converted[key] = int(row[key])
        for key in ("fraction_arrangement", "fraction_baseline", "fraction_third", "D_arrangement_minus_baseline"):
            converted[key] = float(row[key])
        rows.append(converted)
    return sorted(rows, key=lambda r: (str(r["lg"]), int(r["midpoint"])))


def is_boundary(row: dict[str, object]) -> bool:
    return as_bool(row["overlaps_inversion"]) and not as_bool(row["inside_inversion"])


def eligible(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if not is_boundary(row)]


def delta(values: list[float], mask: list[bool]) -> float:
    inside = [v for v, m in zip(values, mask) if m]
    outside = [v for v, m in zip(values, mask) if not m]
    return statistics.fmean(inside) - statistics.fmean(outside)


def rotate(values: list[float], shift: int) -> list[float]:
    shift %= len(values)
    return values[shift:] + values[:shift]


def circular(rows: list[dict[str, object]], key: str = "D_arrangement_minus_baseline") -> tuple[float, float, int, list[float]]:
    rows = eligible(rows)
    values = [float(row[key]) for row in rows]
    mask = [as_bool(row["inside_inversion"]) for row in rows]
    obs = delta(values, mask)
    null = [delta(rotate(values, shift), mask) for shift in range(len(values))]
    rank = sum(v > obs for v in null) + 1
    p = sum(v >= obs for v in null) / len(null)
    return obs, p, rank, null


def bh(pvals: dict[str, float]) -> dict[str, float]:
    ordered = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(ordered)
    out = {}
    running = 1.0
    for rev_i, (lg, p) in enumerate(reversed(ordered), start=1):
        rank = m - rev_i + 1
        running = min(running, p * m / rank)
        out[lg] = min(running, 1.0)
    return out


def reproduce_stage4b(scores: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    pvals = {}
    for lg in LGS:
        lg_rows = [r for r in scores if r["lg"] == lg]
        obs, p, rank, null = circular(lg_rows)
        da, p_a, _, _ = circular(lg_rows, "fraction_arrangement")
        db, p_b, _, _ = circular(lg_rows, "fraction_baseline")
        inside = [float(r["D_arrangement_minus_baseline"]) for r in eligible(lg_rows) if as_bool(r["inside_inversion"])]
        outside = [float(r["D_arrangement_minus_baseline"]) for r in eligible(lg_rows) if not as_bool(r["inside_inversion"])]
        rows.append({
            "lg": lg,
            "n_windows_inside": len(inside),
            "n_windows_outside": len(outside),
            "mean_D_inside": statistics.fmean(inside),
            "mean_D_outside": statistics.fmean(outside),
            "delta_D": obs,
            "delta_A": da,
            "delta_B": db,
            "circular_p": p,
            "circular_rank": rank,
            "null_min": min(null),
            "null_max": max(null),
        })
        pvals[lg] = p
    adjusted = bh(pvals)
    for row in rows:
        row["BH_p"] = adjusted[str(row["lg"])]
    expected = {"LG01": (1.232766, 0.009174, 0.021978), "LG02": (0.892332, 0.010989, 0.021978), "LG07": (-0.005981, 0.555556, 0.555556), "LG12": (1.043170, 0.038462, 0.051282)}
    for row in rows:
        exp = expected[str(row["lg"])]
        for got, target in zip((row["delta_D"], row["circular_p"], row["BH_p"]), exp):
            if abs(float(got) - target) > 1e-5:
                raise ValueError(f"Stage-4B reproduction mismatch for {row['lg']}: got {got}, expected {target}")
    return rows


def load_regions() -> dict[str, dict[str, object]]:
    regions = {}
    for row in read_tsv(REGIONS):
        if row["lg"] in LGS:
            regions[row["lg"]] = {"start": int(row["start"]), "end": int(row["end"]), "region_name": row["region_name"], "coordinate_system": row["coordinate_system"]}
    return regions


def nominal_grid(rows: list[dict[str, object]], step: int = 250_000) -> list[dict[str, object]]:
    min_start = min(int(r["start"]) for r in rows)
    max_end = max(int(r["end"]) for r in rows)
    lookup = {int(r["start"]): r for r in rows}
    grid = []
    start = min_start
    while start <= max_end:
        end = start + step - 1
        row = lookup.get(start)
        grid.append({"start": start, "end": end, "observed": row is not None, "row": row})
        start += step
    return grid


def physical_null_for_lg(lg_rows: list[dict[str, object]], region: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object]]:
    grid = nominal_grid(lg_rows)
    width = int(region["end"]) - int(region["start"])
    observed_inside = [r for r in lg_rows if as_bool(r["inside_inversion"])]
    min_inside = max(3, math.ceil(0.5 * len(observed_inside)))
    obs_rows = eligible(lg_rows)
    observed_delta = delta([float(r["D_arrangement_minus_baseline"]) for r in obs_rows], [as_bool(r["inside_inversion"]) for r in obs_rows])
    out = []
    for idx, cell in enumerate(grid):
        cand_start = int(cell["start"])
        cand_end = cand_start + width
        if cand_end > int(grid[-1]["end"]):
            continue
        observed = [g["row"] for g in grid if g["observed"]]
        inside = [r for r in observed if int(r["start"]) >= cand_start and int(r["end"]) <= cand_end]
        outside = [r for r in observed if int(r["end"]) < cand_start or int(r["start"]) > cand_end]
        if len(inside) < min_inside or not outside:
            continue
        d = statistics.fmean(float(r["D_arrangement_minus_baseline"]) for r in inside) - statistics.fmean(float(r["D_arrangement_minus_baseline"]) for r in outside)
        out.append({"lg": lg_rows[0]["lg"], "candidate_index": idx, "candidate_start": cand_start, "candidate_end": cand_end, "n_observed_inside": len(inside), "n_observed_outside": len(outside), "delta_D": d, "greater_or_equal_observed": d >= observed_delta})
    vals = [float(r["delta_D"]) for r in out]
    summary = {
        "lg": lg_rows[0]["lg"],
        "observed_delta_D": observed_delta,
        "n_valid_positions": len(out),
        "minimum_observed_inside_required": min_inside,
        "physical_null_min": min(vals) if vals else math.nan,
        "physical_null_max": max(vals) if vals else math.nan,
        "observed_rank": sum(v > observed_delta for v in vals) + 1 if vals else "NA",
        "p_one_sided": sum(v >= observed_delta for v in vals) / len(vals) if vals else math.nan,
    }
    return out, summary


def physical_null(scores: list[dict[str, object]], regions: dict[str, dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    all_null = []
    summaries = []
    for lg in LGS:
        rows, summary = physical_null_for_lg([r for r in scores if r["lg"] == lg], regions[lg])
        all_null.extend(rows)
        summaries.append(summary)
    return all_null, summaries


def load_stage4b_tables() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, list[dict[str, str]]], dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    circular = {r["lg"]: r for r in read_tsv(RESULTS / "stage4b_circular_shift_summary.tsv")}
    block = {r["lg"]: r for r in read_tsv(RESULTS / "stage4b_block_placement_summary.tsv")}
    boundary = defaultdict(list)
    for r in read_tsv(RESULTS / "stage4b_boundary_sensitivity.tsv"):
        boundary[r["lg"]].append(r)
    topo = {r["lg"]: r for r in read_tsv(RESULTS / "stage4b_topology_run_summary.tsv")}
    composition = {}
    for r in read_tsv(RESULTS / "stage4a_topology_composition.tsv"):
        composition[f"{r['lg']}|{r['region_class']}"] = r
    qeff = {r["lg"]: r for r in read_tsv(RESULTS / "stage4a_LG_effect_summary.tsv")}
    return circular, block, boundary, topo, composition, qeff


def robustness_summary(primary: list[dict[str, object]], physical: list[dict[str, object]]) -> list[dict[str, object]]:
    circular, block, boundary, _topo, _composition, qeff = load_stage4b_tables()
    phys = {r["lg"]: r for r in physical}
    rows = []
    for row in primary:
        lg = str(row["lg"])
        boundary_ps = [float(r["p_one_sided"]) for r in boundary[lg]]
        boundary_ds = [float(r["observed_delta_D"]) for r in boundary[lg]]
        boundary_robust = all(d > 0 for d in boundary_ds) and max(boundary_ps) <= 0.06 if lg != "LG07" else all(abs(d) < 0.01 for d in boundary_ds)
        loo_rows = [r for r in read_tsv(RESULTS / "stage4b_leave_one_population_out.tsv") if r["lg"] == lg and r["delta_D"] != "NA"]
        loo_vals = [float(r["delta_D"]) for r in loo_rows]
        loo_robust = all(v > 0 for v in loo_vals) if lg != "LG07" else not all(v > 0 for v in loo_vals)
        circ_p = float(circular[lg]["p_one_sided"])
        block_p = float(block[lg]["p_one_sided"])
        phys_p = float(phys[lg]["p_one_sided"])
        bhp = float(row["BH_p"])
        if lg in {"LG01", "LG02"} and circ_p <= 0.02 and bhp <= 0.05 and block_p <= 0.03 and phys_p <= 0.05 and loo_robust:
            interp = "robust_positive"
        elif lg == "LG12" and float(row["delta_D"]) > 0 and circ_p <= 0.05 and block_p <= 0.05 and loo_robust:
            interp = "supportive_positive"
        else:
            interp = "weak_or_null"
        rows.append({
            "lg": lg,
            "delta_D": row["delta_D"],
            "circular_p": circ_p,
            "circular_rank": circular[lg]["observed_rank"],
            "block_p": block_p,
            "block_rank": block[lg]["observed_rank"],
            "physical_coordinate_p": phys_p,
            "physical_coordinate_rank": phys[lg]["observed_rank"],
            "boundary_robust": boundary_robust,
            "leave_one_population_out_robust": loo_robust,
            "final_interpretation": interp,
        })
    return rows


def informative_counts() -> dict[str, int]:
    counts = {}
    for r in read_tsv(RESULTS / "stage4a_LG_effect_summary.tsv"):
        counts[r["lg"]] = int(r["n_informative_quartets"])
    return counts


def make_main_table(primary: list[dict[str, object]], physical: list[dict[str, object]], robust: list[dict[str, object]]) -> list[dict[str, object]]:
    _circular, block, _boundary, _topo, _composition, qeff = load_stage4b_tables()
    counts = informative_counts()
    phys = {r["lg"]: r for r in physical}
    robust_map = {r["lg"]: r for r in robust}
    rows = []
    for r in primary:
        lg = str(r["lg"])
        rows.append({
            "lg": lg,
            "n_windows_inside": r["n_windows_inside"],
            "n_windows_outside": r["n_windows_outside"],
            "n_informative_quartets": counts[lg],
            "mean_D_inside": r["mean_D_inside"],
            "mean_D_outside": r["mean_D_outside"],
            "delta_D": r["delta_D"],
            "delta_A": r["delta_A"],
            "delta_B": r["delta_B"],
            "circular_p": r["circular_p"],
            "BH_p": r["BH_p"],
            "block_p": block[lg]["p_one_sided"],
            "physical_coordinate_p": phys[lg]["p_one_sided"],
            "fraction_quartets_positive": qeff[lg]["fraction_quartets_positive_effect"],
            "final_interpretation": robust_map[lg]["final_interpretation"],
        })
    return rows


def write_md_table(path: Path, rows: list[dict[str, object]]) -> None:
    lines = ["| LG | inside windows | outside windows | informative quartets | mean D inside | mean D outside | ΔD | circular p | BH p | block p | physical p | interpretation |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        lines.append(f"| {r['lg']} | {r['n_windows_inside']} | {r['n_windows_outside']} | {r['n_informative_quartets']} | {float(r['mean_D_inside']):.3f} | {float(r['mean_D_outside']):.3f} | {float(r['delta_D']):.3f} | {float(r['circular_p']):.3f} | {float(r['BH_p']):.3f} | {float(r['block_p']):.3f} | {float(r['physical_coordinate_p']):.3f} | {r['final_interpretation']} |")
    path.write_text("\n".join(lines) + "\n")


def plot_main(scores: list[dict[str, object]], primary: list[dict[str, object]], regions: dict[str, dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 7.8), sharey=True)
    for ax, lg in zip(axes, LGS):
        rows = [r for r in scores if r["lg"] == lg]
        x = [int(r["midpoint"]) / 1_000_000 for r in rows]
        y = [float(r["D_arrangement_minus_baseline"]) for r in rows]
        boundary_x = [int(r["midpoint"]) / 1_000_000 for r in rows if is_boundary(r)]
        boundary_y = [float(r["D_arrangement_minus_baseline"]) for r in rows if is_boundary(r)]
        ax.axvspan(int(regions[lg]["start"]) / 1_000_000, int(regions[lg]["end"]) / 1_000_000, color="#d9d9d9", alpha=0.7)
        ax.axhline(0, color="black", lw=0.8)
        ax.plot(x, y, color="#1f77b4", marker="o", markersize=2.5, lw=1.0)
        ax.scatter(boundary_x, boundary_y, color="#d62728", s=18, zorder=3)
        stat = next(r for r in primary if r["lg"] == lg)
        ax.text(0.99, 0.88, f"{lg}\nΔD={float(stat['delta_D']):.3f}\np={float(stat['circular_p']):.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=8, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
        ax.set_ylabel("D(w)")
        ax.set_ylim(-1.08, 1.08)
    axes[-1].set_xlabel("genomic position (Mb)")
    fig.suptitle("Atlantic cod local genealogy shift toward arrangement partition", y=0.995, fontsize=12)
    fig.tight_layout()
    fig.savefig(MAIN_FIG_PDF)
    fig.savefig(MAIN_FIG_PNG, dpi=300)
    plt.close(fig)


def plot_supp(primary, physical, main_table):
    circular, block, _boundary, _topo, comp, _qeff = load_stage4b_tables()
    loo = read_tsv(RESULTS / "stage4b_leave_one_population_out.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.5))
    ax = axes[0, 0]
    for i, lg in enumerate(LGS):
        nulls = [float(r["delta_D"]) for r in read_tsv(RESULTS / "stage4b_circular_shift_null.tsv") if r["lg"] == lg]
        ax.boxplot(nulls, positions=[i], widths=0.5, patch_artist=True, boxprops={"facecolor": "#d9d9d9"}, medianprops={"color": "black"})
        ax.scatter([i], [float(circular[lg]["observed_delta_D"])], color="#d62728", zorder=3)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(4)); ax.set_xticklabels(LGS); ax.set_ylabel("ΔD"); ax.set_title("A. Circular-shift null")
    ax = axes[0, 1]
    for i, lg in enumerate(LGS):
        vals = [float(r["delta_D"]) for r in loo if r["lg"] == lg and r["delta_D"] != "NA"]
        ax.scatter([i] * len(vals), vals, color="#4daf4a", s=16, alpha=0.8)
        ax.scatter([i], [float(circular[lg]["observed_delta_D"])], color="#d62728", marker="D", zorder=3)
    ax.axhline(0, color="black", lw=0.7); ax.set_xticks(range(4)); ax.set_xticklabels(LGS); ax.set_title("B. Leave-one-population-out"); ax.set_ylabel("ΔD")
    ax = axes[1, 0]
    width = 0.12
    topology_colors = {"arrangement": "#2ca02c", "baseline": "#9467bd", "third": "#8c564b"}
    for i, lg in enumerate(LGS):
        for j, cls in enumerate(["arrangement", "baseline", "third"]):
            inside = float(comp[f"{lg}|inside"][f"mean_fraction_{cls}"])
            outside = float(comp[f"{lg}|outside"][f"mean_fraction_{cls}"])
            ax.bar(i - 0.22 + j * width, inside, width=width, color=topology_colors[cls], alpha=0.9)
            ax.bar(i + 0.18 + j * width, outside, width=width, color=topology_colors[cls], alpha=0.35)
    color_handles = [Patch(facecolor=topology_colors["arrangement"], label="arrangement"),
                     Patch(facecolor=topology_colors["baseline"], label="baseline"),
                     Patch(facecolor=topology_colors["third"], label="third")]
    opacity_handles = [Patch(facecolor="gray", alpha=0.9, label="inside"),
                       Patch(facecolor="gray", alpha=0.35, label="outside")]
    leg1 = ax.legend(handles=color_handles, loc="upper left", fontsize=7, frameon=False, title="topology", title_fontsize=7)
    ax.add_artist(leg1)
    ax.legend(handles=opacity_handles, loc="upper right", fontsize=7, frameon=False, title="region", title_fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels(LGS); ax.set_ylim(0, 1); ax.set_ylabel("mean fraction"); ax.set_title("C. Topology composition")
    ax = axes[1, 1]
    phys_map = {r["lg"]: r for r in physical}
    for i, lg in enumerate(LGS):
        vals = [float(r["delta_D"]) for r in read_tsv(PHYS_NULL) if r["lg"] == lg]
        ax.boxplot(vals, positions=[i], widths=0.5, patch_artist=True, boxprops={"facecolor": "#c6dbef"}, medianprops={"color": "black"})
        ax.scatter([i], [float(phys_map[lg]["observed_delta_D"])], color="#d62728", zorder=3)
    ax.axhline(0, color="black", lw=0.7); ax.set_xticks(range(4)); ax.set_xticklabels(LGS); ax.set_ylabel("ΔD"); ax.set_title("D. Physical-coordinate null")
    fig.tight_layout()
    fig.savefig(SUPP_FIG_PDF)
    fig.savefig(SUPP_FIG_PNG, dpi=300)
    plt.close(fig)


def write_text_outputs(main_table, primary, physical, manifest_sha_placeholder=""):
    mt = {r["lg"]: r for r in main_table}
    MAIN_CAPTION.write_text(textwrap.dedent(f"""
    Figure. Arrangement-associated local genealogy shifts in Atlantic cod supergene regions. The Atlantic cod analysis treats the dataset as a multi-population test of recombination-partitioned genealogy structure. Population arrangement states were frozen from independent structural metadata, and the baseline population-history topology was frozen from collinear sequence before the published local trees were inspected. Each panel shows 250-kb windows on one focal linkage group. The y-axis is D(w)=fA−fB, where fA is the fraction of resolved informative quartets matching the arrangement partition and fB is the fraction matching the independently frozen baseline population-history partition. Gray shading marks the frozen inversion interval; the horizontal line marks D=0; red points indicate boundary-overlap windows. Exact one-sided circular-shift tests keep the ordered D(w) track intact while rotating it relative to the fixed inversion mask. LG01 and LG02 show strong shifts toward the arrangement partition within the inversion intervals (BH-adjusted p={float(mt['LG01']['BH_p']):.3f} and {float(mt['LG02']['BH_p']):.3f}). LG12 shows a supportive positive shift but has fewer informative quartets and weaker orientation evidence. LG07 shows no corresponding spatial shift, providing a contrasting inversion region within the same empirical system. These results establish spatial association, not causality.
    """).strip() + "\n")
    RESULTS_TEXT.write_text(textwrap.dedent(f"""
    We used the Atlantic cod supergene dataset as a multi-population test of whether recombination-partitioned chromosomal backgrounds are associated with persistent local genealogy structure. Arrangement-derived quartet predictions were frozen independently of local trees, and a baseline population-history prediction was frozen from collinear sequence. We then evaluated published 250-kb SNAPP population-tree summaries using a window-level contrast D(w)=fA−fB, where fA is the fraction of informative quartets matching the arrangement partition and fB is the fraction matching the baseline population-history partition.

    Local genealogies within LG01 and LG02 shifted strongly toward the independently defined arrangement partition within the inversion intervals. LG01 had mean inside D={float(mt['LG01']['mean_D_inside']):.3f}, ΔD={float(mt['LG01']['delta_D']):.3f}, exact circular-shift p={float(mt['LG01']['circular_p']):.3f}, and BH-adjusted p={float(mt['LG01']['BH_p']):.3f}. LG02 had ΔD={float(mt['LG02']['delta_D']):.3f}, exact circular-shift p={float(mt['LG02']['circular_p']):.3f}, and BH-adjusted p={float(mt['LG02']['BH_p']):.3f}. These effects were robust to contiguous-block and physical-coordinate-aware nulls and to leave-one-population-out summaries.

    LG12 also showed a positive arrangement-associated shift (ΔD={float(mt['LG12']['delta_D']):.3f}, raw circular-shift p={float(mt['LG12']['circular_p']):.3f}, BH-adjusted p={float(mt['LG12']['BH_p']):.3f}), but we retain the predeclared caveats that LG12 has only 16 informative quartets and that its ancestral/derived orientation evidence is weaker. LG07 showed no corresponding shift (ΔD={float(mt['LG07']['delta_D']):.3f}, p={float(mt['LG07']['circular_p']):.3f}), providing a contrasting inversion region in the same dataset.

    Overall, local genealogies within LG01 and LG02, and supportively within LG12, are unusually aligned with the arrangement partition relative to other positions on the same linkage group. These analyses do not establish causality, do not fit MSRC parameters, and do not make Atlantic cod a direct multi-species persistence-through-speciation test.
    """).strip() + "\n")
    METHODS_TEXT.write_text(textwrap.dedent("""
    We analyzed the Atlantic cod supergene dataset of Matschiner et al. (2022) using the authors' published 250-kb SNAPP population-tree summaries for LG01, LG02, LG07, and LG12. Inversion coordinates, population arrangement states, and the collinear baseline population tree were frozen before local window topologies were used. Informative quartets were predeclared as four-population subsets in which the arrangement-derived split differed from the baseline population-history split.

    For each local tree and informative quartet, we extracted the induced unrooted quartet split and classified it as arrangement-concordant, baseline-concordant, third topology, or unresolved. The primary window-level statistic was D(w)=fA−fB, where fA and fB are the fractions of resolved informative quartets matching the arrangement and baseline splits. For each linkage group, the primary effect was ΔD=mean(Dinside)−mean(Doutside), comparing fully inside inversion windows with fully outside windows and excluding boundary-overlap windows.

    Spatial significance was evaluated with an exact circular-shift null that preserves the ordered D(w) track and the fixed inversion mask while rotating their relative alignment. Exact one-sided p-values include the identity rotation. Benjamini-Hochberg correction was applied across the four primary linkage-group tests. Sensitivity analyses included contiguous-block placement, alternative boundary treatments, leave-one-population-out recalculation, full-tree topology-run summaries, and a physical-coordinate-aware null that shifts the inversion interval along a nominal 250-kb coordinate grid while retaining missing windows as missing.
    """).strip() + "\n")
    BORNHOLM_NOTE.write_text(textwrap.dedent("""
    # Bornholm LG12 validation note

    The targeted Bornholm positive-control window is `LG12_007500001_007750000`. This check was performed only after Stage-2 structural predictions and Stage-3 baseline predictions were frozen. In Stage 4A, Bornholm's closest local relationship in the target window differed from the immediate neighboring windows (`LG12_007250001_007500000` and `LG12_007750001_008000000`), recovering the local-switch behavior expected from the exchanged LG12 segment reported by Matschiner et al. This validation was not used to tune arrangement states, baseline splits, quartet definitions, or any spatial-null statistic.
    """).strip() + "\n")


def write_final_analysis(main_table, robust, physical, checksums, stage4b_checksums, output_paths):
    lines = ["# Atlantic cod final analysis freeze", "", "Atlantic cod empirical analysis is frozen for manuscript use.", "", "## Frozen input checksums", ""]
    for path, value in checksums.items():
        lines.append(f"- `{path}`: `{value}`")
    lines.extend(["", "## Primary Stage-4B results", "", "| LG | ΔD | circular p | BH p | interpretation |", "|---|---:|---:|---:|---|"])
    interp = {r["lg"]: r["final_interpretation"] for r in robust}
    for r in main_table:
        lines.append(f"| {r['lg']} | {float(r['delta_D']):.6f} | {float(r['circular_p']):.6f} | {float(r['BH_p']):.6f} | {interp[r['lg']]} |")
    lines.extend(["", "## Stage-5 physical-coordinate robustness", "", "| LG | p | rank |", "|---|---:|---:|"])
    for r in physical:
        lines.append(f"| {r['lg']} | {float(r['p_one_sided']):.6f} | {r['observed_rank']} |")
    lines.extend(["", "## Manuscript-ready files", ""])
    for path in output_paths:
        lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    FINAL_ANALYSIS.write_text("\n".join(lines) + "\n")


def write_manifest(checksums, stage4b_checksums, output_paths, primary, physical, regions):
    import matplotlib as mpl
    outputs = {str(p.relative_to(REPO_ROOT)): sha256(p) for p in output_paths if p.exists()}
    manifest = {
        "stage": 5,
        "analysis_status": "frozen_for_manuscript_use",
        "frozen_input_checksums": checksums,
        "stage4b_output_checksums": stage4b_checksums,
        "stage5_script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "stage5_script_sha256": sha256(Path(__file__).resolve()),
        "primary_statistic": "Delta_D = mean(D_inside) - mean(D_outside), D = f_A - f_B",
        "primary_null": "exact circular shifts of the ordered Stage-4A D(w) track with fixed inversion mask",
        "sensitivity_nulls": ["contiguous-block placement", "boundary treatment", "leave-one-population-out", "physical-coordinate-aware interval shift"],
        "inversion_coordinates": regions,
        "number_windows_per_lg": {lg: sum(1 for r in load_scores() if r["lg"] == lg) for lg in LGS},
        "informative_quartet_counts": informative_counts(),
        "physical_coordinate_summary": {r["lg"]: r for r in physical},
        "final_output_checksums": outputs,
        "software_versions": {"python": platform.python_version(), "biopython": Bio.__version__, "matplotlib": mpl.__version__},
        "reproducibility_commands": [
            "python3 empirical/atlantic_cod/scripts/05_finalize_cod_analysis.py --run-tests",
            "python3 empirical/atlantic_cod/scripts/05_finalize_cod_analysis.py",
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    checksums = verify_frozen()
    stage4b_checksums = {str(p.relative_to(REPO_ROOT)): sha256(p) for p in STAGE4B_OUTPUTS if p.exists()}
    scores = load_scores()
    regions = load_regions()
    primary = reproduce_stage4b(scores)
    phys_null, phys_summary = physical_null(scores, regions)
    write_tsv(PHYS_NULL, phys_null, ["lg", "candidate_index", "candidate_start", "candidate_end", "n_observed_inside", "n_observed_outside", "delta_D", "greater_or_equal_observed"])
    write_tsv(PHYS_SUMMARY, phys_summary, ["lg", "observed_delta_D", "n_valid_positions", "minimum_observed_inside_required", "physical_null_min", "physical_null_max", "observed_rank", "p_one_sided"])
    robust = robustness_summary(primary, phys_summary)
    write_tsv(ROBUSTNESS, robust, ["lg", "delta_D", "circular_p", "circular_rank", "block_p", "block_rank", "physical_coordinate_p", "physical_coordinate_rank", "boundary_robust", "leave_one_population_out_robust", "final_interpretation"])
    main_table = make_main_table(primary, phys_summary, robust)
    write_tsv(MAIN_TABLE_TSV, main_table, ["lg", "n_windows_inside", "n_windows_outside", "n_informative_quartets", "mean_D_inside", "mean_D_outside", "delta_D", "delta_A", "delta_B", "circular_p", "BH_p", "block_p", "physical_coordinate_p", "fraction_quartets_positive", "final_interpretation"])
    write_md_table(MAIN_TABLE_MD, main_table)
    plot_main(scores, primary, regions)
    plot_supp(primary, phys_summary, main_table)
    write_text_outputs(main_table, primary, phys_summary)
    manuscript_outputs = [MAIN_FIG_PDF, MAIN_FIG_PNG, SUPP_FIG_PDF, SUPP_FIG_PNG, MAIN_CAPTION, MAIN_TABLE_TSV, MAIN_TABLE_MD, RESULTS_TEXT, METHODS_TEXT, BORNHOLM_NOTE, PHYS_NULL, PHYS_SUMMARY, ROBUSTNESS]
    write_final_analysis(main_table, robust, phys_summary, checksums, stage4b_checksums, manuscript_outputs)
    all_outputs = manuscript_outputs + [FINAL_ANALYSIS]
    write_manifest(checksums, stage4b_checksums, all_outputs, primary, phys_summary, regions)


class Stage5Tests(unittest.TestCase):
    def test_checksum_verification(self):
        self.assertIn(str(STAGE2.relative_to(REPO_ROOT)), verify_frozen())

    def test_reproduce_stage4b_delta_and_p(self):
        rows = reproduce_stage4b(load_scores())
        lg01 = next(r for r in rows if r["lg"] == "LG01")
        self.assertAlmostEqual(lg01["delta_D"], 1.232766, places=5)
        self.assertAlmostEqual(lg01["circular_p"], 0.009174, places=5)

    def test_nominal_grid_construction(self):
        rows = [{"start": 1, "end": 250000}, {"start": 500001, "end": 750000}]
        grid = nominal_grid(rows)
        self.assertEqual(len(grid), 3)
        self.assertFalse(grid[1]["observed"])

    def test_missing_windows_remain_na_not_zero(self):
        rows = [{"start": 1, "end": 250000, "D_arrangement_minus_baseline": 1.0}, {"start": 500001, "end": 750000, "D_arrangement_minus_baseline": -1.0}]
        grid = nominal_grid(rows)
        self.assertIsNone(grid[1]["row"])

    def test_physical_coordinate_null_generation(self):
        rows = []
        for i in range(8):
            start = 1 + i * 250000
            rows.append({"lg": "LGX", "start": start, "end": start + 249999, "inside_inversion": "true" if 2 <= i <= 3 else "false", "overlaps_inversion": "true" if 2 <= i <= 3 else "false", "D_arrangement_minus_baseline": float(i)})
        null, summary = physical_null_for_lg(rows, {"start": 500001, "end": 1500000})
        self.assertGreater(len(null), 0)
        self.assertEqual(summary["n_valid_positions"], len(null))

    def test_deterministic_physical_null(self):
        scores = load_scores()
        regions = load_regions()
        a = physical_null(scores, regions)[1]
        b = physical_null(scores, regions)[1]
        self.assertEqual(a, b)

    def test_manuscript_table_consistency(self):
        scores = load_scores(); regions = load_regions(); primary = reproduce_stage4b(scores); phys = physical_null(scores, regions)[1]
        robust = robustness_summary(primary, phys)
        table = make_main_table(primary, phys, robust)
        self.assertEqual(len(table), 4)
        self.assertEqual(next(r for r in table if r["lg"] == "LG07")["final_interpretation"], "weak_or_null")

    def test_figure_input_consistency(self):
        scores = load_scores()
        self.assertEqual({r["lg"] for r in scores}, set(LGS))

    def test_manifest_checksum_determinism(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x"
            p.write_text("abc\n")
            self.assertEqual(sha256(p), sha256(p))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage5Tests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        sys.exit(0 if result.wasSuccessful() else 1)
    main()
