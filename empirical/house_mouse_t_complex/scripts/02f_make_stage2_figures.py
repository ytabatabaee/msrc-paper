#!/usr/bin/env python3
"""Build and validate the House mouse Stage 2 main figure."""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc_mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/msrc_xdgcache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from house_mouse_stage2_utils import FIGURES, RESULTS, read_tsv, rel


PNG = FIGURES / "house_mouse_t_complex_stage2_main.png"
PDF = FIGURES / "house_mouse_t_complex_stage2_main.pdf"
VALIDATION = RESULTS / "stage2_figure_validation.tsv"


def parse_args() -> argparse.Namespace:
    return argparse.ArgumentParser(description=__doc__).parse_args()


def moving_bin(rows: list[dict[str, str]], treatment: str, bin_size: int = 250_000) -> tuple[np.ndarray, np.ndarray]:
    sub = [r for r in rows if r["treatment"] == treatment]
    xs = np.array([float(r["midpoint_bp"]) for r in sub])
    ys = np.array([float(r["delta_species_alt"]) for r in sub])
    if np.nanstd(ys) == 0:
        raise RuntimeError(f"Nonconstant delta series expected for {treatment}")
    bins = np.arange(xs.min(), xs.max() + bin_size, bin_size)
    bx, by = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (xs >= lo) & (xs < hi)
        if mask.any():
            bx.append((lo + hi) / 2)
            by.append(float(np.nanmedian(ys[mask])))
    return np.array(bx), np.array(by)


def run() -> None:
    scan = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    qsum = read_tsv(RESULTS / "stage2_fixed_quartet_summary.tsv")
    balanced = read_tsv(RESULTS / "stage2_balanced_sampling_summary.tsv") if (RESULTS / "stage2_balanced_sampling_summary.tsv").exists() else []
    linkage = read_tsv(RESULTS / "stage2_linkage_sensitivity.tsv")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax = axes[0, 0]
    for treatment, color in [("STANDARD_ONLY", "#4c78a8"), ("ALL_TIPS", "#f58518")]:
        sub = [r for r in scan if r["treatment"] == treatment]
        xs = np.array([float(r["midpoint_bp"]) / 1e6 for r in sub])
        ys = np.array([float(r["delta_species_alt"]) for r in sub])
        ax.plot(xs, ys, color=color, alpha=0.18, linewidth=0.5)
        bx, by = moving_bin(scan, treatment)
        ax.plot(bx / 1e6, by, color=color, linewidth=2.2, label=treatment.replace("_", " "))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("chr17 position (Mb)")
    ax.set_ylabel("q_species - q_t_alt")
    ax.set_title("A. Spatial fixed-quartet support")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    treatments = ["STANDARD_ONLY", "T_ONLY", "ALL_TIPS"]
    comp = ["mean_q_species", "mean_q_t_alt", "mean_q_other", "mean_q_unresolved"]
    labels = ["Q_SPECIES", "Q_T_ALT", "Q_OTHER", "unresolved"]
    bottoms = np.zeros(len(treatments))
    qmap = {r["treatment"]: r for r in qsum}
    colors = ["#4c78a8", "#f58518", "#54a24b", "#bab0ac"]
    for c, label, color in zip(comp, labels, colors):
        vals = np.array([float(qmap[t][c]) for t in treatments])
        ax.bar(treatments, vals, bottom=bottoms, label=label, color=color)
        bottoms += vals
    ax.set_ylim(0, 1)
    ax.set_title("B. Aggregate quartet composition")
    ax.tick_params(axis="x", labelrotation=20)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    if balanced:
        bmap = {r["treatment"]: r for r in balanced}
        btreat = ["B0_STANDARD_MATCHED", "B1_T_MATCHED", "B2_MIXED_BALANCED"]
        bottoms = np.zeros(len(btreat))
        for key, label, color in [
            ("fraction_Q_SPECIES", "Q_SPECIES", "#4c78a8"),
            ("fraction_Q_T_ALT", "Q_T_ALT", "#f58518"),
            ("fraction_Q_OTHER", "Q_OTHER", "#54a24b"),
        ]:
            vals = np.array([float(bmap.get(t, {}).get(key, 0)) for t in btreat])
            ax.bar(btreat, vals, bottom=bottoms, label=label, color=color)
            bottoms += vals
        ax.set_ylim(0, 1)
    ax.set_title("C. Balanced ASTRAL controls")
    ax.tick_params(axis="x", labelrotation=20)

    ax = axes[1, 1]
    by_spacing: dict[int, list[dict[str, str]]] = defaultdict(list)
    for r in linkage:
        by_spacing[int(r["spacing_kb"])].append(r)
    spacings, pp, ntree = [], [], []
    for spacing in sorted(by_spacing):
        rows = by_spacing[spacing]
        spacings.append(spacing)
        pp.append(np.mean([float(r["localPP"]) for r in rows if r["localPP"] != "NA"]))
        ntree.append(np.mean([float(r["n_trees"]) for r in rows]))
    ax.plot(spacings, pp, marker="o", color="#e45756", label="mean localPP")
    ax2 = ax.twinx()
    ax2.plot(spacings, ntree, marker="s", color="#72b7b2", label="mean n trees")
    ax.set_xscale("log")
    ax.set_xlabel("spacing (kb)")
    ax.set_ylabel("localPP")
    ax2.set_ylabel("n trees")
    ax.set_title("D. Spatial thinning sensitivity")

    fig.tight_layout()
    fig.savefig(PNG, dpi=220)
    fig.savefig(PDF)
    validate()
    print(f"Wrote {PNG}")
    print(f"Wrote {PDF}")


def validate() -> None:
    if PNG.stat().st_size <= 20_000:
        raise RuntimeError("Main PNG is too small.")
    if PDF.stat().st_size <= 10_000:
        raise RuntimeError("Main PDF is too small.")
    img = Image.open(PNG).convert("RGB")
    arr = np.asarray(img)
    nonwhite = np.sum(np.any(arr < 245, axis=2))
    distinct = len(np.unique(arr.reshape(-1, 3), axis=0))
    if img.width < 800 or img.height < 600 or nonwhite < 5000 or distinct < 100:
        raise RuntimeError("Main PNG failed nonempty pixel validation.")
    VALIDATION.write_text(
        "file\twidth\theight\tsize_bytes\tnonwhite_pixels\tdistinct_rgb_values\n"
        f"{rel(PNG)}\t{img.width}\t{img.height}\t{PNG.stat().st_size}\t{nonwhite}\t{distinct}\n"
        f"{rel(PDF)}\tNA\tNA\t{PDF.stat().st_size}\tNA\tNA\n"
    )


if __name__ == "__main__":
    raise SystemExit(run())
