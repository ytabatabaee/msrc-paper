#!/usr/bin/env python3
"""Create manuscript figures from the frozen processed quartet table."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np


def load(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _topology_icon(ax, x, y, split, color, emphasized=False):
    """Draw a compact two-pair split icon in axes coordinates."""
    lw = 2.0 if emphasized else 1.2
    left_x, right_x = x - 0.07, x + 0.07
    ax.plot([left_x, left_x, x, right_x, right_x], [y - 0.035, y + 0.035, y, y + 0.035, y - 0.035], color=color, lw=lw, solid_capstyle="round")
    ax.text(x, y - 0.075, split, ha="center", va="top", fontsize=8.5, color=color, fontweight="bold" if emphasized else "normal")


def make_main_panel(data, figdir: Path):
    """Create the compact manuscript Panel A from exact processed values."""
    durations = sorted({float(r["duration"]) for r in data})
    m_values = sorted({float(r["m"]) for r in data})
    lookup = {(float(r["duration"]), float(r["m"])): r for r in data}
    z = np.array([[float(lookup[(duration, m)]["exact_delta_alt_species"]) for duration in durations] for m in m_values])
    strongest = max(data, key=lambda r: float(r["exact_delta_alt_species"]))

    fig = plt.figure(figsize=(11.2, 5.5), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[0.92, 1.28], wspace=0.08, left=0.035, right=0.965, top=0.88, bottom=0.13)
    schematic = fig.add_subplot(gs[0, 0])
    result = fig.add_subplot(gs[0, 1])
    schematic.set_xlim(0, 1); schematic.set_ylim(0, 1); schematic.axis("off")
    schematic.text(0.02, 1.04, "A1  Scenario", fontsize=11, fontweight="bold", transform=schematic.transAxes)
    schematic.text(0.02, 0.985, "Conditional quartet experiment", fontsize=13, fontweight="bold", transform=schematic.transAxes)
    schematic.text(0.02, 0.925, "taxon order 1, 2, 3, 4     configuration 1010", fontsize=9.5, transform=schematic.transAxes, color="#333333")

    schematic.add_patch(Rectangle((0.17, 0.28), 0.64, 0.52, facecolor="#f6f7f8", edgecolor="#a9afb5", lw=1.0, linestyle="--"))
    schematic.text(0.49, 0.83, "structured interval, duration t", ha="center", va="bottom", fontsize=9, color="#4d555d")
    schematic.text(0.07, 0.67, "state 1", ha="left", va="center", fontsize=9, color="#a33b2a", fontweight="bold")
    schematic.text(0.07, 0.42, "state 0", ha="left", va="center", fontsize=9, color="#455b76", fontweight="bold")
    for x, taxon in [(0.29, "1"), (0.59, "3")]:
        schematic.plot(x, 0.67, "o", ms=11, color="#c85a45", mec="#7f281d", mew=1.0)
        schematic.text(x, 0.60, f"taxon {taxon}", ha="center", va="top", fontsize=8)
    for x, taxon in [(0.29, "2"), (0.59, "4")]:
        schematic.plot(x, 0.42, "o", ms=11, mfc="white", mec="#455b76", mew=1.5)
        schematic.text(x, 0.35, f"taxon {taxon}", ha="center", va="top", fontsize=8)
    schematic.add_patch(FancyArrowPatch((0.66, 0.61), (0.66, 0.48), arrowstyle="<->", mutation_scale=12, lw=1.2, color="#555555"))
    schematic.text(0.70, 0.545, "m01 = m10 = m", ha="left", va="center", fontsize=8.5, color="#444444")
    schematic.text(0.49, 0.255, "same-background coalescence only", ha="center", va="top", fontsize=8.5, color="#555555")
    schematic.text(0.49, 0.20, "lambda0 = lambda1 = 1", ha="center", va="top", fontsize=8.5, color="#555555")
    schematic.text(0.49, 0.105, "if unresolved after t, resolve symmetrically\namong the three quartet topologies", ha="center", va="center", fontsize=8.2, color="#555555")
    _topology_icon(schematic, 0.20, 0.08, "12|34  Q_SPECIES", "#1b4f9c")
    _topology_icon(schematic, 0.50, 0.08, "13|24  Q_ALT", "#c23b22", emphasized=True)
    _topology_icon(schematic, 0.80, 0.08, "14|23  Q_OTHER", "#555555")

    result.set_title("A2  Exact q_ALT − q_SPECIES", loc="left", fontsize=11, fontweight="bold", pad=12)
    image = result.imshow(z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-np.max(np.abs(z)), vmax=np.max(np.abs(z)))
    result.set_xticks(range(len(durations)), [f"{x:g}" for x in durations])
    result.set_yticks(range(len(m_values)), [f"{x:g}" for x in m_values])
    result.set_xlabel("structured interval duration t")
    result.set_ylabel("symmetric switching rate m")
    result.tick_params(labelsize=8)
    result.grid(color="white", lw=0.5, alpha=0.45)
    max_col = durations.index(float(strongest["duration"]))
    max_row = m_values.index(float(strongest["m"]))
    result.add_patch(Rectangle((max_col - 0.5, max_row - 0.5), 1, 1, fill=False, edgecolor="#111111", lw=2.0))
    result.annotate("strongest cell\nt=2, m=0\n(0.0061, 0.9878, 0.0061)", xy=(max_col, max_row), xytext=(max_col - 1.5, max_row + 1.15), fontsize=8, ha="center", va="bottom", arrowprops={"arrowstyle": "-", "color": "#222222", "lw": 0.8}, bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#777777", "alpha": 0.9})
    result.text(0.98, 0.03, "Q_ALT > Q_SPECIES throughout tested grid", transform=result.transAxes, ha="right", va="bottom", fontsize=8.5, color="#7d281d", bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "none", "alpha": 0.85})
    cbar = fig.colorbar(image, ax=result, fraction=0.045, pad=0.035)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("q_ALT − q_SPECIES", fontsize=8.5)
    fig.suptitle("Conditional quartet experiment", x=0.035, ha="left", fontsize=15, fontweight="bold")
    fig.savefig(figdir / "conditional_quartet_main_panel.png", dpi=300)
    fig.savefig(figdir / "conditional_quartet_main_panel.pdf")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[4])
    args = ap.parse_args()
    root = args.repo_root.resolve()
    data = load(root / "simulations/quartet/datasets/processed/conditional_quartet_grid.tsv")
    figdir = root / "simulations/quartet/analysis/figures"
    figdir.mkdir(parents=True, exist_ok=True)
    make_main_panel(data, figdir)
    durations = [0.05, 0.10, 0.25, 0.50, 1.00, 2.00]
    colors = {"q_species": "#1b4f9c", "q_alt": "#c23b22", "q_other": "#555555"}
    labels = {"q_species": "Q_SPECIES (12|34)", "q_alt": "Q_ALT (13|24)", "q_other": "Q_OTHER (14|23)"}
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True, constrained_layout=True)
    for ax, duration in zip(axes.flat, durations):
        subset = sorted((r for r in data if float(r["duration"]) == duration), key=lambda r: float(r["m"]))
        x = np.array([float(r["m"]) for r in subset])
        for key in colors:
            ax.plot(x, [float(r[f"exact_{key}"]) for r in subset], color=colors[key], label=labels[key], lw=2)
            ax.scatter(x, [float(r[f"empirical_{key}"]) for r in subset], color=colors[key], s=10, alpha=0.35)
        ax.axhline(1/3, color="black", ls="--", lw=0.8)
        ax.set_title(f"duration = {duration:g}")
        ax.set_xscale("symlog", linthresh=0.0025)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.2)
    for ax in axes[-1]: ax.set_xlabel("symmetric switching rate m")
    for ax in axes[:, 0]: ax.set_ylabel("quartet probability")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(figdir / "conditional_quartet_curves.png", dpi=300)
    fig.savefig(figdir / "conditional_quartet_curves.pdf")
    plt.close(fig)

    durations_all = sorted({float(r["duration"]) for r in data})
    m_all = sorted({float(r["m"]) for r in data})
    z = np.array([[float(next(r for r in data if float(r["duration"]) == d and float(r["m"]) == m)["exact_delta_alt_species"]) for d in durations_all] for m in m_all])
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    im = ax.imshow(z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-np.max(np.abs(z)), vmax=np.max(np.abs(z)), extent=[min(durations_all), max(durations_all), min(m_all), max(m_all)])
    if np.min(z) < 0 < np.max(z): ax.contour(durations_all, m_all, z, levels=[0], colors="black", linewidths=1)
    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=0.0025)
    ax.set_xlabel("structured interval duration"); ax.set_ylabel("symmetric switching rate m")
    ax.set_title("Arrangement-conditioned quartet bias")
    fig.colorbar(im, ax=ax, label="exact Q_ALT − Q_SPECIES")
    fig.savefig(figdir / "conditional_quartet_bias_heatmap.png", dpi=300)
    fig.savefig(figdir / "conditional_quartet_bias_heatmap.pdf")
    plt.close(fig)

    # Barycentric projection: x = q_alt + q_other/2, y = sqrt(3) q_other/2.
    fig, axes = plt.subplots(2, 4, figsize=(10, 5.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for ax, duration in zip(axes.flat, durations_all):
        subset = sorted((r for r in data if float(r["duration"]) == duration), key=lambda r: float(r["m"]))
        xs = [float(r["exact_q_alt"]) + float(r["exact_q_other"]) / 2 for r in subset]
        ys = [np.sqrt(3) * float(r["exact_q_other"]) / 2 for r in subset]
        ax.plot(xs, ys, color="#888888", lw=0.8)
        ax.scatter(xs, ys, c=[float(r["m"]) for r in subset], cmap=cmap, s=24)
        ax.scatter([0.5], [np.sqrt(3)/6], color="black", marker="+", s=55)
        ax.set_title(f"t = {duration:g}", fontsize=9, loc="left")
        ax.set_xlim(0, 1); ax.set_ylim(0, np.sqrt(3)/2); ax.set_aspect("equal"); ax.axis("off")
        ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3)/2, 0], color="black", lw=0.8)
    axes.flat[-1].axis("off")
    axes.flat[0].text(0, -0.08, "Q_SPECIES", transform=axes.flat[0].transAxes, ha="left")
    axes.flat[0].text(1, -0.08, "Q_ALT", transform=axes.flat[0].transAxes, ha="right")
    axes.flat[0].text(0.5, 0.92, "Q_OTHER", transform=axes.flat[0].transAxes, ha="center", va="center", fontsize=8.5, bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.8})
    fig.suptitle("Conditional quartet probability simplex")
    fig.savefig(figdir / "conditional_quartet_simplex.png", dpi=300)
    fig.savefig(figdir / "conditional_quartet_simplex.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
