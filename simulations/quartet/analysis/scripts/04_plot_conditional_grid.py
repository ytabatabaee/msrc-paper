#!/usr/bin/env python3
"""Create manuscript figures from the frozen processed quartet table."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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


def _phylo_tree(ax, x, y, width, height, split, color="#333333", labels=False, lw=1.0):
    """Draw a small rooted quartet tree with the requested split."""
    pair_order = {"12|34": ((0, 1), (2, 3)), "13|24": ((0, 2), (1, 3)), "14|23": ((0, 3), (1, 2))}[split]
    tip_y = [y + height * 0.82, y + height * 0.58, y + height * 0.34, y + height * 0.10]
    root = (x + width * 0.50, y + height * 0.46)
    internal = [(x + width * 0.28, np.mean([tip_y[i] for i in pair_order[0]])), (x + width * 0.72, np.mean([tip_y[i] for i in pair_order[1]]))]
    for node in internal:
        ax.plot([root[0], node[0]], [root[1], node[1]], color=color, lw=lw, solid_capstyle="round")
    for node, pair in zip(internal, pair_order):
        for i in pair:
            ax.plot([node[0], x + width * 0.93], [node[1], tip_y[i]], color=color, lw=lw, solid_capstyle="round")
            if labels:
                ax.text(x + width * 0.96, tip_y[i], str(i + 1), ha="left", va="center", fontsize=8, color=color)


def _draw_arrangement_bar(ax, x, y, width, height, state):
    ax.add_patch(Rectangle((x, y), width, height, facecolor="#f3f4f5", edgecolor="#9ca3aa", lw=0.7))
    segment_color = "#c85a45" if state == 1 else "#57708f"
    hatch = "///" if state == 1 else "\\\\\\"
    ax.add_patch(Rectangle((x + width * 0.62, y + height * 0.12), width * 0.27, height * 0.76, facecolor=segment_color, edgecolor=segment_color, hatch=hatch, alpha=0.9, lw=0.6))


def _draw_schematic(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.00, 1.04, "A  Biological setup", fontsize=11, fontweight="bold", transform=ax.transAxes)
    ax.text(0.00, 0.985, "Species tree Q_SPECIES = 12|34", fontsize=10.5, fontweight="bold", transform=ax.transAxes)
    ax.text(0.57, 0.985, "conditioned arrangement configuration 1010", fontsize=9.2, color="#4b5563", transform=ax.transAxes)
    ax.text(0.57, 0.945, "recombination-suppressed region", fontsize=8.4, color="#4b5563", transform=ax.transAxes)
    _phylo_tree(ax, 0.03, 0.50, 0.42, 0.35, "12|34", color="#263746", labels=True, lw=1.4)
    ax.text(0.03, 0.47, "species tree", fontsize=8.3, color="#4b5563")
    bar_x, bar_w, bar_h = 0.58, 0.35, 0.055
    tip_y = [0.79, 0.68, 0.57, 0.46]
    states = [1, 0, 1, 0]
    for i, (yy, state) in enumerate(zip(tip_y, states), start=1):
        ax.plot([0.49, bar_x], [yy + bar_h / 2, yy + bar_h / 2], color="#9ca3aa", lw=0.8)
        _draw_arrangement_bar(ax, bar_x, yy, bar_w, bar_h, state)
        ax.text(0.95, yy + bar_h / 2, f"{i}: {state}", ha="left", va="center", fontsize=8.2, color="#7f281d" if state else "#455b76", fontweight="bold")
    ax.text(bar_x + bar_w * 0.755, 0.865, "1010", ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#374151")
    ax.text(bar_x + bar_w * 0.755, 0.415, "state 1: reversed/patterned block   state 0: baseline block", ha="center", va="top", fontsize=7.1, color="#4b5563")
    ax.text(0.03, 0.36, "structured interval duration t; symmetric switching m", fontsize=8.5, color="#4b5563")
    ax.text(0.03, 0.315, "conditioning fixes 1010; no forward origin is simulated", fontsize=8.1, color="#7f281d")
    ax.text(0.03, 0.23, "local quartet gene-tree outcomes", fontsize=8.5, color="#4b5563")
    for xx, split, color, emphasis in [(0.05, "12|34", "#1b4f9c", False), (0.36, "13|24", "#c23b22", True), (0.67, "14|23", "#555555", False)]:
        _phylo_tree(ax, xx, 0.035, 0.24, 0.16, split, color=color, lw=1.2 if emphasis else 0.9)
        ax.text(xx + 0.12, 0.005, f"{split}  { {'12|34':'Q_SPECIES','13|24':'Q_ALT','14|23':'Q_OTHER'}[split]}", ha="center", va="bottom", fontsize=7.7, color=color, fontweight="bold" if emphasis else "normal")


def _make_curves_panel(fig, spec, data, selected):
    sub = spec.subgridspec(2, 2, wspace=0.28, hspace=0.35)
    axes = np.array([[fig.add_subplot(sub[i, j]) for j in range(2)] for i in range(2)])
    colors = {"q_species": "#1b4f9c", "q_alt": "#c23b22", "q_other": "#555555"}
    labels = {"q_species": "Q_SPECIES = 12|34", "q_alt": "Q_ALT = 13|24", "q_other": "Q_OTHER = 14|23"}
    for index, (ax, duration) in enumerate(zip(axes.flat, selected)):
        subset = sorted((r for r in data if float(r["duration"]) == duration), key=lambda r: float(r["m"]))
        x = np.array([float(r["m"]) for r in subset])
        for key in ("q_species", "q_alt", "q_other"):
            style = "--" if key == "q_other" else "-"
            marker = "x" if key == "q_other" else "o"
            ax.plot(x, [float(r[f"exact_{key}"]) for r in subset], color=colors[key], lw=1.5, ls=style, label=labels[key])
            ax.scatter(x, [float(r[f"empirical_{key}"]) for r in subset], color=colors[key], s=12, marker=marker, alpha=0.75, zorder=3)
        ax.axhline(1 / 3, color="#777777", ls=":", lw=0.8)
        ax.set_title("B  Quartet probability curves" if index == 0 else f"t = {duration:g}", loc="left", fontsize=9.5 if index == 0 else 9, fontweight="bold")
        ax.set_xscale("symlog", linthresh=0.0025)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.18)
        ax.tick_params(labelsize=7)
        if index >= 2: ax.set_xlabel("symmetric switching rate m", fontsize=8)
        if index % 2 == 0: ax.set_ylabel("probability", fontsize=8)
    axes[0, 0].text(0.03, 0.08, "Q_SPECIES = Q_OTHER by symmetry", transform=axes[0, 0].transAxes, fontsize=6.8, color="#555555", bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.8})
    axes[0, 0].text(0.03, 0.88, "t = 0.05", transform=axes[0, 0].transAxes, fontsize=7.5, color="#333333")
    axes[0, 0].legend(handles=[Line2D([0], [0], color=colors[k], ls="--" if k == "q_other" else "-", marker="x" if k == "q_other" else "o", markersize=4, lw=1.4, label=labels[k]) for k in ("q_species", "q_alt", "q_other")], fontsize=6.8, frameon=False, loc="upper right")


def _make_simplex_panel(fig, spec, data, selected):
    sub = spec.subgridspec(2, 2, wspace=0.16, hspace=0.20)
    axes = np.array([[fig.add_subplot(sub[i, j]) for j in range(2)] for i in range(2)])
    m_values = sorted({float(r["m"]) for r in data})
    norm = plt.Normalize(min(m_values), max(m_values))
    cmap = plt.get_cmap("viridis")
    for index, (ax, duration) in enumerate(zip(axes.flat, selected)):
        subset = sorted((r for r in data if float(r["duration"]) == duration), key=lambda r: float(r["m"]))
        xs = np.array([float(r["exact_q_alt"]) + float(r["exact_q_other"]) / 2 for r in subset])
        ys = np.array([np.sqrt(3) * float(r["exact_q_other"]) / 2 for r in subset])
        ax.plot(xs, ys, color="#8b949e", lw=0.9, zorder=1)
        ax.scatter(xs, ys, c=[float(r["m"]) for r in subset], cmap=cmap, norm=norm, s=23, edgecolor="white", linewidth=0.25, zorder=2)
        ax.scatter([0.5], [np.sqrt(3) / 6], marker="+", color="black", s=50, zorder=3)
        ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3) / 2, 0], color="#30363d", lw=0.8)
        ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.04, np.sqrt(3) / 2 + 0.04); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title("C  Quartet simplex" if index == 0 else f"t = {duration:g}", loc="left", fontsize=9.5 if index == 0 else 9, fontweight="bold", pad=1)
        if index == 0:
            ax.text(0.03, 0.86, "t = 0.05", transform=ax.transAxes, fontsize=7.5, color="#333333")
            ax.text(0.01, -0.04, "Q_SPECIES", transform=ax.transAxes, ha="left", va="top", fontsize=7.3, color="#1b4f9c")
            ax.text(0.99, -0.04, "Q_ALT", transform=ax.transAxes, ha="right", va="top", fontsize=7.3, color="#c23b22")
            ax.text(0.50, 0.93, "Q_OTHER", transform=ax.transAxes, ha="center", va="top", fontsize=7.3, color="#555555")
            ax.text(0.55, 0.14, "low m", transform=ax.transAxes, fontsize=6.8, color="#555555")
            ax.text(0.22, 0.23, "high m", transform=ax.transAxes, fontsize=6.8, color="#555555")
            ax.add_patch(FancyArrowPatch((0.61, 0.17), (0.28, 0.27), transform=ax.transAxes, arrowstyle="->", mutation_scale=8, lw=0.8, color="#555555"))
        if index == 1:
            ax.text(0.52, 0.08, "Q_SPECIES = Q_OTHER symmetry line", transform=ax.transAxes, ha="center", fontsize=6.6, color="#555555")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.035, pad=0.02, shrink=0.78)
    cb.set_label("m", fontsize=8); cb.ax.tick_params(labelsize=7)


def _make_heatmap_panel(fig, spec, data):
    ax = fig.add_subplot(spec)
    durations = sorted({float(r["duration"]) for r in data})
    m_values = sorted({float(r["m"]) for r in data})
    lookup = {(float(r["duration"]), float(r["m"])): r for r in data}
    z = np.array([[float(lookup[(duration, m)]["exact_delta_alt_species"]) for duration in durations] for m in m_values])
    image = ax.imshow(z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-np.max(np.abs(z)), vmax=np.max(np.abs(z)))
    ax.set_xticks(range(len(durations)), [f"{x:g}" for x in durations], fontsize=7)
    ax.set_yticks(range(len(m_values)), [f"{x:g}" for x in m_values], fontsize=7)
    ax.set_xlabel("structured interval duration t", fontsize=8)
    ax.set_ylabel("symmetric switching rate m", fontsize=8)
    ax.set_title("D  Full parameter grid: strength of arrangement-conditioned bias", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    ax.grid(color="white", lw=0.5, alpha=0.45)
    strongest = max(data, key=lambda r: float(r["exact_delta_alt_species"]))
    col = durations.index(float(strongest["duration"])); row = m_values.index(float(strongest["m"]))
    ax.add_patch(Rectangle((col - 0.5, row - 0.5), 1, 1, fill=False, edgecolor="#111111", lw=1.8))
    ax.annotate("strongest\nt=2, m=0", xy=(col, row), xytext=(col - 1.1, row + 1.0), fontsize=7, ha="center", va="bottom", arrowprops={"arrowstyle": "-", "color": "#222222", "lw": 0.7}, bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "#777777", "alpha": 0.9})
    ax.text(0.98, 0.03, "Q_ALT > Q_SPECIES throughout", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#7d281d", bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.85})
    cb = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("exact q_ALT − q_SPECIES", fontsize=7.8); cb.ax.tick_params(labelsize=7)


def make_main_panel(data, figdir: Path):
    """Create the compact manuscript Panel A from exact processed values."""
    selected = [0.05, 0.25, 0.5, 2.0]
    fig = plt.figure(figsize=(14.2, 10.0), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[0.92, 1.08], height_ratios=[1.0, 1.0], wspace=0.16, hspace=0.22, left=0.04, right=0.97, top=0.95, bottom=0.07)
    _draw_schematic(fig.add_subplot(gs[0, 0]))
    _make_curves_panel(fig, gs[0, 1], data, selected)
    _make_simplex_panel(fig, gs[1, 0], data, selected)
    _make_heatmap_panel(fig, gs[1, 1], data)
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
