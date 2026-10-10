#!/usr/bin/env python3
"""Plot the saved mechanistic quartet pilot and matched MSC controls."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np


TOPOLOGY_COLORS = {"Q_S": "#2c7fb8", "Q_A": "#d95f02", "Q_O": "#636363"}
TOPOLOGY_LABELS = {"Q_S": "Q_S = 12|34", "Q_A": "Q_A = 13|24", "Q_O": "Q_O = 14|23"}


def read_tsv(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def tree_icon(ax, x, y, topology, scale=1.0, color="#333333", lw=1.2):
    pairs = {"12|34": ((0, 0.35), (0, -0.35), (0.6, 0)), "13|24": ((0, 0.35), (0, -0.35), (0.6, 0)), "14|23": ((0, 0.35), (0, -0.35), (0.6, 0))}
    split = {"12|34": ((1, 2), (3, 4)), "13|24": ((1, 3), (2, 4)), "14|23": ((1, 4), (2, 3))}[topology]
    ys = {1: y + 0.75 * scale, 2: y + 0.25 * scale, 3: y - 0.25 * scale, 4: y - 0.75 * scale}
    x0, x1, x2 = x, x + 0.35 * scale, x + 0.7 * scale
    for pair, center in zip(split, (0.5, -0.5)):
        ya, yb = ys[pair[0]], ys[pair[1]]
        yc = (ya + yb) / 2
        ax.plot([x0, x1], [ya, ya], color=color, lw=lw)
        ax.plot([x0, x1], [yb, yb], color=color, lw=lw)
        ax.plot([x1, x1], [ya, yb], color=color, lw=lw)
        ax.plot([x1, x2], [yc, yc], color=color, lw=lw)
    ax.plot([x2, x2 + 0.25 * scale], [0, 0], color=color, lw=lw)
    for taxon, yy in ys.items():
        ax.text(x0 - 0.06 * scale, yy, str(taxon), ha="right", va="center", fontsize=7)


def panel_a(ax):
    ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.01, 0.98, "Mechanistic quartet scenario", transform=ax.transAxes, va="top", fontsize=11, fontweight="bold")
    ax.text(0.01, 0.91, "Saved cell Ne=20, τ=.1, origin=.75, copy=.2, cross=.1; one replicate yielded 1010", transform=ax.transAxes, va="top", fontsize=7.2, color="#444444")
    # species tree
    x0, x1, x2 = 0.08, 0.22, 0.36
    y = {1: 0.70, 2: 0.54, 3: 0.30, 4: 0.14}
    for a, b, ya, yb in [(1, 2, y[1], y[2]), (3, 4, y[3], y[4])]:
        yc = (ya + yb) / 2
        ax.plot([x0, x1], [ya, ya], color="#333333", lw=1.5); ax.plot([x0, x1], [yb, yb], color="#333333", lw=1.5)
        ax.plot([x1, x1], [ya, yb], color="#333333", lw=1.5); ax.plot([x1, x2], [yc, yc], color="#333333", lw=1.5)
    ax.plot([x2, x2], [(y[1]+y[2])/2, (y[3]+y[4])/2], color="#333333", lw=1.5)
    ax.plot([x2, 0.43], [0.5, 0.5], color="#333333", lw=1.5)
    for k, yy in y.items(): ax.text(x0 - 0.02, yy, str(k), ha="right", va="center", fontsize=8)
    ax.text(0.18, 0.78, "species tree", fontsize=8, ha="center", fontweight="bold")
    ax.text(0.26, 0.45, "12|34", fontsize=8, color=TOPOLOGY_COLORS["Q_S"], fontweight="bold", ha="center")
    # genomic bars and sampled states from one saved replicate
    ax.text(0.50, 0.78, "arrangement history in a suppressed region", fontsize=8, ha="center", fontweight="bold")
    ax.text(0.50, 0.72, "sampled terminal pattern 1010 (not imposed)", fontsize=7, ha="center", color="#444444")
    for i, (taxon, state) in enumerate(zip((1, 2, 3, 4), "1010")):
        yy = y[taxon]
        ax.text(0.47, yy, f"taxon {taxon}", fontsize=7, ha="right", va="center")
        ax.add_patch(Rectangle((0.49, yy - 0.035), 0.20, 0.07, facecolor="#f4f4f4", edgecolor="#555555", lw=0.7))
        ax.add_patch(Rectangle((0.59, yy - 0.035), 0.07, 0.07, facecolor="#f28e2b" if state == "1" else "#ffffff", edgecolor="#333333", lw=0.8, hatch="///" if state == "1" else ""))
        ax.text(0.675, yy, f"state {state}", fontsize=7, va="center")
    ax.add_patch(FancyArrowPatch((0.43, 0.50), (0.46, 0.64), arrowstyle="->", mutation_scale=10, lw=1, color="#555555"))
    ax.text(0.42, 0.60, "ROOT origin; WF segregation", fontsize=6.5, ha="right", va="bottom", color="#555555")
    ax.annotate("t", xy=(0.77, 0.73), xytext=(0.77, 0.18), arrowprops=dict(arrowstyle="<->", lw=1, color="#777777"), ha="center", va="center", fontsize=8)
    ax.text(0.77, 0.10, "structured interval", fontsize=7, ha="center", color="#555555")
    ax.add_patch(FancyArrowPatch((0.72, 0.55), (0.72, 0.38), arrowstyle="<->", mutation_scale=10, lw=1, color="#777777"))
    ax.text(0.82, 0.47, "cross-arrangement recombination", fontsize=6.5, va="center", color="#555555")
    ax.text(0.50, 0.04, "segregation and transmission at speciation; loci then generate genealogies", fontsize=7, ha="center", color="#555555")
    # outcome icons
    ax.text(0.99, 0.90, "possible local genealogies", fontsize=8, ha="right", fontweight="bold")
    for yy, topo, label in [(0.68, "12|34", "Q_S"), (0.43, "13|24", "Q_A"), (0.18, "14|23", "Q_O")]:
        tree_icon(ax, 0.80, yy, topo, scale=0.22, color=TOPOLOGY_COLORS[label], lw=1.0)
        ax.text(0.98, yy, f"{topo}  {label}", ha="right", va="center", fontsize=7.5, color=TOPOLOGY_COLORS[label], fontweight="bold" if label == "Q_A" else "normal")
    ax.text(0.99, 0.03, "mechanistic primary analysis averages over all terminal histories", ha="right", fontsize=6.5, color="#555555")


def panel_b(ax, rows, controls):
    selected = [0.1, 0.5, 2.0]
    for tau in selected:
        sub = [r for r in rows if float(r["internal_branch_coalescent_units"]) == tau and float(r["origin_fraction_of_root_extension"]) == 0.25 and float(r["initial_copy_fraction"]) == 0.20 and float(r["Ne"]) == 20]
        if not sub: continue
        # use the two cross-fraction points; x is suppression strength 1-fraction
        x = np.array([1-float(r["cross_arrangement_fraction"]) for r in sorted(sub, key=lambda z: float(z["cross_arrangement_fraction"]))])
        for key, col, lab in [("q_species", "Q_S", "Q_S"), ("q_alt", "Q_A", "Q_A"), ("q_other", "Q_O", "Q_O")]:
            y = np.array([float(r[key]) for r in sorted(sub, key=lambda z: float(z["cross_arrangement_fraction"]))])
            lo = np.array([float(r[key+"_ci95_low"]) for r in sorted(sub, key=lambda z: float(z["cross_arrangement_fraction"]))])
            hi = np.array([float(r[key+"_ci95_high"]) for r in sorted(sub, key=lambda z: float(z["cross_arrangement_fraction"]))])
            ax.errorbar(x, y, yerr=[y-lo, hi-y], fmt="o-", color=TOPOLOGY_COLORS[col], label=lab, capsize=2, lw=1.2, ms=4)
        ax.set_title(f"τ = {tau:g}", fontsize=8)
        ax.set_ylim(-0.02, 1.02); ax.axhline(1/3, color="#999999", ls=":", lw=0.8)
        ax.set_xticks([0, 0.9]); ax.set_xticklabels(["0", "0.9"])
        ax.set_xlabel("suppression strength 1 − cross fraction", fontsize=7); ax.set_ylabel("probability", fontsize=7)
        ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.20))


def panel_b_single(ax, rows, controls, tau):
    sub = [r for r in rows if float(r["internal_branch_coalescent_units"]) == tau and float(r["origin_fraction_of_root_extension"]) == 0.25 and float(r["initial_copy_fraction"]) == 0.20 and float(r["Ne"]) == 20]
    sub = sorted(sub, key=lambda z: float(z["cross_arrangement_fraction"]))
    x = np.array([1-float(r["cross_arrangement_fraction"]) for r in sub])
    for key, col, lab in [("q_species", "Q_S", "Q_S"), ("q_alt", "Q_A", "Q_A"), ("q_other", "Q_O", "Q_O")]:
        y = np.array([float(r[key]) for r in sub]); lo = np.array([float(r[key+"_ci95_low"]) for r in sub]); hi = np.array([float(r[key+"_ci95_high"]) for r in sub])
        ax.errorbar(x, y, yerr=[y-lo, hi-y], fmt="o-", color=TOPOLOGY_COLORS[col], label=lab, capsize=2, lw=1.2, ms=4)
    control = next((r for r in controls if float(r["Ne"]) == 20 and float(r["internal_branch_coalescent_units"]) == 2*tau), None)
    if control:
        ax.axhline(float(control["expected_q_species"]), color=TOPOLOGY_COLORS["Q_S"], ls=":", lw=.8)
        ax.axhline(float(control["expected_q_alt"]), color=TOPOLOGY_COLORS["Q_A"], ls=":", lw=.8)
    ax.set_title(f"τ = {tau:g}", fontsize=8); ax.set_ylim(-.02, 1.02); ax.axhline(1/3, color="#999999", ls=":", lw=.8)
    ax.set_xticks([0, .9]); ax.set_xticklabels(["0", ".9"]); ax.set_xlabel("suppression strength", fontsize=7); ax.set_ylabel("probability", fontsize=7); ax.tick_params(labelsize=7)
    ax.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(.5, 1.18))


def panel_c(ax, rows):
    # Show pilot observations in a quartet simplex, grouped by tau.
    for tau, marker in [(0.1, "o"), (0.5, "s"), (2.0, "^")]:
        sub = [r for r in rows if float(r["internal_branch_coalescent_units"]) == tau]
        x = np.array([float(r["q_alt"]) for r in sub]); y = np.array([float(r["q_other"]) for r in sub])
        c = np.array([float(r["cross_arrangement_fraction"]) for r in sub])
        sc = ax.scatter(x, y, c=c, cmap="viridis", vmin=0, vmax=1, s=18, marker=marker, alpha=.65, label=f"τ={tau:g}")
    ax.plot([0, 1], [0, 0], color="#bbbbbb", lw=.7); ax.plot([0, 0], [0, 1], color="#bbbbbb", lw=.7); ax.plot([0, 1], [1, 0], color="#bbbbbb", lw=.7)
    ax.plot([1/3], [1/3], marker="*", color="black", ms=9, label="MSC center")
    ax.set_xlim(-.03, 1.03); ax.set_ylim(-.03, 1.03); ax.set_xlabel("Q_A = 13|24", fontsize=7); ax.set_ylabel("Q_O = 14|23", fontsize=7)
    ax.set_title("Pilot quartet distributions", fontsize=8); ax.tick_params(labelsize=7); ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.text(.04, .94, "Q_S = 1 − Q_A − Q_O", transform=ax.transAxes, fontsize=6.5, color="#555555")


def panel_d(ax, rows):
    sub = [r for r in rows if float(r["Ne"]) == 20 and float(r["origin_fraction_of_root_extension"]) == .25 and float(r["initial_copy_fraction"]) == .20]
    taus = sorted(set(float(r["internal_branch_coalescent_units"]) for r in sub)); cross = sorted(set(float(r["cross_arrangement_fraction"]) for r in sub))
    z = np.full((len(cross), len(taus)), np.nan); sig = []
    for i, cf in enumerate(cross):
        for j, tau in enumerate(taus):
            r = next((r for r in sub if float(r["cross_arrangement_fraction"]) == cf and float(r["internal_branch_coalescent_units"]) == tau), None)
            if r:
                z[i, j] = float(r["delta_alt_species"])
                sig.append((j, i, float(r["delta_ci95_low"]) > 0))
    norm = TwoSlopeNorm(vmin=np.nanmin(z), vcenter=0, vmax=np.nanmax(z)) if np.nanmin(z) < 0 < np.nanmax(z) else None
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", norm=norm, vmin=None if norm else np.nanmin(z), vmax=None if norm else np.nanmax(z))
    ax.set_xticks(range(len(taus)), [f"{v:g}" for v in taus]); ax.set_yticks(range(len(cross)), [f"{v:g}" for v in cross]); ax.set_xlabel("structured branch duration τ", fontsize=7); ax.set_ylabel("cross-arrangement fraction", fontsize=7); ax.set_title("Unconditional quartet bias across a pilot slice", fontsize=8); ax.tick_params(labelsize=7)
    for j, i, supported in sig:
        if supported: ax.plot(j, i, marker="*", ms=10, color="black")
    cb = ax.figure.colorbar(im, ax=ax, fraction=.046, pad=.04); cb.set_label("Δ = P(Q_A) − P(Q_S)", fontsize=7); cb.ax.tick_params(labelsize=6)
    ax.text(.02, -.28, "Selected pilot slice: Ne=20, origin=.25, copy=.20; stars = lower 95% CI > 0", transform=ax.transAxes, fontsize=6.5, color="#555555")


def panel_e(ax, rows, patterns):
    """Compare unconditional Q_A with Q_A conditional on terminal 1010."""
    available = []
    for r in patterns:
        if r["terminal_pattern"] != "1010":
            continue
        p = next((x for x in rows if x["cell_id"] == r["cell_id"]), None)
        if p:
            available.append((p, r))
    available.sort(key=lambda pair: (float(pair[0]["internal_branch_coalescent_units"]), float(pair[0]["terminal_1010_frequency"])))
    if not available:
        ax.text(.5, .5, "No terminal 1010 histories observed", ha="center", va="center")
        ax.axis("off")
        return
    x = np.arange(len(available))
    freq = np.array([float(p["terminal_1010_frequency"]) for p, _ in available])
    unconditional = np.array([float(p["q_alt"]) for p, _ in available])
    conditional = np.array([float(r["q_alt"]) for _, r in available])
    ax2 = ax.twinx()
    ax2.bar(x, freq, color="#bdbdbd", alpha=.45, width=.65, label="frequency of 1010")
    ax.plot(x, unconditional, "o-", color=TOPOLOGY_COLORS["Q_A"], label="unconditional Q_A")
    ax.plot(x, conditional, "D--", color="#7b3294", label="Q_A | terminal 1010")
    ax.set_ylim(0, 1.02); ax2.set_ylim(0, max(.1, freq.max() * 1.35)); ax.set_ylabel("quartet probability", fontsize=7); ax2.set_ylabel("frequency of terminal 1010", fontsize=7)
    ax.set_xlabel("pilot cells with observed 1010", fontsize=7); ax.set_xticks(x); ax.set_xticklabels([f"τ={p['internal_branch_coalescent_units']}\ncf={p['cross_arrangement_fraction']}" for p, _ in available], fontsize=6)
    ax.set_title("Conditional versus unconditional behavior", fontsize=8); ax.tick_params(labelsize=7); ax2.tick_params(labelsize=7)
    handles, labels = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(handles+h2, labels+l2, fontsize=6.5, frameon=False, loc="upper left")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5]); args = ap.parse_args(); root=args.repo_root.resolve()
    out = root / "simulations/quartet/mechanistic/analysis/figures"; out.mkdir(parents=True, exist_ok=True)
    rows=read_tsv(root/"simulations/quartet/mechanistic/datasets/processed/mechanistic_quartet_pilot.tsv")
    patterns=read_tsv(root/"simulations/quartet/mechanistic/datasets/processed/mechanistic_pattern_contributions.tsv")
    controls=read_tsv(root/"simulations/quartet/mechanistic/datasets/processed/msc_control_validation.tsv")
    fig, ax = plt.subplots(figsize=(13, 5.5)); panel_a(ax); ax.text(-.02, 1.03, "A", transform=ax.transAxes, fontsize=14, fontweight="bold"); fig.savefig(out/"mechanistic_scenario.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_scenario.pdf", bbox_inches="tight"); plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True)
    for ax, tau in zip(axs, [0.1, 0.5, 2.0]): panel_b_single(ax, rows, controls, tau)
    fig.text(.01, .98, "B", fontsize=14, fontweight="bold", va="top"); fig.tight_layout(); fig.savefig(out/"mechanistic_unconditional_probabilities.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_unconditional_probabilities.pdf", bbox_inches="tight"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 5)); panel_c(ax, rows); fig.savefig(out/"mechanistic_quartet_simplex.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_quartet_simplex.pdf", bbox_inches="tight"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4.5)); panel_d(ax, rows); ax.text(-.08, 1.04, "C", transform=ax.transAxes, fontsize=14, fontweight="bold"); fig.savefig(out/"mechanistic_quartet_bias_map.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_quartet_bias_map.pdf", bbox_inches="tight"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5)); panel_e(ax, rows, patterns); ax.text(-.08, 1.04, "D", transform=ax.transAxes, fontsize=14, fontweight="bold"); fig.savefig(out/"mechanistic_conditional_vs_unconditional.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_conditional_vs_unconditional.pdf", bbox_inches="tight"); plt.close(fig)
    fig=plt.figure(figsize=(13, 11)); gs=fig.add_gridspec(2,2, width_ratios=[1.25,1], height_ratios=[1.15,1], hspace=.38, wspace=.28)
    axs=[fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1])]
    panel_a(axs[0]); panel_b_single(axs[1], rows, controls, 2.0); panel_d(axs[2], rows); panel_e(axs[3], rows, patterns)
    for label, ax in zip("ABCD", axs): ax.text(-.06, 1.04, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")
    fig.savefig(out/"mechanistic_quartet_pilot.png", dpi=300, bbox_inches="tight"); fig.savefig(out/"mechanistic_quartet_pilot.pdf", bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__": main()
