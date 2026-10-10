#!/usr/bin/env python3
"""Generate diagnostic figures for the completed Experiment 2B pilot."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


COLORS = {0.1: "#d95f02", 0.5: "#1b9e77", 2.0: "#7570b3"}


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, out, name):
    fig.savefig(out / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_a(axs, persistence):
    for ax, tau in zip(axs, [0.5, 2.0]):
        for cross, ls in [(0.1, "-"), (1.0, "--")]:
            for origin, marker in [(0.25, "o"), (0.75, "s")]:
                sub = [r for r in persistence if float(r["tau"]) == tau and float(r["cross_fraction"]) == cross and float(r["origin_fraction"]) == origin and r["metric"] == "daughter_AB_persistence"]
                sub = sorted(sub, key=lambda r: float(r["initial_copy_fraction"]))
                if not sub: continue
                x = [float(r["initial_copy_fraction"]) for r in sub]; y = [float(r["estimate"]) for r in sub]
                lo = [float(r["ci95_low"]) for r in sub]; hi = [float(r["ci95_high"]) for r in sub]
                ax.errorbar(x, y, yerr=[np.array(y)-lo, np.array(hi)-y], color=COLORS[tau], ls=ls, marker=marker, capsize=2, label=f"origin={origin}, cross={cross}")
        ax.set_title(f"τ = {tau:g}"); ax.set_xlabel("initial rearrangement frequency"); ax.set_ylabel("P(persistent in A and B)"); ax.set_ylim(-.03, 1.03); ax.grid(alpha=.2)
    axs[1].legend(fontsize=6.5, frameon=False, ncol=2)


def figure_b(ax, decomposition, summaries):
    choices = ["ne20_tau0p5_orig0p75_copy0p5_cross0p1", "ne100_tau2_orig0p75_copy0p5_cross0p1"]
    tick_x = []; tick_labels = []
    for i, cell in enumerate(choices):
        sub = [r for r in decomposition if r["cell_id"] == cell]
        if not sub: continue
        cats = [r["category"] for r in sub]; vals = [float(r["contribution"]) for r in sub]; lo = [float(r["contribution_bootstrap_ci95_low"]) for r in sub]; hi = [float(r["contribution_bootstrap_ci95_high"]) for r in sub]
        x = np.arange(len(cats)) + i * (len(cats) + 1)
        tick_x.extend(x.tolist()); tick_labels.extend([c.replace("terminal_", "").replace("_", "\n") for c in cats])
        ax.bar(x, vals, color=["#b2182b" if v < 0 else "#2166ac" for v in vals], alpha=.8)
        ax.errorbar(x, vals, yerr=[np.array(vals)-lo, np.array(hi)-vals], fmt="none", ecolor="black", capsize=2, lw=.8)
    ax.axhline(0, color="black", lw=.8); ax.set_ylabel("category contribution C_H = P(H)Δ_H"); ax.set_title("Bias decomposition: contributions sum to unconditional Δ\nLeft: Ne=20, τ=.5, origin=.75, copy=.5, cross=.1 | Right: Ne=100, τ=2, origin=.75, copy=.5, cross=.1", fontsize=10, pad=12)
    ax.set_xticks(tick_x); ax.set_xticklabels(tick_labels, fontsize=5.5, rotation=45, ha="right")
    ax.text(.01, -.18, "Red = negative; blue = positive. Empty categories contribute zero and have no conditional estimate.", transform=ax.transAxes, fontsize=6.5, color="#555555")


def figure_c(ax, summaries, decomposition):
    points = []
    for d in decomposition:
        if d["category"] != "terminal_1010" or int(d["n_histories"]) < 3: continue
        s = next(r for r in summaries if r["cell_id"] == d["cell_id"])
        points.append((d, s))
    x = np.arange(len(points)); width = .25
    q_a_cond = [float(d["conditional_q_alt"]) for d, _ in points]; q_a = [float(s["q_alt"]) for _, s in points]; q_s = [float(s["q_species"]) for _, s in points]
    ax.bar(x-width, q_a_cond, width, label="Q_A | terminal 1010", color="#7b3294")
    ax.bar(x, q_a, width, label="unconditional Q_A", color="#d95f02")
    ax.bar(x+width, q_s, width, label="unconditional Q_S", color="#2c7fb8")
    ax.set_xticks(x, [f"τ={s['tau']}\ncf={s['cross_fraction']}\nn={d['n_histories']}" for d, s in points], fontsize=6); ax.set_ylim(0, 1.05); ax.set_ylabel("quartet probability"); ax.set_title("Conditional versus unconditional bias (supported 1010 cells)"); ax.legend(fontsize=6.5, frameon=False)
    ax.text(.01, -.18, "Only cells with ≥3 independent 1010 histories are shown; the primary pilot remains unconditional.", transform=ax.transAxes, fontsize=6.5, color="#555555")


def figure_d(ax, comparisons):
    for tau in sorted(set(float(r["tau"]) for r in comparisons)):
        sub = [r for r in comparisons if float(r["tau"]) == tau]
        x = np.array([(float(r["suppressed_persistence"]) + float(r["unsuppressed_persistence"])) / 2 for r in sub]); y = np.array([float(r["q_species_difference_suppressed_minus_unsuppressed"]) for r in sub]); lo = np.array([float(r["q_species_difference_ci95_low"]) for r in sub]); hi = np.array([float(r["q_species_difference_ci95_high"]) for r in sub])
        ax.errorbar(x, y, yerr=[y-lo, hi-y], fmt="o", color=COLORS[tau], label=f"τ={tau:g}", alpha=.8, capsize=2)
    ax.axhline(0, color="#777777", ls=":"); ax.set_xlabel("mean daughter-branch persistence"); ax.set_ylabel("P_suppressed(Q_S) − P_unsuppressed(Q_S)"); ax.set_title("Suppression effect on species-topology support"); ax.legend(frameon=False, fontsize=7)


def figure_e(ax, summaries):
    for cross, marker, label in [(0.1, "o", "cross=0.1"), (1.0, "s", "cross=1.0")]:
        sub = [r for r in summaries if float(r["cross_fraction"]) == cross and not r["short_branch_control"]]
        x = np.array([float(r["q_alt"]) for r in sub]); y = np.array([float(r["q_other"]) for r in sub])
        ax.scatter(x, y, marker=marker, alpha=.55, label=label)
    q_o = np.linspace(0, 1, 300); boundary = (1-q_o)/2
    ax.plot(boundary, q_o, color="#555555", ls="--", lw=1, label="Q_A = Q_S")
    ax.fill_betweenx(q_o, boundary, 1-q_o, color="#fddbc7", alpha=.35, label="Q_A > Q_S")
    ax.plot([1/3], [1/3], "*", color="black", ms=9, label="symmetric point")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("Q_A = 13|24"); ax.set_ylabel("Q_O = 14|23"); ax.set_title("Unconditional mechanistic quartet simplex"); ax.legend(fontsize=6.5, frameon=False, loc="upper right")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5]); args = ap.parse_args(); root=args.repo_root.resolve(); out=root/"simulations/quartet/mechanistic/analysis/figures"; out.mkdir(parents=True, exist_ok=True)
    base=root/"simulations/quartet/mechanistic/datasets/processed"; persistence=read(base/"mechanistic_pilot2_persistence.tsv"); summaries=read(base/"mechanistic_pilot2_summary.tsv"); decomposition=read(base/"mechanistic_pilot2_bias_decomposition.tsv"); comparisons=read(base/"mechanistic_pilot2_suppression_comparison.tsv")
    fig, axs=plt.subplots(1,2,figsize=(9,4)); figure_a(axs,persistence); fig.text(.01,.98,"A",fontsize=14,fontweight="bold",va="top"); save(fig,out,"pilot2_persistence")
    fig, ax=plt.subplots(figsize=(9,4.5)); figure_b(ax,decomposition,summaries); fig.text(.01,.98,"B",fontsize=14,fontweight="bold",va="top"); save(fig,out,"pilot2_bias_decomposition")
    fig, ax=plt.subplots(figsize=(8,4.5)); figure_c(ax,summaries,decomposition); fig.text(.01,.98,"C",fontsize=14,fontweight="bold",va="top"); save(fig,out,"pilot2_conditional_unconditional")
    fig, ax=plt.subplots(figsize=(6,4.5)); figure_d(ax,comparisons); fig.text(.01,.98,"D",fontsize=14,fontweight="bold",va="top"); save(fig,out,"pilot2_suppression_effect")
    fig, ax=plt.subplots(figsize=(6,5)); figure_e(ax,summaries); fig.text(.01,.98,"E",fontsize=14,fontweight="bold",va="top"); save(fig,out,"pilot2_mechanistic_simplex")


if __name__ == "__main__": main()
