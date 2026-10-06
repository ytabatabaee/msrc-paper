#!/usr/bin/env python3
"""Create manuscript figures from the frozen processed quartet table."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[4])
    args = ap.parse_args()
    root = args.repo_root.resolve()
    data = load(root / "simulations/quartet/datasets/processed/conditional_quartet_grid.tsv")
    figdir = root / "simulations/quartet/analysis/figures"
    figdir.mkdir(parents=True, exist_ok=True)
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
    fig, axes = plt.subplots(2, 4, figsize=(12, 6), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for ax, duration in zip(axes.flat, durations_all):
        subset = sorted((r for r in data if float(r["duration"]) == duration), key=lambda r: float(r["m"]))
        xs = [float(r["exact_q_alt"]) + float(r["exact_q_other"]) / 2 for r in subset]
        ys = [np.sqrt(3) * float(r["exact_q_other"]) / 2 for r in subset]
        ax.plot(xs, ys, color="#888888", lw=0.8)
        ax.scatter(xs, ys, c=[float(r["m"]) for r in subset], cmap=cmap, s=24)
        ax.scatter([0.5], [np.sqrt(3)/6], color="black", marker="+", s=55)
        ax.set_title(f"duration = {duration:g}")
        ax.set_xlim(0, 1); ax.set_ylim(0, np.sqrt(3)/2); ax.set_aspect("equal"); ax.axis("off")
        ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3)/2, 0], color="black", lw=0.8)
    axes.flat[-1].axis("off")
    axes.flat[0].text(0, -0.08, "Q_SPECIES", transform=axes.flat[0].transAxes, ha="left")
    axes.flat[0].text(1, -0.08, "Q_ALT", transform=axes.flat[0].transAxes, ha="right")
    axes.flat[0].text(0.5, 1.02, "Q_OTHER", transform=axes.flat[0].transAxes, ha="center")
    fig.suptitle("Conditional quartet probability simplex")
    fig.savefig(figdir / "conditional_quartet_simplex.png", dpi=300)
    fig.savefig(figdir / "conditional_quartet_simplex.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
