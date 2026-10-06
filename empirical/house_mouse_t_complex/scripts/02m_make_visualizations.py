#!/usr/bin/env python3
"""Generate tree, quartet-track, recombination-track, and manuscript figures."""

from __future__ import annotations

import io
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from PIL import Image

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    FIGURES,
    RESULTS,
    canonical_split,
    read_tsv,
    tree_bipartitions,
    write_tsv,
)


Q_COLORS = {"q_species": "#4c78a8", "q_t_alt": "#f58518", "q_other": "#54a24b"}
RECOMB_COLORS = {
    "fraction_no_recent_recombination": "#9e9e9e",
    "fraction_recent_or_older": "#756bb1",
    "fraction_very_recent_or_extensive": "#e34a33",
}
PRETTY = {
    "Mus_musculus_domesticus": "M. m. domesticus",
    "Mus_musculus_musculus": "M. m. musculus",
    "Mus_musculus_castaneus": "M. m. castaneus",
    "Mus_spretus": "M. spretus",
    "France": "France",
    "Germany": "Germany",
    "Afghanistan": "Afghanistan",
    "Czech_Republic": "Czech Republic",
    "Kazakhstan": "Kazakhstan",
    "CAST": "CAST",
    "SPRE": "SPRE",
}


def save(fig: plt.Figure, stem: str) -> None:
    png = FIGURES / f"{stem}.png"
    pdf = FIGURES / f"{stem}.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def load_tree(path: Path, outgroup: str) -> object:
    tree = Phylo.read(io.StringIO(path.read_text()), "newick")
    for terminal in tree.get_terminals():
        terminal.name = PRETTY.get(terminal.name, terminal.name)
    for clade in tree.get_nonterminals():
        clade.name = None
    out = next(t for t in tree.get_terminals() if t.name == PRETTY[outgroup])
    tree.root_with_outgroup(out)
    return tree


def tree_label(clade) -> str | None:
    if clade.is_terminal():
        return clade.name
    return None


def draw_tree(ax, tree, title: str, annotation: str, color: str) -> None:
    Phylo.draw(tree, axes=ax, do_show=False, label_func=tree_label, show_confidence=False)
    ax.set_title(title, color=color, fontsize=12, weight="bold")
    ax.text(0.02, 0.02, annotation, transform=ax.transAxes, va="bottom", ha="left", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.9, "boxstyle": "round,pad=0.35"})
    ax.set_xlabel("ASTRAL branch length (CU; visual root only)")


def make_4group_trees() -> None:
    paths = {
        "STANDARD": RESULTS / "stage1_aster" / "T0_STANDARD_subspecies.nwk",
        "ALL": RESULTS / "stage1_aster" / "T1_ALL_WINDOWS_subspecies.nwk",
    }
    annotations = {
        "STANDARD": "focal split: Q_SPECIES\nCU=0.056283; localPP=0.999985\nq_species=0.369910; q_t_alt=0.343656; q_other=0.286434",
        "ALL": "focal split: Q_T_ALT\nCU=0.049097; localPP=1.000000\nq_species=0.365364; q_t_alt=0.323323; q_other=0.311313",
    }
    colors = {"STANDARD": Q_COLORS["q_species"], "ALL": Q_COLORS["q_t_alt"]}
    trees = {key: load_tree(path, "Mus_spretus") for key, path in paths.items()}
    max_x = max(max((clade.branch_length or 0) for clade in tree.find_clades()) for tree in trees.values()) * 1.25
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), sharex=True)
    for ax, key in zip(axes, ("STANDARD", "ALL")):
        draw_tree(ax, trees[key], "STANDARD / noncarrier" if key == "STANDARD" else "ALL standard + pseudo-t", annotations[key], colors[key])
        ax.set_xlim(0, max_x)
    fig.suptitle("Four-group ASTRAL trees; inference is unrooted, SPRE is a visual root", fontsize=14)
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_astral_4group")


def parse_split_id(text: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    left, right = text.split(" || ")
    # Population split tables use source identifiers while plotted trees use
    # publication-facing labels.  Normalize both sides before comparison.
    normalize = lambda token: PRETTY.get(token, token)
    return canonical_split([[normalize(x) for x in left.split("|")],
                            [normalize(x) for x in right.split("|")]])


def set_changed_colors(tree, lost: set, gained: set) -> None:
    all_tips = frozenset(t.name for t in tree.get_terminals())
    from Bio.Phylo.BaseTree import BranchColor
    for clade in tree.find_clades(order="postorder"):
        side = frozenset(t.name for t in clade.get_terminals())
        if not (1 < len(side) < len(all_tips) - 1):
            continue
        split = canonical_split([side, all_tips - side])
        if split in lost:
            clade.color = BranchColor.from_hex("#d62728")
        elif split in gained:
            clade.color = BranchColor.from_hex("#1f77b4")


def make_7population_trees() -> None:
    comparison = read_tsv(RESULTS / "stage2_population_split_comparison.tsv")
    lost = {parse_split_id(row["split_id"]) for row in comparison if row["status"] == "lost"}
    gained = {parse_split_id(row["split_id"]) for row in comparison if row["status"] == "gained"}
    p0 = load_tree(RESULTS / "stage2_population_astral" / "P0_STANDARD.nwk", "SPRE")
    p1 = load_tree(RESULTS / "stage2_population_astral" / "P1_ALL.nwk", "SPRE")
    set_changed_colors(p0, lost, set())
    set_changed_colors(p1, set(), gained)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True)
    draw_tree(axes[0], p0, "P0_STANDARD", "red = split lost in P1\nRF(P0,P1) = 2", "#444444")
    draw_tree(axes[1], p1, "P1_ALL", "blue = split gained in P1\nRF(P0,P1) = 2", "#444444")
    fig.suptitle("Seven-population ASTRAL trees; SPRE is a visual root", fontsize=14)
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_astral_7population")


def bin_q(rows: list[dict[str, str]], treatment: str, bin_size: int = 250_000) -> list[dict[str, float]]:
    subset = [row for row in rows if row["treatment"] == treatment]
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in subset:
        start = int(row["start_bp"])
        b = 5_000_000 + ((start - 5_000_000) // bin_size) * bin_size
        grouped[b].append(row)
    out = []
    for start, block in sorted(grouped.items()):
        out.append({"midpoint": start + bin_size / 2, **{key: float(np.mean([float(row[key]) for row in block])) for key in ("q_species", "q_t_alt", "q_other")}})
    return out


def plot_segments(ax, x: np.ndarray, y: np.ndarray, **kwargs) -> None:
    """Plot contiguous coordinate runs without bridging genomic gaps."""
    if len(x) == 0:
        return
    label = kwargs.pop("label", None)
    breaks = np.flatnonzero(np.diff(x) > 0.375) + 1
    for index, segment in enumerate(np.split(np.arange(len(x)), breaks)):
        if len(segment) >= 2:
            ax.plot(x[segment], y[segment], label=label if index == 0 else None, **kwargs)


def plot_q_axis(ax, rows: list[dict[str, str]], treatment: str, title: str, show_delta: bool = True) -> None:
    raw = [row for row in rows if row["treatment"] == treatment]
    bins = bin_q(rows, treatment)
    xraw = np.array([float(row["midpoint_bp"]) / 1e6 for row in raw])
    for key, color in Q_COLORS.items():
        ax.scatter(xraw, [float(row[key]) for row in raw], s=2, alpha=0.08, color=color)
        xb = np.array([row["midpoint"] / 1e6 for row in bins])
        yb = np.array([row[key] for row in bins])
        plot_segments(ax, xb, yb, color=color, lw=1.7, label=key.replace("q_", ""))
    ax.axhline(1 / 3, color="#777777", ls="--", lw=0.7)
    ax.set_ylim(0, 1)
    ax.set_ylabel("quartet fraction")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, ncol=3)


def make_q_tracks() -> None:
    rows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    plot_q_axis(axes[0], rows, "STANDARD_ONLY", "STANDARD_ONLY; 250-kb mean, raw windows light")
    plot_q_axis(axes[1], rows, "ALL_TIPS", "ALL_TIPS; 250-kb mean, raw windows light")
    axes[1].set_xlabel("chr17 position (Mb)")
    fig.suptitle("Direct local quartet support along chr17", fontsize=14)
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_q_tracks")


def load_recomb() -> list[dict[str, str]]:
    rows = read_tsv(RESULTS / "stage2_recombination_state_500kb.tsv")
    for row in rows:
        total = sum(float(row[key]) for key in RECOMB_COLORS)
        assert abs(total - 1.0) < 1e-8, f"recombination fractions do not sum to 1: {row}"
    return rows


def plot_recomb_axes(axes, rows: list[dict[str, str]], title: str | None = None) -> None:
    for ax, species in zip(axes, ("domesticus", "musculus", "castaneus")):
        sub = [row for row in rows if row["subspecies"] == species]
        x = np.array([(int(row["bin_start"]) + 250_000) / 1e6 for row in sub])
        bottom = np.zeros(len(sub))
        for key, color in RECOMB_COLORS.items():
            values = np.array([float(row[key]) for row in sub])
            ax.bar(x, values, width=0.49, bottom=bottom, color=color, label=key.replace("fraction_", "").replace("_", " "))
            bottom += values
        ax.set_ylim(0, 1)
        ax.set_ylabel(species)
        ax.set_yticks([0, 1])
        ax.grid(axis="y", alpha=0.2)
    axes[-1].set_xlabel("chr17 position (Mb)")
    if title:
        axes[0].set_title(title)
    axes[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")


def make_recombination_track() -> None:
    rows = load_recomb()
    fig, axes = plt.subplots(3, 1, figsize=(12, 6), sharex=True)
    plot_recomb_axes(axes, rows, "Phylogeny-based recombination-state track; 500-kb bins")
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_recombination_track")


def make_q_recombination() -> None:
    qrows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    recomb = load_recomb()
    fig = plt.figure(figsize=(13, 11))
    grid = fig.add_gridspec(5, 1, height_ratios=[1.2, 0.55, 0.7, 0.7, 0.7], hspace=0.25)
    axq = fig.add_subplot(grid[0])
    plot_q_axis(axq, qrows, "ALL_TIPS", "ALL_TIPS quartet support; 250-kb means")
    axd = fig.add_subplot(grid[1], sharex=axq)
    bins = bin_q(qrows, "ALL_TIPS")
    xd = np.array([x["midpoint"] / 1e6 for x in bins])
    yd = np.array([x["q_species"] - x["q_t_alt"] for x in bins])
    plot_segments(axd, xd, yd, color="#333333", lw=1.5)
    axd.axhline(0, color="black", lw=0.7)
    axd.set_ylabel("delta")
    axd.set_title("q_species - q_t_alt")
    axes = [fig.add_subplot(grid[i], sharex=axq) for i in (2, 3, 4)]
    plot_recomb_axes(axes, recomb)
    axes[-1].set_xlabel("chr17 position (Mb)")
    fig.suptitle("Quartet support aligned with source-defined phylogeny-based recombination states", fontsize=14)
    save(fig, "house_mouse_t_complex_q_recombination")


def arrangement_panel(ax) -> None:
    rows = read_tsv(RESULTS / "stage2_arrangement_pattern_summary.tsv")
    patterns = ["SSS", "one_T", "two_T", "TTT"]
    data = {row["status_pattern"]: row for row in rows}
    bottom = np.zeros(4)
    for key, label, color in [("q_species", "Q_SPECIES", Q_COLORS["q_species"]), ("q_t_alt", "Q_T_ALT", Q_COLORS["q_t_alt"]), ("q_other", "Q_OTHER", Q_COLORS["q_other"]), ("q_unresolved", "unresolved", "#bbbbbb")]:
        values = np.array([float(data[p][key]) for p in patterns])
        ax.bar(patterns, values, bottom=bottom, label=label, color=color)
        bottom += values
    ax.set_ylim(0, 1)
    ax.set_title("D. Arrangement-state quartet composition")
    ax.legend(frameon=False, fontsize=8)


def make_controls() -> None:
    balanced = read_tsv(RESULTS / "stage2_balanced_quartet_resampling_summary.tsv")
    linkage = read_tsv(RESULTS / "stage2_linkage_summary.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bmap = {row["treatment"]: row for row in balanced}
    treatments = ["B0_STANDARD_MATCHED", "B1_T_MATCHED", "B2_MIXED_BALANCED"]
    bottom = np.zeros(3)
    for key, label, color in [("fraction_Q_SPECIES", "Q_SPECIES", Q_COLORS["q_species"]), ("fraction_Q_T_ALT", "Q_T_ALT", Q_COLORS["q_t_alt"]), ("fraction_Q_OTHER", "Q_OTHER", Q_COLORS["q_other"])]:
        values = np.array([float(bmap[t][key]) for t in treatments])
        axes[0].bar(treatments, values, bottom=bottom, label=label, color=color)
        bottom += values
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Balanced exact resampling (n=1000)")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].legend(frameon=False)
    x = [int(row["spacing_kb"]) for row in linkage]
    axes[1].plot(x, [float(row["fraction_Q_T_ALT"]) for row in linkage], marker="o", color=Q_COLORS["q_t_alt"], label="Q_T_ALT frequency")
    axes[1].plot(x, [float(row["median_localPP"]) for row in linkage], marker="s", color="#d62728", label="median localPP")
    axes[1].set_xscale("log")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel("spacing (kb)")
    axes[1].set_title("Spatial thinning control")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_controls")


def make_main_v2() -> None:
    qrows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    recomb = load_recomb()
    fig = plt.figure(figsize=(15, 12))
    outer = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.28, wspace=0.24)
    tree_grid = outer[0, 0].subgridspec(1, 2, wspace=0.05)
    trees = {"STANDARD": load_tree(RESULTS / "stage1_aster" / "T0_STANDARD_subspecies.nwk", "Mus_spretus"), "ALL": load_tree(RESULTS / "stage1_aster" / "T1_ALL_WINDOWS_subspecies.nwk", "Mus_spretus")}
    axes_tree = [fig.add_subplot(tree_grid[i]) for i in range(2)]
    draw_tree(axes_tree[0], trees["STANDARD"], "STANDARD", "Q_SPECIES", Q_COLORS["q_species"])
    draw_tree(axes_tree[1], trees["ALL"], "ALL", "Q_T_ALT", Q_COLORS["q_t_alt"])
    axes_tree[0].set_title("A. STANDARD\nQ_SPECIES", color=Q_COLORS["q_species"], fontsize=10)
    axes_tree[1].set_title("ALL TIPS\nQ_T_ALT", color=Q_COLORS["q_t_alt"], fontsize=10)
    axq = fig.add_subplot(outer[0, 1])
    plot_q_axis(axq, qrows, "ALL_TIPS", "B. ALL_TIPS direct quartet support")
    axes_rec = [fig.add_subplot(outer[1, 0].subgridspec(3, 1, hspace=0.05)[i]) for i in range(3)]
    plot_recomb_axes(axes_rec, recomb, "C. Phylogeny-based recombination state")
    axarr = fig.add_subplot(outer[1, 1])
    arrangement_panel(axarr)
    fig.suptitle("House-mouse t-complex: linked genealogy distortion and structural-state context", fontsize=15)
    save(fig, "house_mouse_t_complex_main_v2")


def quartet_by_recombination() -> None:
    qrows = {row["locus_id"]: row for row in read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv") if row["treatment"] == "ALL_TIPS"}
    state_rows = read_tsv(RESULTS / "stage2_recombination_state_5kb.tsv")
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in state_rows:
        if row["locus_id"] in qrows:
            groups[(row["subspecies"], row["recombination_class"])].append(qrows[row["locus_id"]])
    out = []
    for (species, category), subset in sorted(groups.items()):
        out.append({
            "subspecies": species,
            "recombination_class": category,
            "n_windows": len(subset),
            "mean_q_species": np.mean([float(row["q_species"]) for row in subset]),
            "mean_q_t_alt": np.mean([float(row["q_t_alt"]) for row in subset]),
            "mean_q_other": np.mean([float(row["q_other"]) for row in subset]),
            "mean_delta_species_alt": np.mean([float(row["delta_species_alt"]) for row in subset]),
            "median_q_species": np.median([float(row["q_species"]) for row in subset]),
            "median_q_t_alt": np.median([float(row["q_t_alt"]) for row in subset]),
            "median_q_other": np.median([float(row["q_other"]) for row in subset]),
            "median_delta_species_alt": np.median([float(row["delta_species_alt"]) for row in subset]),
        })
    write_tsv(RESULTS / "stage2_quartet_by_recombination_state.tsv", out, list(out[0].keys()))


def validate_figures() -> None:
    stems = [
        "house_mouse_t_complex_astral_4group",
        "house_mouse_t_complex_astral_7population",
        "house_mouse_t_complex_q_tracks",
        "house_mouse_t_complex_recombination_track",
        "house_mouse_t_complex_q_recombination",
        "house_mouse_t_complex_controls",
        "house_mouse_t_complex_main_v2",
    ]
    validation = []
    for stem in stems:
        png = FIGURES / f"{stem}.png"
        pdf = FIGURES / f"{stem}.pdf"
        assert png.stat().st_size > 20_000, f"small PNG: {png}"
        assert pdf.stat().st_size > 10_000, f"small PDF: {pdf}"
        with Image.open(png) as image:
            rgb = image.convert("RGB")
            pixels = np.asarray(rgb)
            nonwhite = int(np.sum(np.any(pixels < 245, axis=2)))
            distinct = int(np.unique(pixels.reshape(-1, 3), axis=0).shape[0])
            assert pixels.shape[1] >= 800 and pixels.shape[0] >= 500
            assert nonwhite > 5_000, f"mostly blank PNG: {png}"
            assert distinct > 100, f"low-variation PNG: {png}"
            validation.append({"figure": stem, "png_bytes": png.stat().st_size,
                               "pdf_bytes": pdf.stat().st_size,
                               "width_px": pixels.shape[1], "height_px": pixels.shape[0],
                               "nonwhite_pixels": nonwhite, "distinct_rgb": distinct})
    write_tsv(RESULTS / "stage2_visualization_validation.tsv", validation,
              list(validation[0].keys()))


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    qrows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    for row in qrows:
        total = sum(float(row[key]) for key in ("q_species", "q_t_alt", "q_other", "q_unresolved"))
        assert abs(total - 1.0) < 1e-8, f"quartet fractions do not sum to 1: {row}"
    make_4group_trees()
    make_7population_trees()
    make_q_tracks()
    make_recombination_track()
    make_q_recombination()
    make_controls()
    quartet_by_recombination()
    make_main_v2()
    validate_figures()
    print("Wrote visualization suite")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
