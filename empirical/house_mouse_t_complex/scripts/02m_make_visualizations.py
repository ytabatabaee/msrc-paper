#!/usr/bin/env python3
"""Generate publication figures from frozen House-mouse Stage-2B tables."""

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
    four_taxon_topology_from_newick,
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
    "France": "France", "Germany": "Germany", "Afghanistan": "Afghanistan",
    "Czech_Republic": "Czech Republic", "Kazakhstan": "Kazakhstan",
    "CAST": "CAST", "SPRE": "SPRE",
}
FOUR_LABELS = ["Mus_musculus_domesticus", "Mus_musculus_musculus", "Mus_musculus_castaneus", "Mus_spretus"]
POP_LABELS = ["France", "Germany", "Afghanistan", "Czech_Republic", "Kazakhstan", "CAST", "SPRE"]


def save(fig: plt.Figure, stem: str | list[str]) -> None:
    stems = [stem] if isinstance(stem, str) else stem
    for name in stems:
        fig.savefig(FIGURES / f"{name}.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
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


def parse_split_id(text: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    left, right = text.split(" || ")
    return canonical_split([[PRETTY.get(x, x) for x in left.split("|")], [PRETTY.get(x, x) for x in right.split("|")]])


def ordered_children(clade, rank: dict[str, int]):
    return sorted(clade.clades, key=lambda child: min(rank[t.name] for t in child.get_terminals()))


def cladogram_layout(tree, label_order: list[str]):
    rank = {label: i for i, label in enumerate(label_order)}
    terminals = list(tree.get_terminals())
    if set(t.name for t in terminals) != set(label_order):
        raise AssertionError(f"terminal labels differ from expected set: {[t.name for t in terminals]}")
    y = {label: float(rank[label]) for label in label_order}
    x = {}
    lines = []

    def visit(clade, depth: int):
        children = ordered_children(clade, rank) if clade.clades else []
        if not children:
            node_y = y[clade.name]
        else:
            child_ys = [visit(child, depth + 1) for child in children]
            node_y = float(np.mean(child_ys))
            for child, child_y in zip(children, child_ys):
                lines.append((depth, x[child], node_y, child_y, child))
            lines.append((depth, depth, min(child_ys), max(child_ys), None))
        x[clade] = float(depth)
        return node_y

    visit(tree.root, 0)
    return x, y, lines


def draw_cladogram(ax, tree, label_order: list[str], title: str, title_color: str = "#333333", focal_split=None, changed_splits=None):
    x, y, lines = cladogram_layout(tree, label_order)
    all_tips = frozenset(label_order)
    changed_splits = changed_splits or {}
    for depth, child_x, parent_y, child_y, child in lines:
        color, lw, ls = "#222222", 1.8, "-"
        if child is not None:
            side = frozenset(t.name for t in child.get_terminals())
            split = canonical_split([side, all_tips - side])
            if focal_split is not None and split == focal_split:
                color, lw = title_color, 3.0
            if split in changed_splits:
                color, lw, ls = changed_splits[split]
        if child is None:
            ax.plot([depth, depth], [parent_y, child_y], color=color, lw=lw, ls=ls, clip_on=True)
        else:
            ax.plot([child_x, depth], [child_y, child_y], color=color, lw=lw, ls=ls, clip_on=True)

    rank = {label: i for i, label in enumerate(label_order)}

    def y_for(clade):
        if not clade.clades:
            return y[clade.name]
        return float(np.mean([y_for(c) for c in ordered_children(clade, rank)]))

    def draw_verticals(clade):
        children = ordered_children(clade, rank) if clade.clades else []
        if children:
            values = [y_for(c) for c in children]
            ax.plot([x[clade], x[clade]], [min(values), max(values)], color="#222222", lw=1.8, clip_on=True)
            for child in children:
                draw_verticals(child)

    draw_verticals(tree.root)
    xmax = max(x.values())
    for label in label_order:
        ax.text(xmax + 0.12, y[label], label, va="center", ha="left", fontsize=10, clip_on=True)
    ax.set_xlim(-0.15, xmax + 3.4)
    ax.set_ylim(-0.8, len(label_order) - 0.2)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, color=title_color, fontsize=12, weight="bold", pad=8)


def focal_split(name: str):
    if name == "Q_SPECIES":
        return canonical_split([["M. m. musculus", "M. m. castaneus"], ["M. m. domesticus", "M. spretus"]])
    if name == "Q_T_ALT":
        return canonical_split([["M. m. domesticus", "M. m. musculus"], ["M. m. castaneus", "M. spretus"]])
    return canonical_split([["M. m. domesticus", "M. m. castaneus"], ["M. m. musculus", "M. spretus"]])


def focal_annotations() -> dict[str, dict[str, float | str]]:
    """Map ASTRAL q1/q2/q3 to biological names using the inferred split."""
    out = {}
    for row in read_tsv(RESULTS / "stage2_astral_summary.tsv"):
        if row["treatment"] not in {"T0_STANDARD", "T1_ALL_WINDOWS"}:
            continue
        key = "STANDARD" if row["treatment"] == "T0_STANDARD" else "ALL"
        inferred = row["topology"]
        if inferred == "Q_SPECIES":
            mapped = {"q_species": float(row["q1"]), "q_t_alt": float(row["q2"]), "q_other": float(row["q3"])}
        elif inferred == "Q_T_ALT":
            mapped = {"q_species": float(row["q2"]), "q_t_alt": float(row["q1"]), "q_other": float(row["q3"])}
        else:
            raise AssertionError(f"unexpected focal topology: {inferred}")
        mapped.update({"inferred": inferred, "CU": float(row["CU_length"]), "localPP": float(row["localPP"])})
        out[key] = mapped
    assert set(out) == {"STANDARD", "ALL"}
    return out


def make_4group_trees() -> None:
    paths = {"STANDARD": RESULTS / "stage1_aster" / "T0_STANDARD_subspecies.nwk", "ALL": RESULTS / "stage1_aster" / "T1_ALL_WINDOWS_subspecies.nwk"}
    annotations = focal_annotations()
    colors = {"STANDARD": Q_COLORS["q_species"], "ALL": Q_COLORS["q_t_alt"]}
    trees = {key: load_tree(path, "Mus_spretus") for key, path in paths.items()}
    labels = [PRETTY[x] for x in FOUR_LABELS]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.8), gridspec_kw={"wspace": 0.55})
    assert len(axes) == 2
    for ax, key, inferred in zip(axes, ("STANDARD", "ALL"), ("Q_SPECIES", "Q_T_ALT")):
        draw_cladogram(ax, trees[key], labels, "STANDARD / noncarrier" if key == "STANDARD" else "ALL standard + pseudo-t", colors[key], focal_split=focal_split(inferred))
        ann = annotations[key]
        text = (f"{ann['inferred']}\nCU = {ann['CU']:.6f}\nlocalPP = {ann['localPP']:.6f}\n"
                f"q_species = {ann['q_species']:.6f}\nq_t_alt = {ann['q_t_alt']:.6f}\nq_other = {ann['q_other']:.6f}")
        ax.text(0.02, -0.04, text, transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"facecolor": "white", "edgecolor": colors[key], "boxstyle": "round,pad=0.35"})
        assert all("q1=" not in item.get_text() and "q2=" not in item.get_text() and "q3=" not in item.get_text() for item in ax.texts)
    fig.suptitle("Four-group ASTRAL cladograms; unrooted inference, SPRE used only for orientation", fontsize=13)
    fig.subplots_adjust(top=0.84, bottom=0.24, left=0.04, right=0.98)
    save(fig, "house_mouse_t_complex_astral_4group")


def make_7population_trees() -> None:
    comparison = read_tsv(RESULTS / "stage2_population_split_comparison.tsv")
    lost = {parse_split_id(row["split_id"]) for row in comparison if row["status"] == "lost"}
    gained = {parse_split_id(row["split_id"]) for row in comparison if row["status"] == "gained"}
    p0 = load_tree(RESULTS / "stage2_population_astral" / "P0_STANDARD.nwk", "SPRE")
    p1 = load_tree(RESULTS / "stage2_population_astral" / "P1_ALL.nwk", "SPRE")
    labels = [PRETTY[x] for x in POP_LABELS]
    changed_p0 = {split: ("#d62728", 3.0, "-") for split in lost}
    changed_p1 = {split: ("#1f77b4", 3.0, "-") for split in gained}
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw={"wspace": 0.62})
    assert len(axes) == 2
    draw_cladogram(axes[0], p0, labels, "P0_STANDARD", changed_splits=changed_p0)
    draw_cladogram(axes[1], p1, labels, "P1_ALL", changed_splits=changed_p1)
    axes[0].text(0.02, -0.04, "red = split lost in P1", transform=axes[0].transAxes, fontsize=9, color="#d62728")
    axes[1].text(0.02, -0.04, "blue = split gained in P1", transform=axes[1].transAxes, fontsize=9, color="#1f77b4")
    fig.suptitle("Seven-population ASTRAL trees; SPRE used only for orientation\nRF(P0,P1) = 2", fontsize=13)
    fig.subplots_adjust(top=0.82, bottom=0.16, left=0.03, right=0.98)
    save(fig, "house_mouse_t_complex_astral_7population")


def bin_q(rows: list[dict[str, str]], treatment: str, bin_size: int = 250_000) -> list[dict[str, float]]:
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["treatment"] == treatment:
            start = int(row["start_bp"])
            grouped[5_000_000 + ((start - 5_000_000) // bin_size) * bin_size].append(row)
    return [{"midpoint": start + bin_size / 2, **{key: float(np.mean([float(row[key]) for row in block])) for key in Q_COLORS}} for start, block in sorted(grouped.items())]


def plot_segments(ax, x: np.ndarray, y: np.ndarray, **kwargs) -> None:
    if len(x) == 0:
        return
    label = kwargs.pop("label", None)
    breaks = np.flatnonzero(np.diff(x) > 0.375) + 1
    for index, segment in enumerate(np.split(np.arange(len(x)), breaks)):
        if len(segment) >= 2:
            ax.plot(x[segment], y[segment], label=label if index == 0 else None, **kwargs)


def plot_q_axis(ax, rows: list[dict[str, str]], treatment: str, title: str) -> None:
    raw = [row for row in rows if row["treatment"] == treatment]
    bins = bin_q(rows, treatment)
    xraw = np.array([float(row["midpoint_bp"]) / 1e6 for row in raw])
    xb = np.array([row["midpoint"] / 1e6 for row in bins])
    for key, color in Q_COLORS.items():
        ax.scatter(xraw, [float(row[key]) for row in raw], s=2, alpha=0.08, color=color)
        plot_segments(ax, xb, np.array([row[key] for row in bins]), color=color, lw=1.7, label=key)
    ax.axhline(1 / 3, color="#777777", ls="--", lw=0.7)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Quartet fraction")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, ncol=3)


def load_recomb() -> list[dict[str, str]]:
    rows = read_tsv(RESULTS / "stage2_recombination_state_500kb.tsv")
    for row in rows:
        total = sum(float(row[key]) for key in RECOMB_COLORS)
        assert abs(total - 1.0) < 1e-8, f"recombination fractions do not sum to 1: {row}"
    return rows


def plot_recomb_axes(axes, rows: list[dict[str, str]], title: str | None = None, all_ylabels: bool = True) -> None:
    for ax, species in zip(axes, ("domesticus", "musculus", "castaneus")):
        sub = [row for row in rows if row["subspecies"] == species]
        x = np.array([(int(row["bin_start"]) + 250_000) / 1e6 for row in sub])
        bottom = np.zeros(len(sub))
        for key, color in RECOMB_COLORS.items():
            values = np.array([float(row[key]) for row in sub])
            ax.bar(x, values, width=0.49, bottom=bottom, color=color, label=key.replace("fraction_", "").replace("_", " "))
            bottom += values
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction of 5-kb windows" if all_ylabels else "", fontsize=8)
        ax.set_yticks([0, 0.5, 1.0])
        ax.text(0.01, 0.82, species, transform=ax.transAxes, fontsize=9, weight="bold")
        ax.grid(axis="y", alpha=0.2)
    axes[-1].set_xlabel("chr17 position (Mb)")
    if title:
        axes[0].set_title(title)
    if all_ylabels:
        axes[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")


def make_q_tracks() -> None:
    rows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    plot_q_axis(axes[0], rows, "STANDARD_ONLY", "STANDARD_ONLY; 250-kb mean, raw windows light")
    plot_q_axis(axes[1], rows, "ALL_TIPS", "ALL_TIPS; 250-kb mean, raw windows light")
    axes[1].set_xlabel("chr17 position (Mb)")
    fig.suptitle("Direct local quartet support along chr17", fontsize=14)
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_q_tracks")


def make_recombination_track() -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 6.5), sharex=True)
    plot_recomb_axes(axes, load_recomb(), "Source-defined phylogeny-based recombination-state track; 500-kb bins")
    fig.tight_layout()
    save(fig, "house_mouse_t_complex_recombination_track")


def make_q_recombination() -> None:
    qrows, recomb = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv"), load_recomb()
    fig = plt.figure(figsize=(13, 11))
    grid = fig.add_gridspec(5, 1, height_ratios=[1.2, 0.55, 0.8, 0.8, 0.8], hspace=0.25)
    axq = fig.add_subplot(grid[0])
    plot_q_axis(axq, qrows, "ALL_TIPS", "ALL_TIPS quartet support; 250-kb means")
    axd = fig.add_subplot(grid[1], sharex=axq)
    bins = bin_q(qrows, "ALL_TIPS")
    plot_segments(axd, np.array([r["midpoint"] / 1e6 for r in bins]), np.array([r["q_species"] - r["q_t_alt"] for r in bins]), color="#333333", lw=1.5)
    axd.axhline(0, color="black", lw=0.7)
    axd.set_ylabel("delta")
    axd.set_title("q_species - q_t_alt")
    axes = [fig.add_subplot(grid[i], sharex=axq) for i in (2, 3, 4)]
    plot_recomb_axes(axes, recomb)
    axes[-1].set_xlabel("chr17 position (Mb)")
    fig.suptitle("Quartet support aligned with source-defined phylogeny-based recombination states", fontsize=14)
    save(fig, "house_mouse_t_complex_q_recombination")


def arrangement_panel(ax) -> None:
    rows, patterns = read_tsv(RESULTS / "stage2_arrangement_pattern_summary.tsv"), ["SSS", "one_T", "two_T", "TTT"]
    data, bottom = {row["status_pattern"]: row for row in rows}, np.zeros(4)
    for key, label, color in [("q_species", "Q_SPECIES", Q_COLORS["q_species"]), ("q_t_alt", "Q_T_ALT", Q_COLORS["q_t_alt"]), ("q_other", "Q_OTHER", Q_COLORS["q_other"]), ("q_unresolved", "unresolved", "#bbbbbb")]:
        values = np.array([float(data[p][key]) for p in patterns])
        ax.bar(patterns, values, bottom=bottom, label=label, color=color)
        bottom += values
    ax.set_ylim(0, 1)
    ax.set_ylabel("Quartet fraction")
    ax.set_title("D. Arrangement-state quartet composition")
    ax.legend(frameon=False, fontsize=8)


def make_controls() -> None:
    balanced, linkage = read_tsv(RESULTS / "stage2_balanced_quartet_resampling_summary.tsv"), read_tsv(RESULTS / "stage2_linkage_summary.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bmap, treatments, bottom = {row["treatment"]: row for row in balanced}, ["B0_STANDARD_MATCHED", "B1_T_MATCHED", "B2_MIXED_BALANCED"], np.zeros(3)
    for key, label, color in [("fraction_Q_SPECIES", "Q_SPECIES", Q_COLORS["q_species"]), ("fraction_Q_T_ALT", "Q_T_ALT", Q_COLORS["q_t_alt"]), ("fraction_Q_OTHER", "Q_OTHER", Q_COLORS["q_other"])]:
        values = np.array([float(bmap[t][key]) for t in treatments])
        axes[0].bar(treatments, values, bottom=bottom, label=label, color=color)
        bottom += values
    axes[0].set_ylim(0, 1); axes[0].set_title("Balanced exact resampling (n=1000)"); axes[0].tick_params(axis="x", rotation=25); axes[0].legend(frameon=False)
    x = [int(row["spacing_kb"]) for row in linkage]
    axes[1].plot(x, [float(row["fraction_Q_T_ALT"]) for row in linkage], marker="o", color=Q_COLORS["q_t_alt"], label="Q_T_ALT frequency")
    axes[1].plot(x, [float(row["median_localPP"]) for row in linkage], marker="s", color="#d62728", label="median localPP")
    axes[1].set_xscale("log"); axes[1].set_ylim(0, 1.05); axes[1].set_xlabel("spacing (kb)"); axes[1].set_title("Spatial thinning control"); axes[1].legend(frameon=False)
    fig.tight_layout(); save(fig, "house_mouse_t_complex_controls")


def draw_main(which: str) -> None:
    qrows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    fig = plt.figure(figsize=(15, 12))
    outer = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], height_ratios=[1.15, 1], hspace=0.28, wspace=0.24)
    tree_grid = outer[0, 0].subgridspec(1, 2, wspace=0.15)
    trees = {"STANDARD": load_tree(RESULTS / "stage1_aster" / "T0_STANDARD_subspecies.nwk", "Mus_spretus"), "ALL": load_tree(RESULTS / "stage1_aster" / "T1_ALL_WINDOWS_subspecies.nwk", "Mus_spretus")}
    labels = [PRETTY[x] for x in FOUR_LABELS]
    axes_tree = [fig.add_subplot(tree_grid[i]) for i in range(2)]
    draw_cladogram(axes_tree[0], trees["STANDARD"], labels, "A. STANDARD\nQ_SPECIES", Q_COLORS["q_species"], focal_split=focal_split("Q_SPECIES"))
    draw_cladogram(axes_tree[1], trees["ALL"], labels, "ALL TIPS\nQ_T_ALT", Q_COLORS["q_t_alt"], focal_split=focal_split("Q_T_ALT"))
    axq = fig.add_subplot(outer[0, 1]); plot_q_axis(axq, qrows, "ALL_TIPS", "B. ALL_TIPS direct quartet support")
    axes_rec = [fig.add_subplot(outer[1, 0].subgridspec(3, 1, hspace=0.05)[i]) for i in range(3)]
    plot_recomb_axes(axes_rec, load_recomb(), "C. Source-defined phylogeny-based recombination classification", all_ylabels=False)
    fig.text(0.015, 0.27, "Fraction of 5-kb windows", rotation=90, va="center", fontsize=8)
    axarr = fig.add_subplot(outer[1, 1]); arrangement_panel(axarr)
    title = "House-mouse t-complex: structural-state mixture, quartet distortion, and summary-tree sensitivity" if which == "v3" else "House-mouse t-complex: linked genealogy distortion and structural-state context"
    fig.suptitle(title, fontsize=15)
    if which == "v2":
        # Preserve both the historical and compatibility filenames.
        save(fig, ["house_mouse_t_complex_stage2_main", "house_mouse_t_complex_main_v2"])
    else:
        save(fig, "house_mouse_t_complex_main_v3")


def quartet_by_recombination() -> None:
    qrows = {row["locus_id"]: row for row in read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv") if row["treatment"] == "ALL_TIPS"}
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(RESULTS / "stage2_recombination_state_5kb.tsv"):
        if row["locus_id"] in qrows:
            groups[(row["subspecies"], row["recombination_class"])].append(qrows[row["locus_id"]])
    out = []
    for (species, category), subset in sorted(groups.items()):
        out.append({"subspecies": species, "recombination_class": category, "n_windows": len(subset), "mean_q_species": np.mean([float(r["q_species"]) for r in subset]), "mean_q_t_alt": np.mean([float(r["q_t_alt"]) for r in subset]), "mean_q_other": np.mean([float(r["q_other"]) for r in subset]), "mean_delta_species_alt": np.mean([float(r["delta_species_alt"]) for r in subset]), "median_q_species": np.median([float(r["q_species"]) for r in subset]), "median_q_t_alt": np.median([float(r["q_t_alt"]) for r in subset]), "median_q_other": np.median([float(r["q_other"]) for r in subset]), "median_delta_species_alt": np.median([float(r["delta_species_alt"]) for r in subset])})
    write_tsv(RESULTS / "stage2_quartet_by_recombination_state.tsv", out, list(out[0].keys()))


def validate_tree_semantics() -> None:
    four = [("T0_STANDARD_subspecies.nwk", "Q_SPECIES"), ("T1_ALL_WINDOWS_subspecies.nwk", "Q_T_ALT")]
    expected = set(PRETTY[x] for x in FOUR_LABELS)
    annotations = focal_annotations()
    for filename, inferred in four:
        tree = load_tree(RESULTS / "stage1_aster" / filename, "Mus_spretus")
        assert set(t.name for t in tree.get_terminals()) == expected
        topology = four_taxon_topology_from_newick((RESULTS / "stage1_aster" / filename).read_text().strip())
        assert topology == inferred
        key = "STANDARD" if "T0_" in filename else "ALL"
        assert abs(sum(float(annotations[key][q]) for q in ("q_species", "q_t_alt", "q_other")) - 1) < 1e-6
        if key == "STANDARD":
            assert annotations[key]["q_species"] > annotations[key]["q_t_alt"]
        else:
            assert annotations[key]["q_t_alt"] > annotations[key]["q_species"]
    p0 = (RESULTS / "stage2_population_astral" / "P0_STANDARD.nwk").read_text().strip()
    p1 = (RESULTS / "stage2_population_astral" / "P1_ALL.nwk").read_text().strip()
    assert len(tree_bipartitions(p0) ^ tree_bipartitions(p1)) == 2
    comparison = read_tsv(RESULTS / "stage2_population_split_comparison.tsv")
    raw_changed = {}
    for row in comparison:
        if row["status"] in {"lost", "gained"}:
            left, right = row["split_id"].split(" || ")
            raw_changed[row["status"]] = canonical_split([left.split("|"), right.split("|")])
    assert raw_changed["lost"] in tree_bipartitions(p0) and raw_changed["lost"] not in tree_bipartitions(p1)
    assert raw_changed["gained"] in tree_bipartitions(p1) and raw_changed["gained"] not in tree_bipartitions(p0)


def validate_figures() -> None:
    stems = ["house_mouse_t_complex_astral_4group", "house_mouse_t_complex_astral_7population", "house_mouse_t_complex_q_tracks", "house_mouse_t_complex_recombination_track", "house_mouse_t_complex_q_recombination", "house_mouse_t_complex_controls", "house_mouse_t_complex_stage2_main", "house_mouse_t_complex_main_v2", "house_mouse_t_complex_main_v3"]
    rows = []
    for stem in stems:
        png, pdf = FIGURES / f"{stem}.png", FIGURES / f"{stem}.pdf"
        assert png.stat().st_size > 20_000 and pdf.stat().st_size > 10_000
        with Image.open(png) as image:
            pixels = np.asarray(image.convert("RGB")); nonwhite = int(np.sum(np.any(pixels < 245, axis=2))); distinct = int(np.unique(pixels.reshape(-1, 3), axis=0).shape[0])
            assert pixels.shape[1] >= 800 and pixels.shape[0] >= 500 and nonwhite > 5_000 and distinct > 100
            rows.append({"figure": stem, "png_bytes": png.stat().st_size, "pdf_bytes": pdf.stat().st_size, "width_px": pixels.shape[1], "height_px": pixels.shape[0], "nonwhite_pixels": nonwhite, "distinct_rgb": distinct, "intended_tree_panels": 2 if "astral" in stem else "NA"})
    write_tsv(RESULTS / "stage2_visualization_validation.tsv", rows, list(rows[0].keys()))


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    qrows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    for row in qrows:
        assert abs(sum(float(row[key]) for key in ("q_species", "q_t_alt", "q_other", "q_unresolved")) - 1.0) < 1e-8
    validate_tree_semantics(); make_4group_trees(); make_7population_trees(); make_q_tracks(); make_recombination_track(); make_q_recombination(); make_controls(); quartet_by_recombination(); draw_main("v2"); draw_main("v3"); validate_figures()
    print("Wrote visualization suite")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
