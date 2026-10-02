#!/usr/bin/env python3
"""Stage 5 fire-ant final manuscript freeze.

This stage verifies the frozen Stage-0 through Stage-4B chain and regenerates
manuscript-ready summary tables, figures, text, FINAL_ANALYSIS.md, and the
final reproducibility manifest. It performs no new primary inference, no
ASTRAL/ASTER run, no boundary optimization, and no MSRC model fitting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import unittest
from collections import Counter
from decimal import Decimal, getcontext
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


getcontext().prec = 60

DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMPIRICAL_ROOT / "results"
FIGURES = EMPIRICAL_ROOT / "figures"

STAGE0_FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
STAGE1_MANIFEST = RESULTS / "stage1_manifest.json"
STAGE2_MANIFEST = RESULTS / "stage2_manifest.json"
STAGE2_FOCAL_PARTITION = DATA_ROOT / "processed" / "stage2_focal_partition.tsv"
STAGE3_MANIFEST = RESULTS / "stage3_manifest.json"
STAGE3_BACKGROUND_PARTITION = DATA_ROOT / "processed" / "stage3_background_partition.tsv"
STAGE4A_MANIFEST = RESULTS / "stage4a_manifest.json"
STAGE4A_WINDOW_SUPPORT = DATA_ROOT / "processed" / "stage4a_window_quartet_support.tsv"
STAGE4A_REGION_SUMMARY = RESULTS / "stage4a_region_summary.tsv"
STAGE4A_PRIMARY_CONTRAST = RESULTS / "stage4a_primary_contrast.tsv"
STAGE4B_MANIFEST = RESULTS / "stage4b_manifest.json"
STAGE4B_PRIMARY_TEST = RESULTS / "stage4b_primary_test.tsv"
STAGE4B_COORDINATE_SUMMARY = RESULTS / "stage4b_coordinate_summary.tsv"
STAGE4B_LENGTH_SUMMARY = RESULTS / "stage4b_coordinate_length_weighted_summary.tsv"
STAGE4B_CIRCULAR_SHIFTS = RESULTS / "stage4b_circular_shifts.tsv"
STAGE4B_COORDINATE_PLACEMENTS = RESULTS / "stage4b_coordinate_placements.tsv"
STAGE4B_LENGTH_WEIGHTED = RESULTS / "stage4b_coordinate_length_weighted.tsv"
STAGE4B_PHYSICAL_ORDER = RESULTS / "stage4b_chr16_physical_order.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"
README = EMPIRICAL_ROOT / "README.md"

ROBUSTNESS_SUMMARY = RESULTS / "stage5_robustness_summary.tsv"
MAIN_TABLE_TSV = RESULTS / "fire_ants_main_table.tsv"
MAIN_TABLE_MD = RESULTS / "fire_ants_main_table.md"
MAIN_CAPTION = RESULTS / "fire_ants_main_figure_caption.md"
RESULTS_TEXT = RESULTS / "fire_ants_results_text.md"
METHODS_TEXT = RESULTS / "fire_ants_methods_text.md"
FINAL_ANALYSIS = RESULTS / "FINAL_ANALYSIS.md"
MANIFEST = RESULTS / "stage5_final_manifest.json"
MAIN_FIG_PDF = FIGURES / "fire_ants_main.pdf"
MAIN_FIG_PNG = FIGURES / "fire_ants_main.png"
SUPP_FIG_PDF = FIGURES / "fire_ants_supplement.pdf"
SUPP_FIG_PNG = FIGURES / "fire_ants_supplement.png"

EXPECTED = {
    "stage0_focal": "5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6",
    "stage1_manifest": "fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b",
    "stage2_manifest": "2b56ac7ea4a4d0a31588b803b30896522af4e90bc1e91bf23db6d0b9db80e427",
    "stage2_focal": "52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002",
    "stage3_manifest": "4ed41d1558dcc055efbb6522efe11b273062a187a7335a2975597fdbe0216fd7",
    "stage3_background": "fa4cc04afcc40deb5e6b3d5ca76dc55b3720dfe070f7d90c3939e253b8006286",
    "stage4a_manifest": "d1fa79815bbe4d1f45a2af7f5a1d5487d39fe3041ccac4f79f961a998ab99bf1",
    "stage4a_window_support": "b8ff3ea95ba31ae5553f56f0aeb45a393bdb0ba2b19081fc1aca822dfdd70758",
    "stage4b_manifest": "3aa375b5b1f52964e338a93aa7cb9330e47fe88986fa3205bfcfd26711df7a2b",
}

DECIMAL_PLACES = Decimal("0.000000000000000")
TSV_LINETERMINATOR = "\n"


class Stage5Error(RuntimeError):
    """Raised when the final fire-ant freeze cannot be verified."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else str(row.get(field, "")) for field in fields})


def dec(value: str | Decimal) -> Decimal:
    return value if isinstance(value, Decimal) else Decimal(value)


def fmt15(value: Decimal) -> str:
    return str(value.quantize(DECIMAL_PLACES))


def fmt3(value: str | Decimal) -> str:
    return f"{float(dec(value)):.3f}"


def check_sha(path: Path, expected: str, label: str) -> None:
    observed = sha256(path)
    if observed != expected:
        raise Stage5Error(f"{label} checksum mismatch: expected {expected}, observed {observed}")


def verify_recorded_checksums(entries: dict[str, str], label: str) -> None:
    for path_text, expected in sorted(entries.items()):
        path = REPO_ROOT / path_text
        if not path.exists():
            raise Stage5Error(f"{label} recorded path is missing: {path_text}")
        check_sha(path, expected, f"{label}: {path_text}")


def verify_chain() -> dict[str, object]:
    check_sha(STAGE0_FOCAL_QUARTET, EXPECTED["stage0_focal"], "Stage-0 focal quartet")
    check_sha(STAGE1_MANIFEST, EXPECTED["stage1_manifest"], "Stage-1 manifest")
    check_sha(STAGE2_MANIFEST, EXPECTED["stage2_manifest"], "Stage-2 manifest")
    check_sha(STAGE2_FOCAL_PARTITION, EXPECTED["stage2_focal"], "Stage-2 focal partition")
    check_sha(STAGE3_MANIFEST, EXPECTED["stage3_manifest"], "Stage-3 manifest")
    check_sha(STAGE3_BACKGROUND_PARTITION, EXPECTED["stage3_background"], "Stage-3 background partition")
    check_sha(STAGE4A_MANIFEST, EXPECTED["stage4a_manifest"], "Stage-4A manifest")
    check_sha(STAGE4A_WINDOW_SUPPORT, EXPECTED["stage4a_window_support"], "Stage-4A window support")
    check_sha(STAGE4B_MANIFEST, EXPECTED["stage4b_manifest"], "Stage-4B manifest")

    stage1 = json.loads(STAGE1_MANIFEST.read_text())
    for key, expected in {
        "coordinate_rows": "213",
        "local_tree_rows": "213",
        "min_tree_tips": "267",
        "max_tree_tips": "267",
        "twisst_topologies": "945",
        "weight_rows": "213",
    }.items():
        if str(stage1["observed_metrics"].get(key)) != expected:
            raise Stage5Error(f"Stage-1 metric {key} changed")
    verify_recorded_checksums(stage1.get("raw_topology_file_checksums", {}), "Stage-1 raw")
    verify_recorded_checksums(stage1.get("processed_output_checksums", {}), "Stage-1 processed")

    stage2 = json.loads(STAGE2_MANIFEST.read_text())
    verify_recorded_checksums(stage2.get("stage2_output_checksums", {}), "Stage-2 output")

    stage3 = json.loads(STAGE3_MANIFEST.read_text())
    if stage3.get("stage3_background_partition_matches_stage0_species_split") is not True:
        raise Stage5Error("Stage-3 background split is not AB|CD")
    if stage3.get("stage3_background_partition_checksum") != EXPECTED["stage3_background"]:
        raise Stage5Error("Stage-3 background partition manifest checksum changed")

    stage4a = json.loads(STAGE4A_MANIFEST.read_text())
    verify_recorded_checksums(stage4a.get("stage4a_output_checksums", {}), "Stage-4A output")
    verify_recorded_checksums(stage4a.get("figure_checksums", {}), "Stage-4A figure")

    stage4b = json.loads(STAGE4B_MANIFEST.read_text())
    stage4b_outputs = {
        path: checksum
        for path, checksum in stage4b.get("output_checksums", {}).items()
        if path != "empirical/fire_ants_chr16/README.md"
    }
    verify_recorded_checksums(stage4b_outputs, "Stage-4B output")
    verify_recorded_checksums(stage4b.get("figure_checksums", {}), "Stage-4B figure")
    if stage4b["physical_chr16_region_sequence"] != ["chr16A", "chr16_supergene", "chr16B"]:
        raise Stage5Error("Stage-4B physical chr16 ordering changed")

    return {"stage1": stage1, "stage2": stage2, "stage3": stage3, "stage4a": stage4a, "stage4b": stage4b}


def load_inputs() -> dict[str, object]:
    _, region_rows = read_tsv(STAGE4A_REGION_SUMMARY)
    regions = {row["region"]: row for row in region_rows}
    _, contrast_rows = read_tsv(STAGE4A_PRIMARY_CONTRAST)
    primary_contrast = next(row for row in contrast_rows if row["comparison"] == "chr16_supergene_vs_chr16_outside")
    _, window_rows = read_tsv(STAGE4A_WINDOW_SUPPORT)
    _, primary_rows = read_tsv(STAGE4B_PRIMARY_TEST)
    _, coordinate_rows = read_tsv(STAGE4B_COORDINATE_SUMMARY)
    _, length_rows = read_tsv(STAGE4B_LENGTH_SUMMARY)
    _, circular_rows = read_tsv(STAGE4B_CIRCULAR_SHIFTS)
    _, coordinate_placements = read_tsv(STAGE4B_COORDINATE_PLACEMENTS)
    _, length_weighted = read_tsv(STAGE4B_LENGTH_WEIGHTED)
    _, region_manifest = read_tsv(REGION_MANIFEST)
    return {
        "regions": regions,
        "primary_contrast": primary_contrast,
        "windows": window_rows,
        "primary_test": primary_rows[0],
        "coordinate_summary": coordinate_rows[0],
        "length_summary": length_rows[0],
        "circular_rows": circular_rows,
        "coordinate_placements": coordinate_placements,
        "length_weighted": length_weighted,
        "region_manifest": {row["region_id"]: row for row in region_manifest},
    }


def verify_final_values(inputs: dict[str, object]) -> None:
    regions = inputs["regions"]
    expected_regions = {
        "chr16_outside": ("44", "0.977638072355247", "0.011196452159259", "0.011165475485494", "-0.966441620195988", "44", "0", "0"),
        "chr16_supergene": ("52", "0.167775556643584", "0.811443799942942", "0.020780643413474", "0.643668243299358", "8", "43", "1"),
        "chr1": ("117", "0.983232450588138", "0.008234827754813", "0.008532721657050", "-0.974997622833325", "117", "0", "0"),
    }
    for region, expected in expected_regions.items():
        row = regions[region]
        observed = (
            row["n_windows"],
            row["mean_q_S"],
            row["mean_q_H"],
            row["mean_q_3"],
            row["mean_D"],
            row["n_species_dominant"],
            row["n_haplotype_dominant"],
            row["n_third_dominant"],
        )
        if observed != expected:
            raise Stage5Error(f"Stage-4A summary changed for {region}: {observed}")
    if inputs["primary_contrast"]["delta_D"] != "1.610109863495346":
        raise Stage5Error("Stage-4A Delta_D changed")
    primary = inputs["primary_test"]
    if (primary["observed_rank"], primary["n_exact_alignments"], primary["p_one_sided"]) != ("2", "96", "0.0208333333"):
        raise Stage5Error("Stage-4B circular result changed")
    unique = inputs["coordinate_summary"]
    if (unique["rank_unique_membership"], unique["n_unique_membership_states"], unique["p_unique_membership"]) != ("3", "97", "0.0309278351"):
        raise Stage5Error("Stage-4B unique-membership result changed")
    length = inputs["length_summary"]
    if (length["total_start_domain_length"], length["extreme_start_domain_length"], length["p_length_weighted"]) != ("12682605", "993416", "0.0783290184"):
        raise Stage5Error("Stage-4B length-weighted result changed")


def manuscript_rows(inputs: dict[str, object]) -> list[dict[str, object]]:
    regions = inputs["regions"]
    return [
        {"region": "chr1", **regions["chr1"]},
        {"region": "chr16_outside", **regions["chr16_outside"]},
        {"region": "chr16_supergene", **regions["chr16_supergene"]},
    ]


def write_summary_tables(inputs: dict[str, object]) -> None:
    robust_rows = [
        {
            "analysis": "exact_circular_shift",
            "observed_delta_D": inputs["primary_test"]["observed_delta_D"],
            "n_or_domain": "96 exact circular alignments",
            "rank": inputs["primary_test"]["observed_rank"],
            "p_value": inputs["primary_test"]["p_one_sided"],
            "role": "primary",
            "preserves": "ordered chr16 D track; 52 inside and 44 outside mask positions",
            "interpretation": "positive_spatial_alignment",
        },
        {
            "analysis": "unique_membership_coordinate",
            "observed_delta_D": inputs["coordinate_summary"]["observed_delta_D"],
            "n_or_domain": "97 unique sampled-window membership states",
            "rank": inputs["coordinate_summary"]["rank_unique_membership"],
            "p_value": inputs["coordinate_summary"]["p_unique_membership"],
            "role": "secondary_sensitivity",
            "preserves": "physical interval width; irregular sampled-window coordinates; equal state weights",
            "interpretation": "positive_spatial_alignment",
        },
        {
            "analysis": "continuous_uniform_physical_start",
            "observed_delta_D": inputs["length_summary"]["observed_delta_D"],
            "n_or_domain": "12682605 bp valid start-coordinate domain",
            "rank": "",
            "p_value": inputs["length_summary"]["p_length_weighted"],
            "role": "secondary_sensitivity",
            "preserves": "physical interval width; uniform physical start coordinate; length-weighted states",
            "interpretation": "positive_spatial_alignment",
        },
    ]
    write_tsv(ROBUSTNESS_SUMMARY, robust_rows, ["analysis", "observed_delta_D", "n_or_domain", "rank", "p_value", "role", "preserves", "interpretation"])

    main_rows = []
    for row in manuscript_rows(inputs):
        main_rows.append(
            {
                "region": row["region"],
                "n_windows": row["n_windows"],
                "mean_q_species": row["mean_q_S"],
                "mean_q_haplotype": row["mean_q_H"],
                "mean_q_third": row["mean_q_3"],
                "mean_D": row["mean_D"],
                "species_dominant": row["n_species_dominant"],
                "haplotype_dominant": row["n_haplotype_dominant"],
                "third_dominant": row["n_third_dominant"],
            }
        )
    fields = ["region", "n_windows", "mean_q_species", "mean_q_haplotype", "mean_q_third", "mean_D", "species_dominant", "haplotype_dominant", "third_dominant"]
    write_tsv(MAIN_TABLE_TSV, main_rows, fields)

    md_lines = [
        "| region | n windows | q_S | q_H | q_3 | D | species dom. | haplotype dom. | third dom. |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in main_rows:
        md_lines.append(
            f"| {row['region']} | {row['n_windows']} | {fmt3(row['mean_q_species'])} | {fmt3(row['mean_q_haplotype'])} | {fmt3(row['mean_q_third'])} | {fmt3(row['mean_D'])} | {row['species_dominant']} | {row['haplotype_dominant']} | {row['third_dominant']} |"
        )
    md_lines.extend(
        [
            "",
            f"`Delta_D = {inputs['primary_contrast']['delta_D']}`.",
            f"Primary exact circular-shift `p = {inputs['primary_test']['p_one_sided']}`.",
            f"Coordinate unique-membership sensitivity `p = {inputs['coordinate_summary']['p_unique_membership']}`.",
            f"Coordinate length-weighted physical-start sensitivity `p = {inputs['length_summary']['p_length_weighted']}`.",
        ]
    )
    MAIN_TABLE_MD.write_text("\n".join(md_lines) + "\n")


def sorted_windows(windows: list[dict[str, str]], chrom: str) -> list[dict[str, str]]:
    return sorted([row for row in windows if row["chrom"] == chrom], key=lambda row: Decimal(row["mid"]))


def region_span(region_manifest: dict[str, dict[str, str]], region: str) -> tuple[float, float]:
    row = region_manifest[region]
    return float(Decimal(row["coordinate_start"])) / 1_000_000, float(Decimal(row["coordinate_end"])) / 1_000_000


def plot_support_by_region(ax, windows: list[dict[str, str]], chrom: str, region_order: list[str], region_manifest: dict[str, dict[str, str]], shade: bool) -> None:
    colors = {"q_S": "#345995", "q_H": "#b23a48", "q_3": "#2f7f4f"}
    labels = {"q_S": "q_S species-history", "q_H": "q_H SB/Sb haplotype", "q_3": "q_3 third resolution"}
    for field in ["q_S", "q_H", "q_3"]:
        first = True
        for region in region_order:
            subset = sorted([row for row in windows if row["chrom"] == chrom and row["region"] == region], key=lambda row: Decimal(row["mid"]))
            if not subset:
                continue
            ax.plot(
                [float(Decimal(row["mid"])) / 1_000_000 for row in subset],
                [float(Decimal(row[field])) for row in subset],
                marker="o",
                markersize=3.0,
                linewidth=0.85,
                color=colors[field],
                label=labels[field] if first else None,
                alpha=0.9,
            )
            first = False
    if shade:
        start, end = region_span(region_manifest, "chr16_supergene")
        ax.axvspan(start, end, color="#d9d9d9", alpha=0.38, label="author-designated supergene analysis region")
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.set_ylabel("TWISST-derived focal quartet support")


def make_main_figure(inputs: dict[str, object]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 7.2), height_ratios=[1.25, 1.0])
    plot_support_by_region(axes[0], inputs["windows"], "chr16", ["chr16A", "chr16_supergene", "chr16B"], inputs["region_manifest"], True)
    axes[0].set_xlabel("chr16 physical midpoint (Mb)")
    axes[0].set_title("A. Chromosome-16 focal quartet support")
    axes[0].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    axes[0].text(0.02, 0.06, "outside: 44/44 species-dominant\nsupergene: 43/52 haplotype-dominant", transform=axes[0].transAxes, fontsize=8, bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.85})

    circular = sorted(inputs["circular_rows"], key=lambda row: Decimal(row["delta_D"]))
    axes[1].scatter(range(1, len(circular) + 1), [float(Decimal(row["delta_D"])) for row in circular], s=24, color="#345995", linewidth=0)
    axes[1].axhline(float(Decimal(inputs["primary_test"]["observed_delta_D"])), color="#b23a48", linewidth=1.2, label="observed Delta_D")
    axes[1].set_title("B. Exact circular-shift null")
    axes[1].set_xlabel("ranked circular alignment")
    axes[1].set_ylabel("Delta_D")
    axes[1].grid(axis="y", color="#eeeeee", linewidth=0.6)
    axes[1].text(0.03, 0.88, "observed rank = 2/96\nexact one-sided p = 0.0208", transform=axes[1].transAxes, fontsize=9, bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.85})
    axes[1].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    with PdfPages(MAIN_FIG_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(MAIN_FIG_PNG, dpi=220)
    plt.close(fig)


def make_supplement_figure(inputs: dict[str, object]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))
    ax = axes[0, 0]
    plot_support_by_region(ax, inputs["windows"], "chr1", ["chr1"], inputs["region_manifest"], False)
    ax.set_xlabel("chr1 physical midpoint (Mb)")
    ax.set_title("A. Chromosome-1 control")
    inset = inset_axes(ax, width="38%", height="42%", loc="lower right", borderpad=1.2)
    chr1 = sorted_windows(inputs["windows"], "chr1")
    for field, color in [("q_H", "#b23a48"), ("q_3", "#2f7f4f")]:
        inset.plot([float(Decimal(row["mid"])) / 1_000_000 for row in chr1], [float(Decimal(row[field])) for row in chr1], marker="o", markersize=1.8, linewidth=0.55, color=color, alpha=0.85)
    inset.set_ylim(0, 0.15)
    inset.set_title("q_H/q_3 zoom", fontsize=8)
    inset.tick_params(labelsize=7)

    ax = axes[0, 1]
    for region in ["chr16A", "chr16_supergene", "chr16B"]:
        subset = sorted([row for row in inputs["windows"] if row["chrom"] == "chr16" and row["region"] == region], key=lambda row: Decimal(row["mid"]))
        ax.plot([float(Decimal(row["mid"])) / 1_000_000 for row in subset], [float(Decimal(row["D"])) for row in subset], marker="o", markersize=3.0, linewidth=0.85, color="#345995")
    start, end = region_span(inputs["region_manifest"], "chr16_supergene")
    ax.axvspan(start, end, color="#d9d9d9", alpha=0.38)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("chr16 physical midpoint (Mb)")
    ax.set_ylabel("D = q_H - q_S")
    ax.set_title("B. Chromosome-16 D track")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)

    ax = axes[1, 0]
    placements = sorted(inputs["coordinate_placements"], key=lambda row: Decimal(row["interval_start"]))
    ax.scatter([(float(Decimal(row["interval_start"])) + float(Decimal(row["interval_end"]))) / 2_000_000 for row in placements], [float(Decimal(row["delta_D"])) for row in placements], s=20, color="#4c78a8", linewidth=0)
    ax.axhline(float(Decimal(inputs["primary_test"]["observed_delta_D"])), color="#b23a48", linewidth=1.1)
    ax.set_xlabel("candidate interval midpoint (Mb)")
    ax.set_ylabel("Delta_D")
    ax.set_title("C. Equal weight per distinct sampled-window membership state")
    ax.text(0.03, 0.88, "p = 0.0309", transform=ax.transAxes, fontsize=9, bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.85})
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)

    ax = axes[1, 1]
    for row in inputs["length_weighted"]:
        x0 = float(Decimal(row["start_min"])) / 1_000_000
        x1 = float(Decimal(row["start_max"])) / 1_000_000
        y = float(Decimal(row["delta_D"]))
        color = "#b23a48" if row["exceeds_or_equals_observed"] == "true" else "#4c78a8"
        ax.hlines(y, x0, x1, color=color, linewidth=1.8, alpha=0.88)
    ax.axhline(float(Decimal(inputs["primary_test"]["observed_delta_D"])), color="#b23a48", linewidth=1.1)
    ax.axvline(float(Decimal(inputs["region_manifest"]["chr16_supergene"]["coordinate_start"])) / 1_000_000, color="#222222", linestyle=":", linewidth=1.0)
    ax.set_xlabel("interval start coordinate (Mb)")
    ax.set_ylabel("Delta_D")
    ax.set_title("D. Uniform physical interval-start sensitivity")
    ax.text(0.03, 0.88, "p = 0.0783", transform=ax.transAxes, fontsize=9, bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.85})
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    fig.tight_layout()
    with PdfPages(SUPP_FIG_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(SUPP_FIG_PNG, dpi=220)
    plt.close(fig)


def write_text_outputs(inputs: dict[str, object]) -> None:
    MAIN_CAPTION.write_text(
        "# Fire-ant main figure caption\n\n"
        "TWISST-derived focal quartet support across the fire-ant chromosome-16 social-supergene analysis. The focal quartet was defined as A = `invicta/macdonaghi_SB`, B = `invicta/macdonaghi_Sb`, C = `richteri_SB`, and D = `richteri_Sb`. `q_S` is support for the independently frozen species-history split AB|CD, `q_H` is support for the SB/Sb haplotype-partition split AC|BD, and `q_3` is support for the third unrooted quartet resolution AD|BC. Support values were obtained by aggregating the authors' published TWISST topology weights; structural/group definitions and the chromosome 1-15 background species split were frozen independently before formal local-support analysis. Each point is one published four-BUSCO-gene window. The shaded region is the author-designated supergene analysis interval / observed BUSCO-window span, not exact inversion breakpoints. Panel B shows the exact 96-alignment circular-shift null for the primary spatial statistic. These results do not rule out introgression.\n"
    )

    RESULTS_TEXT.write_text(
        "# Fire-ant Results text\n\n"
        "Focal species and SB/Sb haplotype partitions were defined before inspecting local quartet support. The focal groups were `invicta/macdonaghi_SB`, `invicta/macdonaghi_Sb`, `richteri_SB`, and `richteri_Sb`, and the chromosome 1-15 ASTRAL species history independently defined the AB|CD species-history split. We then unblinded the published local RAxML/TWISST windows and aggregated the authors' topology weights into focal quartet support.\n\n"
        "Collinear chromosome-16 windows outside the supergene were dominated by the species-history topology: mean `q_S = 0.978`, mean `q_H = 0.011`, and 44/44 windows were species-dominant. In the author-designated supergene region, support shifted toward the SB/Sb haplotype partition: mean `q_S = 0.168`, mean `q_H = 0.811`, and 43/52 windows were haplotype-dominant. The resulting `Delta_D` was 1.610. The chromosome-1 control remained species-dominant in 117/117 windows. Thus, the supergene does not simply increase arbitrary discordance; support shifts specifically toward the independently defined cross-species SB/Sb haplotype partition while the third topology remains low.\n\n"
        "Under the primary exact circular-shift null, the observed supergene alignment ranked second among 96 possible circular alignments of the physically ordered chromosome-16 `D` track, giving a one-sided exact `p = 0.0208`. Coordinate-aware sensitivities gave `p = 0.0309` when equal weight was assigned to each unique sampled-window membership state and `p = 0.0783` when interval starts were distributed uniformly over physical coordinate space. The latter is more conservative because membership states are weighted by the amount of physical start-coordinate space over which they persist.\n\n"
        "The fire-ant supergene therefore defines a coherent alternative genealogy regime across species. The original source study infers recurrent adaptive introgression of `Sb`, and that interpretation remains a central caveat. These results establish association between the recombination-suppressed structural haplotype and local genealogy, but they do not identify whether shared structural ancestry, introgression, or their combination generated the historical pattern. This empirical case is consistent with the manuscript's broader point that structural history and gene flow can be difficult to distinguish from local quartet patterns alone.\n"
    )

    METHODS_TEXT.write_text(
        "# Fire-ant Methods text\n\n"
        "We analyzed the fire-ant chromosome-16 social-supergene dataset from Stolle et al. (2022) using the public repository `wurmlab/2021-fire-ant-social-supergene-introgression`. We used the published local RAxML/TWISST analysis rather than reanalyzing raw sequence data. The normalized dataset contained 213 published four-BUSCO-gene windows: 117 on chr1, 42 in chr16A, 2 in chr16B, and 52 in the author-designated chr16 supergene interval. The TWISST analysis used seven grouped taxa.\n\n"
        "The focal quartet was frozen as A = `invicta/macdonaghi_SB`, B = `invicta/macdonaghi_Sb`, C = `richteri_SB`, and D = `richteri_Sb`. The species-history split was AB|CD, the haplotype split was AC|BD, and the third split was AD|BC. Each grouped TWISST topology was classified by its induced unrooted focal quartet, and topology weights were aggregated per window into `q_S`, `q_H`, and `q_3`. We defined `D = q_H - q_S` and the primary contrast as `Delta_D = mean(D_supergene) - mean(D_chr16 outside)`.\n\n"
        "For spatial inference, chromosome-16 windows were reordered by increasing physical midpoint because the upstream TWISST concatenation order was not physical chromosome order. The primary null was the predeclared exact circular-shift null over all 96 physically ordered chromosome-16 windows, keeping the frozen 52-window supergene mask fixed. The one-sided direction, higher `D` inside the supergene, was predeclared. We also report two coordinate-aware sensitivities using the same frozen interval width: an equal-weight unique-membership-state sensitivity and a continuous length-weighted uniform physical-start sensitivity. No arbitrary coordinate grid was used.\n\n"
        "No raw-sequence reanalysis was performed, no new local trees were inferred, no ASTRAL/ASTER analysis was run in Stage 5, no MSRC parameters were fit, and no supergene boundary was optimized.\n"
    )


def write_final_analysis(inputs: dict[str, object], chain: dict[str, object]) -> None:
    rows = manuscript_rows(inputs)
    table = "\n".join(f"| {row['region']} | {row['n_windows']} | {fmt3(row['mean_q_S'])} | {fmt3(row['mean_q_H'])} | {fmt3(row['mean_q_3'])} | {fmt3(row['mean_D'])} |" for row in rows)
    FINAL_ANALYSIS.write_text(
        "# Fire-ant final analysis freeze\n\n"
        "Fire-ant empirical analysis is frozen for manuscript use.\n\n"
        "## Frozen input checksums\n\n"
        f"- Stage-0 focal quartet: `{sha256(STAGE0_FOCAL_QUARTET)}`\n"
        f"- Stage-1 manifest: `{sha256(STAGE1_MANIFEST)}`\n"
        f"- Stage-2 focal partition: `{sha256(STAGE2_FOCAL_PARTITION)}`\n"
        f"- Stage-2 manifest: `{sha256(STAGE2_MANIFEST)}`\n"
        f"- Stage-3 background partition: `{sha256(STAGE3_BACKGROUND_PARTITION)}`\n"
        f"- Stage-3 manifest: `{sha256(STAGE3_MANIFEST)}`\n"
        f"- Stage-4A manifest: `{sha256(STAGE4A_MANIFEST)}`\n"
        f"- Stage-4A window support: `{sha256(STAGE4A_WINDOW_SUPPORT)}`\n"
        f"- Stage-4B manifest: `{sha256(STAGE4B_MANIFEST)}`\n\n"
        "## Final descriptive result\n\n"
        "| region | n | q_S | q_H | q_3 | D |\n"
        "|---|---:|---:|---:|---:|---:|\n"
        f"{table}\n\n"
        "## Primary spatial result\n\n"
        f"`Delta_D = 1.6101`; circular rank = 2/96; `p = {inputs['primary_test']['p_one_sided']}`.\n\n"
        "## Coordinate sensitivities\n\n"
        f"Unique membership: 3/97, `p = {inputs['coordinate_summary']['p_unique_membership']}`.\n\n"
        f"Uniform physical-start: `p = {inputs['length_summary']['p_length_weighted']}`.\n\n"
        "## Dominant topology counts\n\n"
        "- chr1: 117/117 species\n"
        "- chr16 outside: 44/44 species\n"
        "- supergene: 43/52 haplotype\n\n"
        "## Interpretation\n\n"
        "The chromosome-16 supergene is associated with a pronounced, spatially localized shift from the background species-history quartet toward the cross-species SB/Sb haplotype quartet.\n\n"
        "The source study's recurrent adaptive introgression interpretation remains an explicit caveat. These results do not identify whether shared structural ancestry, introgression, or their combination generated the historical pattern.\n\n"
        "## Manuscript-ready files\n\n"
        "- `empirical/fire_ants_chr16/figures/fire_ants_main.pdf`\n"
        "- `empirical/fire_ants_chr16/figures/fire_ants_main.png`\n"
        "- `empirical/fire_ants_chr16/figures/fire_ants_supplement.pdf`\n"
        "- `empirical/fire_ants_chr16/figures/fire_ants_supplement.png`\n"
        "- `empirical/fire_ants_chr16/results/fire_ants_main_table.tsv`\n"
        "- `empirical/fire_ants_chr16/results/fire_ants_main_table.md`\n"
        "- `empirical/fire_ants_chr16/results/fire_ants_main_figure_caption.md`\n"
        "- `empirical/fire_ants_chr16/results/fire_ants_methods_text.md`\n"
        "- `empirical/fire_ants_chr16/results/fire_ants_results_text.md`\n"
        "- `empirical/fire_ants_chr16/results/stage5_robustness_summary.tsv`\n"
    )


def update_readme() -> None:
    text = README.read_text()
    status_old = "Stage 4B complete — spatial-null inference. It implements the frozen exact circular-shift null and the frozen physical-coordinate-aware same-width interval sensitivity. It does not optimize supergene boundaries, run ASTRAL/ASTER, fit MSRC parameters, or identify a historical mechanism.\n"
    status_new = status_old + "\nStage 5 complete — frozen for manuscript use.\n"
    if status_new not in text:
        if status_old not in text:
            raise Stage5Error("README Stage-4B status block not found")
        text = text.replace(status_old, status_new)
    links = (
        "\n## Stage 5 manuscript freeze\n\n"
        "- `results/FINAL_ANALYSIS.md`\n"
        "- `results/fire_ants_main_table.tsv`\n"
        "- `results/fire_ants_methods_text.md`\n"
        "- `results/fire_ants_results_text.md`\n"
        "- `figures/fire_ants_main.pdf`\n"
        "- `figures/fire_ants_supplement.pdf`\n\n"
        "Stage 6, if performed, is a post-freeze ASTER/ASTRAL4 summary-tree sensitivity analysis and does not alter the Stage-5 primary empirical result.\n"
    )
    if "## Stage 5 manuscript freeze" not in text:
        text += links
    README.write_text(text)


def write_manifest(chain: dict[str, object], inputs: dict[str, object]) -> None:
    final_outputs = [
        ROBUSTNESS_SUMMARY,
        MAIN_TABLE_TSV,
        MAIN_TABLE_MD,
        MAIN_CAPTION,
        RESULTS_TEXT,
        METHODS_TEXT,
        FINAL_ANALYSIS,
        MAIN_FIG_PDF,
        MAIN_FIG_PNG,
        SUPP_FIG_PDF,
        SUPP_FIG_PNG,
        README,
    ]
    data = {
        "stage": 5,
        "analysis_status": "frozen_for_manuscript_use",
        "analysis_frozen": True,
        "frozen_input_checksums": {
            "stage0_focal_quartet": sha256(STAGE0_FOCAL_QUARTET),
            "stage1_manifest": sha256(STAGE1_MANIFEST),
            "stage2_focal_partition": sha256(STAGE2_FOCAL_PARTITION),
            "stage2_manifest": sha256(STAGE2_MANIFEST),
            "stage3_background_partition": sha256(STAGE3_BACKGROUND_PARTITION),
            "stage3_manifest": sha256(STAGE3_MANIFEST),
            "stage4a_manifest": sha256(STAGE4A_MANIFEST),
            "stage4a_window_support": sha256(STAGE4A_WINDOW_SUPPORT),
            "stage4b_manifest": sha256(STAGE4B_MANIFEST),
        },
        "stage4b_central_output_checksums": {
            rel(path): sha256(path)
            for path in [STAGE4B_PRIMARY_TEST, STAGE4B_COORDINATE_SUMMARY, STAGE4B_LENGTH_SUMMARY, STAGE4B_CIRCULAR_SHIFTS, STAGE4B_PHYSICAL_ORDER]
        },
        "final_results": {
            "regions": {
                row["region"]: {
                    "n_windows": row["n_windows"],
                    "mean_q_S": row["mean_q_S"],
                    "mean_q_H": row["mean_q_H"],
                    "mean_q_3": row["mean_q_3"],
                    "mean_D": row["mean_D"],
                    "species_dominant": row["n_species_dominant"],
                    "haplotype_dominant": row["n_haplotype_dominant"],
                    "third_dominant": row["n_third_dominant"],
                }
                for row in manuscript_rows(inputs)
            },
            "Delta_D": inputs["primary_contrast"]["delta_D"],
            "circular_rank": inputs["primary_test"]["observed_rank"],
            "circular_n": inputs["primary_test"]["n_exact_alignments"],
            "circular_p": inputs["primary_test"]["p_one_sided"],
            "unique_coordinate_rank": inputs["coordinate_summary"]["rank_unique_membership"],
            "unique_coordinate_n": inputs["coordinate_summary"]["n_unique_membership_states"],
            "unique_coordinate_p": inputs["coordinate_summary"]["p_unique_membership"],
            "length_weighted_coordinate_p": inputs["length_summary"]["p_length_weighted"],
        },
        "final_output_checksums": {rel(path): sha256(path) for path in final_outputs},
        "reproducibility_command": "python3 empirical/fire_ants_chr16/scripts/05_finalize_fire_ant_analysis.py --run-tests",
        "astral_or_aster_run": False,
        "msrc_model_fit": False,
        "boundaries_optimized": False,
        "new_primary_test_added": False,
        "stage5_script": rel(Path(__file__)),
        "stage5_script_sha256": sha256(Path(__file__)),
    }
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def run_stage5() -> dict[str, object]:
    chain = verify_chain()
    inputs = load_inputs()
    verify_final_values(inputs)
    write_summary_tables(inputs)
    make_main_figure(inputs)
    make_supplement_figure(inputs)
    write_text_outputs(inputs)
    write_final_analysis(inputs, chain)
    update_readme()
    write_manifest(chain, inputs)
    return {
        "analysis_status": "frozen_for_manuscript_use",
        "Delta_D": inputs["primary_contrast"]["delta_D"],
        "circular_p": inputs["primary_test"]["p_one_sided"],
        "unique_coordinate_p": inputs["coordinate_summary"]["p_unique_membership"],
        "length_weighted_coordinate_p": inputs["length_summary"]["p_length_weighted"],
        "manifest_sha": sha256(MANIFEST),
    }


class Stage5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.chain = verify_chain()
        cls.inputs = load_inputs()

    def test_chain_and_values(self) -> None:
        verify_final_values(self.inputs)
        self.assertEqual(sha256(STAGE4B_MANIFEST), EXPECTED["stage4b_manifest"])

    def test_methods_text_contains_frozen_resolutions(self) -> None:
        text = METHODS_TEXT.read_text()
        for phrase in ["AB|CD", "AC|BD", "AD|BC"]:
            self.assertIn(phrase, text)

    def test_results_text_contains_sensitivity_and_introgression(self) -> None:
        text = RESULTS_TEXT.read_text().lower()
        self.assertIn("0.0783", text)
        self.assertIn("introgression", text)
        self.assertIn("do not identify", text)

    def test_figures_and_no_smoothing_claims(self) -> None:
        self.assertTrue(MAIN_FIG_PDF.exists())
        self.assertTrue(SUPP_FIG_PDF.exists())
        script = Path(__file__).read_text().lower()
        self.assertNotIn("." + "rolling", script)
        self.assertNotIn("sav" + "gol", script)
        self.assertNotIn("low" + "ess", script)
        self.assertIn("q_h/q_3 zoom", script)

    def test_no_forbidden_actions(self) -> None:
        script = Path(__file__).read_text().lower()
        self.assertNotIn("import " + "subprocess", script)
        self.assertNotIn("sub" + "process" + ".", script)
        self.assertNotIn("msrc_model_fit" + " = " + "true", script)
        self.assertNotIn("boundaries_optimized" + " = " + "true", script)
        if MANIFEST.exists():
            manifest = json.loads(MANIFEST.read_text())
            self.assertFalse(manifest["astral_or_aster_run"])
            self.assertFalse(manifest["msrc_model_fit"])
            self.assertFalse(manifest["boundaries_optimized"])
            self.assertFalse(manifest["new_primary_test_added"])

    def test_manuscript_tables_match_sources(self) -> None:
        _, rows = read_tsv(MAIN_TABLE_TSV)
        by_region = {row["region"]: row for row in rows}
        self.assertEqual(by_region["chr16_supergene"]["mean_q_haplotype"], "0.811443799942942")
        self.assertEqual(by_region["chr16_outside"]["species_dominant"], "44")

    def test_repeated_outputs_byte_stable_where_practical(self) -> None:
        run_stage5()
        outputs = [ROBUSTNESS_SUMMARY, MAIN_TABLE_TSV, MAIN_TABLE_MD, MAIN_CAPTION, RESULTS_TEXT, METHODS_TEXT, FINAL_ANALYSIS, MAIN_FIG_PDF, MAIN_FIG_PNG, SUPP_FIG_PDF, SUPP_FIG_PNG, MANIFEST]
        first = {path: sha256(path) for path in outputs}
        run_stage5()
        second = {path: sha256(path) for path in outputs}
        self.assertEqual(first, second)


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage5Tests)
    return unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args(argv)
    summary = run_stage5()
    if args.run_tests and not run_tests():
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
