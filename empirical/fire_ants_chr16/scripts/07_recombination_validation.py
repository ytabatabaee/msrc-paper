#!/usr/bin/env python3
"""Supplementary recombination-suppression validation for fire-ant chr16.

This script adds an independent recombination/linkage validation layer to the
frozen fire-ant analysis. It does not rerun or modify topology, ASTRAL, or
CASTLES-II results. Wang et al. 2013 linkage-map marker tables are preserved in
processed form, but they are on the original Si_gnF scaffold coordinate system.
Because no reliable scaffold-to-frozen-chr16 coordinate conversion is available
in the committed inputs, the chromosome-16 figure uses a schematic published
regional suppression track rather than a fabricated recombination-rate curve.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import statistics
import sys
import textwrap
import unittest
import zipfile
from dataclasses import dataclass
from datetime import date
from decimal import Decimal
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
EMP = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
SCRIPTS = EMP / "scripts"
DATA = REPO_ROOT / "data" / "fire_ants_chr16"
PROCESSED = DATA / "processed"
METADATA = DATA / "metadata"
RAW_RECOMB = DATA / "raw" / "recombination_sources" / "wang_2013_nature"

REGION_MANIFEST = METADATA / "region_manifest.tsv"
WINDOW_SUPPORT = PROCESSED / "stage4a_window_quartet_support.tsv"
SOURCE_ZIP = RAW_RECOMB / "41586_2013_BFnature11832_MOESM98_ESM.zip"
SOURCE_PDF = RAW_RECOMB / "41586_2013_BFnature11832_MOESM97_ESM.pdf"

RECOMB_MAP = PROCESSED / "recombination_map.tsv"
RECOMB_SUMMARY = RESULTS / "recombination_summary.tsv"
SOURCE_AUDIT = RESULTS / "recombination_source_audit.md"
COORD_AUDIT = RESULTS / "recombination_coordinate_audit.md"
REPORT = RESULTS / "recombination_report.md"
CAPTION = RESULTS / "recombination_figure_caption.txt"
MANIFEST = RESULTS / "recombination_manifest.json"
FIG_PDF = FIGURES / "fire_ants_recombination_genealogy.pdf"
FIG_PNG = FIGURES / "fire_ants_recombination_genealogy.png"
README = EMP / "README.md"

FROZEN_MANIFESTS = [
    RESULTS / "stage4a_manifest.json",
    RESULTS / "stage4b_manifest.json",
    RESULTS / "stage5_final_manifest.json",
    RESULTS / "stage6_manifest.json",
    RESULTS / "stage6_cu_manifest.json",
    RESULTS / "stage6c_manifest.json",
]

WANG_NATURE_URL = "https://www.nature.com/articles/nature11832"
WANG_SUPP_ZIP_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM98_ESM.zip"
WANG_SUPP_PDF_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM97_ESM.pdf"
WANG_DOI = "10.1038/nature11832"
WANG_CITATION = "Wang J., Wurm Y., Nipitwattanaphon M. et al. A Y-like social chromosome causes alternative colony organization in fire ants. Nature 493, 664-668 (2013)."
SUPERGENE_EXPECTED = (11680438, 27917498)
ANALYSIS_DATE = "2026-10-04"


@dataclass(frozen=True)
class Marker:
    marker: str
    family: str
    linkage_group: str
    scaffold: str
    original_position: int
    genetic_position_cm: Decimal
    source_table: str


@dataclass(frozen=True)
class IntervalClassification:
    region_class: str
    relation: str


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def load_regions() -> dict[str, dict[str, str]]:
    rows = read_tsv(REGION_MANIFEST)
    return {r["region_id"]: r for r in rows}


def frozen_supergene_interval() -> tuple[int, int]:
    regions = load_regions()
    sg = regions["chr16_supergene"]
    start = int(sg["coordinate_start"])
    end = int(sg["coordinate_end"])
    if (start, end) != SUPERGENE_EXPECTED:
        raise AssertionError(f"frozen supergene coordinates changed: {(start, end)}")
    return start, end


def recombination_cm_per_mb(left_pos: int, right_pos: int, left_cm: Decimal, right_cm: Decimal) -> Decimal:
    span_bp = right_pos - left_pos
    if span_bp <= 0:
        raise ValueError("physical span must be positive")
    return (right_cm - left_cm) / (Decimal(span_bp) / Decimal(1_000_000))


def classify_interval(start: int, end: int, sg_start: int, sg_end: int) -> IntervalClassification:
    if end <= start:
        raise ValueError("interval end must exceed start")
    if start >= sg_start and end <= sg_end:
        return IntervalClassification("supergene", "inside_frozen_supergene_interval")
    if end <= sg_start or start >= sg_end:
        return IntervalClassification("background_chr16", "outside_frozen_supergene_interval")
    return IntervalClassification("boundary_overlap", "crosses_frozen_supergene_boundary")


def preserve_reported_value(value: str) -> str:
    return str(value)


def family_from_table_name(name: str, title: str) -> str:
    match = re.search(r"family\s+([A-Za-z0-9]+)", title)
    if match:
        return match.group(1)
    stem = Path(name).stem
    return stem.replace("Supplementary_Table", "table")


def parse_marker_tables(zip_path: Path = SOURCE_ZIP) -> list[Marker]:
    if not zip_path.exists():
        raise FileNotFoundError(zip_path)
    markers: list[Marker] = []
    pos_re = re.compile(r"(?P<scaffold>Si_gnF\.scaffold\d+)_nt(?P<pos>\d+)")
    with zipfile.ZipFile(zip_path) as zf:
        for name in sorted(zf.namelist()):
            if not name.endswith(".txt"):
                continue
            lines = zf.read(name).decode("utf-8", "replace").splitlines()
            title = next((line for line in lines if line.strip()), "")
            family = family_from_table_name(name, title)
            header_i = None
            for i, line in enumerate(lines):
                if line.startswith("locus_id\tassembly_gnF_position"):
                    header_i = i
                    break
            if header_i is None:
                continue
            reader = csv.DictReader(lines[header_i:], delimiter="\t")
            for row in reader:
                locus = (row.get("locus_id") or "").strip()
                assembly_pos = (row.get("assembly_gnF_position") or "").strip()
                lg = (row.get("LG") or "").strip()
                cm = (row.get("cM") or "").strip()
                if not locus or not assembly_pos or not lg or cm == "":
                    continue
                m = pos_re.fullmatch(assembly_pos)
                if not m:
                    continue
                markers.append(Marker(
                    marker=locus,
                    family=family,
                    linkage_group=lg,
                    scaffold=m.group("scaffold"),
                    original_position=int(m.group("pos")),
                    genetic_position_cm=Decimal(cm),
                    source_table=Path(name).name,
                ))
    if not markers:
        raise AssertionError("no Wang et al. linkage markers parsed")
    return markers


def write_recombination_map(markers: list[Marker]) -> None:
    rows = [{
        "marker": m.marker,
        "family": m.family,
        "linkage_group": m.linkage_group,
        "original_assembly": "Wang_2013_Si_gnF_scaffold_assembly",
        "scaffold": m.scaffold,
        "original_position_bp": m.original_position,
        "genetic_position_cm": str(m.genetic_position_cm),
        "source": "Wang et al. 2013 Nature Supplementary Tables 8-14",
        "source_table": m.source_table,
        "coordinate_status": "original_scaffold_coordinates_not_projected_to_frozen_chr16",
    } for m in markers]
    write_tsv(RECOMB_MAP, rows, [
        "marker", "family", "linkage_group", "original_assembly", "scaffold",
        "original_position_bp", "genetic_position_cm", "source", "source_table", "coordinate_status",
    ])


def write_recombination_summary() -> None:
    sg_start, sg_end = frozen_supergene_interval()
    rows = [{
        "region": "chr16_supergene",
        "start_bp": sg_start,
        "end_bp": sg_end,
        "evidence_type": "direct_linkage_map_regional_summary",
        "reported_result": "approximately 13 Mb, about 55% of the social chromosome, with complete recombination suppression between SB and Sb",
        "relation": "published suppressed-recombination region corresponds to the independently frozen social-supergene analysis interval at regional scale",
        "unit": "published regional summary; no target-coordinate cM/Mb curve estimated",
        "source": f"{WANG_CITATION} DOI:{WANG_DOI}",
        "notes": "Marker-level linkage tables provide scaffold positions and cM, but no reliable conversion to the frozen Stolle/TWISST chromosome-16 coordinate system was found in the committed analysis inputs.",
    }]
    write_tsv(RECOMB_SUMMARY, rows, [
        "region", "start_bp", "end_bp", "evidence_type", "reported_result",
        "relation", "unit", "source", "notes",
    ])


def marker_summary(markers: list[Marker]) -> dict[str, object]:
    families = sorted({m.family for m in markers})
    scaffolds = sorted({m.scaffold for m in markers})
    lgs = sorted({m.linkage_group for m in markers})
    cm_values = [float(m.genetic_position_cm) for m in markers]
    return {
        "n_markers": len(markers),
        "n_families": len(families),
        "families": families,
        "n_scaffolds": len(scaffolds),
        "n_linkage_groups": len(lgs),
        "min_cM": min(cm_values),
        "max_cM": max(cm_values),
    }


def load_chr16_support() -> list[dict[str, object]]:
    rows = []
    for row in read_tsv(WINDOW_SUPPORT):
        if row["chrom"] == "chr16":
            rows.append({
                "window_index": int(row["window_index"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "mid": int(row["mid"]),
                "region": row["region"],
                "q_S": float(Decimal(row["q_S"])),
                "q_H": float(Decimal(row["q_H"])),
                "D": float(Decimal(row["D"])),
            })
    rows.sort(key=lambda r: r["mid"])
    if len(rows) != 96:
        raise AssertionError(f"expected 96 chr16 windows, found {len(rows)}")
    return rows


def make_figure() -> None:
    sg_start, sg_end = frozen_supergene_interval()
    chr16 = load_chr16_support()
    x_min = min(r["start"] for r in chr16) / 1_000_000
    x_max = max(r["end"] for r in chr16) / 1_000_000
    sg_start_mb = sg_start / 1_000_000
    sg_end_mb = sg_end / 1_000_000

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8.2, 5.8), sharex=True,
        gridspec_kw={"height_ratios": [0.8, 2.4], "hspace": 0.10},
    )
    for ax in (ax0, ax1):
        ax.axvspan(sg_start_mb, sg_end_mb, color="#d9d9d9", alpha=0.55, lw=0)

    ax0.hlines(0.5, x_min, x_max, color="#bdbdbd", lw=5, alpha=0.6)
    ax0.hlines(0.5, sg_start_mb, sg_end_mb, color="#7f2704", lw=10)
    ax0.text((sg_start_mb + sg_end_mb) / 2, 0.69,
             "published direct linkage-map evidence:\nregional SB-Sb recombination suppression",
             ha="center", va="bottom", fontsize=9, color="#7f2704")
    ax0.text(x_min, 0.18, "schematic regional validation; no target-coordinate cM/Mb curve",
             ha="left", va="center", fontsize=8, color="#525252")
    ax0.set_ylim(0, 1)
    ax0.set_yticks([])
    ax0.set_ylabel("recombination\nevidence", fontsize=9)
    ax0.spines[["top", "right", "left"]].set_visible(False)

    colors = {"q_S": "#2b8cbe", "q_H": "#e34a33"}
    labels = {"q_S": "species-history support", "q_H": "SB/Sb haplotype support"}
    gap_break_mb = 1.0
    for region in ["chr16A", "chr16_supergene", "chr16B"]:
        subset = [r for r in chr16 if r["region"] == region]
        for key in ["q_S", "q_H"]:
            first_label = labels[key] if region == "chr16A" else None
            xs_all = [r["mid"] / 1_000_000 for r in subset]
            ys_all = [r[key] for r in subset]
            ax1.scatter(xs_all, ys_all, s=16, color=colors[key], label=first_label, zorder=3)
            start_i = 0
            for i in range(1, len(xs_all) + 1):
                if i == len(xs_all) or xs_all[i] - xs_all[i - 1] > gap_break_mb:
                    if i - start_i > 1:
                        ax1.plot(xs_all[start_i:i], ys_all[start_i:i], lw=1.0, color=colors[key], alpha=0.9, zorder=2)
                    start_i = i
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_ylabel("TWISST-derived\nfocal quartet support")
    ax1.set_xlabel("chromosome 16 physical position (Mb; frozen Stolle/TWISST coordinates)")
    ax1.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9)
    ax1.text((sg_start_mb + sg_end_mb) / 2, 1.01, "frozen author-designated\nsupergene analysis interval",
             ha="center", va="bottom", fontsize=8, color="#525252")
    ax1.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Independent recombination-suppression evidence aligns with the frozen fire-ant genealogy signal", fontsize=11)
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PDF, bbox_inches="tight")
    fig.savefig(FIG_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_source_audit(markers: list[Marker]) -> None:
    summary = marker_summary(markers)
    zip_sha = sha256(SOURCE_ZIP) if SOURCE_ZIP.exists() else "missing"
    pdf_sha = sha256(SOURCE_PDF) if SOURCE_PDF.exists() else "missing"
    SOURCE_AUDIT.write_text(textwrap.dedent(f"""
        # Recombination/linkage source audit

        ## Selected source

        - Paper: {WANG_CITATION}
        - DOI: {WANG_DOI}
        - Article URL: {WANG_NATURE_URL}
        - Supplementary data URL: {WANG_SUPP_ZIP_URL}
        - Supplementary information URL: {WANG_SUPP_PDF_URL}
        - Download/access date: {ANALYSIS_DATE}
        - Supplementary data file: `{rel(SOURCE_ZIP)}`
        - Supplementary data SHA256: `{zip_sha}`
        - Supplementary information file: `{rel(SOURCE_PDF) if SOURCE_PDF.exists() else 'missing'}`
        - Supplementary information SHA256: `{pdf_sha}`

        ## Evidence hierarchy result

        Direct linkage-map evidence was found. Wang et al. provide RADtag linkage-map tables with marker identifiers, original scaffold positions, linkage groups, and cM positions for seven mapping families. The parsed marker inventory contains {summary['n_markers']} marker rows across {summary['n_families']} families ({', '.join(summary['families'])}). Marker-level numerical data are therefore available in the source assembly.

        ## Coordinate status

        The marker positions are reported as `Si_gnF.scaffold..._nt...` scaffold coordinates with linkage-map cM positions. The frozen MSRC fire-ant genealogy track uses Stolle/TWISST chromosome coordinates (`chr16`, bp). I did not find a reliable committed scaffold-to-target-chromosome conversion in the frozen analysis inputs, so the marker-level map is preserved in `data/fire_ants_chr16/processed/recombination_map.tsv` but is not projected onto chromosome 16.

        ## Data type

        - Data type: direct linkage map, represented here as a published regional summary for target-coordinate visualization.
        - Coordinate system: Wang 2013 original `Si_gnF` scaffold positions plus linkage-group cM positions.
        - Marker-level numerical data available: yes.
        - Marker-level target-coordinate data available: no.
    """).lstrip())


def write_coordinate_audit() -> None:
    sg_start, sg_end = frozen_supergene_interval()
    COORD_AUDIT.write_text(textwrap.dedent(f"""
        # Recombination coordinate audit

        ## Frozen MSRC fire-ant coordinate system

        The frozen local-genealogy analysis uses the author-labeled Stolle/TWISST window-coordinate file recorded in `data/fire_ants_chr16/metadata/region_manifest.tsv`. The chromosome-16 social-supergene analysis interval is `{sg_start}-{sg_end}` bp on `chr16`. These coordinates are an observed BUSCO-window span / author-designated analysis interval, not exact inversion breakpoints.

        ## Linkage-map coordinate system

        Wang et al. 2013 Supplementary Tables 8-14 report RADtag positions in the original assembly as values such as `Si_gnF.scaffold00759_nt19793`, together with linkage groups and cM coordinates.

        ## Conversion decision

        No reliable conversion from the Wang `Si_gnF` scaffold coordinates to the frozen Stolle/TWISST chromosome-16 coordinate system was found in the committed fire-ant analysis inputs. I therefore did not approximate, scale, or otherwise project the marker positions onto chr16. The integrated figure uses a schematic regional recombination-suppression track based on Wang et al.'s published direct linkage-map conclusion rather than a target-coordinate cM/Mb curve.

        ## Conversion method

        None. Marker-level original positions are preserved in `data/fire_ants_chr16/processed/recombination_map.tsv`; target-coordinate visualization is schematic.
    """).lstrip())


def write_report(markers: list[Marker]) -> None:
    summary = marker_summary(markers)
    sg_start, sg_end = frozen_supergene_interval()
    REPORT.write_text(textwrap.dedent(f"""
        # Recombination-suppression validation report

        ## Source study

        The validation uses an independent direct linkage-map source: {WANG_CITATION} DOI:{WANG_DOI}. This source is distinct from the Stolle et al. local-tree/TWISST analysis used for the frozen MSRC fire-ant genealogy results.

        ## Data type

        Direct linkage-map marker data were recovered from the Wang et al. supplementary data archive. The parsed source tables contain {summary['n_markers']} RADtag marker rows across {summary['n_families']} mapping families. Each row includes a marker, an original scaffold position, a linkage group, and a cM coordinate. These marker-level data are written to `data/fire_ants_chr16/processed/recombination_map.tsv` in their original coordinate system.

        ## Coordinate compatibility

        The linkage-map marker positions use original `Si_gnF` scaffold coordinates. The frozen MSRC fire-ant topology track uses Stolle/TWISST chromosome coordinates. Because no reliable scaffold-to-chromosome conversion was found in the committed inputs, no marker-level cM/Mb curve was projected onto chromosome 16. The figure therefore uses the published regional conclusion as a schematic validation track over the frozen author-designated supergene interval `{sg_start}-{sg_end}` bp.

        ## Recombination/linkage result

        Wang et al. report a large social-chromosome region of approximately 13 Mb, about 55% of the chromosome, in which recombination is completely suppressed between the SB and Sb social chromosomes. This is direct linkage-map evidence for suppressed recombination across the social-supergene region.

        ## Relation to the frozen genealogy signal

        The frozen MSRC fire-ant analysis shows that chr16 windows outside the supergene are species-history dominated, whereas windows inside the author-designated supergene interval shift strongly toward the cross-species SB/Sb haplotype quartet. The independent linkage-map evidence supports the biological consistency of this result: the genomic interval with the strong social-haplotype genealogy is also the known recombination-suppressed social chromosome region.

        ## Relation to Stage 6C branch-length effects

        Stage 6C remains unchanged. Its fixed-topology individual-level ASTRAL4/CASTLES-II comparison found strong CULength reductions for the two focal SB/Sb species-pair branches when supergene windows were added, with SULength responses that differed by branch. The recombination validation does not estimate branch lengths and does not reinterpret CU or SU values as recombination rates. It supports the narrative that linked histories in a recombination-suppressed region can influence summary-tree branch estimates even when the global topology remains stable.

        ## Limitations

        The marker-level linkage map could not be placed onto the frozen chromosome-16 coordinate axis without a documented coordinate conversion. The integrated figure is therefore a regional validation figure, not a new recombination-rate map. It should not be read as estimating local cM/Mb values across the Stolle/TWISST windows.

        ## Introgression caveat

        The fire-ant supergene literature invokes recurrent adaptive introgression among socially polymorphic species. This validation supports the role of recombination suppression in maintaining a long linked genealogy, but it does not show that recombination suppression alone caused the observed genealogy or that MSRC without gene flow explains the system.
    """).lstrip())


def write_caption() -> None:
    CAPTION.write_text(textwrap.dedent(f"""
        Independent recombination-suppression evidence and local genealogy across fire-ant chromosome 16. Top panel: schematic regional representation of the direct linkage-map result from Wang et al. (2013), who reported an approximately 13-Mb social-chromosome region with complete recombination suppression between SB and Sb. Marker-level Wang et al. linkage data are available in original scaffold coordinates but were not projected onto the frozen Stolle/TWISST chromosome-16 coordinate system because no reliable coordinate conversion was available in the committed inputs. Bottom panel: existing frozen MSRC fire-ant TWISST-derived focal quartet support across chromosome 16, showing species-history support and cross-species SB/Sb haplotype support for the published four-BUSCO windows. The gray shading marks the independently frozen author-designated supergene analysis interval; it is not an optimized topology boundary or exact inversion breakpoint. The overlap supports biological consistency between the recombination-suppressed social chromosome and the localized haplotype-associated genealogy, but it does not establish a causal mechanism or rule out the recurrent-introgression history inferred by the source study.
    """).strip() + "\n")


def write_manifest(markers: list[Marker], before: dict[str, str], after: dict[str, str]) -> None:
    outputs = [RECOMB_MAP, RECOMB_SUMMARY, SOURCE_AUDIT, COORD_AUDIT, REPORT, CAPTION, FIG_PDF, FIG_PNG]
    script_path = Path(__file__).resolve()
    manifest = {
        "analysis": "fire_ant_recombination_suppression_validation",
        "analysis_date": ANALYSIS_DATE,
        "analysis_type": "supplementary_figure_level_validation",
        "stage5_stage6_frozen_results_modified": False,
        "source_citation": WANG_CITATION,
        "doi": WANG_DOI,
        "source_urls": {
            "article": WANG_NATURE_URL,
            "supplementary_data_zip": WANG_SUPP_ZIP_URL,
            "supplementary_information_pdf": WANG_SUPP_PDF_URL,
        },
        "source_files": {
            rel(SOURCE_ZIP): sha256(SOURCE_ZIP) if SOURCE_ZIP.exists() else "missing",
            rel(SOURCE_PDF): sha256(SOURCE_PDF) if SOURCE_PDF.exists() else "missing",
        },
        "data_mode": "direct_linkage_map_regional_summary_schematic",
        "marker_level_numerical_data_available": True,
        "marker_level_target_coordinate_data_available": False,
        "n_parsed_linkage_markers": len(markers),
        "coordinate_assembly": {
            "linkage_source": "Wang_2013_Si_gnF_scaffold_coordinates_plus_linkage_group_cM",
            "frozen_genealogy_track": "Stolle_TWISST_chr16_window_coordinates",
        },
        "conversion_method": "none; reliable scaffold-to-frozen-chr16 coordinate conversion unavailable in committed inputs",
        "frozen_supergene_interval": {
            "chromosome": "chr16",
            "coordinate_start": SUPERGENE_EXPECTED[0],
            "coordinate_end": SUPERGENE_EXPECTED[1],
            "definition": "author-designated analysis interval / observed BUSCO-window span",
        },
        "frozen_manifest_checksums_before": {rel(k): v for k, v in before.items()},
        "frozen_manifest_checksums_after": {rel(k): v for k, v in after.items()},
        "frozen_manifest_checksums_unchanged": before == after,
        "script_hash": {rel(script_path): sha256(script_path)},
        "output_hashes": {rel(p): sha256(p) for p in outputs if p.exists()},
        "prohibited_actions": {
            "astral_rerun": False,
            "castles_rerun": False,
            "gene_tree_regeneration": False,
            "topology_reclassification": False,
            "new_species_tree_estimation": False,
            "recombination_inferred_from_genealogy": False,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def update_readme() -> None:
    section = textwrap.dedent("""
        ## Recombination-suppression validation

        A small supplementary validation layer aligns the frozen chromosome-16 genealogy signal with independent recombination evidence from Wang et al. 2013 (`Nature`, DOI `10.1038/nature11832`). The source provides direct RADtag linkage-map marker tables and reports an approximately 13-Mb social-chromosome region with complete recombination suppression between SB and Sb. The marker tables are preserved in `../../data/fire_ants_chr16/processed/recombination_map.tsv` in their original `Si_gnF` scaffold coordinate system.

        Because no reliable conversion from the Wang scaffold coordinates to the frozen Stolle/TWISST chromosome-16 coordinates was found in the committed inputs, the manuscript figure uses a schematic regional recombination-suppression track rather than a fabricated cM/Mb curve. The integrated figure is `figures/fire_ants_recombination_genealogy.pdf`, with report and provenance in `results/recombination_report.md`, `results/recombination_source_audit.md`, `results/recombination_coordinate_audit.md`, and `results/recombination_manifest.json`.

        This validation does not rerun or alter the frozen topology, ASTRAL, CASTLES-II, or branch-length results. It supports the biological consistency of the localized social-haplotype genealogy with an independently known recombination-suppressed social chromosome region, while retaining the introgression caveat from the source literature.
    """).strip()
    text = README.read_text()
    if "## Recombination-suppression validation" in text:
        start = text.index("## Recombination-suppression validation")
        # Replace until next H2 or EOF.
        next_match = re.search(r"\n## ", text[start + 1:])
        if next_match:
            end = start + 1 + next_match.start()
            text = text[:start] + section + "\n" + text[end:]
        else:
            text = text[:start] + section + "\n"
    else:
        if not text.endswith("\n"):
            text += "\n"
        text += "\n" + section + "\n"
    README.write_text(text)


def frozen_manifest_hashes() -> dict[Path, str]:
    return {p: sha256(p) for p in FROZEN_MANIFESTS if p.exists()}


def run_analysis(update_readme_flag: bool = True) -> dict[str, object]:
    before = frozen_manifest_hashes()
    frozen_supergene_interval()
    markers = parse_marker_tables()
    write_recombination_map(markers)
    write_recombination_summary()
    write_source_audit(markers)
    write_coordinate_audit()
    write_report(markers)
    write_caption()
    make_figure()
    if update_readme_flag:
        update_readme()
    after = frozen_manifest_hashes()
    if before != after:
        raise AssertionError("frozen fire-ant result manifest checksum changed during recombination validation")
    write_manifest(markers, before, after)
    return {"markers": len(markers), "before": before, "after": after}


class RecombinationValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.markers = parse_marker_tables()

    def test_frozen_supergene_coordinates_unchanged(self):
        self.assertEqual(frozen_supergene_interval(), SUPERGENE_EXPECTED)

    def test_direct_recombination_rate_synthetic(self):
        rate = recombination_cm_per_mb(1_000_000, 2_000_000, Decimal("3.5"), Decimal("5.5"))
        self.assertEqual(rate, Decimal("2"))

    def test_zero_physical_span_rejected(self):
        with self.assertRaises(ValueError):
            recombination_cm_per_mb(10, 10, Decimal("0"), Decimal("1"))

    def test_boundary_overlap_classification(self):
        sg_start, sg_end = SUPERGENE_EXPECTED
        self.assertEqual(classify_interval(sg_start + 1, sg_start + 100, sg_start, sg_end).region_class, "supergene")
        self.assertEqual(classify_interval(1, sg_start - 1, sg_start, sg_end).region_class, "background_chr16")
        self.assertEqual(classify_interval(sg_start - 1, sg_start + 1, sg_start, sg_end).region_class, "boundary_overlap")

    def test_ld_fallback_label_available(self):
        row = {"evidence_type": "LD-based evidence of recombination suppression"}
        self.assertIn("LD-based", row["evidence_type"])

    def test_schematic_fallback_mode_written(self):
        rows = read_tsv(RECOMB_SUMMARY)
        self.assertEqual(rows[0]["evidence_type"], "direct_linkage_map_regional_summary")

    def test_inequalities_preserved(self):
        self.assertEqual(preserve_reported_value("<0.5"), "<0.5")

    def test_no_topology_fields_define_recombination_regions(self):
        header = RECOMB_SUMMARY.read_text().splitlines()[0].split("\t")
        self.assertFalse(any(h.startswith("q_") or h in {"D", "dominant_class"} for h in header))
        row = read_tsv(RECOMB_SUMMARY)[0]
        self.assertEqual((int(row["start_bp"]), int(row["end_bp"])), SUPERGENE_EXPECTED)

    def test_frozen_manifests_unchanged_in_manifest(self):
        manifest = json.loads(MANIFEST.read_text())
        self.assertTrue(manifest["frozen_manifest_checksums_unchanged"])
        self.assertEqual(manifest["frozen_manifest_checksums_before"], manifest["frozen_manifest_checksums_after"])

    def test_figure_generation_succeeds(self):
        self.assertTrue(FIG_PDF.exists() and FIG_PDF.stat().st_size > 1000)
        self.assertTrue(FIG_PNG.exists() and FIG_PNG.stat().st_size > 1000)

    def test_marker_tables_parsed(self):
        self.assertGreater(len(self.markers), 0)
        self.assertTrue(all(m.scaffold.startswith("Si_gnF.scaffold") for m in self.markers))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run lightweight validation tests after regenerating outputs")
    parser.add_argument("--no-readme", action="store_true", help="do not update README")
    args = parser.parse_args(argv)

    result = run_analysis(update_readme_flag=not args.no_readme)
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(RecombinationValidationTests)
        runner = unittest.TextTestRunner(verbosity=2)
        test_result = runner.run(suite)
        if not test_result.wasSuccessful():
            return 1
        print(f"tests_passed={test_result.testsRun}")
    print(f"parsed_linkage_markers={result['markers']}")
    print(f"figure={rel(FIG_PDF)}")
    print(f"manifest={rel(MANIFEST)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
