#!/usr/bin/env python3
"""Supplementary recombination-suppression validation audit for fire-ant chr16.

The preferred deliverable is a quantitative physical-coordinate recombination or
LD track aligned to the frozen chromosome-16 genealogy track. This script first
preserves the recovered Wang et al. 2013 direct linkage-map marker data, then
audits whether either Route A (documented Wang Si_gnF scaffold -> Stolle Si_gnGA
chr16 placement) or Route B (reproducible Yan et al. 2020 physical-coordinate LD
values/genotypes) is available. If neither route is available, it records the
explicit stop condition and archives, but does not promote, the earlier schematic
figure.

It does not rerun or modify frozen topology, TWISST, ASTRAL, CASTLES-II, or
Stage 4A-6C outputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import textwrap
import unittest
import zipfile
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")
import matplotlib
matplotlib.use("Agg")

REPO_ROOT = Path(__file__).resolve().parents[3]
EMP = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
DATA = REPO_ROOT / "data" / "fire_ants_chr16"
PROCESSED = DATA / "processed"
METADATA = DATA / "metadata"
RAW_RECOMB = DATA / "raw" / "recombination_sources"
RAW_WANG = RAW_RECOMB / "wang_2013_nature"
RAW_YAN = RAW_RECOMB / "yan_2020"

REGION_MANIFEST = METADATA / "region_manifest.tsv"
WINDOW_SUPPORT = PROCESSED / "stage4a_window_quartet_support.tsv"
WANG_SOURCE_ZIP = RAW_WANG / "41586_2013_BFnature11832_MOESM98_ESM.zip"
WANG_SOURCE_PDF = RAW_WANG / "41586_2013_BFnature11832_MOESM97_ESM.pdf"
YAN_SUPP_PDF = RAW_YAN / "41559_2019_1081_MOESM1_ESM.pdf"
YAN_SUPP_XLS = RAW_YAN / "41559_2019_1081_MOESM2_ESM.xls"
YAN_BAD_XLSX = RAW_YAN / "41559_2019_1081_MOESM3_ESM.xlsx"

RECOMB_MAP = PROCESSED / "recombination_map.tsv"
WANG_TO_STOLLE_MAP = PROCESSED / "wang_to_stolle_coordinate_map.tsv"
RECOMB_MAP_CHR16 = PROCESSED / "recombination_map_chr16.tsv"
RECOMB_INTERVALS_CHR16 = PROCESSED / "recombination_intervals_chr16.tsv"
RECOMB_BINNED_CHR16 = PROCESSED / "recombination_binned_chr16.tsv"
YAN_LD = PROCESSED / "yan2020_chr16_ld.tsv"
YAN_LD_BINNED = PROCESSED / "yan2020_chr16_ld_binned.tsv"

RECOMB_SUMMARY = RESULTS / "recombination_summary.tsv"
SOURCE_AUDIT = RESULTS / "recombination_source_audit.md"
COORD_AUDIT = RESULTS / "recombination_coordinate_audit.md"
REPORT = RESULTS / "recombination_report.md"
CAPTION = RESULTS / "recombination_figure_caption.txt"
MANIFEST = RESULTS / "recombination_manifest.json"
README = EMP / "README.md"
FIG_PDF = FIGURES / "fire_ants_recombination_genealogy.pdf"
FIG_PNG = FIGURES / "fire_ants_recombination_genealogy.png"
ARCHIVE_PDF = FIGURES / "fire_ants_recombination_genealogy_schematic_archived.pdf"
ARCHIVE_PNG = FIGURES / "fire_ants_recombination_genealogy_schematic_archived.png"

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
YAN_URL = "https://www.nature.com/articles/s41559-019-1081-1"
YAN_SUPP_XLS_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-019-1081-1/MediaObjects/41559_2019_1081_MOESM2_ESM.xls"
YAN_SUPP_PDF_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-019-1081-1/MediaObjects/41559_2019_1081_MOESM1_ESM.pdf"
YAN_DOI = "10.1038/s41559-019-1081-1"
YAN_CITATION = "Yan Z., Martin S.H., Gotzek D. et al. Evolution of a supergene that regulates a trans-species social polymorphism. Nature Ecology & Evolution 4, 240-249 (2020)."
STOLLE_VCF_URL = "https://github.com/wurmlab/wurmlab.github.io/raw/master/data/supergene_introgression/gt.vcf.gz/gt.vcf.gz"
STOLLE_VCF_SIZE_BYTES = 1503782377
SUPERGENE_EXPECTED = (11680438, 27917498)
YAN_INVERSIONS = [
    {"inversion": "In(16)1", "start_bp": 14549064, "end_bp": 24031576},
    {"inversion": "In(16)2", "start_bp": 13705210, "end_bp": 24030990},
    {"inversion": "In(16)3", "start_bp": 12612565, "end_bp": 13683100},
]
ANALYSIS_DATE = "2026-10-04"
STATUS = "QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED"


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
class ScaffoldPlacement:
    scaffold: str
    scaffold_start: int
    scaffold_end: int
    chromosome: str
    chr_start: int
    chr_end: int
    orientation: str


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
    return {r["region_id"]: r for r in read_tsv(REGION_MANIFEST)}


def frozen_supergene_interval() -> tuple[int, int]:
    sg = load_regions()["chr16_supergene"]
    start = int(sg["coordinate_start"])
    end = int(sg["coordinate_end"])
    if (start, end) != SUPERGENE_EXPECTED:
        raise AssertionError(f"frozen supergene coordinates changed: {(start, end)}")
    return start, end


def recombination_cm_per_mb(left_pos: int, right_pos: int, left_cm: Decimal, right_cm: Decimal) -> Decimal:
    span_bp = abs(right_pos - left_pos)
    if span_bp <= 0:
        raise ValueError("physical span must be positive")
    return abs(right_cm - left_cm) / (Decimal(span_bp) / Decimal(1_000_000))


def convert_scaffold_position(position: int, placement: ScaffoldPlacement) -> int:
    if placement.orientation not in {"+", "-", "?"}:
        raise ValueError("orientation must be +, -, or ?")
    if not (placement.scaffold_start <= position <= placement.scaffold_end):
        raise ValueError("position outside scaffold placement")
    if placement.orientation == "-":
        return placement.chr_start + placement.scaffold_end - position
    return placement.chr_start - placement.scaffold_start + position


def calculate_ld_r2(genotypes_a: Iterable[int], genotypes_b: Iterable[int]) -> float:
    a = list(genotypes_a)
    b = list(genotypes_b)
    if len(a) != len(b) or len(a) < 2:
        raise ValueError("equal genotype vectors with at least two samples are required")
    pairs = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    if len(pairs) < 2:
        return math.nan
    xs = [float(x) for x, _ in pairs]
    ys = [float(y) for _, y in pairs]
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    r = cov / math.sqrt(vx * vy)
    return r * r


def family_from_table_name(name: str, title: str) -> str:
    match = re.search(r"family\s+([A-Za-z0-9]+)", title)
    if match:
        return match.group(1)
    return Path(name).stem.replace("Supplementary_Table", "table")


def parse_marker_tables(zip_path: Path = WANG_SOURCE_ZIP) -> list[Marker]:
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
            header_i = next((i for i, line in enumerate(lines) if line.startswith("locus_id\tassembly_gnF_position")), None)
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
                markers.append(Marker(locus, family, lg, m.group("scaffold"), int(m.group("pos")), Decimal(cm), Path(name).name))
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
    write_tsv(RECOMB_MAP, rows, ["marker", "family", "linkage_group", "original_assembly", "scaffold", "original_position_bp", "genetic_position_cm", "source", "source_table", "coordinate_status"])


def marker_summary(markers: list[Marker]) -> dict[str, object]:
    return {
        "n_markers": len(markers),
        "n_families": len({m.family for m in markers}),
        "families": sorted({m.family for m in markers}),
        "n_scaffolds": len({m.scaffold for m in markers}),
        "n_linkage_groups": len({m.linkage_group for m in markers}),
    }


def yan_supplement_contains_ld_track() -> bool:
    if not YAN_SUPP_XLS.exists():
        return False
    # The workbook strings show Supplementary Tables 1-5, gene lists, expression,
    # and metadata, but no numeric Extended Data Fig. 5 LD matrix/track source.
    text = ""
    try:
        import subprocess
        cp = subprocess.run(["strings", str(YAN_SUPP_XLS)], check=False, capture_output=True, text=True, timeout=10)
        text = cp.stdout.lower()
    except Exception:
        return False
    indicators = ["extended data fig. 5", "linkage disequilibrium", "ld dot plot", "r2"]
    return all(ind in text for ind in indicators)


def archive_existing_schematic() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    if FIG_PDF.exists() and not ARCHIVE_PDF.exists():
        shutil.copy2(FIG_PDF, ARCHIVE_PDF)
    if FIG_PNG.exists() and not ARCHIVE_PNG.exists():
        shutil.copy2(FIG_PNG, ARCHIVE_PNG)


def write_recombination_summary(markers: list[Marker]) -> None:
    sg_start, sg_end = frozen_supergene_interval()
    inv_start = min(r["start_bp"] for r in YAN_INVERSIONS)
    inv_end = max(r["end_bp"] for r in YAN_INVERSIONS)
    rows = [
        {
            "analysis_item": "quantitative_track_status",
            "region": "chr16",
            "start_bp": "NA",
            "end_bp": "NA",
            "evidence_type": "route_a_or_b_quantitative_track",
            "metric": STATUS,
            "value": "NA",
            "unit": "NA",
            "source": "route audit",
            "notes": "No documented Wang scaffold-to-Stolle chr16 mapping and no reproducible Yan physical-coordinate LD source values were recovered; no quantitative recombination/LD figure was generated.",
        },
        {
            "analysis_item": "wang_direct_linkage_marker_inventory",
            "region": "source_original_scaffolds",
            "start_bp": "NA",
            "end_bp": "NA",
            "evidence_type": "direct_linkage_map_marker_tables",
            "metric": "parsed_marker_rows",
            "value": len(markers),
            "unit": "markers",
            "source": f"{WANG_CITATION} DOI:{WANG_DOI}",
            "notes": "Markers have Si_gnF scaffold positions and cM coordinates, but are not projected to frozen chr16 without a documented scaffold placement table.",
        },
        {
            "analysis_item": "yan_inversion_union",
            "region": "published_inversion_union",
            "start_bp": inv_start,
            "end_bp": inv_end,
            "evidence_type": "published_physical_inversion_breakpoints",
            "metric": "inversion_union_span",
            "value": inv_end - inv_start,
            "unit": "bp",
            "source": f"{YAN_CITATION} DOI:{YAN_DOI}",
            "notes": "Published SB-reference inversion union only; not a replacement for the frozen MSRC chr16_supergene mask.",
        },
        {
            "analysis_item": "frozen_msrc_supergene_span",
            "region": "chr16_supergene",
            "start_bp": sg_start,
            "end_bp": sg_end,
            "evidence_type": "frozen_genealogy_analysis_interval",
            "metric": "analysis_span",
            "value": sg_end - sg_start,
            "unit": "bp",
            "source": "Stolle/TWISST window coordinates frozen in Stage 0",
            "notes": "Author-defined BUSCO-window analysis span, not exact inversion boundary.",
        },
    ]
    write_tsv(RECOMB_SUMMARY, rows, ["analysis_item", "region", "start_bp", "end_bp", "evidence_type", "metric", "value", "unit", "source", "notes"])


def write_source_audit(markers: list[Marker]) -> None:
    s = marker_summary(markers)
    def maybe_sha(p: Path) -> str:
        return sha256(p) if p.exists() and p.stat().st_size > 0 else "missing"
    SOURCE_AUDIT.write_text(textwrap.dedent(f"""
        # Recombination/linkage source audit

        ## Route A: Wang et al. 2013 direct linkage map

        - Paper: {WANG_CITATION}
        - DOI: {WANG_DOI}
        - Article URL: {WANG_NATURE_URL}
        - Supplementary data URL: {WANG_SUPP_ZIP_URL}
        - Supplementary data file: `{rel(WANG_SOURCE_ZIP)}`
        - Supplementary data SHA256: `{maybe_sha(WANG_SOURCE_ZIP)}`
        - Supplementary information URL: {WANG_SUPP_PDF_URL}
        - Supplementary information file: `{rel(WANG_SOURCE_PDF) if WANG_SOURCE_PDF.exists() else 'missing'}`
        - Supplementary information SHA256: `{maybe_sha(WANG_SOURCE_PDF)}`
        - Data type: direct linkage map.
        - Marker-level numerical data available: yes.
        - Parsed marker rows: {s['n_markers']} across {s['n_families']} families ({', '.join(s['families'])}).
        - Original coordinate system: `Si_gnF.scaffold..._nt...` plus family-specific linkage-group cM.

        A documented Wang `Si_gnF` scaffold to frozen Stolle `Si_gnGA`/`gng20170922wFex.fa` chromosome-16 placement table was not recovered. The upstream Stolle README mentions `linkage_map_supergene.txt` and `2018-05-11-linkage-map/results/linkage_map_supergene.txt`, but that file is not committed in the local frozen snapshot or in the public upstream GitHub tree inspected for this analysis. Therefore Route A did not produce target-coordinate cM/Mb values.

        ## Route B: Yan et al. 2020 physical-coordinate LD

        - Paper: {YAN_CITATION}
        - DOI: {YAN_DOI}
        - Article URL: {YAN_URL}
        - Supplementary table URL: {YAN_SUPP_XLS_URL}
        - Supplementary table file: `{rel(YAN_SUPP_XLS) if YAN_SUPP_XLS.exists() else 'missing'}`
        - Supplementary table SHA256: `{maybe_sha(YAN_SUPP_XLS)}`
        - Supplementary information URL: {YAN_SUPP_PDF_URL}
        - Supplementary information file: `{rel(YAN_SUPP_PDF) if YAN_SUPP_PDF.exists() else 'missing'}`
        - Supplementary information SHA256: `{maybe_sha(YAN_SUPP_PDF)}`
        - Public raw sequence BioProject: PRJNA421367.

        Yan et al. report physical-coordinate LD across chr16 in Extended Data Fig. 5 and provide exact SB-reference inversion breakpoints in the article. The downloadable Supplementary Tables 1-5 contain gene lists, expression tables, and sample metadata, but not the numeric LD matrix or a one-dimensional LD track underlying Extended Data Fig. 5. The BioProject contains large raw sequence data; no small indexed genotype/LD file sufficient to reconstruct chr16 LD was recovered. The Stolle/Wurmlab genome-wide VCF endpoint is approximately {STOLLE_VCF_SIZE_BYTES:,} bytes and is not the Yan Extended Data Fig. 5 source; it was not downloaded blindly for this figure.

        ## Evidence hierarchy conclusion

        `{STATUS}`. Wang 2013 remains direct experimental linkage support in original coordinates, and Yan 2020 remains published physical-coordinate inversion/LD evidence, but neither yielded a reproducible quantitative chr16 recombination/LD track aligned to the frozen MSRC coordinate axis.
    """).lstrip())


def write_coordinate_audit() -> None:
    sg_start, sg_end = frozen_supergene_interval()
    inv_lines = "\n".join(f"- {r['inversion']}: chr16:{r['start_bp']}-{r['end_bp']}" for r in YAN_INVERSIONS)
    content = f"""# Recombination coordinate audit

## Frozen MSRC / Stolle coordinate system

The frozen fire-ant local-genealogy analysis uses Stolle et al. 2022 workflow coordinates from the `Si_gnGA` / `gng20170922wFex.fa` chromosome-level reference. The frozen chromosome-16 supergene analysis interval is `chr16:{sg_start}-{sg_end}`. This is an author-defined BUSCO-window analysis span, not an exact inversion boundary.

## Wang et al. 2013 coordinate system

Wang linkage-map markers are reported on old `Si_gnF` scaffolds, for example `Si_gnF.scaffold00759_nt19793`, with family-specific genetic positions in cM. The Stolle upstream README states that a file named `linkage_map_supergene.txt` / `input/gngs_linkage_map.txt` should contain fields such as `scaffold`, `scaffold_start`, `scaffold_end`, `chr`, `chr_start`, `chr_end`, `orientation`, and possibly `region`, but this placement file was not recovered from the local frozen inputs or public upstream GitHub tree. No Wang marker was accepted as mapped to chr16.

## Yan et al. 2020 coordinate system

Yan et al. report physical chr16 coordinates on the S. invicta SB reference and provide these inversion breakpoints:

{inv_lines}

The union is chr16:{min(r['start_bp'] for r in YAN_INVERSIONS)}-{max(r['end_bp'] for r in YAN_INVERSIONS)}. These coordinates are useful biological annotations, but without the underlying LD values/genotypes or a documented coordinate equivalence/liftover to the Stolle `Si_gnGA` axis, they are not sufficient to generate a quantitative LD track aligned to the frozen MSRC windows.

## Overlay decision

Assembly compatibility checks did not pass for a quantitative overlay. No direct cM/Mb track and no LD r² track were overlaid on the frozen genealogy track. The previous schematic figure was archived and should not be treated as the requested quantitative validation figure.
"""
    COORD_AUDIT.write_text(content)

def write_report(markers: list[Marker]) -> None:
    s = marker_summary(markers)
    inv_start = min(r["start_bp"] for r in YAN_INVERSIONS)
    inv_end = max(r["end_bp"] for r in YAN_INVERSIONS)
    REPORT.write_text(textwrap.dedent(f"""
        # Recombination-suppression validation report

        ## Purpose

        The requested update was to replace the earlier schematic recombination bar with a quantitative chromosome-16 recombination or LD track on a real physical coordinate axis. The target layout was physical chr16 position, quantitative recombination/LD evidence, and the existing frozen `q_species` / `q_haplotype` genealogy track.

        ## Route A: direct Wang linkage map

        Wang et al. 2013 provide direct linkage-map marker tables. I parsed {s['n_markers']} RADtag marker rows across {s['n_families']} mapping families and preserved them in `data/fire_ants_chr16/processed/recombination_map.tsv`. These rows contain original `Si_gnF` scaffold positions and family-specific genetic positions in cM.

        Route A did not produce a chr16 cM/Mb track because no documented `Si_gnF` scaffold-to-Stolle `Si_gnGA`/`gng20170922wFex.fa` chromosome-placement file was recovered. The important candidate file named `linkage_map_supergene.txt` is referenced by the upstream README but was not present in the local frozen inputs or the public upstream repository tree. I did not guess scaffold placements.

        ## Route B: Yan physical-coordinate LD

        Yan et al. 2020 report LD r² across physical chr16 and exact SB-reference inversion breakpoints. The published inversion union is chr16:{inv_start}-{inv_end}. I downloaded the minimal journal supplementary PDF and Supplementary Tables 1-5 and inspected them. They do not contain the numeric Extended Data Fig. 5 LD matrix or a one-dimensional LD track. Reconstructing LD from raw PRJNA421367 reads would require large-scale read/genotype processing, and no small public chr16 genotype/LD file was recovered. I did not digitize the published heatmap or treat pixels as quantitative data.

        ## Result

        `{STATUS}`

        No quantitative recombination or LD track was generated. The previous constant schematic bar was copied to archive filenames for provenance, but it should not be used as the final quantitative validation figure requested here.

        ## Relation to frozen genealogy and Stage 6C

        The frozen topology/TWISST/ASTRAL/CASTLES-II results remain unchanged. The local-genealogy result still shows a pronounced switch from species-history support outside the frozen supergene analysis span to SB/Sb haplotype support inside it. Stage 6C branch-length results remain unchanged. This audit only addresses whether an independent quantitative recombination/LD track can be reproducibly aligned to those frozen coordinates.

        ## Introgression caveat

        The fire-ant supergene is known to have experienced recurrent adaptive introgression. Recombination suppression helps preserve a long linked haplotype after such events, but this audit does not imply that MSRC without gene flow fully explains the system.

        ## Smallest missing objects

        A quantitative figure would require one of the following: (1) a documented Wang `Si_gnF` scaffold-to-Stolle `Si_gnGA` chr16 placement table such as `linkage_map_supergene.txt`, or (2) the numeric Yan Extended Data Fig. 5 SNP/genotype/LD source data in physical chr16 coordinates, preferably as a chr16 VCF/genotype matrix or precomputed r² table.
    """).lstrip())


def write_caption() -> None:
    CAPTION.write_text(textwrap.dedent(f"""
        No manuscript-ready quantitative recombination/LD figure was generated. Direct Wang et al. linkage-map marker tables were recovered in original `Si_gnF` scaffold coordinates, but no documented scaffold-to-Stolle chr16 conversion was recovered. Yan et al. provide published physical-coordinate LD evidence and inversion breakpoints, but the downloadable supplementary tables do not include the numeric LD matrix or track underlying Extended Data Fig. 5. The earlier schematic recombination bar has been archived for provenance only and should not be cited as a quantitative recombination-rate or LD track. The frozen MSRC genealogy statistics are unchanged.
    """).strip() + "\n")


def write_manifest(markers: list[Marker], before: dict[Path, str], after: dict[Path, str]) -> None:
    existing_outputs = [p for p in [RECOMB_MAP, RECOMB_SUMMARY, SOURCE_AUDIT, COORD_AUDIT, REPORT, CAPTION, ARCHIVE_PDF, ARCHIVE_PNG] if p.exists()]
    source_files = [p for p in [WANG_SOURCE_ZIP, WANG_SOURCE_PDF, YAN_SUPP_PDF, YAN_SUPP_XLS] if p.exists()]
    manifest = {
        "analysis": "fire_ant_recombination_suppression_validation",
        "analysis_date": ANALYSIS_DATE,
        "status": STATUS,
        "analysis_type": "supplementary_quantitative_track_audit",
        "stage4a_to_stage6c_frozen_results_modified": False,
        "route_a_wang_scaffold_to_chr16_mapping_found": False,
        "route_a_wang_markers_parsed": len(markers),
        "route_a_wang_markers_mapped_to_chr16": 0,
        "route_a_cm_per_mb_calculated": False,
        "route_b_needed": True,
        "route_b_yan_source_used": {
            "article": YAN_URL,
            "supplementary_pdf": YAN_SUPP_PDF_URL,
            "supplementary_xls": YAN_SUPP_XLS_URL,
            "bioproject": "PRJNA421367",
        },
        "route_b_ld_reconstructed_or_obtained": False,
        "quantitative_track_generated": False,
        "final_quantitative_figure": None,
        "archived_schematic_figures": [rel(p) for p in [ARCHIVE_PDF, ARCHIVE_PNG] if p.exists()],
        "stolle_vcf_checked_not_downloaded": {"url": STOLLE_VCF_URL, "content_length_bytes": STOLLE_VCF_SIZE_BYTES},
        "yan_inversion_coordinates": YAN_INVERSIONS,
        "published_inversion_union": {
            "chromosome": "chr16",
            "start_bp": min(r["start_bp"] for r in YAN_INVERSIONS),
            "end_bp": max(r["end_bp"] for r in YAN_INVERSIONS),
        },
        "frozen_supergene_interval": {"chromosome": "chr16", "coordinate_start": SUPERGENE_EXPECTED[0], "coordinate_end": SUPERGENE_EXPECTED[1]},
        "assembly_compatibility_for_quantitative_overlay": False,
        "source_files": {rel(p): sha256(p) for p in source_files},
        "frozen_manifest_checksums_before": {rel(k): v for k, v in before.items()},
        "frozen_manifest_checksums_after": {rel(k): v for k, v in after.items()},
        "frozen_manifest_checksums_unchanged": before == after,
        "script_hash": {rel(Path(__file__).resolve()): sha256(Path(__file__).resolve())},
        "output_hashes": {rel(p): sha256(p) for p in existing_outputs},
        "prohibited_actions": {
            "astral_rerun": False,
            "castles_rerun": False,
            "gene_tree_regeneration": False,
            "topology_reclassification": False,
            "new_species_tree_estimation": False,
            "recombination_inferred_from_genealogy": False,
            "heatmap_digitized": False,
            "large_raw_sequence_download": False,
            "stolle_vcf_downloaded": False,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def update_readme() -> None:
    section = textwrap.dedent(f"""
        ## Recombination-suppression validation

        A supplementary audit attempted to replace the earlier schematic recombination bar with a quantitative chromosome-16 recombination or LD track on a real physical coordinate axis. The audit followed two routes: Wang et al. 2013 direct linkage-map markers and Yan et al. 2020 physical-coordinate LD evidence.

        Result: `{STATUS}`. Wang et al. marker-level linkage data were recovered and preserved in `../../data/fire_ants_chr16/processed/recombination_map.tsv`, but no documented `Si_gnF` scaffold-to-Stolle `Si_gnGA` chr16 placement table was recovered. Yan et al. provide exact inversion breakpoints and published chr16 LD figures, but the downloadable supplementary tables do not contain the numeric LD matrix or one-dimensional LD track, and the available raw/genome-wide genotype resources are too large or not source-specific enough to download blindly for this figure.

        The previous schematic figure has been archived as `figures/fire_ants_recombination_genealogy_schematic_archived.pdf` / `.png` for provenance only. It should not be cited as a quantitative recombination-rate or LD track. See `results/recombination_report.md`, `results/recombination_source_audit.md`, `results/recombination_coordinate_audit.md`, and `results/recombination_manifest.json`.

        This validation audit did not rerun or alter the frozen topology, TWISST, ASTRAL, CASTLES-II, or branch-length results.
    """).strip()
    text = README.read_text()
    if "## Recombination-suppression validation" in text:
        start = text.index("## Recombination-suppression validation")
        m = re.search(r"\n## ", text[start + 1:])
        if m:
            end = start + 1 + m.start()
            text = text[:start] + section + "\n" + text[end:]
        else:
            text = text[:start] + section + "\n"
    else:
        text = text.rstrip() + "\n\n" + section + "\n"
    README.write_text(text)


def frozen_manifest_hashes() -> dict[Path, str]:
    return {p: sha256(p) for p in FROZEN_MANIFESTS if p.exists()}


def run_analysis(update_readme_flag: bool = True) -> dict[str, object]:
    before = frozen_manifest_hashes()
    frozen_supergene_interval()
    markers = parse_marker_tables()
    write_recombination_map(markers)
    archive_existing_schematic()
    write_recombination_summary(markers)
    write_source_audit(markers)
    write_coordinate_audit()
    write_report(markers)
    write_caption()
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
        cls.manifest = json.loads(MANIFEST.read_text())

    def test_frozen_manifests_unchanged(self):
        self.assertTrue(self.manifest["frozen_manifest_checksums_unchanged"])
        self.assertEqual(self.manifest["frozen_manifest_checksums_before"], self.manifest["frozen_manifest_checksums_after"])

    def test_frozen_supergene_coordinates_unchanged(self):
        self.assertEqual(frozen_supergene_interval(), SUPERGENE_EXPECTED)

    def test_yan_inversion_coordinates_stored_exactly(self):
        self.assertEqual(self.manifest["yan_inversion_coordinates"], YAN_INVERSIONS)

    def test_no_wang_mapping_accepted_without_documentation(self):
        self.assertFalse(self.manifest["route_a_wang_scaffold_to_chr16_mapping_found"])
        self.assertEqual(self.manifest["route_a_wang_markers_mapped_to_chr16"], 0)
        self.assertFalse(WANG_TO_STOLLE_MAP.exists())
        self.assertFalse(RECOMB_MAP_CHR16.exists())

    def test_orientation_conversion_synthetic(self):
        plus = ScaffoldPlacement("s1", 1, 100, "chr16", 1000, 1099, "+")
        minus = ScaffoldPlacement("s1", 1, 100, "chr16", 1000, 1099, "-")
        self.assertEqual(convert_scaffold_position(25, plus), 1024)
        self.assertEqual(convert_scaffold_position(25, minus), 1075)

    def test_cm_per_mb_synthetic(self):
        self.assertEqual(recombination_cm_per_mb(1_000_000, 2_000_000, Decimal("3.5"), Decimal("5.5")), Decimal("2"))

    def test_zero_physical_distance_rejected(self):
        with self.assertRaises(ValueError):
            recombination_cm_per_mb(10, 10, Decimal("0"), Decimal("1"))

    def test_each_mapping_family_available_separately(self):
        families = {m.family for m in self.markers}
        self.assertEqual(families, {"M013", "M047", "M173", "P008", "P016", "P033", "P034"})

    def test_ld_r2_synthetic_haploid(self):
        self.assertAlmostEqual(calculate_ld_r2([0, 0, 1, 1], [0, 0, 1, 1]), 1.0)
        self.assertAlmostEqual(calculate_ld_r2([0, 0, 1, 1], [0, 1, 0, 1]), 0.0)

    def test_ld_never_labeled_cm_per_mb(self):
        if YAN_LD.exists():
            header = YAN_LD.read_text().splitlines()[0].lower()
            self.assertNotIn("cm_per_mb", header)
        self.assertFalse(self.manifest["route_b_ld_reconstructed_or_obtained"])

    def test_assembly_compatibility_required_before_overlay(self):
        self.assertFalse(self.manifest["assembly_compatibility_for_quantitative_overlay"])
        self.assertFalse(self.manifest["quantitative_track_generated"])

    def test_no_constant_schematic_bar_final_result(self):
        self.assertEqual(self.manifest["status"], STATUS)
        self.assertIsNone(self.manifest["final_quantitative_figure"])
        self.assertGreaterEqual(len(self.manifest["archived_schematic_figures"]), 0)

    def test_marker_tables_parsed(self):
        self.assertEqual(len(self.markers), self.manifest["route_a_wang_markers_parsed"])
        self.assertGreater(len(self.markers), 0)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true")
    parser.add_argument("--no-readme", action="store_true")
    args = parser.parse_args(argv)
    result = run_analysis(update_readme_flag=not args.no_readme)
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(RecombinationValidationTests)
        runner = unittest.TextTestRunner(verbosity=2)
        test_result = runner.run(suite)
        if not test_result.wasSuccessful():
            return 1
        print(f"tests_passed={test_result.testsRun}")
    print(f"status={STATUS}")
    print(f"parsed_linkage_markers={result['markers']}")
    print(f"manifest={rel(MANIFEST)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
