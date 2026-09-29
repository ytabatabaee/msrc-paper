#!/usr/bin/env python3
"""Build the Atlantic cod Stage-1 window-tree table from published source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import tempfile
import unittest
from dataclasses import dataclass
from datetime import date
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"

FIG4_SOURCE = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM6_ESM_source_data_fig4.txt"
FIG2_SOURCE = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM4_ESM_source_data_fig2.txt"
ZENODO_METADATA = DATA_ROOT / "raw" / "zenodo" / "zenodo_4560275_record.json"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"
WINDOW_TABLE = DATA_ROOT / "processed" / "cod_window_trees.tsv"
POP_MAPPING = DATA_ROOT / "metadata" / "population_name_mapping.tsv"
BASELINE_SOURCE = DATA_ROOT / "metadata" / "baseline_tree_source.tsv"
POSTERIOR_INVENTORY = EMPIRICAL_ROOT / "results" / "stage1_posterior_inventory.tsv"
WINDOW_INVENTORY_TSV = EMPIRICAL_ROOT / "results" / "stage1_window_inventory.tsv"
WINDOW_INVENTORY_MD = EMPIRICAL_ROOT / "results" / "stage1_window_inventory.md"
BORNHOLM_MD = EMPIRICAL_ROOT / "results" / "stage1_bornholm_window.md"
PROVENANCE_MD = EMPIRICAL_ROOT / "results" / "stage1_data_provenance.md"
REPORT_MD = EMPIRICAL_ROOT / "results" / "stage1_report.md"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"

FOCAL_LGS = ("LG01", "LG02", "LG07", "LG12")
WINDOW_SIZE = 250_000
MINIMUM_SITES = 500


@dataclass(frozen=True)
class Region:
    lg: str
    region_name: str
    start: int
    end: int


@dataclass
class ParsedTree:
    lg: str
    window_id: str
    start: int
    end: int
    tree_newick: str
    taxa: tuple[str, ...]


def require_tree_parser():
    try:
        from Bio import Phylo  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "Biopython is required for Stage-1 Newick validation. "
            "Install biopython in the configured Python environment before rerunning."
        ) from exc
    return Phylo


def parse_newick_taxa(newick: str) -> tuple[str, ...]:
    phylo = require_tree_parser()
    try:
        tree = phylo.read(io.StringIO(newick), "newick")
    except Exception as exc:  # pragma: no cover - exact parser exception varies
        raise ValueError(f"Malformed Newick tree: {exc}") from exc
    taxa = [terminal.name for terminal in tree.get_terminals()]
    if any(label is None or label == "" for label in taxa):
        raise ValueError("Tree contains an empty terminal label.")
    if len(taxa) != len(set(taxa)):
        raise ValueError("Tree contains duplicate taxon labels.")
    return tuple(sorted(taxa))


def parse_window_id(window_id: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"(LG\d{2})_(\d{9})_(\d{9})", window_id)
    if not match:
        raise ValueError(f"Malformed window ID: {window_id}")
    lg, start, end = match.groups()
    start_int = int(start)
    end_int = int(end)
    if end_int < start_int:
        raise ValueError(f"Window end precedes start: {window_id}")
    return lg, start_int, end_int


def read_regions(path: Path = REGION_MANIFEST) -> dict[str, Region]:
    regions: dict[str, Region] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["lg"] in FOCAL_LGS:
                regions[row["lg"]] = Region(
                    lg=row["lg"],
                    region_name=row["region_name"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                )
    missing = set(FOCAL_LGS) - set(regions)
    if missing:
        raise ValueError(f"Missing focal regions: {', '.join(sorted(missing))}")
    return regions


def classify_window(start: int, end: int, region: Region) -> tuple[bool, bool, int, int]:
    inside = start >= region.start and end <= region.end
    overlaps = start <= region.end and end >= region.start
    return inside, overlaps, start - region.start, region.end - end


def iter_fig4_trees(path: Path = FIG4_SOURCE) -> list[ParsedTree]:
    parsed: list[ParsedTree] = []
    current_lg: str | None = None
    seen: set[str] = set()
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            heading = re.fullmatch(r"### Population-trees \((LG\d{2})\) \([a-d]\)", line)
            if heading:
                current_lg = heading.group(1)
                continue
            if line == "window\ttree":
                continue
            if current_lg is None:
                continue
            if "\t" not in line:
                continue
            window_id, tree_newick = line.split("\t", 1)
            lg, start, end = parse_window_id(window_id)
            if lg != current_lg:
                raise ValueError(f"Window {window_id} appears under {current_lg}.")
            if window_id in seen:
                raise ValueError(f"Duplicate window ID: {window_id}")
            taxa = parse_newick_taxa(tree_newick)
            parsed.append(ParsedTree(lg=lg, window_id=window_id, start=start, end=end, tree_newick=tree_newick, taxa=taxa))
            seen.add(window_id)
    return [tree for tree in parsed if tree.lg in FOCAL_LGS]


def expected_window_ids(lg: str, observed: list[ParsedTree]) -> list[str]:
    if not observed:
        return []
    max_end = max(tree.end for tree in observed)
    ids = []
    start = 1
    while start <= max_end:
        end = start + WINDOW_SIZE - 1
        ids.append(f"{lg}_{start:09d}_{end:09d}")
        start = end + 1
    return ids


def write_window_table(trees: list[ParsedTree], regions: dict[str, Region]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    fieldnames = [
        "lg",
        "window_id",
        "start",
        "end",
        "midpoint",
        "window_size_bp",
        "n_sites",
        "tree_newick",
        "n_taxa",
        "inside_inversion",
        "overlaps_inversion",
        "inversion_id",
        "distance_to_left_boundary",
        "distance_to_right_boundary",
        "source_file",
        "tree_type",
        "posterior_support_available",
        "posterior_tree_file",
        "mcc_tree_file",
        "analysis_reference",
        "notes",
    ]
    for tree in sorted(trees, key=lambda item: (item.lg, item.start, item.end)):
        region = regions[tree.lg]
        inside, overlaps, left_dist, right_dist = classify_window(tree.start, tree.end, region)
        rows.append(
            {
                "lg": tree.lg,
                "window_id": tree.window_id,
                "start": str(tree.start),
                "end": str(tree.end),
                "midpoint": str((tree.start + tree.end) // 2),
                "window_size_bp": str(tree.end - tree.start + 1),
                "n_sites": "",
                "tree_newick": tree.tree_newick,
                "n_taxa": str(len(tree.taxa)),
                "inside_inversion": str(inside).lower(),
                "overlaps_inversion": str(overlaps).lower(),
                "inversion_id": region.region_name,
                "distance_to_left_boundary": str(left_dist),
                "distance_to_right_boundary": str(right_dist),
                "source_file": str(FIG4_SOURCE.relative_to(REPO_ROOT)),
                "tree_type": "published_source_data_mcc_window_tree",
                "posterior_support_available": "not_in_source_data_fig4",
                "posterior_tree_file": "",
                "mcc_tree_file": "",
                "analysis_reference": "Matschiner et al. 2022 Source Data Fig. 4",
                "notes": "n_sites is unavailable in Source Data Fig. 4; windows reflect published tree rows that passed the authors' filtering.",
            }
        )
    WINDOW_TABLE.parent.mkdir(parents=True, exist_ok=True)
    with WINDOW_TABLE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_population_mapping(trees: list[ParsedTree]) -> list[str]:
    labels = sorted({taxon for tree in trees for taxon in tree.taxa})
    with POP_MAPPING.open("w", newline="") as handle:
        fieldnames = ["source_label", "canonical_label", "description", "source", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for label in labels:
            writer.writerow(
                {
                    "source_label": label,
                    "canonical_label": label,
                    "description": "unresolved",
                    "source": "Matschiner et al. 2022 Source Data Fig. 4",
                    "notes": "Canonical label preserves source label; biological population identity is not inferred in Stage 1.",
                }
            )
    return labels


def summarize_inventory(trees: list[ParsedTree], regions: dict[str, Region]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    by_lg = {lg: [tree for tree in trees if tree.lg == lg] for lg in FOCAL_LGS}
    for lg, lg_trees in by_lg.items():
        observed_ids = {tree.window_id for tree in lg_trees}
        expected_ids = expected_window_ids(lg, lg_trees)
        inside = outside = boundary = 0
        for tree in lg_trees:
            full, overlaps, _, _ = classify_window(tree.start, tree.end, regions[lg])
            if full:
                inside += 1
            elif overlaps:
                boundary += 1
            else:
                outside += 1
        first = min(lg_trees, key=lambda tree: tree.start).window_id if lg_trees else ""
        last = max(lg_trees, key=lambda tree: tree.start).window_id if lg_trees else ""
        rows.append(
            {
                "lg": lg,
                "windows_present": str(len(lg_trees)),
                "first_window": first,
                "last_window": last,
                "fully_inside_inversion": str(inside),
                "outside_inversion": str(outside),
                "boundary_overlapping": str(boundary),
                "expected_nonoverlapping_windows_1_to_last_observed": str(len(expected_ids)),
                "missing_or_excluded_windows": str(len(set(expected_ids) - observed_ids)),
            }
        )
    with WINDOW_INVENTORY_TSV.open("w", newline="") as handle:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# Atlantic cod Stage-1 window inventory",
        "",
        "The inventory is based on Matschiner et al. (2022) Source Data Fig. 4. Coordinates are parsed from author window IDs and stored as 1-based inclusive intervals.",
        "",
        "| LG | windows present | first window | last window | fully inside | outside | boundary-overlapping | missing/excluded |",
        "|---|---:|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['lg']} | {row['windows_present']} | {row['first_window']} | {row['last_window']} | "
            f"{row['fully_inside_inversion']} | {row['outside_inversion']} | {row['boundary_overlapping']} | {row['missing_or_excluded_windows']} |"
        )
    WINDOW_INVENTORY_MD.write_text("\n".join(lines) + "\n")
    return rows


def first_newick_from_source(path: Path) -> str:
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("(") and line.endswith(";"):
                return line
    return ""


def write_baseline_source() -> bool:
    tree = first_newick_from_source(FIG2_SOURCE)
    located = bool(tree)
    if located:
        parse_newick_taxa(tree)
    with BASELINE_SOURCE.open("w", newline="") as handle:
        fieldnames = ["source_name", "source_type", "path", "tree_present", "tree_type", "used_for_stage3", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "source_name": "Nature Source Data Fig. 2",
                "source_type": "published_source_data",
                "path": str(FIG2_SOURCE.relative_to(REPO_ROOT)),
                "tree_present": str(located).lower(),
                "tree_type": "maximum_credibility_population_tree_outside_supergenes",
                "used_for_stage3": "false",
                "notes": "Located for source discovery only; Stage 3 baseline classification has not been performed.",
            }
        )
    return located


def load_zenodo_files() -> dict[str, dict[str, str]]:
    if not ZENODO_METADATA.exists():
        return {}
    record = json.loads(ZENODO_METADATA.read_text())
    files: dict[str, dict[str, str]] = {}
    for entry in record.get("files", []):
        files[entry["key"]] = {
            "size": str(entry.get("size", "")),
            "checksum": entry.get("checksum", ""),
            "url": entry.get("links", {}).get("self", ""),
        }
    return files


def write_posterior_inventory(trees: list[ParsedTree]) -> None:
    files = load_zenodo_files()
    with POSTERIOR_INVENTORY.open("w", newline="") as handle:
        fieldnames = [
            "lg",
            "window_id",
            "mcc_available",
            "posterior_available",
            "posterior_tree_count_if_known",
            "source",
            "estimated_size",
            "downloaded",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for tree in sorted(trees, key=lambda item: (item.lg, item.start)):
            writer.writerow(
                {
                    "lg": tree.lg,
                    "window_id": tree.window_id,
                    "mcc_available": "true",
                    "posterior_available": "not_publicly_listed_per_window_in_zenodo_record",
                    "posterior_tree_count_if_known": "",
                    "source": "Source Data Fig. 4 contains MCC-like window trees; Zenodo record lists supergene-wide SNAPP .trees files but no per-window posterior files.",
                    "estimated_size": "",
                    "downloaded": "false",
                }
            )
    # Record supergene-wide posterior files as comments in the provenance/report, not as window rows.
    _ = files


def write_bornholm(trees: list[ParsedTree]) -> dict[str, str]:
    target = [tree for tree in trees if tree.lg == "LG12" and tree.start <= 7_500_000 and tree.end >= 7_750_000]
    if not target:
        target = [tree for tree in trees if tree.lg == "LG12" and tree.start == 7_500_001 and tree.end == 7_750_000]
    row: dict[str, str]
    if target:
        tree = target[0]
        row = {
            "LG": tree.lg,
            "window coordinates": f"{tree.start}-{tree.end}",
            "window ID": tree.window_id,
            "source file": str(FIG4_SOURCE.relative_to(REPO_ROOT)),
            "tree present?": "true",
            "taxa present?": ",".join(tree.taxa),
        }
    else:
        row = {
            "LG": "LG12",
            "window coordinates": "7,500,001-7,750,000 expected",
            "window ID": "",
            "source file": str(FIG4_SOURCE.relative_to(REPO_ROOT)),
            "tree present?": "false",
            "taxa present?": "",
        }
    lines = ["# Atlantic cod Stage-1 Bornholm LG12 window check", ""]
    for key, value in row.items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("No topology interpretation or biological classification was performed in Stage 1.")
    BORNHOLM_MD.write_text("\n".join(lines) + "\n")
    return row


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_download_manifest() -> list[dict[str, str]]:
    files = [
        (
            "GitHub mmatschiner/supergenes cod_phylogenomics source",
            "split_vcf.sh",
            "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/src/split_vcf.sh",
            DATA_ROOT / "raw" / "github_supergenes" / "split_vcf.sh",
            "Cod phylogenomics inversion-coordinate provenance.",
        ),
        (
            "GitHub mmatschiner/supergenes cod_phylogenomics source",
            "make_snapp_xmls_windows.sh",
            "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/src/make_snapp_xmls_windows.sh",
            DATA_ROOT / "raw" / "github_supergenes" / "make_snapp_xmls_windows.sh",
            "Window-design provenance.",
        ),
        (
            "GitHub mmatschiner/supergenes cod_phylogenomics source",
            "make_snapp_xmls_windows.slurm",
            "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/src/make_snapp_xmls_windows.slurm",
            DATA_ROOT / "raw" / "github_supergenes" / "make_snapp_xmls_windows.slurm",
            "Window naming, coordinate, sample, and filtering provenance.",
        ),
        (
            "GitHub mmatschiner/supergenes cod_phylogenomics source",
            "combine_snapp_results_windows.sh",
            "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/src/combine_snapp_results_windows.sh",
            DATA_ROOT / "raw" / "github_supergenes" / "combine_snapp_results_windows.sh",
            "MCC and posterior tree-file provenance.",
        ),
        (
            "Nature Source Data Fig. 4",
            "41559_2022_1661_MOESM6_ESM.txt",
            "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM6_ESM.txt",
            FIG4_SOURCE,
            "Canonical 250-kb published population-window trees.",
        ),
        (
            "Nature Source Data Fig. 2",
            "41559_2022_1661_MOESM4_ESM.txt",
            "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM4_ESM.txt",
            FIG2_SOURCE,
            "Baseline population-tree source discovery only.",
        ),
        (
            "Zenodo record metadata",
            "zenodo_4560275_record.json",
            "https://zenodo.org/api/records/4560275",
            ZENODO_METADATA,
            "Inventory of public MCC/posterior files without archive download.",
        ),
    ]
    rows = []
    for source, original, url, path, purpose in files:
        if not path.exists():
            continue
        rows.append(
            {
                "source": source,
                "original_filename": original,
                "url_or_doi": url,
                "download_date": date.today().isoformat(),
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256(path),
                "purpose": purpose,
                "local_path": str(path.relative_to(REPO_ROOT)),
            }
        )
    with DOWNLOAD_MANIFEST.open("w", newline="") as handle:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_provenance(downloads: list[dict[str, str]]) -> None:
    lines = [
        "# Atlantic cod Stage-1 data provenance",
        "",
        "Stage 1 uses already-computed published outputs only. No raw reads, whole-genome VCFs, BEAST, or SNAPP reruns were used.",
        "",
        "## Upstream workflow observations",
        "",
        "- `split_vcf.sh` freezes the cod phylogenomics inversion intervals as `LG01:9114741-26192386`, `LG02:18489307-24048607`, `LG07:13610591-23019113`, and `LG12:638100-14327837`.",
        "- Other source-repository sub-workflows, especially demography, may use different inversion limits; Stage 1 does not harmonize them.",
        "- The paper also notes assembly-placement uncertainty around part of the LG02 region when comparing gadMor2 and gadMor3. Stage 1 does not resolve this; it uses the `cod_phylogenomics` coordinate set because that is the coordinate system used for the published window-tree workflow normalized here.",
        "- `make_snapp_xmls_windows.sh` sets `window_size=250000`, `min_n_sites=500`, `max_n_sites=1000`, and `min_site_dist=50`.",
        "- `make_snapp_xmls_windows.slurm` encodes window IDs as `LG##_000000001_000250000`, using 1-based inclusive coordinates with zero-padded starts and ends.",
        "- Windows are non-overlapping: each next start is the previous end plus one.",
        "- Incomplete terminal behavior follows the authors' loop over starts while `window_start < last_pos`; Source Data Fig. 4 contains only windows that produced published tree rows.",
        "- For `LG12`, the XML-generation script excludes `Gadmor_lfc1` and `Gadmor_lfc2` before population-tree inference.",
        "- `combine_snapp_results_windows.sh` combines replicate `.trees` files and writes maximum-clade-credibility `.tre` files with TreeAnnotator.",
        "- Source Data Fig. 4 supplies the published window tree Newick strings directly; per-window `n_sites` is not included there.",
        "",
        "## Downloaded files",
        "",
        "| local path | source | size bytes | SHA256 | purpose |",
        "|---|---|---:|---|---|",
    ]
    for row in downloads:
        lines.append(f"| `{row['local_path']}` | {row['source']} | {row['size_bytes']} | `{row['sha256']}` | {row['purpose']} |")
    PROVENANCE_MD.write_text("\n".join(lines) + "\n")


def write_report(
    downloads: list[dict[str, str]],
    inventory: list[dict[str, str]],
    taxa: list[str],
    baseline_located: bool,
    bornholm: dict[str, str],
) -> None:
    total_size = sum(int(row["size_bytes"]) for row in downloads)
    trees_parsed = sum(int(row["windows_present"]) for row in inventory)
    zenodo_files = load_zenodo_files()
    posterior_available = all(
        f"gadus_morhua_supergene_lg{lg[-2:]}_snapp.trees" in zenodo_files for lg in FOCAL_LGS
    )
    lines = [
        "# Atlantic cod Stage-1 report",
        "",
        "Stage 1 converted published processed population-tree source data into a canonical 250-kb window table. No topology enrichment, quartet classification, arrangement-state assignment, SNAPP rerun, or raw-data processing was performed.",
        "",
        "## Data acquired",
        "",
        f"Total downloaded disk size recorded in `download_manifest.tsv`: {total_size} bytes.",
        "",
        "| file | size bytes | purpose |",
        "|---|---:|---|",
    ]
    for row in downloads:
        lines.append(f"| `{row['local_path']}` | {row['size_bytes']} | {row['purpose']} |")
    lines.extend(["", "## Provenance", "", "Primary tree source: Nature Source Data Fig. 4. Baseline-tree discovery source: Nature Source Data Fig. 2. Zenodo was inspected through record metadata only; no archive or large genotype file was downloaded.", ""])
    lines.extend(["## Window inventory", "", "| LG | windows | fully inside | outside | boundary-overlapping | missing/excluded |", "|---|---:|---:|---:|---:|---:|"])
    for row in inventory:
        lines.append(f"| {row['lg']} | {row['windows_present']} | {row['fully_inside_inversion']} | {row['outside_inversion']} | {row['boundary_overlapping']} | {row['missing_or_excluded_windows']} |")
    lines.extend(
        [
            "",
            "## Inversion membership",
            "",
            "Canonical coordinates are 1-based inclusive. `inside_inversion=true` only when the full window lies within the frozen inversion interval. Boundary-overlapping windows are tracked separately with `overlaps_inversion=true` and `inside_inversion=false`.",
            "",
            "## Tree parsing",
            "",
            f"Parsed successfully: {trees_parsed} trees. Failures: 0.",
            "",
            "## Populations",
            "",
            "Canonical labels preserve source labels in Stage 1:",
            "",
            ", ".join(taxa),
            "",
            "## Baseline tree",
            "",
            f"Located: {str(baseline_located).lower()}. The tree is recorded for Stage-3 source discovery only.",
            "",
            "## Posterior trees",
            "",
            f"Supergene-wide SNAPP posterior `.trees` files are listed in Zenodo: {str(posterior_available).lower()}. Per-window posterior files are not publicly listed in the Zenodo record inspected during Stage 1.",
            "",
            "## Bornholm check",
            "",
            f"LG12 7.50-7.75 Mb window located: {bornholm.get('tree present?', 'false')}. Window ID: {bornholm.get('window ID', '')}.",
            "",
            "## Blockers",
            "",
            "- Source Data Fig. 4 does not include per-window `n_sites`, so that field is blank in the canonical table.",
            "- Biological population descriptions for abbreviated source labels remain unresolved in Stage 1.",
            "- Per-window posterior SNAPP samples were not identified in the public Zenodo record; only published source-data trees and supergene-wide posterior files were inventoried.",
            "- Coordinate caution: Stage 1 freezes the `cod_phylogenomics/src/split_vcf.sh` intervals for window-tree analysis. Other source-repository workflows and paper notes, especially around LG02 assembly placement, must not be mixed into this coordinate set silently.",
        ]
    )
    REPORT_MD.write_text("\n".join(lines) + "\n")


def build_all() -> None:
    regions = read_regions()
    trees = iter_fig4_trees()
    if not trees:
        raise ValueError("No trees were parsed from Source Data Fig. 4.")
    downloads = write_download_manifest()
    write_window_table(trees, regions)
    taxa = write_population_mapping(trees)
    inventory = summarize_inventory(trees, regions)
    baseline_located = write_baseline_source()
    write_posterior_inventory(trees)
    bornholm = write_bornholm(trees)
    write_provenance(downloads)
    write_report(downloads, inventory, taxa, baseline_located, bornholm)


class BuildWindowTableTests(unittest.TestCase):
    def test_newick_parsing(self) -> None:
        self.assertEqual(parse_newick_taxa("(A:1,(B:1,C:1):1);"), ("A", "B", "C"))

    def test_window_coordinate_parsing(self) -> None:
        self.assertEqual(parse_window_id("LG12_007500001_007750000"), ("LG12", 7_500_001, 7_750_000))

    def test_full_containment_classification(self) -> None:
        region = Region("LG01", "test", 10, 20)
        self.assertEqual(classify_window(10, 20, region)[:2], (True, True))

    def test_partial_overlap_classification(self) -> None:
        region = Region("LG01", "test", 10, 20)
        self.assertEqual(classify_window(5, 10, region)[:2], (False, True))

    def test_boundary_coordinates(self) -> None:
        region = Region("LG01", "test", 10, 20)
        self.assertEqual(classify_window(1, 9, region)[:2], (False, False))
        self.assertEqual(classify_window(21, 30, region)[:2], (False, False))

    def test_duplicate_window_ids(self) -> None:
        content = "### Population-trees (LG01) (a)\nwindow\ttree\nLG01_000000001_000250000\t(A:1,B:1);\nLG01_000000001_000250000\t(A:1,B:1);\n"
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write(content)
            tmp = Path(handle.name)
        self.addCleanup(lambda: tmp.unlink(missing_ok=True))
        with self.assertRaisesRegex(ValueError, "Duplicate window ID"):
            iter_fig4_trees(tmp)

    def test_duplicate_taxon_labels(self) -> None:
        with self.assertRaisesRegex(ValueError, "duplicate taxon"):
            parse_newick_taxa("(A:1,A:1);")

    def test_malformed_trees(self) -> None:
        with self.assertRaises(ValueError):
            parse_newick_taxa("(A:1,B:1")

    def test_canonical_population_name_mapping(self) -> None:
        tree = ParsedTree("LG01", "LG01_000000001_000250000", 1, 250000, "(Gadmor_a_spc:1,Gadmor_b_spc:1);", ("Gadmor_a_spc", "Gadmor_b_spc"))
        labels = sorted({taxon for taxon in tree.taxa})
        self.assertEqual(labels, ["Gadmor_a_spc", "Gadmor_b_spc"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="Run built-in unit tests and exit.")
    args = parser.parse_args(argv)
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(BuildWindowTableTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        return 0 if result.wasSuccessful() else 1
    build_all()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
