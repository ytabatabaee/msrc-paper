#!/usr/bin/env python3
"""Freeze Atlantic cod arrangement states without reading local tree topology."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import sys
import tempfile
import unittest
from dataclasses import dataclass
from datetime import date
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"

STAGE1_MAPPING = DATA_ROOT / "metadata" / "population_name_mapping.tsv"
POPULATION_ARRANGEMENTS = DATA_ROOT / "processed" / "population_arrangements.tsv"
ARRANGEMENT_ORIENTATION = DATA_ROOT / "metadata" / "arrangement_orientation.tsv"
STRUCTURAL_CANDIDATES = DATA_ROOT / "processed" / "stage2_structural_quartet_candidates.tsv"
STRUCTURAL_CANDIDATES_SHA = DATA_ROOT / "processed" / "stage2_structural_quartet_candidates.sha256"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"
SOURCE_MANIFEST = DATA_ROOT / "metadata" / "source_manifest.tsv"
STAGE2_MANIFEST = EMPIRICAL_ROOT / "results" / "stage2_manifest.json"
STATE_SUMMARY_TSV = EMPIRICAL_ROOT / "results" / "stage2_population_state_summary.tsv"
STATE_SUMMARY_MD = EMPIRICAL_ROOT / "results" / "stage2_population_state_summary.md"
STAGE2_REPORT = EMPIRICAL_ROOT / "results" / "stage2_report.md"

SUPPLEMENT_PDF = DATA_ROOT / "raw" / "nature_supplement" / "41559_2022_1661_MOESM1_ESM_supplementary_information.pdf"
SOURCE_DATA_FIG3 = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM5_ESM_source_data_fig3.txt"
GITHUB_TABLES = DATA_ROOT / "raw" / "github_supergenes" / "tables"

LGS = ("LG01", "LG02", "LG07", "LG12")
STATE_VALUES = {"ancestral", "derived", "unknown", "not_applicable"}
ECOTYPE_VALUES = {"migratory", "stationary", "unspecified", "outgroup"}
FORBIDDEN_GENEALOGY_TOKENS = (
    "tree_newick",
    "local_topology",
    "quartet_support",
    "q1",
    "q2",
    "q3",
    "topology",
    "window_tree",
)


@dataclass(frozen=True)
class Population:
    tree_label: str
    canonical_population: str
    species: str
    location: str
    ecotype: str
    notes: str


POPULATIONS: tuple[Population, ...] = (
    Population("Gadmor_avc_spc", "More_stationary", "Gadus_morhua", "Møre", "stationary", "AVE/Møre coastal-stationary source label."),
    Population("Gadmor_avo_spc", "More_migratory", "Gadus_morhua", "Møre", "migratory", "AVE/Møre oceanic-migratory source label."),
    Population("Gadmor_bat_spc", "Labrador", "Gadus_morhua", "Labrador", "stationary", "Other sampling localities were considered stationary by the authors."),
    Population("Gadmor_bor_spc", "Bornholm_Basin", "Gadus_morhua", "Bornholm Basin", "stationary", "Other sampling localities were considered stationary by the authors."),
    Population("Gadmor_icc_spc", "Iceland_stationary", "Gadus_morhua", "Iceland", "stationary", "ICC/Iceland coastal-stationary source label."),
    Population("Gadmor_ico_spc", "Iceland_migratory", "Gadus_morhua", "Iceland", "migratory", "ICO/Iceland oceanic-migratory source label."),
    Population("Gadmor_kie_spc", "Kiel_Bight", "Gadus_morhua", "Kiel Bight", "stationary", "Other sampling localities were considered stationary by the authors."),
    Population("Gadmor_lfc_spc", "Lofoten_stationary", "Gadus_morhua", "Lofoten", "stationary", "LOF/Lofoten coastal-stationary source label."),
    Population("Gadmor_lfo_spc", "Lofoten_migratory", "Gadus_morhua", "Lofoten", "migratory", "LOF/Lofoten oceanic-migratory source label."),
    Population("Gadmor_low_spc", "Suffolk", "Gadus_morhua", "Suffolk", "stationary", "LOW/Suffolk; other sampling localities were considered stationary by the authors."),
    Population("Gadmor_twc_spc", "Newfoundland_stationary", "Gadus_morhua", "Newfoundland", "stationary", "TWI/Newfoundland coastal-stationary source label."),
    Population("Gadmor_two_spc", "Newfoundland_migratory", "Gadus_morhua", "Newfoundland", "migratory", "TWI/Newfoundland oceanic-migratory source label."),
)

ORIENTATION = {
    "LG01": {
        "gadMor2_state": "derived",
        "gadMor_Stat_state": "ancestral",
        "evidence_type": "outgroup_and_contig_alignment",
        "source": "Matschiner et al. 2022 Table 1 and main-text orientation statement",
        "confidence": "direct_outgroup_alignment",
        "notes": "Haddock outgroup/assembly evidence indicates gadMor2 carries the derived arrangement for LG01.",
    },
    "LG02": {
        "gadMor2_state": "ancestral",
        "gadMor_Stat_state": "derived",
        "evidence_type": "outgroup_colinearity_and_assembly_comparison",
        "source": "Matschiner et al. 2022 Table 1 and main-text orientation statement",
        "confidence": "direct_outgroup_alignment",
        "notes": "A melAeg contig colinear with gadMor2 near an LG02 end indicates the derived arrangement is carried by gadMor_Stat.",
    },
    "LG07": {
        "gadMor2_state": "derived",
        "gadMor_Stat_state": "ancestral",
        "evidence_type": "outgroup_and_contig_alignment",
        "source": "Matschiner et al. 2022 Table 1 and main-text orientation statement",
        "confidence": "direct_outgroup_alignment",
        "notes": "Haddock outgroup/assembly evidence indicates gadMor2 carries the derived arrangement for LG07.",
    },
    "LG12": {
        "gadMor2_state": "ancestral",
        "gadMor_Stat_state": "derived",
        "evidence_type": "demographic_inference_after_contig_mapping_was_uninformative",
        "source": "Matschiner et al. 2022 Table 1, Supplementary Note 2, and Supplementary Table 5 footnote",
        "confidence": "demographic_inference",
        "notes": "Contig mapping did not identify orientation; demographic inference suggested the Lofoten stationary/Suffolk/Kiel Bight haplotype is derived.",
    },
}

DERIVED_LABELS = {
    "LG01": {"Gadmor_two_spc", "Gadmor_ico_spc", "Gadmor_lfo_spc", "Gadmor_avo_spc"},
    "LG02": {"Gadmor_lfc_spc", "Gadmor_avc_spc", "Gadmor_low_spc", "Gadmor_kie_spc"},
    "LG07": {
        "Gadmor_two_spc",
        "Gadmor_twc_spc",
        "Gadmor_bat_spc",
        "Gadmor_ico_spc",
        "Gadmor_icc_spc",
        "Gadmor_lfo_spc",
        "Gadmor_lfc_spc",
    },
    "LG12": {"Gadmor_lfc_spc", "Gadmor_low_spc", "Gadmor_kie_spc"},
}

ASSIGNMENT_SOURCE = {
    "LG01": "Supplementary Table 5 derived-arrangement underlines; Table 1 orientation",
    "LG02": "Supplementary Table 5 derived-arrangement underlines; Table 1 orientation",
    "LG07": "Supplementary Table 5 derived-arrangement underlines; Table 1 orientation",
    "LG12": "Supplementary Table 5 derived-arrangement underlines and footnote; demographic orientation inference",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reject_genealogy_inputs(paths: list[Path]) -> None:
    for path in paths:
        lower = path.name.lower()
        if path.name == "cod_window_trees.tsv":
            raise ValueError("Stage 2 must not read cod_window_trees.tsv.")
        if any(token in lower for token in FORBIDDEN_GENEALOGY_TOKENS):
            raise ValueError(f"Forbidden genealogy/topology input filename: {path}")
        if not path.exists() or not path.is_file():
            continue
        try:
            with path.open(newline="") as handle:
                first = handle.readline().strip()
        except UnicodeDecodeError:
            continue
        columns = first.split("\t")
        for column in columns:
            lowered = column.lower()
            if any(token == lowered or token in lowered for token in FORBIDDEN_GENEALOGY_TOKENS):
                raise ValueError(f"Forbidden genealogy/topology column {column!r} in {path}.")


def read_stage1_labels(path: Path = STAGE1_MAPPING) -> list[str]:
    reject_genealogy_inputs([path])
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "source_label" not in (reader.fieldnames or []):
            raise ValueError("population_name_mapping.tsv lacks source_label.")
        labels = [row["source_label"] for row in reader]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate source labels in population_name_mapping.tsv.")
    expected = sorted(pop.tree_label for pop in POPULATIONS)
    if sorted(labels) != expected:
        raise ValueError(f"Stage-1 labels do not match curated Stage-2 labels: {sorted(labels)}")
    return sorted(labels)


def validate_vocab(state: str, ecotype: str | None = None) -> None:
    if state not in STATE_VALUES:
        raise ValueError(f"Invalid state: {state}")
    if ecotype is not None and ecotype not in ECOTYPE_VALUES:
        raise ValueError(f"Invalid ecotype: {ecotype}")


def state_for(label: str, lg: str) -> str:
    if label in DERIVED_LABELS[lg]:
        return "derived"
    return "ancestral"


def write_population_mapping() -> None:
    with STAGE1_MAPPING.open("w", newline="") as handle:
        fieldnames = ["source_label", "canonical_label", "species", "location", "ecotype", "source", "confidence", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for pop in sorted(POPULATIONS, key=lambda item: item.tree_label):
            validate_vocab("unknown", pop.ecotype)
            writer.writerow(
                {
                    "source_label": pop.tree_label,
                    "canonical_label": pop.canonical_population,
                    "species": pop.species,
                    "location": pop.location,
                    "ecotype": pop.ecotype,
                    "source": "Matschiner et al. 2022 Supplementary Table 4; cod_phylogenomics include_ids source-table convention",
                    "confidence": "explicit_or_source_label_supported",
                    "notes": pop.notes,
                }
            )


def write_orientation() -> None:
    with ARRANGEMENT_ORIENTATION.open("w", newline="") as handle:
        fieldnames = ["lg", "gadMor2_state", "gadMor_Stat_state", "evidence_type", "source", "confidence", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for lg in LGS:
            row = {"lg": lg, **ORIENTATION[lg]}
            writer.writerow(row)


def write_population_arrangements() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pop in sorted(POPULATIONS, key=lambda item: item.tree_label):
        row = {
            "tree_label": pop.tree_label,
            "canonical_population": pop.canonical_population,
            "location": pop.location,
            "ecotype": pop.ecotype,
            "notes": pop.notes,
        }
        for lg in LGS:
            state = state_for(pop.tree_label, lg)
            validate_vocab(state, pop.ecotype)
            row[f"{lg}_state"] = state
            row[f"{lg}_confidence"] = ORIENTATION[lg]["confidence"] if state != "unknown" else "uncertain"
            row[f"{lg}_source"] = ASSIGNMENT_SOURCE[lg]
        rows.append(row)
    fieldnames = [
        "tree_label",
        "canonical_population",
        "location",
        "ecotype",
        "LG01_state",
        "LG01_confidence",
        "LG01_source",
        "LG02_state",
        "LG02_confidence",
        "LG02_source",
        "LG07_state",
        "LG07_confidence",
        "LG07_source",
        "LG12_state",
        "LG12_confidence",
        "LG12_source",
        "notes",
    ]
    with POPULATION_ARRANGEMENTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def is_strict_2_2(states: list[str]) -> bool:
    return states.count("ancestral") == 2 and states.count("derived") == 2 and "unknown" not in states and "not_applicable" not in states


def arrangement_split(taxa: list[str], states: list[str]) -> str:
    if not is_strict_2_2(states):
        raise ValueError("Arrangement split requires exactly 2 ancestral and 2 derived states.")
    ancestral = [taxon for taxon, state in zip(taxa, states) if state == "ancestral"]
    derived = [taxon for taxon, state in zip(taxa, states) if state == "derived"]
    return f"{','.join(ancestral)}|{','.join(derived)}"


def enumerate_quartets(pop_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    by_label = {row["tree_label"]: row for row in pop_rows}
    labels = sorted(by_label)
    for lg in LGS:
        quartet_index = 1
        for combo in itertools.combinations(labels, 4):
            states = [by_label[label][f"{lg}_state"] for label in combo]
            if not is_strict_2_2(states):
                continue
            rows.append(
                {
                    "lg": lg,
                    "quartet_id": f"{lg}_SQ{quartet_index:04d}",
                    "taxon1": combo[0],
                    "taxon2": combo[1],
                    "taxon3": combo[2],
                    "taxon4": combo[3],
                    "taxon1_state": states[0],
                    "taxon2_state": states[1],
                    "taxon3_state": states[2],
                    "taxon4_state": states[3],
                    "arrangement_split": arrangement_split(list(combo), states),
                    "strict_2_2": "true",
                    "state_source_complete": "true",
                    "notes": "Structural metadata only; no local tree topology used.",
                }
            )
            quartet_index += 1
    with STRUCTURAL_CANDIDATES.open("w", newline="") as handle:
        fieldnames = [
            "lg",
            "quartet_id",
            "taxon1",
            "taxon2",
            "taxon3",
            "taxon4",
            "taxon1_state",
            "taxon2_state",
            "taxon3_state",
            "taxon4_state",
            "arrangement_split",
            "strict_2_2",
            "state_source_complete",
            "notes",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    STRUCTURAL_CANDIDATES_SHA.write_text(sha256(STRUCTURAL_CANDIDATES) + "\n")
    return rows


def append_download_manifest() -> None:
    additions = [
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "include_ids.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/include_ids.txt", GITHUB_TABLES / "include_ids.txt", "Stage-2 source-label/ecotype convention evidence."),
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "dsuite_order_colinear.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/dsuite_order_colinear.txt", GITHUB_TABLES / "dsuite_order_colinear.txt", "Stage-2 population-label ordering evidence; no topology data."),
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "dsuite_order_inversion_lg01.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/dsuite_order_inversion_lg01.txt", GITHUB_TABLES / "dsuite_order_inversion_lg01.txt", "Stage-2 population-label ordering evidence; no topology data."),
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "dsuite_order_inversion_lg02.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/dsuite_order_inversion_lg02.txt", GITHUB_TABLES / "dsuite_order_inversion_lg02.txt", "Stage-2 population-label ordering evidence; no topology data."),
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "dsuite_order_inversion_lg07.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/dsuite_order_inversion_lg07.txt", GITHUB_TABLES / "dsuite_order_inversion_lg07.txt", "Stage-2 population-label ordering evidence; no topology data."),
        ("GitHub mmatschiner/supergenes cod_phylogenomics table", "dsuite_order_inversion_lg12.txt", "https://raw.githubusercontent.com/mmatschiner/supergenes/main/cod_phylogenomics/data/tables/dsuite_order_inversion_lg12.txt", GITHUB_TABLES / "dsuite_order_inversion_lg12.txt", "Stage-2 population-label ordering evidence; no topology data."),
        ("Nature Supplementary Information", "41559_2022_1661_MOESM1_ESM.pdf", "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM1_ESM.pdf", SUPPLEMENT_PDF, "Stage-2 Supplementary Tables 4 and 5 plus orientation notes."),
        ("Nature Source Data Fig. 3", "41559_2022_1661_MOESM5_ESM.txt", "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM5_ESM.txt", SOURCE_DATA_FIG3, "Downloaded during Stage-2 discovery but excluded from structural inputs because it contains topology."),
    ]
    existing: list[dict[str, str]] = []
    if DOWNLOAD_MANIFEST.exists():
        with DOWNLOAD_MANIFEST.open(newline="") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    existing_paths = {row["local_path"] for row in existing}
    for source, original, url, path, purpose in additions:
        if not path.exists():
            continue
        local_path = str(path.relative_to(REPO_ROOT))
        if local_path in existing_paths:
            continue
        existing.append(
            {
                "source": source,
                "original_filename": original,
                "url_or_doi": url,
                "download_date": date.today().isoformat(),
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256(path),
                "purpose": purpose,
                "local_path": local_path,
            }
        )
    with DOWNLOAD_MANIFEST.open("w", newline="") as handle:
        fieldnames = ["source", "original_filename", "url_or_doi", "download_date", "size_bytes", "sha256", "purpose", "local_path"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(existing)


def update_source_manifest() -> None:
    if not SOURCE_MANIFEST.exists():
        return
    rows: list[dict[str, str]]
    with SOURCE_MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        fieldnames = handle.readline()
    by_name = {row["source_name"]: row for row in rows}
    if "Nature Source Data Fig. 3" in by_name:
        by_name["Nature Source Data Fig. 3"].update(
            {
                "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM5_ESM.txt",
                "expected_size": str(SOURCE_DATA_FIG3.stat().st_size) + " bytes" if SOURCE_DATA_FIG3.exists() else "327.9 KB",
                "downloaded": "true" if SOURCE_DATA_FIG3.exists() else "false",
                "notes": "Downloaded during Stage-2 discovery but excluded from Stage-2 structural inputs because it contains topology.",
            }
        )
    if "Nature Supplementary Information" not in by_name and SUPPLEMENT_PDF.exists():
        rows.append(
            {
                "source_name": "Nature Supplementary Information",
                "source_type": "supplementary_metadata",
                "citation": "Matschiner et al. 2022 Supplementary Information",
                "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM1_ESM.pdf",
                "data_role": "Stage-2 population labels, ecotypes, orientation notes, and derived-arrangement underlines",
                "expected_format": "PDF",
                "expected_size": str(SUPPLEMENT_PDF.stat().st_size) + " bytes",
                "downloaded": "true",
                "notes": "Used for Supplementary Tables 4 and 5 and Supplementary Note 2; no local window trees used.",
            }
        )
    with SOURCE_MANIFEST.open("w", newline="") as handle:
        fieldnames = ["source_name", "source_type", "citation", "url", "data_role", "expected_format", "expected_size", "downloaded", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summaries(pop_rows: list[dict[str, str]], quartet_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary: list[dict[str, str]] = []
    atlantic = [row for row in pop_rows if row["ecotype"] != "outgroup"]
    for lg in LGS:
        states = [row[f"{lg}_state"] for row in atlantic]
        q_count = sum(1 for row in quartet_rows if row["lg"] == lg)
        summary.append(
            {
                "lg": lg,
                "n_atlantic_cod_populations": str(len(atlantic)),
                "n_ancestral": str(states.count("ancestral")),
                "n_derived": str(states.count("derived")),
                "n_unknown": str(states.count("unknown")),
                "n_strict_2_2_quartets": str(q_count),
                "orientation_evidence": ORIENTATION[lg]["confidence"],
            }
        )
    with STATE_SUMMARY_TSV.open("w", newline="") as handle:
        fieldnames = list(summary[0].keys())
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    lines = [
        "# Atlantic cod Stage-2 population-state summary",
        "",
        "Assignments are structural metadata only. No local window-tree topology was used.",
        "",
        "| LG | ancestral | derived | unknown | strict 2:2 quartets | orientation evidence |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary:
        lines.append(f"| {row['lg']} | {row['n_ancestral']} | {row['n_derived']} | {row['n_unknown']} | {row['n_strict_2_2_quartets']} | {row['orientation_evidence']} |")
    lines.extend(["", "## Assignments", ""])
    for pop in pop_rows:
        lines.append(
            f"- `{pop['tree_label']}` = {pop['canonical_population']} ({pop['ecotype']}): "
            f"LG01 {pop['LG01_state']}, LG02 {pop['LG02_state']}, LG07 {pop['LG07_state']}, LG12 {pop['LG12_state']}."
        )
    STATE_SUMMARY_MD.write_text("\n".join(lines) + "\n")
    return summary


def write_manifest(summary: list[dict[str, str]], quartet_rows: list[dict[str, str]]) -> None:
    input_paths = [
        STAGE1_MAPPING,
        SUPPLEMENT_PDF,
        GITHUB_TABLES / "include_ids.txt",
        GITHUB_TABLES / "dsuite_order_colinear.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg01.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg02.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg07.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg12.txt",
    ]
    manifest = {
        "stage": 2,
        "script": str(Path(__file__).relative_to(REPO_ROOT)),
        "script_sha256": sha256(Path(__file__)),
        "local_tree_topology_data_read": False,
        "forbidden_inputs": ["data/atlantic_cod/processed/cod_window_trees.tsv", "Source Data Fig. 4 tree_newick", "quartet support", "window topology"],
        "structural_metadata_inputs": [
            {"path": str(path.relative_to(REPO_ROOT)), "sha256": sha256(path)}
            for path in input_paths
            if path.exists()
        ],
        "excluded_downloaded_files": [
            {
                "path": str(SOURCE_DATA_FIG3.relative_to(REPO_ROOT)),
                "reason": "Downloaded during discovery but excluded because it contains topology.",
                "sha256": sha256(SOURCE_DATA_FIG3) if SOURCE_DATA_FIG3.exists() else "",
            }
        ],
        "number_of_populations": len(POPULATIONS),
        "known_states_per_lg": {
            row["lg"]: int(row["n_ancestral"]) + int(row["n_derived"])
            for row in summary
        },
        "unknown_states_per_lg": {row["lg"]: int(row["n_unknown"]) for row in summary},
        "strict_2_2_quartet_candidates_per_lg": {
            row["lg"]: int(row["n_strict_2_2_quartets"])
            for row in summary
        },
        "candidate_table": str(STRUCTURAL_CANDIDATES.relative_to(REPO_ROOT)),
        "candidate_table_sha256": sha256(STRUCTURAL_CANDIDATES),
        "anti_circularity_declaration": "No local tree topology, Source Data Fig. 4 tree grouping, quartet support, or window topology information was used to construct the Stage-2 arrangement predictions.",
    }
    STAGE2_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_report(pop_rows: list[dict[str, str]], summary: list[dict[str, str]]) -> None:
    lines = [
        "# Atlantic cod Stage-2 report",
        "",
        "Stage 2 freezes population identities, ecotype labels, and chromosomal arrangement states from structural metadata only.",
        "",
        "## Population-label resolution",
        "",
        "| tree label | canonical population | location | ecotype |",
        "|---|---|---|---|",
    ]
    for pop in pop_rows:
        lines.append(f"| `{pop['tree_label']}` | {pop['canonical_population']} | {pop['location']} | {pop['ecotype']} |")
    lines.extend(["", "All 12 Stage-1 tree labels are Atlantic cod population/ecotype labels; no outgroup labels occur among them.", "", "## Arrangement orientation", "", "| LG | gadMor2 | gadMor_Stat | evidence | confidence |", "|---|---|---|---|---|"])
    for lg in LGS:
        orientation = ORIENTATION[lg]
        lines.append(f"| {lg} | {orientation['gadMor2_state']} | {orientation['gadMor_Stat_state']} | {orientation['evidence_type']} | {orientation['confidence']} |")
    lines.extend(["", "## Population arrangement states", "", "| population | LG01 | LG02 | LG07 | LG12 |", "|---|---|---|---|---|"])
    for pop in pop_rows:
        lines.append(f"| {pop['canonical_population']} | {pop['LG01_state']} | {pop['LG02_state']} | {pop['LG07_state']} | {pop['LG12_state']} |")
    lines.extend(["", "## Strict 2:2 candidates", "", "| LG | candidates |", "|---|---:|"])
    for row in summary:
        lines.append(f"| {row['lg']} | {row['n_strict_2_2_quartets']} |")
    lines.extend(
        [
            "",
            "No topology support, local tree grouping, or inside/outside enrichment was calculated.",
            "",
            "## Ambiguities",
            "",
            "- LG12 orientation is weaker than LG01, LG02, and LG07 because contig mapping was uninformative; the Stage-2 orientation follows the authors' demographic inference and Supplementary Table 5 footnote.",
            "- Source Data Fig. 3 was downloaded during discovery but excluded from Stage-2 structural inputs because it contains topology.",
            "- No arrangement state is inferred from local tree clustering or from the Bornholm LG12 topology switch.",
            "",
            "## Anti-circularity declaration",
            "",
            "No local tree topology, Source Data Fig. 4 tree grouping, quartet support, or window topology information was used to construct the Stage-2 arrangement predictions.",
        ]
    )
    STAGE2_REPORT.write_text("\n".join(lines) + "\n")


def build_all() -> None:
    approved_inputs = [
        STAGE1_MAPPING,
        SUPPLEMENT_PDF,
        GITHUB_TABLES / "include_ids.txt",
        GITHUB_TABLES / "dsuite_order_colinear.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg01.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg02.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg07.txt",
        GITHUB_TABLES / "dsuite_order_inversion_lg12.txt",
    ]
    reject_genealogy_inputs(approved_inputs)
    read_stage1_labels()
    write_population_mapping()
    write_orientation()
    pop_rows = write_population_arrangements()
    quartet_rows = enumerate_quartets(pop_rows)
    append_download_manifest()
    update_source_manifest()
    summary = write_summaries(pop_rows, quartet_rows)
    write_manifest(summary, quartet_rows)
    write_report(pop_rows, summary)


class Stage2Tests(unittest.TestCase):
    def test_vocab_validation(self) -> None:
        validate_vocab("ancestral", "stationary")
        with self.assertRaises(ValueError):
            validate_vocab("0", "stationary")

    def test_unknown_handling(self) -> None:
        self.assertFalse(is_strict_2_2(["ancestral", "derived", "derived", "unknown"]))

    def test_exact_2_2_detection(self) -> None:
        self.assertTrue(is_strict_2_2(["ancestral", "ancestral", "derived", "derived"]))
        self.assertFalse(is_strict_2_2(["ancestral", "derived", "derived", "derived"]))

    def test_arrangement_split(self) -> None:
        self.assertEqual(arrangement_split(["A", "B", "C", "D"], ["ancestral", "derived", "ancestral", "derived"]), "A,C|B,D")

    def test_deterministic_quartet_ordering(self) -> None:
        rows = [
            {"tree_label": "B", "LG01_state": "derived"},
            {"tree_label": "A", "LG01_state": "ancestral"},
            {"tree_label": "D", "LG01_state": "derived"},
            {"tree_label": "C", "LG01_state": "ancestral"},
        ]
        for row in rows:
            for lg in ("LG02", "LG07", "LG12"):
                row[f"{lg}_state"] = "unknown"
        with tempfile.TemporaryDirectory() as tmp:
            old = STRUCTURAL_CANDIDATES
            old_sha = STRUCTURAL_CANDIDATES_SHA
            try:
                globals()["STRUCTURAL_CANDIDATES"] = Path(tmp) / "q.tsv"
                globals()["STRUCTURAL_CANDIDATES_SHA"] = Path(tmp) / "q.sha256"
                out = enumerate_quartets(rows)
            finally:
                globals()["STRUCTURAL_CANDIDATES"] = old
                globals()["STRUCTURAL_CANDIDATES_SHA"] = old_sha
        self.assertEqual(out[0]["taxon1"], "A")
        self.assertEqual(out[0]["arrangement_split"], "A,C|B,D")

    def test_duplicate_population_rejection(self) -> None:
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write("source_label\nA\nA\n")
            path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        with self.assertRaises(ValueError):
            read_stage1_labels(path)

    def test_outgroup_vocab(self) -> None:
        validate_vocab("not_applicable", "outgroup")

    def test_incomplete_state_exclusion(self) -> None:
        self.assertFalse(is_strict_2_2(["ancestral", "ancestral", "derived", "not_applicable"]))

    def test_sha256_reproducibility(self) -> None:
        with tempfile.NamedTemporaryFile("wb", delete=False) as handle:
            handle.write(b"abc")
            path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        self.assertEqual(sha256(path), sha256(path))

    def test_genealogy_column_rejection(self) -> None:
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write("source_label\ttree_newick\nA\t(A,B);\n")
            path = Path(handle.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        with self.assertRaises(ValueError):
            reject_genealogy_inputs([path])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="Run built-in Stage-2 tests and exit.")
    args = parser.parse_args(argv)
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage2Tests))
        return 0 if result.wasSuccessful() else 1
    build_all()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
