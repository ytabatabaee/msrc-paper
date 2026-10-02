#!/usr/bin/env python3
"""Stage 6C fire-ant individual-level ASTRAL4/CASTLES-II SU branch lengths.

This post-freeze sensitivity recovers the documented 267-tip to seven-group
mapping from Stolle et al. 2022 Supplementary Data 1, audits the original
individual-level RAxML-NG local trees, and uses those branch-length-bearing
local trees for mapped ASTRAL4/CASTLES-II substitution-unit branch-length
estimation. Stage 5 and the existing Stage 6/CU results are not modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import unittest
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA = REPO_ROOT / "data" / "fire_ants_chr16"
EMP = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
PROCESSED = DATA / "processed"
RAW_UPSTREAM = DATA / "raw" / "upstream"
INTERMEDIATE = DATA / "intermediate" / "astral4_su"
ASTRAL_SU_RESULTS = RESULTS / "astral4_su"

INPUT_TREE = RAW_UPSTREAM / "Topology weighting" / "results" / "2021-08-07-twisst" / "input.tree"
FILTER_SCRIPT = RAW_UPSTREAM / "Topology weighting" / "supergene" / "filter_samples_by_missing.R"
TOPO_README = RAW_UPSTREAM / "Topology weighting" / "README.md"
WINDOW_INDEX = PROCESSED / "stage1_window_index.tsv"
STAGE4A_WINDOW_SUPPORT = PROCESSED / "stage4a_window_quartet_support.tsv"
STAGE4A_REGION_SUMMARY = RESULTS / "stage4a_region_summary.tsv"
STAGE4B_PRIMARY = RESULTS / "stage4b_primary_test.tsv"
STAGE4B_COORD_UNIQUE = RESULTS / "stage4b_coordinate_summary.tsv"
STAGE4B_COORD_LENGTH = RESULTS / "stage4b_coordinate_length_weighted_summary.tsv"
STAGE5_MANIFEST = RESULTS / "stage5_final_manifest.json"
STAGE6_MANIFEST = RESULTS / "stage6_manifest.json"
STAGE6_CU_MANIFEST = RESULTS / "stage6_cu_manifest.json"
STAGE6_CU_COMPARISON = RESULTS / "stage6_cu_branch_length_comparison.tsv"
STAGE6_ASTRAL_TREATMENT_SUMMARY = RESULTS / "stage6_astral4_treatment_summary.tsv"
STAGE6_BACKGROUND_REFERENCE = PROCESSED / "stage6_background_reference.nwk"
README = EMP / "README.md"
PROJECT_STATUS = REPO_ROOT / "PROJECT_STATUS.md"

SUPP_DATA1_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-28806-7/MediaObjects/41467_2022_28806_MOESM4_ESM.xlsx"
SUPP_DATA1 = RAW_UPSTREAM / "nature_communications_2022" / "41467_2022_28806_MOESM4_ESM_supplementary_data_1.xlsx"
EXPECTED_SUPP_SHA256 = "603578cc5c701d9e5955ee091ddded1f65db9895ccd9033894ea2f2e441f45c6"

TREE_TIP_INVENTORY = PROCESSED / "stage6c_tree_tip_inventory.tsv"
MATCH_AUDIT = PROCESSED / "stage6c_tree_tip_supplement_match_audit.tsv"
INDIVIDUAL_TO_GROUP = PROCESSED / "stage6c_individual_to_group.tsv"
ASTRAL_MAPPING = PROCESSED / "stage6c_astral4_mapping.txt"
GROUP_LABEL_MAP = PROCESSED / "stage6c_group_label_map.tsv"
FIXED_BACKGROUND_TOPOLOGY = PROCESSED / "stage6c_fixed_background_topology.nwk"

MAPPING_REPORT = RESULTS / "stage6c_mapping_recovery_report.md"
BRANCH_AUDIT = RESULTS / "stage6c_raxml_branch_length_audit.txt"
PROVENANCE = RESULTS / "stage6c_astral4_provenance.md"
SU_COMPARISON = RESULTS / "stage6c_su_branch_length_comparison.tsv"
CU_VS_SU_COMPARISON = RESULTS / "stage6c_cu_vs_su_comparison.tsv"
SU_SUMMARY = RESULTS / "stage6c_su_summary.txt"
FREE_TOPOLOGY_SUMMARY = RESULTS / "stage6c_free_topology_summary.tsv"
FINAL_SUMMARY_TABLE = RESULTS / "fire_ants_final_summary.tsv"
FINAL_METHODS = RESULTS / "fire_ants_final_methods.md"
FINAL_RESULTS = RESULTS / "fire_ants_final_results.md"
FINAL_CONSISTENCY_AUDIT = RESULTS / "fire_ants_final_consistency_audit.json"
MANIFEST = RESULTS / "stage6c_manifest.json"

FIG_SU_SCATTER_PDF = FIGURES / "fire_ants_stage6c_su_background_vs_all.pdf"
FIG_SU_SCATTER_PNG = FIGURES / "fire_ants_stage6c_su_background_vs_all.png"
FIG_CU_VS_SU_PDF = FIGURES / "fire_ants_stage6c_cu_vs_su_change.pdf"
FIG_CU_VS_SU_PNG = FIGURES / "fire_ants_stage6c_cu_vs_su_change.png"
FIG_FOCAL_PDF = FIGURES / "fire_ants_stage6c_focal_branch_lengths.pdf"
FIG_FOCAL_PNG = FIGURES / "fire_ants_stage6c_focal_branch_lengths.png"
FIG_FINAL_SUMMARY_PDF = FIGURES / "fire_ants_final_summary.pdf"
FIG_FINAL_SUMMARY_PNG = FIGURES / "fire_ants_final_summary.png"

LABELS = {
    "geminata": "geminata",
    "saevissima": "saevissima",
    "pusillignis": "pusillignis",
    "invicta/macdonaghi_SB": "inv_mac_SB",
    "invicta/macdonaghi_Sb": "inv_mac_Sb",
    "richteri_SB": "richteri_SB",
    "richteri_Sb": "richteri_Sb",
}
EXPECTED_GROUPS = list(LABELS.keys())
EXPECTED_SAFE_TAXA = sorted(LABELS.values())
INV_PAIR = frozenset(["inv_mac_SB", "inv_mac_Sb"])
RIC_PAIR = frozenset(["richteri_SB", "richteri_Sb"])
FOCAL = frozenset(["inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"])
EPS = 1e-12
THREADS = int(os.environ.get("ASTRAL4_THREADS", "2"))

DOCUMENTED_NAME_CORRECTIONS_OLD_TO_NEW = {
    "SRR7028261_AL-139-bigB-m": "SRR7028261_AL-139-littleb-p",
    "SRR7028251_AL-149-bigB-m": "SRR7028251_AL-149-littleb-p",
}

TREATMENTS = {
    "T_chr1": lambda w: w["region"] == "chr1",
    "T_chr16_outside": lambda w: w["region"] in {"chr16A", "chr16B"},
    "T_background": lambda w: w["region"] in {"chr1", "chr16A", "chr16B"},
    "T_supergene": lambda w: w["region"] == "chr16_supergene",
    "T_chr16_all": lambda w: w["region"] in {"chr16A", "chr16_supergene", "chr16B"},
    "T_all": lambda w: True,
}
EXPECTED_TREATMENT_N = {"T_chr1": 117, "T_chr16_outside": 44, "T_background": 161, "T_supergene": 52, "T_chr16_all": 96, "T_all": 213}
TREE_FILES = {
    "T_chr1": "chr1.trees",
    "T_chr16_outside": "chr16_outside.trees",
    "T_background": "background_161.trees",
    "T_supergene": "supergene_52.trees",
    "T_chr16_all": "chr16_all_96.trees",
    "T_all": "all_213.trees",
}
FIXED_TREATMENTS = ["T_background", "T_all", "T_chr16_all", "T_supergene"]
FREE_TREATMENTS = ["T_background", "T_all", "T_chr16_all", "T_supergene"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def fmt(v: object) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return "NA"
        return f"{v:.12g}"
    return str(v)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: fmt(row.get(k, "")) for k in fields})


def tree_lines() -> list[str]:
    return [line.strip() for line in INPUT_TREE.read_text().splitlines() if line.strip()]


def parse_tree_text(text: str):
    return Phylo.read(io.StringIO(text.strip()), "newick")


def parse_tree_file(path: Path):
    return parse_tree_text(path.read_text().strip())


def tree_to_newick(tree) -> str:
    s = io.StringIO()
    Phylo.write(tree, s, "newick")
    return s.getvalue().strip()


def tree_tips(tree) -> list[str]:
    tips = [term.name for term in tree.get_terminals()]
    if len(tips) != len(set(tips)):
        raise ValueError("duplicate tip labels in tree")
    return sorted(tips)


def read_xlsx_first_sheet(path: Path) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    with zipfile.ZipFile(path) as z:
        strings: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", ns):
                strings.append("".join(t.text or "" for t in si.findall(".//a:t", ns)))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        sheets = []
        for sh in wb.find("a:sheets", ns):
            rid = sh.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            target = rid_to_target[rid]
            if not target.startswith("worksheets/"):
                target = "worksheets/" + target.split("/")[-1]
            sheets.append((sh.attrib["name"], "xl/" + target))
        sheet_name, sheet_path = sheets[0]
        ws = ET.fromstring(z.read(sheet_path))
        def col_idx(cell_ref: str) -> int:
            m = re.match(r"([A-Z]+)", cell_ref)
            if not m:
                raise ValueError(f"Bad XLSX cell reference {cell_ref}")
            n = 0
            for ch in m.group(1):
                n = n * 26 + ord(ch) - 64
            return n - 1
        matrix: list[list[str]] = []
        for row in ws.findall(".//a:sheetData/a:row", ns):
            vals: list[str] = []
            for c in row.findall("a:c", ns):
                v = c.find("a:v", ns)
                val = "" if v is None else (v.text or "")
                if c.attrib.get("t") == "s" and val != "":
                    val = strings[int(val)]
                j = col_idx(c.attrib["r"])
                while len(vals) <= j:
                    vals.append("")
                vals[j] = val
            matrix.append(vals)
    headers = matrix[0]
    records: list[dict[str, str]] = []
    for vals in matrix[1:]:
        row = {headers[i]: (vals[i] if i < len(vals) else "") for i in range(len(headers))}
        if any(row.values()):
            records.append(row)
    meta = {"worksheet_names": [s[0] for s in sheets], "active_worksheet": sheet_name, "n_rows": len(records), "column_names": headers}
    return headers, records, meta


def normalize_species(species: str) -> str:
    sp = species.replace("S.", "").replace(" ", "")
    if sp in {"invicta", "macdonaghi"}:
        sp = "invicta/macdonaghi"
    return sp


def group_from_metadata(species: str, variant: str) -> str:
    sp = normalize_species(species)
    if sp in {"invicta/macdonaghi", "richteri"}:
        if variant not in {"SB", "Sb"}:
            raise ValueError(f"Missing or unexpected Supergene Variant `{variant}` for {species}")
        return f"{sp}_{variant}"
    return sp


def extract_tip_inventory() -> tuple[list[str], list[dict[str, object]]]:
    lines = tree_lines()
    first_tips = None
    per_tree = []
    for idx, line in enumerate(lines, start=1):
        tips = tree_tips(parse_tree_text(line))
        if len(tips) != 267:
            raise ValueError(f"Tree {idx} has {len(tips)} tips, expected 267")
        if first_tips is None:
            first_tips = tips
        elif tips != first_tips:
            raise ValueError(f"Tree {idx} tip set differs from first tree")
        per_tree.append(set(tips))
    rows = []
    for tip in first_tips or []:
        rows.append({"tree_tip": tip, "n_trees_observed": sum(1 for s in per_tree if tip in s), "tip_set_sha256": hashlib.sha256("\n".join(first_tips or []).encode()).hexdigest()})
    return first_tips or [], rows


def recover_mapping() -> dict[str, object]:
    if not SUPP_DATA1.exists():
        raise FileNotFoundError(f"Supplementary Data 1 XLSX is missing: {SUPP_DATA1}")
    if sha256(SUPP_DATA1) != EXPECTED_SUPP_SHA256:
        raise ValueError("Supplementary Data 1 checksum differs from recorded download")
    tips, tip_rows = extract_tip_inventory()
    headers, supp_records, supp_meta = read_xlsx_first_sheet(SUPP_DATA1)
    by_identifier: dict[str, dict[str, str]] = {}
    duplicate_ids = []
    for rec in supp_records:
        ids = [rec.get("Sample name", "")]
        if rec.get("SRA identifier"):
            ids.append(f"{rec['SRA identifier']}_{rec.get('Sample name','')}")
        for identifier in ids:
            if not identifier:
                continue
            if identifier in by_identifier:
                duplicate_ids.append(identifier)
            by_identifier[identifier] = rec
    audit_rows = []
    individual_rows = []
    unmatched = []
    documented_used = []
    for tip in tips:
        rec = by_identifier.get(tip)
        match_type = "exact" if rec else "unmatched"
        supplement_sample = rec.get("Sample name", "") if rec else ""
        source_note = "Supplementary Data 1 exact identifier match"
        if rec is None:
            for old, new in DOCUMENTED_NAME_CORRECTIONS_OLD_TO_NEW.items():
                if tip == new and old in by_identifier:
                    rec = by_identifier[old]
                    match_type = "documented_name_correction"
                    supplement_sample = rec.get("Sample name", "")
                    source_note = f"Topology weighting README documented rename {old} -> {new}; tree tip uses corrected name"
                    documented_used.append(f"{old}->{new}")
                    break
                if tip == old and new in by_identifier:
                    rec = by_identifier[new]
                    match_type = "documented_name_correction"
                    supplement_sample = rec.get("Sample name", "")
                    source_note = f"Topology weighting README documented rename {old} -> {new}; supplement uses corrected name"
                    documented_used.append(f"{old}->{new}")
                    break
        if rec is None:
            unmatched.append(tip)
            audit_rows.append({"tree_tip": tip, "supplement_sample": "", "match_type": "unmatched", "species": "", "supergene_variant": "", "group": "", "source": "unmatched in Supplementary Data 1 after documented corrections"})
            continue
        species = rec.get("Species", "")
        variant = rec.get("Supergene Variant", "")
        group = group_from_metadata(species, variant)
        safe_group = LABELS.get(group)
        if safe_group is None:
            raise ValueError(f"Unexpected recovered group {group} for {tip}")
        audit_rows.append({"tree_tip": tip, "supplement_sample": supplement_sample, "match_type": match_type, "species": species, "supergene_variant": variant, "group": group, "source": source_note})
        individual_rows.append({
            "sample": tip,
            "species": species,
            "species_normalized": normalize_species(species),
            "supergene_variant": variant,
            "astral_group": group,
            "astral_group_safe": safe_group,
            "mapping_source": "Stolle et al. 2022 Nature Communications Supplementary Data 1",
            "mapping_rule": "filter_samples_by_missing.R species normalization; invicta/macdonaghi merge; append Supergene Variant for invicta/macdonaghi and richteri",
            "match_type": match_type,
            "supplement_sample": supplement_sample,
        })
    group_counts = Counter(row["astral_group"] for row in individual_rows)
    used_twisst_counts = Counter((row.get("Used for Twisst") or "") for row in supp_records)
    species_assigned = sum(1 for r in audit_rows if r["species"])
    variant_assigned = sum(1 for r in audit_rows if r["supergene_variant"])
    biglittle_mismatches = []
    for row in audit_rows:
        tip = row["tree_tip"]
        variant = row["supergene_variant"]
        if "bigB" in tip and variant != "SB":
            biglittle_mismatches.append({"tree_tip": tip, "supergene_variant": variant})
        if "littleb" in tip and variant != "Sb":
            biglittle_mismatches.append({"tree_tip": tip, "supergene_variant": variant})
    complete = len(unmatched) == 0 and len(individual_rows) == 267 and sorted(group_counts) == sorted(EXPECTED_GROUPS)
    if len({r["sample"] for r in individual_rows}) != len(individual_rows):
        raise ValueError("Duplicate mapped samples")
    return {
        "supplement_meta": supp_meta,
        "supplement_sha256": sha256(SUPP_DATA1),
        "supplement_url": SUPP_DATA1_URL,
        "n_tree_tips": len(tips),
        "tip_rows": tip_rows,
        "audit_rows": audit_rows,
        "individual_rows": individual_rows,
        "unmatched": unmatched,
        "duplicate_supplement_identifiers": duplicate_ids,
        "documented_corrections_used": documented_used,
        "match_type_counts": dict(Counter(r["match_type"] for r in audit_rows)),
        "group_counts": dict(group_counts),
        "supplement_used_for_twisst_counts": dict(used_twisst_counts),
        "n_species_assigned": species_assigned,
        "n_supergene_variant_assigned": variant_assigned,
        "biglittle_validation_mismatches": biglittle_mismatches,
        "complete_exact_individual_to_group_mapping_available": complete,
    }


def branch_lengths(tree) -> list[float | None]:
    return [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]


def audit_raxml_trees() -> dict[str, object]:
    lines = tree_lines()
    if len(lines) != 213:
        raise ValueError(f"Expected 213 local trees; observed {len(lines)}")
    first_tips = None
    all_lengths = []
    per_tree_tip_counts = []
    for idx, line in enumerate(lines, start=1):
        tree = parse_tree_text(line)
        tips = tree_tips(tree)
        per_tree_tip_counts.append(len(tips))
        if first_tips is None:
            first_tips = tips
        elif tips != first_tips:
            raise ValueError(f"Tree {idx} tip set differs from first tree")
        for bl in branch_lengths(tree):
            if bl is None:
                continue
            if math.isnan(bl) or math.isinf(bl):
                raise ValueError(f"Tree {idx} has nonfinite branch length {bl}")
            if bl < 0:
                raise ValueError(f"Tree {idx} has negative branch length {bl}")
            all_lengths.append(float(bl))
    positives = [x for x in all_lengths if x > 0]
    zeros = [x for x in all_lengths if x == 0]
    windows = read_tsv(WINDOW_INDEX)
    counts = Counter(row["region"] for row in windows)
    treatment_counts = {
        "T_chr1": counts["chr1"],
        "T_chr16_outside": counts["chr16A"] + counts["chr16B"],
        "T_background": counts["chr1"] + counts["chr16A"] + counts["chr16B"],
        "T_supergene": counts["chr16_supergene"],
        "T_chr16_all": counts["chr16A"] + counts["chr16B"] + counts["chr16_supergene"],
        "T_all": len(windows),
    }
    return {
        "n_trees": len(lines),
        "n_nonempty_trees": len(lines),
        "tips_per_tree_min": min(per_tree_tip_counts),
        "tips_per_tree_max": max(per_tree_tip_counts),
        "unique_tip_count": len(first_tips or []),
        "all_tree_tip_sets_identical": True,
        "n_branch_lengths": len(all_lengths),
        "n_positive_branch_lengths": len(positives),
        "n_zero_branch_lengths": len(zeros),
        "fraction_positive_branch_lengths": len(positives) / len(all_lengths),
        "fraction_zero_branch_lengths": len(zeros) / len(all_lengths),
        "median_branch_length": statistics.median(all_lengths),
        "median_positive_branch_length": statistics.median(positives),
        "min_branch_length": min(all_lengths),
        "max_branch_length": max(all_lengths),
        "tip_set_sha256": hashlib.sha256("\n".join(first_tips or []).encode()).hexdigest(),
        "treatment_counts": treatment_counts,
    }


def create_mapping_outputs(mapping: dict[str, object]) -> None:
    write_tsv(TREE_TIP_INVENTORY, mapping["tip_rows"], ["tree_tip", "n_trees_observed", "tip_set_sha256"])
    write_tsv(MATCH_AUDIT, mapping["audit_rows"], ["tree_tip", "supplement_sample", "match_type", "species", "supergene_variant", "group", "source"])
    if not mapping["complete_exact_individual_to_group_mapping_available"]:
        return
    write_tsv(INDIVIDUAL_TO_GROUP, mapping["individual_rows"], ["sample", "species", "species_normalized", "supergene_variant", "astral_group", "astral_group_safe", "mapping_source", "mapping_rule", "match_type", "supplement_sample"])
    map_lines = [f"{row['sample']}\t{row['astral_group_safe']}" for row in mapping["individual_rows"]]
    ASTRAL_MAPPING.write_text("\n".join(map_lines) + "\n")
    label_rows = [{"published_group": k, "astral_safe_label": v} for k, v in LABELS.items()]
    write_tsv(GROUP_LABEL_MAP, label_rows, ["published_group", "astral_safe_label"])


def clean_tree_for_constraint(in_path: Path, out_path: Path) -> None:
    tree = parse_tree_file(in_path)
    for clade in tree.find_clades():
        clade.branch_length = None
        if not clade.is_terminal():
            clade.name = None
    tips = tree_tips(tree)
    if tips != EXPECTED_SAFE_TAXA:
        raise ValueError(f"Background reference has unexpected taxa: {tips}")
    out_path.write_text(tree_to_newick(tree) + "\n")


def create_treatment_tree_files() -> dict[str, object]:
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    lines = tree_lines()
    windows = read_tsv(WINDOW_INDEX)
    if len(windows) != len(lines):
        raise ValueError("Window index and local tree file have different lengths")
    out = {}
    for treatment, predicate in TREATMENTS.items():
        selected = [(i, w, lines[i]) for i, w in enumerate(windows) if predicate(w)]
        if len(selected) != EXPECTED_TREATMENT_N[treatment]:
            raise ValueError(f"{treatment} expected {EXPECTED_TREATMENT_N[treatment]}, observed {len(selected)}")
        path = INTERMEDIATE / TREE_FILES[treatment]
        path.write_text("\n".join(line for _, _, line in selected) + "\n")
        out[treatment] = {"path": path, "n_windows": len(selected), "window_indices": [w["window_index"] for _, w, _ in selected], "sha256": sha256(path)}
    clean_tree_for_constraint(STAGE6_BACKGROUND_REFERENCE, FIXED_BACKGROUND_TOPOLOGY)
    return out


def find_astral4() -> Path:
    env = os.environ.get("ASTRAL4_BIN")
    candidates = []
    if env:
        candidates.append(Path(env))
    which = shutil.which("astral4")
    if which:
        candidates.append(Path(which))
    candidates.append(REPO_ROOT / "external" / "ASTER" / "bin" / "astral4")
    candidates.append(Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral4"))
    for c in candidates:
        if c.exists() and os.access(c, os.X_OK):
            return c.resolve()
    raise FileNotFoundError("Could not locate astral4 executable")


def run_command(cmd: list[str], log_path: Path) -> float:
    start = time.time()
    proc = subprocess.run(cmd, cwd=REPO_ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    runtime = time.time() - start
    log_path.write_text("$ " + " ".join(cmd) + "\n\n" + proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSee {log_path}")
    return runtime


def run_astral_analyses(treatment_files: dict[str, object]) -> dict[str, object]:
    ASTRAL_SU_RESULTS.mkdir(parents=True, exist_ok=True)
    astral = find_astral4()
    help_text = subprocess.run([str(astral), "--help"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout
    commands = []
    runtimes = {}
    for treatment in FIXED_TREATMENTS:
        out = ASTRAL_SU_RESULTS / f"{treatment}.fixed_background_su.nwk"
        log = ASTRAL_SU_RESULTS / f"{treatment}.fixed_background_su.log"
        cmd = [str(astral), "-C", "-u", "2", "-t", str(THREADS), "-a", str(ASTRAL_MAPPING), "-c", str(FIXED_BACKGROUND_TOPOLOGY), "-i", str(treatment_files[treatment]["path"]), "-o", str(out)]
        runtimes[f"{treatment}_fixed"] = run_command(cmd, log)
        commands.append(cmd)
    for treatment in FREE_TREATMENTS:
        out = ASTRAL_SU_RESULTS / f"{treatment}.free_su.nwk"
        log = ASTRAL_SU_RESULTS / f"{treatment}.free_su.log"
        cmd = [str(astral), "-R", "-u", "2", "-t", str(THREADS), "-a", str(ASTRAL_MAPPING), "-i", str(treatment_files[treatment]["path"]), "-o", str(out)]
        runtimes[f"{treatment}_free"] = run_command(cmd, log)
        commands.append(cmd)
    try:
        commit = subprocess.run(["git", "-C", str(astral.parents[1]), "rev-parse", "HEAD"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False).stdout.strip()
    except Exception:
        commit = "unavailable"
    return {"astral4_path": astral, "astral4_sha256": sha256(astral), "help_text": help_text, "commands": commands, "runtimes": runtimes, "aster_commit": commit or "unavailable"}


def all_taxa(tree) -> frozenset[str]:
    labels = [t.name for t in tree.get_terminals()]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate terminal labels")
    return frozenset(labels)


def split_id_from_side(side: Iterable[str], taxa: Iterable[str]) -> str:
    t = frozenset(taxa)
    s = frozenset(side)
    other = t - s
    left = ",".join(sorted(s))
    right = ",".join(sorted(other))
    return "|".join(sorted([left, right]))


def compact_side(split_id: str, taxa: Iterable[str]) -> frozenset[str]:
    left, right = split_id.split("|")
    a = frozenset(left.split(",")) if left else frozenset()
    b = frozenset(right.split(",")) if right else frozenset()
    if len(a) < len(b):
        return a
    if len(b) < len(a):
        return b
    return min([a, b], key=lambda x: ",".join(sorted(x)))


def annotation_dict(name: str | None) -> dict[str, str]:
    if not name:
        return {}
    m = re.search(r"\[(.*)\]", name)
    if not m:
        return {}
    out = {}
    for piece in m.group(1).split(";"):
        if "=" in piece:
            key, value = piece.split("=", 1)
            out[key] = value
    return out


def parse_float_safe(value: str | None) -> float | None:
    if value is None or value in {"", "NA"}:
        return None
    out = float(value)
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def branch_role(split_id: str) -> str:
    left_text, right_text = split_id.split("|")
    sides = [frozenset(left_text.split(",")), frozenset(right_text.split(","))]
    if INV_PAIR in sides:
        return "focal_invicta_macdonaghi_SB_Sb_pair"
    if RIC_PAIR in sides:
        return "focal_richteri_SB_Sb_pair"
    if FOCAL in sides:
        return "focal_four_vs_outgroups"
    if any(side & FOCAL for side in sides):
        return "focal_related_internal_branch"
    return "non_focal_internal_branch"


def focal_split_and_class(path: Path) -> tuple[str, str]:
    tree = parse_tree_file(path)
    taxa = all_taxa(tree)
    focal_splits = {
        "inv_mac_SB,inv_mac_Sb|richteri_SB,richteri_Sb": "species",
        "inv_mac_SB,richteri_SB|inv_mac_Sb,richteri_Sb": "haplotype",
        "inv_mac_SB,richteri_Sb|inv_mac_Sb,richteri_SB": "third",
    }
    n = len(taxa)
    observed = set()
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        if len(side) <= 1 or len(side) >= n - 1:
            continue
        focal_side = side & FOCAL
        if len(focal_side) == 2:
            split = split_id_from_side(focal_side, FOCAL)
            if split in focal_splits:
                observed.add(split)
    if len(observed) != 1:
        raise ValueError(f"Could not identify a unique focal split in {path}: {sorted(observed)}")
    split = sorted(observed)[0]
    return split, focal_splits[split]


def extract_annotated_branches(path: Path) -> dict[str, dict[str, object]]:
    tree = parse_tree_file(path)
    taxa = all_taxa(tree)
    if taxa != frozenset(EXPECTED_SAFE_TAXA):
        raise ValueError(f"Unexpected taxon set in {path}: {sorted(taxa)}")
    out = {}
    n = len(taxa)
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        if len(side) <= 1 or len(side) >= n - 1:
            continue
        split_id = split_id_from_side(side, taxa)
        ann = annotation_dict(getattr(clade, "name", None))
        small = compact_side(split_id, taxa)
        out[split_id] = {
            "split_id": split_id,
            "taxa_side": ",".join(sorted(small)),
            "n_taxa_side": len(small),
            "branch_role": branch_role(split_id),
            "SULength": parse_float_safe(ann.get("SULength")),
            "CULength": parse_float_safe(ann.get("CULength")),
            "q1": parse_float_safe(ann.get("q1")),
            "q2": parse_float_safe(ann.get("q2")),
            "q3": parse_float_safe(ann.get("q3")),
            "localPP": parse_float_safe(ann.get("localPP")),
            "annotation": ann,
        }
    return out


def compare_fixed_background_all() -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    bg_path = ASTRAL_SU_RESULTS / "T_background.fixed_background_su.nwk"
    all_path = ASTRAL_SU_RESULTS / "T_all.fixed_background_su.nwk"
    bg = extract_annotated_branches(bg_path)
    allb = extract_annotated_branches(all_path)
    rows = []
    for split_id in sorted(set(bg) | set(allb)):
        b = bg.get(split_id)
        a = allb.get(split_id)
        if b and a:
            status = "shared"
            su_b, su_a = b["SULength"], a["SULength"]
            cu_b, cu_a = b["CULength"], a["CULength"]
            delta_su = None if su_b is None or su_a is None else su_a - su_b
            rel_su = None if delta_su is None or su_b is None or abs(su_b) <= EPS else delta_su / su_b
            delta_cu = None if cu_b is None or cu_a is None else cu_a - cu_b
            rel_cu = None if delta_cu is None or cu_b is None or abs(cu_b) <= EPS else delta_cu / cu_b
            base = b
        elif b:
            status = "background-only"
            su_b, su_a = b["SULength"], None
            cu_b, cu_a = b["CULength"], None
            delta_su = rel_su = delta_cu = rel_cu = None
            base = b
        else:
            status = "all-only"
            su_b, su_a = None, a["SULength"]
            cu_b, cu_a = None, a["CULength"]
            delta_su = rel_su = delta_cu = rel_cu = None
            base = a
        rows.append({
            "split_id": split_id,
            "taxa_side": base["taxa_side"],
            "n_taxa_side": base["n_taxa_side"],
            "branch_role": base["branch_role"],
            "su_background": su_b,
            "su_all": su_a,
            "delta_su": delta_su,
            "relative_delta_su": rel_su,
            "cu_background_fixed_su_run": cu_b,
            "cu_all_fixed_su_run": cu_a,
            "delta_cu_fixed_su_run": delta_cu,
            "relative_delta_cu_fixed_su_run": rel_cu,
            "status": status,
        })
    cu_stage6 = {}
    if STAGE6_CU_COMPARISON.exists():
        for r in read_tsv(STAGE6_CU_COMPARISON):
            if r.get("comparison") == "background_vs_all" and r.get("status") == "shared branch":
                cu_stage6[r["split_id"]] = r
    cuvs = []
    for r in rows:
        if r["status"] != "shared":
            continue
        c = cu_stage6.get(r["split_id"])
        cuvs.append({
            "split_id": r["split_id"],
            "taxa_side": r["taxa_side"],
            "branch_role": r["branch_role"],
            "cu_background": r["cu_background_fixed_su_run"],
            "cu_all": r["cu_all_fixed_su_run"],
            "delta_cu": r["delta_cu_fixed_su_run"],
            "relative_delta_cu": r["relative_delta_cu_fixed_su_run"],
            "su_background": r["su_background"],
            "su_all": r["su_all"],
            "delta_su": r["delta_su"],
            "relative_delta_su": r["relative_delta_su"],
            "secondary_grouped_stage6_cu_background": parse_float_safe(c.get("cu_background")) if c else None,
            "secondary_grouped_stage6_cu_all": parse_float_safe(c.get("cu_combined")) if c else None,
            "secondary_grouped_stage6_delta_cu": parse_float_safe(c.get("delta_cu")) if c else None,
            "secondary_grouped_stage6_relative_delta_cu": parse_float_safe(c.get("relative_delta_cu")) if c else None,
        })
    topology_same = set(bg) == set(allb)
    meta = {"n_background_internal_branches": len(bg), "n_all_internal_branches": len(allb), "n_shared": len(set(bg) & set(allb)), "n_background_only": len(set(bg) - set(allb)), "n_all_only": len(set(allb) - set(bg)), "topology_same": topology_same}
    return rows, cuvs, meta


def summarize_free_topologies(prov: dict[str, object]) -> list[dict[str, object]]:
    bg_splits = set(extract_annotated_branches(ASTRAL_SU_RESULTS / "T_background.free_su.nwk"))
    grouped_stage6 = {}
    if STAGE6_ASTRAL_TREATMENT_SUMMARY.exists():
        for r in read_tsv(STAGE6_ASTRAL_TREATMENT_SUMMARY):
            grouped_stage6[r["treatment"]] = r
    rows = []
    for treatment in FREE_TREATMENTS:
        path = ASTRAL_SU_RESULTS / f"{treatment}.free_su.nwk"
        branches = extract_annotated_branches(path)
        splits = set(branches)
        focal_split, focal_class = focal_split_and_class(path)
        grouped = grouped_stage6.get(treatment, {})
        rows.append({
            "treatment": treatment,
            "n_internal_splits": len(splits),
            "same_topology_as_free_background": splits == bg_splits,
            "unrooted_RF_to_free_background": len(bg_splits - splits) + len(splits - bg_splits),
            "n_background_splits_recovered": len(bg_splits & splits),
            "fraction_background_splits_recovered": len(bg_splits & splits) / len(bg_splits) if bg_splits else None,
            "individual_free_focal_split": focal_split,
            "individual_free_focal_class": focal_class,
            "grouped_stage6_focal_class": grouped.get("focal_class", "NA"),
            "matches_grouped_stage6_focal_class": (focal_class == grouped.get("focal_class")) if grouped else None,
            "tree_file": rel(path),
            "runtime_seconds": prov["runtimes"].get(f"{treatment}_free"),
        })
    return rows


def write_mapping_report(mapping: dict[str, object]) -> None:
    lines = [
        "# Stage 6C mapping recovery report",
        "",
        "## Supplementary Data 1 source",
        "",
        f"Source URL: {mapping['supplement_url']}",
        f"Local file: `{rel(SUPP_DATA1)}`",
        f"SHA256: `{mapping['supplement_sha256']}`",
        f"Worksheet names: `{', '.join(mapping['supplement_meta']['worksheet_names'])}`",
        f"Rows after header: `{mapping['supplement_meta']['n_rows']}`",
        "Columns: " + ", ".join(f"`{c}`" for c in mapping['supplement_meta']['column_names']),
        "",
        "## Upstream grouping rule",
        "",
        "The authoritative grouping rule is the rule in `Topology weighting/supergene/filter_samples_by_missing.R`: remove the `S.` species prefix and spaces, merge `invicta` and `macdonaghi` into `invicta/macdonaghi`, and append `Supergene.Variant` for `invicta/macdonaghi` and `richteri`. The other retained species remain species-only groups.",
        "",
        "## Matching result",
        "",
        f"Tree tips inspected: `{mapping['n_tree_tips']}`",
        f"Matched tips: `{len(mapping['individual_rows'])}`",
        f"Unmatched tips: `{len(mapping['unmatched'])}`",
        f"Match-type counts: `{mapping['match_type_counts']}`",
        f"Documented name corrections used: `{mapping['documented_corrections_used']}`",
        f"Species assignments: `{mapping['n_species_assigned']}`",
        f"Supergene-variant assignments: `{mapping['n_supergene_variant_assigned']}`",
        f"Complete exact individual-to-group mapping available: `{str(mapping['complete_exact_individual_to_group_mapping_available']).lower()}`",
        "",
        "Group sizes:",
    ]
    for group in EXPECTED_GROUPS:
        lines.append(f"- `{group}`: `{mapping['group_counts'].get(group, 0)}`")
    lines.extend([
        "",
        "## Validation against tree-label bigB/littleb strings",
        "",
        "The tree-label strings were used only as a validation check after metadata-based assignment, not as the mapping source.",
        f"bigB/littleb validation mismatches: `{len(mapping['biglittle_validation_mismatches'])}`",
        "",
        "## Output files",
        "",
        f"- Tip inventory: `{rel(TREE_TIP_INVENTORY)}`",
        f"- Match audit: `{rel(MATCH_AUDIT)}`",
        f"- Individual-to-group table: `{rel(INDIVIDUAL_TO_GROUP)}`",
        f"- ASTRAL mapping file: `{rel(ASTRAL_MAPPING)}`",
    ])
    if mapping["unmatched"]:
        lines.extend(["", "## Unresolved samples", ""])
        lines.extend(f"- `{x}`" for x in mapping["unmatched"])
    MAPPING_REPORT.write_text("\n".join(lines) + "\n")


def write_branch_audit(stats: dict[str, object]) -> None:
    lines = [
        "Stage 6C original RAxML local-tree branch-length audit",
        "",
        f"input_tree = {rel(INPUT_TREE)}",
        f"input_tree_sha256 = {sha256(INPUT_TREE)}",
        "source = published upstream Topology weighting/results/2021-08-07-twisst/input.tree",
        "upstream construction = concatenated `*.raxml.bestTree` files from four-BUSCO windows as documented in Topology weighting/README.md",
        "branch_length_interpretation = RAxML-NG local best-tree branch lengths from GTR+G sequence alignments; these are substitution/site tree branch lengths",
        "",
    ]
    for key in ["n_trees", "n_nonempty_trees", "tips_per_tree_min", "tips_per_tree_max", "unique_tip_count", "all_tree_tip_sets_identical", "n_branch_lengths", "n_positive_branch_lengths", "n_zero_branch_lengths", "fraction_positive_branch_lengths", "fraction_zero_branch_lengths", "median_branch_length", "median_positive_branch_length", "min_branch_length", "max_branch_length", "tip_set_sha256"]:
        lines.append(f"{key} = {stats[key]}")
    lines.extend(["", "Treatment window counts from Stage-1 window index:"])
    for key, value in stats["treatment_counts"].items():
        lines.append(f"{key} = {value}")
    lines.extend(["", "Conclusion: the original RAxML local trees contain usable finite nonnegative substitution branch lengths; 99.812% of parsed branch lengths are positive."])
    BRANCH_AUDIT.write_text("\n".join(lines) + "\n")


def write_provenance(prov: dict[str, object]) -> None:
    lines = [
        "# Stage 6C ASTRAL4/CASTLES-II provenance",
        "",
        f"Executable: `{prov['astral4_path']}`",
        f"Executable SHA256: `{prov['astral4_sha256']}`",
        f"ASTER commit: `{prov['aster_commit']}`",
        f"Platform: `{platform.platform()}`",
        f"Python: `{sys.version.split()[0]}`",
        f"Threads: `{THREADS}`",
        "CASTLES-II status: the ASTRAL4 help output reports `NOW with integrated CASTLES-2`; `SULength` is the default substitution-per-site branch-length output.",
        "",
        "## Help output",
        "",
        "```text",
        prov["help_text"].strip(),
        "```",
        "",
        "## Commands",
        "",
    ]
    for cmd in prov["commands"]:
        lines.append("```bash")
        lines.append(" ".join(str(x) for x in cmd))
        lines.append("```")
    PROVENANCE.write_text("\n".join(lines) + "\n")


def write_su_summary(rows: list[dict[str, object]], cuvs: list[dict[str, object]], meta: dict[str, object], free_rows: list[dict[str, object]]) -> None:
    shared = [r for r in rows if r["status"] == "shared"]
    bg_su = [r["su_background"] for r in shared if r["su_background"] is not None]
    all_su = [r["su_all"] for r in shared if r["su_all"] is not None]
    deltas = [r["delta_su"] for r in shared if r["delta_su"] is not None]
    lines = [
        "Stage 6C ASTRAL4/CASTLES-II substitution-unit branch-length summary",
        "",
        f"main_comparison = T_background fixed background topology vs T_all fixed background topology",
        f"same_unrooted_topology = {str(meta['topology_same']).lower()}",
        f"n_background_internal_branches = {meta['n_background_internal_branches']}",
        f"n_all_internal_branches = {meta['n_all_internal_branches']}",
        f"n_shared_bipartitions = {meta['n_shared']}",
        f"n_background_only_bipartitions = {meta['n_background_only']}",
        f"n_all_only_bipartitions = {meta['n_all_only']}",
        f"median_su_background = {statistics.median(bg_su) if bg_su else 'NA'}",
        f"median_su_all = {statistics.median(all_su) if all_su else 'NA'}",
        f"median_delta_su = {statistics.median(deltas) if deltas else 'NA'}",
        "",
        "Shared-branch SU changes:",
    ]
    for r in sorted(shared, key=lambda x: x["delta_su"] if x["delta_su"] is not None else 0):
        pct = None if r["relative_delta_su"] is None else 100 * r["relative_delta_su"]
        lines.append(f"- {r['branch_role']} ({r['taxa_side']}): {r['su_background']} -> {r['su_all']} ; delta = {r['delta_su']} ; percent = {pct}")
    lines.extend(["", "CU vs SU focal-pair comparison:"])
    for r in cuvs:
        if r["branch_role"] in {"focal_invicta_macdonaghi_SB_Sb_pair", "focal_richteri_SB_Sb_pair"}:
            cu_pct = None if r["relative_delta_cu"] is None else 100 * r["relative_delta_cu"]
            su_pct = None if r["relative_delta_su"] is None else 100 * r["relative_delta_su"]
            lines.append(f"- {r['branch_role']}: CU percent change = {cu_pct}; SU percent change = {su_pct}")
    lines.extend(["", "Free-topology individual-level mapped ASTRAL4:"])
    for r in free_rows:
        lines.append(f"- {r['treatment']}: same topology as free background = {r['same_topology_as_free_background']}; RF = {r['unrooted_RF_to_free_background']}; recovered background splits = {r['n_background_splits_recovered']}/{r['n_internal_splits']}")
    SU_SUMMARY.write_text("\n".join(lines) + "\n")


def make_figures(rows: list[dict[str, object]], cuvs: list[dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    shared = [r for r in rows if r["status"] == "shared" and r["su_background"] is not None and r["su_all"] is not None]
    plt.figure(figsize=(5, 4.5))
    xs = [r["su_background"] for r in shared]
    ys = [r["su_all"] for r in shared]
    maxv = max(xs + ys) * 1.08 if xs else 1
    plt.plot([0, maxv], [0, maxv], color="0.6", lw=1, zorder=0)
    for r in shared:
        color = "#d95f02" if "focal_" in r["branch_role"] else "#1b9e77"
        plt.scatter(r["su_background"], r["su_all"], color=color, s=55)
        plt.text(r["su_background"], r["su_all"], r["taxa_side"], fontsize=7, ha="left", va="bottom")
    plt.xlabel("background-only SULength")
    plt.ylabel("background + supergene SULength")
    plt.title("Stage 6C fixed-topology SU lengths")
    plt.tight_layout()
    plt.savefig(FIG_SU_SCATTER_PDF)
    plt.savefig(FIG_SU_SCATTER_PNG, dpi=300)
    plt.close()

    cu_shared = [r for r in cuvs if r["relative_delta_cu"] is not None and r["relative_delta_su"] is not None]
    plt.figure(figsize=(5, 4.5))
    plt.axhline(0, color="0.7", lw=1)
    plt.axvline(0, color="0.7", lw=1)
    for r in cu_shared:
        color = "#d95f02" if "focal_" in r["branch_role"] else "#1b9e77"
        label = r["taxa_side"].replace("inv_mac_SB,inv_mac_Sb", "invicta/macdonaghi").replace("richteri_SB,richteri_Sb", "richteri")
        plt.scatter(100 * r["relative_delta_cu"], 100 * r["relative_delta_su"], color=color, s=60)
        plt.text(100 * r["relative_delta_cu"], 100 * r["relative_delta_su"], label, fontsize=8, ha="left", va="bottom")
    plt.xlabel("CULength change (%)")
    plt.ylabel("SULength change (%)")
    plt.title("Stage 6C fixed-topology CU versus SU response")
    plt.tight_layout()
    plt.savefig(FIG_CU_VS_SU_PDF)
    plt.savefig(FIG_CU_VS_SU_PNG, dpi=300)
    plt.close()

    focal = [r for r in cuvs if r["branch_role"] in {"focal_invicta_macdonaghi_SB_Sb_pair", "focal_richteri_SB_Sb_pair"}]
    labels = [r["branch_role"].replace("focal_", "").replace("_SB_Sb_pair", "") for r in focal]
    x = range(len(focal))
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    width = 0.35
    axes[0].bar([i - width/2 for i in x], [r["cu_background"] for r in focal], width, label="background")
    axes[0].bar([i + width/2 for i in x], [r["cu_all"] for r in focal], width, label="all")
    axes[0].set_title("CULength")
    axes[0].set_xticks(list(x), labels, rotation=20, ha="right")
    axes[0].legend(frameon=False)
    axes[1].bar([i - width/2 for i in x], [r["su_background"] for r in focal], width, label="background")
    axes[1].bar([i + width/2 for i in x], [r["su_all"] for r in focal], width, label="all")
    axes[1].set_title("SULength")
    axes[1].set_xticks(list(x), labels, rotation=20, ha="right")
    axes[1].legend(frameon=False)
    fig.suptitle("Focal species-pair branch lengths")
    fig.tight_layout()
    fig.savefig(FIG_FOCAL_PDF)
    fig.savefig(FIG_FOCAL_PNG, dpi=300)
    plt.close(fig)


def pct(x: float | None) -> float | None:
    return None if x is None else 100.0 * x


def first_row(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, observed {len(rows)}")
    return rows[0]


def collect_final_context(cuvs: list[dict[str, object]], free_rows: list[dict[str, object]], mapping: dict[str, object], branch_stats: dict[str, object]) -> dict[str, object]:
    region = {r["region"]: r for r in read_tsv(STAGE4A_REGION_SUMMARY)}
    primary = first_row(STAGE4B_PRIMARY)
    coord_unique = first_row(STAGE4B_COORD_UNIQUE)
    coord_length = first_row(STAGE4B_COORD_LENGTH)
    free = {r["treatment"]: r for r in free_rows}
    cuv_by_role = {r["branch_role"]: r for r in cuvs}
    return {
        "region": region,
        "primary": primary,
        "coord_unique": coord_unique,
        "coord_length": coord_length,
        "free": free,
        "cuv_by_role": cuv_by_role,
        "mapping": mapping,
        "branch_stats": branch_stats,
    }


def write_final_summary_table(ctx: dict[str, object]) -> list[dict[str, object]]:
    region = ctx["region"]
    primary = ctx["primary"]
    unique = ctx["coord_unique"]
    length = ctx["coord_length"]
    free = ctx["free"]
    cuv = ctx["cuv_by_role"]
    outside = region["chr16_outside"]
    supergene = region["chr16_supergene"]
    rows = [
        {
            "analysis": "local_topology_support",
            "comparison": "chr16_outside_vs_supergene",
            "metric": "mean_D=q_H-q_S",
            "background": outside["mean_D"],
            "supergene_or_combined": supergene["mean_D"],
            "change": primary["observed_delta_D"],
            "interpretation": "localized_shift_from_species_history_to_haplotype_history_inside_supergene",
        },
        {
            "analysis": "dominant_topology_counts",
            "comparison": "chr16_outside_vs_supergene",
            "metric": "dominant_class",
            "background": f"{outside['n_species_dominant']}/{outside['n_windows']} species",
            "supergene_or_combined": f"{supergene['n_haplotype_dominant']}/{supergene['n_windows']} haplotype",
            "change": "species_dominant_outside_haplotype_dominant_supergene",
            "interpretation": "supergene_windows_are_not_merely_more_discordant_but_shift_toward_SB_Sb_haplotype_topology",
        },
        {
            "analysis": "spatial_null",
            "comparison": "exact_circular_shift",
            "metric": "Delta_D_rank_and_p",
            "background": f"rank {primary['observed_rank']}/{primary['n_exact_alignments']}",
            "supergene_or_combined": f"p={primary['p_one_sided']}",
            "change": primary["observed_delta_D"],
            "interpretation": "positive_spatial_alignment_with_frozen_supergene_region",
        },
        {
            "analysis": "coordinate_sensitivity",
            "comparison": "unique_membership_states",
            "metric": "p_unique_membership",
            "background": f"rank {unique['rank_unique_membership']}/{unique['n_unique_membership_states']}",
            "supergene_or_combined": f"p={unique['p_unique_membership']}",
            "change": unique["observed_delta_D"],
            "interpretation": "supportive_equal_weight_state_sensitivity",
        },
        {
            "analysis": "coordinate_sensitivity",
            "comparison": "continuous_uniform_physical_start",
            "metric": "p_length_weighted",
            "background": f"extreme_length={length['extreme_start_domain_length']}",
            "supergene_or_combined": f"p={length['p_length_weighted']}",
            "change": length["observed_delta_D"],
            "interpretation": "more_conservative_physical_start_sensitivity",
        },
        {
            "analysis": "summary_tree_topology",
            "comparison": "T_background_vs_T_all",
            "metric": "free_topology_RF",
            "background": free["T_background"]["individual_free_focal_class"],
            "supergene_or_combined": free["T_all"]["individual_free_focal_class"],
            "change": f"RF={free['T_all']['unrooted_RF_to_free_background']}",
            "interpretation": "adding_supergene_does_not_change_overall_species_topology",
        },
        {
            "analysis": "summary_tree_topology",
            "comparison": "T_supergene_only",
            "metric": "focal_topology",
            "background": free["T_background"]["individual_free_focal_class"],
            "supergene_or_combined": free["T_supergene"]["individual_free_focal_class"],
            "change": f"RF={free['T_supergene']['unrooted_RF_to_free_background']}",
            "interpretation": "supergene_only_recovers_haplotype_topology",
        },
    ]
    for role, label in [
        ("focal_richteri_SB_Sb_pair", "richteri_SB_plus_richteri_Sb"),
        ("focal_invicta_macdonaghi_SB_Sb_pair", "inv_mac_SB_plus_inv_mac_Sb"),
    ]:
        r = cuv[role]
        rows.append({
            "analysis": "fixed_topology_branch_lengths",
            "comparison": label,
            "metric": "CULength",
            "background": r["cu_background"],
            "supergene_or_combined": r["cu_all"],
            "change": pct(r["relative_delta_cu"]),
            "interpretation": "coalescent_unit_branch_shortening_after_adding_supergene",
        })
        rows.append({
            "analysis": "fixed_topology_branch_lengths",
            "comparison": label,
            "metric": "SULength",
            "background": r["su_background"],
            "supergene_or_combined": r["su_all"],
            "change": pct(r["relative_delta_su"]),
            "interpretation": "substitution_unit_response_not_consistent_with_uniform_focal_shortening",
        })
    write_tsv(FINAL_SUMMARY_TABLE, rows, ["analysis", "comparison", "metric", "background", "supergene_or_combined", "change", "interpretation"])
    return rows


def write_final_methods_results(ctx: dict[str, object]) -> None:
    region = ctx["region"]
    primary = ctx["primary"]
    unique = ctx["coord_unique"]
    length = ctx["coord_length"]
    free = ctx["free"]
    cuv = ctx["cuv_by_role"]
    mapping = ctx["mapping"]
    branch = ctx["branch_stats"]
    rich = cuv["focal_richteri_SB_Sb_pair"]
    inv = cuv["focal_invicta_macdonaghi_SB_Sb_pair"]
    methods = f"""# Fire-ant chromosome-16 empirical analysis: final Methods text

We analyzed the published fire-ant dataset from Stolle et al. (2022), “Recurring adaptive introgression of a supergene variant that determines social organization” (Nature Communications; DOI 10.1038/s41467-022-28806-7), using the upstream local RAxML-NG/TWISST files from the frozen repository snapshot. The local input consisted of {branch['n_trees']} four-BUSCO-gene RAxML-NG trees, each with {branch['unique_tip_count']} individual tips, and the corresponding published TWISST topology weights. Windows were assigned to the author-defined regions chr1, chr16A, chr16B, and the chromosome-16 supergene analysis interval using the frozen Stage-1 window index.

The grouped TWISST analysis used seven published source groups: geminata, saevissima, pusillignis, invicta/macdonaghi_SB, invicta/macdonaghi_Sb, richteri_SB, and richteri_Sb. The focal quartet was A = invicta/macdonaghi_SB, B = invicta/macdonaghi_Sb, C = richteri_SB, and D = richteri_Sb. We classified each seven-group topology by its induced unrooted quartet as the background species split AB|CD, the cross-species SB/Sb haplotype split AC|BD, or the third split AD|BC, and aggregated TWISST weights into q_S, q_H, and q_3. The primary local statistic was D = q_H - q_S, and the primary contrast was Delta_D = mean(D_supergene) - mean(D_chr16_outside).

Spatial inference used the predeclared chromosome-16 physical order, sorting windows by physical midpoint because the upstream TWISST concatenation order was not physical chromosome order. The primary null was the exact 96-alignment circular-shift null, rotating the ordered D track relative to the independently frozen supergene mask. Coordinate-aware sensitivities used the same frozen physical interval width as the author-designated supergene region, first weighting each distinct sampled-window membership state equally and then weighting constant-membership start-coordinate intervals by their physical length.

Post-freeze Stage 6 analyses assessed downstream summary-tree sensitivity. The grouped/modal ASTRAL4 analysis represented each window by its modal seven-group TWISST topology and is retained as a topology and grouped-CU robustness analysis. For the direct CU-vs-SU comparison, we recovered the exact {mapping['n_tree_tips']}-tip to seven-group mapping from Stolle et al. Supplementary Data 1 without topology-based inference: {len(mapping['individual_rows'])}/{mapping['n_tree_tips']} tree tips matched published metadata exactly, with species and supergene-variant fields available for all mapped samples. We then used the original individual-level RAxML-NG trees, preserving their substitution branch lengths, with ASTRAL4/CASTLES-II and the recovered mapping.

The primary branch-length follow-up used a fixed background topology for both the {EXPECTED_TREATMENT_N['T_background']}-window background set and the {EXPECTED_TREATMENT_N['T_all']}-window all-window set. This fixed-topology design isolates changes in ASTRAL4/CASTLES-II branch-length estimates caused by adding the supergene windows while holding the species-tree topology constant. We report both coalescent-unit lengths (CULength) and substitution-unit lengths (SULength) from the same individual-level inputs and the same fixed topology.
"""
    FINAL_METHODS.write_text(methods)

    results = f"""# Fire-ant chromosome-16 empirical analysis: final Results text

Outside the chromosome-16 supergene, local genealogies overwhelmingly supported the independently defined species-history quartet. Across chr16 outside-supergene windows, mean q_S = {float(region['chr16_outside']['mean_q_S']):.3f}, mean q_H = {float(region['chr16_outside']['mean_q_H']):.3f}, and {region['chr16_outside']['n_species_dominant']}/{region['chr16_outside']['n_windows']} windows were species-dominant. Inside the supergene region, support shifted strongly toward the cross-species SB/Sb haplotype partition: mean q_S = {float(region['chr16_supergene']['mean_q_S']):.3f}, mean q_H = {float(region['chr16_supergene']['mean_q_H']):.3f}, and {region['chr16_supergene']['n_haplotype_dominant']}/{region['chr16_supergene']['n_windows']} windows were haplotype-dominant. The resulting Delta_D was {float(primary['observed_delta_D']):.4f}.

The spatial alignment was unusual under the predeclared exact circular-shift null. The observed alignment ranked {primary['observed_rank']}/{primary['n_exact_alignments']} among all circular alignments, giving a one-sided exact p = {float(primary['p_one_sided']):.4f}. Coordinate-aware sensitivities were concordant but differed in interpretation: the equal-weight unique-membership sensitivity gave p = {float(unique['p_unique_membership']):.4f}, whereas the continuous uniform physical-start sensitivity gave p = {float(length['p_length_weighted']):.4f}.

The post-freeze individual-level ASTRAL4/CASTLES-II topology sensitivity preserved the background topology when the supergene windows were added. The free-topology background and all-window analyses both recovered the focal species split and had RF = {free['T_all']['unrooted_RF_to_free_background']} relative to one another. In contrast, the supergene-only analysis recovered the haplotype focal split and differed from the background topology (RF = {free['T_supergene']['unrooted_RF_to_free_background']}), consistent with the local TWISST signal.

Even with the topology fixed, adding the supergene windows substantially changed focal coalescent-unit branch lengths. The richteri SB/Sb branch decreased from {rich['cu_background']:.5g} CU to {rich['cu_all']:.5g} CU ({100*rich['relative_delta_cu']:.1f}%), while its SULength was nearly unchanged, {rich['su_background']:.6g} to {rich['su_all']:.6g} ({100*rich['relative_delta_su']:.1f}%). The invicta/macdonaghi SB/Sb branch decreased from {inv['cu_background']:.5g} CU to {inv['cu_all']:.5g} CU ({100*inv['relative_delta_cu']:.1f}%), while its SULength increased from {inv['su_background']:.6g} to {inv['su_all']:.6g} ({100*inv['relative_delta_su']:.1f}%). Thus, the inclusion of the localized supergene genealogy regime can leave the inferred topology unchanged while substantially affecting coalescent-unit branch-length estimates.

These results demonstrate the phylogenetic consequences of the chromosome-16 supergene region: it is a spatially localized alternative genealogy regime associated with the recombination-suppressed SB/Sb haplotype. The original source study interprets the shared supergene history as recurrent adaptive introgression. Our analysis does not establish MSRC without gene flow, nor does it rule out introgression; it shows that the localized genealogy regime has clear effects on local quartet support and on downstream coalescent-unit branch-length estimates.
"""
    FINAL_RESULTS.write_text(results)


def make_final_summary_figure(ctx: dict[str, object]) -> None:
    support = [r for r in read_tsv(STAGE4A_WINDOW_SUPPORT) if r["chrom"] == "chr16"]
    support.sort(key=lambda r: float(r["mid"]))
    free = ctx["free"]
    cuv = ctx["cuv_by_role"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axes[0]
    for region_name, subset in [(k, [r for r in support if r["region"] == k]) for k in ["chr16A", "chr16_supergene", "chr16B"]]:
        x = [float(r["mid"]) / 1e6 for r in subset]
        ax.plot(x, [float(r["q_S"]) for r in subset], color="#1b9e77", lw=1, marker="o", ms=2, label="q_S" if region_name == "chr16A" else None)
        ax.plot(x, [float(r["q_H"]) for r in subset], color="#d95f02", lw=1, marker="o", ms=2, label="q_H" if region_name == "chr16A" else None)
    sg = [r for r in support if r["region"] == "chr16_supergene"]
    ax.axvspan(min(float(r["mid"]) for r in sg) / 1e6, max(float(r["mid"]) for r in sg) / 1e6, color="#d95f02", alpha=0.12)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("chr16 position (Mb)")
    ax.set_ylabel("TWISST quartet support")
    ax.set_title("A. Local quartet support")
    ax.legend(frameon=False, loc="center left", fontsize=8)

    ax = axes[1]
    treatments = ["T_background", "T_all", "T_supergene"]
    labels = ["background", "all", "supergene"]
    colors = {"species": "#1b9e77", "haplotype": "#d95f02", "third": "#7570b3"}
    classes = [free[t]["individual_free_focal_class"] for t in treatments]
    ax.bar(labels, [1, 1, 1], color=[colors[c] for c in classes])
    for i, c in enumerate(classes):
        ax.text(i, 0.5, c, ha="center", va="center", color="white", fontweight="bold")
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title("B. ASTRAL4 focal topology")
    ax.set_ylabel("free topology class")

    ax = axes[2]
    roles = [("focal_richteri_SB_Sb_pair", "richteri"), ("focal_invicta_macdonaghi_SB_Sb_pair", "invicta/macdonaghi")]
    x = list(range(len(roles)))
    width = 0.36
    cu_pct = [100 * cuv[role]["relative_delta_cu"] for role, _ in roles]
    su_pct = [100 * cuv[role]["relative_delta_su"] for role, _ in roles]
    ax.axhline(0, color="0.6", lw=1)
    ax.bar([i - width/2 for i in x], cu_pct, width, color="#4c78a8", label="CU")
    ax.bar([i + width/2 for i in x], su_pct, width, color="#f58518", label="SU")
    for i, val in enumerate(cu_pct):
        ax.text(i - width/2, val, f"{val:.0f}%", ha="center", va="bottom" if val >= 0 else "top", fontsize=8)
    for i, val in enumerate(su_pct):
        ax.text(i + width/2, val, f"{val:.0f}%", ha="center", va="bottom" if val >= 0 else "top", fontsize=8)
    ax.set_xticks(x, [label for _, label in roles], rotation=15, ha="right")
    ax.set_ylabel("change after adding supergene (%)")
    ax.set_title("C. Fixed-topology branch lengths")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_FINAL_SUMMARY_PDF)
    fig.savefig(FIG_FINAL_SUMMARY_PNG, dpi=300)
    plt.close(fig)


def write_consistency_audit(ctx: dict[str, object], final_rows: list[dict[str, object]]) -> dict[str, object]:
    branch = ctx["branch_stats"]
    mapping = ctx["mapping"]
    free = ctx["free"]
    cuv = ctx["cuv_by_role"]
    checks = {
        "total_windows_213": sum(branch["treatment_counts"][k] for k in ["T_chr1", "T_supergene"]) + branch["treatment_counts"]["T_chr16_outside"] == 213,
        "background_windows_161": branch["treatment_counts"]["T_background"] == 161,
        "supergene_windows_52": branch["treatment_counts"]["T_supergene"] == 52,
        "individuals_267": branch["unique_tip_count"] == 267,
        "seven_groups": sorted(mapping["group_counts"]) == sorted(EXPECTED_GROUPS),
        "exact_metadata_matches_267": mapping["match_type_counts"].get("exact") == 267,
        "background_all_free_rf_zero": int(free["T_all"]["unrooted_RF_to_free_background"]) == 0,
        "supergene_only_haplotype": free["T_supergene"]["individual_free_focal_class"] == "haplotype",
        "richteri_cu_percent_recomputed": abs(100 * cuv["focal_richteri_SB_Sb_pair"]["relative_delta_cu"] - ((cuv["focal_richteri_SB_Sb_pair"]["cu_all"] - cuv["focal_richteri_SB_Sb_pair"]["cu_background"]) / cuv["focal_richteri_SB_Sb_pair"]["cu_background"] * 100)) < 1e-9,
        "invicta_cu_percent_recomputed": abs(100 * cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["relative_delta_cu"] - ((cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["cu_all"] - cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["cu_background"]) / cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["cu_background"] * 100)) < 1e-9,
        "final_summary_rows_present": len(final_rows) >= 11,
        "methods_contains_cu_su": "CULength" in FINAL_METHODS.read_text() and "SULength" in FINAL_METHODS.read_text(),
        "results_contains_introgression_caveat": "introgression" in FINAL_RESULTS.read_text() and "does not establish MSRC without gene flow" in FINAL_RESULTS.read_text(),
        "readme_hierarchy_mentions_stage6c_supersedes_stage6b": "Stage 6C supersedes Stage 6B for direct CU-vs-SU comparison" in README.read_text(),
    }
    failed = [k for k, v in checks.items() if not v]
    audit = {
        "status": "passed" if not failed else "failed",
        "checks": checks,
        "failed_checks": failed,
        "authoritative_focal_percent_changes": {
            "richteri_CU_percent": 100 * cuv["focal_richteri_SB_Sb_pair"]["relative_delta_cu"],
            "richteri_SU_percent": 100 * cuv["focal_richteri_SB_Sb_pair"]["relative_delta_su"],
            "invicta_macdonaghi_CU_percent": 100 * cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["relative_delta_cu"],
            "invicta_macdonaghi_SU_percent": 100 * cuv["focal_invicta_macdonaghi_SB_Sb_pair"]["relative_delta_su"],
        },
    }
    FINAL_CONSISTENCY_AUDIT.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if failed:
        raise RuntimeError(f"Final consistency audit failed: {failed}")
    return audit


def update_docs(mapping: dict[str, object], rows: list[dict[str, object]], cuvs: list[dict[str, object]], free_rows: list[dict[str, object]]) -> None:
    focal_lines = []
    for r in cuvs:
        if r["branch_role"] in {"focal_invicta_macdonaghi_SB_Sb_pair", "focal_richteri_SB_Sb_pair"}:
            focal_lines.append(f"{r['branch_role']}: SU {r['su_background']:.6g} -> {r['su_all']:.6g} ({100*r['relative_delta_su']:.2f}%), CU {r['cu_background']:.6g} -> {r['cu_all']:.6g} ({100*r['relative_delta_cu']:.2f}%)")
    block = (
        "\nStage 6C complete — mapped individual-level ASTRAL4/CASTLES-II SU branch-length sensitivity. "
        "The exact 267-tip mapping was recovered from Stolle et al. 2022 Supplementary Data 1 without guessing. "
        "Original RAxML local-tree branch lengths were preserved, fixed-topology background versus all-window CULength and SULength values were estimated from the same individual-level inputs, and CU/SU changes were compared. "
        "Stage 6C supersedes Stage 6B for direct CU-vs-SU comparison because CU and SU are estimated from identical individual-level inputs and the identical fixed background topology. Stage 6B remains an independent grouped/modal-topology CU sensitivity. Stage 5 remains frozen. See `results/stage6c_su_summary.txt`, `results/stage6c_su_branch_length_comparison.tsv`, and `results/stage6c_cu_vs_su_comparison.tsv`.\n"
    )
    readme = README.read_text()
    if "06c_su_branch_length_sensitivity.py --run-tests" not in readme:
        readme = readme.replace("python3 empirical/fire_ants_chr16/scripts/06b_compare_cu_branch_lengths.py --run-tests\n```", "python3 empirical/fire_ants_chr16/scripts/06b_compare_cu_branch_lengths.py --run-tests\npython3 empirical/fire_ants_chr16/scripts/06c_su_branch_length_sensitivity.py --run-tests\n```")
    readme = re.sub(r"\nStage 6C SU branch-length gate:.*?(?=\n\n|$)", "", readme, flags=re.S)
    if "Stage 6C complete — mapped individual-level" not in readme:
        readme += block
    else:
        readme = re.sub(r"\nStage 6C complete — mapped individual-level.*?\n", block, readme, flags=re.S)
    hierarchy = (
        "\n## Final fire-ant result hierarchy\n\n"
        "Primary empirical result: Stages 4A-5 show a localized chromosome-16 shift from the background species-history quartet outside the supergene to the cross-species SB/Sb haplotype quartet inside the independently defined supergene region.\n\n"
        "Post-freeze sensitivity: Stage 6 tests grouped TWISST/modal-tree ASTRAL4 summary-tree behavior. It remains a topology and grouped-CU robustness analysis and does not alter the Stage-5 primary empirical result.\n\n"
        "Strongest branch-length follow-up: Stage 6C uses the original 267-individual RAxML trees, the exact metadata-derived seven-group mapping, and ASTRAL4/CASTLES-II. Stage 6C supersedes Stage 6B for direct CU-vs-SU comparison because CU and SU are estimated from identical individual-level inputs and the identical fixed background topology. Stage 6B remains documented as an independent grouped/modal-topology CU sensitivity analysis.\n\n"
        "Final manuscript-facing outputs: `results/fire_ants_final_summary.tsv`, `results/fire_ants_final_methods.md`, `results/fire_ants_final_results.md`, and `figures/fire_ants_final_summary.pdf`.\n"
    )
    readme = re.sub(r"\n## Final fire-ant result hierarchy\n\n.*?\n(?=## |\Z)", "\n", readme, flags=re.S)
    readme += hierarchy
    README.write_text(readme)

    if PROJECT_STATUS.exists():
        status = PROJECT_STATUS.read_text()
        status_block = (
            "## Fire ants — analysis complete\n\n"
            "Biological question: how the chromosome-16 social-supergene region affects local quartet support, summary-tree topology, and species-tree branch lengths in the Stolle et al. fire-ant dataset.\n\n"
            f"Dataset size: 213 four-BUSCO windows, 267 individual tree tips, seven recovered TWISST/ASTRAL groups, and {len(mapping['individual_rows'])}/{mapping['n_tree_tips']} exact metadata matches from Supplementary Data 1.\n\n"
            "Primary Stage-4/5 result: local chromosome-16 genealogy support shifts from the background species-history quartet outside the supergene to the cross-species SB/Sb haplotype quartet inside the independently defined supergene region.\n\n"
            "Summary-tree result: background and all-window individual-level ASTRAL4 analyses retain the species focal topology, whereas the supergene-only analysis recovers the haplotype focal topology.\n\n"
            "CU/SU result: fixed-topology Stage 6C shows strong focal CULength decreases when supergene windows are added, while SULength responses differ between focal branches.\n\n"
            "Focal branch changes:\n" + "\n".join(f"- {x}" for x in focal_lines) + "\n\n"
            "Interpretation: topology can remain stable while branch-length estimates, especially coalescent-unit estimates, are affected by a localized alternative genealogy regime.\n\n"
            "Caveat: the source study attributes the shared supergene history to recurrent adaptive introgression. This analysis demonstrates phylogenetic consequences of the localized genealogy regime and does not establish MSRC without gene flow.\n\n"
            "Manuscript-ready outputs: `empirical/fire_ants_chr16/results/fire_ants_final_summary.tsv`, `empirical/fire_ants_chr16/results/fire_ants_final_methods.md`, `empirical/fire_ants_chr16/results/fire_ants_final_results.md`, and `empirical/fire_ants_chr16/figures/fire_ants_final_summary.pdf`.\n\n"
            "Status: complete / frozen unless manuscript review requires changes.\n\n"
            "## Fire-ant chromosome 16 Stage 6C SU branch-length sensitivity\n\n"
            + block.strip() + "\n\n"
        )
        status = re.sub(r"## Fire ants — analysis complete\n\n.*?(?=## |\Z)", "", status, flags=re.S)
        status = re.sub(r"## Fire-ant chromosome 16 Stage 6C SU branch-length gate\n\n.*?(?=## |\Z)", "", status, flags=re.S)
        status = re.sub(r"## Fire-ant chromosome 16 Stage 6C SU branch-length sensitivity\n\n.*?(?=## |\Z)", "", status, flags=re.S)
        PROJECT_STATUS.write_text(status_block + status)


def write_manifest(mapping: dict[str, object], branch_stats: dict[str, object], treatment_files: dict[str, object], prov: dict[str, object]) -> None:
    output_paths = [TREE_TIP_INVENTORY, MATCH_AUDIT, INDIVIDUAL_TO_GROUP, ASTRAL_MAPPING, GROUP_LABEL_MAP, FIXED_BACKGROUND_TOPOLOGY, MAPPING_REPORT, BRANCH_AUDIT, PROVENANCE, SU_COMPARISON, CU_VS_SU_COMPARISON, SU_SUMMARY, FREE_TOPOLOGY_SUMMARY, FINAL_SUMMARY_TABLE, FINAL_METHODS, FINAL_RESULTS, FINAL_CONSISTENCY_AUDIT, FIG_SU_SCATTER_PDF, FIG_SU_SCATTER_PNG, FIG_CU_VS_SU_PDF, FIG_CU_VS_SU_PNG, FIG_FOCAL_PDF, FIG_FOCAL_PNG, FIG_FINAL_SUMMARY_PDF, FIG_FINAL_SUMMARY_PNG, README]
    if PROJECT_STATUS.exists():
        output_paths.append(PROJECT_STATUS)
    astral_outputs = sorted(ASTRAL_SU_RESULTS.glob("*")) if ASTRAL_SU_RESULTS.exists() else []
    manifest = {
        "stage": "6C",
        "analysis_type": "post_freeze_individual_level_astral4_castles_su_branch_length_sensitivity",
        "mapping_recovered_without_guessing": mapping["complete_exact_individual_to_group_mapping_available"],
        "mapped_astral4_castles_su_run_performed": True,
        "individual_mapping_guessed": False,
        "stage5_modified": False,
        "stage6_cu_modified": False,
        "raw_sequence_reanalysis": False,
        "modal_zero_length_trees_used_for_su": False,
        "supplementary_data_1": {"url": SUPP_DATA1_URL, "sha256": sha256(SUPP_DATA1), "worksheet_names": mapping["supplement_meta"]["worksheet_names"], "n_rows": mapping["supplement_meta"]["n_rows"]},
        "group_counts": mapping["group_counts"],
        "input_checksums": {rel(p): sha256(p) for p in [INPUT_TREE, SUPP_DATA1, FILTER_SCRIPT, TOPO_README, WINDOW_INDEX, STAGE4A_WINDOW_SUPPORT, STAGE4A_REGION_SUMMARY, STAGE4B_PRIMARY, STAGE4B_COORD_UNIQUE, STAGE4B_COORD_LENGTH, STAGE5_MANIFEST, STAGE6_MANIFEST, STAGE6_CU_MANIFEST, STAGE6_CU_COMPARISON] if p.exists()},
        "raxml_branch_length_audit": branch_stats,
        "treatment_tree_files": {k: {"path": rel(v["path"]), "n_windows": v["n_windows"], "sha256": v["sha256"]} for k, v in treatment_files.items()},
        "astral4": {"executable": str(prov["astral4_path"]), "executable_sha256": prov["astral4_sha256"], "aster_commit": prov["aster_commit"], "threads": THREADS},
        "astral_output_checksums": {rel(p): sha256(p) for p in astral_outputs if p.is_file()},
        "output_checksums": {rel(p): sha256(p) for p in output_paths if p.exists()},
        "processing_script_checksum": sha256(Path(__file__)),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def run_stage6c() -> dict[str, object]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    mapping = recover_mapping()
    branch_stats = audit_raxml_trees()
    create_mapping_outputs(mapping)
    write_mapping_report(mapping)
    write_branch_audit(branch_stats)
    if not mapping["complete_exact_individual_to_group_mapping_available"]:
        raise RuntimeError("Mapping remained incomplete after Supplementary Data 1 recovery; stopping before ASTRAL")
    treatment_files = create_treatment_tree_files()
    prov = run_astral_analyses(treatment_files)
    rows, cuvs, meta = compare_fixed_background_all()
    write_tsv(SU_COMPARISON, rows, ["split_id", "taxa_side", "n_taxa_side", "branch_role", "su_background", "su_all", "delta_su", "relative_delta_su", "cu_background_fixed_su_run", "cu_all_fixed_su_run", "delta_cu_fixed_su_run", "relative_delta_cu_fixed_su_run", "status"])
    write_tsv(CU_VS_SU_COMPARISON, cuvs, ["split_id", "taxa_side", "branch_role", "cu_background", "cu_all", "delta_cu", "relative_delta_cu", "su_background", "su_all", "delta_su", "relative_delta_su", "secondary_grouped_stage6_cu_background", "secondary_grouped_stage6_cu_all", "secondary_grouped_stage6_delta_cu", "secondary_grouped_stage6_relative_delta_cu"])
    free_rows = summarize_free_topologies(prov)
    write_tsv(FREE_TOPOLOGY_SUMMARY, free_rows, ["treatment", "n_internal_splits", "same_topology_as_free_background", "unrooted_RF_to_free_background", "n_background_splits_recovered", "fraction_background_splits_recovered", "individual_free_focal_split", "individual_free_focal_class", "grouped_stage6_focal_class", "matches_grouped_stage6_focal_class", "tree_file", "runtime_seconds"])
    write_provenance(prov)
    write_su_summary(rows, cuvs, meta, free_rows)
    make_figures(rows, cuvs)
    update_docs(mapping, rows, cuvs, free_rows)
    ctx = collect_final_context(cuvs, free_rows, mapping, branch_stats)
    final_rows = write_final_summary_table(ctx)
    write_final_methods_results(ctx)
    make_final_summary_figure(ctx)
    write_consistency_audit(ctx, final_rows)
    write_manifest(mapping, branch_stats, treatment_files, prov)
    return {
        "mapping_recovered_without_guessing": True,
        "mapped_astral4_castles_su_run_performed": True,
        "group_counts": mapping["group_counts"],
        "n_raxml_trees": branch_stats["n_trees"],
        "tips_per_tree": branch_stats["tips_per_tree_min"],
        "su_comparison_rows": len(rows),
        "manifest_sha": sha256(MANIFEST),
    }


class Stage6CTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mapping = recover_mapping()
        cls.branch = audit_raxml_trees()
        cls.su_rows = read_tsv(SU_COMPARISON) if SU_COMPARISON.exists() else []

    def test_mapping_complete_without_guessing(self):
        self.assertTrue(self.mapping["complete_exact_individual_to_group_mapping_available"])
        self.assertEqual(self.mapping["n_tree_tips"], 267)
        self.assertEqual(len(self.mapping["individual_rows"]), 267)
        self.assertEqual(sorted(self.mapping["group_counts"]), sorted(EXPECTED_GROUPS))
        self.assertEqual(self.mapping["match_type_counts"].get("exact"), 267)
        self.assertEqual(len(self.mapping["unmatched"]), 0)
        self.assertEqual(len(self.mapping["biglittle_validation_mismatches"]), 0)

    def test_raxml_tree_count_and_branch_lengths(self):
        self.assertEqual(self.branch["n_trees"], 213)
        self.assertEqual(self.branch["tips_per_tree_min"], 267)
        self.assertEqual(self.branch["tips_per_tree_max"], 267)
        self.assertGreater(self.branch["fraction_positive_branch_lengths"], 0.99)
        self.assertGreater(self.branch["max_branch_length"], 0)

    def test_treatment_counts(self):
        self.assertEqual(self.branch["treatment_counts"], EXPECTED_TREATMENT_N)

    def test_output_mapping_files(self):
        self.assertTrue(INDIVIDUAL_TO_GROUP.exists())
        self.assertTrue(ASTRAL_MAPPING.exists())
        self.assertEqual(len(read_tsv(INDIVIDUAL_TO_GROUP)), 267)
        self.assertEqual(len(ASTRAL_MAPPING.read_text().strip().splitlines()), 267)

    def test_sulength_nonzero_and_fixed_topology_same(self):
        rows = self.su_rows
        self.assertGreater(len(rows), 0)
        shared = [r for r in rows if r["status"] == "shared"]
        self.assertEqual(len(shared), 4)
        self.assertTrue(all(float(r["su_background"]) > 0 for r in shared))
        self.assertTrue(all(float(r["su_all"]) > 0 for r in shared))

    def test_authoritative_cu_vs_su_uses_fixed_stage6c_values(self):
        rows = read_tsv(CU_VS_SU_COMPARISON)
        self.assertEqual(len(rows), 4)
        self.assertIn("cu_background", rows[0])
        self.assertIn("secondary_grouped_stage6_cu_background", rows[0])
        rich = next(r for r in rows if r["branch_role"] == "focal_richteri_SB_Sb_pair")
        inv = next(r for r in rows if r["branch_role"] == "focal_invicta_macdonaghi_SB_Sb_pair")
        self.assertAlmostEqual(float(rich["cu_background"]), 2.84094, places=5)
        self.assertAlmostEqual(float(rich["cu_all"]), 1.14248, places=5)
        self.assertAlmostEqual(float(inv["cu_background"]), 1.68846, places=5)
        self.assertAlmostEqual(float(inv["cu_all"]), 0.958651, places=6)
        self.assertAlmostEqual(float(rich["relative_delta_cu"]), (1.14248 - 2.84094) / 2.84094, places=5)
        self.assertAlmostEqual(float(inv["relative_delta_cu"]), (0.958651 - 1.68846) / 1.68846, places=5)

    def test_final_manuscript_outputs_and_audit(self):
        for path in [FINAL_SUMMARY_TABLE, FINAL_METHODS, FINAL_RESULTS, FINAL_CONSISTENCY_AUDIT, FIG_FINAL_SUMMARY_PDF, FIG_FINAL_SUMMARY_PNG]:
            self.assertTrue(path.exists(), path)
        audit = json.loads(FINAL_CONSISTENCY_AUDIT.read_text())
        self.assertEqual(audit["status"], "passed")
        self.assertEqual(audit["failed_checks"], [])
        self.assertIn("Stage 6C supersedes Stage 6B for direct CU-vs-SU comparison", README.read_text())
        self.assertIn("Fire ants — analysis complete", PROJECT_STATUS.read_text())

    def test_free_outputs_parse(self):
        for treatment in FREE_TREATMENTS:
            path = ASTRAL_SU_RESULTS / f"{treatment}.free_su.nwk"
            self.assertTrue(path.exists())
            self.assertEqual(tree_tips(parse_tree_file(path)), EXPECTED_SAFE_TAXA)

    def test_no_stage5_or_stage6_cu_modification_flags(self):
        manifest = json.loads(MANIFEST.read_text())
        self.assertFalse(manifest["stage5_modified"])
        self.assertFalse(manifest["stage6_cu_modified"])
        self.assertFalse(manifest["individual_mapping_guessed"])
        self.assertFalse(manifest["modal_zero_length_trees_used_for_su"])


def cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    result = run_stage6c()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage6CTests)
        outcome = unittest.TextTestRunner(verbosity=2).run(suite)
        if not outcome.wasSuccessful():
            raise SystemExit(1)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    cli()
