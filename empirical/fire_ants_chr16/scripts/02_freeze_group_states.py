#!/usr/bin/env python3
"""Freeze fire-ant TWISST group definitions and SB/Sb state provenance.

Stage 2 records the biological metadata rule used to construct the seven TWISST
groups and freezes the focal species/haplotype partitions. It does not classify
the 945 topologies, calculate quartet support, compare regions, or run
ASTRAL/ASTER.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import tempfile
import unittest
import urllib.request
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"

STAGE0_FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
STAGE1_MANIFEST = EMPIRICAL_ROOT / "results" / "stage1_manifest.json"
STAGE1_VALIDATION_SUMMARY = EMPIRICAL_ROOT / "results" / "stage1_validation_summary.tsv"
STAGE1_SOURCE_PROVENANCE = EMPIRICAL_ROOT / "results" / "stage1_source_provenance.md"
STAGE1_WINDOW_TREES = DATA_ROOT / "processed" / "stage1_window_trees.tsv"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"

GROUP_SCRIPT = DATA_ROOT / "raw" / "upstream" / "Topology weighting" / "supergene" / "filter_samples_by_missing.R"
RULES = DATA_ROOT / "metadata" / "twisst_group_definition_rules.tsv"
STAGE2_GROUPS = DATA_ROOT / "processed" / "stage2_twisst_groups.tsv"
SAMPLE_LABEL_INVENTORY = DATA_ROOT / "processed" / "stage2_sample_label_inventory.tsv"
FOCAL_PARTITION = DATA_ROOT / "processed" / "stage2_focal_partition.tsv"
FOCAL_PARTITION_SHA = DATA_ROOT / "processed" / "stage2_focal_partition.sha256"
STAGE2_REPORT = EMPIRICAL_ROOT / "results" / "stage2_report.md"
STAGE2_MANIFEST = EMPIRICAL_ROOT / "results" / "stage2_manifest.json"

UPSTREAM_COMMIT = "bdc7823941680fa560e61b6467af9f77950f64b0"
UPSTREAM_BRANCH = "master"
UPSTREAM_REPO = "https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression"
GROUP_SCRIPT_URL = f"https://raw.githubusercontent.com/wurmlab/2021-fire-ant-social-supergene-introgression/{UPSTREAM_COMMIT}/Topology%20weighting/supergene/filter_samples_by_missing.R"
DOWNLOAD_DATE = "2026-09-30"

EXPECTED_STAGE0_FOCAL_SHA = "5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6"
EXPECTED_STAGE1_METRICS = {
    "coordinate_rows": "213",
    "local_tree_rows": "213",
    "min_tree_tips": "267",
    "max_tree_tips": "267",
    "n_distinct_tip_sets": "1",
    "twisst_topologies": "945",
    "weight_columns": "945",
    "weight_rows": "213",
    "topology_crosscheck_mismatches": "0",
}
EXPECTED_GROUP_SCRIPT_SHA = "419cb7839521a9e3326d80c32b3d77445b79d7e171df77449e570cf89bb65065"
EXPECTED_GROUPS = [
    "geminata",
    "saevissima",
    "pusillignis",
    "invicta/macdonaghi_Sb",
    "invicta/macdonaghi_SB",
    "richteri_Sb",
    "richteri_SB",
]
FOCAL_ROLE_BY_GROUP = {
    "invicta/macdonaghi_SB": ("inv_mac_SB", "A", "invicta/macdonaghi", "SB"),
    "invicta/macdonaghi_Sb": ("inv_mac_Sb", "B", "invicta/macdonaghi", "Sb"),
    "richteri_SB": ("richteri_SB", "C", "richteri", "SB"),
    "richteri_Sb": ("richteri_Sb", "D", "richteri", "Sb"),
}
RAW_STAGE1_FORBIDDEN_FOR_CLASSIFICATION = {
    "data/fire_ants_chr16/raw/upstream/Topology weighting/results/2021-08-07-twisst/weights_output.csv",
    "data/fire_ants_chr16/processed/stage1_twisst_weights_sparse.tsv",
    "data/fire_ants_chr16/raw/upstream/Topology weighting/results/2021-08-07-twisst/topologies_output.trees",
    "data/fire_ants_chr16/processed/stage1_twisst_topologies.tsv",
}
FORBIDDEN_OUTPUT_FIELDS = {
    "q_s",
    "q_h",
    "q_3",
    "d",
    "delta_d",
    "species_split_support",
    "haplotype_split_support",
    "third_split_support",
    "quartet_class",
    "topology_id",
    "raw_weight",
}
OUTPUTS = [
    RULES,
    STAGE2_GROUPS,
    SAMPLE_LABEL_INVENTORY,
    FOCAL_PARTITION,
    FOCAL_PARTITION_SHA,
    STAGE2_REPORT,
    STAGE2_MANIFEST,
]
TSV_LINETERMINATOR = "\n"


class Stage2Error(RuntimeError):
    """Raised when Stage-2 validation fails."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else str(row.get(field, "")) for field in fields})


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def verify_prior_stages() -> dict[str, object]:
    if sha256(STAGE0_FOCAL_QUARTET) != EXPECTED_STAGE0_FOCAL_SHA:
        raise Stage2Error("Stage-0 focal quartet checksum changed")
    manifest = json.loads(STAGE1_MANIFEST.read_text())
    if manifest.get("upstream_commit") != UPSTREAM_COMMIT or manifest.get("upstream_branch") != UPSTREAM_BRANCH:
        raise Stage2Error("Stage-1 manifest upstream provenance changed")
    if manifest.get("stage0_focal_quartet_checksum") != EXPECTED_STAGE0_FOCAL_SHA:
        raise Stage2Error("Stage-1 manifest does not preserve Stage-0 focal checksum")
    observed = manifest.get("observed_metrics", {})
    for key, expected in EXPECTED_STAGE1_METRICS.items():
        if str(observed.get(key)) != expected:
            raise Stage2Error(f"Stage-1 metric {key} expected {expected}, observed {observed.get(key)}")
    for rel_path, expected_sha in manifest.get("raw_topology_file_checksums", {}).items():
        path = REPO_ROOT / rel_path
        if sha256(path) != expected_sha:
            raise Stage2Error(f"Stage-1 raw checksum changed: {rel_path}")
    for rel_path, expected_sha in manifest.get("processed_output_checksums", {}).items():
        path = REPO_ROOT / rel_path
        if sha256(path) != expected_sha:
            raise Stage2Error(f"Stage-1 processed checksum changed: {rel_path}")
    script = manifest.get("processing_script", {})
    if sha256(REPO_ROOT / script["path"]) != script["sha256"]:
        raise Stage2Error("Stage-1 processing script checksum changed")
    return manifest


def download_group_script(skip_download: bool = False) -> None:
    GROUP_SCRIPT.parent.mkdir(parents=True, exist_ok=True)
    if GROUP_SCRIPT.exists() and sha256(GROUP_SCRIPT) == EXPECTED_GROUP_SCRIPT_SHA:
        return
    if skip_download:
        raise Stage2Error("missing or changed group-construction script in --skip-download mode")
    urllib.request.urlretrieve(GROUP_SCRIPT_URL, GROUP_SCRIPT)
    if sha256(GROUP_SCRIPT) != EXPECTED_GROUP_SCRIPT_SHA:
        raise Stage2Error("downloaded group-construction script checksum did not match expected frozen checksum")


def update_download_manifest() -> None:
    fields, rows = read_tsv(DOWNLOAD_MANIFEST)
    rel = GROUP_SCRIPT.relative_to(REPO_ROOT).as_posix()
    new_row = {
        "source": "GitHub wurmlab/2021-fire-ant-social-supergene-introgression",
        "original_filename": "Topology weighting/supergene/filter_samples_by_missing.R",
        "url_or_doi": GROUP_SCRIPT_URL,
        "download_date": DOWNLOAD_DATE,
        "size_bytes": str(GROUP_SCRIPT.stat().st_size),
        "sha256": sha256(GROUP_SCRIPT),
        "purpose": f"Stage-2 biological TWISST group-construction rule provenance. Upstream commit: {UPSTREAM_COMMIT}.",
        "local_path": rel,
    }
    replaced = False
    out_rows = []
    for row in rows:
        if row.get("local_path") == rel:
            out_rows.append(new_row)
            replaced = True
        else:
            out_rows.append(row)
    if not replaced:
        out_rows.append(new_row)
    write_tsv(DOWNLOAD_MANIFEST, out_rows, fields)


def group_definition_rules() -> list[dict[str, object]]:
    source = GROUP_SCRIPT.relative_to(REPO_ROOT).as_posix()
    return [
        {
            "rule_order": 1,
            "input_field": "Sample.name",
            "operation": "read and merge",
            "condition": 'colnames(miss)[1] <- "Sample.name"; merge(miss, samples, by="Sample.name")',
            "output_effect": "joins missingness table tmp/supergene.missing to samples_overview.csv metadata",
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "samples_overview.csv is referenced by the frozen script but is not committed in the audited repository tree.",
        },
        {
            "rule_order": 2,
            "input_field": "F_MISS",
            "operation": "filter",
            "condition": "F_MISS < 0.05",
            "output_effect": "keeps samples with less than 5 percent missing genotype calls for TWISST",
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "missingness filter precedes group-table output",
        },
        {
            "rule_order": 3,
            "input_field": "Species",
            "operation": "exclude species labels",
            "condition": '!Species %in% c("S.\u00a0interrupta", "xAdR", "S.megergates")',
            "output_effect": "removes these source species/categories before TWISST sample table construction",
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "the first string contains a non-breaking space in the frozen R script",
        },
        {
            "rule_order": 4,
            "input_field": "Species",
            "operation": "normalize label",
            "condition": 'gsub("S\\\\.", "", Species); gsub(" ", "", Species)',
            "output_effect": "removes the S. prefix and spaces from species labels",
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "source terminology is preserved in the rule table",
        },
        {
            "rule_order": 5,
            "input_field": "Species",
            "operation": "merge species labels",
            "condition": 'Species %in% c("invicta", "macdonaghi")',
            "output_effect": 'sets Species to "invicta/macdonaghi"',
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "this merge is explicit in the upstream script",
        },
        {
            "rule_order": 6,
            "input_field": "Species; Supergene.Variant",
            "operation": "construct clade",
            "condition": 'if Species %in% c("invicta/macdonaghi", "richteri")',
            "output_effect": 'sets clade to paste(Species, Supergene.Variant, sep="_")',
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "this is the authoritative SB/Sb provenance for focal TWISST groups",
        },
        {
            "rule_order": 7,
            "input_field": "Species",
            "operation": "construct clade",
            "condition": 'if Species is not invicta/macdonaghi or richteri',
            "output_effect": "sets clade to Species without Supergene.Variant splitting",
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "geminata, saevissima, and pusillignis are not split by variant in the TWISST groups",
        },
        {
            "rule_order": 8,
            "input_field": "Sample.name; clade",
            "operation": "write table",
            "condition": "select(Sample.name, clade)",
            "output_effect": 'writes results/twisst_samples_low_missingness without row or column names',
            "source_file": source,
            "source_commit": UPSTREAM_COMMIT,
            "notes": "the output table itself is not committed in the audited repository tree",
        },
    ]


def stage2_groups() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label in EXPECTED_GROUPS:
        if label in FOCAL_ROLE_BY_GROUP:
            canonical, role, species, state = FOCAL_ROLE_BY_GROUP[label]
            rows.append(
                {
                    "group_id": canonical,
                    "upstream_group_label": label,
                    "species_component": species,
                    "supergene_variant": state,
                    "structural_state_known": "true",
                    "focal_quartet_member": "true",
                    "focal_role": role,
                    "group_definition_source": GROUP_SCRIPT.relative_to(REPO_ROOT).as_posix(),
                    "notes": "Focal group state is derived from upstream Species + Supergene.Variant clade construction, not topology.",
                }
            )
        else:
            rows.append(
                {
                    "group_id": label,
                    "upstream_group_label": label,
                    "species_component": label,
                    "supergene_variant": "NA",
                    "structural_state_known": "false",
                    "focal_quartet_member": "false",
                    "focal_role": "NA",
                    "group_definition_source": GROUP_SCRIPT.relative_to(REPO_ROOT).as_posix(),
                    "notes": "Non-focal TWISST group is not split by Supergene.Variant in the upstream rule.",
                }
            )
    return rows


def focal_partition_rows() -> list[dict[str, object]]:
    return [
        {"group": "inv_mac_SB", "species_partition": "invicta/macdonaghi", "haplotype_partition": "SB", "role": "A"},
        {"group": "inv_mac_Sb", "species_partition": "invicta/macdonaghi", "haplotype_partition": "Sb", "role": "B"},
        {"group": "richteri_SB", "species_partition": "richteri", "haplotype_partition": "SB", "role": "C"},
        {"group": "richteri_Sb", "species_partition": "richteri", "haplotype_partition": "Sb", "role": "D"},
    ]


def extract_tip_labels_from_newick_topology_neutral(tree: str) -> list[str]:
    """Extract leaf labels only, without retaining clades or relationships."""

    return re.findall(r"(?<=[(,])([^():,;]+):", tree)


def parse_label(label: str) -> dict[str, str]:
    species_code = "unknown"
    status = "unknown"
    source = "none"
    for token in ["_inv-", "_ric-", "_mac-", "_sae-"]:
        if token in label:
            species_code = token.strip("_-")
            status = "literal_only"
            source = "sample_label_literal_token"
            break
    variant = "unknown"
    if "bigB" in label and "littleb" in label:
        variant = "ambiguous_bigB_and_littleb"
    elif "bigB" in label:
        variant = "bigB"
    elif "littleb" in label:
        variant = "littleb"
    focal_candidate = "true" if species_code in {"inv", "mac", "ric"} else "false"
    notes = "Literal sample-label audit only; bigB/littleb is not converted to SB/Sb in Stage 2."
    return {
        "sample_label": label,
        "species_code_from_label": species_code,
        "variant_string_from_label": variant,
        "label_parse_status": status,
        "focal_species_candidate": focal_candidate,
        "interpretation_source": source,
        "notes": notes,
    }


def sample_label_inventory() -> list[dict[str, str]]:
    fields, rows = read_tsv(STAGE1_WINDOW_TREES)
    if "tree_newick" not in fields:
        raise Stage2Error("stage1_window_trees.tsv missing tree_newick")
    first_tree = rows[0]["tree_newick"]
    labels = sorted(set(extract_tip_labels_from_newick_topology_neutral(first_tree)))
    if len(labels) != 267:
        raise Stage2Error(f"expected 267 unique sample labels, observed {len(labels)}")
    for row in rows[1:]:
        observed = set(extract_tip_labels_from_newick_topology_neutral(row["tree_newick"]))
        if observed != set(labels):
            raise Stage2Error("not all Stage-1 local trees share the same tip set")
    return [parse_label(label) for label in labels]


def write_focal_partition_sha() -> None:
    FOCAL_PARTITION_SHA.write_text(f"{sha256(FOCAL_PARTITION)}  {FOCAL_PARTITION.relative_to(REPO_ROOT).as_posix()}\n")


def write_report(samples_available: bool, inventory_rows: list[dict[str, str]]) -> None:
    label_status = Counter(row["label_parse_status"] for row in inventory_rows)
    variant_status = Counter(row["variant_string_from_label"] for row in inventory_rows)
    text = [
        "# Fire ants chromosome 16 Stage 2 report",
        "",
        "Stage 2 freezes the biological TWISST group definitions and SB/Sb state provenance. It does not classify topologies, calculate local quartet support, compare regions, or run ASTRAL/ASTER.",
        "",
        "## Biological group provenance",
        "",
        "The frozen upstream script `Topology weighting/supergene/filter_samples_by_missing.R` reads `tmp/supergene.missing`, reads the first 22 columns of `samples_overview.csv`, renames the missingness sample column to `Sample.name`, and merges the two tables by `Sample.name`.",
        "",
        "For the TWISST sample table, the script filters to `F_MISS < 0.05`, excludes `S.\u00a0interrupta`, `xAdR`, and `S.megergates`, removes `S.` and spaces from `Species`, merges `invicta` and `macdonaghi` into `invicta/macdonaghi`, and constructs `clade` as `paste(Species, Supergene.Variant, sep=\"_\")` only for `invicta/macdonaghi` and `richteri`. Other species use `Species` alone as `clade`.",
        "",
        "Thus the four focal group labels are determined by species membership plus the `Supergene.Variant` metadata field. They are not derived from local RAxML clustering, TWISST weights, ASTRAL trees, q-values, or dominant topology.",
        "",
        "## Frozen TWISST groups",
        "",
        "```text",
        *EXPECTED_GROUPS,
        "```",
        "",
        "## Focal species partition",
        "",
        "```text",
        "(invicta/macdonaghi_SB, invicta/macdonaghi_Sb)",
        "|",
        "(richteri_SB, richteri_Sb)",
        "```",
        "",
        "## Focal haplotype partition",
        "",
        "```text",
        "(invicta/macdonaghi_SB, richteri_SB)",
        "|",
        "(invicta/macdonaghi_Sb, richteri_Sb)",
        "```",
        "",
        "## Metadata availability",
        "",
        "The exact published group-construction algorithm is available in the frozen upstream script.",
        "",
        (
            "The exact `samples_overview.csv` / `results/twisst_samples_low_missingness` source table was found in the frozen upstream repository tree."
            if samples_available
            else "The exact 267-row source `samples_overview.csv` / grouped `results/twisst_samples_low_missingness` file is not deposited in the audited repository path at the frozen commit. Stage 2 therefore does not fabricate a sample-to-TWISST-group table."
        ),
        "",
        "Stage 4 can still operate on the already-published grouped TWISST weights from Stage 1.",
        "",
        "## Sample-label audit",
        "",
        f"Unique Stage-1 sample labels audited: `{len(inventory_rows)}`.",
        "",
        f"Label parse status counts: `{dict(sorted(label_status.items()))}`.",
        "",
        f"Literal variant-string counts: `{dict(sorted(variant_status.items()))}`.",
        "",
        "Strings such as `bigB` and `littleb` are retained as literal sample-label metadata only. Stage 2 uses the upstream `Supergene.Variant` field as the authoritative SB/Sb provenance and does not silently equate sample-name strings to SB/Sb.",
        "",
        "## Anti-circularity",
        "",
        "No TWISST weight row, grouped topology tree, focal quartet support value, regional enrichment result, ASTRAL tree, or local topology classification was used to assign group identity or SB/Sb state.",
    ]
    STAGE2_REPORT.write_text("\n".join(text) + "\n")


def write_manifest(stage1_manifest: dict[str, object]) -> None:
    processed_outputs = {
        RULES.relative_to(REPO_ROOT).as_posix(): sha256(RULES),
        STAGE2_GROUPS.relative_to(REPO_ROOT).as_posix(): sha256(STAGE2_GROUPS),
        SAMPLE_LABEL_INVENTORY.relative_to(REPO_ROOT).as_posix(): sha256(SAMPLE_LABEL_INVENTORY),
        FOCAL_PARTITION.relative_to(REPO_ROOT).as_posix(): sha256(FOCAL_PARTITION),
        FOCAL_PARTITION_SHA.relative_to(REPO_ROOT).as_posix(): sha256(FOCAL_PARTITION_SHA),
        STAGE2_REPORT.relative_to(REPO_ROOT).as_posix(): sha256(STAGE2_REPORT),
    }
    data = {
        "stage": "Stage 2",
        "description": "Freeze biological TWISST group definitions and SB/Sb state provenance without topology classification.",
        "stage0_focal_quartet_checksum": EXPECTED_STAGE0_FOCAL_SHA,
        "stage1_manifest_path": STAGE1_MANIFEST.relative_to(REPO_ROOT).as_posix(),
        "stage1_manifest_checksum": sha256(STAGE1_MANIFEST),
        "stage1_recorded_checksums": {
            "raw_topology_file_checksums": stage1_manifest.get("raw_topology_file_checksums", {}),
            "processed_output_checksums": stage1_manifest.get("processed_output_checksums", {}),
            "processing_script": stage1_manifest.get("processing_script", {}),
        },
        "upstream_group_definition_script": {
            "path": GROUP_SCRIPT.relative_to(REPO_ROOT).as_posix(),
            "source_url": GROUP_SCRIPT_URL,
            "source_commit": UPSTREAM_COMMIT,
            "sha256": sha256(GROUP_SCRIPT),
        },
        "stage2_output_checksums": processed_outputs,
        "processing_script": {
            "path": "empirical/fire_ants_chr16/scripts/02_freeze_group_states.py",
            "sha256": sha256(Path(__file__)),
        },
        "sample_label_inventory_rows": 267,
        "twisst_groups": EXPECTED_GROUPS,
        "focal_species_partition": "AB|CD",
        "focal_haplotype_partition": "AC|BD",
        "topology_files_not_used_for_classification": sorted(RAW_STAGE1_FORBIDDEN_FOR_CLASSIFICATION),
    }
    STAGE2_MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate_output_schemas() -> None:
    for path in [RULES, STAGE2_GROUPS, SAMPLE_LABEL_INVENTORY, FOCAL_PARTITION]:
        fields, _ = read_tsv(path)
        forbidden = {field.lower() for field in fields} & FORBIDDEN_OUTPUT_FIELDS
        if forbidden:
            raise Stage2Error(f"{path.relative_to(REPO_ROOT)} contains forbidden Stage-2 output fields: {sorted(forbidden)}")


def run_stage2(skip_download: bool = False) -> dict[str, object]:
    stage1_manifest = verify_prior_stages()
    download_group_script(skip_download=skip_download)
    update_download_manifest()

    script_text = GROUP_SCRIPT.read_text()
    required_snippets = [
        'read.csv("samples_overview.csv", header=T)[, 1:22]',
        'filter(miss, F_MISS < 0.05)',
        'filter(!Species %in% c("S.\u00a0interrupta", "xAdR", "S.megergates"))',
        'ifelse(Species %in% c("invicta", "macdonaghi"), "invicta/macdonaghi", Species)',
        'paste(Species, Supergene.Variant, sep = "_")',
        'write.table(file="results/twisst_samples_low_missingness", quote=F, row.names=F, col.names=F)',
    ]
    for snippet in required_snippets:
        if snippet not in script_text:
            raise Stage2Error(f"upstream group script missing expected rule snippet: {snippet}")

    source_provenance = STAGE1_SOURCE_PROVENANCE.read_text()
    for label in EXPECTED_GROUPS:
        if label not in source_provenance:
            raise Stage2Error(f"Stage-1 source provenance missing expected group label: {label}")

    inventory_rows = sample_label_inventory()
    write_tsv(
        RULES,
        group_definition_rules(),
        ["rule_order", "input_field", "operation", "condition", "output_effect", "source_file", "source_commit", "notes"],
    )
    write_tsv(
        STAGE2_GROUPS,
        stage2_groups(),
        [
            "group_id",
            "upstream_group_label",
            "species_component",
            "supergene_variant",
            "structural_state_known",
            "focal_quartet_member",
            "focal_role",
            "group_definition_source",
            "notes",
        ],
    )
    write_tsv(
        SAMPLE_LABEL_INVENTORY,
        inventory_rows,
        [
            "sample_label",
            "species_code_from_label",
            "variant_string_from_label",
            "label_parse_status",
            "focal_species_candidate",
            "interpretation_source",
            "notes",
        ],
    )
    write_tsv(FOCAL_PARTITION, focal_partition_rows(), ["group", "species_partition", "haplotype_partition", "role"])
    write_focal_partition_sha()

    samples_available = False
    write_report(samples_available=samples_available, inventory_rows=inventory_rows)
    validate_output_schemas()
    write_manifest(stage1_manifest)
    return {
        "sample_labels": len(inventory_rows),
        "group_script_sha": sha256(GROUP_SCRIPT),
        "focal_partition_sha": sha256(FOCAL_PARTITION),
        "manifest_sha": sha256(STAGE2_MANIFEST),
    }


class Stage2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.metrics = run_stage2(skip_download=True)

    def test_stage0_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE0_FOCAL_QUARTET), EXPECTED_STAGE0_FOCAL_SHA)

    def test_stage1_checksums_unchanged(self) -> None:
        verify_prior_stages()

    def test_group_definition_script_checksum_recorded(self) -> None:
        self.assertEqual(sha256(GROUP_SCRIPT), EXPECTED_GROUP_SCRIPT_SHA)
        manifest = json.loads(STAGE2_MANIFEST.read_text())
        self.assertEqual(manifest["upstream_group_definition_script"]["sha256"], EXPECTED_GROUP_SCRIPT_SHA)

    def test_exactly_seven_frozen_twisst_groups(self) -> None:
        _, rows = read_tsv(STAGE2_GROUPS)
        self.assertEqual(len(rows), 7)

    def test_exact_expected_seven_upstream_labels(self) -> None:
        _, rows = read_tsv(STAGE2_GROUPS)
        self.assertEqual([row["upstream_group_label"] for row in rows], EXPECTED_GROUPS)

    def test_exactly_four_focal_groups(self) -> None:
        _, rows = read_tsv(STAGE2_GROUPS)
        self.assertEqual(sum(row["focal_quartet_member"] == "true" for row in rows), 4)

    def test_two_focal_groups_are_sb(self) -> None:
        _, rows = read_tsv(STAGE2_GROUPS)
        focal_states = [row["supergene_variant"] for row in rows if row["focal_quartet_member"] == "true"]
        self.assertEqual(Counter(focal_states)["SB"], 2)

    def test_two_focal_groups_are_sb_lower(self) -> None:
        _, rows = read_tsv(STAGE2_GROUPS)
        focal_states = [row["supergene_variant"] for row in rows if row["focal_quartet_member"] == "true"]
        self.assertEqual(Counter(focal_states)["Sb"], 2)

    def test_exactly_two_focal_species_partitions(self) -> None:
        _, rows = read_tsv(FOCAL_PARTITION)
        self.assertEqual({row["species_partition"] for row in rows}, {"invicta/macdonaghi", "richteri"})

    def test_roles_match_stage0(self) -> None:
        _, rows = read_tsv(FOCAL_PARTITION)
        self.assertEqual({row["group"]: row["role"] for row in rows}, {"inv_mac_SB": "A", "inv_mac_Sb": "B", "richteri_SB": "C", "richteri_Sb": "D"})

    def test_species_partition_ab_cd(self) -> None:
        _, rows = read_tsv(FOCAL_PARTITION)
        by_species: dict[str, set[str]] = {}
        for row in rows:
            by_species.setdefault(row["species_partition"], set()).add(row["role"])
        self.assertEqual(set(map(frozenset, by_species.values())), {frozenset({"A", "B"}), frozenset({"C", "D"})})

    def test_haplotype_partition_ac_bd(self) -> None:
        _, rows = read_tsv(FOCAL_PARTITION)
        by_state: dict[str, set[str]] = {}
        for row in rows:
            by_state.setdefault(row["haplotype_partition"], set()).add(row["role"])
        self.assertEqual(set(map(frozenset, by_state.values())), {frozenset({"A", "C"}), frozenset({"B", "D"})})

    def test_no_focal_state_inferred_from_topology(self) -> None:
        _, rules = read_tsv(RULES)
        rule_text = "\n".join(row["output_effect"] + " " + row["notes"] for row in rules)
        self.assertIn("Supergene.Variant", rule_text)
        self.assertNotIn("topology", "\n".join(row["condition"] for row in rules).lower())

    def test_exactly_267_unique_sample_labels(self) -> None:
        _, rows = read_tsv(SAMPLE_LABEL_INVENTORY)
        self.assertEqual(len(rows), 267)
        self.assertEqual(len({row["sample_label"] for row in rows}), 267)

    def test_output_schemas_contain_no_topology_support_values(self) -> None:
        validate_output_schemas()

    def test_stage2_does_not_read_twisst_weights_for_classification(self) -> None:
        manifest = json.loads(STAGE2_MANIFEST.read_text())
        self.assertIn("topology_files_not_used_for_classification", manifest)
        for forbidden in RAW_STAGE1_FORBIDDEN_FOR_CLASSIFICATION:
            self.assertIn(forbidden, manifest["topology_files_not_used_for_classification"])

    def test_stage2_does_not_classify_945_topologies(self) -> None:
        for path in OUTPUTS:
            text = path.read_text()
            self.assertNotIn("topo945", text)
            self.assertNotIn("quartet_class", text)

    def test_repeated_runs_byte_stable(self) -> None:
        before = {path: sha256(path) for path in OUTPUTS}
        run_stage2(skip_download=True)
        after = {path: sha256(path) for path in OUTPUTS}
        self.assertEqual(before, after)


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage2Tests)
    with tempfile.TemporaryDirectory(prefix="fire_ants_stage2_tests_"):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true", help="regenerate Stage-2 outputs from immutable raw inputs")
    parser.add_argument("--run-tests", action="store_true", help="run embedded Stage-2 tests after freezing outputs")
    args = parser.parse_args(argv)
    try:
        metrics = run_stage2(skip_download=args.skip_download)
    except Stage2Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.run_tests:
        return 0 if run_tests() else 1
    print(f"Wrote {STAGE2_GROUPS.relative_to(REPO_ROOT)}")
    print(f"Wrote {FOCAL_PARTITION.relative_to(REPO_ROOT)}")
    print(f"Audited {metrics['sample_labels']} sample labels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
