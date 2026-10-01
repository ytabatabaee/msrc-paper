#!/usr/bin/env python3
"""Freeze the independent fire-ant chromosome 1-15 background history.

Stage 3 records the published chromosome 1-15 ASTRAL species-history baseline
before the formal local TWISST support test. It does not inspect the supergene
ASTRAL tree, read TWISST weights, classify TWISST topologies, calculate
quartet-support statistics, or compare genomic regions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
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
STAGE2_MANIFEST = EMPIRICAL_ROOT / "results" / "stage2_manifest.json"
STAGE2_FOCAL_PARTITION = DATA_ROOT / "processed" / "stage2_focal_partition.tsv"
STAGE0_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"
COALESCENCE_README = DATA_ROOT / "raw" / "upstream" / "Coalescence-based phylogenies" / "README.md"

RAW_BACKGROUND_TREE = DATA_ROOT / "raw" / "upstream" / "Trees in newick format" / "Astral.10SNP-genes.chr1-15.nwk"
PROCESSED_BACKGROUND_TREE = DATA_ROOT / "processed" / "stage3_background_tree.nwk"
PROCESSED_BACKGROUND_TREE_SHA = DATA_ROOT / "processed" / "stage3_background_tree.sha256"
BACKGROUND_PROVENANCE = DATA_ROOT / "metadata" / "background_tree_provenance.tsv"
BACKGROUND_PARTITION = DATA_ROOT / "processed" / "stage3_background_partition.tsv"
BACKGROUND_PARTITION_SHA = DATA_ROOT / "processed" / "stage3_background_partition.sha256"
STAGE3_REPORT = EMPIRICAL_ROOT / "results" / "stage3_report.md"
STAGE3_MANIFEST = EMPIRICAL_ROOT / "results" / "stage3_manifest.json"

UPSTREAM_REPO = "https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression"
UPSTREAM_COMMIT = "bdc7823941680fa560e61b6467af9f77950f64b0"
UPSTREAM_BRANCH = "master"
BACKGROUND_TREE_URL = f"https://raw.githubusercontent.com/wurmlab/2021-fire-ant-social-supergene-introgression/{UPSTREAM_COMMIT}/Trees%20in%20newick%20format/Astral.10SNP-genes.chr1-15.nwk"
DOWNLOAD_DATE = "2026-10-01"

EXPECTED_STAGE0_FOCAL_SHA = "5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6"
EXPECTED_STAGE1_MANIFEST_SHA = "fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b"
EXPECTED_STAGE2_MANIFEST_SHA = "2b56ac7ea4a4d0a31588b803b30896522af4e90bc1e91bf23db6d0b9db80e427"
EXPECTED_STAGE2_FOCAL_PARTITION_SHA = "52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002"
EXPECTED_BACKGROUND_TREE_SHA = "64c032bf646c1bdef4699d419bd19686f0d9713a57755f0c77b30aea27668638"
EXPECTED_STAGE1_METRICS = {
    "local_tree_rows": "213",
    "min_tree_tips": "267",
    "max_tree_tips": "267",
    "n_distinct_tip_sets": "1",
    "twisst_topologies": "945",
    "weight_rows": "213",
    "weight_columns": "945",
    "topology_crosscheck_mismatches": "0",
}
EXPECTED_BACKGROUND_RELATIONSHIP = "invicta/macdonaghi|richteri"
EXPECTED_STAGE3_SPLIT = "AB|CD"
FORBIDDEN_INPUT_PATTERNS = (
    re.compile(r"Astral[.]10SNP-genes[.]supergene[.]nwk$", re.IGNORECASE),
    re.compile(r"weights_output[.]csv$", re.IGNORECASE),
    re.compile(r"stage1_twisst_weights_sparse[.]tsv$", re.IGNORECASE),
    re.compile(r"topologies_output[.]trees$", re.IGNORECASE),
    re.compile(r"stage1_twisst_topologies[.]tsv$", re.IGNORECASE),
)
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
    PROCESSED_BACKGROUND_TREE,
    PROCESSED_BACKGROUND_TREE_SHA,
    BACKGROUND_PROVENANCE,
    BACKGROUND_PARTITION,
    BACKGROUND_PARTITION_SHA,
    STAGE3_REPORT,
    STAGE3_MANIFEST,
]
TSV_LINETERMINATOR = "\n"


class Stage3Error(RuntimeError):
    """Raised when Stage-3 validation fails."""


class ParsedTree:
    def __init__(self, tips: tuple[str, ...], splits: tuple[str, ...]) -> None:
        self.tips = tips
        self.splits = splits


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def assert_allowed_input(path: Path) -> None:
    normalized = path.as_posix()
    for pattern in FORBIDDEN_INPUT_PATTERNS:
        if pattern.search(normalized):
            raise Stage3Error(f"Stage 3 must not read forbidden input: {path}")


def read_text_checked(path: Path) -> str:
    assert_allowed_input(path)
    return path.read_text()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    assert_allowed_input(path)
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


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.i = 0

    def parse(self) -> object:
        if not self.text:
            raise Stage3Error("empty Newick tree")
        node = self.parse_subtree()
        self.skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ";":
            raise Stage3Error("Newick tree missing terminal semicolon")
        self.i += 1
        self.skip_ws()
        if self.i != len(self.text):
            raise Stage3Error("unexpected content after Newick semicolon")
        return node

    def skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def parse_subtree(self) -> object:
        self.skip_ws()
        if self.i >= len(self.text):
            raise Stage3Error("unexpected end of Newick tree")
        if self.text[self.i] == "(":
            self.i += 1
            children = [self.parse_subtree()]
            while True:
                self.skip_ws()
                if self.i >= len(self.text):
                    raise Stage3Error("unterminated internal node")
                if self.text[self.i] == ",":
                    self.i += 1
                    children.append(self.parse_subtree())
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise Stage3Error(f"unexpected Newick character {self.text[self.i]!r}")
            self.consume_optional_label()
            self.consume_optional_branch_length()
            if len(children) < 2:
                raise Stage3Error("internal node has fewer than two children")
            return tuple(children)
        label = self.consume_tip_label()
        if not label:
            raise Stage3Error("empty tip label")
        self.consume_optional_branch_length()
        return label

    def consume_tip_label(self) -> str:
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in "(),:;":
            self.i += 1
        return self.text[start:self.i].strip()

    def consume_optional_label(self) -> None:
        self.skip_ws()
        while self.i < len(self.text) and self.text[self.i] not in ",:();":
            self.i += 1

    def consume_optional_branch_length(self) -> None:
        self.skip_ws()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            start = self.i
            while self.i < len(self.text) and self.text[self.i] not in ",();":
                self.i += 1
            if self.i == start:
                raise Stage3Error("empty branch length")


def collect_tips(node: object, tips: list[str]) -> None:
    if isinstance(node, str):
        tips.append(node)
    else:
        for child in node:
            collect_tips(child, tips)


def graph_and_leaves(node: object) -> tuple[dict[int, set[int]], dict[int, str]]:
    graph: dict[int, set[int]] = {}
    leaves: dict[int, str] = {}
    next_id = 0

    def add(current: object, parent: int | None = None) -> int:
        nonlocal next_id
        node_id = next_id
        next_id += 1
        graph.setdefault(node_id, set())
        if parent is not None:
            graph[node_id].add(parent)
            graph[parent].add(node_id)
        if isinstance(current, str):
            leaves[node_id] = current
        else:
            for child in current:
                add(child, node_id)
        return node_id

    add(node)
    return graph, leaves


def unrooted_splits(node: object) -> tuple[str, ...]:
    graph, leaves_by_id = graph_and_leaves(node)
    all_labels = frozenset(leaves_by_id.values())
    leaf_ids = set(leaves_by_id)
    splits: set[str] = set()
    for a, neighbors in graph.items():
        for b in neighbors:
            if a > b:
                continue
            stack = [a]
            seen = {b}
            component: set[int] = set()
            while stack:
                current = stack.pop()
                if current in seen:
                    continue
                seen.add(current)
                component.add(current)
                stack.extend(graph[current] - seen)
            left = frozenset(leaves_by_id[i] for i in component & leaf_ids)
            if len(left) <= 1 or len(left) >= len(all_labels) - 1:
                continue
            right = all_labels - left
            splits.add("|".join(sorted([",".join(sorted(left)), ",".join(sorted(right))])))
    return tuple(sorted(splits))


def parse_newick(text: str) -> ParsedTree:
    root = NewickParser(text).parse()
    tips: list[str] = []
    collect_tips(root, tips)
    if len(tips) != len(set(tips)):
        duplicates = sorted(label for label, count in Counter(tips).items() if count > 1)
        raise Stage3Error(f"duplicate tip labels: {', '.join(duplicates)}")
    return ParsedTree(tips=tuple(tips), splits=unrooted_splits(root))


def verify_prior_stages() -> tuple[dict[str, object], dict[str, object]]:
    if sha256(STAGE0_FOCAL_QUARTET) != EXPECTED_STAGE0_FOCAL_SHA:
        raise Stage3Error("Stage-0 focal quartet checksum changed")
    if sha256(STAGE1_MANIFEST) != EXPECTED_STAGE1_MANIFEST_SHA:
        raise Stage3Error("Stage-1 manifest checksum changed")
    if sha256(STAGE2_MANIFEST) != EXPECTED_STAGE2_MANIFEST_SHA:
        raise Stage3Error("Stage-2 manifest checksum changed")
    if sha256(STAGE2_FOCAL_PARTITION) != EXPECTED_STAGE2_FOCAL_PARTITION_SHA:
        raise Stage3Error("Stage-2 focal partition checksum changed")

    stage1 = json.loads(read_text_checked(STAGE1_MANIFEST))
    stage2 = json.loads(read_text_checked(STAGE2_MANIFEST))
    if stage1.get("upstream_commit") != UPSTREAM_COMMIT:
        raise Stage3Error("Stage-1 upstream commit changed")
    for key, expected in EXPECTED_STAGE1_METRICS.items():
        if str(stage1.get("observed_metrics", {}).get(key)) != expected:
            raise Stage3Error(f"Stage-1 metric {key} expected {expected}, observed {stage1.get('observed_metrics', {}).get(key)}")
    for rel_path, expected_sha in stage1.get("raw_topology_file_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage3Error(f"Stage-1 raw checksum changed: {rel_path}")
    for rel_path, expected_sha in stage1.get("processed_output_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage3Error(f"Stage-1 processed checksum changed: {rel_path}")
    if stage2.get("stage1_manifest_checksum") != EXPECTED_STAGE1_MANIFEST_SHA:
        raise Stage3Error("Stage-2 manifest no longer points at the frozen Stage-1 manifest")
    for rel_path, expected_sha in stage2.get("stage2_output_checksums", {}).items():
        if sha256(REPO_ROOT / rel_path) != expected_sha:
            raise Stage3Error(f"Stage-2 output checksum changed: {rel_path}")
    return stage1, stage2


def download_background_tree(skip_download: bool = False) -> None:
    RAW_BACKGROUND_TREE.parent.mkdir(parents=True, exist_ok=True)
    if RAW_BACKGROUND_TREE.exists() and sha256(RAW_BACKGROUND_TREE) == EXPECTED_BACKGROUND_TREE_SHA:
        return
    if skip_download:
        raise Stage3Error("missing or changed chr1-15 background tree in --skip-download mode")
    urllib.request.urlretrieve(BACKGROUND_TREE_URL, RAW_BACKGROUND_TREE)
    if sha256(RAW_BACKGROUND_TREE) != EXPECTED_BACKGROUND_TREE_SHA:
        raise Stage3Error("downloaded chr1-15 background tree checksum did not match expected frozen checksum")


def update_download_manifest() -> None:
    fields, rows = read_tsv(DOWNLOAD_MANIFEST)
    rel = RAW_BACKGROUND_TREE.relative_to(REPO_ROOT).as_posix()
    new_row = {
        "source": "GitHub wurmlab/2021-fire-ant-social-supergene-introgression",
        "original_filename": "Trees in newick format/Astral.10SNP-genes.chr1-15.nwk",
        "url_or_doi": BACKGROUND_TREE_URL,
        "download_date": DOWNLOAD_DATE,
        "size_bytes": str(RAW_BACKGROUND_TREE.stat().st_size),
        "sha256": sha256(RAW_BACKGROUND_TREE),
        "purpose": f"Published chromosome 1-15 ASTRAL species-history baseline. Upstream commit: {UPSTREAM_COMMIT}.",
        "local_path": rel,
    }
    out_rows: list[dict[str, str]] = []
    replaced = False
    for row in rows:
        if row.get("local_path") == rel:
            out_rows.append(new_row)
            replaced = True
        else:
            out_rows.append(row)
    if not replaced:
        out_rows.append(new_row)
    write_tsv(DOWNLOAD_MANIFEST, out_rows, fields)


def parse_method_provenance() -> dict[str, str]:
    text = read_text_checked(COALESCENCE_README)
    required = [
        "single-copy genes",
        "raxml-ng --all",
        "astral.5.14.3.jar",
        "--reps 100",
        'INPUTTREES="astral.chr1-15/input.10SNPs.AR142-AR66-pruned.tre"',
    ]
    for snippet in required:
        if snippet not in text:
            raise Stage3Error(f"coalescence README missing expected method snippet: {snippet}")
    return {
        "tree_id": "Astral.10SNP-genes.chr1-15",
        "genomic_source": "chromosomes 1-15",
        "n_genes": "1631 (reported in source paper/upstream context; not stated in frozen Coalescence-based phylogenies/README.md)",
        "gene_tree_method": "RAxML-NG per single-copy gene using best-fit partitions; retained 10-SNP informative gene trees",
        "summary_method": "ASTRAL consensus tree with bootstrap support",
        "astral_version": "ASTRAL-III 5.14.3",
        "bootstrap_replicates": "100",
        "rooting": "published tree rooted using S. geminata for visualization/interpretation",
        "biological_role": "species/background phylogeny from chromosomes 1-15",
        "source_file": COALESCENCE_README.relative_to(REPO_ROOT).as_posix(),
        "source_commit": UPSTREAM_COMMIT,
        "notes": "Frozen README states the aim was to construct a supergene-region tree and a species phylogeny from chromosomes 1 to 15; Stage 3 uses only the chr1-15 tree.",
    }


def background_tree_metrics() -> dict[str, object]:
    text = read_text_checked(RAW_BACKGROUND_TREE).strip()
    parsed = parse_newick(text)
    return {
        "n_tips": len(parsed.tips),
        "tip_set_sha256": sha256_text("\n".join(sorted(parsed.tips)) + "\n"),
        "tree_sha256": sha256(RAW_BACKGROUND_TREE),
        "has_duplicate_tips": "false",
    }


def infer_species_relationship_from_background_tree() -> dict[str, object]:
    text = read_text_checked(RAW_BACKGROUND_TREE).strip()
    parsed = parse_newick(text)
    tips = set(parsed.tips)
    richteri = frozenset(t for t in tips if "_ric-" in t)
    inv_mac = frozenset(t for t in tips if "_inv-" in t or "_mac-" in t)
    if not richteri or not inv_mac:
        raise Stage3Error("background tree lacks expected richteri or invicta/macdonaghi labels")
    return {
        "relationship": EXPECTED_BACKGROUND_RELATIONSHIP,
        "richteri_tips": len(richteri),
        "invicta_macdonaghi_tips": len(inv_mac),
        "relationship_basis": "published chromosome 1-15 species/background interpretation plus frozen Stage-2 group membership",
        "individual_tree_literal_split_check": "not_used_for_stage3_partition",
        "notes": "Literal _ric-, _inv-, and _mac- labels were counted as a sanity check only. Stage 3 does not derive a new individual-level quartet or require SB/Sb subclades in the chr1-15 tree.",
    }


def write_background_tree_outputs() -> None:
    shutil.copyfile(RAW_BACKGROUND_TREE, PROCESSED_BACKGROUND_TREE)
    PROCESSED_BACKGROUND_TREE_SHA.write_text(f"{sha256(PROCESSED_BACKGROUND_TREE)}  {PROCESSED_BACKGROUND_TREE.relative_to(REPO_ROOT).as_posix()}\n")


def background_partition_rows() -> list[dict[str, object]]:
    source = PROCESSED_BACKGROUND_TREE.relative_to(REPO_ROOT).as_posix()
    return [
        {
            "group": "inv_mac_SB",
            "species_partition": "invicta/macdonaghi",
            "role": "A",
            "background_side": "side_1",
            "background_tree_source": source,
            "evidence_type": "published chr1-15 species/background relationship plus frozen Stage-2 group membership",
            "notes": "SB/Sb state not inferred from tree topology.",
        },
        {
            "group": "inv_mac_Sb",
            "species_partition": "invicta/macdonaghi",
            "role": "B",
            "background_side": "side_1",
            "background_tree_source": source,
            "evidence_type": "published chr1-15 species/background relationship plus frozen Stage-2 group membership",
            "notes": "SB/Sb state not inferred from tree topology.",
        },
        {
            "group": "richteri_SB",
            "species_partition": "richteri",
            "role": "C",
            "background_side": "side_2",
            "background_tree_source": source,
            "evidence_type": "published chr1-15 species/background relationship plus frozen Stage-2 group membership",
            "notes": "SB/Sb state not inferred from tree topology.",
        },
        {
            "group": "richteri_Sb",
            "species_partition": "richteri",
            "role": "D",
            "background_side": "side_2",
            "background_tree_source": source,
            "evidence_type": "published chr1-15 species/background relationship plus frozen Stage-2 group membership",
            "notes": "SB/Sb state not inferred from tree topology.",
        },
    ]


def stage0_species_split() -> str:
    _, rows = read_tsv(STAGE0_QUARTET)
    for row in rows:
        if row.get("resolution_id") == "species_split":
            return row.get("split", "")
    raise Stage3Error("Stage-0 focal quartet lacks species_split row")


def write_partition_sha() -> None:
    BACKGROUND_PARTITION_SHA.write_text(f"{sha256(BACKGROUND_PARTITION)}  {BACKGROUND_PARTITION.relative_to(REPO_ROOT).as_posix()}\n")


def write_report(metrics: dict[str, object], relationship: dict[str, object], stage0_match: bool) -> None:
    text = [
        "# Fire ants chromosome 16 Stage 3 report",
        "",
        "## Purpose",
        "",
        "Stage 3 freezes the independent chromosome 1-15 background species-history relationship before the first formal local TWISST support test.",
        "",
        "## Prior-stage checksum validation",
        "",
        f"- Stage 0 focal quartet: `{EXPECTED_STAGE0_FOCAL_SHA}`",
        f"- Stage 1 manifest: `{EXPECTED_STAGE1_MANIFEST_SHA}`",
        f"- Stage 2 manifest: `{EXPECTED_STAGE2_MANIFEST_SHA}`",
        f"- Stage 2 focal partition: `{EXPECTED_STAGE2_FOCAL_PARTITION_SHA}`",
        "",
        "All recorded prior-stage raw and processed checksums were verified before Stage 3 outputs were written.",
        "",
        "## Published background-tree provenance",
        "",
        "The frozen upstream `Coalescence-based phylogenies/README.md` states that the authors constructed a tree representing the supergene region and a tree representing the species phylogeny from chromosomes 1 to 15. It records single-copy-gene alignments, RAxML-NG per-gene phylogenies, ASTRAL `5.14.3`, and 100 ASTRAL bootstrap replicates.",
        "",
        "Stage 3 uses only the published chromosome 1-15 ASTRAL tree:",
        "",
        f"```text\n{RAW_BACKGROUND_TREE.relative_to(REPO_ROOT).as_posix()}\n```",
        "",
        f"Tree tips: `{metrics['n_tips']}`.",
        "",
        "## Background species relationship",
        "",
        "The independent background relationship frozen for the focal comparison is:",
        "",
        "```text\ninvicta/macdonaghi | richteri\n```",
        "",
        f"The chr1-15 tree contains `{relationship['invicta_macdonaghi_tips']}` literal `_inv-`/`_mac-` tips and `{relationship['richteri_tips']}` literal `_ric-` tips. These counts are a sanity check only; the Stage-3 partition is frozen from the published chromosome 1-15 species/background interpretation plus Stage-2 group membership.",
        "",
        "This is a species-level background statement, not a new representative-individual quartet analysis.",
        "",
        "## Mapping to frozen focal quartet",
        "",
        "Combining the background species relationship with frozen Stage-2 group membership gives:",
        "",
        "```text\n(inv_mac_SB, inv_mac_Sb) | (richteri_SB, richteri_Sb)\n```",
        "",
        f"Resulting focal split: `{EXPECTED_STAGE3_SPLIT}`.",
        "",
        f"Stage-3 background split matches frozen Stage-0 `species_split`: `{str(stage0_match).lower()}`.",
        "",
        "The independently published chromosome 1-15 species history supports the predeclared focal species partition AB|CD.",
        "",
        "## Anti-circularity",
        "",
        "Background species tree: derived from the published chromosome 1-15 ASTRAL analysis.",
        "",
        "Supergene/local genealogy: not inspected in Stage 3.",
        "",
        "Stage 3 did not read the supergene ASTRAL tree, read TWISST weights, classify grouped TWISST topologies, calculate q-support statistics, select representatives by tree position, compare regions, or run ASTRAL/ASTER.",
        "",
        "## Stage boundary",
        "",
        "Stage 4A is the first formal unblinded local-topology stage.",
    ]
    STAGE3_REPORT.write_text("\n".join(text) + "\n")


def write_manifest(stage1: dict[str, object], stage2: dict[str, object], metrics: dict[str, object], relationship: dict[str, object], stage0_match: bool) -> None:
    output_checksums = {path.relative_to(REPO_ROOT).as_posix(): sha256(path) for path in OUTPUTS if path != STAGE3_MANIFEST}
    data = {
        "stage": "Stage 3",
        "description": "Freeze independent chromosome 1-15 background species-history relationship without local TWISST support analysis.",
        "stage0_focal_quartet_checksum": EXPECTED_STAGE0_FOCAL_SHA,
        "stage1_manifest_checksum": sha256(STAGE1_MANIFEST),
        "stage2_manifest_checksum": sha256(STAGE2_MANIFEST),
        "stage2_focal_partition_checksum": sha256(STAGE2_FOCAL_PARTITION),
        "upstream_repository": UPSTREAM_REPO,
        "upstream_branch": UPSTREAM_BRANCH,
        "upstream_commit": UPSTREAM_COMMIT,
        "raw_chr1_15_astral_tree": {
            "path": RAW_BACKGROUND_TREE.relative_to(REPO_ROOT).as_posix(),
            "source_url": BACKGROUND_TREE_URL,
            "sha256": sha256(RAW_BACKGROUND_TREE),
        },
        "normalized_background_tree_checksum": sha256(PROCESSED_BACKGROUND_TREE),
        "background_tree_provenance_checksum": sha256(BACKGROUND_PROVENANCE),
        "stage3_background_partition_checksum": sha256(BACKGROUND_PARTITION),
        "stage3_background_partition_matches_stage0_species_split": stage0_match,
        "background_tree_metrics": metrics,
        "background_species_relationship": relationship,
        "stage1_prior_manifest_reference": {
            "path": STAGE1_MANIFEST.relative_to(REPO_ROOT).as_posix(),
            "recorded_metrics": stage1.get("observed_metrics", {}),
        },
        "stage2_prior_manifest_reference": {
            "path": STAGE2_MANIFEST.relative_to(REPO_ROOT).as_posix(),
            "focal_species_partition": stage2.get("focal_species_partition"),
            "focal_haplotype_partition": stage2.get("focal_haplotype_partition"),
        },
        "processing_script": {
            "path": "empirical/fire_ants_chr16/scripts/03_freeze_background_history.py",
            "sha256": sha256(Path(__file__)),
        },
        "supergene_astral_tree_read": False,
        "twisst_topology_files_read": False,
        "twisst_weight_files_read": False,
        "q_support_statistics_calculated": False,
        "outputs": [path.relative_to(REPO_ROOT).as_posix() for path in OUTPUTS],
    }
    STAGE3_MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate_output_schemas() -> None:
    for path in [BACKGROUND_PROVENANCE, BACKGROUND_PARTITION]:
        fields, _ = read_tsv(path)
        forbidden = {field.lower() for field in fields} & FORBIDDEN_OUTPUT_FIELDS
        if forbidden:
            raise Stage3Error(f"{path.relative_to(REPO_ROOT)} contains forbidden Stage-3 output fields: {sorted(forbidden)}")


def run_stage3(skip_download: bool = False) -> dict[str, object]:
    stage1, stage2 = verify_prior_stages()
    download_background_tree(skip_download=skip_download)
    update_download_manifest()
    if sha256(RAW_BACKGROUND_TREE) != EXPECTED_BACKGROUND_TREE_SHA:
        raise Stage3Error("raw background tree checksum mismatch")
    write_background_tree_outputs()
    metrics = background_tree_metrics()
    provenance = parse_method_provenance()
    relationship = infer_species_relationship_from_background_tree()
    write_tsv(BACKGROUND_PROVENANCE, [provenance], ["tree_id", "genomic_source", "n_genes", "gene_tree_method", "summary_method", "astral_version", "bootstrap_replicates", "rooting", "biological_role", "source_file", "source_commit", "notes"])
    write_tsv(BACKGROUND_PARTITION, background_partition_rows(), ["group", "species_partition", "role", "background_side", "background_tree_source", "evidence_type", "notes"])
    write_partition_sha()
    if relationship["relationship"] != EXPECTED_BACKGROUND_RELATIONSHIP:
        raise Stage3Error("background species relationship did not match expected invicta/macdonaghi|richteri")
    split = stage0_species_split()
    stage0_match = split == EXPECTED_STAGE3_SPLIT
    if not stage0_match:
        raise Stage3Error(f"Stage-3 background split {EXPECTED_STAGE3_SPLIT} does not match Stage-0 species_split {split}")
    validate_output_schemas()
    write_report(metrics, relationship, stage0_match)
    write_manifest(stage1, stage2, metrics, relationship, stage0_match)
    return {
        "n_tips": metrics["n_tips"],
        "background_tree_sha": sha256(RAW_BACKGROUND_TREE),
        "background_partition_sha": sha256(BACKGROUND_PARTITION),
        "manifest_sha": sha256(STAGE3_MANIFEST),
    }


class Stage3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.metrics = run_stage3(skip_download=True)

    def test_stage0_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE0_FOCAL_QUARTET), EXPECTED_STAGE0_FOCAL_SHA)

    def test_stage1_manifest_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE1_MANIFEST), EXPECTED_STAGE1_MANIFEST_SHA)

    def test_stage2_manifest_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE2_MANIFEST), EXPECTED_STAGE2_MANIFEST_SHA)

    def test_stage2_focal_partition_checksum_unchanged(self) -> None:
        self.assertEqual(sha256(STAGE2_FOCAL_PARTITION), EXPECTED_STAGE2_FOCAL_PARTITION_SHA)

    def test_frozen_upstream_commit_unchanged(self) -> None:
        stage1 = json.loads(read_text_checked(STAGE1_MANIFEST))
        self.assertEqual(stage1["upstream_commit"], UPSTREAM_COMMIT)

    def test_background_astral_tree_hash_recorded(self) -> None:
        manifest = json.loads(read_text_checked(STAGE3_MANIFEST))
        self.assertEqual(manifest["raw_chr1_15_astral_tree"]["sha256"], EXPECTED_BACKGROUND_TREE_SHA)

    def test_background_tree_parses(self) -> None:
        parse_newick(read_text_checked(RAW_BACKGROUND_TREE))

    def test_no_duplicate_tips(self) -> None:
        parsed = parse_newick(read_text_checked(RAW_BACKGROUND_TREE))
        self.assertEqual(len(parsed.tips), len(set(parsed.tips)))

    def test_tip_count_deterministic(self) -> None:
        parsed = parse_newick(read_text_checked(RAW_BACKGROUND_TREE))
        self.assertEqual(len(parsed.tips), self.metrics["n_tips"])

    def test_provenance_records_1631_chr1_15_genes_if_confirmed(self) -> None:
        _, rows = read_tsv(BACKGROUND_PROVENANCE)
        self.assertIn("1631", rows[0]["n_genes"])
        self.assertIn("chromosomes 1-15", rows[0]["genomic_source"])

    def test_astral_version_recorded(self) -> None:
        _, rows = read_tsv(BACKGROUND_PROVENANCE)
        self.assertEqual(rows[0]["astral_version"], "ASTRAL-III 5.14.3")

    def test_published_background_species_relationship(self) -> None:
        relationship = infer_species_relationship_from_background_tree()
        self.assertEqual(relationship["relationship"], EXPECTED_BACKGROUND_RELATIONSHIP)
        self.assertEqual(relationship["relationship_basis"], "published chromosome 1-15 species/background interpretation plus frozen Stage-2 group membership")

    def test_stage3_mapping_implies_ab_cd(self) -> None:
        _, rows = read_tsv(BACKGROUND_PARTITION)
        sides: dict[str, set[str]] = {}
        for row in rows:
            sides.setdefault(row["background_side"], set()).add(row["role"])
        self.assertEqual(set(map(frozenset, sides.values())), {frozenset({"A", "B"}), frozenset({"C", "D"})})

    def test_stage3_split_matches_stage0_species_split(self) -> None:
        self.assertEqual(stage0_species_split(), EXPECTED_STAGE3_SPLIT)

    def test_no_stage1_twisst_weight_file_is_read(self) -> None:
        with self.assertRaises(Stage3Error):
            read_text_checked(DATA_ROOT / "raw" / "upstream" / "Topology weighting" / "results" / "2021-08-07-twisst" / "weights_output.csv")
        with self.assertRaises(Stage3Error):
            read_text_checked(DATA_ROOT / "processed" / "stage1_twisst_weights_sparse.tsv")

    def test_no_stage1_grouped_topology_file_is_read(self) -> None:
        with self.assertRaises(Stage3Error):
            read_text_checked(DATA_ROOT / "raw" / "upstream" / "Topology weighting" / "results" / "2021-08-07-twisst" / "topologies_output.trees")
        with self.assertRaises(Stage3Error):
            read_text_checked(DATA_ROOT / "processed" / "stage1_twisst_topologies.tsv")

    def test_no_supergene_astral_tree_is_read(self) -> None:
        with self.assertRaises(Stage3Error):
            read_text_checked(DATA_ROOT / "raw" / "upstream" / "Trees in newick format" / "Astral.10SNP-genes.supergene.nwk")

    def test_no_q_support_statistic_is_produced(self) -> None:
        validate_output_schemas()
        manifest = json.loads(read_text_checked(STAGE3_MANIFEST))
        self.assertFalse(manifest["q_support_statistics_calculated"])

    def test_repeated_runs_byte_stable(self) -> None:
        before = {path: sha256(path) for path in OUTPUTS}
        run_stage3(skip_download=True)
        after = {path: sha256(path) for path in OUTPUTS}
        self.assertEqual(before, after)


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage3Tests)
    with tempfile.TemporaryDirectory(prefix="fire_ants_stage3_tests_"):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true", help="regenerate Stage-3 outputs from immutable raw inputs")
    parser.add_argument("--run-tests", action="store_true", help="run embedded Stage-3 tests after freezing outputs")
    args = parser.parse_args(argv)
    try:
        metrics = run_stage3(skip_download=args.skip_download)
    except Stage3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.run_tests:
        return 0 if run_tests() else 1
    print(f"Wrote {BACKGROUND_PARTITION.relative_to(REPO_ROOT)}")
    print(f"Parsed {metrics['n_tips']} background-tree tips")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
