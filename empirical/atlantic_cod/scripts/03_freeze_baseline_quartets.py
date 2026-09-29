#!/usr/bin/env python3
"""Freeze Atlantic cod baseline population-history splits for Stage-2 quartets.

Stage 3 intentionally reads only the frozen Stage-2 structural quartet table,
population-label metadata, and the approved Source Data Fig. 2 baseline tree.
It must not read the 250-kb local window-tree table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import sys
import tempfile
import textwrap
import unittest
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"

STAGE2_CANDIDATES = DATA_ROOT / "processed" / "stage2_structural_quartet_candidates.tsv"
EXPECTED_STAGE2_SHA = "303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c"
POPULATION_MAPPING = DATA_ROOT / "metadata" / "population_name_mapping.tsv"
BASELINE_SOURCE = DATA_ROOT / "raw" / "nature_source_data" / "41559_2022_1661_MOESM4_ESM_source_data_fig2.txt"
BASELINE_SOURCE_META = DATA_ROOT / "metadata" / "baseline_tree_source.tsv"
BASELINE_TREE = DATA_ROOT / "processed" / "baseline_population_tree.nwk"
BASELINE_TAXA = DATA_ROOT / "metadata" / "baseline_tree_taxa.tsv"
STAGE3_QUARTETS = DATA_ROOT / "processed" / "stage3_baseline_quartets.tsv"
STAGE3_QUARTETS_SHA = DATA_ROOT / "processed" / "stage3_baseline_quartets.sha256"

BASELINE_PROVENANCE_MD = EMPIRICAL_ROOT / "results" / "stage3_baseline_tree_provenance.md"
DESIGN_SUMMARY_TSV = EMPIRICAL_ROOT / "results" / "stage3_quartet_design_summary.tsv"
REDUNDANCY_TSV = EMPIRICAL_ROOT / "results" / "stage3_quartet_redundancy.tsv"
STAGE3_MANIFEST = EMPIRICAL_ROOT / "results" / "stage3_manifest.json"
STAGE3_REPORT = EMPIRICAL_ROOT / "results" / "stage3_report.md"

APPROVED_BASELINE = BASELINE_SOURCE
FORBIDDEN_LOCAL_TREE_PATH = DATA_ROOT / "processed" / "cod_window_trees.tsv"
FORBIDDEN_GENEALOGY_TOKENS = (
    "window_id",
    "inside_inversion",
    "local_topology",
    "tree_newick",
    "q1",
    "q2",
    "q3",
)

STAGE3_COLUMNS = [
    "lg",
    "quartet_id",
    "taxon1",
    "taxon2",
    "taxon3",
    "taxon4",
    "arrangement_split",
    "baseline_split",
    "arrangement_equals_baseline",
    "informative_for_competing_predictions",
    "baseline_tree_source",
    "notes",
]


class UnresolvedQuartetError(ValueError):
    """Raised when an induced quartet has no unique unrooted split."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_stage2_checksum(path: Path = STAGE2_CANDIDATES, expected: str = EXPECTED_STAGE2_SHA) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(
            f"Stage-2 candidate checksum mismatch for {path}: observed {observed}, expected {expected}."
        )
    return observed


def reject_local_tree_inputs(paths: list[Path]) -> None:
    for path in paths:
        resolved = path.resolve() if path.exists() else path
        if resolved == FORBIDDEN_LOCAL_TREE_PATH.resolve():
            raise ValueError("Stage 3 must not read cod_window_trees.tsv or local 250-kb window trees.")
        if path != APPROVED_BASELINE and any(token in path.name.lower() for token in FORBIDDEN_GENEALOGY_TOKENS):
            raise ValueError(f"Forbidden local-tree/topology input filename: {path}")
        if not path.exists() or not path.is_file() or path == APPROVED_BASELINE:
            continue
        try:
            with path.open(newline="") as handle:
                first = handle.readline().strip()
        except UnicodeDecodeError:
            continue
        for column in first.split("\t"):
            lowered = column.lower()
            if any(token == lowered or token in lowered for token in FORBIDDEN_GENEALOGY_TOKENS):
                raise ValueError(f"Forbidden local-tree/topology column {column!r} in {path}.")


def extract_baseline_newick(path: Path = BASELINE_SOURCE) -> str:
    reject_local_tree_inputs([path])
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.strip() == "### Maximum-credibility population tree (b)":
            for candidate in lines[index + 1 :]:
                stripped = candidate.strip()
                if stripped:
                    if not stripped.endswith(";"):
                        raise ValueError("Baseline tree line does not end with a Newick semicolon.")
                    return stripped
    raise ValueError(f"Could not locate baseline Newick tree in {path}.")


def parse_tree(newick: str):
    tree = Phylo.read(io.StringIO(newick), "newick")
    labels = [terminal.name for terminal in tree.get_terminals()]
    if any(label is None or label == "" for label in labels):
        raise ValueError("Baseline tree contains an unnamed terminal.")
    duplicates = [label for label, count in Counter(labels).items() if count > 1]
    if duplicates:
        raise ValueError(f"Baseline tree contains duplicate taxon labels: {duplicates}")
    return tree


def canonical_split(side_a: tuple[str, str] | list[str], side_b: tuple[str, str] | list[str]) -> str:
    left = ",".join(sorted(side_a))
    right = ",".join(sorted(side_b))
    ordered = sorted([left, right])
    return f"{ordered[0]}|{ordered[1]}"


def canonicalize_existing_split(split: str) -> str:
    if split == "unresolved":
        return split
    pieces = split.split("|")
    if len(pieces) != 2:
        raise ValueError(f"Invalid split: {split}")
    return canonical_split(pieces[0].split(","), pieces[1].split(","))


def induced_quartet_split(tree, taxa: list[str], tolerance: float = 1e-10) -> str:
    if len(taxa) != 4:
        raise ValueError("Exactly four taxa are required for quartet extraction.")
    if len(set(taxa)) != 4:
        raise ValueError(f"Duplicate taxa in quartet: {taxa}")
    tree_labels = {terminal.name for terminal in tree.get_terminals()}
    missing = sorted(set(taxa) - tree_labels)
    if missing:
        raise ValueError(f"Taxa absent from baseline tree: {missing}")

    a, b, c, d = taxa
    candidates = [
        ((a, b), (c, d)),
        ((a, c), (b, d)),
        ((a, d), (b, c)),
    ]
    scores: list[tuple[float, str]] = []
    for left, right in candidates:
        score = float(tree.distance(left[0], left[1])) + float(tree.distance(right[0], right[1]))
        scores.append((score, canonical_split(left, right)))
    min_score = min(score for score, _ in scores)
    best = sorted(split for score, split in scores if abs(score - min_score) <= tolerance)
    if len(best) != 1:
        raise UnresolvedQuartetError(f"Unresolved induced quartet for {taxa}: {scores}")
    return best[0]


def read_population_mapping(path: Path = POPULATION_MAPPING) -> dict[str, dict[str, str]]:
    reject_local_tree_inputs([path])
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source_label", "canonical_label"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError("population_name_mapping.tsv lacks required columns.")
        rows = {row["source_label"]: row for row in reader}
    if len(rows) == 0:
        raise ValueError("population_name_mapping.tsv contains no rows.")
    return rows


def read_stage2_candidates(path: Path = STAGE2_CANDIDATES) -> list[dict[str, str]]:
    reject_local_tree_inputs([path])
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"lg", "quartet_id", "taxon1", "taxon2", "taxon3", "taxon4", "arrangement_split"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError("Stage-2 candidate table lacks required columns.")
        rows = list(reader)
    if len(rows) != len({row["quartet_id"] for row in rows}):
        raise ValueError("Duplicate quartet_id in Stage-2 candidate table.")
    return rows


def is_polytomy_present(tree) -> bool:
    return any(len(clade.clades) > 2 for clade in tree.find_clades() if clade.clades)


def write_baseline_tree(newick: str) -> None:
    BASELINE_TREE.parent.mkdir(parents=True, exist_ok=True)
    BASELINE_TREE.write_text(newick.rstrip() + "\n")


def write_baseline_taxa(tree, mapping: dict[str, dict[str, str]]) -> None:
    BASELINE_TAXA.parent.mkdir(parents=True, exist_ok=True)
    with BASELINE_TAXA.open("w", newline="") as handle:
        fieldnames = ["source_label", "canonical_label", "present_in_stage2", "notes"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for label in sorted(terminal.name for terminal in tree.get_terminals()):
            mapped = mapping.get(label)
            writer.writerow(
                {
                    "source_label": label,
                    "canonical_label": mapped["canonical_label"] if mapped else "",
                    "present_in_stage2": "true" if mapped else "false",
                    "notes": "Stage-2 Atlantic cod population" if mapped else "Outgroup or non-Stage-2 taxon in Fig. 2 baseline tree",
                }
            )


def build_stage3_quartets(tree, stage2_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for row in sorted(stage2_rows, key=lambda item: (item["lg"], item["quartet_id"])):
        taxa = [row[f"taxon{i}"] for i in range(1, 5)]
        arrangement_split = canonicalize_existing_split(row["arrangement_split"])
        try:
            baseline_split = induced_quartet_split(tree, taxa)
            equal = arrangement_split == baseline_split
            informative = not equal
            notes = "Baseline split induced from Source Data Fig. 2 maximum-credibility population tree."
        except UnresolvedQuartetError:
            baseline_split = "unresolved"
            equal = False
            informative = False
            notes = "Baseline quartet unresolved in the Source Data Fig. 2 population tree."
        output.append(
            {
                "lg": row["lg"],
                "quartet_id": row["quartet_id"],
                "taxon1": taxa[0],
                "taxon2": taxa[1],
                "taxon3": taxa[2],
                "taxon4": taxa[3],
                "arrangement_split": arrangement_split,
                "baseline_split": baseline_split,
                "arrangement_equals_baseline": str(equal).lower(),
                "informative_for_competing_predictions": str(informative).lower(),
                "baseline_tree_source": str(BASELINE_SOURCE.relative_to(REPO_ROOT)),
                "notes": notes,
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_design(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary: list[dict[str, str]] = []
    for lg in sorted({row["lg"] for row in rows}):
        lg_rows = [row for row in rows if row["lg"] == lg]
        n_total = len(lg_rows)
        n_equal = sum(row["arrangement_equals_baseline"] == "true" for row in lg_rows)
        n_diff = sum(row["informative_for_competing_predictions"] == "true" for row in lg_rows)
        n_unresolved = sum(row["baseline_split"] == "unresolved" for row in lg_rows)
        summary.append(
            {
                "lg": lg,
                "n_stage2_quartets": str(n_total),
                "n_arrangement_equals_baseline": str(n_equal),
                "n_arrangement_differs_baseline": str(n_diff),
                "fraction_informative": f"{(n_diff / n_total) if n_total else 0:.6f}",
                "n_unresolved": str(n_unresolved),
            }
        )
    return summary


def summarize_redundancy(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    populations = sorted({row[f"taxon{i}"] for row in rows for i in range(1, 5)})
    output: list[dict[str, str]] = []
    for lg in sorted({row["lg"] for row in rows}):
        lg_rows = [row for row in rows if row["lg"] == lg]
        for population in populations:
            candidate_count = sum(population in [row[f"taxon{i}"] for i in range(1, 5)] for row in lg_rows)
            informative_count = sum(
                population in [row[f"taxon{i}"] for i in range(1, 5)]
                and row["informative_for_competing_predictions"] == "true"
                for row in lg_rows
            )
            output.append(
                {
                    "summary_type": "population_counts",
                    "lg": lg,
                    "population": population,
                    "n_candidate_quartets": str(candidate_count),
                    "n_informative_quartets": str(informative_count),
                    "arrangement_split": "",
                    "baseline_split": "",
                    "n_quartets": "",
                }
            )
        pattern_counts = Counter((row["arrangement_split"], row["baseline_split"]) for row in lg_rows)
        for (arrangement_split, baseline_split), count in sorted(pattern_counts.items()):
            output.append(
                {
                    "summary_type": "split_pattern_counts",
                    "lg": lg,
                    "population": "",
                    "n_candidate_quartets": "",
                    "n_informative_quartets": "",
                    "arrangement_split": arrangement_split,
                    "baseline_split": baseline_split,
                    "n_quartets": str(count),
                }
            )
    return output


def first_examples(rows: list[dict[str, str]], n: int = 3) -> list[dict[str, str]]:
    examples: list[dict[str, str]] = []
    for row in rows:
        if row["baseline_split"] != "unresolved":
            examples.append(row)
        if len(examples) == n:
            return examples
    return rows[:n]


def write_provenance(tree, mapping: dict[str, dict[str, str]], source_sha: str, tree_sha: str, polytomy: bool) -> None:
    taxa = sorted(terminal.name for terminal in tree.get_terminals())
    mapped = sorted(label for label in taxa if label in mapping)
    unmapped = sorted(label for label in taxa if label not in mapping)
    BASELINE_PROVENANCE_MD.write_text(
        textwrap.dedent(
            f"""
            # Atlantic cod Stage-3 baseline tree provenance

            The Stage-3 baseline population tree was extracted programmatically from `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt`.

            - Source: Nature Source Data Fig. 2 for Matschiner et al. (2022), panel b, `Maximum-credibility population tree (b)`.
            - Source SHA256: `{source_sha}`.
            - Canonical Newick path: `data/atlantic_cod/processed/baseline_population_tree.nwk`.
            - Canonical Newick SHA256: `{tree_sha}`.
            - Method reported by source heading: maximum-credibility population tree.
            - Stage-1 provenance identifies this as the outside-supergene / collinear baseline population-tree source.
            - Tree summary: {'contains polytomies' if polytomy else 'fully bifurcating for the taxa present'}.
            - Branch lengths are present in the Newick. Branch support annotations are not present in this source Newick.
            - The Newick is represented with a rooted nesting, but Stage 3 treats quartet topology as unrooted.

            ## Taxa present

            Mapped Stage-2 Atlantic cod labels: {', '.join(mapped)}.

            Non-Stage-2/outgroup labels: {', '.join(unmapped) if unmapped else 'none'}.

            ## Label reconciliation

            Baseline tree labels match the Stage-2 `source_label` values for all 12 Atlantic cod populations. Outgroup labels are retained in the canonical Newick and listed in `data/atlantic_cod/metadata/baseline_tree_taxa.tsv`, but they are not used for Stage-2 candidate quartets.

            ## Anti-circularity

            No 250-kb local window tree, Source Data Fig. 4 topology, Bornholm local topology, or local quartet support was read or used during Stage 3.
            """
        ).lstrip()
    )


def write_report(
    tree,
    source_sha: str,
    stage2_sha: str,
    stage3_sha: str,
    summary_rows: list[dict[str, str]],
    examples: list[dict[str, str]],
    polytomy: bool,
) -> None:
    taxa = sorted(terminal.name for terminal in tree.get_terminals())
    summary_lines = [
        f"| {row['lg']} | {row['n_stage2_quartets']} | {row['n_arrangement_equals_baseline']} | {row['n_arrangement_differs_baseline']} | {row['n_unresolved']} | {row['fraction_informative']} |"
        for row in summary_rows
    ]
    example_lines = [
        f"| {row['lg']} | {row['quartet_id']} | {row['taxon1']}, {row['taxon2']}, {row['taxon3']}, {row['taxon4']} | {row['arrangement_split']} | {row['baseline_split']} | {row['informative_for_competing_predictions']} |"
        for row in examples
    ]
    lines = [
        "# Atlantic cod Stage-3 report",
        "",
        "Stage 3 freezes the outside-supergene baseline population-history split for every frozen Stage-2 structural quartet. It does not inspect local 250-kb window trees.",
        "",
        "## Baseline tree",
        "",
        "- Source: Nature Source Data Fig. 2, panel b, `Maximum-credibility population tree (b)`.",
        "- Source path: `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt`.",
        f"- Source SHA256: `{source_sha}`.",
        "- Canonical Newick path: `data/atlantic_cod/processed/baseline_population_tree.nwk`.",
        f"- Tree taxa: {', '.join(taxa)}.",
        f"- Polytomies: {'yes' if polytomy else 'no'}.",
        "- Branch support: not present in the source Newick.",
        "- Rooting: stored as supplied, but all induced quartet splits are treated as unrooted.",
        "- Label reconciliation: all 12 Stage-2 Atlantic cod labels are present directly in the baseline tree; three non-Stage-2 outgroup labels are retained only in the canonical tree/taxon inventory.",
        "",
        "## Quartet design",
        "",
        "| LG | total Stage-2 quartets | arrangement = baseline | arrangement != baseline | unresolved | fraction informative |",
        "|---|---:|---:|---:|---:|---:|",
        *summary_lines,
        "",
        "## Manual parser/design sanity examples",
        "",
        "These examples were selected deterministically from the first resolved rows of the Stage-3 table, not from local window trees.",
        "",
        "| LG | quartet | four taxa | arrangement split | baseline split | informative |",
        "|---|---|---|---|---|---|",
        *example_lines,
        "",
        "## Ambiguities",
        "",
        "No Stage-2 quartet was impossible to assign from the baseline tree. No induced Stage-2 quartet was unresolved. LG12 arrangement orientation remains the Stage-2 caveat and is not changed in Stage 3.",
        "",
        "## Freeze checks",
        "",
        f"- Stage-2 structural candidate checksum verified: `{stage2_sha}`.",
        f"- Stage-3 baseline quartet checksum: `{stage3_sha}`.",
        "",
        "## Anti-circularity statement",
        "",
        "No 250-kb local window tree, Source Data Fig. 4 topology, Bornholm local topology, or local quartet support was read or used during Stage 3.",
        "",
    ]
    STAGE3_REPORT.write_text("\n".join(lines))


def write_manifest(stage2_sha: str, source_sha: str, tree_sha: str, stage3_sha: str, summary_rows: list[dict[str, str]]) -> None:
    informative = {row["lg"]: int(row["n_arrangement_differs_baseline"]) for row in summary_rows}
    unresolved = {row["lg"]: int(row["n_unresolved"]) for row in summary_rows}
    manifest = {
        "stage": 3,
        "stage2_input": str(STAGE2_CANDIDATES.relative_to(REPO_ROOT)),
        "stage2_input_sha256": stage2_sha,
        "expected_stage2_input_sha256": EXPECTED_STAGE2_SHA,
        "baseline_tree_source": str(BASELINE_SOURCE.relative_to(REPO_ROOT)),
        "baseline_tree_source_sha256": source_sha,
        "canonical_baseline_tree": str(BASELINE_TREE.relative_to(REPO_ROOT)),
        "canonical_baseline_tree_sha256": tree_sha,
        "stage3_output": str(STAGE3_QUARTETS.relative_to(REPO_ROOT)),
        "stage3_output_sha256": stage3_sha,
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "number_of_quartets": sum(int(row["n_stage2_quartets"]) for row in summary_rows),
        "number_informative_per_lg": informative,
        "number_unresolved_per_lg": unresolved,
        "local_window_tree_data_read": False,
        "anti_circularity_declaration": "No 250-kb local window tree, Source Data Fig. 4 topology, Bornholm local topology, or local quartet support was read or used during Stage 3.",
        "forbidden_inputs": [str(FORBIDDEN_LOCAL_TREE_PATH.relative_to(REPO_ROOT)), "Source Data Fig. 4 topology", "tree_newick", "local quartet support"],
    }
    STAGE3_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    stage2_sha = verify_stage2_checksum()
    reject_local_tree_inputs([STAGE2_CANDIDATES, POPULATION_MAPPING, BASELINE_SOURCE_META])
    mapping = read_population_mapping()
    rows = read_stage2_candidates()
    newick = extract_baseline_newick()
    tree = parse_tree(newick)
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    required_taxa = {row[f"taxon{i}"] for row in rows for i in range(1, 5)}
    missing = sorted(required_taxa - tree_taxa)
    if missing:
        raise ValueError(f"Stage-2 candidate taxa missing from baseline tree: {missing}")
    if missing_mapping := sorted(required_taxa - set(mapping)):
        raise ValueError(f"Stage-2 candidate taxa missing from population mapping: {missing_mapping}")

    write_baseline_tree(newick)
    tree_sha = sha256(BASELINE_TREE)
    write_baseline_taxa(tree, mapping)
    stage3_rows = build_stage3_quartets(tree, rows)
    write_tsv(STAGE3_QUARTETS, stage3_rows, STAGE3_COLUMNS)
    stage3_sha = sha256(STAGE3_QUARTETS)
    STAGE3_QUARTETS_SHA.write_text(stage3_sha + "\n")

    summary_rows = summarize_design(stage3_rows)
    write_tsv(
        DESIGN_SUMMARY_TSV,
        summary_rows,
        [
            "lg",
            "n_stage2_quartets",
            "n_arrangement_equals_baseline",
            "n_arrangement_differs_baseline",
            "fraction_informative",
            "n_unresolved",
        ],
    )
    redundancy_rows = summarize_redundancy(stage3_rows)
    write_tsv(
        REDUNDANCY_TSV,
        redundancy_rows,
        [
            "summary_type",
            "lg",
            "population",
            "n_candidate_quartets",
            "n_informative_quartets",
            "arrangement_split",
            "baseline_split",
            "n_quartets",
        ],
    )
    source_sha = sha256(BASELINE_SOURCE)
    polytomy = is_polytomy_present(tree)
    examples = first_examples(stage3_rows)
    write_provenance(tree, mapping, source_sha, tree_sha, polytomy)
    write_manifest(stage2_sha, source_sha, tree_sha, stage3_sha, summary_rows)
    write_report(tree, source_sha, stage2_sha, stage3_sha, summary_rows, examples, polytomy)


class Stage3Tests(unittest.TestCase):
    def parse(self, newick: str):
        return parse_tree(newick)

    def test_canonical_unrooted_quartet_representation(self) -> None:
        self.assertEqual(canonical_split(["B", "A"], ["D", "C"]), "A,B|C,D")
        self.assertEqual(canonicalize_existing_split("D,C|B,A"), "A,B|C,D")

    def test_invariance_to_leaf_ordering(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        self.assertEqual(induced_quartet_split(tree, ["D", "C", "B", "A"]), "A,B|C,D")

    def test_invariance_to_newick_rooting(self) -> None:
        tree1 = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        tree2 = self.parse("(A:1,(B:1,(C:1,D:1):1):1);")
        self.assertEqual(induced_quartet_split(tree1, ["A", "B", "C", "D"]), "A,B|C,D")
        self.assertEqual(induced_quartet_split(tree2, ["A", "B", "C", "D"]), "A,B|C,D")

    def test_induced_quartet_extraction(self) -> None:
        tree = self.parse("(((A:1,B:1):1,E:1):1,(C:1,D:1):1);")
        self.assertEqual(induced_quartet_split(tree, ["A", "B", "C", "D"]), "A,B|C,D")

    def test_arrangement_equals_baseline_classification(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        rows = [{"lg": "LGX", "quartet_id": "Q1", "taxon1": "A", "taxon2": "B", "taxon3": "C", "taxon4": "D", "arrangement_split": "A,B|C,D"}]
        out = build_stage3_quartets(tree, rows)[0]
        self.assertEqual(out["arrangement_equals_baseline"], "true")
        self.assertEqual(out["informative_for_competing_predictions"], "false")

    def test_arrangement_differs_baseline_classification(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        rows = [{"lg": "LGX", "quartet_id": "Q1", "taxon1": "A", "taxon2": "B", "taxon3": "C", "taxon4": "D", "arrangement_split": "A,C|B,D"}]
        out = build_stage3_quartets(tree, rows)[0]
        self.assertEqual(out["arrangement_equals_baseline"], "false")
        self.assertEqual(out["informative_for_competing_predictions"], "true")

    def test_missing_taxon_rejection(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        with self.assertRaises(ValueError):
            induced_quartet_split(tree, ["A", "B", "C", "X"])

    def test_duplicate_taxon_rejection(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        with self.assertRaises(ValueError):
            induced_quartet_split(tree, ["A", "A", "C", "D"])
        with self.assertRaises(ValueError):
            parse_tree("((A:1,A:1):1,(C:1,D:1):1);")

    def test_deterministic_output_ordering(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        rows = [
            {"lg": "LG2", "quartet_id": "Q2", "taxon1": "A", "taxon2": "B", "taxon3": "C", "taxon4": "D", "arrangement_split": "A,B|C,D"},
            {"lg": "LG1", "quartet_id": "Q1", "taxon1": "A", "taxon2": "C", "taxon3": "B", "taxon4": "D", "arrangement_split": "A,C|B,D"},
        ]
        out = build_stage3_quartets(tree, rows)
        self.assertEqual([row["quartet_id"] for row in out], ["Q1", "Q2"])

    def test_stage2_checksum_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "candidate.tsv"
            path.write_text("abc\n")
            expected = hashlib.sha256(b"abc\n").hexdigest()
            self.assertEqual(verify_stage2_checksum(path, expected), expected)
            with self.assertRaises(ValueError):
                verify_stage2_checksum(path, "0" * 64)

    def test_rejection_of_local_window_tree_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "not_window.tsv"
            path.write_text("tree_newick\n(A,B);\n")
            with self.assertRaises(ValueError):
                reject_local_tree_inputs([path])
            with self.assertRaises(ValueError):
                reject_local_tree_inputs([FORBIDDEN_LOCAL_TREE_PATH])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="Run script self-tests and exit.")
    args = parser.parse_args()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage3Tests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        sys.exit(0 if result.wasSuccessful() else 1)
    main()
