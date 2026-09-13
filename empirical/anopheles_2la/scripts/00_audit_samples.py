#!/usr/bin/env python3
"""Audit Anopheles 2La sample metadata and freeze structural-only quartet candidates.

This Stage-0 script uses only sample, geography, species, and 2La karyotype
metadata. It must not read local-tree, gene-tree, topology, or quartet-support
outputs.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import sys
import unittest
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "anopheles_2la"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "anopheles_2la"
SAMPLE_MANIFEST = DATA_ROOT / "metadata" / "sample_manifest.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"
AVAILABILITY = DATA_ROOT / "metadata" / "ag3_10_fontaine_rebuild_availability.tsv"
QUARTETS_OUT = DATA_ROOT / "processed" / "candidate_strict_quartets.tsv"
REPORT_OUT = EMPIRICAL_ROOT / "results" / "stage0_data_audit.md"

ALLOWED_STATES = {"A0_homozygous", "A1_homozygous", "heterokaryotype", "unknown"}
STRICT_STATES = {"A0_homozygous", "A1_homozygous"}
TOPOLOGY_LEAK_PATTERNS = (
    "topology",
    "gene_tree",
    "genetree",
    "local_tree",
    "quartet_support",
    "qqs",
    "q1",
    "q2",
    "q3",
    "fontaine_topology",
)
SAMPLE_FIELDS = [
    "sample_id",
    "species",
    "population",
    "country",
    "locality",
    "latitude_if_public",
    "longitude_if_public",
    "sex",
    "dataset",
    "dataset_version",
    "reference_build",
    "2La_karyotype",
    "2La_state",
    "karyotype_source",
    "karyotype_confidence",
    "phased_data_available",
    "sequence_data_available",
    "include_strict_test",
    "exclusion_reason",
    "notes",
]
QUARTET_FIELDS = [
    "quartet_id",
    "design_class",
    "sample_1",
    "sample_2",
    "sample_3",
    "sample_4",
    "species_1",
    "species_2",
    "species_3",
    "species_4",
    "population_1",
    "population_2",
    "population_3",
    "population_4",
    "arrangement_state_1",
    "arrangement_state_2",
    "arrangement_state_3",
    "arrangement_state_4",
    "implied_arrangement_split",
    "same_geography_possible",
    "one_sample_per_species",
    "repeated_species_present",
    "notes",
]


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def forbidden_columns(fields: list[str]) -> list[str]:
    bad = []
    for field in fields:
        lowered = field.lower()
        if any(pattern == lowered or pattern in lowered for pattern in TOPOLOGY_LEAK_PATTERNS):
            bad.append(field)
    return bad


def validate_samples(fields: list[str], rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    missing = [field for field in SAMPLE_FIELDS if field not in fields]
    if missing:
        errors.append(f"sample_manifest.tsv missing required fields: {', '.join(missing)}")
    bad_cols = forbidden_columns(fields)
    if bad_cols:
        errors.append(f"sample_manifest.tsv contains forbidden topology/genealogy columns: {', '.join(bad_cols)}")

    sample_ids = [row.get("sample_id", "").strip() for row in rows]
    counts = Counter(sample_ids)
    duplicates = sorted(s for s, n in counts.items() if s and n > 1)
    if duplicates:
        errors.append(f"duplicate sample_id values: {', '.join(duplicates)}")

    for i, row in enumerate(rows, 2):
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            errors.append(f"row {i} has empty sample_id")
        if not row.get("species", "").strip():
            errors.append(f"row {i} ({sample_id or 'unknown sample'}) has missing species")
        state = row.get("2La_state", "").strip()
        if state not in ALLOWED_STATES:
            errors.append(f"row {i} ({sample_id}) has invalid 2La_state: {state!r}")
        include = truthy(row.get("include_strict_test", ""))
        if include and state not in STRICT_STATES:
            errors.append(f"row {i} ({sample_id}) is included for strict test but has state {state!r}")
        if include and not row.get("karyotype_source", "").strip():
            errors.append(f"row {i} ({sample_id}) is included for strict test without karyotype_source")
    return errors


def validate_regions(path: Path) -> list[str]:
    fields, rows = read_tsv(path)
    errors: list[str] = []
    required = {"region_name", "reference_build", "chromosome", "start", "end", "source", "notes"}
    missing = sorted(required - set(fields))
    if missing:
        errors.append(f"region_manifest.tsv missing required fields: {', '.join(missing)}")
    for i, row in enumerate(rows, 2):
        start_raw = row.get("start", "").strip()
        end_raw = row.get("end", "").strip()
        if start_raw:
            try:
                start = int(start_raw)
                if start < 1:
                    errors.append(f"row {i} has start < 1")
            except ValueError:
                errors.append(f"row {i} has non-integer start: {start_raw!r}")
                start = None
        else:
            start = None
        if end_raw:
            try:
                end = int(end_raw)
                if end < 1:
                    errors.append(f"row {i} has end < 1")
            except ValueError:
                errors.append(f"row {i} has non-integer end: {end_raw!r}")
                end = None
        else:
            end = None
        if start is not None and end is not None and start > end:
            errors.append(f"row {i} has start > end")
    return errors


def same_geography(rows: tuple[dict[str, str], ...]) -> bool:
    countries = {r.get("country", "").strip() for r in rows if r.get("country", "").strip()}
    localities = {r.get("locality", "").strip() for r in rows if r.get("locality", "").strip()}
    populations = {r.get("population", "").strip() for r in rows if r.get("population", "").strip()}
    return bool((len(countries) == 1 and countries) or (len(localities) == 1 and localities) or (len(populations) == 1 and populations))


def classify_design(rows: tuple[dict[str, str], ...]) -> str:
    species = [r["species"] for r in rows]
    states_by_species: dict[str, set[str]] = {}
    pop_state_pairs: set[tuple[str, str, str]] = set()
    for row in rows:
        states_by_species.setdefault(row["species"], set()).add(row["2La_state"])
        pop_state_pairs.add((row.get("species", ""), row.get("population", ""), row["2La_state"]))

    one_per_species = len(set(species)) == 4
    has_within_species_replacement = any({"A0_homozygous", "A1_homozygous"} <= states for states in states_by_species.values())
    geo = same_geography(rows)
    labels = []
    if one_per_species:
        labels.append("Design A - cross-species strict quartet")
    if has_within_species_replacement:
        labels.append("Design B - within-species arrangement replacement")
    if geo:
        labels.append("Design C - population/geography-controlled quartet")
    if not labels:
        labels.append("strict 2:2 arrangement quartet")
    return "; ".join(labels)


def generate_quartets(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    eligible = [
        row
        for row in rows
        if truthy(row.get("include_strict_test", "")) and row.get("2La_state") in STRICT_STATES
    ]
    eligible.sort(key=lambda r: (r.get("species", ""), r.get("population", ""), r.get("sample_id", "")))
    out: list[dict[str, object]] = []
    for combo in itertools.combinations(eligible, 4):
        state_counts = Counter(row["2La_state"] for row in combo)
        if state_counts["A0_homozygous"] != 2 or state_counts["A1_homozygous"] != 2:
            continue
        ordered = sorted(combo, key=lambda r: (r["2La_state"], r.get("species", ""), r.get("sample_id", "")))
        a0 = [r["sample_id"] for r in ordered if r["2La_state"] == "A0_homozygous"]
        a1 = [r["sample_id"] for r in ordered if r["2La_state"] == "A1_homozygous"]
        species = [r.get("species", "") for r in ordered]
        quartet_id = f"Q{len(out) + 1:06d}"
        out.append(
            {
                "quartet_id": quartet_id,
                "design_class": classify_design(tuple(ordered)),
                "sample_1": ordered[0]["sample_id"],
                "sample_2": ordered[1]["sample_id"],
                "sample_3": ordered[2]["sample_id"],
                "sample_4": ordered[3]["sample_id"],
                "species_1": ordered[0].get("species", ""),
                "species_2": ordered[1].get("species", ""),
                "species_3": ordered[2].get("species", ""),
                "species_4": ordered[3].get("species", ""),
                "population_1": ordered[0].get("population", ""),
                "population_2": ordered[1].get("population", ""),
                "population_3": ordered[2].get("population", ""),
                "population_4": ordered[3].get("population", ""),
                "arrangement_state_1": ordered[0]["2La_state"],
                "arrangement_state_2": ordered[1]["2La_state"],
                "arrangement_state_3": ordered[2]["2La_state"],
                "arrangement_state_4": ordered[3]["2La_state"],
                "implied_arrangement_split": f"{a0[0]},{a0[1]}|{a1[0]},{a1[1]}",
                "same_geography_possible": same_geography(tuple(ordered)),
                "one_sample_per_species": len(set(species)) == 4,
                "repeated_species_present": len(set(species)) < 4,
                "notes": "Frozen from 2La_state only; no topology, tree, or quartet-support columns consulted.",
            }
        )
    return out


def aggregate_availability(path: Path = AVAILABILITY) -> tuple[Counter[str], Counter[str], int]:
    if not path.exists():
        return Counter(), Counter(), 0
    _, rows = read_tsv(path)
    by_species: Counter[str] = Counter()
    by_state: Counter[str] = Counter()
    total = 0
    for row in rows:
        count = int(row["sample_count"])
        total += count
        by_species[row["species"]] += count
        state = row.get("2La_state_if_fixed", "") or "unknown"
        by_state[state] += count
    return by_species, by_state, total


def summarize(rows: list[dict[str, str]], quartets: list[dict[str, object]]) -> str:
    species_counts = Counter(row.get("species", "missing") or "missing" for row in rows)
    state_counts = Counter(row.get("2La_state", "missing") or "missing" for row in rows)
    aggregate_species, aggregate_states, aggregate_total = aggregate_availability()
    usable_sequence = sum(1 for row in rows if truthy(row.get("sequence_data_available", "")))
    usable_phased = sum(1 for row in rows if truthy(row.get("phased_data_available", "")))
    one_species = sum(1 for row in quartets if truthy(str(row.get("one_sample_per_species", ""))))
    geo = sum(1 for row in quartets if truthy(str(row.get("same_geography_possible", ""))))
    design_b = sum(1 for row in quartets if "Design B" in str(row.get("design_class", "")))

    lines = [
        "# Stage 0 Anopheles 2La data audit",
        "",
        "Generated by `python3 empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests`.",
        "",
        "## Scope",
        "",
        "Stage 0 freezes only data/sample inventory and structural metadata-derived candidate quartet predictions. It does not read local trees, gene trees, topology labels, quartet support, or Fontaine topology results.",
        "",
        "## Source data located",
        "",
        "- Local repository search found no pre-existing Anopheles files before this scaffold.",
        "- Dryad Fontaine et al. 2015 dataset DOI `10.5061/dryad.f4114` was located; listed files total 2.17 GB, so no large archives were downloaded.",
        "- MalariaGEN Ag3.10 was located as a current distribution containing the `fontaine-2015-rebuild` sample set with 72 samples.",
        "- MalariaGEN Ag3 metadata/karyotype access was identified via `malariagen_data.Ag3().sample_metadata()` and `ag3.karyotype(\"2La\")`.",
        "- `ag3_10_fontaine_rebuild_availability.tsv` records aggregate public sample counts from the Ag3.10 documentation; it is not used to emit sample-level quartets.",
        "",
        "## Reference and regions",
        "",
        "- Reference build recorded for Stage 0: MalariaGEN Ag3 uses AgamP4-aligned data; exact 2La breakpoints trace back to PEST/AgamP3/anoGam3-style coordinates and require equivalence confirmation.",
        "- Chromosome/contig: `2L`.",
        "- 2La interval recorded provisionally: `2L:20524058-42165532`, from White et al. 2007/Fontaine et al. 2015 coordinate lineage and reuse in Ag1000G examples.",
        "- Left/right flanks are represented as configurable candidates only. No final flank length is frozen.",
        "",
        "## Sample manifest status",
        "",
        f"- Sample rows with sample-level provenance: {len(rows)}.",
        f"- Samples with sequence data available: {usable_sequence}.",
        f"- Samples with phased data available: {usable_phased}.",
        "",
        "### Samples by species",
        "",
    ]
    if species_counts:
        for species, count in sorted(species_counts.items()):
            lines.append(f"- {species}: {count}")
    else:
        lines.append("- none: sample-level manifest not yet populated")

    lines.extend(["", "### Samples by 2La state", ""])
    for state in ["A0_homozygous", "A1_homozygous", "heterokaryotype", "unknown"]:
        lines.append(f"- {state}: {state_counts.get(state, 0)}")

    lines.extend(["", "### Aggregate Ag3.10 fontaine-2015-rebuild availability", ""])
    lines.append(f"- Public aggregate sample count: {aggregate_total}.")
    if aggregate_species:
        for species, count in sorted(aggregate_species.items()):
            lines.append(f"- {species}: {count}")
    lines.append("")
    lines.append("Aggregate arrangement-state availability, using only species-level fixed-state literature for non-polymorphic taxa:")
    for state in ["A0_homozygous", "A1_homozygous", "heterokaryotype", "unknown"]:
        lines.append(f"- {state}: {aggregate_states.get(state, 0)}")
    lines.append("")
    lines.append("These aggregate counts suggest cross-species A0/A1 contrasts should be available, but they are not frozen as theorem-test quartets until sample IDs and per-sample provenance are present.")

    lines.extend(
        [
            "",
            "## Candidate strict quartets",
            "",
            f"- Strict 2:2 A0/A1 quartets: {len(quartets)}.",
            f"- One-sample-per-species strict quartets: {one_species}.",
            f"- Geography/population-controlled possibilities: {geo}.",
            f"- Within-species A0-vs-A1 replacement possibilities: {design_b}.",
            "",
            "## Blockers before Stage 1",
            "",
            "- Retrieve sample-level metadata and 2La karyotype calls for the proof-of-principle dataset.",
            "- Confirm exact 2La breakpoint coordinates directly from Fontaine et al. table S11 or an equivalent machine-readable source, not just rounded text summaries.",
            "- Confirm whether the Fontaine rebuild and Ag3 karyotypes use identical AgamP4 coordinates and chromosome naming.",
            "- Do not inspect topology/genealogy outputs until `candidate_strict_quartets.tsv` is populated and frozen.",
        ]
    )
    return "\n".join(lines) + "\n"


def run(sample_manifest: Path, region_manifest: Path, quartets_out: Path, report_out: Path) -> int:
    fields, rows = read_tsv(sample_manifest)
    errors = []
    errors.extend(validate_samples(fields, rows))
    errors.extend(validate_regions(region_manifest))
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2
    quartets = generate_quartets(rows)
    write_tsv(quartets_out, quartets, QUARTET_FIELDS)
    report_out.parent.mkdir(parents=True, exist_ok=True)
    report_out.write_text(summarize(rows, quartets))
    return 0


class AuditTests(unittest.TestCase):
    def test_allowed_states(self) -> None:
        fields = SAMPLE_FIELDS[:]
        rows = [
            {"sample_id": "s1", "species": "gambiae", "2La_state": "A0_homozygous", "include_strict_test": "true", "karyotype_source": "test"},
            {"sample_id": "s2", "species": "gambiae", "2La_state": "bad", "include_strict_test": "false"},
        ]
        errors = validate_samples(fields, rows)
        self.assertTrue(any("invalid 2La_state" in error for error in errors))

    def test_duplicate_sample_ids(self) -> None:
        rows = [
            {"sample_id": "s1", "species": "gambiae", "2La_state": "unknown"},
            {"sample_id": "s1", "species": "coluzzii", "2La_state": "unknown"},
        ]
        errors = validate_samples(SAMPLE_FIELDS[:], rows)
        self.assertTrue(any("duplicate sample_id" in error for error in errors))

    def test_missing_species(self) -> None:
        rows = [{"sample_id": "s1", "species": "", "2La_state": "unknown"}]
        errors = validate_samples(SAMPLE_FIELDS[:], rows)
        self.assertTrue(any("missing species" in error for error in errors))

    def test_forbidden_topology_columns(self) -> None:
        fields = SAMPLE_FIELDS[:] + ["local_tree_topology"]
        errors = validate_samples(fields, [])
        self.assertTrue(any("forbidden topology" in error for error in errors))

    def test_invalid_coordinates(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "regions.tsv"
            write_tsv(
                path,
                [{"region_name": "bad", "reference_build": "AgamP4", "chromosome": "2L", "start": 20, "end": 10, "source": "test", "notes": ""}],
                ["region_name", "reference_build", "chromosome", "start", "end", "source", "notes"],
            )
            errors = validate_regions(path)
        self.assertTrue(any("start > end" in error for error in errors))

    def test_strict_quartet_classification(self) -> None:
        rows = [
            {"sample_id": "a", "species": "sp1", "population": "p", "country": "x", "locality": "l", "2La_state": "A0_homozygous", "include_strict_test": "true"},
            {"sample_id": "b", "species": "sp2", "population": "p", "country": "x", "locality": "l", "2La_state": "A0_homozygous", "include_strict_test": "true"},
            {"sample_id": "c", "species": "sp3", "population": "p", "country": "x", "locality": "l", "2La_state": "A1_homozygous", "include_strict_test": "true"},
            {"sample_id": "d", "species": "sp4", "population": "p", "country": "x", "locality": "l", "2La_state": "A1_homozygous", "include_strict_test": "true"},
        ]
        quartets = generate_quartets(rows)
        self.assertEqual(len(quartets), 1)
        self.assertIn("Design A", quartets[0]["design_class"])
        self.assertIn("Design C", quartets[0]["design_class"])
        self.assertEqual(quartets[0]["implied_arrangement_split"], "a,b|c,d")

    def test_impossible_2_to_2_not_emitted(self) -> None:
        rows = [
            {"sample_id": "a", "species": "sp1", "2La_state": "A0_homozygous", "include_strict_test": "true"},
            {"sample_id": "b", "species": "sp2", "2La_state": "A0_homozygous", "include_strict_test": "true"},
            {"sample_id": "c", "species": "sp3", "2La_state": "A0_homozygous", "include_strict_test": "true"},
            {"sample_id": "d", "species": "sp4", "2La_state": "A1_homozygous", "include_strict_test": "true"},
        ]
        self.assertEqual(generate_quartets(rows), [])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate Anopheles 2La Stage-0 manifests and generate structural-only strict 2:2 candidate quartets."
    )
    parser.add_argument("--sample-manifest", type=Path, default=SAMPLE_MANIFEST, help="Input sample manifest TSV.")
    parser.add_argument("--region-manifest", type=Path, default=REGION_MANIFEST, help="Input region manifest TSV.")
    parser.add_argument("--quartets-out", type=Path, default=QUARTETS_OUT, help="Output candidate strict quartet TSV.")
    parser.add_argument("--report-out", type=Path, default=REPORT_OUT, help="Output Stage-0 markdown report.")
    parser.add_argument("--run-tests", action="store_true", help="Run embedded validation tests before generating outputs.")
    args = parser.parse_args()

    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(AuditTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            return 1
    return run(args.sample_manifest, args.region_manifest, args.quartets_out, args.report_out)


if __name__ == "__main__":
    raise SystemExit(main())
