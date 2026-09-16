#!/usr/bin/env python3
"""Audit Anopheles 2La sample metadata and freeze structural-only quartet candidates.

This Stage-0 script uses only sample, geography, species, and 2La karyotype
metadata. It must not read local-tree, gene-tree, topology, or quartet-support
outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
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
FROZEN_QUARTETS = DATA_ROOT / "processed" / "frozen_strict_quartets_stage1a.tsv"
FROZEN_CHECKSUM = DATA_ROOT / "processed" / "frozen_strict_quartets_stage1a.sha256"
DESIGN_SUMMARY_TSV = EMPIRICAL_ROOT / "results" / "stage1a_quartet_design_summary.tsv"
DESIGN_SUMMARY_MD = EMPIRICAL_ROOT / "results" / "stage1a_quartet_design_summary.md"
STAGE1A_REPORT = EMPIRICAL_ROOT / "results" / "stage1a_metadata_karyotype_freeze.md"
API_PROVENANCE = DATA_ROOT / "metadata" / "stage1a_malariagen_api_provenance.json"
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
    "sample_set",
    "species",
    "taxon",
    "population",
    "cohort",
    "country",
    "admin1_iso",
    "admin1_name",
    "locality",
    "year",
    "month",
    "latitude_if_public",
    "longitude_if_public",
    "sex",
    "dataset",
    "dataset_version",
    "reference_build",
    "2La_karyotype",
    "2La_karyotype_raw",
    "2La_state",
    "karyotype_source",
    "karyotype_method",
    "karyotype_confidence",
    "karyotype_confidence_or_support",
    "state_basis",
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
    "species_pattern",
    "populations_used",
    "repeated_species",
    "state_evidence_type",
    "geography_matched",
    "recommended_primary",
    "independence_group",
    "same_geography_possible",
    "one_sample_per_species",
    "repeated_species_present",
    "notes",
]
DESIGN_SUMMARY_FIELDS = ["metric", "value", "notes"]


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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
        if include and not row.get("state_basis", "").strip():
            errors.append(f"row {i} ({sample_id}) is included for strict test without state_basis")
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


def evidence_type(rows: tuple[dict[str, str], ...]) -> str:
    bases = sorted({r.get("state_basis", "") or "unrecorded" for r in rows})
    if bases == ["direct_sample_karyotype"]:
        return "direct_sample_karyotype_only"
    if "species_fixed" in bases:
        return "includes_species_fixed"
    return ";".join(bases)


def independence_group(rows: tuple[dict[str, str], ...]) -> str:
    species = ",".join(sorted({r.get("species", "") for r in rows}))
    countries = ",".join(sorted({r.get("country", "") for r in rows if r.get("country", "")}))
    return f"species={species};countries={countries or 'unrecorded'}"


def is_recommended_primary(rows: tuple[dict[str, str], ...]) -> bool:
    return (
        len({r.get("species", "") for r in rows}) == 4
        and all(r.get("2La_state") in STRICT_STATES for r in rows)
        and all((r.get("state_basis") or "") in {"direct_sample_karyotype", "species_fixed"} for r in rows)
    )


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
        populations = [r.get("population", "") for r in ordered if r.get("population", "")]
        evidence = evidence_type(tuple(ordered))
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
                "species_pattern": ";".join(f"{r.get('species', '')}:{r.get('2La_state', '')}" for r in ordered),
                "populations_used": ";".join(sorted(set(populations))),
                "repeated_species": ";".join(sorted(sp for sp, n in Counter(species).items() if n > 1)),
                "state_evidence_type": evidence,
                "geography_matched": same_geography(tuple(ordered)),
                "recommended_primary": is_recommended_primary(tuple(ordered)) and evidence in {"direct_sample_karyotype_only", "includes_species_fixed"},
                "independence_group": independence_group(tuple(ordered)),
                "same_geography_possible": same_geography(tuple(ordered)),
                "one_sample_per_species": len(set(species)) == 4,
                "repeated_species_present": len(set(species)) < 4,
                "notes": "Frozen from 2La_state only; no topology, tree, or quartet-support columns consulted.",
            }
        )
    return out


def design_summary_rows(rows: list[dict[str, str]], quartets: list[dict[str, object]]) -> list[dict[str, object]]:
    state_counts = Counter(row.get("2La_state", "unknown") or "unknown" for row in rows)
    species_counts = Counter(row.get("species", "missing") or "missing" for row in rows)
    basis_counts = Counter(row.get("state_basis", "missing") or "missing" for row in rows if row.get("2La_state") in STRICT_STATES)
    independence_groups = {row.get("independence_group", "") for row in quartets}
    summary = [
        {"metric": "total_samples", "value": len(rows), "notes": "sample_manifest.tsv rows"},
        {"metric": "A0_homozygotes", "value": state_counts.get("A0_homozygous", 0), "notes": ""},
        {"metric": "A1_homozygotes", "value": state_counts.get("A1_homozygous", 0), "notes": ""},
        {"metric": "heterokaryotypes", "value": state_counts.get("heterokaryotype", 0), "notes": ""},
        {"metric": "unknown_or_unresolved", "value": state_counts.get("unknown", 0), "notes": ""},
        {"metric": "strict_2_to_2_candidate_quartets", "value": len(quartets), "notes": ""},
        {"metric": "four_distinct_species_quartets", "value": sum(1 for q in quartets if truthy(str(q.get("one_sample_per_species", "")))), "notes": ""},
        {"metric": "direct_sample_level_only_quartets", "value": sum(1 for q in quartets if q.get("state_evidence_type") == "direct_sample_karyotype_only"), "notes": ""},
        {"metric": "quartets_relying_partly_on_species_fixed_state", "value": sum(1 for q in quartets if q.get("state_evidence_type") == "includes_species_fixed"), "notes": ""},
        {"metric": "geography_matched_quartets", "value": sum(1 for q in quartets if truthy(str(q.get("geography_matched", "")))), "notes": ""},
        {"metric": "recommended_primary_quartets", "value": sum(1 for q in quartets if truthy(str(q.get("recommended_primary", "")))), "notes": ""},
        {"metric": "independence_groups", "value": len(independence_groups), "notes": ""},
        {"metric": "direct_sample_karyotype_states", "value": basis_counts.get("direct_sample_karyotype", 0), "notes": "resolved homokaryotypic samples"},
        {"metric": "species_fixed_states", "value": basis_counts.get("species_fixed", 0), "notes": "resolved homokaryotypic samples"},
    ]
    for species, count in sorted(species_counts.items()):
        summary.append({"metric": f"samples_species_{species}", "value": count, "notes": ""})
    return summary


def write_design_md(path: Path, rows: list[dict[str, str]], quartets: list[dict[str, object]], checksum: str) -> None:
    summary = design_summary_rows(rows, quartets)
    top = sorted(
        quartets,
        key=lambda q: (
            not truthy(str(q.get("recommended_primary", ""))),
            not truthy(str(q.get("one_sample_per_species", ""))),
            str(q.get("state_evidence_type", "")),
            str(q.get("quartet_id", "")),
        ),
    )[:20]
    lines = [
        "# Stage 1A quartet design summary",
        "",
        "Generated from structural/karyotype metadata only. No local tree, topology, q1/q2/q3, QQS/BQS, distance-derived phylogeny, or sequence-derived topology input was read.",
        "",
        "## Metrics",
        "",
    ]
    for row in summary:
        lines.append(f"- {row['metric']}: {row['value']}")
    lines.extend(["", "## Top candidate designs", ""])
    if top:
        for q in top:
            lines.append(f"- {q['quartet_id']}: {q['implied_arrangement_split']} ({q['design_class']}; {q['state_evidence_type']})")
    else:
        lines.append("- none: no sample-level strict 2:2 quartets are available")
    lines.extend(["", "## Frozen table", "", f"- `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`", f"- sha256: `{checksum}`", ""])
    path.write_text("\n".join(lines))


def write_stage1a_report(path: Path, rows: list[dict[str, str]], quartets: list[dict[str, object]], checksum: str) -> None:
    provenance = {}
    if API_PROVENANCE.exists():
        provenance = json.loads(API_PROVENANCE.read_text())
    state_counts = Counter(row.get("2La_state", "unknown") or "unknown" for row in rows)
    species_counts = Counter(row.get("species", "missing") or "missing" for row in rows)
    basis_counts = Counter(row.get("state_basis", "missing") or "missing" for row in rows if row.get("2La_state") in STRICT_STATES)
    unresolved = [row for row in rows if row.get("2La_state") == "unknown" or not truthy(row.get("include_strict_test", ""))]
    top = [q for q in quartets if truthy(str(q.get("recommended_primary", "")))][:10]
    lines = [
        "# Stage 1A metadata/karyotype freeze",
        "",
        "Stage 1A stops after structural/karyotype prediction freezing. It does not infer or inspect local genealogies.",
        "",
        "## API retrieval",
        "",
        f"- MalariaGEN package version: {provenance.get('malariagen_data_package_version', 'not recorded')}",
        f"- Requested release: {provenance.get('release_requested', '3.10')}",
        f"- Sample-set query: {provenance.get('sample_set_query', 'fontaine')}",
        f"- Documented sample set identifier: {provenance.get('documented_sample_set_identifier', 'fontaine-2015-rebuild')}",
        f"- Programmatically discovered sample set identifier: {provenance.get('sample_set_identifier', 'not discovered')}",
        f"- Retrieval status: {provenance.get('retrieval_status', 'not run')}",
        f"- Samples retrieved: {provenance.get('n_samples_retrieved', len(rows))}",
        "",
        "## Reconciliation",
        "",
        "- Expected Stage-0 aggregate total: 72.",
        f"- Retrieved sample-manifest rows: {len(rows)}.",
        "- See `empirical/anopheles_2la/results/stage1a_sample_reconciliation.tsv`.",
        "",
        "## Karyotype provenance",
        "",
        "- Preferred method: `malariagen_data.Ag3().karyotype(\"2La\", sample_sets=<fontaine sample set>)`.",
        "- Species-fixed states are allowed only with `state_basis=species_fixed` and are not mislabeled as direct sample-level calls.",
        "",
        "## State counts",
        "",
    ]
    for state in ["A0_homozygous", "A1_homozygous", "heterokaryotype", "unknown"]:
        lines.append(f"- {state}: {state_counts.get(state, 0)}")
    lines.extend(["", "## Counts by species", ""])
    if species_counts:
        for species, count in sorted(species_counts.items()):
            lines.append(f"- {species}: {count}")
    else:
        lines.append("- none: API retrieval blocked before sample-level metadata were returned")
    lines.extend(
        [
            "",
            "## Evidence basis",
            "",
            f"- direct_sample_karyotype states: {basis_counts.get('direct_sample_karyotype', 0)}",
            f"- species_fixed states: {basis_counts.get('species_fixed', 0)}",
            f"- unresolved/excluded samples: {len(unresolved)}",
            "",
            "## Strict quartet predictions",
            "",
            f"- strict 2:2 candidate quartets: {len(quartets)}",
            f"- four-distinct-species strict quartets: {sum(1 for q in quartets if truthy(str(q.get('one_sample_per_species', ''))))}",
            f"- geography-matched strict quartets: {sum(1 for q in quartets if truthy(str(q.get('geography_matched', ''))))}",
            f"- frozen file: `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`",
            f"- sha256: `{checksum}`",
            "",
            "## Best Design A candidates",
            "",
        ]
    )
    if top:
        for q in top:
            lines.append(f"- {q['quartet_id']}: {q['implied_arrangement_split']} ({q['state_evidence_type']})")
    else:
        lines.append("- none available from authoritative sample-level manifest")
    lines.extend(
        [
            "",
            "## Design B arrangement-replacement contrasts",
            "",
            "- none available until polymorphic gambiae/coluzzii sample-level 2La states are retrieved.",
            "",
            "## Design C geography-matched contrasts",
            "",
            "- none available until sample-level geography and 2La states are retrieved.",
            "",
            "## Proceed to Stage 1B?",
            "",
            "No. Authoritative sample-level metadata/karyotype retrieval remains blocked in this environment, so no clean strict 2:2 predictions have been populated beyond the deterministic empty freeze.",
        ]
    )
    if provenance.get("error"):
        lines.extend(["", "## Retrieval blocker", "", f"- {provenance.get('error_type')}: {provenance.get('error')}"])
    path.write_text("\n".join(lines) + "\n")


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


def run(sample_manifest: Path, region_manifest: Path, quartets_out: Path, report_out: Path, freeze_stage1a: bool = False) -> int:
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
    if freeze_stage1a:
        write_tsv(FROZEN_QUARTETS, quartets, QUARTET_FIELDS)
        checksum = sha256(FROZEN_QUARTETS)
        FROZEN_CHECKSUM.write_text(f"{checksum}  {FROZEN_QUARTETS.name}\n")
        write_tsv(DESIGN_SUMMARY_TSV, design_summary_rows(rows, quartets), DESIGN_SUMMARY_FIELDS)
        write_design_md(DESIGN_SUMMARY_MD, rows, quartets, checksum)
        write_stage1a_report(STAGE1A_REPORT, rows, quartets, checksum)
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

    def test_included_state_requires_provenance(self) -> None:
        rows = [{"sample_id": "s1", "species": "gambiae", "2La_state": "A0_homozygous", "include_strict_test": "true", "karyotype_source": "source", "state_basis": ""}]
        errors = validate_samples(SAMPLE_FIELDS[:], rows)
        self.assertTrue(any("without state_basis" in error for error in errors))

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

    def test_deterministic_ordering_and_checksum(self) -> None:
        rows = [
            {"sample_id": "d", "species": "sp4", "population": "", "country": "", "2La_state": "A1_homozygous", "include_strict_test": "true", "state_basis": "direct_sample_karyotype"},
            {"sample_id": "a", "species": "sp1", "population": "", "country": "", "2La_state": "A0_homozygous", "include_strict_test": "true", "state_basis": "direct_sample_karyotype"},
            {"sample_id": "c", "species": "sp3", "population": "", "country": "", "2La_state": "A1_homozygous", "include_strict_test": "true", "state_basis": "direct_sample_karyotype"},
            {"sample_id": "b", "species": "sp2", "population": "", "country": "", "2La_state": "A0_homozygous", "include_strict_test": "true", "state_basis": "direct_sample_karyotype"},
        ]
        first = generate_quartets(rows)
        second = generate_quartets(list(reversed(rows)))
        self.assertEqual(first, second)
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "q.tsv"
            write_tsv(path, first, QUARTET_FIELDS)
            self.assertEqual(sha256(path), sha256(path))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate Anopheles 2La Stage-0 manifests and generate structural-only strict 2:2 candidate quartets."
    )
    parser.add_argument("--sample-manifest", type=Path, default=SAMPLE_MANIFEST, help="Input sample manifest TSV.")
    parser.add_argument("--region-manifest", type=Path, default=REGION_MANIFEST, help="Input region manifest TSV.")
    parser.add_argument("--quartets-out", type=Path, default=QUARTETS_OUT, help="Output candidate strict quartet TSV.")
    parser.add_argument("--report-out", type=Path, default=REPORT_OUT, help="Output Stage-0 markdown report.")
    parser.add_argument("--freeze-stage1a", action="store_true", help="Also write Stage-1A frozen quartet copy, checksum, and design reports.")
    parser.add_argument("--run-tests", action="store_true", help="Run embedded validation tests before generating outputs.")
    args = parser.parse_args()

    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(AuditTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            return 1
    return run(args.sample_manifest, args.region_manifest, args.quartets_out, args.report_out, args.freeze_stage1a)


if __name__ == "__main__":
    raise SystemExit(main())
