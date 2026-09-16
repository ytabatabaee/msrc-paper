#!/usr/bin/env python3
"""Fetch Anopheles 2La sample metadata and karyotypes via the MalariaGEN API.

Stage 1A is metadata-only. This script must not call haplotypes, SNP calls,
distance trees, NJ trees, or any observed topology API.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import importlib.metadata
import inspect
import json
import math
import os
import re
import tempfile
import sys
import unittest
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "anopheles_2la"
RESULTS = REPO_ROOT / "empirical" / "anopheles_2la" / "results"
SAMPLE_MANIFEST = DATA_ROOT / "metadata" / "sample_manifest.tsv"
RECONCILIATION = RESULTS / "stage1a_sample_reconciliation.tsv"
PROVENANCE_JSON = DATA_ROOT / "metadata" / "stage1a_malariagen_api_provenance.json"
PROVENANCE_DIR = RESULTS / "provenance"
RAW_KARYOTYPES = DATA_ROOT / "metadata" / "stage1a_raw_2la_karyotypes.tsv"
EXPECTED = DATA_ROOT / "metadata" / "ag3_10_fontaine_rebuild_availability.tsv"
AUDIT_SCRIPT = REPO_ROOT / "empirical" / "anopheles_2la" / "scripts" / "00_audit_samples.py"

TARGET_RELEASE = "3.10"
TARGET_SAMPLE_SET_SUBSTRING = "fontaine"
EXPECTED_TOTAL = 72

FORBIDDEN_API_NAMES = {
    "haplotypes",
    "snp_calls",
    "njt",
    "plot_njt",
    "biallelic_snp_calls",
    "site_allele_frequencies",
    "plot_haplotype_clustering",
    "haplotype_pairwise_distances",
}

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

RECON_FIELDS = ["taxon", "expected_stage0_count", "retrieved_count", "difference", "status", "notes"]

STATUS_AUTH = "AUTHENTICATION_REQUIRED"
STATUS_API = "API_UNAVAILABLE"
STATUS_MISSING = "SAMPLE_SET_NOT_FOUND"
STATUS_SUCCESS = "SUCCESS"


class Stage1AError(RuntimeError):
    status = STATUS_API


class AuthenticationRequired(Stage1AError):
    status = STATUS_AUTH


class SampleSetNotFound(Stage1AError):
    status = STATUS_MISSING


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def atomic_replace_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=path.parent, newline="") as handle:
        handle.write(text)
        tmp = Path(handle.name)
    os.replace(tmp, path)


def atomic_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with src.open("rb") as in_handle, tempfile.NamedTemporaryFile("wb", delete=False, dir=dest.parent) as out_handle:
        for chunk in iter(lambda: in_handle.read(1 << 20), b""):
            out_handle.write(chunk)
        tmp = Path(out_handle.name)
    os.replace(tmp, dest)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def timestamp_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def sanitize_text(text: str) -> str:
    text = re.sub(r"Bearer\s+[A-Za-z0-9._~+/=-]+", "Bearer <redacted>", text, flags=re.I)
    text = re.sub(r"access[_ -]?token[=:]\s*[^,\s]+", "access_token=<redacted>", text, flags=re.I)
    text = re.sub(r"refresh[_ -]?token[=:]\s*[^,\s]+", "refresh_token=<redacted>", text, flags=re.I)
    text = re.sub(r"token=([^&\s]+)", "token=<redacted>", text, flags=re.I)
    text = re.sub(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", "<redacted-account>", text)
    text = re.sub(r"(/[^\s,;:]+)+/(?:[^/\s,;:]+\.json)", "<redacted-credential-path>", text)
    return text


def sanitize_error(exc: BaseException) -> str:
    return sanitize_text(str(exc))


def is_auth_error(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    markers = (
        "401",
        "403",
        "unauthorized",
        "forbidden",
        "permission",
        "anonymous caller",
        "credentials",
        "authentication",
        "storage.objects.get access",
    )
    return any(marker in text for marker in markers)


def print_auth_message() -> None:
    print(
        "MalariaGEN GCS authentication is required.\n"
        "Outside Google Colab, configure Google Application Default Credentials.\n"
        "Recommended command:\n"
        "    gcloud auth application-default login",
        file=sys.stderr,
    )


def write_attempt_provenance(provenance: dict[str, object]) -> Path:
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)
    status = str(provenance.get("retrieval_status", "UNKNOWN")).lower()
    path = PROVENANCE_DIR / f"stage1a_{timestamp_utc()}_{status}.json"
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str) + "\n")
    return path


def load_audit_module() -> Any:
    spec = importlib.util.spec_from_file_location("anopheles_stage0_audit", AUDIT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load audit script from {AUDIT_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install_topology_api_guard(ag3: Any) -> None:
    def blocked(*args: Any, **kwargs: Any) -> None:
        raise RuntimeError("Stage 1A topology/sequence API embargo violated")

    for name in FORBIDDEN_API_NAMES:
        if hasattr(ag3, name):
            setattr(ag3, name, blocked)


def package_version() -> str:
    try:
        return importlib.metadata.version("malariagen-data")
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def normalise_taxon(value: str) -> str:
    value = (value or "").strip()
    mapping = {
        "Anopheles arabiensis": "arabiensis",
        "Anopheles coluzzii": "coluzzii",
        "Anopheles gambiae": "gambiae",
        "Anopheles melas": "melas",
        "Anopheles merus": "merus",
        "Anopheles quadriannulatus": "quadriannulatus",
    }
    return mapping.get(value, value)


def state_from_karyotype(raw: str) -> tuple[str, str]:
    value = (raw or "").strip()
    compact = value.replace(" ", "").replace("^", "")
    if not value or value.lower() in {"nan", "none", ".", "unknown", "unassigned"}:
        return "unknown", "unresolved"
    standard_tokens = {"2L+a/2L+a", "2L+a|2L+a", "+/+", "standard", "standard_hom"}
    inverted_tokens = {"2La/2La", "2La|2La", "inverted", "inverted_hom"}
    het_tokens = {"2L+a/2La", "2La/2L+a", "+/2La", "2La/+", "heterozygote", "heterokaryotype"}
    if compact in standard_tokens or compact.lower() in standard_tokens:
        return "A0_homozygous", "direct_sample_karyotype"
    if compact in inverted_tokens or compact.lower() in inverted_tokens:
        return "A1_homozygous", "direct_sample_karyotype"
    if compact in het_tokens or compact.lower() in het_tokens:
        return "heterokaryotype", "direct_sample_karyotype"
    return "unknown", "unresolved"


def species_fixed_state(taxon: str) -> str:
    fixed = {
        "melas": "A0_homozygous",
        "quadriannulatus": "A0_homozygous",
        "arabiensis": "A1_homozygous",
        "merus": "A1_homozygous",
    }
    return fixed.get(taxon, "")


def find_column(columns: list[str], candidates: list[str]) -> str:
    by_lower = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return ""


def discover_fontaine_sample_set(ag3: Any, release: str, substring: str) -> tuple[str, Any]:
    sample_sets = ag3.sample_sets(release=release)
    id_col = find_column(list(sample_sets.columns), ["sample_set", "id"])
    if not id_col:
        raise RuntimeError(f"Could not find sample-set identifier column in {list(sample_sets.columns)!r}")
    matches = [
        str(value)
        for value in sample_sets[id_col].tolist()
        if substring.lower() in str(value).lower()
    ]
    if len(matches) != 1:
        raise SampleSetNotFound(f"Expected exactly one sample set containing {substring!r}; found {matches!r}")
    return matches[0], sample_sets


def dataframe_to_tsv_rows(df: Any) -> list[dict[str, str]]:
    rows = [{str(k): fmt(v) for k, v in row.items()} for row in df.to_dict(orient="records")]
    if rows and "sample_id" in rows[0]:
        rows.sort(key=lambda row: row.get("sample_id", ""))
    return rows


def build_sample_rows(metadata: Any, karyotypes: Any, sample_set: str, release: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    meta = metadata.copy()
    kdf = karyotypes.copy()
    if "sample_id" not in meta.columns:
        raise RuntimeError("sample_metadata() returned no sample_id column")
    if "sample_id" not in kdf.columns:
        kdf = kdf.reset_index().rename(columns={"index": "sample_id"})
    karyotype_col = find_column(list(kdf.columns), ["karyotype", "call", "2La_karyotype", "genotype"])
    if not karyotype_col:
        non_id = [c for c in kdf.columns if c != "sample_id"]
        if len(non_id) == 1:
            karyotype_col = non_id[0]
        else:
            raise RuntimeError(f"Could not identify karyotype column in {list(kdf.columns)!r}")

    joined = meta.merge(kdf[["sample_id", karyotype_col]], on="sample_id", how="left")
    columns = list(joined.columns)
    taxon_col = find_column(columns, ["taxon", "aim_species", "species", "taxon_call"])
    country_col = find_column(columns, ["country"])
    location_col = find_column(columns, ["location", "locality"])
    year_col = find_column(columns, ["year"])
    month_col = find_column(columns, ["month"])
    lat_col = find_column(columns, ["latitude", "lat"])
    lon_col = find_column(columns, ["longitude", "lon", "long"])
    sex_col = find_column(columns, ["sex_call", "sex"])
    cohort_col = find_column(columns, ["cohort_admin1_year", "cohort", "sample_set"])

    sample_rows: list[dict[str, object]] = []
    raw_rows = dataframe_to_tsv_rows(kdf)
    for record in joined.to_dict(orient="records"):
        sample_id = fmt(record.get("sample_id", ""))
        taxon = normalise_taxon(fmt(record.get(taxon_col, ""))) if taxon_col else ""
        raw = fmt(record.get(karyotype_col, ""))
        state, basis = state_from_karyotype(raw)
        source = "MalariaGEN Ag3 karyotype('2La')"
        method = "API inversion karyotype call"
        confidence = ""
        include = state in {"A0_homozygous", "A1_homozygous"} and taxon in {"gambiae", "coluzzii", "arabiensis"}
        exclusion = "" if include else "not included: unresolved state or classifier validity not established for this taxon"
        if state == "unknown":
            fixed = species_fixed_state(taxon)
            if fixed:
                state = fixed
                basis = "species_fixed"
                source = "literature species-level fixed 2La arrangement"
                method = "species-level fixed arrangement assignment"
                include = True
                exclusion = ""
        sample_rows.append(
            {
                "sample_id": sample_id,
                "sample_set": sample_set,
                "species": taxon,
                "taxon": taxon,
                "population": fmt(record.get(cohort_col, "")) if cohort_col else "",
                "cohort": fmt(record.get(cohort_col, "")) if cohort_col else "",
                "country": fmt(record.get(country_col, "")) if country_col else "",
                "admin1_iso": fmt(record.get("admin1_iso", "")),
                "admin1_name": fmt(record.get("admin1_name", "")),
                "locality": fmt(record.get(location_col, "")) if location_col else "",
                "year": fmt(record.get(year_col, "")) if year_col else "",
                "month": fmt(record.get(month_col, "")) if month_col else "",
                "latitude_if_public": fmt(record.get(lat_col, "")) if lat_col else "",
                "longitude_if_public": fmt(record.get(lon_col, "")) if lon_col else "",
                "sex": fmt(record.get(sex_col, "")) if sex_col else "",
                "dataset": "MalariaGEN Ag3",
                "dataset_version": f"Ag3.{release}",
                "reference_build": "AgamP4",
                "2La_karyotype": raw,
                "2La_karyotype_raw": raw,
                "2La_state": state,
                "karyotype_source": source,
                "karyotype_method": method,
                "karyotype_confidence": confidence,
                "karyotype_confidence_or_support": confidence,
                "state_basis": basis,
                "phased_data_available": "true",
                "sequence_data_available": "true",
                "include_strict_test": include,
                "exclusion_reason": exclusion,
                "notes": "Stage 1A metadata-only freeze; no local genealogy or sequence-derived topology accessed.",
            }
        )
    sample_rows.sort(key=lambda r: (str(r["species"]), str(r["country"]), str(r["sample_id"])))
    return sample_rows, raw_rows


def validate_sample_rows(sample_rows: list[dict[str, object]]) -> None:
    sample_ids = [str(row.get("sample_id", "")) for row in sample_rows]
    missing = [idx for idx, sample_id in enumerate(sample_ids, start=1) if not sample_id]
    if missing:
        raise RuntimeError(f"Stage 1A validation failed: missing sample_id in rows {missing[:5]}")
    duplicates = sorted(sample_id for sample_id, count in Counter(sample_ids).items() if count > 1)
    if duplicates:
        raise RuntimeError(f"Stage 1A validation failed: duplicate sample_id values: {duplicates[:10]}")


def run_fetch_tests() -> None:
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FetchTests))
    if not result.wasSuccessful():
        raise RuntimeError("Stage 1A fetch validation tests failed")


def expected_counts(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in read_tsv(path):
        counts[row["species"]] += int(row["sample_count"])
    return counts


def reconciliation_rows(expected: Counter[str], retrieved: Counter[str], api_status: str, allow_discrepancy: bool = False) -> list[dict[str, object]]:
    taxa = sorted(set(expected) | set(retrieved))
    rows = []
    for taxon in taxa:
        exp = expected.get(taxon, 0)
        got = retrieved.get(taxon, 0)
        diff = got - exp
        status = "ok" if diff == 0 and api_status == STATUS_SUCCESS else "blocked" if api_status != STATUS_SUCCESS else "documented_mismatch" if allow_discrepancy else "mismatch"
        rows.append(
            {
                "taxon": taxon,
                "expected_stage0_count": exp,
                "retrieved_count": got,
                "difference": diff,
                "status": status,
                "notes": "" if status == "ok" else "Mismatch explicitly allowed after authenticated sample-level retrieval." if status == "documented_mismatch" else "Official API retrieval did not complete; Stage-0 aggregate count retained for expectation only." if status == "blocked" else "Authenticated retrieval count differs from Stage-0 aggregate expectation; rerun with --allow-count-discrepancy only after documenting the source discrepancy.",
            }
        )
    if not taxa:
        rows.append({"taxon": "all", "expected_stage0_count": EXPECTED_TOTAL, "retrieved_count": 0, "difference": -EXPECTED_TOTAL, "status": "blocked", "notes": "No sample-level metadata retrieved."})
    return rows


def validate_reconciliation(rows: list[dict[str, object]], allow_discrepancy: bool) -> None:
    mismatches = [row for row in rows if row["status"] == "mismatch"]
    if mismatches and not allow_discrepancy:
        summary = ", ".join(f"{row['taxon']} expected {row['expected_stage0_count']} got {row['retrieved_count']}" for row in mismatches)
        raise RuntimeError(f"Sample reconciliation failed: {summary}")


def run_audit_tests(audit: Any) -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(audit.AuditTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise RuntimeError("Stage 0/1A validation tests failed")


def stage_success_outputs(
    tmpdir: Path,
    audit: Any,
    sample_rows: list[dict[str, object]],
    raw_rows: list[dict[str, str]],
    recon_rows: list[dict[str, object]],
    provenance: dict[str, object],
) -> dict[Path, Path]:
    temp_sample = tmpdir / "sample_manifest.tsv"
    temp_raw = tmpdir / "stage1a_raw_2la_karyotypes.tsv"
    temp_recon = tmpdir / "stage1a_sample_reconciliation.tsv"
    temp_prov = tmpdir / "stage1a_malariagen_api_provenance.json"
    temp_candidate = tmpdir / "candidate_strict_quartets.tsv"
    temp_frozen = tmpdir / "frozen_strict_quartets_stage1a.tsv"
    temp_checksum = tmpdir / "frozen_strict_quartets_stage1a.sha256"
    temp_summary_tsv = tmpdir / "stage1a_quartet_design_summary.tsv"
    temp_summary_md = tmpdir / "stage1a_quartet_design_summary.md"
    temp_report = tmpdir / "stage1a_metadata_karyotype_freeze.md"
    temp_stage0 = tmpdir / "stage0_data_audit.md"

    write_tsv(temp_sample, sample_rows, SAMPLE_FIELDS)
    raw_fields = sorted({field for row in raw_rows for field in row}) if raw_rows else ["sample_id"]
    write_tsv(temp_raw, raw_rows, raw_fields)
    write_tsv(temp_recon, recon_rows, RECON_FIELDS)

    fields, rows = audit.read_tsv(temp_sample)
    errors = []
    errors.extend(audit.validate_samples(fields, rows))
    errors.extend(audit.validate_regions(audit.REGION_MANIFEST))
    if errors:
        raise RuntimeError("Stage 1A validation failed: " + "; ".join(errors))
    quartets = audit.generate_quartets(rows)
    audit.write_tsv(temp_candidate, quartets, audit.QUARTET_FIELDS)
    temp_stage0.write_text(audit.summarize(rows, quartets))
    audit.write_tsv(temp_frozen, quartets, audit.QUARTET_FIELDS)
    checksum = audit.sha256(temp_frozen)
    temp_checksum.write_text(f"{checksum}  {temp_frozen.name}\n")
    audit.write_tsv(temp_summary_tsv, audit.design_summary_rows(rows, quartets), audit.DESIGN_SUMMARY_FIELDS)

    provenance["retrieval_status"] = STATUS_SUCCESS
    provenance["n_candidate_strict_quartets"] = len(quartets)
    provenance["frozen_strict_quartets_sha256"] = checksum
    temp_prov.write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str) + "\n")
    canonical_before = audit.API_PROVENANCE
    audit.API_PROVENANCE = temp_prov
    try:
        audit.write_design_md(temp_summary_md, rows, quartets, checksum)
        audit.write_stage1a_report(temp_report, rows, quartets, checksum)
    finally:
        audit.API_PROVENANCE = canonical_before

    return {
        temp_sample: SAMPLE_MANIFEST,
        temp_raw: RAW_KARYOTYPES,
        temp_recon: RECONCILIATION,
        temp_prov: PROVENANCE_JSON,
        temp_candidate: audit.QUARTETS_OUT,
        temp_stage0: audit.REPORT_OUT,
        temp_frozen: audit.FROZEN_QUARTETS,
        temp_checksum: audit.FROZEN_CHECKSUM,
        temp_summary_tsv: audit.DESIGN_SUMMARY_TSV,
        temp_summary_md: audit.DESIGN_SUMMARY_MD,
        temp_report: audit.STAGE1A_REPORT,
    }


def run(args: argparse.Namespace) -> int:
    provenance: dict[str, object] = {
        "script": str(Path(__file__).relative_to(REPO_ROOT)),
        "command": sanitize_text(" ".join(sys.argv)),
        "attempt_timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "malariagen_data_package_version": package_version(),
        "release_requested": args.release,
        "sample_set_query": args.sample_set_substring,
        "documented_sample_set_identifier": "fontaine-2015-rebuild",
        "retrieval_status": "STARTED",
        "topology_data_embargo": "haplotypes/snp_calls/njt/topology methods not called",
    }
    expected = expected_counts(EXPECTED)
    try:
        import google.auth
        from google.auth.credentials import AnonymousCredentials
        import malariagen_data

        if args.anonymous_credentials:
            google.auth.default = lambda scopes=None, **kwargs: (AnonymousCredentials(), None)

        ag3_kwargs = {"check_location": False, "show_progress": False}
        if args.url:
            ag3_kwargs["url"] = args.url
        ag3 = malariagen_data.Ag3(**ag3_kwargs)
        install_topology_api_guard(ag3)
        sample_set, sample_sets = discover_fontaine_sample_set(ag3, args.release, args.sample_set_substring)
        provenance["sample_set_identifier"] = sample_set
        provenance["available_sample_sets_columns"] = list(sample_sets.columns)
        provenance["available_sample_sets_sha256"] = hashlib.sha256(sample_sets.to_csv(index=False).encode()).hexdigest()
        metadata = ag3.sample_metadata(sample_sets=sample_set)
        karyotypes = ag3.karyotype("2La", sample_sets=sample_set)
        sample_rows, raw_rows = build_sample_rows(metadata, karyotypes, sample_set, args.release)
        validate_sample_rows(sample_rows)
        retrieved = Counter(row["species"] for row in sample_rows)
        recon_rows = reconciliation_rows(expected, retrieved, STATUS_SUCCESS, args.allow_count_discrepancy)
        validate_reconciliation(recon_rows, args.allow_count_discrepancy)
        run_fetch_tests()
        provenance.update(
            {
                "retrieval_status": STATUS_SUCCESS,
                "n_samples_retrieved": len(sample_rows),
                "raw_karyotype_value_counts": Counter(row["2La_karyotype_raw"] for row in sample_rows),
                "state_counts": Counter(row["2La_state"] for row in sample_rows),
                "species_counts": retrieved,
                "api_methods_used": ["sample_sets", "sample_metadata", "karyotype"],
                "api_method_signatures": {
                    name: str(inspect.signature(getattr(ag3, name)))
                    for name in ["sample_sets", "sample_metadata", "karyotype"]
                },
            }
        )
        audit = load_audit_module()
        run_audit_tests(audit)
        with tempfile.TemporaryDirectory(prefix="stage1a_fetch_", dir=str(RESULTS)) as tmp:
            outputs = stage_success_outputs(Path(tmp), audit, sample_rows, raw_rows, recon_rows, provenance)
            temp_sample = next(src for src, dst in outputs.items() if dst == SAMPLE_MANIFEST)
            temp_prov = next(src for src, dst in outputs.items() if dst == PROVENANCE_JSON)
            provenance["sample_manifest"] = str(SAMPLE_MANIFEST.relative_to(REPO_ROOT))
            provenance["sample_manifest_sha256"] = sha256(temp_sample)
            temp_prov.write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str) + "\n")
            for src, dest in outputs.items():
                atomic_copy(src, dest)
        write_attempt_provenance(provenance)
        print(f"Stage 1A retrieval status: {STATUS_SUCCESS}")
        exit_code = 0
    except Exception as exc:
        status = getattr(exc, "status", STATUS_API)
        if status == STATUS_API and is_auth_error(exc):
            status = STATUS_AUTH
        provenance.update(
            {
                "retrieval_status": status,
                "error_type": type(exc).__name__,
                "error": sanitize_error(exc),
                "n_samples_retrieved": 0,
                "expected_total": EXPECTED_TOTAL,
                "note": "No canonical sample manifest or frozen quartet output was replaced after this failed attempt.",
            }
        )
        path = write_attempt_provenance(provenance)
        print(f"Stage 1A retrieval status: {status}", file=sys.stderr)
        print(f"Attempt provenance written to {path.relative_to(REPO_ROOT)}", file=sys.stderr)
        if status == STATUS_AUTH:
            print_auth_message()
        elif status == STATUS_MISSING:
            print("The authenticated API was reachable, but the requested Fontaine sample set was not located.", file=sys.stderr)
        else:
            print("MalariaGEN API access failed before Stage 1A could be safely frozen.", file=sys.stderr)
        exit_code = 3 if status == STATUS_MISSING else 2
    return exit_code


class FetchTests(unittest.TestCase):
    def test_state_normalisation(self) -> None:
        self.assertEqual(state_from_karyotype("2L+a/2L+a")[0], "A0_homozygous")
        self.assertEqual(state_from_karyotype("2La/2La")[0], "A1_homozygous")
        self.assertEqual(state_from_karyotype("2La/2L+a")[0], "heterokaryotype")
        self.assertEqual(state_from_karyotype("")[0], "unknown")

    def test_species_fixed_not_direct(self) -> None:
        self.assertEqual(species_fixed_state("quadriannulatus"), "A0_homozygous")
        self.assertEqual(species_fixed_state("merus"), "A1_homozygous")
        self.assertEqual(species_fixed_state("gambiae"), "")

    def test_forbidden_api_guard(self) -> None:
        class Dummy:
            def haplotypes(self) -> None:
                return None

            def sample_metadata(self) -> None:
                return None

        dummy = Dummy()
        install_topology_api_guard(dummy)
        with self.assertRaises(RuntimeError):
            dummy.haplotypes()
        self.assertIsNone(dummy.sample_metadata())

    def test_auth_error_detection(self) -> None:
        self.assertTrue(is_auth_error(RuntimeError("Anonymous caller denied, 401")))
        self.assertTrue(is_auth_error(RuntimeError("storage.objects.get access forbidden, 403")))
        self.assertFalse(is_auth_error(RuntimeError("sample set missing")))


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch MalariaGEN Ag3 Fontaine rebuild metadata and 2La karyotypes for Stage 1A.")
    parser.add_argument("--release", default=TARGET_RELEASE, help="Ag3 release number to inspect, default 3.10.")
    parser.add_argument("--sample-set-substring", default=TARGET_SAMPLE_SET_SUBSTRING, help="Substring used to discover the Fontaine rebuild sample set.")
    parser.add_argument("--url", default="", help="Optional MalariaGEN API storage URL. Leave empty to use the package default.")
    parser.add_argument("--anonymous-credentials", action="store_true", help="Use anonymous Google credentials for public GCS access when supported.")
    parser.add_argument("--allow-count-discrepancy", action="store_true", help="Permit an authenticated sample-count mismatch after documenting it in reconciliation output.")
    parser.add_argument("--run-tests", action="store_true", help="Run embedded tests before fetching.")
    args = parser.parse_args()
    if args.run_tests:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FetchTests))
        if not result.wasSuccessful():
            return 1
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
