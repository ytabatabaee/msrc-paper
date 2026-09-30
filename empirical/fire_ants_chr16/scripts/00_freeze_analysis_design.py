#!/usr/bin/env python3
"""Freeze the fire-ant chromosome-16 Stage-0 analysis design.

Stage 0 uses only provenance, structural labels, author-designated coordinate
regions, and biologically defined SB/Sb group labels. It must not read TWISST
weights, local trees, Newick files, topology outputs, or quartet-support files.
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
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "fire_ants_chr16"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "fire_ants_chr16"

RAW_UPSTREAM = DATA_ROOT / "raw" / "upstream"
COORDINATES = RAW_UPSTREAM / "Topology weighting" / "results" / "window_coordinates"
SOURCE_MANIFEST = DATA_ROOT / "metadata" / "source_manifest.tsv"
DOWNLOAD_MANIFEST = DATA_ROOT / "metadata" / "download_manifest.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"
FOCAL_GROUP_MANIFEST = DATA_ROOT / "metadata" / "focal_group_manifest.tsv"
UPSTREAM_PROVENANCE = DATA_ROOT / "metadata" / "upstream_repository_provenance.json"
FOCAL_QUARTET = DATA_ROOT / "processed" / "stage0_focal_quartet.tsv"
FOCAL_QUARTET_SHA = DATA_ROOT / "processed" / "stage0_focal_quartet.sha256"
REGION_INVENTORY = EMPIRICAL_ROOT / "results" / "stage0_region_inventory.tsv"
STAGE0_REPORT = EMPIRICAL_ROOT / "results" / "stage0_report.md"

EXPECTED_DIRS = [
    EMPIRICAL_ROOT,
    EMPIRICAL_ROOT / "config",
    EMPIRICAL_ROOT / "figures",
    EMPIRICAL_ROOT / "results",
    EMPIRICAL_ROOT / "scripts",
    EMPIRICAL_ROOT / "tests",
    DATA_ROOT,
    DATA_ROOT / "raw",
    DATA_ROOT / "metadata",
    DATA_ROOT / "intermediate",
    DATA_ROOT / "processed",
]

EXPECTED_REGIONS = ["chr1", "chr16A", "chr16B", "chr16_supergene"]
EXPECTED_COUNTS = {
    "chr1": 117,
    "chr16A": 42,
    "chr16B": 2,
    "chr16_supergene": 52,
}
EXPECTED_GROUPS = {
    "inv_mac_SB": ("invicta/macdonaghi_SB", "invicta/macdonaghi", "SB", "A"),
    "inv_mac_Sb": ("invicta/macdonaghi_Sb", "invicta/macdonaghi", "Sb", "B"),
    "richteri_SB": ("richteri_SB", "richteri", "SB", "C"),
    "richteri_Sb": ("richteri_Sb", "richteri", "Sb", "D"),
}
EXPECTED_SPLITS = {
    "species_split": "AB|CD",
    "haplotype_split": "AC|BD",
    "third_split": "AD|BC",
}
RAW_FILES = [
    RAW_UPSTREAM / "README.md",
    RAW_UPSTREAM / "Topology weighting" / "README.md",
    RAW_UPSTREAM / "Topology weighting" / "supergene" / "README.md",
    COORDINATES,
    RAW_UPSTREAM / "Coalescence-based phylogenies" / "README.md",
]

FORBIDDEN_INPUT_PATTERNS = (
    re.compile(r"(^|[/_.-])weights?([/_.-]|$)", re.IGNORECASE),
    re.compile(r"(^|[/_.-])topolog(?:y|ies)([/_.-]|$)", re.IGNORECASE),
    re.compile(r"(^|/)input[.]tree$", re.IGNORECASE),
    re.compile(r"(^|[/_.-])newick([/_.-]|$)", re.IGNORECASE),
    re.compile(r"[.](?:nwk|tre|tree)$", re.IGNORECASE),
    re.compile(r"(^|[/_.-])quartet[_-]?support([/_.-]|$)", re.IGNORECASE),
    re.compile(r"(^|[/_.-])qqs([/_.-]|$)", re.IGNORECASE),
    re.compile(r"(^|[/_.-])q[123]([/_.-]|$)", re.IGNORECASE),
)

TSV_LINETERMINATOR = "\n"


class Stage0Error(RuntimeError):
    """Raised when Stage-0 validation fails."""


def fmt(value: object) -> str:
    return "" if value is None else str(value)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_stage0_input_allowed(path: Path, purpose: str = "analysis_input") -> None:
    """Reject topology-bearing analysis inputs while allowing documentation text."""

    normalized = path.as_posix()
    if purpose == "documentation":
        return
    if normalized.endswith("/Topology weighting/results/window_coordinates"):
        return
    if path.name == "window_coordinates":
        return
    for pattern in FORBIDDEN_INPUT_PATTERNS:
        if pattern.search(normalized):
            raise Stage0Error(f"Stage 0 may not read topology-bearing analysis input: {path}")


def read_tsv(path: Path, *, purpose: str = "analysis_input") -> tuple[list[str], list[dict[str, str]]]:
    assert_stage0_input_allowed(path, purpose=purpose)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def parse_coord(value: str) -> float:
    return float(value)


def fmt_coord(value: object) -> str:
    number = float(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:.12g}"


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def validate_required_dirs() -> list[str]:
    return [f"missing required directory: {path.relative_to(REPO_ROOT)}" for path in EXPECTED_DIRS if not path.is_dir()]


def validate_manifest_schema() -> list[str]:
    errors: list[str] = []
    schemas = {
        SOURCE_MANIFEST: [
            "source_name",
            "source_type",
            "citation",
            "url_or_doi",
            "data_role",
            "expected_format",
            "expected_size",
            "downloaded",
            "notes",
        ],
        DOWNLOAD_MANIFEST: [
            "source",
            "original_filename",
            "url_or_doi",
            "download_date",
            "size_bytes",
            "sha256",
            "purpose",
            "local_path",
        ],
        REGION_MANIFEST: [
            "region_id",
            "chromosome",
            "source_region_label",
            "role",
            "coordinate_start",
            "coordinate_end",
            "coordinate_definition",
            "source",
            "notes",
        ],
        FOCAL_GROUP_MANIFEST: [
            "canonical_group",
            "upstream_group_label",
            "species_group",
            "haplotype",
            "role",
            "source",
            "notes",
        ],
    }
    for path, required in schemas.items():
        if not path.exists():
            errors.append(f"missing manifest: {path.relative_to(REPO_ROOT)}")
            continue
        fields, _ = read_tsv(path, purpose="metadata")
        missing = [field for field in required if field not in fields]
        if missing:
            errors.append(f"{path.relative_to(REPO_ROOT)} missing fields: {', '.join(missing)}")
    return errors


def validate_upstream_provenance() -> list[str]:
    errors: list[str] = []
    if not UPSTREAM_PROVENANCE.exists():
        return ["missing upstream_repository_provenance.json"]
    data = json.loads(UPSTREAM_PROVENANCE.read_text())
    expected = {
        "repository_url": "https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression",
        "default_branch": "master",
        "head_commit_sha": "bdc7823941680fa560e61b6467af9f77950f64b0",
        "retrieval_date": "2026-09-30",
    }
    for key, value in expected.items():
        require(data.get(key) == value, f"upstream provenance {key} expected {value!r}, observed {data.get(key)!r}", errors)
    require(len(data.get("head_commit_sha", "")) == 40, "upstream HEAD commit SHA must be 40 characters", errors)
    return errors


def validate_downloads() -> list[str]:
    errors: list[str] = []
    fields, rows = read_tsv(DOWNLOAD_MANIFEST, purpose="metadata")
    require("sha256" in fields, "download_manifest.tsv must include sha256", errors)
    by_path = {row.get("local_path", ""): row for row in rows}
    for raw_path in RAW_FILES:
        rel = raw_path.relative_to(REPO_ROOT).as_posix()
        row = by_path.get(rel)
        if row is None:
            errors.append(f"download_manifest.tsv missing raw file: {rel}")
            continue
        if not raw_path.exists():
            errors.append(f"downloaded raw file missing on disk: {rel}")
            continue
        observed_size = raw_path.stat().st_size
        observed_sha = sha256(raw_path)
        require(row.get("size_bytes") == str(observed_size), f"{rel} size mismatch in download_manifest.tsv", errors)
        require(row.get("sha256") == observed_sha, f"{rel} sha256 mismatch in download_manifest.tsv", errors)
    return errors


def coordinate_rows() -> tuple[list[str], list[dict[str, str]]]:
    fields, rows = read_tsv(COORDINATES)
    required = ["chrom", "start", "end", "mid", "region", "window", "gene1", "gene2", "gene3", "gene4"]
    missing = [field for field in required if field not in fields]
    if missing:
        raise Stage0Error(f"window_coordinates missing fields: {', '.join(missing)}")
    return fields, rows


def region_inventory_rows() -> list[dict[str, object]]:
    _, rows = coordinate_rows()
    inventory: dict[str, dict[str, object]] = {}
    for i, row in enumerate(rows, 2):
        region = row["region"]
        if region not in EXPECTED_REGIONS:
            raise Stage0Error(f"unrecognized region label in window_coordinates row {i}: {region}")
        start = parse_coord(row["start"])
        end = parse_coord(row["end"])
        if start > end:
            raise Stage0Error(f"window_coordinates row {i} has start > end")
        item = inventory.setdefault(
            region,
            {
                "region_id": region,
                "chromosome": row["chrom"],
                "window_count": 0,
                "observed_window_span_start": start,
                "observed_window_span_end": end,
                "first_window_id": row["window"],
                "last_window_id": row["window"],
                "coordinate_definition": "author-designated analysis region / observed BUSCO-window span",
            },
        )
        if item["chromosome"] != row["chrom"]:
            raise Stage0Error(f"region {region} appears on multiple chromosomes")
        item["window_count"] = int(item["window_count"]) + 1
        item["observed_window_span_start"] = min(float(item["observed_window_span_start"]), start)
        item["observed_window_span_end"] = max(float(item["observed_window_span_end"]), end)
        item["last_window_id"] = row["window"]

    missing = [region for region in EXPECTED_REGIONS if region not in inventory]
    if missing:
        raise Stage0Error(f"window_coordinates missing expected regions: {', '.join(missing)}")
    for item in inventory.values():
        item["observed_window_span_start"] = fmt_coord(item["observed_window_span_start"])
        item["observed_window_span_end"] = fmt_coord(item["observed_window_span_end"])
    return [inventory[region] for region in EXPECTED_REGIONS]


def validate_region_manifest() -> list[str]:
    errors: list[str] = []
    _, rows = read_tsv(REGION_MANIFEST, purpose="metadata")
    by_region = {row["region_id"]: row for row in rows}
    inventory = {row["region_id"]: row for row in region_inventory_rows()}
    require(set(by_region) == set(EXPECTED_REGIONS), "region_manifest.tsv must contain exactly the expected four regions", errors)
    for region in EXPECTED_REGIONS:
        manifest = by_region.get(region, {})
        inv = inventory[region]
        require(manifest.get("coordinate_start") == fmt_coord(inv["observed_window_span_start"]), f"{region} coordinate_start does not match window_coordinates span", errors)
        require(manifest.get("coordinate_end") == fmt_coord(inv["observed_window_span_end"]), f"{region} coordinate_end does not match window_coordinates span", errors)
        require(
            manifest.get("coordinate_definition") == "author-designated analysis region / observed BUSCO-window span",
            f"{region} coordinate_definition must distinguish BUSCO-window span from breakpoints",
            errors,
        )
    observed_counts = {row["region_id"]: int(row["window_count"]) for row in inventory.values()}
    require(observed_counts == EXPECTED_COUNTS, f"region counts changed: {observed_counts}", errors)
    require(sum(observed_counts.values()) == 213, "region counts must sum to 213 windows", errors)
    return errors


def validate_groups() -> list[str]:
    errors: list[str] = []
    _, rows = read_tsv(FOCAL_GROUP_MANIFEST, purpose="metadata")
    by_group = {row["canonical_group"]: row for row in rows}
    require(set(by_group) == set(EXPECTED_GROUPS), "focal_group_manifest.tsv must contain exactly the four focal groups", errors)
    roles = set()
    for group, (upstream, species_group, haplotype, role) in EXPECTED_GROUPS.items():
        row = by_group.get(group, {})
        require(row.get("upstream_group_label") == upstream, f"{group} upstream label mismatch", errors)
        require(row.get("species_group") == species_group, f"{group} species_group mismatch", errors)
        require(row.get("haplotype") == haplotype, f"{group} haplotype mismatch", errors)
        require(row.get("haplotype") in {"SB", "Sb"}, f"{group} must have exactly one SB or Sb state", errors)
        require(row.get("role") == role, f"{group} role mismatch", errors)
        roles.add(row.get("role", ""))
    require(roles == {"A", "B", "C", "D"}, "focal groups must map to roles A, B, C, D", errors)
    return errors


def validate_quartet() -> list[str]:
    errors: list[str] = []
    fields, rows = read_tsv(FOCAL_QUARTET, purpose="metadata")
    required = [
        "quartet_id",
        "role_A",
        "role_B",
        "role_C",
        "role_D",
        "group_A",
        "group_B",
        "group_C",
        "group_D",
        "resolution_id",
        "split",
        "role_split",
        "group_split",
        "interpretation",
        "primary_future_quantity",
    ]
    missing = [field for field in required if field not in fields]
    require(not missing, f"stage0_focal_quartet.tsv missing fields: {', '.join(missing)}", errors)
    require(len(rows) == 3, "stage0_focal_quartet.tsv must contain exactly three unrooted resolutions", errors)
    quartet_ids = {row.get("quartet_id") for row in rows}
    require(quartet_ids == {"fire_ants_chr16_focal"}, "stage0_focal_quartet.tsv must contain exactly one focal quartet", errors)
    for row in rows:
        groups = [row.get(f"group_{role}") for role in "ABCD"]
        require(len(set(groups)) == 4, "focal quartet must have four distinct groups", errors)
        require(groups == ["inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"], "focal quartet group order must match roles A-D", errors)
    splits = {row.get("resolution_id"): row.get("split") for row in rows}
    require(splits == EXPECTED_SPLITS, f"unexpected quartet splits: {splits}", errors)
    split_values = set(splits.values())
    require(split_values == {"AB|CD", "AC|BD", "AD|BC"}, "quartet resolutions must be the three distinct unrooted splits", errors)
    return errors


def write_region_inventory() -> None:
    fields = [
        "region_id",
        "chromosome",
        "window_count",
        "observed_window_span_start",
        "observed_window_span_end",
        "first_window_id",
        "last_window_id",
        "coordinate_definition",
    ]
    write_tsv(REGION_INVENTORY, region_inventory_rows(), fields)


def write_focal_quartet_sha() -> None:
    FOCAL_QUARTET_SHA.write_text(f"{sha256(FOCAL_QUARTET)}  {FOCAL_QUARTET.relative_to(REPO_ROOT).as_posix()}\n")


def write_stage0_report() -> None:
    provenance = json.loads(UPSTREAM_PROVENANCE.read_text())
    inventory = region_inventory_rows()
    raw_lines = []
    _, downloads = read_tsv(DOWNLOAD_MANIFEST, purpose="metadata")
    for row in downloads:
        raw_lines.append(f"- `{row['local_path']}`: `{row['sha256']}`")

    lines = [
        "# Fire ants chromosome 16 Stage 0 report",
        "",
        "Stage 0 scaffolds the fire-ant chromosome-16 social-supergene analysis and freezes the design before local topology support is reproduced in this repository.",
        "",
        "## Upstream repository",
        "",
        f"- Repository: `{provenance['repository_url']}`",
        f"- Default branch: `{provenance['default_branch']}`",
        f"- HEAD commit: `{provenance['head_commit_sha']}`",
        f"- Retrieval date: `{provenance['retrieval_date']}`",
        "",
        "## Region inventory",
        "",
        "| region | chromosome | windows | observed BUSCO-window span |",
        "| --- | --- | ---: | --- |",
    ]
    for row in inventory:
        lines.append(
            f"| `{row['region_id']}` | `{row['chromosome']}` | {row['window_count']} | "
            f"{row['observed_window_span_start']}-{row['observed_window_span_end']} |"
        )
    lines.extend(
        [
            "",
            "Coordinates are author-designated analysis regions and observed BUSCO-window spans. They are not treated as nucleotide-resolved inversion breakpoints.",
            "",
            "## Frozen focal quartet",
            "",
            "- A: `invicta/macdonaghi_SB` (`inv_mac_SB`)",
            "- B: `invicta/macdonaghi_Sb` (`inv_mac_Sb`)",
            "- C: `richteri_SB` (`richteri_SB`)",
            "- D: `richteri_Sb` (`richteri_Sb`)",
            "",
            "Frozen resolutions:",
            "",
            "- `species_split`: `AB|CD`",
            "- `haplotype_split`: `AC|BD`",
            "- `third_split`: `AD|BC`",
            "",
            f"Stage-0 focal-quartet checksum: `{sha256(FOCAL_QUARTET)}`",
            "",
            "## Raw source-file checksums",
            "",
            *raw_lines,
            "",
            "## Guardrails",
            "",
            "The Stage-0 validator rejects topology-bearing analysis inputs such as TWISST weights, topology-output trees, input trees, Newick files, quartet-support files, QQS files, and precomputed q1/q2/q3 support tables.",
            "",
            "No TWISST weights, local trees, Newick files, ASTRAL/ASTER outputs, or quartet-support values were parsed in Stage 0.",
        ]
    )
    STAGE0_REPORT.write_text("\n".join(lines) + "\n")


def validate_all() -> list[str]:
    errors: list[str] = []
    errors.extend(validate_required_dirs())
    errors.extend(validate_manifest_schema())
    errors.extend(validate_upstream_provenance())
    errors.extend(validate_downloads())
    errors.extend(validate_region_manifest())
    errors.extend(validate_groups())
    errors.extend(validate_quartet())
    return errors


def freeze() -> None:
    write_region_inventory()
    write_focal_quartet_sha()
    write_stage0_report()
    errors = validate_all()
    if errors:
        raise Stage0Error("\n".join(errors))


class Stage0Tests(unittest.TestCase):
    def test_required_directories_exist(self) -> None:
        self.assertEqual(validate_required_dirs(), [])

    def test_source_manifest_schema(self) -> None:
        fields, _ = read_tsv(SOURCE_MANIFEST, purpose="metadata")
        self.assertIn("source_name", fields)
        self.assertIn("citation", fields)
        self.assertIn("downloaded", fields)

    def test_upstream_head_provenance_captured(self) -> None:
        self.assertEqual(validate_upstream_provenance(), [])

    def test_copied_raw_files_have_sha256_values(self) -> None:
        self.assertEqual(validate_downloads(), [])

    def test_coordinate_region_labels_are_recognized(self) -> None:
        labels = {row["region_id"] for row in region_inventory_rows()}
        self.assertEqual(labels, set(EXPECTED_REGIONS))

    def test_region_counts_sum_correctly(self) -> None:
        counts = {row["region_id"]: int(row["window_count"]) for row in region_inventory_rows()}
        self.assertEqual(counts, EXPECTED_COUNTS)
        self.assertEqual(sum(counts.values()), 213)

    def test_focal_quartet_has_four_distinct_groups(self) -> None:
        _, rows = read_tsv(FOCAL_QUARTET, purpose="metadata")
        for row in rows:
            self.assertEqual(len({row[f"group_{role}"] for role in "ABCD"}), 4)

    def test_each_group_has_one_sb_state(self) -> None:
        _, rows = read_tsv(FOCAL_GROUP_MANIFEST, purpose="metadata")
        self.assertEqual(Counter(row["haplotype"] for row in rows), Counter({"SB": 2, "Sb": 2}))
        for row in rows:
            self.assertIn(row["haplotype"], {"SB", "Sb"})

    def test_three_distinct_unrooted_quartet_resolutions(self) -> None:
        _, rows = read_tsv(FOCAL_QUARTET, purpose="metadata")
        splits = {row["resolution_id"]: row["split"] for row in rows}
        self.assertEqual(splits, EXPECTED_SPLITS)
        self.assertEqual(set(splits.values()), {"AB|CD", "AC|BD", "AD|BC"})

    def test_stage0_script_rejects_topology_bearing_files(self) -> None:
        bad_names = [
            "weights_output.csv",
            "topologies_output.trees",
            "input.tree",
            "Trees in newick format/example.nwk",
            "quartet_support.tsv",
            "QQS.tsv",
            "q1.tsv",
            "q2.tsv",
            "q3.tsv",
        ]
        for name in bad_names:
            with self.subTest(name=name):
                with self.assertRaises(Stage0Error):
                    assert_stage0_input_allowed(Path(name))
        assert_stage0_input_allowed(Path("Topology weighting/README.md"), purpose="documentation")

    def test_repeated_runs_are_deterministic(self) -> None:
        freeze()
        before = {
            REGION_INVENTORY: sha256(REGION_INVENTORY),
            FOCAL_QUARTET_SHA: sha256(FOCAL_QUARTET_SHA),
            STAGE0_REPORT: sha256(STAGE0_REPORT),
        }
        freeze()
        after = {
            REGION_INVENTORY: sha256(REGION_INVENTORY),
            FOCAL_QUARTET_SHA: sha256(FOCAL_QUARTET_SHA),
            STAGE0_REPORT: sha256(STAGE0_REPORT),
        }
        self.assertEqual(before, after)

    def test_stage0_checksum_is_unchanged_on_rerun(self) -> None:
        freeze()
        first = FOCAL_QUARTET_SHA.read_text()
        freeze()
        self.assertEqual(first, FOCAL_QUARTET_SHA.read_text())


def run_tests() -> bool:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage0Tests)
    with tempfile.TemporaryDirectory(prefix="fire_ants_stage0_tests_"):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return result.wasSuccessful()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run embedded Stage-0 validation tests after freezing outputs")
    args = parser.parse_args(argv)

    try:
        freeze()
    except Stage0Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.run_tests:
        return 0 if run_tests() else 1
    print(f"Wrote {REGION_INVENTORY.relative_to(REPO_ROOT)}")
    print(f"Wrote {FOCAL_QUARTET_SHA.relative_to(REPO_ROOT)}")
    print(f"Wrote {STAGE0_REPORT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
