#!/usr/bin/env python3
"""Design-agnostic real-data gate for Anopheles 2La Stage 1B."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path


NO_ELIGIBLE_STATUS = "NO_ELIGIBLE_PROSPECTIVE_DESIGNS"
NO_ELIGIBLE_MESSAGE = (
    "Stage 1A completed successfully, but no prospectively eligible "
    "structural-state analysis designs were available."
)

DESIGN_A_TABLES = ("frozen_strict_quartets_stage1a.tsv",)
DESIGN_B_TABLES = (
    "frozen_design_b_contrasts_stage1a.tsv",
    "frozen_arrangement_replacement_contrasts_stage1a.tsv",
)
DESIGN_C_TABLES = (
    "frozen_design_c_geography_controls_stage1a.tsv",
    "frozen_geography_population_controls_stage1a.tsv",
)


@dataclass(frozen=True)
class Stage1BGateResult:
    status: str
    design_A_enabled: bool
    design_B_enabled: bool
    design_C_enabled: bool
    design_A_count: int
    design_B_count: int
    design_C_count: int
    message: str

    @property
    def enabled_designs(self) -> set[str]:
        enabled = set()
        if self.design_A_enabled:
            enabled.add("A")
        if self.design_B_enabled:
            enabled.add("B")
        if self.design_C_enabled:
            enabled.add("C")
        return enabled

    def as_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "design_A_enabled": self.design_A_enabled,
            "design_B_enabled": self.design_B_enabled,
            "design_C_enabled": self.design_C_enabled,
            "design_A_count": self.design_A_count,
            "design_B_count": self.design_B_count,
            "design_C_count": self.design_C_count,
            "message": self.message,
        }


class Stage1BGateError(RuntimeError):
    pass


class NoEligibleProspectiveDesigns(RuntimeError):
    def __init__(self, result: Stage1BGateResult):
        self.result = result
        super().__init__(f"{result.status}: {result.message}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _nonempty_rows(path: Path) -> int:
    return len(_read_rows(path))


def _count_design_rows(path: Path, token: str, fallback_all_rows: bool = False) -> int:
    rows = _read_rows(path)
    if not rows:
        return 0
    if "design_class" not in rows[0]:
        return len(rows) if fallback_all_rows else 0
    return sum(1 for row in rows if token in row.get("design_class", ""))


def _validate_marker(processed_root: Path) -> None:
    marker = processed_root / "stage1a_freeze_complete.json"
    if not marker.exists():
        raise Stage1BGateError(f"real Stage 1B is blocked until Stage 1A gate exists: {marker}")
    try:
        payload = json.loads(marker.read_text())
    except json.JSONDecodeError as exc:
        raise Stage1BGateError(f"Stage 1A completion marker is not valid JSON: {marker}") from exc
    if not isinstance(payload, dict):
        raise Stage1BGateError(f"Stage 1A completion marker must be a JSON object: {marker}")
    blocked_statuses = {"AUTHENTICATION_REQUIRED", "API_ERROR", "MISSING_SAMPLE_SET", "FAILED"}
    status = str(payload.get("status", payload.get("retrieval_status", ""))).upper()
    if status in blocked_statuses:
        raise Stage1BGateError(f"Stage 1A completion marker is not a successful scientific completion: {status}")


def _validate_checksum(table: Path) -> None:
    checksum_file = table.with_suffix(".sha256")
    if not table.exists() or not checksum_file.exists():
        return
    expected = checksum_file.read_text().split()[0]
    import hashlib

    h = hashlib.sha256()
    with table.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    observed = h.hexdigest()
    if observed != expected:
        raise Stage1BGateError(f"checksum validation failed for {table}: expected {expected}, observed {observed}")


def _validate_provenance(data_root: Path) -> None:
    provenance = data_root / "metadata" / "stage1a_malariagen_api_provenance.json"
    if not provenance.exists():
        return
    try:
        payload = json.loads(provenance.read_text())
    except json.JSONDecodeError as exc:
        raise Stage1BGateError(f"Stage 1A provenance is not valid JSON: {provenance}") from exc
    if not isinstance(payload, dict):
        raise Stage1BGateError(f"Stage 1A provenance must be a JSON object: {provenance}")


def evaluate_stage1b_gate(data_root: Path) -> Stage1BGateResult:
    processed = data_root / "processed"
    _validate_marker(processed)
    _validate_provenance(data_root)

    combined = processed / "frozen_strict_quartets_stage1a.tsv"
    for table in [combined, *(processed / name for name in DESIGN_B_TABLES), *(processed / name for name in DESIGN_C_TABLES)]:
        _validate_checksum(table)

    design_a_count = _count_design_rows(combined, "Design_A", fallback_all_rows=True)
    design_b_count = _count_design_rows(combined, "Design_B") + sum(_nonempty_rows(processed / name) for name in DESIGN_B_TABLES)
    design_c_count = _count_design_rows(combined, "Design_C") + sum(_nonempty_rows(processed / name) for name in DESIGN_C_TABLES)

    result = Stage1BGateResult(
        status="READY" if any([design_a_count, design_b_count, design_c_count]) else NO_ELIGIBLE_STATUS,
        design_A_enabled=design_a_count > 0,
        design_B_enabled=design_b_count > 0,
        design_C_enabled=design_c_count > 0,
        design_A_count=design_a_count,
        design_B_count=design_b_count,
        design_C_count=design_c_count,
        message=(
            "Stage 1B real-data gate passed with at least one prospectively frozen design."
            if any([design_a_count, design_b_count, design_c_count])
            else NO_ELIGIBLE_MESSAGE
        ),
    )
    return result


def require_stage1b_designs(data_root: Path) -> Stage1BGateResult:
    result = evaluate_stage1b_gate(data_root)
    if result.status == NO_ELIGIBLE_STATUS:
        raise NoEligibleProspectiveDesigns(result)
    return result
