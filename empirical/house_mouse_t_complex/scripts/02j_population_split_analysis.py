#!/usr/bin/env python3
"""Extract and compare split-specific seven-population ASTRAL annotations."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import RESULTS, annotated_split_records, canonical_split, write_tsv  # noqa: E402


TREE_DIR = RESULTS / "stage2_population_astral"
OUT = RESULTS / "stage2_population_split_table.tsv"
COMPARE = RESULTS / "stage2_population_split_comparison.tsv"
TREATMENTS = ["P0_STANDARD", "P1_ALL", "P2_THIN_50kb", "P3_THIN_100kb"]


def split_id(split: tuple[tuple[str, ...], tuple[str, ...]]) -> str:
    return "|".join(split[0]) + " || " + "|".join(split[1])


def run() -> None:
    rows = []
    by_treatment = {}
    for treatment in TREATMENTS:
        records = annotated_split_records((TREE_DIR / f"{treatment}.nwk").read_text().strip())
        by_treatment[treatment] = {record["split"]: record for record in records}
        for split, record in sorted(by_treatment[treatment].items()):
            rows.append({
                "treatment": treatment,
                "split_id": split_id(split),
                "side_a": "|".join(split[0]),
                "side_b": "|".join(split[1]),
                "CU": record.get("CULength", "NA"),
                "SU": record.get("SULength", "NA"),
                "localPP": record.get("localPP", "NA"),
                "q1": record.get("q1", "NA"),
                "q2": record.get("q2", "NA"),
                "q3": record.get("q3", "NA"),
            })
    write_tsv(OUT, rows, ["treatment", "split_id", "side_a", "side_b", "CU", "SU", "localPP", "q1", "q2", "q3"])
    p0 = by_treatment["P0_STANDARD"]
    p1 = by_treatment["P1_ALL"]
    comparison = []
    for split in sorted(set(p0) | set(p1)):
        in0 = split in p0
        in1 = split in p1
        status = "shared" if in0 and in1 else "lost" if in0 else "gained"
        record0 = p0.get(split, {})
        record1 = p1.get(split, {})
        comparison.append({
            "split_id": split_id(split),
            "status": status,
            "CU_P0": record0.get("CULength", "NA"),
            "CU_P1": record1.get("CULength", "NA"),
            "delta_CU": float(record1["CULength"]) - float(record0["CULength"]) if in0 and in1 and "CULength" in record0 and "CULength" in record1 else "NA",
            "localPP_P0": record0.get("localPP", "NA"),
            "localPP_P1": record1.get("localPP", "NA"),
        })
    write_tsv(COMPARE, comparison, ["split_id", "status", "CU_P0", "CU_P1", "delta_CU", "localPP_P0", "localPP_P1"])
    print(f"Wrote {OUT}\nWrote {COMPARE}")


if __name__ == "__main__":
    run()
