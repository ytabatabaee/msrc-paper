#!/usr/bin/env python3
"""Compare frozen de novo structural breakpoints with published outlier bounds."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
MANIFEST = RESULTS_DIR / "denovo_breakpoint_manifest.json"
CONSENSUS = RESULTS_DIR / "denovo_consensus_structural_breakpoints.tsv"
PUBLISHED = REPO_ROOT / "rearrangements" / "outlier-regions-by-maf2synteny.tsv"
OUT = RESULTS_DIR / "denovo_vs_published_boundaries.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def numeric_values(row: dict[str, str]) -> list[int]:
    values = []
    for key, value in row.items():
        if not isinstance(key, str):
            continue
        if not value:
            continue
        low = key.lower()
        if any(token in low for token in ["start", "end", "bound", "break", "coord", "position", "pos"]):
            for hit in re.findall(r"\d+", value.replace(",", "")):
                val = int(hit)
                if 1 <= val <= 100_000_000:
                    values.append(val)
    return values


def read_published_boundaries(path: Path) -> list[dict[str, object]]:
    out = []
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        for idx, fields in enumerate(reader, 1):
            cleaned = [f.strip().strip('"') for f in fields]
            bounds = []
            # The published table has an unlabeled leading row index. The final
            # triplet is the GalGal6 axis: chr4, Gal_start, Gal_end.
            if len(cleaned) >= 3 and cleaned[-3] == "chr4":
                bounds = [int(cleaned[-2]), int(cleaned[-1])]
            else:
                row = dict(zip(header or [], cleaned))
                vals = sorted(set(numeric_values(row)))
                if len(vals) >= 2:
                    bounds = [min(vals), max(vals)]
                elif len(vals) == 1:
                    bounds = [vals[0]]
            if not bounds:
                continue
            label = cleaned[1] if len(cleaned) > 1 else f"published_row_{idx}"
            for side, pos in zip(["left", "right"], sorted(bounds)):
                out.append({"published_boundary_id": f"PUB_{idx:03d}_{side}", "published_region": label, "published_position": pos})
    dedup = {}
    for row in out:
        dedup[(row["published_region"], row["published_position"])] = row
    return sorted(dedup.values(), key=lambda r: int(r["published_position"]))


def published_boundaries(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out = []
    for idx, row in enumerate(rows, 1):
        vals = sorted(set(numeric_values(row)))
        if len(vals) >= 2:
            bounds = [min(vals), max(vals)]
        elif len(vals) == 1:
            bounds = [vals[0]]
        else:
            continue
        label = row.get("region", row.get("outlier_region", row.get("id", f"published_row_{idx}")))
        for side, pos in zip(["left", "right"], bounds):
            out.append({"published_boundary_id": f"PUB_{idx:03d}_{side}", "published_region": label, "published_position": pos})
    dedup = {}
    for row in out:
        dedup[(row["published_region"], row["published_position"])] = row
    return sorted(dedup.values(), key=lambda r: int(r["published_position"]))


def compare() -> list[dict[str, object]]:
    if not MANIFEST.exists():
        raise SystemExit(f"Refusing comparison before freeze manifest exists: {MANIFEST}")
    consensus = read_tsv(CONSENSUS)
    published = read_published_boundaries(PUBLISHED)
    rows = []
    for bp in consensus:
        pos = int(float(bp["reference_position"]))
        nearest = min(published, key=lambda p: abs(pos - int(p["published_position"]))) if published else None
        inside = [
            p["published_region"]
            for p in published
            if abs(pos - int(p["published_position"])) == 0
        ]
        rows.append(
            {
                "consensus_breakpoint_id": bp["consensus_breakpoint_id"],
                "denovo_reference_position": pos,
                "confidence_class": bp["confidence_class"],
                "n_species_support": bp["n_species_support"],
                "nearest_published_boundary_id": nearest["published_boundary_id"] if nearest else "",
                "nearest_published_region": nearest["published_region"] if nearest else "",
                "nearest_published_position": nearest["published_position"] if nearest else "",
                "distance_bp": abs(pos - int(nearest["published_position"])) if nearest else "",
                "exact_published_boundary_overlap": int(bool(inside)),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    rows = compare()
    write_tsv(
        OUT,
        rows,
        [
            "consensus_breakpoint_id", "denovo_reference_position", "confidence_class",
            "n_species_support", "nearest_published_boundary_id", "nearest_published_region",
            "nearest_published_position", "distance_bp", "exact_published_boundary_overlap",
        ],
    )
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
