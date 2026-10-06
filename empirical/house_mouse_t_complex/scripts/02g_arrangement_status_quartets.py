#!/usr/bin/env python3
"""Decompose four-group quartet support by standard/pseudo-t arrangement."""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    Q_OTHER,
    Q_SPECIES,
    Q_T_ALT,
    Q_UNRESOLVED,
    RESULTS,
    load_stage1,
    read_tsv,
    score_to_row,
    topological_distance_matrix,
    write_tsv,
)


PATTERNS = ["SSS", "SST", "STS", "STT", "TSS", "TST", "TTS", "TTT"]
GROUPS = [
    "Mus musculus domesticus",
    "Mus musculus musculus",
    "Mus musculus castaneus",
]
OUT = RESULTS / "stage2_arrangement_pattern_quartets.tsv"
SUMMARY = RESULTS / "stage2_arrangement_pattern_summary.tsv"
SPATIAL = RESULTS / "stage2_arrangement_pattern_spatial_summary.tsv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-tests", action="store_true")
    return p.parse_args()


def count_pattern(newick: str, pattern: str, mapping: list[dict[str, str]]) -> dict[str, int]:
    selected: list[str] = []
    pattern_groups: dict[str, list[str]] = {}
    for group, state in zip(GROUPS, pattern):
        tips = sorted(
            row["tree_tip"]
            for row in mapping
            if row["subspecies"] == group
            and ((state == "T" and row["t_status"] == "pseudo-t_haplotype") or
                 (state == "S" and row["t_status"] == "standard_noncarrier"))
        )
        pattern_groups[group] = tips
        selected.extend(tips)
    outgroup = sorted(row["tree_tip"] for row in mapping if row["subspecies"] == "Mus spretus")
    selected.extend(outgroup)
    if not all(pattern_groups.values()) or not outgroup:
        return {Q_SPECIES: 0, Q_T_ALT: 0, Q_OTHER: 0, Q_UNRESOLVED: 0}
    tips = sorted(set(selected))
    index = {tip: i for i, tip in enumerate(tips)}
    dist = topological_distance_matrix(newick, tips)
    d = [index[x] for x in pattern_groups[GROUPS[0]]]
    m = [index[x] for x in pattern_groups[GROUPS[1]]]
    c = [index[x] for x in pattern_groups[GROUPS[2]]]
    s = [index[x] for x in outgroup]
    import numpy as np
    dd, mm, cc, ss = np.meshgrid(d, m, c, s, indexing="ij")
    species = dist[mm, cc] + dist[dd, ss]
    t_alt = dist[dd, mm] + dist[cc, ss]
    other = dist[dd, cc] + dist[mm, ss]
    stacked = np.stack([species, t_alt, other], axis=0)
    minimum = stacked.min(axis=0)
    ties = (stacked == minimum).sum(axis=0)
    return {
        Q_SPECIES: int(((species == minimum) & (ties == 1)).sum()),
        Q_T_ALT: int(((t_alt == minimum) & (ties == 1)).sum()),
        Q_OTHER: int(((other == minimum) & (ties == 1)).sum()),
        Q_UNRESOLVED: int((ties != 1).sum()),
    }


def pattern_pool(pattern: str) -> str:
    if pattern in {"SSS", "TTT"}:
        return "all_same_state"
    if pattern.count("T") == 1:
        return "one_T"
    if pattern.count("T") == 2:
        return "two_T"
    return pattern


def run() -> None:
    meta, trees, mapping = load_stage1()
    rows: list[dict[str, object]] = []
    pooled: dict[str, Counter[str]] = defaultdict(Counter)
    windows: dict[str, set[str]] = defaultdict(set)
    spatial: list[dict[str, object]] = []
    scan = {row["locus_id"]: row for row in read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv") if row["treatment"] == "ALL_TIPS"}
    for record, tree in zip(meta, trees):
        locus = record["locus_id"]
        spatial_row: dict[str, object] = {key: record[key] for key in ("locus_id", "start_bp", "end_bp", "midpoint_bp")}
        for pattern in PATTERNS:
            counts = count_pattern(tree, pattern, mapping)
            scored = score_to_row(counts)
            row = {
                "locus_id": locus,
                "start_bp": record["start_bp"],
                "end_bp": record["end_bp"],
                "midpoint_bp": record["midpoint_bp"],
                "status_pattern": pattern,
                "n_quartets": sum(counts.values()),
                "q_species": scored["q_species"],
                "q_t_alt": scored["q_t_alt"],
                "q_other": scored["q_other"],
                "q_unresolved": scored["q_unresolved"],
                "delta_species_alt": scored["delta_species_alt"],
            }
            rows.append(row)
            pool = pattern_pool(pattern)
            for key, value in counts.items():
                pooled[pattern][key] += value
                pooled[pool][key] += value
            if counts[Q_SPECIES] + counts[Q_T_ALT] + counts[Q_OTHER] + counts[Q_UNRESOLVED] > 0:
                windows[pattern].add(locus)
                windows[pool].add(locus)
            if pattern in {"SSS", "TTT"}:
                spatial_row[f"{pattern}_delta"] = scored["delta_species_alt"]
        all_tip = scan[locus]
        spatial_row["one_T_delta"] = "NA"
        spatial_row["two_T_delta"] = "NA"
        spatial_row["ALL_TIPS_delta"] = all_tip["delta_species_alt"]
        pattern_scores = {r["status_pattern"]: r for r in rows[-8:]}
        for pool in ("one_T", "two_T"):
            subset = [pattern_scores[p] for p in PATTERNS if pattern_pool(p) == pool]
            totals = sum(float(x["n_quartets"]) for x in subset)
            if totals:
                spatial_row[f"{pool}_delta"] = sum(float(x["n_quartets"]) * float(x["delta_species_alt"]) for x in subset) / totals
        spatial.append(spatial_row)
    fields = ["locus_id", "start_bp", "end_bp", "midpoint_bp", "status_pattern", "n_quartets", "q_species", "q_t_alt", "q_other", "q_unresolved", "delta_species_alt"]
    write_tsv(OUT, rows, fields)
    summary_rows = []
    for pattern in PATTERNS + ["all_same_state", "one_T", "two_T"]:
        counts = pooled[pattern]
        scored = score_to_row(counts)
        summary_rows.append({
            "status_pattern": pattern,
            "n_induced_quartets": sum(counts.values()),
            "n_windows_with_pattern": len(windows[pattern]),
            "q_species": scored["q_species"],
            "q_t_alt": scored["q_t_alt"],
            "q_other": scored["q_other"],
            "q_unresolved": scored["q_unresolved"],
            "delta_species_alt": scored["delta_species_alt"],
        })
    write_tsv(SUMMARY, summary_rows, ["status_pattern", "n_induced_quartets", "n_windows_with_pattern", "q_species", "q_t_alt", "q_other", "q_unresolved", "delta_species_alt"])
    write_tsv(SPATIAL, spatial, ["locus_id", "start_bp", "end_bp", "midpoint_bp", "SSS_delta", "TTT_delta", "one_T_delta", "two_T_delta", "ALL_TIPS_delta"])
    print(f"Wrote {OUT}\nWrote {SUMMARY}\nWrote {SPATIAL}")


def main() -> int:
    args = parse_args()
    if args.run_tests:
        for pattern in PATTERNS:
            assert len(pattern) == 3 and set(pattern) <= {"S", "T"}
        return 0
    run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
