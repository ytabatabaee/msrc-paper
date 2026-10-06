#!/usr/bin/env python3
"""Stage 2E spatial thinning and pseudo-replication sensitivity."""

from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    STAGE2_DATA,
    find_astral,
    load_stage1,
    parse_astral_result,
    read_tsv,
    rel,
    run_astral,
    spatial_thin_records,
    write_tsv,
    write_tree_file,
)


OUT = RESULTS / "stage2_linkage_sensitivity.tsv"
BLOCKS = RESULTS / "stage2_source_defined_blocks.tsv"
OUT_DIR = RESULTS / "stage2_linkage_astral"
DATA_DIR = STAGE2_DATA / "linkage_astral"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path("IST-2017-78-v1+1_Data.zip"))
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--max-phases", type=int, default=None, help="Deprecated phase limit; omit for the complete feasible sweep.")
    p.add_argument("--astral-phases", type=int, default=10, help="Run/cache ASTRAL for this many initial phases per spacing; exact quartet winners are computed for every phase.")
    return p.parse_args()


def run(args: argparse.Namespace) -> None:
    meta, trees, _mapping = load_stage1()
    astral = find_astral()
    map_file = PROCESSED / "house_mouse_t_complex_astral_subspecies.map"
    scan = {row["locus_id"]: row for row in read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv") if row["treatment"] == "ALL_TIPS"}
    rows = []
    for spacing in [5_000, 10_000, 20_000, 50_000, 100_000, 250_000, 500_000]:
        phases = list(range(0, spacing, 5_000))
        if args.max_phases is not None and len(phases) > args.max_phases:
            phases = phases[: args.max_phases]
        for phase in phases:
            idx = spatial_thin_records(meta, spacing, phase)
            name = f"spacing_{spacing // 1000}kb_phase_{phase // 1000}kb"
            selected_ids = {meta[i]["locus_id"] for i in idx}
            exact_counts = {"Q_SPECIES": 0.0, "Q_T_ALT": 0.0, "Q_OTHER": 0.0}
            exact_total = 0.0
            for locus_id in selected_ids:
                record = scan[locus_id]
                denominator = float(record["n_resolved_quartets"]) + float(record["n_unresolved_quartets"])
                exact_total += denominator
                for key, column in (("Q_SPECIES", "q_species"), ("Q_T_ALT", "q_t_alt"), ("Q_OTHER", "q_other")):
                    exact_counts[key] += float(record[column]) * denominator
            exact_topology = max(exact_counts, key=exact_counts.get)
            result = {"topology": exact_topology, "CU_length": "NA", "localPP": "NA", "q1": "NA", "q2": "NA", "q3": "NA"}
            tree_file = DATA_DIR / f"{name}.tre"
            if phase // 5_000 < args.astral_phases:
                if not tree_file.exists():
                    write_tree_file(tree_file, [trees[i] for i in idx])
                astral_output = OUT_DIR / f"{name}.nwk"
                astral_log = astral_output.with_suffix(".log")
                if astral_output.exists() and astral_log.exists():
                    result = parse_astral_result(astral_output, astral_log)
                else:
                    result = run_astral(astral, tree_file, map_file, astral_output, threads=args.threads)
            rows.append(
                {
                    "spacing_kb": spacing // 1000,
                    "phase_kb": phase // 1000,
                    "n_trees": len(idx),
                    "topology": result.get("topology", exact_topology),
                    "exact_topology": exact_topology,
                    "matches_q_species": exact_topology == "Q_SPECIES",
                    "matches_q_t_alt": exact_topology == "Q_T_ALT",
                    "matches_q_other": exact_topology == "Q_OTHER",
                    "exact_delta_species_alt": (exact_counts["Q_SPECIES"] - exact_counts["Q_T_ALT"]) / exact_total,
                    "CU": result.get("CU_length", "NA"),
                    "localPP": result.get("localPP", "NA"),
                    "q1": result.get("q1", "NA"),
                    "q2": result.get("q2", "NA"),
                    "q3": result.get("q3", "NA"),
                    "tree_file": rel(tree_file) if phase // 5_000 < args.astral_phases else "NA",
                }
            )
    write_tsv(OUT, rows, ["spacing_kb", "phase_kb", "n_trees", "topology", "exact_topology", "matches_q_species", "matches_q_t_alt", "matches_q_other", "exact_delta_species_alt", "CU", "localPP", "q1", "q2", "q3", "tree_file"])
    summary = []
    for spacing in sorted({int(row["spacing_kb"]) for row in rows}):
        subset = [row for row in rows if int(row["spacing_kb"]) == spacing]
        astral_rows = [row for row in subset if row["localPP"] != "NA"]
        summary.append({
            "spacing_kb": spacing,
            "n_phases": len(subset),
            "mean_n_trees": sum(int(row["n_trees"]) for row in subset) / len(subset),
            "fraction_Q_SPECIES": sum(row["exact_topology"] == "Q_SPECIES" for row in subset) / len(subset),
            "fraction_Q_T_ALT": sum(row["exact_topology"] == "Q_T_ALT" for row in subset) / len(subset),
            "fraction_Q_OTHER": sum(row["exact_topology"] == "Q_OTHER" for row in subset) / len(subset),
            "median_delta_species_alt": float(np.median([float(row["exact_delta_species_alt"]) for row in subset])),
            "median_localPP": float(np.median([float(row["localPP"]) for row in astral_rows])) if astral_rows else "NA",
            "min_localPP": min(float(row["localPP"]) for row in astral_rows) if astral_rows else "NA",
            "max_localPP": max(float(row["localPP"]) for row in astral_rows) if astral_rows else "NA",
            "median_CU": float(np.median([float(row["CU"]) for row in astral_rows if row["CU"] != "NA"])) if astral_rows else "NA",
        })
    write_tsv(RESULTS / "stage2_linkage_summary.tsv", summary, ["spacing_kb", "n_phases", "mean_n_trees", "fraction_Q_SPECIES", "fraction_Q_T_ALT", "fraction_Q_OTHER", "median_delta_species_alt", "median_localPP", "min_localPP", "max_localPP", "median_CU"])
    write_tsv(BLOCKS, inspect_blocks(args.archive), ["source_defined_blocks_available", "source_path", "notes"])
    print(f"Wrote {OUT}")


def inspect_blocks(archive: Path) -> list[dict[str, object]]:
    if not archive.exists():
        return [{"source_defined_blocks_available": False, "source_path": "NA", "notes": "Archive not found for block inspection."}]
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
    candidates = [n for n in names if "Non-Recombined" in n or "NonRecomb" in n]
    # The archive contains concatenated non-recombined-region tree files, but no interval table.
    return [
        {
            "source_defined_blocks_available": False,
            "source_path": ",".join(candidates),
            "notes": "Source archive names concatenated non-recombined-region tree files but no machine-readable genomic interval definitions were recovered; no post-hoc blocks inferred.",
        }
    ]


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
