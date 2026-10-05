#!/usr/bin/env python3
"""Stage 2C seven-population ASTRAL sensitivity."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    STAGE2_DATA,
    find_astral,
    load_stage1,
    prune_newick,
    rf_distance,
    run_astral,
    tree_bipartitions,
    write_map,
    write_tree_file,
    write_tsv,
)


OUT = RESULTS / "stage2_population_astral.tsv"
OUT_DIR = RESULTS / "stage2_population_astral"
DATA_DIR = STAGE2_DATA / "population_astral"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--threads", type=int, default=2)
    return p.parse_args()


def run(args: argparse.Namespace) -> None:
    meta, trees, mapping = load_stage1()
    astral = find_astral()
    standard = {r["tree_tip"] for r in mapping if r["t_status"] in {"standard_noncarrier", "outgroup_not_t_haplotype"}}
    all_tips = {r["tree_tip"] for r in mapping}
    treatments = {
        "P0_STANDARD": (trees, standard),
        "P1_ALL": (trees, all_tips),
        "P2_THIN_50kb": ([trees[i] for i in thin_indices(meta, 50_000)], all_tips),
        "P3_THIN_100kb": ([trees[i] for i in thin_indices(meta, 100_000)], all_tips),
    }
    rows = []
    newicks = {}
    for name, (newick_list, tips) in treatments.items():
        tree_file = DATA_DIR / f"{name}.tre"
        map_file = DATA_DIR / f"{name}.map"
        if not tree_file.exists():
            write_tree_file(tree_file, [prune_newick(nw, tips) if tips != all_tips else nw for nw in newick_list])
        write_map(map_file, mapping, tips, "population")
        result = run_astral(astral, tree_file, map_file, OUT_DIR / f"{name}.nwk", threads=args.threads)
        newicks[name] = str(result["newick"])
        rows.append(
            {
                "treatment": name,
                "n_gene_trees": result.get("n_gene_trees", len(newick_list)),
                "n_populations": result.get("n_species", "NA"),
                "newick": result.get("newick", ""),
                "CU_length": result.get("CU_length", "NA"),
                "localPP": result.get("localPP", "NA"),
                "q1": result.get("q1", "NA"),
                "q2": result.get("q2", "NA"),
                "q3": result.get("q3", "NA"),
            }
        )
    p0 = newicks["P0_STANDARD"]
    p1 = newicks["P1_ALL"]
    changed = sorted(tree_bipartitions(p0) ^ tree_bipartitions(p1))
    for row in rows:
        row["rf_distance_vs_P0"] = 0 if row["treatment"] == "P0_STANDARD" else rf_distance(p0, newicks[row["treatment"]])
        row["changed_splits_vs_P0"] = "NA" if row["treatment"] == "P0_STANDARD" else len(tree_bipartitions(p0) ^ tree_bipartitions(newicks[row["treatment"]]))
    write_tsv(
        OUT,
        rows,
        ["treatment", "n_gene_trees", "n_populations", "rf_distance_vs_P0", "changed_splits_vs_P0", "CU_length", "localPP", "q1", "q2", "q3", "newick"],
    )
    print(f"Wrote {OUT}")


def thin_indices(meta: list[dict[str, str]], spacing: int) -> list[int]:
    bins = {}
    for i, row in enumerate(meta):
        start = int(row["start_bp"])
        b = (start - 5_000_000) // spacing
        if b not in bins or start < bins[b][0]:
            bins[b] = (start, i)
    return [i for _start, i in sorted(bins.values())]


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))

