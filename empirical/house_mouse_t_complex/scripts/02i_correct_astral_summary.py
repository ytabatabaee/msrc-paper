#!/usr/bin/env python3
"""Regenerate the compact Stage-2 ASTRAL summary from split-specific annotations."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    Q_SPECIES,
    Q_T_ALT,
    RESULTS,
    annotation_for_split,
    assert_q_normalized,
    canonical_split,
    read_tsv,
    write_tsv,
)


OUT = RESULTS / "stage2_astral_summary.tsv"
FOCAL = {
    "T0_STANDARD": (RESULTS / "stage1_aster" / "T0_STANDARD_subspecies.nwk", Q_SPECIES),
    "T1_ALL_WINDOWS": (RESULTS / "stage1_aster" / "T1_ALL_WINDOWS_subspecies.nwk", Q_T_ALT),
}
FOCAL_SPLITS = {
    Q_SPECIES: canonical_split([
        ["Mus_musculus_musculus", "Mus_musculus_castaneus"],
        ["Mus_musculus_domesticus", "Mus_spretus"],
    ]),
    Q_T_ALT: canonical_split([
        ["Mus_musculus_domesticus", "Mus_musculus_musculus"],
        ["Mus_musculus_castaneus", "Mus_spretus"],
    ]),
}


def run() -> None:
    rows = []
    for treatment, (path, topology) in FOCAL.items():
        annotation = annotation_for_split(path.read_text().strip(), FOCAL_SPLITS[topology])
        if annotation is None:
            raise RuntimeError(f"No focal annotation found for {treatment}")
        assert_q_normalized(annotation)
        rows.append({
            "treatment": treatment,
            "n_gene_trees": 4046,
            "topology": topology,
            "rf_distance_vs_P0": "NA",
            "changed_splits_vs_P0": "NA",
            "CU_length": annotation["CULength"],
            "localPP": annotation["localPP"],
            "q1": annotation["q1"],
            "q2": annotation["q2"],
            "q3": annotation["q3"],
            "status": "focal_split_annotation",
        })
    balanced = {row["treatment"]: row for row in read_tsv(RESULTS / "stage2_balanced_sampling_summary.tsv")}
    for treatment in ["B0_STANDARD_MATCHED", "B1_T_MATCHED", "B2_MIXED_BALANCED"]:
        row = balanced[treatment]
        rows.append({
            "treatment": treatment,
            "n_gene_trees": 4046,
            "topology": "see_stage2_balanced_sampling_summary.tsv",
            "rf_distance_vs_P0": "NA",
            "changed_splits_vs_P0": "NA",
            "CU_length": row["mean_CU_length"],
            "localPP": row["mean_localPP"],
            "q1": "NA",
            "q2": "NA",
            "q3": "NA",
            "status": "ASTRAL_replicate_summary_no_single_focal_split",
        })
    population = read_tsv(RESULTS / "stage2_population_astral.tsv")
    for row in population:
        rows.append({
            "treatment": row["treatment"],
            "n_gene_trees": row["n_gene_trees"],
            "topology": row["treatment"],
            "rf_distance_vs_P0": row["rf_distance_vs_P0"],
            "changed_splits_vs_P0": row["changed_splits_vs_P0"],
            "CU_length": "NA",
            "localPP": "NA",
            "q1": "NA",
            "q2": "NA",
            "q3": "NA",
            "status": "population_split_table_required",
        })
    write_tsv(OUT, rows, ["treatment", "n_gene_trees", "topology", "rf_distance_vs_P0", "changed_splits_vs_P0", "CU_length", "localPP", "q1", "q2", "q3", "status"])
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    run()
