#!/usr/bin/env python3
"""Stage 2A exact fixed-quartet scoring for house-mouse local trees."""

from __future__ import annotations

import argparse
import sys
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    Q_OTHER,
    Q_SPECIES,
    Q_T_ALT,
    RESULTS,
    aggregate_quartet_rows,
    load_stage1,
    quartet_counts_for_tree,
    score_to_row,
    tips_by_species,
    write_tsv,
)


SCAN = RESULTS / "stage2_fixed_quartet_scan.tsv"
SUMMARY = RESULTS / "stage2_fixed_quartet_summary.tsv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-tests", action="store_true")
    return p.parse_args()


def run() -> None:
    meta, trees, mapping = load_stage1()
    treatments = ["STANDARD_ONLY", "T_ONLY", "ALL_TIPS"]
    rows = []
    for treatment in treatments:
        groups = tips_by_species(mapping, treatment)
        for locus, newick in zip(meta, trees, strict=True):
            counts = quartet_counts_for_tree(newick, groups)
            row = {
                "locus_id": locus["locus_id"],
                "start_bp": locus["start_bp"],
                "end_bp": locus["end_bp"],
                "midpoint_bp": locus["midpoint_bp"],
                "treatment": treatment,
            }
            row.update(score_to_row(counts))
            rows.append(row)
    fields = [
        "locus_id",
        "start_bp",
        "end_bp",
        "midpoint_bp",
        "treatment",
        "n_resolved_quartets",
        "n_unresolved_quartets",
        "q_species",
        "q_t_alt",
        "q_other",
        "q_unresolved",
        "delta_species_alt",
    ]
    write_tsv(SCAN, rows, fields)
    write_tsv(
        SUMMARY,
        aggregate_quartet_rows([{k: str(v) for k, v in row.items()} for row in rows]),
        [
            "treatment",
            "n_loci",
            "mean_q_species",
            "mean_q_t_alt",
            "mean_q_other",
            "mean_q_unresolved",
            "mean_delta_species_alt",
        ],
    )
    for row in rows:
        total = sum(float(row[k]) for k in ["q_species", "q_t_alt", "q_other", "q_unresolved"])
        if abs(total - 1.0) > 1e-9:
            raise RuntimeError(f"Quartet proportions do not sum to 1: {row['locus_id']} {row['treatment']} {total}")
    print(f"Wrote {SCAN}")
    print(f"Wrote {SUMMARY}")


class QuartetTests(unittest.TestCase):
    groups = {
        "Mus musculus domesticus": ["d"],
        "Mus musculus musculus": ["m"],
        "Mus musculus castaneus": ["c"],
        "Mus spretus": ["s"],
    }

    def assert_topology(self, newick: str, expected: str) -> None:
        counts = quartet_counts_for_tree(newick, self.groups)
        for key in [Q_SPECIES, Q_T_ALT, Q_OTHER]:
            self.assertEqual(counts[key], 1 if key == expected else 0)

    def test_species(self) -> None:
        self.assert_topology("((m,c),(d,s));", Q_SPECIES)

    def test_t_alt(self) -> None:
        self.assert_topology("((d,m),(c,s));", Q_T_ALT)

    def test_other(self) -> None:
        self.assert_topology("((d,c),(m,s));", Q_OTHER)

    def test_polytomy(self) -> None:
        counts = quartet_counts_for_tree("(d,m,c,s);", self.groups)
        self.assertEqual(counts["Q_UNRESOLVED"], 1)


def main() -> int:
    args = parse_args()
    if args.run_tests:
        return 0 if unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QuartetTests)).wasSuccessful() else 1
    run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

