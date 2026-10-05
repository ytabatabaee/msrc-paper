#!/usr/bin/env python3
"""Stage 2B deterministic balanced-sampling ASTRAL controls."""

from __future__ import annotations

import argparse
import json
import sys
import unittest
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    Q_OTHER,
    Q_SPECIES,
    Q_T_ALT,
    RESULTS,
    STAGE2_DATA,
    deterministic_balanced_samples,
    find_astral,
    load_stage1,
    mapping_by_tip,
    prune_newick,
    rel,
    run_astral,
    parse_astral_result,
    write_map,
    write_tree_file,
    write_tsv,
)


OUT_DIR = RESULTS / "stage2_balanced_astral"
DATA_DIR = STAGE2_DATA / "balanced_astral"
DETAIL = RESULTS / "stage2_balanced_sampling_astral.tsv"
SUMMARY = RESULTS / "stage2_balanced_sampling_summary.tsv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-replicates", type=int, default=200)
    p.add_argument("--seed", type=int, default=20261005)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--run-tests", action="store_true")
    return p.parse_args()


def groups_for_selected(selected: list[str], by_tip: dict[str, dict[str, str]]) -> dict[str, int]:
    c = Counter(by_tip[t]["subspecies"] for t in selected)
    return dict(c)


def run(args: argparse.Namespace) -> None:
    meta, trees, mapping = load_stage1()
    by_tip = mapping_by_tip(mapping)
    astral = find_astral()
    samples = deterministic_balanced_samples(mapping, args.n_replicates, args.seed)
    rows = []
    for sample in samples:
        rep = int(sample["replicate"])
        treatment = str(sample["treatment"])
        selected = list(sample["selected_tips"])
        rep_dir = DATA_DIR / f"replicate_{rep:04d}"
        tree_file = rep_dir / f"{treatment}.tre"
        map_file = rep_dir / f"{treatment}.map"
        output = OUT_DIR / f"replicate_{rep:04d}_{treatment}.nwk"
        if not tree_file.exists():
            write_tree_file(tree_file, [prune_newick(nw, set(selected)) for nw in trees])
        if not map_file.exists():
            write_map(map_file, mapping, selected, "subspecies")
        log = output.with_suffix(".log")
        if output.exists() and log.exists():
            result = parse_astral_result(output, log) | {"tree_file": rel(tree_file), "map_file": rel(map_file)}
        else:
            result = run_astral(astral, tree_file, map_file, output, threads=args.threads, seed=args.seed + rep)
        topology = str(result.get("topology", "NA"))
        rows.append(
            {
                "replicate": rep,
                "treatment": treatment,
                "selected_tips": ",".join(selected),
                "group_counts": json.dumps(groups_for_selected(selected, by_tip), sort_keys=True),
                "n_gene_trees": result.get("n_gene_trees", len(trees)),
                "inferred_unrooted_split": topology,
                "matches_q_species": topology == Q_SPECIES,
                "matches_q_t_alt": topology == Q_T_ALT,
                "matches_q_other": topology == Q_OTHER,
                "CU_length": result.get("CU_length", "NA"),
                "SU_length": result.get("SU_length", "NA"),
                "localPP": result.get("localPP", "NA"),
                "q1": result.get("q1", "NA"),
                "q2": result.get("q2", "NA"),
                "q3": result.get("q3", "NA"),
                "astral_score": result.get("astral_score", "NA"),
                "tree_file": result.get("tree_file", rel(tree_file)),
                "map_file": result.get("map_file", rel(map_file)),
                "output_tree": result.get("tree_output", rel(output)),
            }
        )
    fields = [
        "replicate",
        "treatment",
        "selected_tips",
        "group_counts",
        "n_gene_trees",
        "inferred_unrooted_split",
        "matches_q_species",
        "matches_q_t_alt",
        "matches_q_other",
        "CU_length",
        "SU_length",
        "localPP",
        "q1",
        "q2",
        "q3",
        "astral_score",
        "tree_file",
        "map_file",
        "output_tree",
    ]
    rows.sort(key=lambda r: (int(r["replicate"]), str(r["treatment"])))
    write_tsv(DETAIL, rows, fields)
    by_treatment: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_treatment[str(row["treatment"])].append(row)
    summary = []
    for treatment, subset in sorted(by_treatment.items()):
        n = len(subset)
        counts = Counter(str(r["inferred_unrooted_split"]) for r in subset)
        summary.append(
            {
                "treatment": treatment,
                "n_replicates": n,
                "fraction_Q_SPECIES": counts[Q_SPECIES] / n,
                "fraction_Q_T_ALT": counts[Q_T_ALT] / n,
                "fraction_Q_OTHER": counts[Q_OTHER] / n,
                "other_or_unresolved": (n - counts[Q_SPECIES] - counts[Q_T_ALT] - counts[Q_OTHER]) / n,
                "mean_localPP": mean_float([r["localPP"] for r in subset]),
                "mean_CU_length": mean_float([r["CU_length"] for r in subset]),
            }
        )
    write_tsv(
        SUMMARY,
        summary,
        [
            "treatment",
            "n_replicates",
            "fraction_Q_SPECIES",
            "fraction_Q_T_ALT",
            "fraction_Q_OTHER",
            "other_or_unresolved",
            "mean_localPP",
            "mean_CU_length",
        ],
    )
    print(f"Wrote {DETAIL}")
    print(f"Wrote {SUMMARY}")


def mean_float(values: list[object]) -> float:
    xs = [float(v) for v in values if str(v) != "NA"]
    return sum(xs) / len(xs) if xs else float("nan")


class BalancedTests(unittest.TestCase):
    def test_determinism_and_counts(self) -> None:
        mapping = [
            {"tree_tip": f"{sp}_{status}_{i}", "subspecies": sp, "t_status": t_status}
            for sp in ["Mus musculus castaneus", "Mus musculus domesticus", "Mus musculus musculus"]
            for status, t_status in [("std", "standard_noncarrier"), ("t", "pseudo-t_haplotype")]
            for i in range(5)
        ]
        mapping.extend({"tree_tip": f"s_{i}", "subspecies": "Mus spretus", "t_status": "outgroup_not_t_haplotype"} for i in range(8))
        a = deterministic_balanced_samples(mapping, 3, 17)
        b = deterministic_balanced_samples(mapping, 3, 17)
        self.assertEqual(a, b)
        by_tip = mapping_by_tip(mapping)
        for row in a:
            counts = groups_for_selected(row["selected_tips"], by_tip)
            self.assertEqual(counts["Mus spretus"], 3)
            self.assertEqual(counts["Mus musculus castaneus"], 6 if row["treatment"] == "B2_MIXED_BALANCED" else 3)


def main() -> int:
    args = parse_args()
    if args.run_tests:
        return 0 if unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BalancedTests)).wasSuccessful() else 1
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
