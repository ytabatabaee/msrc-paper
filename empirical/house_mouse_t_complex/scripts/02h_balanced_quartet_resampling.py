#!/usr/bin/env python3
"""Fast exact balanced quartet resampling across the frozen local trees."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    Q_OTHER,
    Q_SPECIES,
    Q_T_ALT,
    RESULTS,
    deterministic_balanced_samples,
    load_stage1,
    mapping_by_tip,
    topological_distance_matrix,
    read_tsv,
    write_tsv,
)


OUT = RESULTS / "stage2_balanced_quartet_resampling.tsv"
SUMMARY = RESULTS / "stage2_balanced_quartet_resampling_summary.tsv"
VALIDATION = RESULTS / "stage2_balanced_quartet_astral_validation.tsv"
SPECIES = ["Mus musculus domesticus", "Mus musculus musculus", "Mus musculus castaneus"]
ALL_GROUPS = SPECIES + ["Mus spretus"]
Q_NAMES = [Q_SPECIES, Q_T_ALT, Q_OTHER]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-replicates", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261005)
    p.add_argument("--run-tests", action="store_true")
    p.add_argument("--validate-astral", action="store_true")
    return p.parse_args()


def group_tips(mapping: list[dict[str, str]]) -> dict[str, list[str]]:
    out = {group: [] for group in ALL_GROUPS}
    for row in mapping:
        if row["subspecies"] not in out:
            continue
        if row["subspecies"] == "Mus spretus" or row["mapping_confidence"] in {"strong", "exact"}:
            out[row["subspecies"]].append(row["tree_tip"])
    for group in out:
        out[group].sort()
    return out


def build_total_counts(trees: list[str], mapping: list[dict[str, str]]) -> tuple[np.ndarray, dict[str, list[str]]]:
    """Count each possible individual quartet once across all windows.

    The resulting tensor has axes domesticus, musculus, castaneus, spretus
    followed by Q_SPECIES/Q_T_ALT/Q_OTHER. It is small (55-tip universe) and
    makes 1,000 resamples a table lookup instead of 3,000 ASTRAL runs.
    """
    groups = group_tips(mapping)
    shape = tuple(len(groups[group]) for group in ALL_GROUPS) + (3,)
    total = np.zeros(shape, dtype=np.int64)
    all_tips = [tip for group in ALL_GROUPS for tip in groups[group]]
    index = {tip: i for i, tip in enumerate(all_tips)}
    offsets = []
    for group in ALL_GROUPS:
        offsets.append(np.array([index[tip] for tip in groups[group]], dtype=np.int32))
    for n, tree in enumerate(trees, start=1):
        dist = topological_distance_matrix(tree, all_tips)
        dd, mm, cc, ss = np.meshgrid(*offsets, indexing="ij")
        species = dist[mm, cc] + dist[dd, ss]
        t_alt = dist[dd, mm] + dist[cc, ss]
        other = dist[dd, cc] + dist[mm, ss]
        stacked = np.stack([species, t_alt, other], axis=0)
        minimum = stacked.min(axis=0)
        ties = (stacked == minimum).sum(axis=0)
        for q, values in enumerate((species, t_alt, other)):
            total[..., q] += ((values == minimum) & (ties == 1)).astype(np.int64)
        if n % 500 == 0:
            print(f"precomputed {n}/{len(trees)} windows", flush=True)
    return total, groups


def counts_for_selection(total: np.ndarray, groups: dict[str, list[str]], selected: list[str], mapping_by_tip_: dict[str, dict[str, str]]) -> np.ndarray:
    indices = []
    for group in ALL_GROUPS:
        chosen = [tip for tip in selected if mapping_by_tip_[tip]["subspecies"] == group]
        indices.append([groups[group].index(tip) for tip in chosen])
    return total[np.ix_(*indices)].sum(axis=tuple(range(4)))


def winner(counts: np.ndarray) -> str:
    maximum = int(counts.max())
    winners = [name for name, value in zip(Q_NAMES, counts) if int(value) == maximum]
    return winners[0] if len(winners) == 1 else "Q_UNRESOLVED"


def run(args: argparse.Namespace) -> None:
    _meta, trees, mapping = load_stage1()
    total, groups = build_total_counts(trees, mapping)
    by_tip = mapping_by_tip(mapping)
    samples = deterministic_balanced_samples(mapping, args.n_replicates, args.seed)
    rows = []
    for sample in samples:
        counts = counts_for_selection(total, groups, list(sample["selected_tips"]), by_tip)
        n = int(counts.sum())
        rows.append({
            "replicate": sample["replicate"],
            "treatment": sample["treatment"],
            "selected_tips": ",".join(sample["selected_tips"]),
            "n_gene_trees": len(trees),
            "n_resolved_quartets": n,
            "q_species": int(counts[0]) / n,
            "q_t_alt": int(counts[1]) / n,
            "q_other": int(counts[2]) / n,
            "winner": winner(counts),
        })
    fields = ["replicate", "treatment", "selected_tips", "n_gene_trees", "n_resolved_quartets", "q_species", "q_t_alt", "q_other", "winner"]
    write_tsv(OUT, rows, fields)
    summaries = []
    for treatment in sorted({str(row["treatment"]) for row in rows}):
        subset = [row for row in rows if row["treatment"] == treatment]
        counts = Counter(str(row["winner"]) for row in subset)
        summaries.append({
            "treatment": treatment,
            "n_replicates": len(subset),
            "fraction_Q_SPECIES": counts[Q_SPECIES] / len(subset),
            "fraction_Q_T_ALT": counts[Q_T_ALT] / len(subset),
            "fraction_Q_OTHER": counts[Q_OTHER] / len(subset),
            "fraction_unresolved": counts["Q_UNRESOLVED"] / len(subset),
            "mean_q_species": np.mean([float(x["q_species"]) for x in subset]),
            "mean_q_t_alt": np.mean([float(x["q_t_alt"]) for x in subset]),
            "mean_q_other": np.mean([float(x["q_other"]) for x in subset]),
        })
    write_tsv(SUMMARY, summaries, ["treatment", "n_replicates", "fraction_Q_SPECIES", "fraction_Q_T_ALT", "fraction_Q_OTHER", "fraction_unresolved", "mean_q_species", "mean_q_t_alt", "mean_q_other"])
    print(f"Wrote {OUT}\nWrote {SUMMARY}")


def validate_astral() -> None:
    exact = {(int(row["replicate"]), row["treatment"]): row for row in read_tsv(OUT)}
    astral_path = RESULTS / "stage2_balanced_sampling_astral.tsv"
    rows = []
    for row in read_tsv(astral_path):
        key = (int(row["replicate"]), row["treatment"])
        if key not in exact:
            continue
        exact_winner = exact[key]["winner"]
        astral_winner = row["inferred_unrooted_split"]
        rows.append({
            "replicate": row["replicate"],
            "treatment": row["treatment"],
            "exact_winner": exact_winner,
            "astral_topology": astral_winner,
            "concordant": exact_winner == astral_winner,
            "exact_q_species": exact[key]["q_species"],
            "exact_q_t_alt": exact[key]["q_t_alt"],
            "exact_q_other": exact[key]["q_other"],
        })
    write_tsv(VALIDATION, rows, ["replicate", "treatment", "exact_winner", "astral_topology", "concordant", "exact_q_species", "exact_q_t_alt", "exact_q_other"])
    print(f"Wrote {VALIDATION} ({len(rows)} matched replicates)")


def main() -> int:
    args = parse_args()
    if args.run_tests:
        assert winner(np.array([4, 2, 1])) == Q_SPECIES
        assert winner(np.array([1, 4, 2])) == Q_T_ALT
        assert winner(np.array([1, 2, 4])) == Q_OTHER
        return 0
    if args.validate_astral:
        validate_astral()
        return 0
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
