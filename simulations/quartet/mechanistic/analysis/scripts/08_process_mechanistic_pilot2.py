#!/usr/bin/env python3
"""Process Experiment 2B with history-level uncertainty and bias decomposition."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import yaml


BOOTSTRAPS = 2000
BOOT_SEED = 20261021


def read_rows(path: Path, delimiter=","):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def mean_ci(values):
    a = np.asarray(values, dtype=float)
    if len(a) == 0:
        return (float("nan"),) * 4
    mean = float(a.mean())
    if len(a) < 2:
        return mean, float("nan"), float("nan"), float("nan")
    sd = float(a.std(ddof=1)); se = sd / math.sqrt(len(a))
    return mean, sd, mean - 1.96 * se, mean + 1.96 * se


def category(record: dict) -> str:
    status = record["root_end_status"]
    pattern = record["terminal_pattern"]
    if status == "lost": return "lost"
    if status == "fixed": return "fixed"
    if pattern == "1010": return "terminal_1010"
    if pattern == "0101": return "terminal_0101"
    if record["is_2_2_pattern"].lower() == "true": return "other_2_2"
    if status == "segregating": return "segregating_non_2_2"
    return "remaining"


def bootstrap_indices(n: int, count: int, rng: np.random.Generator):
    return rng.integers(0, n, size=(count, n))


def process_cell(config_path: Path, root: Path, seed_offset: int):
    config = yaml.safe_load(config_path.read_text())
    raw = root / config["output"]["directory"]
    records = [r for r in read_rows(raw / "replicate_summary.csv") if r["accepted"].lower() == "true"]
    if not records:
        raise RuntimeError(f"No accepted records in {raw}")
    meta = config["experiment_metadata"]
    deltas = np.array([float(r["q2"]) - float(r["q1"]) for r in records])
    q_s = np.array([float(r["q1"]) for r in records]); q_a = np.array([float(r["q2"]) for r in records]); q_o = np.array([float(r["q3"]) for r in records])
    categories = np.array([category(r) for r in records])
    rng = np.random.default_rng(BOOT_SEED + seed_offset)
    indices = bootstrap_indices(len(records), BOOTSTRAPS, rng)
    boot_delta = deltas[indices].mean(axis=1)
    delta_mean, delta_sd, delta_low, delta_high = mean_ci(deltas)
    persistence = np.array([r["persisted_to_target"].lower() == "true" for r in records], dtype=float)
    root_persist = np.array([r["root_end_status"] == "segregating" for r in records], dtype=float)
    two_two = np.array([r["is_2_2_pattern"].lower() == "true" for r in records], dtype=float)
    pattern_1010 = np.array([r["terminal_pattern"] == "1010" for r in records], dtype=float)
    pattern_0101 = np.array([r["terminal_pattern"] == "0101" for r in records], dtype=float)
    summary = {"cell_id": meta["cell_id"], "Ne": config["species_tree"]["default_effective_population_size"], "tau": meta["internal_branch_coalescent_units"], "origin_fraction": meta["origin_fraction_of_root_extension"], "initial_copy_fraction": meta["initial_copy_fraction"], "cross_fraction": meta["cross_arrangement_fraction"], "short_branch_control": meta["short_branch_control"], "n_histories": len(records), "loci_per_history": config["experiment"]["loci_per_replicate"], "root_persistence": float(root_persist.mean()), "daughter_AB_persistence": float(persistence.mean()), "terminal_2_2_frequency": float(two_two.mean()), "terminal_1010_frequency": float(pattern_1010.mean()), "terminal_0101_frequency": float(pattern_0101.mean()), "terminal_unordered_2_2_frequency": float((two_two).mean()), "q_species": float(q_s.mean()), "q_species_sd_history": float(q_s.std(ddof=1)), "q_alt": float(q_a.mean()), "q_alt_sd_history": float(q_a.std(ddof=1)), "q_other": float(q_o.mean()), "q_other_sd_history": float(q_o.std(ddof=1)), "delta": delta_mean, "delta_sd_history": delta_sd, "delta_ci95_low": delta_low, "delta_ci95_high": delta_high, "delta_bootstrap_ci95_low": float(np.quantile(boot_delta, .025)), "delta_bootstrap_ci95_high": float(np.quantile(boot_delta, .975)), "positive_delta_history_fraction": float((deltas > 0).mean()), "config": str(config_path.relative_to(root)), "raw_directory": str(raw.relative_to(root))}
    persistence_rows = {"cell_id": meta["cell_id"], "Ne": config["species_tree"]["default_effective_population_size"], "tau": meta["internal_branch_coalescent_units"], "origin_fraction": meta["origin_fraction_of_root_extension"], "initial_copy_fraction": meta["initial_copy_fraction"], "cross_fraction": meta["cross_arrangement_fraction"], "metric": "root_persistence", "estimate": float(root_persist.mean()), "ci95_low": float(np.quantile(root_persist[indices].mean(axis=1), .025)), "ci95_high": float(np.quantile(root_persist[indices].mean(axis=1), .975)), "n_histories": len(records)}
    persistence_rows_all = [persistence_rows]
    for name, values in [("daughter_AB_persistence", persistence), ("terminal_2_2_frequency", two_two), ("terminal_1010_frequency", pattern_1010), ("terminal_0101_frequency", pattern_0101)]:
        boot = values[indices].mean(axis=1)
        persistence_rows_all.append({**persistence_rows, "metric": name, "estimate": float(values.mean()), "ci95_low": float(np.quantile(boot, .025)), "ci95_high": float(np.quantile(boot, .975))})
    decomposition = []
    for name in ["lost", "fixed", "terminal_1010", "terminal_0101", "other_2_2", "segregating_non_2_2", "remaining"]:
        mask = categories == name; weight = float(mask.mean())
        conditional_delta = float(deltas[mask].mean()) if mask.any() else float("nan")
        contribution = weight * conditional_delta if mask.any() else 0.0
        boot_contribution = []
        for sample in indices:
            sampled_cat = categories[sample]; sampled_delta = deltas[sample]
            present = sampled_cat == name
            boot_contribution.append(float(present.mean() * sampled_delta[present].mean()) if present.any() else 0.0)
        decomposition.append({"cell_id": meta["cell_id"], "category": name, "n_histories": int(mask.sum()), "weight": weight, "conditional_delta": conditional_delta, "conditional_q_species": float(q_s[mask].mean()) if mask.any() else float("nan"), "conditional_q_alt": float(q_a[mask].mean()) if mask.any() else float("nan"), "contribution": contribution, "contribution_bootstrap_ci95_low": float(np.quantile(boot_contribution, .025)), "contribution_bootstrap_ci95_high": float(np.quantile(boot_contribution, .975)), "config": str(config_path.relative_to(root))})
    summary["decomposition_sum"] = float(sum(r["contribution"] for r in decomposition))
    summary["decomposition_difference"] = summary["decomposition_sum"] - summary["delta"]
    return summary, persistence_rows_all, decomposition


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5]); args = ap.parse_args(); root = args.repo_root.resolve()
    configs = sorted((root / "simulations/quartet/mechanistic/datasets/configs/pilot2").glob("*.yaml"))
    summaries = []; persistence = []; decomposition = []
    for i, path in enumerate(configs):
        s, p, d = process_cell(path, root, i * 1009); summaries.append(s); persistence.extend(p); decomposition.extend(d)
    processed = root / "simulations/quartet/mechanistic/datasets/processed"
    write_tsv(processed / "mechanistic_pilot2_summary.tsv", summaries)
    write_tsv(processed / "mechanistic_pilot2_persistence.tsv", persistence)
    write_tsv(processed / "mechanistic_pilot2_bias_decomposition.tsv", decomposition)
    # Suppression comparisons and MSC expectation differences.
    by_key = {}
    for row in summaries:
        key = (row["Ne"], row["tau"], row["origin_fraction"], row["initial_copy_fraction"], row["short_branch_control"])
        by_key.setdefault(key, {})[float(row["cross_fraction"])] = row
    comparison = []
    for key, vals in by_key.items():
        if .1 not in vals or 1.0 not in vals: continue
        s, u = vals[.1], vals[1.0]
        tau = float(s["tau"]); expected_s = 1 - (2/3) * math.exp(-2*tau); expected_a = (1/3) * math.exp(-2*tau)
        delta_diff = s["delta"] - u["delta"]; delta_se = math.sqrt((s["delta_sd_history"] ** 2 / s["n_histories"]) + (u["delta_sd_history"] ** 2 / u["n_histories"])); qs_diff = s["q_species"] - u["q_species"]; qs_se = math.sqrt((s["q_species_sd_history"] ** 2 / s["n_histories"]) + (u["q_species_sd_history"] ** 2 / u["n_histories"]))
        comparison.append({"Ne": key[0], "tau": key[1], "origin_fraction": key[2], "initial_copy_fraction": key[3], "short_branch_control": key[4], "suppressed_cell": s["cell_id"], "unsuppressed_cell": u["cell_id"], "delta_suppressed": s["delta"], "delta_unsuppressed": u["delta"], "delta_difference_suppressed_minus_unsuppressed": delta_diff, "delta_difference_ci95_low": delta_diff - 1.96 * delta_se, "delta_difference_ci95_high": delta_diff + 1.96 * delta_se, "q_species_difference_suppressed_minus_unsuppressed": qs_diff, "q_species_difference_ci95_low": qs_diff - 1.96 * qs_se, "q_species_difference_ci95_high": qs_diff + 1.96 * qs_se, "suppressed_persistence": s["daughter_AB_persistence"], "unsuppressed_persistence": u["daughter_AB_persistence"], "msc_expected_q_species": expected_s, "msc_expected_q_alt": expected_a, "unsuppressed_q_species_minus_msc": u["q_species"] - expected_s, "unsuppressed_q_alt_minus_msc": u["q_alt"] - expected_a})
    write_tsv(processed / "mechanistic_pilot2_suppression_comparison.tsv", comparison)
    (root / "simulations/quartet/mechanistic/analysis/results/pilot2_processing_summary.json").write_text(json.dumps({"cells": len(summaries), "expected_primary_cells": 32, "short_branch_controls": sum(bool(x["short_branch_control"]) for x in summaries), "bootstrap_replicates": BOOTSTRAPS, "max_decomposition_difference": max(abs(float(x["decomposition_difference"])) for x in summaries), "all_decompositions_reconstruct": all(abs(float(x["decomposition_difference"])) < 1e-12 for x in summaries)}, indent=2) + "\n")
    print(f"processed {len(summaries)} pilot2 cells")


if __name__ == "__main__": main()
