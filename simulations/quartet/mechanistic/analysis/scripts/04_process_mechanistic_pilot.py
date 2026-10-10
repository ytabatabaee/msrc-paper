#!/usr/bin/env python3
"""Process independent-history mechanistic pilot outputs and MSC controls."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import yaml


def mean_ci(values):
    values = [float(x) for x in values]
    n = len(values)
    mean = sum(values) / n if n else float("nan")
    if n < 2:
        return mean, float("nan"), float("nan"), float("nan")
    sd = math.sqrt(sum((x - mean) ** 2 for x in values) / (n - 1))
    se = sd / math.sqrt(n)
    return mean, sd, mean - 1.96 * se, mean + 1.96 * se


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def process_pilot(root: Path):
    rows = []
    pattern_rows = []
    config_root = root / "simulations/quartet/mechanistic/datasets/configs/pilot"
    for config_path in sorted(config_root.glob("*.yaml")):
        config = yaml.safe_load(config_path.read_text())
        raw = root / config["output"]["directory"]
        with (raw / "replicate_summary.csv").open() as handle:
            records = list(csv.DictReader(handle))
        accepted = [r for r in records if r["accepted"].lower() == "true"]
        if not accepted:
            continue
        meta = config["experiment_metadata"]
        q_values = {key: [float(r[key]) for r in accepted] for key in ("q1", "q2", "q3")}
        deltas = [q2 - q1 for q1, q2 in zip(q_values["q1"], q_values["q2"])]
        q_stats = {key: mean_ci(values) for key, values in q_values.items()}
        delta_mean, delta_sd, delta_low, delta_high = mean_ci(deltas)
        patterns = [r["terminal_pattern"] for r in accepted]
        pattern_counts = {pattern: patterns.count(pattern) for pattern in sorted(set(patterns))}
        pattern_1010 = [r for r in accepted if r["terminal_pattern"] == "1010"]
        conditional_q_alt = mean_ci([float(r["q2"]) for r in pattern_1010]) if pattern_1010 else (float("nan"), float("nan"), float("nan"), float("nan"))
        persistent = sum(r["persisted_to_target"].lower() == "true" for r in accepted) / len(accepted)
        two_two = sum(r["is_2_2_pattern"].lower() == "true" for r in accepted) / len(accepted)
        row = {
            "cell_id": meta["cell_id"], "seed": config["seed"], "Ne": config["species_tree"]["default_effective_population_size"], "internal_branch_coalescent_units": meta["internal_branch_coalescent_units"], "internal_branch_generations": meta["internal_branch_generations"], "origin_fraction_of_root_extension": meta["origin_fraction_of_root_extension"], "origin_time_from_branch_start": config["rearrangement"]["origin_time_from_branch_start"], "initial_copy_count": config["rearrangement"]["initial_copy_count"], "initial_copy_fraction": meta["initial_copy_fraction"], "cross_arrangement_fraction": meta["cross_arrangement_fraction"], "baseline_rate": config["recombination"]["baseline_rate"], "independent_replicates": len(accepted), "loci_per_replicate": config["experiment"]["loci_per_replicate"], "persistence_probability": persistent, "terminal_1010_frequency": patterns.count("1010") / len(patterns), "terminal_two_two_frequency": two_two, "affected_locus_fraction": "NA", "affected_locus_note": "not emitted by the mechanistic replicate API; terminal pattern and history persistence are reported instead", "q_species": q_stats["q1"][0], "q_species_sd_between_replicates": q_stats["q1"][1], "q_species_ci95_low": q_stats["q1"][2], "q_species_ci95_high": q_stats["q1"][3], "q_alt": q_stats["q2"][0], "q_alt_sd_between_replicates": q_stats["q2"][1], "q_alt_ci95_low": q_stats["q2"][2], "q_alt_ci95_high": q_stats["q2"][3], "q_other": q_stats["q3"][0], "q_other_sd_between_replicates": q_stats["q3"][1], "q_other_ci95_low": q_stats["q3"][2], "q_other_ci95_high": q_stats["q3"][3], "delta_alt_species": delta_mean, "delta_sd_between_replicates": delta_sd, "delta_ci95_low": delta_low, "delta_ci95_high": delta_high, "positive_delta_replicate_fraction": sum(x > 0 for x in deltas) / len(deltas), "positive_delta_supported": delta_low > 0, "q_alt_given_1010": conditional_q_alt[0], "q_alt_given_1010_ci95_low": conditional_q_alt[2], "q_alt_given_1010_ci95_high": conditional_q_alt[3], "terminal_pattern_counts": json.dumps(pattern_counts, sort_keys=True), "config": str(config_path.relative_to(root)), "raw_directory": str(raw.relative_to(root)),
        }
        rows.append(row)
        for pattern, count in pattern_counts.items():
            subset = [r for r in accepted if r["terminal_pattern"] == pattern]
            pattern_rows.append({"cell_id": meta["cell_id"], "terminal_pattern": pattern, "count": count, "frequency": count / len(accepted), "q_species": mean_ci([float(r["q1"]) for r in subset])[0], "q_alt": mean_ci([float(r["q2"]) for r in subset])[0], "q_other": mean_ci([float(r["q3"]) for r in subset])[0], "n_replicates": len(subset)})
    fields = list(rows[0].keys()) if rows else []
    write_tsv(root / "simulations/quartet/mechanistic/datasets/processed/mechanistic_quartet_pilot.tsv", fields, rows)
    write_tsv(root / "simulations/quartet/mechanistic/datasets/processed/mechanistic_pattern_contributions.tsv", ["cell_id", "terminal_pattern", "count", "frequency", "q_species", "q_alt", "q_other", "n_replicates"], pattern_rows)
    return rows


def process_controls(root: Path):
    rows = []
    for path in sorted((root / "simulations/quartet/mechanistic/datasets/raw/controls").glob("*/summary.json")):
        summary = json.loads(path.read_text())
        config = yaml.safe_load((path.parent / "config.resolved.yaml").read_text())
        tau = float(config["control_metadata"].get("unrooted_internal_branch_coalescent_units", 2 * float(config["control_metadata"].get("internal_branch_coalescent_units", 0.0)))); n = int(summary["n_loci"])
        observed = summary["probabilities"]; expected = summary["expected_probabilities"]
        se = [math.sqrt(p * (1 - p) / n) for p in observed]
        rows.append({"cell_id": config["control_metadata"]["cell_id"], "Ne": config["species_tree"]["default_effective_population_size"], "internal_branch_coalescent_units": tau, "n_loci": n, "q_species": observed[0], "q_alt": observed[1], "q_other": observed[2], "expected_q_species": expected[0], "expected_q_alt": expected[1], "expected_q_other": expected[2], "z_species": (observed[0] - expected[0]) / se[0], "z_alt": (observed[1] - expected[1]) / se[1], "z_other": (observed[2] - expected[2]) / se[2], "max_abs_error": max(abs(o - e) for o, e in zip(observed, expected)), "raw_directory": str(path.parent.relative_to(root))})
    write_tsv(root / "simulations/quartet/mechanistic/datasets/processed/msc_control_validation.tsv", list(rows[0].keys()), rows)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5]); args = ap.parse_args()
    root = args.repo_root.resolve(); pilot = process_pilot(root); controls = process_controls(root)
    candidates = [r for r in pilot if r["positive_delta_supported"]]
    (root / "simulations/quartet/mechanistic/analysis/results/processing_summary.json").write_text(json.dumps({"pilot_cells": len(pilot), "msc_controls": len(controls), "supported_positive_delta_cells": len(candidates), "pilot_cells_with_positive_point_estimate": sum(r["delta_alt_species"] > 0 for r in pilot)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
