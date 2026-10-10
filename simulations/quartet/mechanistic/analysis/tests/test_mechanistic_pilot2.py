#!/usr/bin/env python3
"""Regression tests for Experiment 2B outputs and identities."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import yaml
from msrcsim.structured_coalescent import quartet_index


ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "simulations/quartet/mechanistic"


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    assert quartet_index([("1", "2"), ("3", "4")]) == 0
    assert quartet_index([("1", "3"), ("2", "4")]) == 1
    assert quartet_index([("1", "4"), ("2", "3")]) == 2
    configs = sorted((BASE / "datasets/configs/pilot2").glob("*.yaml"))
    assert len(configs) == 36
    seeds = []
    for path in configs:
        cfg = yaml.safe_load(path.read_text()); seeds.append(int(cfg["seed"]))
        assert cfg["conditioning"]["mode"] == "none"
        assert cfg["experiment"]["persistence_target_branches"] == ["A", "B"]
        assert cfg["species_tree"]["newick"].count("1:") == 1 and all(f"{t}:" in cfg["species_tree"]["newick"] for t in ("1", "2", "3", "4"))
        tau = float(cfg["experiment_metadata"]["internal_branch_coalescent_units"])
        assert abs(float(cfg["experiment_metadata"]["internal_branch_generations"]) / (2 * cfg["species_tree"]["default_effective_population_size"]) - tau) < 1e-12
        raw = ROOT / cfg["output"]["directory"]
        for name in ("replicate_summary.csv", "prevalence_summary.json", "config.resolved.yaml", "run_metadata.json"):
            assert (raw / name).is_file(), (path, name)
        metadata = json.loads((raw / "run_metadata.json").read_text())
        assert metadata["config_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert int(metadata["return_code"]) == 0
    assert len(seeds) == len(set(seeds))

    summary = rows(BASE / "datasets/processed/mechanistic_pilot2_summary.tsv")
    persistence = rows(BASE / "datasets/processed/mechanistic_pilot2_persistence.tsv")
    decomposition = rows(BASE / "datasets/processed/mechanistic_pilot2_bias_decomposition.tsv")
    comparison = rows(BASE / "datasets/processed/mechanistic_pilot2_suppression_comparison.tsv")
    assert len(summary) == 36 and len(comparison) == 16
    for row in summary:
        assert int(row["n_histories"]) == 100
        assert abs(sum(float(row[k]) for k in ("q_species", "q_alt", "q_other")) - 1) < 1e-12
        assert abs(float(row["decomposition_sum"]) - float(row["delta"])) < 1e-12
        assert abs(float(row["decomposition_difference"])) < 1e-12
    by_cell = {}
    for row in decomposition:
        by_cell.setdefault(row["cell_id"], []).append(row)
    expected_categories = {"lost", "fixed", "terminal_1010", "terminal_0101", "other_2_2", "segregating_non_2_2", "remaining"}
    for cell, parts in by_cell.items():
        assert {p["category"] for p in parts} == expected_categories
        assert sum(int(p["n_histories"]) for p in parts) == 100
    assert all(float(row["delta"]) < 0 for row in summary)
    assert any(float(row["unsuppressed_q_species_minus_msc"]) != 0 for row in comparison)
    assert all(float(row["ci95_low"]) <= float(row["ci95_high"]) for row in persistence)
    print("mechanistic pilot2 regression tests passed")


if __name__ == "__main__":
    main()
