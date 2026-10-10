#!/usr/bin/env python3
"""Offline integrity tests for the saved mechanistic pilot."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import hashlib

import yaml


ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "simulations/quartet/mechanistic"


def rows(path, delimiter="\t"):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def main():
    configs = sorted((BASE / "datasets/configs/pilot").glob("*.yaml"))
    assert len(configs) == 72
    seeds = []
    for p in configs:
        cfg = yaml.safe_load(p.read_text())
        assert cfg["mode"] == "replicate_experiment"
        assert cfg["conditioning"]["mode"] == "none"
        tree_text = cfg["species_tree"]["newick"]
        assert all(f"{taxon}:" in tree_text for taxon in ("1", "2", "3", "4"))
        assert cfg["rearrangement"]["origin_branch"] == "ROOT"
        assert cfg["rearrangement"]["initial_copy_count"] >= 1
        seeds.append(int(cfg["seed"]))
        raw = BASE / "datasets/raw/pilot" / p.stem
        for required in ("replicate_summary.csv", "prevalence_summary.json", "config.resolved.yaml", "run_metadata.json"):
            assert (raw / required).is_file(), (p, required)
        meta = json.loads((raw / "run_metadata.json").read_text())
        assert meta["config_sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
    assert len(set(seeds)) == len(seeds)

    processed = rows(BASE / "datasets/processed/mechanistic_quartet_pilot.tsv")
    assert len(processed) == 72
    for r in processed:
        qs = [float(r[k]) for k in ("q_species", "q_alt", "q_other")]
        assert abs(sum(qs) - 1) < 1e-10
        assert abs(float(r["delta_alt_species"]) - (qs[1] - qs[0])) < 1e-10
        assert float(r["delta_ci95_low"]) <= float(r["delta_ci95_high"])

    controls = rows(BASE / "datasets/processed/msc_control_validation.tsv")
    assert len(controls) == 6
    for r in controls:
        assert max(abs(float(r[k])) for k in ("z_species", "z_alt", "z_other")) < 4

    for name in ("mechanistic_scenario", "mechanistic_unconditional_probabilities", "mechanistic_quartet_simplex", "mechanistic_quartet_bias_map", "mechanistic_conditional_vs_unconditional"):
        assert (BASE / "analysis/figures" / f"{name}.png").stat().st_size > 1000
        assert (BASE / "analysis/figures" / f"{name}.pdf").stat().st_size > 1000
    assert all(float(r["delta_alt_species"]) < 0 for r in processed)
    print("mechanistic pilot integrity tests passed")


if __name__ == "__main__":
    main()
