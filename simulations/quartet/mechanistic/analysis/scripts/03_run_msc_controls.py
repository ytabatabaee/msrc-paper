#!/usr/bin/env python3
"""Run ordinary MSC controls through msrc-sim's public genealogy API."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import yaml

from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_msc_genealogy


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5])
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    root = args.repo_root.resolve()
    for config_path in sorted((root / "simulations/quartet/mechanistic/datasets/configs/controls").glob("*.yaml")):
        config = yaml.safe_load(config_path.read_text())
        out = root / "simulations/quartet/mechanistic/datasets/raw/controls" / config["control_metadata"]["cell_id"]
        config_sha256 = hashlib.sha256(config_path.read_bytes()).hexdigest()
        if (out / "summary.json").is_file() and (out / "config.resolved.yaml").is_file():
            existing = json.loads((out / "summary.json").read_text())
            metadata = out / "run_metadata.json"
            existing_sha256 = ""
            if metadata.is_file():
                existing_sha256 = json.loads(metadata.read_text()).get("config_sha256", "")
            if (int(existing.get("n_loci", -1)) == int(config["n_loci"])
                    and existing_sha256 == config_sha256):
                continue
        out.mkdir(parents=True, exist_ok=True)
        st = config["species_tree"]
        tree = SpeciesTree(st["newick"], st["default_effective_population_size"], st["root_extension"])
        rng = np.random.default_rng(int(config["seed"]))
        counts = np.zeros(3, dtype=int)
        start = time.perf_counter()
        for locus_id in range(int(config["n_loci"])):
            counts[simulate_msc_genealogy(locus_id, tree, rng, record_events=False).topology_index] += 1
        n = int(config["n_loci"])
        tau = float(config["control_metadata"]["unrooted_internal_branch_coalescent_units"])
        expected = [1 - (2 / 3) * np.exp(-tau), (1 / 3) * np.exp(-tau), (1 / 3) * np.exp(-tau)]
        with (out / "quartet_counts.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t"); writer.writerow(["topology", "count", "probability", "expected_probability"])
            for name, count, obs, exp in zip(("12|34", "13|24", "14|23"), counts, counts / n, expected): writer.writerow([name, int(count), float(obs), float(exp)])
        (out / "summary.json").write_text(json.dumps({"mode": "ordinary_msc_control", "seed": int(config["seed"]), "n_loci": n, "counts": counts.tolist(), "probabilities": (counts / n).tolist(), "expected_probabilities": expected, "runtime_seconds": time.perf_counter() - start}, indent=2) + "\n")
        (out / "config.resolved.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        (out / "run_metadata.json").write_text(json.dumps({"config": str(config_path.relative_to(root)), "config_sha256": config_sha256, "api": "msrcsim.structured_coalescent.simulate_msc_genealogy"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
