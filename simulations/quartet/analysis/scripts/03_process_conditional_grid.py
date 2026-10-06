#!/usr/bin/env python3
"""Normalize msrc-sim conditional outputs and create validation summaries."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import yaml


TOPOLOGIES = [("q_species", "12|34"), ("q_alt", "13|24"), ("q_other", "14|23")]
FIELDS = ["duration", "m", "m01", "m10", "lambda0", "lambda1", "configuration", "seed", "n_loci", "empirical_q_species", "empirical_q_alt", "empirical_q_other", "exact_q_species", "exact_q_alt", "exact_q_other", "empirical_delta_alt_species", "exact_delta_alt_species", "abs_error_q_species", "abs_error_q_alt", "abs_error_q_other", "mc_se_q_species", "mc_se_q_alt", "mc_se_q_other", "alt_dominant_empirical", "alt_dominant_exact"]


def se(p: float, n: int) -> float:
    return math.sqrt(max(p * (1.0 - p) / n, 0.0))


def rows(root: Path, smoke: bool = False):
    raw_root = root / "simulations/quartet/datasets/raw/conditional_grid" / ("smoke" if smoke else "")
    config_root = root / "simulations/quartet/datasets/configs/conditional_grid" / ("smoke" if smoke else "")
    for config_path in sorted(config_root.glob("*.yaml")):
        c = yaml.safe_load(config_path.read_text())
        raw_dir = root / c["output"]["directory"]
        summary = json.loads((raw_dir / "summary.json").read_text())
        counts = summary["counts"]
        n = int(c["num_loci"])
        empirical = summary.get("empirical_probabilities", summary.get("empirical"))
        exact = summary["exact_probabilities"] if "exact_probabilities" in summary else summary["exact"]
        duration = float(c["structured_interval"]["duration"])
        m = float(c["structured_interval"]["migration"]["m01"])
        out = {"duration": duration, "m": m, "m01": m, "m10": float(c["structured_interval"]["migration"]["m10"]), "lambda0": 1.0, "lambda1": 1.0, "configuration": "1010", "seed": int(c["seed"]), "n_loci": n}
        for i, (name, _) in enumerate(TOPOLOGIES):
            out[f"empirical_{name}"] = float(empirical[i])
            out[f"exact_{name}"] = float(exact[i])
            out[f"abs_error_{name}"] = abs(float(empirical[i]) - float(exact[i]))
            out[f"mc_se_{name}"] = se(float(empirical[i]), n)
        out["empirical_delta_alt_species"] = out["empirical_q_alt"] - out["empirical_q_species"]
        out["exact_delta_alt_species"] = out["exact_q_alt"] - out["exact_q_species"]
        out["alt_dominant_empirical"] = out["empirical_q_alt"] > out["empirical_q_species"]
        out["alt_dominant_exact"] = out["exact_q_alt"] > out["exact_q_species"]
        yield out


def write_tsv(path: Path, fieldnames, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(data)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[4])
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    data = list(rows(args.repo_root.resolve(), args.smoke))
    if args.smoke:
        return
    result_root = args.repo_root / "simulations/quartet/analysis/results"
    processed = args.repo_root / "simulations/quartet/datasets/processed/conditional_quartet_grid.tsv"
    write_tsv(processed, FIELDS, data)
    validation = []
    for row in data:
        for name, topology in TOPOLOGIES:
            error = row[f"empirical_{name}"] - row[f"exact_{name}"]
            mcse = row[f"mc_se_{name}"]
            validation.append({"duration": row["duration"], "m": row["m"], "topology": topology, "empirical_probability": row[f"empirical_{name}"], "exact_probability": row[f"exact_{name}"], "absolute_error": abs(error), "monte_carlo_se": mcse, "z_error": error / mcse if mcse else 0.0})
    write_tsv(result_root / "conditional_exact_validation.tsv", ["duration", "m", "topology", "empirical_probability", "exact_probability", "absolute_error", "monte_carlo_se", "z_error"], validation)
    boundary = []
    for duration in sorted({r["duration"] for r in data}):
        subset = [r for r in data if r["duration"] == duration]
        positive = [r for r in subset if r["exact_dominant_alt"]] if "exact_dominant_alt" in subset[0] else [r for r in subset if r["alt_dominant_exact"]]
        maximum = max(subset, key=lambda r: r["exact_delta_alt_species"])
        boundary.append({"duration": duration, "smallest_tested_m_with_q_alt_gt_q_species": min((r["m"] for r in positive), default=""), "largest_tested_m_with_q_alt_gt_q_species": max((r["m"] for r in positive), default=""), "maximum_delta_alt_species": maximum["exact_delta_alt_species"], "m_at_maximum_delta": maximum["m"], "q_species_at_maximum": maximum["exact_q_species"], "q_alt_at_maximum": maximum["exact_q_alt"], "q_other_at_maximum": maximum["exact_q_other"]})
    write_tsv(result_root / "conditional_dominance_boundary.tsv", ["duration", "smallest_tested_m_with_q_alt_gt_q_species", "largest_tested_m_with_q_alt_gt_q_species", "maximum_delta_alt_species", "m_at_maximum_delta", "q_species_at_maximum", "q_alt_at_maximum", "q_other_at_maximum"], boundary)
    stats = {"max_abs_error": max(r["absolute_error"] for r in validation), "median_abs_error": sorted(r["absolute_error"] for r in validation)[len(validation)//2], "fraction_within_2_se": sum(abs(r["z_error"]) <= 2 for r in validation) / len(validation), "fraction_within_3_se": sum(abs(r["z_error"]) <= 3 for r in validation) / len(validation)}
    (result_root / "conditional_validation_summary.json").write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    main()
