#!/usr/bin/env python3
"""Synthetic power analysis using the predeclared inside-vs-flank statistic."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
import random
import statistics
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "anopheles_2la"
RESULTS = EMPIRICAL_ROOT / "results" / "stage1b_synthetic"


def load_ingest():
    path = EMPIRICAL_ROOT / "scripts" / "10_ingest_local_trees.py"
    spec = importlib.util.spec_from_file_location("stage1b_ingest", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["stage1b_ingest"] = module
    sys.modules["stage1b_ingest_power"] = module
    spec.loader.exec_module(module)
    return module


def simulate_rows(m, rng: random.Random, n_quartets: int, blocks: int, f_out: float, f_in: float) -> list[dict[str, object]]:
    rows = []
    pos = m.LEFT_BP - 1_000_000
    for q in range(n_quartets):
        for region, p in [("left_flank", f_out), ("inversion_2La", f_in), ("right_flank", f_out)]:
            for b in range(blocks):
                rows.append(
                    {
                        "quartet_id": f"q{q}",
                        "start": pos,
                        "end": pos + 1,
                        "midpoint": pos,
                        "region": region,
                        "topology_class": "arrangement" if rng.random() < p else "species",
                        "usable": True,
                    }
                )
                pos += 2
    return rows


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_power(tsv: Path, pdf: Path) -> None:
    try:
        os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc_mplconfig")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    rows = []
    with tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fig, ax = plt.subplots(figsize=(7, 4))
    for f_in in sorted({r["f_in"] for r in rows}):
        subset = [r for r in rows if r["f_in"] == f_in and r["independent_quartets"] == "5" and r["f_out"] == "0.33"]
        subset.sort(key=lambda r: int(r["effective_blocks_per_region"]))
        ax.plot([int(r["effective_blocks_per_region"]) for r in subset], [float(r["power"]) for r in subset], marker="o", label=f"f_in={f_in}")
    ax.set_xlabel("effective blocks per region")
    ax.set_ylabel("power")
    ax.set_ylim(0, 1.02)
    ax.set_title("SYNTHETIC VALIDATION power analysis")
    ax.legend()
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synthetic-validation", action="store_true", required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=99)
    parser.add_argument("--out", type=Path, default=RESULTS / "power_analysis.tsv")
    parser.add_argument("--pdf", type=Path, default=RESULTS / "power_analysis.pdf")
    args = parser.parse_args()
    m = load_ingest()
    rng = random.Random(m.SEED)
    rows = []
    for n_quartets in [1, 3, 5, 10]:
        for blocks in [10, 25, 50, 100]:
            for f_out in [0.2, 0.33]:
                for f_in in [0.4, 0.6, 0.8]:
                    hits = 0
                    deltas = []
                    for _ in range(args.replicates):
                        sim = simulate_rows(m, rng, n_quartets, blocks, f_out, f_in)
                        test = m.circular_shift_test(sim, args.permutations, rng.randrange(10**9))
                        if float(test["observed_delta"]) > 0 and float(test["p_value"]) < 0.05:
                            hits += 1
                        deltas.append(float(test["observed_delta"]))
                    rows.append(
                        {
                            "SYNTHETIC_ONLY": "TRUE",
                            "independent_quartets": n_quartets,
                            "effective_blocks_per_region": blocks,
                            "f_out": f_out,
                            "f_in": f_in,
                            "replicates": args.replicates,
                            "permutations_per_replicate": args.permutations,
                            "power": hits / args.replicates,
                            "mean_delta": statistics.mean(deltas),
                        }
                    )
    write_rows(args.out, rows)
    plot_power(args.out, args.pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
