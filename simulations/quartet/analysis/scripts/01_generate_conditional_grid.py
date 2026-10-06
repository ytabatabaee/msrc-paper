#!/usr/bin/env python3
"""Generate frozen conditional quartet configurations."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import yaml

DURATIONS = [0.02, 0.05, 0.10, 0.25, 0.50, 1.00, 2.00]
M_VALUES = [0.000, 0.0025, 0.005, 0.010, 0.025, 0.050, 0.100, 0.250, 0.500, 1.000]
MASTER_SEED = 20261006
SMOKE_DURATIONS = [0.05, 0.50, 1.00]
SMOKE_M_VALUES = [0.000, 0.050, 1.000]


def cell_seed(duration_index: int, migration_index: int, master_seed: int = MASTER_SEED) -> int:
    payload = f"{master_seed}:duration:{duration_index}:migration:{migration_index}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % (2**32)


def cell_id(duration_index: int, migration_index: int) -> str:
    return f"d{duration_index:02d}_m{migration_index:02d}"


def write_config(path: Path, duration: float, m: float, seed: int, n_loci: int, raw_root: str) -> None:
    config = {
        "mode": "conditional",
        "seed": int(seed),
        "num_loci": int(n_loci),
        "structured_interval": {
            "duration": float(duration),
            "configuration": "1010",
            "migration": {"m01": float(m), "m10": float(m)},
            "coalescence": {"lambda0": 1.0, "lambda1": 1.0},
        },
        "output": {
            "directory": f"{raw_root}/{path.stem}",
            "record_resolved_config": True,
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[4])
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    durations = SMOKE_DURATIONS if args.smoke else DURATIONS
    m_values = SMOKE_M_VALUES if args.smoke else M_VALUES
    n_loci = 5000 if args.smoke else 100000
    config_root = args.repo_root / "simulations/quartet/datasets/configs/conditional_grid"
    raw_root = "simulations/quartet/datasets/raw/conditional_grid/smoke" if args.smoke else "simulations/quartet/datasets/raw/conditional_grid"
    for di, duration in enumerate(durations):
        for mi, m in enumerate(m_values):
            write_config(config_root / ("smoke/" if args.smoke else "") / f"{cell_id(di, mi)}.yaml", duration, m, cell_seed(di, mi), n_loci, raw_root)


if __name__ == "__main__":
    main()
