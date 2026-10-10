#!/usr/bin/env python3
"""Generate frozen Experiment 2B persistence/refinement configurations."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import yaml


MASTER_SEED = 20261020
NE_VALUES = [20, 100]
TAU_VALUES = [0.5, 2.0]
ORIGIN_FRACTIONS = [0.25, 0.75]
COPY_FRACTIONS = [0.20, 0.50]
CROSS_FRACTIONS = [0.1, 1.0]
SHORT_BRANCH_CELLS = [(20, 0.25, 0.20, 0.1), (20, 0.75, 0.50, 0.1), (100, 0.25, 0.20, 0.1), (100, 0.75, 0.50, 0.1)]
PILOT2_REPLICATES = 100
PILOT2_LOCI = 20
SMOKE_REPLICATES = 3
SMOKE_LOCI = 10


def seed_for(*parts: object) -> int:
    key = ":".join([str(MASTER_SEED), *map(str, parts)]).encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "big") % (2**32)


def make_tree(ne: int, tau: float) -> tuple[str, int, int, int]:
    tip = 2 * ne
    internal = int(round(2 * ne * tau))
    root_extension = 2 * ne
    return f"((1:{tip},2:{tip})A:{internal},(3:{tip},4:{tip})B:{internal})ROOT;", tip, internal, root_extension


def cell_id(ne: int, tau: float, origin: float, copy: float, cross: float, short: bool = False) -> str:
    prefix = "short_" if short else ""
    return f"{prefix}ne{ne}_tau{tau:g}_orig{origin:g}_copy{copy:g}_cross{cross:g}".replace(".", "p")


def write_config(path: Path, *, ne: int, tau: float, origin: float, copy: float, cross: float, reps: int, loci: int, smoke: bool, short: bool) -> None:
    name = cell_id(ne, tau, origin, copy, cross, short)
    _newick, _tip, internal, root_extension = make_tree(ne, tau)
    origin_time = max(1, int(round(root_extension * origin)))
    total = 2 * ne
    copies = max(1, min(total, int(round(total * copy))))
    raw_dir = f"simulations/quartet/mechanistic/datasets/raw/{'pilot2_smoke' if smoke else 'pilot2'}/{name}"
    config = {
        "mode": "replicate_experiment",
        "seed": seed_for("pilot2", name),
        "population_process": {"model": "wright_fisher"},
        "species_tree": {
            "newick": _newick,
            "time_units": "generations",
            "root_extension": root_extension,
            "default_effective_population_size": ne,
        },
        "rearrangement": {
            "id": "inv_1",
            "type": "inversion",
            "origin_branch": "ROOT",
            "origin_time_from_branch_start": origin_time,
            "initial_copy_count": copies,
            "selection": {"model": "genic", "coefficient": 0.0},
        },
        "recombination": {"baseline_rate": 0.02, "effective_cross_arrangement_fraction": cross},
        "sampling": {"samples_per_species": 1},
        "conditioning": {"mode": "none"},
        "experiment": {
            "replicates": reps,
            "loci_per_replicate": loci,
            "asymmetry_threshold": 0.05,
            "persistence_target_branches": ["A", "B"],
        },
        "output": {"directory": raw_dir, "record_resolved_config": True},
        "experiment_metadata": {
            "experiment": "mechanistic_pilot2",
            "cell_id": name,
            "internal_branch_coalescent_units": tau,
            "internal_branch_generations": internal,
            "origin_fraction_of_root_extension": origin,
            "initial_copy_fraction": copy,
            "cross_arrangement_fraction": cross,
            "persistence_target_definition": "both daughter ancestral branches A and B end segregating",
            "unconditional": True,
            "short_branch_control": short,
            "smoke": smoke,
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5])
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    root = args.repo_root.resolve()
    directory = root / "simulations/quartet/mechanistic/datasets/configs" / ("pilot2_smoke" if args.smoke else "pilot2")
    reps, loci = (SMOKE_REPLICATES, SMOKE_LOCI) if args.smoke else (PILOT2_REPLICATES, PILOT2_LOCI)
    cells = [(ne, tau, origin, copy, cross, False) for ne in NE_VALUES for tau in TAU_VALUES for origin in ORIGIN_FRACTIONS for copy in COPY_FRACTIONS for cross in CROSS_FRACTIONS]
    cells.extend((ne, 0.1, origin, copy, cross, True) for ne, origin, copy, cross in SHORT_BRANCH_CELLS)
    for ne, tau, origin, copy, cross, short in cells:
        name = cell_id(ne, tau, origin, copy, cross, short)
        write_config(directory / f"{name}.yaml", ne=ne, tau=tau, origin=origin, copy=copy, cross=cross, reps=reps, loci=loci, smoke=args.smoke, short=short)
    print(f"generated {len(cells)} {'smoke' if args.smoke else 'pilot2'} configurations")


if __name__ == "__main__":
    main()
