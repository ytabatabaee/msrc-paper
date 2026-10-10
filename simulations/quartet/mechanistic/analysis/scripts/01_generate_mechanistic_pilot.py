#!/usr/bin/env python3
"""Generate frozen unconditional mechanistic WF pilot configs and MSC controls."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import yaml

NE_VALUES = [20, 100]
TAU_VALUES = [0.1, 0.5, 2.0]
ORIGIN_FRACTIONS = [0.25, 0.75]
COPY_FRACTIONS = [0.025, 0.20, 0.50]
CROSS_FRACTIONS = [1.0, 0.1]
MASTER_SEED = 20261010
PILOT_REPLICATES = 12
PILOT_LOCI = 20
SMOKE_REPLICATES = 3
SMOKE_LOCI = 10
CONTROL_LOCI = 20000


def seed_for(*parts: object) -> int:
    key = ":".join([str(MASTER_SEED), *map(str, parts)]).encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "big") % (2**32)


def make_tree(ne: int, tau: float) -> tuple[str, int, int, int]:
    tip = 2 * ne
    internal = int(round(2 * ne * tau))
    root_extension = 2 * ne
    newick = f"((1:{tip},2:{tip})A:{internal},(3:{tip},4:{tip})B:{internal})ROOT;"
    return newick, tip, internal, root_extension


def write_mechanistic(path: Path, *, cell_id: str, ne: int, tau: float, origin_fraction: float, copy_fraction: float, cross_fraction: float, reps: int, loci: int, smoke: bool) -> None:
    newick, _tip, internal, root_extension = make_tree(ne, tau)
    origin_time = max(1, int(round(root_extension * origin_fraction)))
    total = 2 * ne
    copies = max(1, min(total, int(round(total * copy_fraction))))
    raw_dir = f"simulations/quartet/mechanistic/datasets/raw/{'smoke' if smoke else 'pilot'}/{cell_id}"
    config = {
        "mode": "replicate_experiment",
        "seed": seed_for("mechanistic", cell_id),
        "population_process": {"model": "wright_fisher"},
        "species_tree": {"newick": newick, "time_units": "generations", "root_extension": root_extension, "default_effective_population_size": ne},
        "rearrangement": {"id": "inv_1", "type": "inversion", "origin_branch": "ROOT", "origin_time_from_branch_start": origin_time, "initial_copy_count": copies, "selection": {"model": "genic", "coefficient": 0.0}},
        "recombination": {"baseline_rate": 0.02, "effective_cross_arrangement_fraction": cross_fraction},
        "sampling": {"samples_per_species": 1},
        "conditioning": {"mode": "none"},
        "experiment": {"replicates": reps, "loci_per_replicate": loci, "asymmetry_threshold": 0.05, "persistence_target_branches": ["ROOT"]},
        "output": {"directory": raw_dir, "record_resolved_config": True},
        "experiment_metadata": {"cell_id": cell_id, "internal_branch_coalescent_units": tau, "internal_branch_generations": internal, "origin_fraction_of_root_extension": origin_fraction, "initial_copy_fraction": copy_fraction, "cross_arrangement_fraction": cross_fraction, "unconditional": True, "smoke": smoke},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def write_control(path: Path, *, cell_id: str, ne: int, tau: float, n_loci: int) -> None:
    newick, _tip, internal, root_extension = make_tree(ne, tau)
    config = {"mode": "ordinary_msc_control", "seed": seed_for("msc", cell_id), "n_loci": n_loci, "species_tree": {"newick": newick, "time_units": "generations", "root_extension": root_extension, "default_effective_population_size": ne}, "control_metadata": {"cell_id": cell_id, "daughter_branch_coalescent_units": tau, "unrooted_internal_branch_coalescent_units": 2 * tau, "internal_branch_generations_per_side": internal, "expected_formula": "qS=1-(2/3)exp(-t_unrooted), qA=qO=(1/3)exp(-t_unrooted)"}}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5])
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    root = args.repo_root.resolve()
    config_root = root / "simulations/quartet/mechanistic/datasets/configs"
    nreps, nloci = (SMOKE_REPLICATES, SMOKE_LOCI) if args.smoke else (PILOT_REPLICATES, PILOT_LOCI)
    for ne in NE_VALUES:
        for tau in TAU_VALUES:
            for origin_fraction in ORIGIN_FRACTIONS:
                for copy_fraction in COPY_FRACTIONS:
                    for cross_fraction in CROSS_FRACTIONS:
                        cell_id = f"ne{ne}_tau{tau:g}_orig{origin_fraction:g}_copy{copy_fraction:g}_cross{cross_fraction:g}".replace(".", "p")
                        write_mechanistic(config_root / ("smoke/" if args.smoke else "pilot/") / f"{cell_id}.yaml", cell_id=cell_id, ne=ne, tau=tau, origin_fraction=origin_fraction, copy_fraction=copy_fraction, cross_fraction=cross_fraction, reps=nreps, loci=nloci, smoke=args.smoke)
    for ne in NE_VALUES:
        for tau in TAU_VALUES:
            cell_id = f"ne{ne}_tau{tau:g}".replace(".", "p")
            out = f"simulations/quartet/mechanistic/datasets/raw/controls/{cell_id}"
            path = config_root / "controls" / f"{cell_id}.yaml"
            write_control(path, cell_id=cell_id, ne=ne, tau=tau, n_loci=CONTROL_LOCI if not args.smoke else 500)


if __name__ == "__main__":
    main()
