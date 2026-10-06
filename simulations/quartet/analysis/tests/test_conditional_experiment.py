from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import yaml
from msrcsim.analytic import TOPOLOGY_NAMES, TOPOLOGY_PAIRS

ROOT = Path(__file__).resolve().parents[4] / "simulations/quartet"
CONFIG_ROOT = ROOT / "datasets/configs/conditional_grid"


def test_topology_mapping_and_arrangement():
    assert tuple(TOPOLOGY_NAMES) == ("12|34", "13|24", "14|23")
    assert TOPOLOGY_PAIRS == {0: ((0, 1), (2, 3)), 1: ((0, 2), (1, 3)), 2: ((0, 3), (1, 2))}
    config = "1010"
    assert [i + 1 for i, x in enumerate(config) if x == "1"] == [1, 3]
    assert [i + 1 for i, x in enumerate(config) if x == "0"] == [2, 4]


def test_configs_seeds_and_schema():
    configs = sorted(CONFIG_ROOT.glob("*.yaml"))
    if not configs: pytest.skip("full grid configs not generated")
    assert len(configs) == 70
    seeds = []
    for path in configs:
        c = yaml.safe_load(path.read_text())
        seeds.append(c["seed"])
        assert c["mode"] == "conditional"
        assert c["structured_interval"]["configuration"] == "1010"
        assert c["structured_interval"]["migration"]["m01"] == c["structured_interval"]["migration"]["m10"]
        assert c["structured_interval"]["coalescence"] == {"lambda0": 1.0, "lambda1": 1.0}
        assert c["num_loci"] == 100000
    assert len(set(seeds)) == 70


def test_processed_probabilities_and_raw_match():
    processed = ROOT / "datasets/processed/conditional_quartet_grid.tsv"
    if not processed.exists(): pytest.skip("full grid not run")
    with processed.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 70
    for row in rows:
        assert abs(sum(float(row[f"exact_q_{x}"]) for x in ("species", "alt", "other")) - 1) < 1e-12
        assert abs(sum(float(row[f"empirical_q_{x}"]) for x in ("species", "alt", "other")) - 1) < 1e-12
    first = rows[0]
    raw = ROOT / "datasets/raw/conditional_grid/d00_m00/summary.json"
    summary = json.loads(raw.read_text())
    assert summary["counts"] == [32111, 36078, 31811]
    assert summary["exact_probabilities"] == [float(first[f"exact_q_{x}"]) for x in ("species", "alt", "other")]


def test_exact_delta_is_monotone_in_switching_rate():
    processed = ROOT / "datasets/processed/conditional_quartet_grid.tsv"
    if not processed.exists(): pytest.skip("full grid not run")
    with processed.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    for duration in sorted({float(r["duration"]) for r in rows}):
        values = [float(r["exact_delta_alt_species"]) for r in sorted((r for r in rows if float(r["duration"]) == duration), key=lambda r: float(r["m"]))]
        assert all(a >= b for a, b in zip(values, values[1:]))


def test_validation_errors_are_mc_consistent():
    validation = ROOT / "analysis/results/conditional_exact_validation.tsv"
    if not validation.exists(): pytest.skip("full grid not run")
    with validation.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 210
    assert sum(abs(float(r["z_error"])) <= 3 for r in rows) / len(rows) > 0.98


def test_repeated_raw_cells_are_complete():
    raw = ROOT / "datasets/raw/conditional_grid"
    if not raw.exists(): pytest.skip("full grid not run")
    cells = list(raw.glob("d*_m*"))
    if not cells: pytest.skip("full grid not run")
    assert len(cells) == 70
    for cell in cells:
        assert all((cell / x).exists() for x in ("summary.json", "quartet_probabilities.csv", "config.resolved.yaml", "run_metadata.json"))


def test_plot_source_is_processed_table():
    script = (ROOT / "analysis/scripts/04_plot_conditional_grid.py").read_text()
    assert "datasets/processed/conditional_quartet_grid.tsv" in script
