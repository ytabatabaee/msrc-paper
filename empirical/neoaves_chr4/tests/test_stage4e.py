import importlib.util
import os
from pathlib import Path

BASE = Path(__file__).parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = load("stage4e_prepare", BASE / "scripts/05_stage4e_prepare_aster.py")
runner = load("stage4e_runner", BASE / "scripts/05_stage4e_run_aster.py")


def test_prepare_script_has_no_aster_subprocess_invocation():
    source = (BASE / "scripts/05_stage4e_prepare_aster.py").read_text()
    assert "ASTRAL4 inference was run" in source
    assert "subprocess.run" not in source
    assert "05_stage4e_run_aster.py" in source


def test_runner_defaults_to_dry_run(tmp_path, monkeypatch):
    input_path = tmp_path / "input.tre"
    input_path.write_text("((A,B),(C,D));\n")
    aster = tmp_path / "astral4"
    aster.write_text("#!/bin/sh\nexit 0\n")
    aster.chmod(0o755)
    out = tmp_path / "out"
    args = type(
        "Args",
        (),
        {
            "execute": False,
            "allow_cluster_execution": False,
            "output_dir": out,
            "aster_bin": aster,
            "input": input_path,
            "treatment": "TINY",
            "threads": 2,
        },
    )()
    runner.run_treatment(args)
    assert (out / "exit_code.txt").read_text() == "DRY_RUN\n"
    assert not (out / "tree.nwk").exists()


def test_real_execution_outside_slurm_raises(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(runner, "recognized_cluster_hostname", lambda: False)
    try:
        runner.require_cluster_execution(True, False)
    except RuntimeError as exc:
        assert "SLURM job" in str(exc)
    else:
        raise AssertionError("execution outside SLURM should fail")


def test_manifest_counts_if_available():
    path = BASE / "results/stage4d_locus_manifest.tsv"
    if not path.exists():
        return
    rows = prepare.read_tsv(path)
    counts = prepare.validate_manifest(rows)
    assert counts["n_loci"] == 63430
    assert counts["n_nonchr4_loci"] == 57168
    assert counts["n_pnas_outlier_loci"] == 1431
    assert counts["n_struct_stringent_loci"] == 495
    assert counts["n_struct_primary_loci"] == 835
    assert counts["n_struct_inclusive_loci"] == 1308


def test_block_representatives_and_random_scaffolds_are_deterministic():
    rows = [
        {"locus_id": "z", "chromosome": "chr4", "midpoint": "250000"},
        {"locus_id": "a", "chromosome": "chr4", "midpoint": "250000"},
        {"locus_id": "b", "chromosome": "chr4", "midpoint": "750000"},
        {"locus_id": "r", "chromosome": "chr4_random", "midpoint": "250000"},
    ]
    expected = {("chr4", 0): "a", ("chr4", 1): "b", ("chr4_random", 0): "r"}
    assert prepare.representatives(rows) == expected
    assert prepare.representatives(list(reversed(rows))) == expected


def test_random_replicate_seed_and_chr4_count(tmp_path):
    rows = [
        {"locus_id": "n", "chromosome": "chr5", "midpoint": "1", "is_chr4": "False"},
        {"locus_id": "a", "chromosome": "chr4", "midpoint": "1", "is_chr4": "True"},
        {"locus_id": "b", "chromosome": "chr4", "midpoint": "2", "is_chr4": "True"},
        {"locus_id": "c", "chromosome": "chr4", "midpoint": "3", "is_chr4": "True"},
    ]
    seed = prepare.random_seed("T_RANDOM_MATCHED_PNAS", 1)
    ids = prepare.random_replicate_ids(rows, 2, seed)
    assert seed == prepare.random_seed("T_RANDOM_MATCHED_PNAS", 1)
    assert len(set(ids) & {"a", "b", "c"}) == 2
    assert "n" in ids


def test_selection_code_avoids_topology_columns():
    source = (BASE / "scripts/05_stage4e_prepare_aster.py").read_text()
    for forbidden in ("q1", "q2", "q3", "topology", "support"):
        assert forbidden not in source.lower().split("def write_analysis_plan", 1)[0]


def test_stage4d_historical_paths_are_not_removed():
    assert (BASE / "results/stage4d").exists()
    assert (BASE / "results/stage4d_report.md").exists()
    assert (BASE / "results/stage4d_manifest.json").exists()
