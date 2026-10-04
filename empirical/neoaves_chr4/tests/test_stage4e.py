import importlib.util
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


def fake_aster(path, body):
    path.write_text(body)
    path.chmod(0o755)
    return path


def write_tiny_manifest(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    tiny = data / "TINY.tre"
    tiny.write_text("((A,B),(C,D));\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "treatment\tscope\trelative_path\tn_loci\tn_chr4_loci\ttrees_with_fewer_than_four_taxa_omitted\tsha256\tdescription\n"
        f"TINY\tunit\tstage4e/TINY.tre\t1\t0\t\t{runner.sha256(tiny)}\ttiny fixture\n"
    )
    return data, manifest


def test_preflight_smoke_test_uses_tiny_synthetic_input(tmp_path):
    aster = fake_aster(
        tmp_path / "astral4",
        """#!/bin/sh
if [ "$1" = "-h" ]; then
  echo "fake astral4 help"
  exit 0
fi
while [ "$#" -gt 0 ]; do
  case "$1" in
    -i) input="$2"; shift 2 ;;
    -o) output="$2"; shift 2 ;;
    -u|-t|--length) shift 2 ;;
    *) shift ;;
  esac
done
case "$input" in
  *stage4e*.tre) echo "real Stage-4E input was used" >&2; exit 13 ;;
esac
grep -q "((A,B),(C,D));" "$input" || exit 14
printf '((A,B),(C,D));\\n' > "$output"
""",
    )
    output_root = tmp_path / "out"
    report = runner.run_smoke_test(aster, output_root)
    assert report["smoke_test_passed"] is True
    assert report["smoke_test_returncode"] == 0
    tiny_input = output_root / "aster_smoke_test/tiny_input.tre"
    assert tiny_input.read_text() == "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
    assert " -t 1 " in f" {' '.join(report['smoke_test_command'])} "
    assert "stage4e" not in str(tiny_input)


def test_preflight_smoke_test_failure_causes_preflight_failure(tmp_path):
    data, manifest = write_tiny_manifest(tmp_path)
    aster = fake_aster(
        tmp_path / "astral4",
        """#!/bin/sh
if [ "$1" = "-h" ]; then
  echo "fake astral4 help"
  exit 0
fi
exit 9
""",
    )
    args = type(
        "Args",
        (),
        {"aster_bin": aster, "output_root": tmp_path / "out", "input_manifest": manifest, "data_dir": data},
    )()
    try:
        runner.run_preflight(args)
    except RuntimeError as exc:
        assert "smoke test failed" in str(exc).lower()
    else:
        raise AssertionError("preflight should fail when the smoke test fails")
    payload = (tmp_path / "out/stage4e_preflight_runtime.json").read_text()
    assert '"smoke_test_passed": false' in payload


def test_preflight_smoke_test_records_runtime_provenance(tmp_path):
    data, manifest = write_tiny_manifest(tmp_path)
    aster = fake_aster(
        tmp_path / "astral4",
        """#!/bin/sh
if [ "$1" = "-h" ]; then
  echo "fake astral4 help"
  exit 0
fi
while [ "$#" -gt 0 ]; do
  case "$1" in
    -o) output="$2"; shift 2 ;;
    -i|-u|-t|--length) shift 2 ;;
    *) shift ;;
  esac
done
printf '((A,B),(C,D));\\n' > "$output"
""",
    )
    args = type(
        "Args",
        (),
        {"aster_bin": aster, "output_root": tmp_path / "out", "input_manifest": manifest, "data_dir": data},
    )()
    runner.run_preflight(args)
    payload = (tmp_path / "out/stage4e_preflight_runtime.json").read_text()
    assert '"smoke_test_passed": true' in payload
    assert '"smoke_test_returncode": 0' in payload
    assert "tiny_input.tre" in payload
    assert "stage4e/TINY.tre" not in payload


def test_stage4e_slurm_bridges2_headers():
    scripts = sorted((BASE / "cluster/stage4e_aster").glob("*.sbatch"))
    assert scripts
    for script in scripts:
        text = script.read_text()
        assert "#SBATCH --mem" not in text
        assert "#SBATCH -p RM-shared" in text
        assert "#SBATCH --cpus-per-task=32" in text
        assert "# #SBATCH -A YOUR_ALLOCATION" in text


def test_python_code_does_not_submit_slurm_or_run_large_inputs():
    for path in [
        BASE / "scripts/05_stage4e_prepare_aster.py",
        BASE / "scripts/05_stage4e_run_aster.py",
        BASE / "scripts/05_stage4e_collect_aster.py",
        BASE / "scripts/05_stage4e_analyze_aster.py",
    ]:
        text = path.read_text()
        assert "subprocess.run(['sbatch'" not in text
        assert 'subprocess.run(["sbatch"' not in text
    source = (BASE / "scripts/05_stage4e_run_aster.py").read_text()
    smoke_source = source.split("def run_smoke_test", 1)[1].split("def run_preflight", 1)[0]
    assert "T0.tre" not in smoke_source
    assert "stage4e" not in smoke_source
