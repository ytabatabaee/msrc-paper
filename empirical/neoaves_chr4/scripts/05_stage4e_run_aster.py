#!/usr/bin/env python3
"""Dry-run-safe ASTER/ASTRAL4 runner for Stage 4E cluster jobs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def astra_command(aster_bin: Path, input_path: Path, output_path: Path, threads: int) -> list[str]:
    return [
        str(aster_bin),
        "-i",
        str(input_path),
        "-o",
        str(output_path),
        "-u",
        "2",
        "-t",
        str(threads),
        "--length",
        "CULength",
    ]


def running_under_slurm() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID"))


def recognized_cluster_hostname() -> bool:
    host = platform.node().lower()
    return any(token in host for token in ("bridges", "psc.edu", "login", "compute"))


def require_cluster_execution(execute: bool, allow_cluster_execution: bool) -> None:
    if not execute:
        return
    if running_under_slurm():
        return
    if allow_cluster_execution and recognized_cluster_hostname():
        return
    raise RuntimeError("Large Stage-4E ASTER analysis may only be executed inside a SLURM job.")


def help_text(aster_bin: Path) -> str:
    result = subprocess.run([str(aster_bin), "-h"], capture_output=True, text=True, check=False)
    return (result.stdout + result.stderr).strip()


def validate_four_taxon_newick(path: Path) -> None:
    import dendropy

    text = path.read_text().strip()
    if not text:
        raise RuntimeError("Smoke-test output tree is empty")
    if not text.endswith(";"):
        raise RuntimeError("Smoke-test output tree does not end in ';'")
    tree = dendropy.Tree.get(data=text, schema="newick", rooting="force-unrooted", preserve_underscores=True)
    taxa = {node.taxon.label for node in tree.leaf_node_iter()}
    if taxa != {"A", "B", "C", "D"}:
        raise RuntimeError(f"Smoke-test output has unexpected taxa: {sorted(taxa)}")


def run_smoke_test(aster_bin: Path, output_root: Path) -> dict[str, object]:
    smoke_dir = output_root / "aster_smoke_test"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    tiny_input = smoke_dir / "tiny_input.tre"
    tiny_output = smoke_dir / "tiny_output.nwk"
    stdout_path = smoke_dir / "stdout.log"
    stderr_path = smoke_dir / "stderr.log"
    tiny_input.write_text("((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n")
    cmd = astra_command(aster_bin, tiny_input, tiny_output, 1)
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        result = subprocess.run(cmd, stdout=stdout, stderr=stderr, check=False)
    passed = False
    error = ""
    try:
        if result.returncode != 0:
            raise RuntimeError(f"Smoke-test return code was {result.returncode}")
        if not tiny_output.exists():
            raise RuntimeError("Smoke-test output tree was not created")
        validate_four_taxon_newick(tiny_output)
        passed = True
    except Exception as exc:
        error = str(exc)
    return {
        "smoke_test_command": cmd,
        "smoke_test_returncode": result.returncode,
        "smoke_test_passed": passed,
        "smoke_test_tree": str(tiny_output),
        "smoke_test_stdout": str(stdout_path),
        "smoke_test_stderr": str(stderr_path),
        "smoke_test_error": error,
    }


def run_preflight(args: argparse.Namespace) -> None:
    if not args.aster_bin.exists():
        raise RuntimeError(f"ASTER_BIN does not exist: {args.aster_bin}")
    if not os.access(args.aster_bin, os.X_OK):
        raise RuntimeError(f"ASTER_BIN is not executable: {args.aster_bin}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(args.output_root)
    rows = []
    with args.input_manifest.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            path = args.data_dir / Path(row["relative_path"]).name
            if not path.exists():
                raise RuntimeError(f"Missing input: {path}")
            if row["sha256"] and sha256(path) != row["sha256"]:
                raise RuntimeError(f"Input checksum mismatch: {path}")
            rows.append(row)
    smoke = run_smoke_test(args.aster_bin, args.output_root)
    report = {
        "status": "preflight_passed",
        "large_inference_run": False,
        "aster_bin": str(args.aster_bin.resolve()),
        "aster_bin_sha256": sha256(args.aster_bin),
        "aster_help": help_text(args.aster_bin),
        "hostname": platform.node(),
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "memory_request": os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU"),
        "output_root": str(args.output_root),
        "output_free_bytes": usage.free,
        "n_inputs_checked": len(rows),
        **smoke,
    }
    if not smoke["smoke_test_passed"]:
        report["status"] = "preflight_failed"
    out = args.output_root / "stage4e_preflight_runtime.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    if not smoke["smoke_test_passed"]:
        raise RuntimeError(f"ASTER smoke test failed; see {smoke['smoke_test_stderr']}: {smoke['smoke_test_error']}")
    print(json.dumps({"status": "preflight_passed", "report": str(out)}, indent=2))


def run_treatment(args: argparse.Namespace) -> None:
    require_cluster_execution(args.execute, args.allow_cluster_execution)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK") or args.threads)
    tree = args.output_dir / "tree.nwk"
    cmd = astra_command(args.aster_bin, args.input, tree, threads)
    (args.output_dir / "command.txt").write_text(" ".join(cmd) + "\n")
    provenance = {
        "treatment": args.treatment,
        "dry_run": not args.execute,
        "command": cmd,
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "aster_bin": str(args.aster_bin.resolve()),
        "aster_bin_sha256": sha256(args.aster_bin) if args.aster_bin.exists() else None,
        "hostname": platform.node(),
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "cpu_count": threads,
        "memory_request": os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU"),
    }
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if not args.execute:
        (args.output_dir / "exit_code.txt").write_text("DRY_RUN\n")
        print("DRY RUN: " + " ".join(cmd))
        return
    with (args.output_dir / "stdout.log").open("w") as stdout, (args.output_dir / "stderr.log").open("w") as stderr:
        result = subprocess.run(cmd, stdout=stdout, stderr=stderr, check=False)
    (args.output_dir / "exit_code.txt").write_text(f"{result.returncode}\n")
    raise SystemExit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true", help="Actually launch ASTER/ASTRAL4; default is dry-run.")
    parser.add_argument("--allow-cluster-execution", action="store_true")
    parser.add_argument("--aster-bin", type=Path, required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--treatment")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--input-manifest", type=Path)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    if args.preflight:
        for required in ("input_manifest", "data_dir", "output_root"):
            if getattr(args, required) is None:
                parser.error(f"--preflight requires --{required.replace('_', '-')}")
        run_preflight(args)
        return
    if args.input is None or args.output_dir is None or not args.treatment:
        parser.error("--input, --output-dir, and --treatment are required outside --preflight")
    run_treatment(args)


if __name__ == "__main__":
    main()
