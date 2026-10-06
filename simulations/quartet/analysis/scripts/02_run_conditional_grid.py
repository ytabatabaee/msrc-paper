#!/usr/bin/env python3
"""Run frozen conditional configs through the public msrc-sim CLI."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complete(raw_dir: Path) -> bool:
    return all((raw_dir / name).is_file() for name in ("summary.json", "quartet_probabilities.csv", "config.resolved.yaml"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[4])
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--msrc-sim", default=None, help="Executable; defaults to msrc-sim on PATH")
    args = ap.parse_args()
    root = args.repo_root.resolve()
    config_root = root / "simulations/quartet/datasets/configs/conditional_grid" / ("smoke" if args.smoke else "")
    configs = sorted(config_root.glob("*.yaml"))
    if not configs:
        raise SystemExit("No frozen configs found; run 01_generate_conditional_grid.py first")
    executable = args.msrc_sim or shutil.which("msrc-sim")
    if not executable:
        raise SystemExit("Could not locate msrc-sim executable")
    for config_path in configs:
        config = yaml.safe_load(config_path.read_text())
        raw_dir = root / config["output"]["directory"]
        metadata_path = raw_dir / "run_metadata.json"
        if complete(raw_dir) and metadata_path.is_file():
            continue
        raw_dir.mkdir(parents=True, exist_ok=True)
        command = [executable, "--config", str(config_path)]
        started = time.perf_counter()
        proc = subprocess.run(command, cwd=root, text=True, capture_output=True)
        elapsed = time.perf_counter() - started
        metadata = {
            "command": command,
            "config": str(config_path.relative_to(root)),
            "config_sha256": sha256(config_path),
            "seed": int(config["seed"]),
            "runtime_seconds": elapsed,
            "return_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "msrc_sim_executable": executable,
        }
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        if proc.returncode != 0:
            raise RuntimeError(f"msrc-sim failed for {config_path}:\n{proc.stderr}")


if __name__ == "__main__":
    main()
