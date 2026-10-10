#!/usr/bin/env python3
"""Run frozen replicate-experiment configs through msrc-sim-replicates."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_one(root: Path, config_path: Path, executable: str) -> None:
    config = yaml.safe_load(config_path.read_text())
    raw = root / config["output"]["directory"]
    metadata_path = raw / "run_metadata.json"
    complete = all((raw / name).is_file() for name in ("replicate_summary.csv", "prevalence_summary.json", "config.resolved.yaml"))
    if complete and metadata_path.is_file():
        return
    raw.mkdir(parents=True, exist_ok=True)
    command = [executable, "--config", str(config_path)]
    start = time.perf_counter()
    proc = subprocess.run(command, cwd=root, text=True, capture_output=True)
    metadata = {"command": command, "config": str(config_path.relative_to(root)), "config_sha256": sha256(config_path), "seed": int(config["seed"]), "runtime_seconds": time.perf_counter() - start, "return_code": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr, "executable": executable}
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    if proc.returncode != 0:
        raise RuntimeError(f"msrc-sim-replicates failed for {config_path}:\n{proc.stderr}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[5])
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--executable", default=None)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    root = args.repo_root.resolve()
    config_root = root / "simulations/quartet/mechanistic/datasets/configs" / ("smoke" if args.smoke else "pilot")
    configs = sorted(config_root.glob("*.yaml"))
    if not configs:
        raise SystemExit("No mechanistic configs found; run 01_generate_mechanistic_pilot.py")
    executable = args.executable or shutil.which("msrc-sim-replicates")
    if not executable:
        raise SystemExit("Could not locate msrc-sim-replicates")
    if args.workers <= 1:
        for config_path in configs:
            run_one(root, config_path, executable)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, root, config_path, executable) for config_path in configs]
            for future in as_completed(futures):
                future.result()


if __name__ == "__main__":
    main()
