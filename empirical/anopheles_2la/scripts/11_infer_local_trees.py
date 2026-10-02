#!/usr/bin/env python3
"""Prepare deterministic local-tree inference commands without topology feedback."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE_ROOT = REPO_ROOT / "empirical" / "anopheles_2la" / "tests" / "fixtures"
DATA_ROOT = REPO_ROOT / "data" / "anopheles_2la"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage1b_gate import NoEligibleProspectiveDesigns, require_stage1b_designs


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def guard(synthetic_validation: bool, real_mode: bool, paths: list[Path]) -> None:
    if synthetic_validation:
        root = FIXTURE_ROOT.resolve()
        for path in paths:
            resolved = path.resolve()
            if root not in [resolved, *resolved.parents]:
                raise RuntimeError(f"synthetic validation input must live under {root}: {path}")
    if real_mode:
        require_stage1b_designs(DATA_ROOT)


def iqtree_version(executable: str) -> str:
    if not shutil.which(executable):
        return "not_found"
    try:
        out = subprocess.run([executable, "--version"], check=False, text=True, capture_output=True, timeout=20)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception as exc:
        return f"version_unavailable: {exc}"


def build_commands(fasta_dir: Path, windows: list[dict[str, str]], outdir: Path, executable: str, seed: int, model_policy: str, threads: str) -> list[dict[str, object]]:
    rows = []
    for window in windows:
        fasta = fasta_dir / f"{window['window_id']}.fasta"
        prefix = outdir / window["window_id"]
        cmd = [executable, "-s", str(fasta), "-seed", str(seed), "-nt", threads, "-pre", str(prefix)]
        if model_policy == "iqtree_model_selection":
            cmd.extend(["-m", "MFP"])
        else:
            cmd.extend(["-m", model_policy])
        rows.append(
            {
                "window_id": window["window_id"],
                "input_alignment": str(fasta),
                "command": " ".join(cmd),
                "software": "IQ-TREE",
                "software_version": iqtree_version(executable),
                "model_policy": model_policy,
                "seed": seed,
                "thread_policy": threads,
                "bootstrap_support_policy": "none_primary_analysis",
                "exit_status": "not_run",
                "resulting_tree": str(prefix) + ".treefile",
                "expected_outputs": ";".join(
                    [
                        str(prefix) + ".treefile",
                        str(prefix) + ".iqtree",
                        str(prefix) + ".log",
                        str(prefix) + ".ckp.gz",
                        str(prefix) + ".model.gz",
                    ]
                ),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fields = list(rows[0].keys()) if rows else ["window_id", "input_alignment", "command", "software", "software_version", "model_policy", "seed", "thread_policy", "bootstrap_support_policy", "exit_status", "resulting_tree", "expected_outputs"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", type=Path)
    parser.add_argument("--windows", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--estimator", default="iqtree2")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--model-policy", default="iqtree_model_selection")
    parser.add_argument("--threads", default="AUTO")
    parser.add_argument("--synthetic-validation", action="store_true")
    parser.add_argument("--real-mode", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        return run_tests()
    if not (args.fasta_dir and args.windows and args.outdir and args.manifest):
        parser.error("--fasta-dir, --windows, --outdir, and --manifest are required unless running tests")
    try:
        guard(args.synthetic_validation, args.real_mode, [args.fasta_dir, args.windows])
    except NoEligibleProspectiveDesigns as exc:
        print(exc.result.status)
        print(exc.result.message)
        return 0
    rows = build_commands(args.fasta_dir, read_tsv(args.windows), args.outdir, args.estimator, args.seed, args.model_policy, args.threads)
    if args.execute:
        for row in rows:
            cmd = row["command"].split()
            result = subprocess.run(cmd, check=False)
            row["exit_status"] = result.returncode
    write_tsv(args.manifest, rows)
    args.manifest.with_suffix(".json").write_text(
        json.dumps(
            {
                "tree_estimation_separate_from_quartet_interpretation": True,
                "prediction_table_used": False,
                "structural_prediction_inputs_forbidden": [
                    "predicted arrangement split",
                    "Design-A/B/C topology prediction",
                    "q1/q2/q3 expectation",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    return 0


class TreeInferenceTests(unittest.TestCase):
    def test_command_uses_only_alignment_window_and_fixed_policy(self) -> None:
        rows = build_commands(
            Path("alignments"),
            [{"window_id": "win_0001"}],
            Path("trees"),
            "iqtree2",
            1729,
            "iqtree_model_selection",
            "AUTO",
        )
        command = rows[0]["command"]
        self.assertEqual(command, "iqtree2 -s alignments/win_0001.fasta -seed 1729 -nt AUTO -pre trees/win_0001 -m MFP")
        forbidden = ["predicted", "arrangement", "Design_A", "Design_B", "Design_C", "q1", "q2", "q3"]
        self.assertFalse(any(token in command for token in forbidden))
        self.assertNotIn("prediction_table", rows[0])


def run_tests() -> int:
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TreeInferenceTests))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
