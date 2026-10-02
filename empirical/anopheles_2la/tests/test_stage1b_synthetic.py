#!/usr/bin/env python3
"""Regression tests for the synthetic-only Stage 1B/2 pipeline."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "anopheles_2la"


def load_ingest():
    path = EMPIRICAL_ROOT / "scripts" / "10_ingest_local_trees.py"
    spec = importlib.util.spec_from_file_location("stage1b_ingest_test", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["stage1b_ingest_test"] = module
    spec.loader.exec_module(module)
    return module


def load_script(name: str, filename: str):
    path = EMPIRICAL_ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Stage1BScenarioTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.m = load_ingest()
        cls.root = EMPIRICAL_ROOT / "tests" / "fixtures"
        cls.m.generate_synthetic_fixtures(cls.root)
        cls.windows = cls.m.read_tsv(cls.root / "synthetic_genomic_windows.tsv")
        cls.quartets = cls.m.read_tsv(cls.root / "synthetic_frozen_strict_quartets_stage1a.tsv")

    def scenario(self, name: str):
        return self.m.ingest_local_trees(self.windows, self.quartets, self.root / "synthetic_gene_trees" / name)

    def test_guard_blocks_real_mode_without_stage1a_marker(self) -> None:
        with self.assertRaises(RuntimeError):
            self.m.guard_inputs([], False, True)

    def test_synthetic_manifest_is_labeled_and_fake(self) -> None:
        rows = self.m.read_tsv(self.root / "synthetic_stage1a_manifest.tsv")
        self.assertEqual(len(rows), 72)
        self.assertTrue(all(r["synthetic_fixture"] == "true" for r in rows))
        self.assertTrue(all(r["sample_id"].startswith("SYN_") for r in rows))

    def test_null_not_systematically_positive(self) -> None:
        test = self.m.circular_shift_test(self.scenario("null_msc_like"), 99)
        self.assertLess(float(test["observed_delta"]), 0.25)

    def test_msrc_positive_detected(self) -> None:
        test = self.m.circular_shift_test(self.scenario("msrc_positive"), 99)
        self.assertGreater(float(test["observed_delta"]), 0.35)
        self.assertLess(float(test["p_value"]), 0.05)

    def test_broad_alternative_not_2la_specific(self) -> None:
        summary = self.m.summarize_inside_outside(self.scenario("broad_alternative"))
        eligible = [r for r in summary if "B_NULL" not in r["quartet_id"]]
        self.assertTrue(all(float(r["f_out"]) > 0.45 for r in eligible))
        self.assertTrue(all(abs(float(r["delta_arr"])) < 0.25 for r in eligible))

    def test_design_b_and_geography_controls(self) -> None:
        rows = self.scenario("msrc_positive")
        contrasts = self.m.design_b_contrasts(rows, self.root / "synthetic_design_b_pairs.tsv")
        labels = {r["contrast_id"]: r["classification"] for r in contrasts}
        self.assertEqual(labels["SYN_DESIGN_B_POSITIVE"], "changed_in_predicted_direction")
        self.assertEqual(labels["SYN_DESIGN_B_NULL"], "null_like")
        geography = self.m.geography_summary(self.quartets)
        self.assertTrue(any(r["country_matched"] is True for r in geography))


class Stage1BGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.gate = load_script("stage1b_gate_test", "stage1b_gate.py")
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "processed").mkdir()
        (self.root / "metadata").mkdir()
        (self.root / "processed" / "stage1a_freeze_complete.json").write_text('{"status": "SUCCESS"}\n')
        (self.root / "metadata" / "stage1a_malariagen_api_provenance.json").write_text('{"retrieval_status": "SUCCESS"}\n')

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def write_rows(self, filename: str, rows: list[dict[str, str]]) -> None:
        fields = ["quartet_id", "design_class", "sample_1", "sample_2", "sample_3", "sample_4"]
        with (self.root / "processed" / filename).open("w", newline="") as handle:
            import csv

            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})

    def test_case_a_design_a_only_passes(self) -> None:
        self.write_rows("frozen_strict_quartets_stage1a.tsv", [{"quartet_id": "a", "design_class": "Design_A_four_species_strict_2x2"}])
        result = self.gate.evaluate_stage1b_gate(self.root)
        self.assertEqual(result.status, "READY")
        self.assertTrue(result.design_A_enabled)
        self.assertFalse(result.design_B_enabled)
        self.assertFalse(result.design_C_enabled)

    def test_case_b_design_b_only_passes(self) -> None:
        self.write_rows("frozen_strict_quartets_stage1a.tsv", [])
        self.write_rows("frozen_design_b_contrasts_stage1a.tsv", [{"quartet_id": "b", "design_class": "Design_B_arrangement_replacement"}])
        result = self.gate.evaluate_stage1b_gate(self.root)
        self.assertEqual(result.status, "READY")
        self.assertFalse(result.design_A_enabled)
        self.assertTrue(result.design_B_enabled)

    def test_case_c_design_c_only_passes(self) -> None:
        self.write_rows("frozen_strict_quartets_stage1a.tsv", [])
        self.write_rows("frozen_design_c_geography_controls_stage1a.tsv", [{"quartet_id": "c", "design_class": "Design_C_geography_population_matched"}])
        result = self.gate.evaluate_stage1b_gate(self.root)
        self.assertEqual(result.status, "READY")
        self.assertTrue(result.design_C_enabled)

    def test_case_d_all_empty_clean_stop_before_topology_access(self) -> None:
        self.write_rows("frozen_strict_quartets_stage1a.tsv", [])
        result = self.gate.evaluate_stage1b_gate(self.root)
        self.assertEqual(result.status, self.gate.NO_ELIGIBLE_STATUS)
        ingest = load_ingest()
        original_root = ingest.DATA_ROOT
        ingest.DATA_ROOT = self.root
        try:
            with mock.patch.object(ingest, "ingest_local_trees", side_effect=AssertionError("topology accessed")):
                code = ingest.main([
                    "--real-mode",
                    "--windows",
                    str(self.root / "windows.tsv"),
                    "--quartets",
                    str(self.root / "processed" / "frozen_strict_quartets_stage1a.tsv"),
                    "--tree-dir",
                    str(self.root / "trees"),
                    "--out",
                    str(self.root / "out.tsv"),
                ])
        finally:
            ingest.DATA_ROOT = original_root
        self.assertEqual(code, 0)

    def test_case_e_missing_completion_marker_fails(self) -> None:
        (self.root / "processed" / "stage1a_freeze_complete.json").unlink()
        with self.assertRaises(self.gate.Stage1BGateError):
            self.gate.evaluate_stage1b_gate(self.root)

    def test_case_f_fetch_run_tests_does_not_call_live_api_path(self) -> None:
        fetch = load_script("stage1a_fetch_test", "01_fetch_sample_karyotypes.py")

        class DummyResult:
            def wasSuccessful(self) -> bool:
                return True

        with mock.patch.object(fetch, "run", side_effect=AssertionError("live API path invoked")):
            with mock.patch.object(fetch.unittest.TextTestRunner, "run", return_value=DummyResult()):
                self.assertEqual(fetch.main(["--run-tests"]), 0)


if __name__ == "__main__":
    raise SystemExit(0 if unittest.main(verbosity=2, exit=False).result.wasSuccessful() else 1)
