import importlib.util
import json
import math
import unittest
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts/04c_independent_validation.py"
SPEC = importlib.util.spec_from_file_location("stage4c", SCRIPT)
stage4c = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(stage4c)


class Stage4CTests(unittest.TestCase):
    def test_reference_selection_hierarchy_and_chr4_exclusion(self):
        candidates = stage4c.candidate_rows()
        selected = next(r for r in candidates if r["reference_id"] == stage4c.REFERENCE_ID)
        self.assertTrue(selected["available"])
        self.assertIn("genome-wide quartet", selected["hierarchy_rank"])
        self.assertTrue(stage4c.is_chr4("chr4"))
        self.assertTrue(stage4c.is_chr4("chr4_AADN_random"))
        self.assertFalse(stage4c.is_chr4("chr14"))

    def test_q_mapping_formula(self):
        q = stage4c.compute_qqs({"1": 10, "2": 14, "3": 14, "4": 8})
        self.assertIsNotNone(q)
        self.assertAlmostEqual(sum(q.values()), 1)
        self.assertEqual(max(q, key=q.get), "q1")

    def test_reference_margin(self):
        self.assertAlmostEqual(stage4c.reference_margin({"q1": .5, "q2": .3, "q3": .2}, "q1"), .2)
        self.assertAlmostEqual(stage4c.reference_margin({"q1": .2, "q2": .5, "q3": .3}, "q2"), .2)

    def test_stage4b_reproduced_and_treatments_unchanged(self):
        manifest = json.loads(stage4c.STAGE4B_MANIFEST.read_text())
        stage4c.validate_frozen()
        self.assertEqual(manifest["treatments"]["T3"], "frozen-mask windows weight 0.5; background weight 1")
        rows = stage4c.read_tsv(stage4c.STAGE4B)
        n61 = next(r for r in rows if r["structural_mask"] == "primary" and r["clade"] == "N61" and r["treatment"] == "T0")
        self.assertAlmostEqual(float(n61["q2"]), .453848662753)

    def test_controls_deterministic(self):
        synthetic = [
            {"start": 1, "end": 2, "gene": "1", "topology": "q1"},
            {"start": 3, "end": 4, "gene": "2", "topology": "q1"},
            {"start": 5, "end": 6, "gene": "3", "topology": "q2"},
        ]
        self.assertEqual(stage4c.make_runs(synthetic), stage4c.make_runs(list(reversed(synthetic))))
        self.assertEqual([len(x) for x in stage4c.make_runs(synthetic)], [2, 1])

    def test_resampling_unit_is_chromosome(self):
        source = SCRIPT.read_text()
        self.assertIn("chromosome_bootstrap", source)
        self.assertNotIn("bootstrap individual windows", source)
        manifest = json.loads((stage4c.RESULTS / "stage4c_manifest.json").read_text()) if (stage4c.RESULTS / "stage4c_manifest.json").exists() else None
        if manifest:
            self.assertIn("chromosome", manifest["uncertainty_unit"])

    def test_reference_tree_extraction_not_applicable(self):
        candidates = stage4c.candidate_rows()
        self.assertFalse(next(r for r in candidates if r["reference_id"] == "resolved_genetrees_absent")["available"])


if __name__ == "__main__":
    unittest.main()
