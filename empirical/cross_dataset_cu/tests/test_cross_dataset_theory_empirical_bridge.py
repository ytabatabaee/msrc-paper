#!/usr/bin/env python3
import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_cross_dataset_theory_empirical_bridge.py"
SPEC = importlib.util.spec_from_file_location("bridge", SCRIPT)
bridge = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(bridge)


class BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.mouse, cls.excluded = bridge.load_threshold_rows()
        cls.mouse["patterns"] = bridge.load_mouse()[1]

    def test_threshold_formula_and_invalid_crossing(self):
        self.assertAlmostEqual(bridge.threshold(0.8, 0.2, 0.3, 0.5), 0.75)
        self.assertIsNone(bridge.threshold(0.8, 0.2, 0.5, 0.3))

    def test_mouse_weights_and_reconstruction(self):
        patterns = self.mouse["patterns"]
        self.assertAlmostEqual(sum(float(patterns[x]["pattern_weight"]) for x in ("SSS", "TTT")), 0.35357142857)
        self.assertAlmostEqual(sum(float(patterns[x]["pattern_weight"]) for x in ("SST", "STS", "STT", "TSS", "TST", "TTS")), 0.64642857143)
        self.assertLess(self.mouse["max_error"], 2e-8)

    def test_empirical_sides_of_threshold(self):
        an_fire = [r for r in self.rows if r["dataset"] != "House mouse chr17 t-complex"]
        self.assertTrue(all(r["observed_below_or_above_threshold"] == "below" for r in an_fire))
        self.assertEqual(self.rows[-1]["observed_below_or_above_threshold"], "above")

    def test_topology_labels(self):
        self.assertEqual(self.rows[-1]["background_topology"], "Q_SPECIES")
        self.assertEqual(self.rows[-1]["pooled_topology"], "Q_T_ALT")
        self.assertEqual(self.rows[-1]["affected_only_topology"], "Q_T_ALT")


if __name__ == "__main__":
    unittest.main()
