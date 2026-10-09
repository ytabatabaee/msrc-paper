import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_stage2c_effective_msrc_validation.py"
SPEC = importlib.util.spec_from_file_location("stage2c", SCRIPT)
stage2c = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(stage2c)


class Stage2CEffectiveMSRCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aggregate = stage2c.load_aggregate()
        cls.rows = stage2c.prediction_rows(cls.aggregate)

    def test_symmetric_formula(self):
        self.assertAlmostEqual(stage2c.delta_from_q(0.6856), 0.5284, places=4)
        self.assertAlmostEqual(sum(stage2c.predicted_q(0.5284)), 1.0)

    def test_frozen_aggregate_values(self):
        self.assertAlmostEqual(self.aggregate["STT"]["q_arrangement"], 0.6856, places=4)
        self.assertAlmostEqual(self.aggregate["TST"]["q_arrangement"], 0.6937, places=4)
        self.assertAlmostEqual(self.aggregate["TTS"]["q_arrangement"], 0.7068, places=4)

    def test_leave_one_out_and_stt_only(self):
        loo = [r for r in self.rows if r["analysis"] == "leave_one_configuration_out"]
        stt = [r for r in self.rows if r["analysis"] == "STT_only_prediction"]
        self.assertEqual(len(loo), 3)
        self.assertEqual(len(stt), 2)
        self.assertTrue(all(float(r["max_abs_error"]) < 0.03 for r in loo))

    def test_block_assignment_is_complete(self):
        blocks = stage2c.load_block_rows(stage2c.PRIMARY_BLOCK_SIZE)
        self.assertGreater(len(blocks), 5)
        self.assertTrue(all(set(v) == set(stage2c.PATTERNS) for v in blocks.values()))


if __name__ == "__main__":
    unittest.main()
