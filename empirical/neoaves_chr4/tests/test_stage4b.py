import importlib.util
import json
import math
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts/04b_inference_treatments.py"
SPEC = importlib.util.spec_from_file_location("stage4b", SCRIPT)
stage4b = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(stage4b)


class Stage4BTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bp = stage4b.read_tsv(stage4b.BREAKPOINTS)
        cls.masks = stage4b.build_masks(cls.bp)
        cls.tracks = stage4b.load_tracks()

    def test_t0_t1_reproduce_stage4a_exactly(self):
        prior = stage4b.read_tsv(stage4b.STAGE4A_WB)
        for clade, track in self.tracks.items():
            for treatment, scheme in (("T0", "window-weighted"), ("T1", "block-normalized")):
                got = stage4b.summarize(track, stage4b.weights_for(track, treatment, self.masks["primary"]))
                want = next(r for r in prior if r["clade"] == clade and r["weighting_scheme"] == scheme)
                for i, q in enumerate(stage4b.TOPOLOGIES, 1):
                    self.assertAlmostEqual(got[q], float(want[f"q_{i}"]), places=11)

    def test_t2_and_masks_use_only_frozen_structural_data(self):
        altered = [dict(r, invented_topology="q3") for r in self.bp]
        self.assertEqual(self.masks, stage4b.build_masks(altered))
        manifest = json.loads(stage4b.STAGE4A_MANIFEST.read_text())
        for rel, digest in manifest["inputs"].items():
            self.assertEqual(stage4b.sha256(stage4b.ROOT / rel), digest)

    def test_t3_is_topology_neutral_and_weights_sum(self):
        track = self.tracks["N61"]
        weights = stage4b.weights_for(track, "T3", self.masks["primary"])
        changed = [dict(r, topology="q1" if r["topology"] != "q1" else "q3") for r in track]
        self.assertEqual(weights, stage4b.weights_for(changed, "T3", self.masks["primary"]))
        expected = sum(stage4b.T3_ASSOCIATED_WEIGHT if stage4b.in_mask(float(r["midpoint"]), self.masks["primary"]) else 1 for r in track)
        self.assertAlmostEqual(sum(weights), expected)

    def test_frequencies_and_margin(self):
        for track in self.tracks.values():
            for treatment in stage4b.TREATMENTS:
                result = stage4b.summarize(track, stage4b.weights_for(track, treatment, self.masks["primary"]))
                self.assertAlmostEqual(sum(result[q] for q in stage4b.TOPOLOGIES), 1)
                self.assertAlmostEqual(result["M_species"], result["q1"] - max(result["q2"], result["q3"]))

    def test_deterministic_blocks(self):
        track = [{"topology": q} for q in ("q1", "q1", "q2", "q3", "q3")]
        self.assertEqual([[r["topology"] for r in run] for run in stage4b.make_runs(track)], [["q1", "q1"], ["q2"], ["q3", "q3"]])
        self.assertEqual(stage4b.make_runs(track), stage4b.make_runs(track))

    def test_leave_one_event_logic(self):
        events = stage4b.frozen_events()
        event = events[0]
        removed = {event["left_breakpoint_id"], event["right_breakpoint_id"]} - {""}
        loo = stage4b.build_masks(self.bp, removed)["primary"]
        self.assertNotEqual(loo, self.masks["primary"])
        self.assertEqual(len(events), 16)

    def test_stage3b_remains_underpowered(self):
        stage4b.validate_frozen_inputs()


if __name__ == "__main__":
    unittest.main()
