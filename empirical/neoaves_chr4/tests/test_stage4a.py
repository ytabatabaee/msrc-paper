import csv
import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts/04a_structural_topology_support.py"
SPEC = importlib.util.spec_from_file_location("stage4a", SCRIPT)
stage4a = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(stage4a)


class Stage4ATests(unittest.TestCase):
    def test_structural_masks_use_only_structural_columns(self):
        rows = [{"reference_position": "1000000", "in_stringent": "true", "in_primary": "true", "in_inclusive": "true", "dominant_topology": "q1"}]
        a = stage4a.build_masks(rows)
        rows[0]["dominant_topology"] = "q3"
        self.assertEqual(a, stage4a.build_masks(rows))

    def test_block_collapsing(self):
        track = [{"dominant_topology": q, "start": i * 10, "end": i * 10 + 9} for i, q in enumerate(["q1", "q1", "q2", "q2", "q1"])]
        runs = stage4a.collapse_topology_runs(track)
        self.assertEqual([r["topology"] for r in runs], ["q1", "q2", "q1"])
        self.assertEqual([r["n_windows"] for r in runs], [2, 2, 1])

    def test_topology_count_conservation(self):
        tracks = stage4a.load_tracks()
        for track in tracks.values():
            counts = stage4a.topology_summary(track, "q1")["counts"]
            self.assertEqual(sum(counts.values()), len(track))

    def test_deterministic_circular_shifts(self):
        track = [{"midpoint": i, "dominant_topology": q} for i, q in enumerate(["q1", "q2", "q2", "q3", "q1"])]
        mask = [(0, 1)]
        self.assertEqual(stage4a.circular_enrichment(track, mask, "q2", 100, 7), stage4a.circular_enrichment(track, mask, "q2", 100, 7))

    def test_event_unit_assignment_and_count(self):
        with stage4a.EVENTS.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        primary = [r for r in rows if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == 250000]
        self.assertEqual(len(primary), 16)
        self.assertEqual(len({r["event_id"] for r in primary}), 16)

    def test_no_stage3b_eligibility_relaxation(self):
        stage4a.validate_stage3b_unchanged()


if __name__ == "__main__":
    unittest.main()
