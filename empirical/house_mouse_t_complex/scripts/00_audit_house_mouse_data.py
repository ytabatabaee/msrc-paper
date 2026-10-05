#!/usr/bin/env python3
"""Stage 0 audit for the Kelemen & Vicoso house-mouse t-complex archive."""

from __future__ import annotations

import argparse
import sys
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import (  # noqa: E402
    PRIMARY_ARCHIVE,
    RESULTS,
    archive_audit,
    parse_coordinates_from_name,
    parse_newick,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Path to IST-2017-78-v1+1_Data.zip")
    parser.add_argument("--run-tests", action="store_true", help="Run lightweight synthetic tests.")
    return parser.parse_args()


class Stage0SyntheticTests(unittest.TestCase):
    def test_newick_inventory_and_malformed_handling(self) -> None:
        parsed = parse_newick("(dom:0.1,(mus:0.2,cas:0.3)99:0.4,spretus:0.5);")
        self.assertEqual(tuple(parsed["tips"]), ("dom", "mus", "cas", "spretus"))
        self.assertTrue(parsed["branch_lengths_present"])
        self.assertTrue(parsed["internal_labels_present"])
        with self.assertRaises(Exception):
            parse_newick("(a:0.1,b:0.2)")

    def test_coordinate_parsing(self) -> None:
        self.assertEqual(parse_coordinates_from_name("0-5000.fa.contree")[:2], (5000000, 5004999))
        self.assertEqual(parse_coordinates_from_name("chr17_5000000_5004999.treefile")[:2], (10000000, 10004999))
        self.assertEqual(parse_coordinates_from_name("window_start_5000000.treefile")[:2], (10000000, 10004999))
        self.assertEqual(parse_coordinates_from_name("treefile_001.treefile")[:2], (None, None))

    def test_nested_zip_audit_deterministic_outputs(self) -> None:
        with TemporaryDirectory() as tmp:
            archive = Path(tmp) / "IST-2017-78-v1+1_Data.zip"
            nested = Path(tmp) / "ML_trees.zip"
            with zipfile.ZipFile(nested, "w") as z:
                z.writestr("chr17_5000000_5004999.treefile", "(a:1,b:1,c:1,d:1);")
                z.writestr("chr17_5005000_5009999.treefile", "(a:1,b:1,c:1,d:1);")
                z.writestr("chr17_5010000_5014999.treefile", "(a:1,b:1,c:1,d:1)")
            with zipfile.ZipFile(archive, "w") as z:
                z.write(nested, PRIMARY_ARCHIVE)
                z.writestr(
                    "Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/"
                    "2-Tree_topologies_for_all_5kb_windows/1-ML_IQtree/Tree_results_test.txt",
                    "topology\tcount\n",
                )
                z.writestr("README.txt", "synthetic fixture")
            first = archive_audit(archive, write_outputs=False)
            second = archive_audit(archive, write_outputs=False)
            self.assertEqual(first["tree_inventory"]["n_final_tree_files"], 3)
            self.assertEqual(first["tree_inventory"]["n_valid_newicks"], 2)
            self.assertEqual(first["tree_inventory"]["n_malformed_newicks"], 1)
            self.assertEqual(first["stage0_gate"]["overall_status"], "FAIL")
            self.assertEqual(first["coordinate_summary"], second["coordinate_summary"])
            self.assertEqual(first["stage0_gate"]["overall_status"], second["stage0_gate"]["overall_status"])


def run_tests() -> int:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage0SyntheticTests)
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


def main() -> int:
    args = parse_args()
    if args.run_tests:
        return run_tests()
    if args.archive is None:
        raise SystemExit("--archive is required unless --run-tests is used")
    if not args.archive.exists():
        raise SystemExit(f"Archive not found: {args.archive}")
    manifest = archive_audit(args.archive)
    print(f"Wrote {RESULTS / 'stage0_source_audit.md'}")
    print(f"Stage 0 status: {manifest['stage0_gate']['overall_status']}")
    print(f"Final ML tree files: {manifest['tree_inventory'].get('n_final_tree_files', 0)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
