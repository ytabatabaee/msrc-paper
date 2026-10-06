#!/usr/bin/env python3
"""Audit deterministic label reconciliation for the coverage-only archive."""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

from Bio import Phylo

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import ML_ARCHIVES, make_tree_records, nested_zip_bytes  # noqa: E402
from house_mouse_stage2_utils import RESULTS, load_stage1, write_tsv  # noqa: E402
from importlib.util import module_from_spec, spec_from_file_location

spec = spec_from_file_location("filtering", SCRIPT_DIR / "02d_filtering_robustness.py")
filtering = module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(filtering)


OUT = RESULTS / "stage2_third_filter_tip_label_audit.tsv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    return p.parse_args()


def run(args: argparse.Namespace) -> None:
    _meta, _trees, mapping = load_stage1()
    primary = {row["tree_tip"] for row in mapping}
    with zipfile.ZipFile(args.archive) as outer:
        data, _ = nested_zip_bytes(outer, ML_ARCHIVES["Coverage_RAW_SNPs"])
    if data is None:
        raise RuntimeError("Coverage_RAW_SNPs ML archive is missing")
    records, _ = make_tree_records(data)
    first = next(record for record in records if record.tree_valid)
    tree = Phylo.read(io.StringIO(first.newick), "newick")
    rows = []
    for terminal in sorted(tree.get_terminals(), key=lambda x: x.name):
        raw = terminal.name
        normalized = filtering.normalize_filter_tip(raw)
        if normalized in primary:
            match_type = "exact_deterministic_normalization"
            matched = normalized
            notes = "Only OverallCovFiltered, DelRemoved, and the tHaplSubset .fa suffix were normalized."
        else:
            match_type = "unmatched"
            matched = "NA"
            notes = "No exact primary mapping label after deterministic normalization; no fuzzy matching attempted."
        rows.append({"raw_tip": raw, "normalized_tip": normalized, "matched_primary_tip": matched, "match_type": match_type, "notes": notes})
    write_tsv(OUT, rows, ["raw_tip", "normalized_tip", "matched_primary_tip", "match_type", "notes"])
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    run(parse_args())
