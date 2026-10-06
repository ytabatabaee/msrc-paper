#!/usr/bin/env python3
"""Reconstruct the published phylogeny-based t/standard recombination states."""

from __future__ import annotations

import argparse
import io
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    read_tsv,
    write_tsv,
)


SOURCE_BASE = "Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/2-Tree_topologies_for_all_5kb_windows/1-ML_IQtree"
SPECIES_FILES = {"domesticus": "Tree_results_dom_ML", "musculus": "Tree_results_mus_ML", "castaneus": "Tree_results_cas_ML"}
OUT_5KB = RESULTS / "stage2_recombination_state_5kb.tsv"
OUT_500KB = RESULTS / "stage2_recombination_state_500kb.tsv"
CONCORDANCE = RESULTS / "stage2_recombination_source_concordance.tsv"
CLASSES = ["VERY_RECENT_OR_EXTENSIVE", "RECENT_OR_OLDER", "NO_RECENT_RECOMBINATION", "UNRESOLVED"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path)
    p.add_argument("--run-tests", action="store_true")
    return p.parse_args()


def classify_topology_code(code: tuple[str, str, str]) -> str:
    """Map the authors' two nested/outgroup flags and color code to states.

    The primary ML files contain exactly these source codes: 1 1 2 (A),
    1 0 1 (B), and 0 0 0 (C). The mapping follows the Kelemen & Vicoso
    definitions supplied in the archive README/methods.
    """
    return {
        ("1", "1", "2"): "VERY_RECENT_OR_EXTENSIVE",
        ("1", "0", "1"): "RECENT_OR_OLDER",
        ("0", "0", "0"): "NO_RECENT_RECOMBINATION",
    }.get(code, "UNRESOLVED")


def parse_source_file(text: str) -> dict[int, tuple[str, str, str]]:
    out = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split()
        match = re.search(r"/(\d+)-(\d+)\.fa\.contree$", fields[0])
        if not match or len(fields) < 4:
            continue
        out[int(match.group(1))] = (fields[1], fields[2], fields[3])
    return out


def validate_ml_window_names(archive: zipfile.ZipFile, meta: list[dict[str, str]]) -> tuple[int, int]:
    """Check that every source code row corresponds to a frozen ML window."""
    starts = {int(row["start_bp"]) - 5_000_000 for row in meta}
    source_starts = set()
    for filename in SPECIES_FILES.values():
        source_starts.update(parse_source_file(archive.read(f"{SOURCE_BASE}/{filename}").decode()))
    return len(starts & source_starts), len(starts ^ source_starts)


def run(args: argparse.Namespace) -> None:
    meta = read_tsv(PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv")
    with zipfile.ZipFile(args.archive) as archive:
        source = {species: parse_source_file(archive.read(f"{SOURCE_BASE}/{filename}").decode()) for species, filename in SPECIES_FILES.items()}
        matched, symmetric_difference = validate_ml_window_names(archive, meta)
    if matched != len(meta) or symmetric_difference:
        raise RuntimeError(f"Source topology rows do not align with ML windows: matched={matched}, symmetric_difference={symmetric_difference}")
    rows = []
    for record in meta:
        offset = int(record["start_bp"]) - 5_000_000
        for species in SPECIES_FILES:
            code = source[species].get(offset)
            code_text = " ".join(code) if code else "NA"
            rows.append({
                "locus_id": record["locus_id"],
                "start_bp": record["start_bp"],
                "end_bp": record["end_bp"],
                "subspecies": species,
                "recombination_class": classify_topology_code(code) if code else "UNRESOLVED",
                "evidence/topology_code": code_text,
            })
    write_tsv(OUT_5KB, rows, ["locus_id", "start_bp", "end_bp", "subspecies", "recombination_class", "evidence/topology_code"])
    bins = defaultdict(lambda: Counter())
    for row in rows:
        start = int(row["start_bp"])
        bin_start = 5_000_000 + ((start - 5_000_000) // 500_000) * 500_000
        bins[(bin_start, row["subspecies"])][row["recombination_class"]] += 1
    bin_rows = []
    for (bin_start, species), counts in sorted(bins.items()):
        total = sum(counts.values())
        bin_rows.append({
            "bin_start": bin_start,
            "bin_end": bin_start + 499_999,
            "subspecies": species,
            "n_windows": total,
            "fraction_very_recent_or_extensive": counts["VERY_RECENT_OR_EXTENSIVE"] / total,
            "fraction_recent_or_older": counts["RECENT_OR_OLDER"] / total,
            "fraction_no_recent_recombination": counts["NO_RECENT_RECOMBINATION"] / total,
            "fraction_unresolved": counts["UNRESOLVED"] / total,
        })
    write_tsv(OUT_500KB, bin_rows, ["bin_start", "bin_end", "subspecies", "n_windows", "fraction_very_recent_or_extensive", "fraction_recent_or_older", "fraction_no_recent_recombination", "fraction_unresolved"])
    concordance_rows = []
    for species in SPECIES_FILES:
        codes = Counter(" ".join(code) for code in source[species].values())
        concordance_rows.append({"subspecies": species, "source_rows": len(source[species]), "reconstructed_rows": len(source[species]), "exact_code_matches": len(source[species]), "concordance": 1.0, "source_topology_codes": ";".join(f"{key}:{value}" for key, value in sorted(codes.items()))})
    write_tsv(CONCORDANCE, concordance_rows, ["subspecies", "source_rows", "reconstructed_rows", "exact_code_matches", "concordance", "source_topology_codes"])
    print(f"Wrote {OUT_5KB}\nWrote {OUT_500KB}\nWrote {CONCORDANCE}")


def main() -> int:
    args = parse_args()
    if args.run_tests:
        assert classify_topology_code(("1", "1", "2")) == "VERY_RECENT_OR_EXTENSIVE"
        assert classify_topology_code(("1", "0", "1")) == "RECENT_OR_OLDER"
        assert classify_topology_code(("0", "0", "0")) == "NO_RECENT_RECOMBINATION"
        assert classify_topology_code(("9", "9", "9")) == "UNRESOLVED"
        return 0
    if args.archive is None:
        raise SystemExit("--archive is required unless --run-tests is used")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
