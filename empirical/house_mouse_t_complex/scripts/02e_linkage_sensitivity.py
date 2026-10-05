#!/usr/bin/env python3
"""Stage 2E spatial thinning and pseudo-replication sensitivity."""

from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    STAGE2_DATA,
    find_astral,
    load_stage1,
    rel,
    run_astral,
    spatial_thin_records,
    write_tsv,
    write_tree_file,
)


OUT = RESULTS / "stage2_linkage_sensitivity.tsv"
BLOCKS = RESULTS / "stage2_source_defined_blocks.tsv"
OUT_DIR = RESULTS / "stage2_linkage_astral"
DATA_DIR = STAGE2_DATA / "linkage_astral"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path("IST-2017-78-v1+1_Data.zip"))
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--max-phases", type=int, default=10, help="Bound local runtime; script remains deterministic and resumable.")
    return p.parse_args()


def run(args: argparse.Namespace) -> None:
    meta, trees, _mapping = load_stage1()
    astral = find_astral()
    map_file = PROCESSED / "house_mouse_t_complex_astral_subspecies.map"
    rows = []
    for spacing in [5_000, 10_000, 20_000, 50_000, 100_000, 250_000, 500_000]:
        phases = list(range(0, spacing, 5_000))
        if len(phases) > args.max_phases:
            phases = phases[: args.max_phases]
        for phase in phases:
            idx = spatial_thin_records(meta, spacing, phase)
            name = f"spacing_{spacing // 1000}kb_phase_{phase // 1000}kb"
            tree_file = DATA_DIR / f"{name}.tre"
            if not tree_file.exists():
                write_tree_file(tree_file, [trees[i] for i in idx])
            result = run_astral(astral, tree_file, map_file, OUT_DIR / f"{name}.nwk", threads=args.threads)
            rows.append(
                {
                    "spacing_kb": spacing // 1000,
                    "phase_kb": phase // 1000,
                    "n_trees": len(idx),
                    "topology": result.get("topology", "NA"),
                    "matches_q_species": result.get("topology") == "Q_SPECIES",
                    "matches_q_t_alt": result.get("topology") == "Q_T_ALT",
                    "CU": result.get("CU_length", "NA"),
                    "localPP": result.get("localPP", "NA"),
                    "q1": result.get("q1", "NA"),
                    "q2": result.get("q2", "NA"),
                    "q3": result.get("q3", "NA"),
                    "tree_file": rel(tree_file),
                }
            )
    write_tsv(OUT, rows, ["spacing_kb", "phase_kb", "n_trees", "topology", "matches_q_species", "matches_q_t_alt", "CU", "localPP", "q1", "q2", "q3", "tree_file"])
    write_tsv(BLOCKS, inspect_blocks(args.archive), ["source_defined_blocks_available", "source_path", "notes"])
    print(f"Wrote {OUT}")


def inspect_blocks(archive: Path) -> list[dict[str, object]]:
    if not archive.exists():
        return [{"source_defined_blocks_available": False, "source_path": "NA", "notes": "Archive not found for block inspection."}]
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
    candidates = [n for n in names if "Non-Recombined" in n or "NonRecomb" in n]
    # The archive contains concatenated non-recombined-region tree files, but no interval table.
    return [
        {
            "source_defined_blocks_available": False,
            "source_path": ",".join(candidates),
            "notes": "Source archive names concatenated non-recombined-region tree files but no machine-readable genomic interval definitions were recovered; no post-hoc blocks inferred.",
        }
    ]


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))

