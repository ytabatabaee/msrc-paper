#!/usr/bin/env python3
"""Initial descriptive summaries for frozen house-mouse local trees."""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import (  # noqa: E402
    FIGURES,
    METADATA,
    PROCESSED,
    RESULTS,
    load_mapping,
    parse_newick,
    read_tsv,
    write_tsv,
)


TREE_FILE = PROCESSED / "house_mouse_t_complex_ml_5kb.tre"
TREE_META = PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv"
SUMMARY_TSV = RESULTS / "stage1_local_tree_summary.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()


def split_for_four(tips: list[str], groups: dict[str, str]) -> str:
    target = ["Mus musculus domesticus", "Mus musculus musculus", "Mus musculus castaneus", "Mus spretus"]
    present = {g: [t for t in tips if groups.get(t) == g] for g in target}
    if any(not present[g] for g in target):
        return "NA"
    return "quartet_available_mapping_not_topology_scored"


def write_inventory_figures(rows: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc_mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/msrc_xdgcache")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    FIGURES.mkdir(parents=True, exist_ok=True)
    starts = np.array([float(r["start_bp"]) for r in rows])
    bins = np.arange(starts.min(), starts.max() + 250_000, 250_000)
    counts, edges = np.histogram(starts, bins=bins)
    if counts.max() == counts.min():
        raise RuntimeError("Inventory figure would be uninformative: 250-kb bin counts are constant.")
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.bar(edges[:-1] / 1e6, counts, width=np.diff(edges) / 1e6, align="edge", color="#4c78a8", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("chr17 position (Mb)")
    ax.set_ylabel("5-kb trees per 250-kb bin")
    ax.set_title("House mouse t-complex window density")
    fig.tight_layout()
    fig.savefig(FIGURES / "house_mouse_t_complex_tree_inventory.png", dpi=220)
    fig.savefig(FIGURES / "house_mouse_t_complex_tree_inventory.pdf")


def main() -> int:
    parse_args()
    if not TREE_FILE.exists() or not TREE_META.exists():
        raise SystemExit("Stage 1 gene-tree file and metadata table are required before 01b.")
    lines = [line.strip() for line in TREE_FILE.read_text().splitlines() if line.strip()]
    meta = read_tsv(TREE_META)
    mapping = {row["tree_tip"]: row.get("subspecies", "NA") for row in load_mapping(METADATA / "tip_mapping.tsv")}
    parsed = [parse_newick(line) for line in lines]
    all_tips = sorted({tip for item in parsed for tip in item["tips"]})
    topology_strings = Counter(line for line in lines)
    rows = []
    for row, item, newick in zip(meta, parsed, lines, strict=True):
        tips = set(item["tips"])
        rows.append(
            {
                "locus_id": row["locus_id"],
                "start_bp": row["start_bp"],
                "end_bp": row["end_bp"],
                "midpoint_bp": row["midpoint_bp"],
                "n_tips": len(tips),
                "n_missing_tips": len(set(all_tips) - tips),
                "fraction_tips_present": f"{len(tips) / len(all_tips):.8f}" if all_tips else "NA",
                "topology_exact_newick_count": topology_strings[newick],
                "domesticus_musculus_castaneus_spretus_quartet_status": split_for_four(list(tips), mapping),
                "q1": "NA",
                "q2": "NA",
                "q3": "NA",
            }
        )
    write_tsv(
        SUMMARY_TSV,
        rows,
        [
            "locus_id",
            "start_bp",
            "end_bp",
            "midpoint_bp",
            "n_tips",
            "n_missing_tips",
            "fraction_tips_present",
            "topology_exact_newick_count",
            "domesticus_musculus_castaneus_spretus_quartet_status",
            "q1",
            "q2",
            "q3",
        ],
    )
    write_inventory_figures(rows)
    print(f"Wrote {SUMMARY_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
