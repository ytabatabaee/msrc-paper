#!/usr/bin/env python3
"""Stage 2D robustness across the three published ML-tree filtering datasets."""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import ML_ARCHIVES, make_tree_records, nested_zip_bytes  # noqa: E402
from house_mouse_stage2_utils import (  # noqa: E402
    PROCESSED,
    RESULTS,
    STAGE2_DATA,
    aggregate_quartet_rows,
    find_astral,
    load_stage1,
    prune_newick,
    quartet_counts_for_tree,
    run_astral,
    score_to_row,
    tips_by_species,
    write_tree_file,
    write_tsv,
)


OUT = RESULTS / "stage2_filtering_robustness.tsv"
OUT_DIR = RESULTS / "stage2_filtering_astral"
DATA_DIR = STAGE2_DATA / "filtering_astral"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    return p.parse_args()


def run(args: argparse.Namespace) -> None:
    _meta, _trees, mapping = load_stage1()
    groups_standard = tips_by_species(mapping, "STANDARD_ONLY")
    groups_all = tips_by_species(mapping, "ALL_TIPS")
    astral = find_astral()
    rows = []
    with zipfile.ZipFile(args.archive) as outer:
        for dataset, expected in ML_ARCHIVES.items():
            data, _info = nested_zip_bytes(outer, expected)
            if data is None:
                rows.append({"dataset": dataset, "status": "missing"})
                continue
            records, coord = make_tree_records(data)
            valid = sorted([r for r in records if r.tree_valid], key=lambda r: (r.start_bp or 0, r.source_filename))
            normalized = [normalize_filter_newick(r.newick) for r in valid]
            for treatment, groups, map_file, tips in [
                ("STANDARD_ONLY", groups_standard, PROCESSED / "house_mouse_t_complex_astral_standard_only_subspecies.map", {t for g in groups_standard.values() for t in g}),
                ("ALL_TIPS", groups_all, PROCESSED / "house_mouse_t_complex_astral_subspecies.map", {t for g in groups_all.values() for t in g}),
            ]:
                qrows = []
                compatible = True
                try:
                    for nw in normalized:
                        qr = score_to_row(quartet_counts_for_tree(nw, groups))
                        qr["treatment"] = treatment
                        qrows.append({k: str(v) for k, v in qr.items()})
                except Exception:
                    compatible = False
                    qrows = []
                if not compatible:
                    rows.append(
                        {
                            "dataset": dataset,
                            "treatment": treatment,
                            "status": "tip_label_incompatible",
                            "n_usable_trees": len(valid),
                            "coordinates_recoverable": coord.get("coordinates_recoverable_for_all_final_trees"),
                            "tip_labels_compatible": False,
                        }
                    )
                    continue
                agg = aggregate_quartet_rows(qrows)[0]
                tree_file = DATA_DIR / f"{dataset}_{treatment}.tre"
                write_tree_file(tree_file, [prune_newick(nw, tips) if treatment == "STANDARD_ONLY" else nw for nw in normalized])
                result = run_astral(astral, tree_file, map_file, OUT_DIR / f"{dataset}_{treatment}.nwk", threads=args.threads)
                rows.append(
                    {
                        "dataset": dataset,
                        "treatment": treatment,
                        "status": "ok",
                        "n_usable_trees": len(valid),
                        "coordinates_recoverable": coord.get("coordinates_recoverable_for_all_final_trees"),
                        "tip_labels_compatible": True,
                        "mean_q_species": agg["mean_q_species"],
                        "mean_q_t_alt": agg["mean_q_t_alt"],
                        "mean_q_other": agg["mean_q_other"],
                        "mean_q_unresolved": agg["mean_q_unresolved"],
                        "astral_topology": result.get("topology", "NA"),
                        "CU": result.get("CU_length", "NA"),
                        "localPP": result.get("localPP", "NA"),
                        "q1": result.get("q1", "NA"),
                        "q2": result.get("q2", "NA"),
                        "q3": result.get("q3", "NA"),
                    }
                )
    write_tsv(
        OUT,
        rows,
        [
            "dataset",
            "treatment",
            "status",
            "n_usable_trees",
            "coordinates_recoverable",
            "tip_labels_compatible",
            "mean_q_species",
            "mean_q_t_alt",
            "mean_q_other",
            "mean_q_unresolved",
            "astral_topology",
            "CU",
            "localPP",
            "q1",
            "q2",
            "q3",
        ],
    )
    print(f"Wrote {OUT}")


def normalize_filter_newick(newick: str) -> str:
    out = newick.replace(".._", "")
    out = out.replace("OverallCovFiltered", "")
    # Stage-1 canonical pseudo-t labels retain the original .fa marker.
    out = re_add_t_hapl_fa(out)
    return out


def re_add_t_hapl_fa(text: str) -> str:
    import re

    def repl(match: re.Match[str]) -> str:
        label = match.group(1)
        if label.endswith(".fa"):
            return label
        return label + ".fa"

    return re.sub(r"([A-Za-z0-9_.-]+_tHaplSubset)(?=[:),])", repl, text)


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
