#!/usr/bin/env python3
"""Analyze completed Stage 4E ASTER/ASTRAL4 species trees only."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage4d_core import (  # noqa: E402
    ROOT,
    RESULTS,
    named_splits,
    rf_comparison,
    focal_split_states,
    frozen_focal_groups,
    read_tsv,
)

BASE = ROOT / "empirical/neoaves_chr4"
TREE_DIR = BASE / "results/stage4e/trees"
OUT = BASE / "results/stage4e"


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def manifest_counts() -> dict[str, dict[str, str]]:
    p = RESULTS / "stage4e_input_manifest.tsv"
    if not p.exists():
        return {}
    return {row["treatment"]: row for row in read_tsv(p)}


def completed_trees() -> dict[str, str]:
    return {p.stem: p.read_text() for p in sorted(TREE_DIR.glob("*.nwk"))}


def split_set(newick: str) -> set[tuple[str, ...]]:
    return set(named_splits(newick)[1])


def parse_label(label: object) -> dict[str, object]:
    text = "" if label is None else str(label)
    out = {"localPP": "", "q1": "", "q2": "", "q3": "", "CULength": ""}
    for key in out:
        match = re.search(rf"{key}=([^,;:)\\]]+)", text)
        if match:
            out[key] = match.group(1)
    if text and not out["localPP"]:
        try:
            float(text)
            out["localPP"] = text
        except ValueError:
            pass
    return out


def comparison_rows(trees: dict[str, str]) -> list[dict[str, object]]:
    counts = manifest_counts()
    rows = []
    t0 = trees.get("T0")
    non = trees.get("NONCHR4")
    for treatment, newick in sorted(trees.items()):
        if treatment.startswith("J48_"):
            continue
        row = {"treatment": treatment}
        row["n_loci"] = counts.get(treatment, {}).get("n_loci", "")
        row["n_chr4_loci"] = counts.get(treatment, {}).get("n_chr4_loci", "")
        if t0:
            c = rf_comparison(newick, t0)
            row.update(rf_to_T0=c["rf"], nrf_to_T0=c["normalized_rf"], same_as_T0=c["rf"] == 0)
        if non:
            c = rf_comparison(newick, non)
            row.update(rf_to_NONCHR4=c["rf"], nrf_to_NONCHR4=c["normalized_rf"], same_as_NONCHR4=c["rf"] == 0)
        rows.append(row)
    return rows


def changed_split_rows(trees: dict[str, str]) -> list[dict[str, object]]:
    rows = []
    t0 = split_set(trees["T0"]) if "T0" in trees else set()
    non = split_set(trees["NONCHR4"]) if "NONCHR4" in trees else set()
    for treatment, newick in sorted(trees.items()):
        if treatment.startswith("J48_"):
            continue
        s = split_set(newick)
        rows.append(
            {
                "treatment": treatment,
                "splits_lost_from_T0": len(t0 - s) if t0 else "PENDING_CLUSTER_OUTPUT",
                "splits_gained_vs_T0": len(s - t0) if t0 else "PENDING_CLUSTER_OUTPUT",
                "splits_matching_NONCHR4": len(s & non) if non else "PENDING_CLUSTER_OUTPUT",
            }
        )
    return rows


def focal_rows(trees: dict[str, str]) -> list[dict[str, object]]:
    rows = []
    groups = frozen_focal_groups()
    baseline_t0 = {r["clade"]: r for r in focal_split_states(trees["T0"], groups)} if "T0" in trees else {}
    baseline_non = {r["clade"]: r for r in focal_split_states(trees["NONCHR4"], groups)} if "NONCHR4" in trees else {}
    for treatment, newick in sorted(trees.items()):
        for row in focal_split_states(newick, groups):
            row.update(
                treatment=treatment,
                CULength=parse_label(row.get("branch_support"))["CULength"],
                matches_T0=(row["state"] == baseline_t0.get(row["clade"], {}).get("state")) if baseline_t0 else "",
                matches_NONCHR4=(row["state"] == baseline_non.get(row["clade"], {}).get("state")) if baseline_non else "",
            )
            rows.append(row)
    return rows


def branch_rows(trees: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if "T0" not in trees:
        return [], []
    _, t0_splits = named_splits(trees["T0"])
    support_rows = []
    length_rows = []
    for treatment, newick in sorted(trees.items()):
        taxa, splits = named_splits(newick)
        for split, label in splits.items():
            if split not in t0_splits:
                continue
            parsed = parse_label(label)
            support_rows.append({"treatment": treatment, "split": "|".join(split), **parsed})
            length_rows.append(
                {
                    "treatment": treatment,
                    "split": "|".join(split),
                    "CULength": parsed["CULength"],
                    "t0_CULength": parse_label(t0_splits[split])["CULength"],
                }
            )
    return support_rows, length_rows


def random_status_rows() -> list[dict[str, object]]:
    p = RESULTS / "stage4e_random_matched_manifest.tsv"
    if not p.exists():
        return [{"regime": "ALL", "status": "PENDING_CLUSTER_OUTPUT"}]
    regimes = sorted({row["regime"] for row in read_tsv(p)})
    return [{"regime": regime, "status": "PENDING_CLUSTER_OUTPUT"} for regime in regimes]


def compare_published(trees: dict[str, str]) -> None:
    if "T0" not in trees:
        return
    published = BASE / "external/stiller2024/63K.tre"
    if not published.exists():
        return
    c = rf_comparison(trees["T0"], published.read_text())
    write_tsv(OUT / "stage4e_t0_vs_published63k.tsv", [{"comparison": "ASTER_T0_vs_published_63K", **c, "gate": "no"}], ["comparison", "rf", "normalized_rf", "branches_only_a", "branches_only_b", "gate"])


def write_figure_scaffold(trees: dict[str, str]) -> None:
    required = {"T0", "NONCHR4", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K"}
    script = OUT / "stage4e_figure_scaffold.R"
    script.write_text(
        """#!/usr/bin/env Rscript
required <- c("T0","NONCHR4","T_PNAS","T_STRUCT_PRIMARY","T_BLOCK_500K")
tree_dir <- "empirical/neoaves_chr4/results/stage4e/trees"
missing <- required[!file.exists(file.path(tree_dir, paste0(required, ".nwk")))]
if (length(missing)) {
  stop(paste("PENDING_CLUSTER_OUTPUT:", paste(missing, collapse=", ")))
}
message("Stage 4E figure scaffolding is ready; implement final panels after tree validation.")
"""
    )
    if not required.issubset(trees):
        print("PENDING_CLUSTER_OUTPUT: figure scaffold written but final figure refused until required trees exist")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-pending", action="store_true")
    args = parser.parse_args()
    trees = completed_trees()
    if not trees and not args.allow_pending:
        raise SystemExit("PENDING_CLUSTER_OUTPUT: no completed Stage 4E trees are available")
    write_tsv(OUT / "stage4e_tree_comparison.tsv", comparison_rows(trees), ["treatment", "n_loci", "n_chr4_loci", "rf_to_T0", "nrf_to_T0", "rf_to_NONCHR4", "nrf_to_NONCHR4", "same_as_T0", "same_as_NONCHR4"])
    write_tsv(OUT / "stage4e_changed_splits.tsv", changed_split_rows(trees), ["treatment", "splits_lost_from_T0", "splits_gained_vs_T0", "splits_matching_NONCHR4"])
    write_tsv(OUT / "stage4e_focal_relationships.tsv", focal_rows(trees) if trees else [], ["treatment", "clade", "state", "branch_support", "CULength", "matches_T0", "matches_NONCHR4", "n_C1", "n_C2", "n_S", "n_O"])
    support, lengths = branch_rows(trees)
    write_tsv(OUT / "stage4e_branch_support.tsv", support, ["treatment", "split", "localPP", "q1", "q2", "q3", "CULength"])
    write_tsv(OUT / "stage4e_branch_length_sensitivity.tsv", lengths, ["treatment", "split", "CULength", "t0_CULength"])
    write_tsv(OUT / "stage4e_random_matched_status.tsv", random_status_rows(), ["regime", "status"])
    compare_published(trees)
    write_figure_scaffold(trees)
    print("Stage 4E analysis tables updated from completed trees only.")


if __name__ == "__main__":
    main()
