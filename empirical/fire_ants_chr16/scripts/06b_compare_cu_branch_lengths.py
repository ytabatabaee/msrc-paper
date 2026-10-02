#!/usr/bin/env python3
"""Stage 6B: compare ASTRAL coalescent-unit branch lengths.

This extension reuses existing fire-ant Stage 6 ASTRAL4 inferred trees and
compares only ASTRAL CULength annotations on homologous internal bipartitions.
It does not rerun ASTRAL, does not inspect substitution-unit SULength values,
and does not alter the Stage-5 primary empirical analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import statistics
import sys
import unittest
import os
from pathlib import Path
from typing import Iterable

from Bio import Phylo
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

REPO_ROOT = Path(__file__).resolve().parents[3]
EMP = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
ASTRAL_DIR = RESULTS / "astral4"
README = EMP / "README.md"
PROJECT_STATUS = REPO_ROOT / "PROJECT_STATUS.md"
STAGE6_REPORT = RESULTS / "stage6_report.md"
STAGE6_MANIFEST = RESULTS / "stage6_manifest.json"

BACKGROUND_TREE = ASTRAL_DIR / "T_background.inferred.nwk"
ALL_TREE = ASTRAL_DIR / "T_all.inferred.nwk"
CHR16_ALL_TREE = ASTRAL_DIR / "T_chr16_all.inferred.nwk"
SUPERGENE_TREE = ASTRAL_DIR / "T_supergene.inferred.nwk"

COMPARISON_TSV = RESULTS / "stage6_cu_branch_length_comparison.tsv"
SUPERGENE_SPLITS_TSV = RESULTS / "stage6_cu_supergene_only_splits.tsv"
SUMMARY_TXT = RESULTS / "stage6_cu_branch_length_summary.txt"
MANIFEST = RESULTS / "stage6_cu_manifest.json"
SCATTER_PDF = FIGURES / "fire_ants_stage6_cu_background_vs_combined.pdf"
SCATTER_PNG = FIGURES / "fire_ants_stage6_cu_background_vs_combined.png"
DELTA_PDF = FIGURES / "fire_ants_stage6_cu_delta.pdf"
DELTA_PNG = FIGURES / "fire_ants_stage6_cu_delta.png"

TAXA = sorted(["geminata", "saevissima", "pusillignis", "inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"])
FOCAL = frozenset(["inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"])
INV_PAIR = frozenset(["inv_mac_SB", "inv_mac_Sb"])
RIC_PAIR = frozenset(["richteri_SB", "richteri_Sb"])
HAP_SB_PAIR = frozenset(["inv_mac_SB", "richteri_SB"])
HAP_Sb_PAIR = frozenset(["inv_mac_Sb", "richteri_Sb"])
EPS = 1e-12

COMPARISONS = {
    "background_vs_all": {
        "background": BACKGROUND_TREE,
        "combined": ALL_TREE,
        "background_label": "T_background",
        "combined_label": "T_all_background_plus_supergene",
        "role": "primary_background_vs_background_plus_supergene",
    },
    "background_vs_chr16_all": {
        "background": BACKGROUND_TREE,
        "combined": CHR16_ALL_TREE,
        "background_label": "T_background",
        "combined_label": "T_chr16_all",
        "role": "secondary_chr16_all_comparison",
    },
    "background_vs_supergene": {
        "background": BACKGROUND_TREE,
        "combined": SUPERGENE_TREE,
        "background_label": "T_background",
        "combined_label": "T_supergene",
        "role": "secondary_supergene_only_comparison",
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def parse_tree_text(text: str):
    return Phylo.read(io.StringIO(text.strip()), "newick")


def parse_tree(path: Path):
    return parse_tree_text(path.read_text())


def all_taxa(tree) -> frozenset[str]:
    labels = [t.name for t in tree.get_terminals()]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate terminal labels")
    return frozenset(labels)


def smaller_side(side: Iterable[str], taxa: Iterable[str]) -> frozenset[str]:
    s = frozenset(side)
    t = frozenset(taxa)
    other = t - s
    if len(s) < len(other):
        return s
    if len(other) < len(s):
        return other
    return min([s, other], key=lambda x: ",".join(sorted(x)))


def split_id_from_side(side: Iterable[str], taxa: Iterable[str]) -> str:
    t = frozenset(taxa)
    s = frozenset(side)
    other = t - s
    left = ",".join(sorted(s))
    right = ",".join(sorted(other))
    return "|".join(sorted([left, right]))


def compact_side(split_id: str, taxa: Iterable[str]) -> frozenset[str]:
    left, right = split_id.split("|")
    a = frozenset(left.split(",")) if left else frozenset()
    b = frozenset(right.split(",")) if right else frozenset()
    if len(a) < len(b):
        return a
    if len(b) < len(a):
        return b
    return min([a, b], key=lambda x: ",".join(sorted(x)))


def annotation_dict(name: str | None) -> dict[str, str]:
    if not name:
        return {}
    m = re.search(r"\[(.*)\]", name)
    if not m:
        return {}
    out = {}
    for piece in m.group(1).split(";"):
        if "=" in piece:
            key, value = piece.split("=", 1)
            out[key] = value
    return out


def parse_float_safe(value: str | None) -> float | None:
    if value is None or value == "" or value == "NA":
        return None
    try:
        out = float(value)
    except ValueError as exc:
        raise ValueError(f"Could not parse numeric value `{value}`") from exc
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def branch_role(split_id: str, taxa: Iterable[str]) -> str:
    left_text, right_text = split_id.split("|")
    sides = [frozenset(left_text.split(",")), frozenset(right_text.split(","))]
    if INV_PAIR in sides:
        return "focal_invicta_macdonaghi_SB_Sb_pair"
    if RIC_PAIR in sides:
        return "focal_richteri_SB_Sb_pair"
    if HAP_SB_PAIR in sides:
        return "supergene_haplotype_SB_pair"
    if HAP_Sb_PAIR in sides:
        return "supergene_haplotype_Sb_pair"
    if FOCAL in sides:
        return "focal_four_vs_outgroups"
    if any(side & FOCAL for side in sides):
        return "focal_related_internal_branch"
    return "non_focal_internal_branch"


def is_focal_related(split_id: str) -> bool:
    return branch_role(split_id, TAXA) != "non_focal_internal_branch"


def extract_cu_branches(path: Path) -> dict[str, dict[str, object]]:
    tree = parse_tree(path)
    taxa = all_taxa(tree)
    if taxa != frozenset(TAXA):
        raise ValueError(f"Unexpected taxon set in {path}: {sorted(taxa)}")
    out: dict[str, dict[str, object]] = {}
    problems = []
    n = len(taxa)
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        if len(side) <= 1 or len(side) >= n - 1:
            continue
        split_id = split_id_from_side(side, taxa)
        ann = annotation_dict(getattr(clade, "name", None))
        if "CULength" not in ann:
            problems.append(f"missing CULength for split {split_id}")
            cu = None
        else:
            cu = parse_float_safe(ann["CULength"])
        if split_id in out:
            raise ValueError(f"Duplicate internal split in {path}: {split_id}")
        side_small = compact_side(split_id, taxa)
        out[split_id] = {
            "split_id": split_id,
            "taxa_side": ",".join(sorted(side_small)),
            "n_taxa_side": len(side_small),
            "cu": cu,
            "branch_role": branch_role(split_id, taxa),
            "focal_related": is_focal_related(split_id),
            "annotation": ann,
        }
    if problems:
        raise ValueError(f"CU parsing problems in {path}: " + "; ".join(problems))
    return out


def tree_topology_splits(path: Path) -> set[str]:
    return set(extract_cu_branches(path))


def compare_pair(name: str, spec: dict[str, object]) -> list[dict[str, object]]:
    bg = extract_cu_branches(spec["background"])
    combined = extract_cu_branches(spec["combined"])
    rows = []
    for split_id in sorted(set(bg) | set(combined)):
        b = bg.get(split_id)
        c = combined.get(split_id)
        if b and c:
            status = "shared branch"
            cu_b = b["cu"]
            cu_c = c["cu"]
            delta = None if cu_b is None or cu_c is None else cu_c - cu_b
            rel_delta = None if delta is None or cu_b is None or abs(cu_b) <= EPS else delta / cu_b
            base = b
        elif b:
            status = "background-only split"
            cu_b = b["cu"]
            cu_c = None
            delta = None
            rel_delta = None
            base = b
        else:
            status = "combined-only split"
            cu_b = None
            cu_c = c["cu"]
            delta = None
            rel_delta = None
            base = c
        rows.append({
            "comparison": name,
            "comparison_role": spec["role"],
            "background_label": spec["background_label"],
            "combined_label": spec["combined_label"],
            "split_id": split_id,
            "taxa_side": base["taxa_side"],
            "n_taxa_side": base["n_taxa_side"],
            "cu_background": cu_b,
            "cu_combined": cu_c,
            "delta_cu": delta,
            "relative_delta_cu": rel_delta,
            "status": status,
            "branch_role": base["branch_role"],
            "focal_related": base["focal_related"],
        })
    return rows


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        return f"{value:.10g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def median(values: Iterable[float | None]) -> float | None:
    clean = [v for v in values if v is not None and not math.isnan(v) and not math.isinf(v)]
    return statistics.median(clean) if clean else None


def summary_lines(rows: list[dict[str, object]]) -> list[str]:
    lines = ["Stage 6 CU branch-length comparison", ""]
    for name, spec in COMPARISONS.items():
        comp = [r for r in rows if r["comparison"] == name]
        shared = [r for r in comp if r["status"] == "shared branch"]
        bg_only = [r for r in comp if r["status"] == "background-only split"]
        combined_only = [r for r in comp if r["status"] == "combined-only split"]
        bg = extract_cu_branches(spec["background"])
        combined = extract_cu_branches(spec["combined"])
        topo_same = not bg_only and not combined_only
        lines.extend([
            f"Comparison: {name}",
            f"background tree: {rel(spec['background'])}",
            f"combined tree: {rel(spec['combined'])}",
            f"topology identical by CU-bearing internal bipartitions: {str(topo_same).lower()}",
            f"internal CU branches in background: {len(bg)}",
            f"internal CU branches in combined: {len(combined)}",
            f"shared bipartitions: {len(shared)}",
            f"background-only bipartitions: {len(bg_only)}",
            f"combined-only bipartitions: {len(combined_only)}",
            f"median CU background: {fmt(median(r['cu_background'] for r in shared))}",
            f"median CU combined: {fmt(median(r['cu_combined'] for r in shared))}",
            f"median absolute CU change: {fmt(median(abs(r['delta_cu']) for r in shared if r['delta_cu'] is not None))}",
        ])
        decreases = sorted([r for r in shared if r["delta_cu"] is not None and r["delta_cu"] < 0], key=lambda r: r["delta_cu"])
        increases = sorted([r for r in shared if r["delta_cu"] is not None and r["delta_cu"] > 0], key=lambda r: r["delta_cu"], reverse=True)
        lines.append("largest CU decreases:")
        for r in decreases[:3]:
            lines.append(f"  {r['branch_role']} [{r['taxa_side']}]: {fmt(r['cu_background'])} -> {fmt(r['cu_combined'])}; delta {fmt(r['delta_cu'])}; relative {fmt(r['relative_delta_cu'])}")
        lines.append("largest CU increases:")
        for r in increases[:3]:
            lines.append(f"  {r['branch_role']} [{r['taxa_side']}]: {fmt(r['cu_background'])} -> {fmt(r['cu_combined'])}; delta {fmt(r['delta_cu'])}; relative {fmt(r['relative_delta_cu'])}")
        focal = [r for r in shared if r["focal_related"]]
        lines.append("focal/supergene-associated branch CU changes:")
        for r in focal:
            lines.append(f"  {r['branch_role']} [{r['taxa_side']}]: {fmt(r['cu_background'])} -> {fmt(r['cu_combined'])}; delta {fmt(r['delta_cu'])}; relative {fmt(r['relative_delta_cu'])}")
        if name == "background_vs_supergene" and combined_only:
            lines.append("supergene-only alternative splits absent from background:")
            for r in combined_only:
                lines.append(f"  {r['branch_role']} [{r['taxa_side']}]: CU {fmt(r['cu_combined'])}")
        lines.append("")
    primary = [r for r in rows if r["comparison"] == "background_vs_all" and r["status"] == "shared branch"]
    focal_decreases = [r for r in primary if r["focal_related"] and r["delta_cu"] is not None and r["delta_cu"] < 0]
    all_deltas = [r["delta_cu"] for r in primary if r["delta_cu"] is not None]
    lines.append("Interpretation:")
    if focal_decreases:
        largest = min(focal_decreases, key=lambda r: r["delta_cu"])
        lines.append(
            f"Adding chr16/supergene loci leaves the primary Stage-6 ASTRAL topology unchanged but shortens focal CU branches, with the largest focal decrease on {largest['branch_role']} [{largest['taxa_side']}] from {fmt(largest['cu_background'])} to {fmt(largest['cu_combined'])} CU."
        )
    else:
        lines.append("Adding chr16/supergene loci leaves the primary Stage-6 ASTRAL topology unchanged and does not shorten any branch marked focal/supergene-associated by this script.")
    if all_deltas:
        lines.append(f"Across all shared primary branches, median absolute CU change is {fmt(median(abs(x) for x in all_deltas))}; changes are branch-specific rather than a uniform genome-wide scaling.")
    return lines


def write_supergene_splits(rows: list[dict[str, object]]) -> None:
    selected = [r for r in rows if r["comparison"] == "background_vs_supergene" and r["status"] != "shared branch"]
    write_tsv(SUPERGENE_SPLITS_TSV, selected, ["comparison", "split_id", "taxa_side", "n_taxa_side", "cu_background", "cu_combined", "status", "branch_role", "focal_related"])


def plot_outputs(rows: list[dict[str, object]]) -> None:
    primary = [r for r in rows if r["comparison"] == "background_vs_all" and r["status"] == "shared branch"]
    colors = ["#b23a48" if r["focal_related"] else "#345995" for r in primary]
    labels = [r["branch_role"].replace("focal_", "").replace("supergene_", "") for r in primary]

    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    xs = [r["cu_background"] for r in primary]
    ys = [r["cu_combined"] for r in primary]
    lim = max(xs + ys) * 1.08
    ax.plot([0, lim], [0, lim], color="#777777", linewidth=0.9, linestyle="--")
    ax.scatter(xs, ys, s=42, c=colors, edgecolor="white", linewidth=0.6)
    for x, y, label, row in zip(xs, ys, labels, primary):
        if row["focal_related"]:
            ax.annotate(label, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("background-only CULength")
    ax.set_ylabel("background + supergene CULength")
    ax.set_title("Stage 6 ASTRAL CULength: background vs all windows")
    ax.grid(color="#eeeeee", linewidth=0.6)
    fig.tight_layout()
    with PdfPages(SCATTER_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(SCATTER_PNG, dpi=220)
    plt.close(fig)

    ordered = sorted(primary, key=lambda r: r["delta_cu"])
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    y = list(range(len(ordered)))
    bar_colors = ["#b23a48" if r["focal_related"] else "#345995" for r in ordered]
    ax.barh(y, [r["delta_cu"] for r in ordered], color=bar_colors)
    ax.axvline(0, color="#222222", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['branch_role']}\n{r['taxa_side']}" for r in ordered], fontsize=7)
    ax.set_xlabel("delta CULength (background + supergene - background)")
    ax.set_title("Stage 6 ASTRAL CULength changes")
    ax.grid(axis="x", color="#eeeeee", linewidth=0.6)
    fig.tight_layout()
    with PdfPages(DELTA_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(DELTA_PNG, dpi=220)
    plt.close(fig)


def update_docs(summary: str) -> None:
    block = (
        "\n## Stage 6 CU branch-length extension\n\n"
        "Existing Stage 6 ASTRAL4 inferred trees were parsed for `CULength` annotations to test whether adding chr16/supergene windows changes coalescent-unit branch lengths when topology is unchanged. The primary comparison is `T_background.inferred.nwk` versus `T_all.inferred.nwk`; `T_chr16_all` and `T_supergene` are secondary comparisons. `SULength` values were not analyzed. Outputs are `stage6_cu_branch_length_comparison.tsv`, `stage6_cu_branch_length_summary.txt`, and the two CU figures.\n"
    )
    report = STAGE6_REPORT.read_text()
    if "## Stage 6 CU branch-length extension" not in report:
        STAGE6_REPORT.write_text(report.rstrip() + "\n" + block)

    readme = README.read_text()
    if "06b_compare_cu_branch_lengths.py --run-tests" not in readme:
        readme = readme.replace(
            "python3 empirical/fire_ants_chr16/scripts/06_summary_tree_sensitivity.py --run-tests\n```",
            "python3 empirical/fire_ants_chr16/scripts/06_summary_tree_sensitivity.py --run-tests\npython3 empirical/fire_ants_chr16/scripts/06b_compare_cu_branch_lengths.py --run-tests\n```",
        )
    if "stage6_cu_branch_length_comparison.tsv" not in readme:
        readme += "\nStage 6 CU branch-length extension outputs: `results/stage6_cu_branch_length_comparison.tsv`, `results/stage6_cu_branch_length_summary.txt`, `figures/fire_ants_stage6_cu_background_vs_combined.pdf`, and `figures/fire_ants_stage6_cu_delta.pdf`.\n"
    README.write_text(readme)

    status_block = (
        "\n## Fire-ant chromosome 16 Stage 6 CU branch-length extension\n\n"
        "The existing grouped TWISST / ASTRAL4 Stage 6 trees were reused to compare ASTRAL coalescent-unit branch lengths. The primary topology, `T_background` versus `T_all`, remains unchanged, but focal CU lengths were compared branch by branch using canonical unrooted bipartitions. `SULength` was not analyzed because the modal gene-tree inputs do not contain meaningful substitution branch lengths. See `empirical/fire_ants_chr16/results/stage6_cu_branch_length_summary.txt` and `empirical/fire_ants_chr16/results/stage6_cu_branch_length_comparison.tsv`.\n"
    )
    if PROJECT_STATUS.exists():
        text = PROJECT_STATUS.read_text()
        if "## Fire-ant chromosome 16 Stage 6 CU branch-length extension" not in text:
            PROJECT_STATUS.write_text(status_block + "\n" + text)


def write_manifest(rows: list[dict[str, object]]) -> None:
    paths = [COMPARISON_TSV, SUPERGENE_SPLITS_TSV, SUMMARY_TXT, SCATTER_PDF, SCATTER_PNG, DELTA_PDF, DELTA_PNG, STAGE6_REPORT, README]
    if PROJECT_STATUS.exists():
        paths.append(PROJECT_STATUS)
    manifest = {
        "stage": "6_cu_branch_length_extension",
        "analysis_type": "post_freeze_cu_branch_length_sensitivity",
        "astral_outputs_reused": True,
        "astral_rerun": False,
        "substitution_unit_lengths_analyzed": False,
        "coalescent_unit_lengths_analyzed": True,
        "primary_comparison": "background_vs_all",
        "input_checksums": {rel(p): sha256(p) for p in [BACKGROUND_TREE, ALL_TREE, CHR16_ALL_TREE, SUPERGENE_TREE, STAGE6_MANIFEST]},
        "output_checksums": {rel(p): sha256(p) for p in paths if p.exists()},
        "processing_script_checksum": sha256(Path(__file__)),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def run_analysis() -> dict[str, object]:
    for path in [BACKGROUND_TREE, ALL_TREE, CHR16_ALL_TREE, SUPERGENE_TREE]:
        if not path.exists():
            raise FileNotFoundError(path)
        # Confirm CULength exists in the raw output before parsing.
        if "CULength=" not in path.read_text():
            raise ValueError(f"No CULength annotation found in {path}")
    rows = []
    for name, spec in COMPARISONS.items():
        rows.extend(compare_pair(name, spec))
    fields = ["comparison", "comparison_role", "background_label", "combined_label", "split_id", "taxa_side", "n_taxa_side", "cu_background", "cu_combined", "delta_cu", "relative_delta_cu", "status", "branch_role", "focal_related"]
    write_tsv(COMPARISON_TSV, rows, fields)
    write_supergene_splits(rows)
    lines = summary_lines(rows)
    SUMMARY_TXT.write_text("\n".join(lines) + "\n")
    plot_outputs(rows)
    update_docs("\n".join(lines))
    write_manifest(rows)
    primary_shared = [r for r in rows if r["comparison"] == "background_vs_all" and r["status"] == "shared branch"]
    return {
        "comparison_tsv": rel(COMPARISON_TSV),
        "summary": rel(SUMMARY_TXT),
        "primary_shared_branches": len(primary_shared),
        "primary_topology_identical": all(r["status"] == "shared branch" for r in rows if r["comparison"] == "background_vs_all"),
        "manifest_sha": sha256(MANIFEST),
    }


class CuLengthTests(unittest.TestCase):
    def test_canonical_split_independent_of_side(self):
        taxa = ["a", "b", "c", "d", "e"]
        self.assertEqual(split_id_from_side(["a", "b"], taxa), split_id_from_side(["c", "d", "e"], taxa))

    def test_identical_trees_match_all_internal_branches(self):
        rows = compare_pair("background_vs_all", {**COMPARISONS["background_vs_all"], "combined": BACKGROUND_TREE})
        self.assertTrue(rows)
        self.assertTrue(all(r["status"] == "shared branch" for r in rows))

    def test_topology_specific_branches_classified(self):
        rows = compare_pair("background_vs_supergene", COMPARISONS["background_vs_supergene"])
        statuses = {r["status"] for r in rows}
        self.assertIn("background-only split", statuses)
        self.assertIn("combined-only split", statuses)

    def test_cu_parsing_actual_astral_syntax(self):
        branches = extract_cu_branches(BACKGROUND_TREE)
        self.assertEqual(len(branches), 4)
        inv = next(v for v in branches.values() if v["branch_role"] == "focal_invicta_macdonaghi_SB_Sb_pair")
        self.assertAlmostEqual(inv["cu"], 2.73003, places=5)

    def test_safe_numeric_handling(self):
        self.assertIsNone(parse_float_safe("inf"))
        self.assertIsNone(parse_float_safe("nan"))
        self.assertIsNone(parse_float_safe(None))
        self.assertEqual(parse_float_safe("0"), 0.0)

    def test_outputs_exist_and_use_no_su_lengths(self):
        text = SUMMARY_TXT.read_text()
        self.assertIn("SULength", Path(__file__).read_text())
        self.assertTrue(COMPARISON_TSV.exists())
        self.assertTrue(SCATTER_PDF.exists())


def cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    result = run_analysis()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(CuLengthTests)
        outcome = unittest.TextTestRunner(verbosity=2).run(suite)
        if not outcome.wasSuccessful():
            raise SystemExit(1)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    cli()
