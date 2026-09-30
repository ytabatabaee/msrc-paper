#!/usr/bin/env python3
"""Stage 6: Atlantic cod quartet-support profiles and ASTRAL4 sensitivity.

Post-freeze analysis only. Verifies Stage-5 frozen inputs, reports local
quartet-support tracks derived from Stage 4A, and runs ASTER/ASTRAL4 population
summary-tree sensitivity analyses on the published 250-kb local population-tree
windows.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import Bio
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMP_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP_ROOT / "results"
FIGURES = EMP_ROOT / "figures"
ASTRAL_DIR = RESULTS / "astral4"
INTERMEDIATE = DATA_ROOT / "intermediate" / "astral4"

STAGE5_MANIFEST = RESULTS / "stage5_final_manifest.json"
COD_WINDOWS = DATA_ROOT / "processed" / "cod_window_trees.tsv"
STAGE4A_SCORES = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"
BASELINE_TREE = DATA_ROOT / "processed" / "baseline_population_tree.nwk"
BASELINE_12 = DATA_ROOT / "processed" / "baseline_population_tree_12taxa.nwk"
REGIONS = DATA_ROOT / "metadata" / "region_manifest.tsv"

LOCAL_SUPPORT = DATA_ROOT / "processed" / "stage6_local_quartet_support.tsv"
LOCAL_SUPPORT_SUMMARY = RESULTS / "stage6_local_quartet_support_summary.tsv"
QUARTET_SUPPORT_FIG_PDF = FIGURES / "atlantic_cod_quartet_support_tracks.pdf"
QUARTET_SUPPORT_FIG_PNG = FIGURES / "atlantic_cod_quartet_support_tracks.png"
ASTER_PROVENANCE = RESULTS / "stage6_aster_provenance.md"
BRANCH_SUPPORT = RESULTS / "stage6_astral4_baseline_branch_support.tsv"
TREATMENT_SUMMARY = RESULTS / "stage6_astral4_treatment_summary.tsv"
DOWNWEIGHTING = RESULTS / "stage6_astral4_downweighting.tsv"
SPLIT_CHANGES = RESULTS / "stage6_astral4_split_changes.tsv"
ASTRAL_FIG_PDF = FIGURES / "atlantic_cod_astral4_sensitivity.pdf"
ASTRAL_FIG_PNG = FIGURES / "atlantic_cod_astral4_sensitivity.png"
BRANCH_FIG_PDF = FIGURES / "atlantic_cod_astral4_branch_support.pdf"
BRANCH_FIG_PNG = FIGURES / "atlantic_cod_astral4_branch_support.png"
REPORT = RESULTS / "stage6_report.md"
METHODS_TEXT = RESULTS / "stage6_methods_text.md"
RESULTS_TEXT = RESULTS / "stage6_results_text.md"
MANIFEST = RESULTS / "stage6_manifest.json"

LGS = ("LG01", "LG02", "LG07", "LG12")
THREADS = int(os.environ.get("ASTRAL4_THREADS", "2"))
M_VALUES = [1, 2, 5, 10, 20, 40, "ALL"]
FINITE_REPS = int(os.environ.get("STAGE6_DOWNWEIGHT_REPS", "100"))
SEED_BASE = 20260929

FROZEN_EXPECTED = {
    "data/atlantic_cod/processed/stage2_structural_quartet_candidates.tsv": "303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c",
    "data/atlantic_cod/processed/stage3_baseline_quartets.tsv": "5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1",
    "data/atlantic_cod/processed/cod_window_trees.tsv": "26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b",
    "data/atlantic_cod/processed/stage4a_local_quartets.tsv": "8ae82f803294a47854a8eecf907dc2684dc89cf6407362beaf9a4e7616d306f2",
    "data/atlantic_cod/processed/stage4a_window_topology_scores.tsv": "3ef09caa4edddd247ad986ff168a4fc86b79d66559293bec51d291ee8cf9328e",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: fmt(row.get(k, "")) for k in fields})


def fmt(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        if math.isnan(v): return "NA"
        return f"{v:.6f}"
    return str(v)


def as_bool(s: str | bool) -> bool:
    if isinstance(s, bool): return s
    if s.lower() == "true": return True
    if s.lower() == "false": return False
    raise ValueError(s)


def verify_frozen() -> dict[str, str]:
    manifest = json.loads(STAGE5_MANIFEST.read_text())
    observed = {}
    for rel, expected in FROZEN_EXPECTED.items():
        path = REPO_ROOT / rel
        got = sha256(path)
        if got != expected:
            raise ValueError(f"Frozen checksum mismatch for {rel}: {got} != {expected}")
        if manifest.get("frozen_input_checksums", {}).get(rel) != expected:
            raise ValueError(f"Stage-5 manifest checksum mismatch for {rel}")
        observed[rel] = got
    return observed


def parse_tree(newick: str):
    return Phylo.read(io.StringIO(newick.strip()), "newick")


def taxa(tree) -> list[str]:
    labels = [t.name for t in tree.get_terminals()]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate taxon labels")
    return sorted(labels)


def canonical_split(side: Iterable[str], all_taxa: Iterable[str]) -> str:
    s = frozenset(side)
    allset = frozenset(all_taxa)
    other = allset - s
    a = ",".join(sorted(s)); b = ",".join(sorted(other))
    return "|".join(sorted([a, b]))


def splits(tree, include_trivial: bool = False) -> set[str]:
    all_taxa = taxa(tree)
    n = len(all_taxa)
    out = set()
    for clade in tree.find_clades():
        side = [t.name for t in clade.get_terminals()]
        if not include_trivial and (len(side) <= 1 or len(side) >= n - 1):
            continue
        if len(side) == n:
            continue
        out.add(canonical_split(side, all_taxa))
    return out


def rf_to_baseline(tree_newick: str, baseline_newick: str) -> dict[str, object]:
    t = parse_tree(tree_newick); b = parse_tree(baseline_newick)
    st = splits(t); sb = splits(b)
    rf = len(st - sb) + len(sb - st)
    denom = len(st) + len(sb)
    return {"unrooted_RF": rf, "normalized_RF": rf / denom if denom else 0.0, "n_baseline_internal_splits": len(sb), "n_recovered_baseline_internal_splits": len(st & sb), "fraction_baseline_splits_recovered": len(st & sb) / len(sb) if sb else 1.0, "exact_topology_match": st == sb}


def annotation_dict(name: str | None) -> dict[str, float]:
    if not name:
        return {}
    m = re.search(r"\[(.*)\]", name)
    if not m:
        return {}
    out = {}
    for piece in m.group(1).split(";"):
        if "=" in piece:
            k, v = piece.split("=", 1)
            try:
                out[k] = float(v)
            except ValueError:
                pass
    return out


def annotated_branch_support(newick: str) -> dict[str, dict[str, float]]:
    tree = parse_tree(newick)
    all_taxa = taxa(tree)
    n = len(all_taxa)
    out = {}
    for clade in tree.find_clades():
        side = [t.name for t in clade.get_terminals()]
        if len(side) <= 1 or len(side) >= n - 1:
            continue
        ann = annotation_dict(getattr(clade, "name", None))
        if ann:
            out[canonical_split(side, all_taxa)] = ann
    return out


def find_astral4() -> Path:
    env = os.environ.get("ASTRAL4_BIN")
    candidates = []
    if env: candidates.append(Path(env))
    p = shutil.which("astral4")
    if p: candidates.append(Path(p))
    candidates.append(REPO_ROOT / "external" / "ASTER" / "bin" / "astral4")
    for c in candidates:
        if c.exists() and os.access(c, os.X_OK):
            return c.resolve()
    raise RuntimeError("No working astral4 found in $ASTRAL4_BIN, PATH, or external/ASTER/bin/astral4")


def run_cmd(cmd: list[str], log: Path) -> float:
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    with log.open("w") as err:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=err, check=True)
    return time.time() - start


def astral_infer(astral: Path, input_file: Path, output: Path, log: Path) -> tuple[str, float, str]:
    cmd = [str(astral), "-R", "-u", "2", "-t", str(THREADS), "-i", str(input_file), "-o", str(output)]
    if output.exists() and output.read_text().strip():
        return " ".join(cmd), 0.0, output.read_text().strip()
    rt = run_cmd(cmd, log)
    return " ".join(cmd), rt, output.read_text().strip()


def astral_score(astral: Path, input_file: Path, baseline: Path, output: Path, log: Path) -> tuple[str, float, str]:
    cmd = [str(astral), "-C", "-u", "2", "-t", str(THREADS), "-c", str(baseline), "-i", str(input_file), "-o", str(output)]
    if output.exists() and output.read_text().strip():
        return " ".join(cmd), 0.0, output.read_text().strip()
    rt = run_cmd(cmd, log)
    return " ".join(cmd), rt, output.read_text().strip()


def load_windows() -> list[dict[str, object]]:
    rows = []
    for r in read_tsv(COD_WINDOWS):
        row = dict(r)
        for k in ["start", "end", "midpoint"]: row[k] = int(r[k])
        rows.append(row)
    return sorted(rows, key=lambda x: (x["lg"], x["start"]))


def write_treefile(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for r in sorted(rows, key=lambda x: (x["lg"], int(x["start"]))):
            f.write(str(r["tree_newick"]).strip() + "\n")


def validate_treefile(path: Path, expected_taxa: list[str]) -> int:
    count = 0
    for line in path.read_text().splitlines():
        if not line.strip(): continue
        if taxa(parse_tree(line)) != expected_taxa:
            raise ValueError(f"Taxon mismatch in {path}")
        count += 1
    return count


def prune_baseline_12(expected_taxa: list[str]) -> str:
    tree = Phylo.read(BASELINE_TREE, "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in expected_taxa:
            tree.prune(terminal)
    if taxa(tree) != expected_taxa:
        raise ValueError("Pruned baseline taxa mismatch")
    s = io.StringIO(); Phylo.write(tree, s, "newick")
    newick = s.getvalue().strip()
    BASELINE_12.write_text(newick + "\n")
    # fully resolved check: n-3 unrooted splits
    if len(splits(parse_tree(newick))) != len(expected_taxa) - 3:
        raise ValueError("Pruned baseline tree is not fully resolved")
    return newick


def create_treatments(windows: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    inside = [r for r in windows if as_bool(r["inside_inversion"])]
    outside = [r for r in windows if not as_bool(r["inside_inversion"]) and not as_bool(r["overlaps_inversion"])]
    treatments = {
        "T0_all": windows,
        "T1_outside": outside,
        "T2_inside": inside,
    }
    for lg in LGS:
        treatments[f"LG{lg[-2:]}_inside"] = [r for r in inside if r["lg"] == lg]
        treatments[f"T_drop_{lg}"] = [r for r in windows if not (r["lg"] == lg and as_bool(r["inside_inversion"]))]
    return treatments


def write_treatment_files(treatments: dict[str, list[dict[str, object]]], expected_taxa: list[str]) -> dict[str, Path]:
    paths = {}
    aliases = {"T0_all": "all_426", "T1_outside": "fully_outside", "T2_inside": "fully_inside", "LGLG01_inside": "LG01_inside"}
    for name, rows in treatments.items():
        if name.startswith("LG") and name.endswith("_inside"):
            fname = name.replace("LG01", "LG01")
        standard = {"T0_all": "all_426.trees", "T1_outside": "fully_outside.trees", "T2_inside": "fully_inside.trees"}.get(name, f"{name}.trees")
        if name in [f"LG{lg[-2:]}_inside" for lg in LGS]:
            standard = name.replace("LG01", "LG01") + ".trees"
        path = INTERMEDIATE / standard
        write_treefile(path, rows)
        validate_treefile(path, expected_taxa)
        paths[name] = path
    # friendly copies for required exact names
    for lg in LGS:
        key = f"LG{lg[-2:]}_inside"
        if key in paths:
            required = INTERMEDIATE / f"{lg}_inside.trees"
            if required != paths[key]:
                required.write_text(paths[key].read_text())
                paths[key] = required
    return paths


def local_support_table() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = []
    for r in read_tsv(STAGE4A_SCORES):
        qA = float(r["fraction_arrangement"]); qB = float(r["fraction_baseline"]); q3 = float(r["fraction_third"])
        if int(r["n_resolved"]) > 0 and abs((qA + qB + q3) - 1.0) > 1e-5:
            raise ValueError(f"Local support fractions do not sum to 1 for {r['window_id']}")
        rows.append({"lg": r["lg"], "window_id": r["window_id"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "inside_inversion": r["inside_inversion"], "overlaps_inversion": r["overlaps_inversion"], "n_informative_quartets": r["n_informative_quartets"], "n_resolved": r["n_resolved"], "q_arrangement": qA, "q_baseline": qB, "q_third": q3, "D_arrangement_minus_baseline": r["D_arrangement_minus_baseline"]})
    summary = []
    for lg in LGS:
        for cls, pred in [("inside", lambda x: as_bool(x["inside_inversion"])), ("outside", lambda x: (not as_bool(x["inside_inversion"]) and not as_bool(x["overlaps_inversion"])))]:
            sel = [r for r in rows if r["lg"] == lg and pred(r)]
            summary.append({"lg": lg, "region_class": cls,
                            "mean_q_arrangement": statistics.fmean(float(r["q_arrangement"]) for r in sel),
                            "mean_q_baseline": statistics.fmean(float(r["q_baseline"]) for r in sel),
                            "mean_q_third": statistics.fmean(float(r["q_third"]) for r in sel),
                            "median_q_arrangement": statistics.median(float(r["q_arrangement"]) for r in sel),
                            "median_q_baseline": statistics.median(float(r["q_baseline"]) for r in sel),
                            "median_q_third": statistics.median(float(r["q_third"]) for r in sel)})
    return rows, summary


def plot_local_support(rows: list[dict[str, object]], regions: dict[str, dict[str, int]]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 8.0), sharey=True)
    colors = {"q_arrangement": "#2ca02c", "q_baseline": "#9467bd", "q_third": "#8c564b"}
    labels = {"q_arrangement": "arrangement", "q_baseline": "baseline", "q_third": "third"}
    for ax, lg in zip(axes, LGS):
        lgrows = sorted([r for r in rows if r["lg"] == lg], key=lambda r: int(r["midpoint"]))
        x = [int(r["midpoint"]) / 1_000_000 for r in lgrows]
        ax.axvspan(regions[lg]["start"] / 1_000_000, regions[lg]["end"] / 1_000_000, color="#d9d9d9", alpha=0.65)
        for key in ["q_arrangement", "q_baseline", "q_third"]:
            ax.plot(x, [float(r[key]) for r in lgrows], lw=1.0, marker="o", ms=2, color=colors[key], label=labels[key])
        ax.set_ylim(-0.03, 1.03); ax.set_ylabel(lg)
        ax.axhline(0.5, color="black", lw=0.4, alpha=0.5)
    axes[0].legend(loc="upper right", fontsize=8, ncol=3)
    axes[-1].set_xlabel("genomic position (Mb)")
    fig.suptitle("Atlantic cod local quartet-support fractions", y=0.995)
    fig.tight_layout()
    fig.savefig(QUARTET_SUPPORT_FIG_PDF)
    fig.savefig(QUARTET_SUPPORT_FIG_PNG, dpi=300)
    plt.close(fig)


def regions() -> dict[str, dict[str, int]]:
    out = {}
    for r in read_tsv(REGIONS):
        if r["lg"] in LGS:
            out[r["lg"]] = {"start": int(r["start"]), "end": int(r["end"])}
    return out


def treatment_count_expected(name: str, rows: list[dict[str, object]]) -> None:
    expected = {"T0_all": 426, "T1_outside": 245, "T2_inside": 176, "LG01_inside": 66, "LG02_inside": 21, "LG07_inside": 36, "LG12_inside": 53}
    if name in expected and len(rows) != expected[name]:
        raise ValueError(f"Unexpected count for {name}: {len(rows)} != {expected[name]}")


def summarize_local_dominance(summary: list[dict[str, object]]) -> dict[str, str]:
    out = {}
    for lg in LGS:
        inside = next(r for r in summary if r["lg"] == lg and r["region_class"] == "inside")
        vals = {"arrangement-dominant": float(inside["mean_q_arrangement"]), "baseline-dominant": float(inside["mean_q_baseline"]), "mixed": float(inside["mean_q_third"])}
        best = max(vals, key=vals.get)
        if vals[best] < 0.5:
            best = "mixed"
        out[lg] = best
    return out


def baseline_branch_ids(baseline_newick: str) -> dict[str, str]:
    return {sp: f"B{i:02d}" for i, sp in enumerate(sorted(splits(parse_tree(baseline_newick))), start=1)}


def mean_support_for_baseline(scored_newick: str, branch_ids: dict[str, str]) -> tuple[float, float, float]:
    supp = annotated_branch_support(scored_newick)
    q1s = [supp.get(split, {}).get("q1", math.nan) for split in branch_ids]
    q1s = [x for x in q1s if not math.isnan(x)]
    if not q1s:
        return math.nan, math.nan, math.nan
    return statistics.fmean(q1s), statistics.median(q1s), min(q1s)


def downweight_rows(windows: list[dict[str, object]], m, rep: int) -> tuple[list[dict[str, object]], int]:
    outside = [r for r in windows if not as_bool(r["inside_inversion"]) and not as_bool(r["overlaps_inversion"])]
    if m == "ALL":
        return outside + [r for r in windows if as_bool(r["inside_inversion"])], 0
    seed = SEED_BASE + 1000 * int(m) + rep
    rng = random.Random(seed)
    chosen = list(outside)
    for lg in LGS:
        inside = [r for r in windows if r["lg"] == lg and as_bool(r["inside_inversion"])]
        k = min(int(m), len(inside))
        chosen.extend(rng.sample(inside, k))
    return chosen, seed


def render_tree_label(newick: str, ax, title: str) -> None:
    tree = parse_tree(newick)
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("["):
            clade.name = ""
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(""); ax.set_ylabel("")


def plot_astral_figures(treatment_rows, down_rows, branch_rows, inferred_newicks, baseline_newick):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    render_tree_label(baseline_newick, axes[0, 0], "A. frozen collinear baseline")
    render_tree_label(inferred_newicks["T0_all"], axes[0, 1], "B. ASTRAL4 all windows")
    render_tree_label(inferred_newicks["T1_outside"], axes[1, 0], "C. ASTRAL4 outside-only")
    ax = axes[1, 1]
    finite = [r for r in down_rows if r["m"] != "ALL"]
    positions = []
    labels = []
    data = []
    for i, m in enumerate([1, 2, 5, 10, 20, 40]):
        vals = [float(r["fraction_baseline_splits_recovered"]) for r in finite if int(r["m"]) == m]
        if vals:
            positions.append(i); labels.append(str(m)); data.append(vals)
    allv = [float(r["fraction_baseline_splits_recovered"]) for r in down_rows if r["m"] == "ALL"]
    if allv:
        positions.append(len(positions)); labels.append("ALL"); data.append(allv)
    ax.boxplot(data, positions=positions, patch_artist=True, boxprops={"facecolor":"#c6dbef"}, medianprops={"color":"black"})
    ax.set_xticks(positions); ax.set_xticklabels(labels); ax.set_ylim(0,1.05)
    ax.set_xlabel("max inversion windows per LG (m)"); ax.set_ylabel("baseline splits recovered")
    ax.set_title("D. inversion-window downweighting", fontsize=9)
    fig.tight_layout(); fig.savefig(ASTRAL_FIG_PDF); fig.savefig(ASTRAL_FIG_PNG, dpi=300); plt.close(fig)

    # branch support figure fixed baseline branches
    selected = [r for r in branch_rows if r["treatment"] in {"T0_all", "T1_outside", "T2_inside", "downweight_m1_mean"}]
    branches = sorted({r["branch_id"] for r in selected})
    treatments = ["T0_all", "T1_outside", "T2_inside", "downweight_m1_mean"]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    width = 0.18
    colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    for j, tr in enumerate(treatments):
        vals = []
        for b in branches:
            row = next((r for r in selected if r["treatment"] == tr and r["branch_id"] == b), None)
            vals.append(float(row["q1"]) if row else math.nan)
        ax.bar([i + (j-1.5)*width for i in range(len(branches))], vals, width=width, color=colors[j], label=tr)
    ax.set_xticks(range(len(branches))); ax.set_xticklabels(branches, rotation=45, ha="right")
    ax.set_ylabel("ASTRAL4 q1 for fixed baseline resolution"); ax.set_ylim(0,1.05); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(BRANCH_FIG_PDF); fig.savefig(BRANCH_FIG_PNG, dpi=300); plt.close(fig)


def main(run_astral: bool = True) -> None:
    frozen = verify_frozen()
    astral = find_astral4()
    astral_help = subprocess.run([str(astral)], capture_output=True, text=True).stdout + subprocess.run([str(astral)], capture_output=True, text=True).stderr
    aster_root = astral.parents[1]
    commit = "unknown"
    status = "not a git checkout"
    if (aster_root / ".git").exists():
        commit = subprocess.check_output(["git", "-C", str(aster_root), "rev-parse", "HEAD"], text=True).strip()
        status = subprocess.check_output(["git", "-C", str(aster_root), "status", "--short"], text=True).strip() or "clean"

    local_rows, local_summary = local_support_table()
    write_tsv(LOCAL_SUPPORT, local_rows, ["lg","window_id","start","end","midpoint","inside_inversion","overlaps_inversion","n_informative_quartets","n_resolved","q_arrangement","q_baseline","q_third","D_arrangement_minus_baseline"])
    write_tsv(LOCAL_SUPPORT_SUMMARY, local_summary, ["lg","region_class","mean_q_arrangement","mean_q_baseline","mean_q_third","median_q_arrangement","median_q_baseline","median_q_third"])
    plot_local_support(local_rows, regions())

    windows = load_windows()
    expected_taxa = taxa(parse_tree(str(windows[0]["tree_newick"])))
    baseline_newick = prune_baseline_12(expected_taxa)
    branch_ids = baseline_branch_ids(baseline_newick)
    treatments = create_treatments(windows)
    # Rename per-LG keys to exact requested names in outputs
    renamed = {}
    for k, v in treatments.items():
        if k.startswith("LG") and k.endswith("_inside"):
            lg = "LG" + k[2:4]
            renamed[f"{lg}_inside"] = v
        else:
            renamed[k] = v
    treatments = renamed
    for name, rows in treatments.items():
        treatment_count_expected(name, rows)
    treatment_files = write_treatment_files(treatments, expected_taxa)

    ASTRAL_DIR.mkdir(parents=True, exist_ok=True)
    commands = []
    treatment_summary = []
    branch_rows = []
    inferred_newicks = {}
    if run_astral:
        for name, input_path in treatment_files.items():
            safe = name
            inf_out = ASTRAL_DIR / f"{safe}.inferred.nwk"; inf_log = ASTRAL_DIR / f"{safe}.inferred.log"
            score_out = ASTRAL_DIR / f"{safe}.baseline_scored.nwk"; score_log = ASTRAL_DIR / f"{safe}.baseline_scored.log"
            cmd, rt, inf_newick = astral_infer(astral, input_path, inf_out, inf_log)
            commands.append({"treatment": name, "command_type": "infer", "command": cmd, "runtime_seconds": rt})
            scmd, srt, scored_newick = astral_score(astral, input_path, BASELINE_12, score_out, score_log)
            commands.append({"treatment": name, "command_type": "score_baseline", "command": scmd, "runtime_seconds": srt})
            inferred_newicks[name] = inf_newick
            rf = rf_to_baseline(inf_newick, baseline_newick)
            treatment_summary.append({"treatment": name, "n_input_trees": validate_treefile(input_path, expected_taxa), "ASTER_commit": commit, "command": cmd, "runtime_seconds": rt, **rf})
            supp = annotated_branch_support(scored_newick)
            for split, bid in branch_ids.items():
                ann = supp.get(split, {})
                branch_rows.append({"treatment": name, "branch_id": bid, "canonical_split": split, "q1": ann.get("q1", math.nan), "q2": ann.get("q2", math.nan), "q3": ann.get("q3", math.nan), "localPP1": ann.get("pp1", ann.get("localPP", math.nan)), "localPP2": ann.get("pp2", math.nan), "localPP3": ann.get("pp3", math.nan), "n_input_trees": validate_treefile(input_path, expected_taxa)})

        # Downweighting replicates
        down_rows = []
        m1_q1_by_branch = defaultdict(list)
        outside = [r for r in windows if not as_bool(r["inside_inversion"]) and not as_bool(r["overlaps_inversion"])]
        for m in M_VALUES:
            reps = [0] if m == "ALL" else list(range(1, FINITE_REPS + 1))
            for rep in reps:
                chosen, seed = downweight_rows(windows, m, rep)
                label = f"downweight_m{m}_rep{rep:03d}" if m != "ALL" else "downweight_mALL"
                path = INTERMEDIATE / "downweighting" / f"{label}.trees"
                write_treefile(path, chosen); validate_treefile(path, expected_taxa)
                inf_out = ASTRAL_DIR / "downweighting" / f"{label}.inferred.nwk"; inf_log = ASTRAL_DIR / "downweighting" / f"{label}.inferred.log"
                score_out = ASTRAL_DIR / "downweighting" / f"{label}.baseline_scored.nwk"; score_log = ASTRAL_DIR / "downweighting" / f"{label}.baseline_scored.log"
                cmd, rt, inf_newick = astral_infer(astral, path, inf_out, inf_log)
                scmd, srt, scored_newick = astral_score(astral, path, BASELINE_12, score_out, score_log)
                rf = rf_to_baseline(inf_newick, baseline_newick)
                meanq, medq, minq = mean_support_for_baseline(scored_newick, branch_ids)
                down_rows.append({"m": m, "replicate": rep, "seed": seed, "n_input_trees": len(chosen), **rf, "mean_baseline_q1": meanq, "median_baseline_q1": medq, "minimum_baseline_q1": minq})
                if m == 1:
                    supp = annotated_branch_support(scored_newick)
                    for split, bid in branch_ids.items():
                        if split in supp and "q1" in supp[split]:
                            m1_q1_by_branch[bid].append(supp[split]["q1"])
        write_tsv(DOWNWEIGHTING, down_rows, ["m","replicate","seed","n_input_trees","unrooted_RF","normalized_RF","n_baseline_internal_splits","n_recovered_baseline_internal_splits","fraction_baseline_splits_recovered","exact_topology_match","mean_baseline_q1","median_baseline_q1","minimum_baseline_q1"])
        # Add m=1 mean branch rows
        base_row_by_id = {r["branch_id"]: r for r in branch_rows if r["treatment"] == "T1_outside"}
        for split, bid in branch_ids.items():
            vals = m1_q1_by_branch.get(bid, [])
            branch_rows.append({"treatment": "downweight_m1_mean", "branch_id": bid, "canonical_split": split, "q1": statistics.fmean(vals) if vals else math.nan, "q2": math.nan, "q3": math.nan, "localPP1": math.nan, "localPP2": math.nan, "localPP3": math.nan, "n_input_trees": 249})
        # Split changes among baseline, T0, T1
        t0_splits = splits(parse_tree(inferred_newicks["T0_all"])); t1_splits = splits(parse_tree(inferred_newicks["T1_outside"])); base_splits = set(branch_ids)
        all_s = sorted(t0_splits | t1_splits | base_splits)
        t0_supp = annotated_branch_support(inferred_newicks["T0_all"]); t1_supp = annotated_branch_support(inferred_newicks["T1_outside"])
        split_rows = [{"canonical_split": sp, "present_in_baseline": sp in base_splits, "present_in_T0_all": sp in t0_splits, "present_in_T1_outside": sp in t1_splits, "support_T0": t0_supp.get(sp, {}).get("q1", math.nan), "support_T1": t1_supp.get(sp, {}).get("q1", math.nan)} for sp in all_s]
        write_tsv(SPLIT_CHANGES, split_rows, ["canonical_split","present_in_baseline","present_in_T0_all","present_in_T1_outside","support_T0","support_T1"])
        write_tsv(TREATMENT_SUMMARY, treatment_summary, ["treatment","n_input_trees","ASTER_commit","command","runtime_seconds","unrooted_RF","normalized_RF","n_baseline_internal_splits","n_recovered_baseline_internal_splits","fraction_baseline_splits_recovered","exact_topology_match"])
        write_tsv(BRANCH_SUPPORT, branch_rows, ["treatment","branch_id","canonical_split","q1","q2","q3","localPP1","localPP2","localPP3","n_input_trees"])
        plot_astral_figures(treatment_summary, down_rows, branch_rows, inferred_newicks, baseline_newick)
    else:
        treatment_summary = []
        down_rows = []
        branch_rows = []
        commands = []

    # provenance, report, methods/results
    help_excerpt = "\n".join(astral_help.splitlines()[:40])
    ASTER_PROVENANCE.write_text(textwrap.dedent(f"""
    # Stage-6 ASTER/ASTRAL4 provenance

    - ASTER repository: https://github.com/chaoszhang/ASTER
    - ASTRAL4 executable: `{astral}`
    - Executable SHA256: `{sha256(astral)}`
    - ASTER Git commit: `{commit}`
    - ASTER Git status: `{status}`
    - Platform: `{platform.platform()}`
    - Threads: `{THREADS}`
    - Compile command: not run in this repository; existing PATH executable used.
    - Date run: generated by Stage-6 script at execution time.

    ## ASTRAL4 help/version excerpt

    ```text
    {help_excerpt}
    ```
    """).strip() + "\n")

    write_stage6_texts(local_summary, treatment_summary, down_rows, commit, astral)
    outputs = [LOCAL_SUPPORT, LOCAL_SUPPORT_SUMMARY, QUARTET_SUPPORT_FIG_PDF, QUARTET_SUPPORT_FIG_PNG, BASELINE_12, ASTER_PROVENANCE, BRANCH_SUPPORT, TREATMENT_SUMMARY, DOWNWEIGHTING, SPLIT_CHANGES, ASTRAL_FIG_PDF, ASTRAL_FIG_PNG, BRANCH_FIG_PDF, BRANCH_FIG_PNG, REPORT, METHODS_TEXT, RESULTS_TEXT]
    manifest = {"stage": 6, "frozen_input_checksums": frozen, "ASTER_commit": commit, "ASTRAL4_executable": str(astral), "ASTRAL4_executable_sha256": sha256(astral), "threads": THREADS, "commands": commands, "treatment_definitions": {k: len(v) for k, v in treatments.items()}, "downweighting_m_values": M_VALUES, "downweighting_replicates_per_finite_m": FINITE_REPS, "seed_rule": "seed = 20260929 + 1000*m + replicate_id", "output_checksums": {str(p.relative_to(REPO_ROOT)): sha256(p) for p in outputs if p.exists()}}
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_stage6_texts(local_summary, treatment_summary, down_rows, commit, astral):
    dominance = summarize_local_dominance(local_summary)
    ts = {r["treatment"]: r for r in treatment_summary}
    def get(t, k, default="NA"):
        return ts.get(t, {}).get(k, default)
    if down_rows:
        m1 = [r for r in down_rows if str(r["m"]) == "1"]
        m1_exact = sum(str(r["exact_topology_match"]).lower() == "true" for r in m1) / len(m1)
        m1_medrf = statistics.median(float(r["normalized_RF"]) for r in m1)
        m1_range = (min(float(r["normalized_RF"]) for r in m1), max(float(r["normalized_RF"]) for r in m1))
        m1_rec = statistics.fmean(float(r["fraction_baseline_splits_recovered"]) for r in m1)
    else:
        m1_exact = m1_medrf = m1_rec = math.nan; m1_range=(math.nan, math.nan)
    REPORT.write_text(textwrap.dedent(f"""
    # Atlantic cod Stage-6 report

    ## A. Local quartet-support profiles

    Local support fractions are calculated within each published 250-kb population tree from the predeclared informative quartets. They are not ASTRAL4 branch q1/q2/q3 values. Inside-inversion local support is summarized as: LG01 {dominance.get('LG01')}, LG02 {dominance.get('LG02')}, LG07 {dominance.get('LG07')}, and LG12 {dominance.get('LG12')}.

    ## B. ASTRAL4 all-windows result

    The all-window ASTRAL4 population summary tree uses all 426 published 250-kb local population-tree windows. Its unrooted RF distance to the frozen 12-population baseline is {get('T0_all','unrooted_RF')} and it recovers {get('T0_all','n_recovered_baseline_internal_splits')}/{get('T0_all','n_baseline_internal_splits')} baseline internal splits.

    ## C. Collinear-only result

    The fully outside treatment uses 245 windows outside all frozen inversion intervals and excludes boundary-overlap windows. Its RF distance to the frozen baseline is {get('T1_outside','unrooted_RF')} and it recovers {get('T1_outside','n_recovered_baseline_internal_splits')}/{get('T1_outside','n_baseline_internal_splits')} baseline internal splits. This directly tests whether filtering inversion regions changes the inferred population-tree topology.

    ## D. Inversion-only result

    The inversion-only treatment uses 176 fully inside windows as a diagnostic linked-supergene tree set. Its RF distance to the frozen baseline is {get('T2_inside','unrooted_RF')} and it recovers {get('T2_inside','n_recovered_baseline_internal_splits')}/{get('T2_inside','n_baseline_internal_splits')} baseline internal splits.

    ## E. Leave-one-inversion-out

    Leave-one-inversion-out treatments are reported in `stage6_astral4_treatment_summary.tsv`. They identify which fully inside inversion-window set most changes topology or baseline-split recovery when removed.

    ## F. Downweighting

    Progressive inversion-window subsampling keeps all 245 fully outside windows and retains at most m fully inside windows per LG. The results in `stage6_astral4_downweighting.tsv` describe whether reducing effective inversion representation moves ASTRAL4 population-tree inference toward the collinear baseline topology.

    ## G. Block-balanced m=1

    For m=1, each inversion contributes at most one sampled inside window per replicate, plus all outside windows. Across 100 deterministic replicates, exact baseline matches occurred in {m1_exact:.3f} of replicates, median normalized RF was {m1_medrf:.3f} with range {m1_range[0]:.3f}-{m1_range[1]:.3f}, and mean baseline-split recovery was {m1_rec:.3f}.

    ## H. Interpretation

    Treating many linked inversion windows as separate input trees changes the quartet support available to a summary-tree estimator. The downweighting analysis tests whether reducing their effective representation shifts inference toward the collinear population-history signal. This empirical stress test does not show that ASTRAL4 is wrong, does not prove MSRC, and does not fit MSRC parameters.
    """).strip()+"\n")
    METHODS_TEXT.write_text(textwrap.dedent(f"""
    Stage 6 used ASTER/ASTRAL4 from `{astral}` at ASTER commit `{commit}`. ASTRAL4 population summary trees were inferred with `astral4 -R -u 2 -t {THREADS} -i INPUT -o OUTPUT`, where `-R` requests a more thorough search and `-u 2` reports detailed branch support annotations. The fixed frozen 12-population baseline tree was scored for each treatment with `astral4 -C -u 2 -t {THREADS} -c baseline_population_tree_12taxa.nwk -i INPUT -o OUTPUT`.

    Treatments were defined from the 426 published 250-kb local population-tree windows: all windows, fully outside windows, fully inside windows, leave-one-inversion-out sets, and per-LG inside-only sets. Boundary-overlap windows were excluded from clean inside/outside treatments. ASTRAL4 outputs are described as population summary trees because the taxa are Atlantic cod populations, not separate species.

    Downweighting kept all 245 fully outside windows and subsampled fully inside inversion windows independently within each LG with m in {{1,2,5,10,20,40,ALL}}. Finite m settings used 100 deterministic replicates with seed `20260929 + 1000*m + replicate_id`; `ALL` was deterministic. Inferred trees were compared with the frozen baseline using unrooted RF distance and fraction of recovered baseline splits. Fixed-baseline scoring was parsed from the observed ASTER `-u 2` annotation fields (`q1`, `q2`, `q3`, `pp1`, `pp2`, `pp3`, and `localPP`).
    """).strip()+"\n")
    RESULTS_TEXT.write_text(textwrap.dedent(f"""
    We first summarized local quartet-support profiles across the four Atlantic cod supergene linkage groups. These local support fractions, q_arrangement, q_baseline, and q_third, are calculated from predeclared informative quartets within each 250-kb local population tree and are distinct from ASTRAL4 branch q1/q2/q3 annotations. The support tracks show whether the dominant local quartet resolution changes within the frozen inversion intervals.

    We then used ASTER/ASTRAL4 as a population summary-tree sensitivity analysis. All published 250-kb windows, fully outside windows, fully inside windows, leave-one-inversion-out treatments, and per-LG inside-only treatments were analyzed as predeclared input sets. Each inferred population summary tree was compared with the independently frozen 12-population collinear baseline topology, and the same frozen baseline tree was also scored under every treatment using ASTRAL4 `-C -u 2`.

    Finally, we reduced the effective contribution of linked inversion windows by subsampling at most m inside windows per LG while keeping all fully outside windows. The m=1 treatment approximates each inversion region contributing one linked genealogy rather than many adjacent 250-kb windows. These results test whether replicated linked inversion windows materially influence quartet-summary population-tree inference without fitting MSRC parameters or rerunning phylogenetic inference.
    """).strip()+"\n")


class Tests(unittest.TestCase):
    def test_checksum_verification(self):
        self.assertIn("data/atlantic_cod/processed/cod_window_trees.tsv", verify_frozen())
    def test_extract_all_newicks_and_counts(self):
        self.assertEqual(len(load_windows()), 426)
        tr = create_treatments(load_windows())
        self.assertEqual(len(tr["T0_all"]), 426); self.assertEqual(len(tr["T1_outside"]), 245); self.assertEqual(len(tr["T2_inside"]), 176)
    def test_exact_taxa(self):
        labels = taxa(parse_tree(str(load_windows()[0]["tree_newick"])))
        self.assertEqual(len(labels), 12)
    def test_prune_baseline(self):
        labels = taxa(parse_tree(str(load_windows()[0]["tree_newick"])))
        nwk = prune_baseline_12(labels)
        self.assertEqual(taxa(parse_tree(nwk)), labels)
    def test_rf(self):
        a = "((A,B),(C,D));"; b = "((A,C),(B,D));"
        self.assertEqual(rf_to_baseline(a, a)["unrooted_RF"], 0)
        self.assertGreater(rf_to_baseline(b, a)["unrooted_RF"], 0)
    def test_canonical_split(self):
        self.assertEqual(canonical_split(["B","A"], ["A","B","C","D"]), "A,B|C,D")
    def test_local_normalization(self):
        rows, _ = local_support_table()
        for r in rows:
            if int(r["n_resolved"]):
                self.assertAlmostEqual(float(r["q_arrangement"])+float(r["q_baseline"])+float(r["q_third"]), 1.0, places=5)
    def test_downweighting_deterministic(self):
        windows = load_windows(); a, sa = downweight_rows(windows, 2, 7); b, sb = downweight_rows(windows, 2, 7)
        self.assertEqual(sa, sb); self.assertEqual([r["window_id"] for r in a], [r["window_id"] for r in b])
    def test_parser_saved_output(self):
        s = "(((C:0,D:0)'[q1=0.6;q2=0.3;q3=0.1;localPP=0.9;pp1=0.9;pp2=0.05;pp3=0.05]':0,B:0):0,A:0);"
        supp = annotated_branch_support(s)
        self.assertTrue(any(abs(v.get("q1",0)-0.6)<1e-9 for v in supp.values()))
    def test_branch_id_mapping(self):
        ids = baseline_branch_ids("(((A,B),C),(D,E));")
        self.assertTrue(ids)
    def test_deterministic_ordering(self):
        ids1 = baseline_branch_ids("(((A,B),C),(D,E));"); ids2 = baseline_branch_ids("(((A,B),C),(D,E));")
        self.assertEqual(ids1, ids2)


def tiny_integration(astral: Path) -> None:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp = td / "tiny.trees"; out = td / "out.nwk"; log = td / "out.log"
        inp.write_text("((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n")
        astral_infer(astral, inp, out, log)
        if not out.exists() or not out.read_text().strip():
            raise AssertionError("ASTRAL4 integration output missing")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tests", action="store_true")
    ap.add_argument("--run-integration-tests", action="store_true")
    ap.add_argument("--skip-astral", action="store_true", help="Generate local-support outputs only; for development")
    args = ap.parse_args()
    if args.run_tests:
        res = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        sys.exit(0 if res.wasSuccessful() else 1)
    if args.run_integration_tests:
        tiny_integration(find_astral4()); sys.exit(0)
    main(run_astral=not args.skip_astral)
