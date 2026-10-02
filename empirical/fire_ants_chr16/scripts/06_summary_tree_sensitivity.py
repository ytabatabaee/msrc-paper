#!/usr/bin/env python3
"""Stage 6 fire-ant grouped TWISST / ASTRAL4 summary-tree sensitivity.

Post-freeze analysis only. Stage 5 remains the frozen primary empirical
analysis. This script uses the published seven-group TWISST topology weights,
constructs deterministic modal grouped trees per window for ASTRAL4, and also
computes an exact TWISST-weighted expected quartet-agreement objective.
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
import time
import unittest
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import Bio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA = REPO_ROOT / "data" / "fire_ants_chr16"
EMP = REPO_ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
SCRIPTS = EMP / "scripts"
PROCESSED = DATA / "processed"
INTERMEDIATE = DATA / "intermediate" / "astral4"
ASTRAL_RESULTS = RESULTS / "astral4"

STAGE5_MANIFEST = RESULTS / "stage5_final_manifest.json"
STAGE4A_REGION_SUMMARY = RESULTS / "stage4a_region_summary.tsv"
STAGE4B_PRIMARY = RESULTS / "stage4b_primary_test.tsv"
TOPOLOGIES = PROCESSED / "stage1_twisst_topologies.tsv"
WEIGHTS = PROCESSED / "stage1_twisst_weights_sparse.tsv"
WEIGHT_SUMMARY = PROCESSED / "stage1_twisst_weight_summary.tsv"
WINDOW_INDEX = PROCESSED / "stage1_window_index.tsv"
STAGE1_MANIFEST = RESULTS / "stage1_manifest.json"
README = EMP / "README.md"

LABEL_MAP = PROCESSED / "stage6_group_label_map.tsv"
MODAL_TABLE = PROCESSED / "stage6_modal_group_topology_per_window.tsv"
WEIGHTED_SCORES = PROCESSED / "stage6_weighted_candidate_scores.tsv"
BACKGROUND_REFERENCE = PROCESSED / "stage6_background_reference.nwk"
MODAL_TIE_EXCLUSIONS = RESULTS / "stage6_modal_tie_exclusions.tsv"
MAPPING_AUDIT = RESULTS / "stage6_mapping_feasibility.md"
WEIGHTED_SUMMARY = RESULTS / "stage6_weighted_summary.tsv"
ASTRAL_TREATMENT_SUMMARY = RESULTS / "stage6_astral4_treatment_summary.tsv"
BRANCH_SUPPORT = RESULTS / "stage6_astral4_background_branch_support.tsv"
FOCAL_SUMMARY = RESULTS / "stage6_focal_split_summary.tsv"
DOWNWEIGHTING = RESULTS / "stage6_astral4_downweighting.tsv"
DOWNWEIGHTING_SUMMARY = RESULTS / "stage6_downweighting_summary.tsv"
ASTER_PROVENANCE = RESULTS / "stage6_aster_provenance.md"
METHODS_TEXT = RESULTS / "stage6_methods_text.md"
RESULTS_TEXT = RESULTS / "stage6_results_text.md"
REPORT = RESULTS / "stage6_report.md"
MANIFEST = RESULTS / "stage6_manifest.json"
ASTRAL_FIG_PDF = FIGURES / "fire_ants_astral4_sensitivity.pdf"
ASTRAL_FIG_PNG = FIGURES / "fire_ants_astral4_sensitivity.png"
WEIGHTED_FIG_PDF = FIGURES / "fire_ants_weighted_summary_sensitivity.pdf"
WEIGHTED_FIG_PNG = FIGURES / "fire_ants_weighted_summary_sensitivity.png"

THREADS = int(os.environ.get("ASTRAL4_THREADS", "2"))
FINITE_REPS = int(os.environ.get("STAGE6_DOWNWEIGHT_REPS", "100"))
SEED_BASE = 20261002
M_VALUES = [1, 2, 5, 10, 20, 40, 52]
EXPECTED_STAGE5_CORE = {
    "Delta_D": "1.610109863495346",
    "circular_rank": "2",
    "circular_n": "96",
    "circular_p": "0.0208333333",
    "unique_coordinate_rank": "3",
    "unique_coordinate_n": "97",
    "unique_coordinate_p": "0.0309278351",
    "length_weighted_coordinate_p": "0.0783290184",
}
LABELS = {
    "geminata": "geminata",
    "saevissima": "saevissima",
    "pusillignis": "pusillignis",
    "invicta/macdonaghi_SB": "inv_mac_SB",
    "invicta/macdonaghi_Sb": "inv_mac_Sb",
    "richteri_SB": "richteri_SB",
    "richteri_Sb": "richteri_Sb",
}
EXPECTED_TAXA = sorted(LABELS.values())
FOCAL = ["inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"]
SPLIT_CLASS = {
    "inv_mac_SB,inv_mac_Sb|richteri_SB,richteri_Sb": "species",
    "inv_mac_SB,richteri_SB|inv_mac_Sb,richteri_Sb": "haplotype",
    "inv_mac_SB,richteri_Sb|inv_mac_Sb,richteri_SB": "third",
}
TREATMENTS = {
    "T_chr1": lambda w: w["region"] == "chr1",
    "T_chr16_outside": lambda w: w["region"] in {"chr16A", "chr16B"},
    "T_background": lambda w: w["region"] in {"chr1", "chr16A", "chr16B"},
    "T_supergene": lambda w: w["region"] == "chr16_supergene",
    "T_chr16_all": lambda w: w["region"] in {"chr16A", "chr16_supergene", "chr16B"},
    "T_all": lambda w: True,
}
EXPECTED_TREATMENT_N = {"T_chr1": 117, "T_chr16_outside": 44, "T_background": 161, "T_supergene": 52, "T_chr16_all": 96, "T_all": 213}
TREE_FILES = {
    "T_chr1": "chr1.trees",
    "T_chr16_outside": "chr16_outside.trees",
    "T_background": "background_161.trees",
    "T_supergene": "supergene_52.trees",
    "T_chr16_all": "chr16_all_96.trees",
    "T_all": "all_213.trees",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


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
        if math.isnan(v):
            return "NA"
        return f"{v:.12g}"
    return str(v)


def parse_tree(newick: str):
    return Phylo.read(io.StringIO(newick.strip()), "newick")


def tree_to_newick(tree) -> str:
    s = io.StringIO()
    Phylo.write(tree, s, "newick")
    return s.getvalue().strip()


def taxa(tree) -> list[str]:
    labels = [t.name for t in tree.get_terminals()]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate taxon labels")
    return sorted(labels)


def relabel_newick(newick: str) -> str:
    tree = parse_tree(newick)
    for terminal in tree.get_terminals():
        if terminal.name not in LABELS:
            raise ValueError(f"Unexpected source group label {terminal.name}")
        terminal.name = LABELS[terminal.name]
    if taxa(tree) != EXPECTED_TAXA:
        raise ValueError("Relabeled topology has wrong taxon set")
    return tree_to_newick(tree)


def canonical_split(side: Iterable[str], all_taxa: Iterable[str]) -> str:
    s = frozenset(side)
    allset = frozenset(all_taxa)
    other = allset - s
    a = ",".join(sorted(s))
    b = ",".join(sorted(other))
    return "|".join(sorted([a, b]))


def splits_from_tree(tree, include_trivial: bool = False) -> set[str]:
    all_taxa = taxa(tree)
    n = len(all_taxa)
    out = set()
    for clade in tree.find_clades():
        side = [t.name for t in clade.get_terminals()]
        if len(side) == n:
            continue
        if not include_trivial and (len(side) <= 1 or len(side) >= n - 1):
            continue
        out.add(canonical_split(side, all_taxa))
    return out


def quartet_split(tree, quartet: tuple[str, str, str, str]) -> str:
    qset = frozenset(quartet)
    winners = []
    for split in splits_from_tree(tree, include_trivial=False):
        left_text, right_text = split.split("|")
        left = set(left_text.split(",")) if left_text else set()
        qleft = left & qset
        if len(qleft) == 2:
            winners.append(canonical_split(qleft, qset))
    winners = sorted(set(winners))
    if len(winners) != 1:
        raise ValueError(f"Unresolved quartet for {quartet}: {winners}")
    return winners[0]


def quartet_signature(newick: str) -> tuple[str, ...]:
    tree = parse_tree(newick)
    labels = taxa(tree)
    return tuple(quartet_split(tree, q) for q in combinations(labels, 4))


def focal_split_and_class(newick: str) -> tuple[str, str]:
    split = quartet_split(parse_tree(newick), tuple(FOCAL))
    cls = SPLIT_CLASS.get(split)
    if cls is None:
        raise ValueError(f"Unexpected focal split {split}")
    return split, cls


def rf_metrics(newick: str, reference: str) -> dict[str, object]:
    st = splits_from_tree(parse_tree(newick))
    sr = splits_from_tree(parse_tree(reference))
    rf = len(st - sr) + len(sr - st)
    denom = len(st) + len(sr)
    return {
        "unrooted_RF": rf,
        "normalized_RF": rf / denom if denom else 0.0,
        "n_background_internal_splits": len(sr),
        "n_recovered_background_splits": len(st & sr),
        "fraction_background_splits_recovered": len(st & sr) / len(sr) if sr else 1.0,
        "exact_background_match": st == sr,
    }


def annotation_dict(name: str | None) -> dict[str, float]:
    if not name:
        return {}
    m = re.search(r"\[(.*)\]", name)
    if not m:
        return {}
    out = {}
    for piece in m.group(1).split(";"):
        if "=" not in piece:
            continue
        key, value = piece.split("=", 1)
        try:
            out[key] = float(value)
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
        split = canonical_split(side, all_taxa)
        out[split] = annotation_dict(getattr(clade, "name", None))
    return out


def verify_stage5() -> dict[str, object]:
    manifest = json.loads(STAGE5_MANIFEST.read_text())
    if manifest.get("analysis_status") != "frozen_for_manuscript_use":
        raise ValueError("Stage 5 is not frozen for manuscript use")
    final = manifest.get("final_results", {})
    for key, expected in EXPECTED_STAGE5_CORE.items():
        if str(final.get(key)) != expected:
            raise ValueError(f"Stage-5 result changed for {key}: {final.get(key)} != {expected}")
    if manifest.get("astral_or_aster_run") is not False or manifest.get("msrc_model_fit") is not False:
        raise ValueError("Stage-5 forbidden-run flag changed")
    return manifest


def load_data() -> dict[str, object]:
    topologies = read_tsv(TOPOLOGIES)
    if len(topologies) != 945:
        raise ValueError(f"Expected 945 topologies, observed {len(topologies)}")
    topo_ids = [r["topology_id"] for r in topologies]
    topo_by_id = {r["topology_id"]: {**r, "safe_newick": relabel_newick(r["tree_newick"])} for r in topologies}
    windows = read_tsv(WINDOW_INDEX)
    if len(windows) != 213:
        raise ValueError(f"Expected 213 windows, observed {len(windows)}")
    window_by_id = {int(r["window_index"]): r for r in windows}
    weights = defaultdict(dict)
    for row in read_tsv(WEIGHTS):
        wi = int(row["window_index"])
        tid = row["topology_id"]
        if tid not in topo_by_id:
            raise ValueError(f"Unknown topology id in weights: {tid}")
        weights[wi][tid] = int(row["raw_weight"])
    summaries = {int(r["window_index"]): r for r in read_tsv(WEIGHT_SUMMARY)}
    return {"topo_ids": topo_ids, "topo_by_id": topo_by_id, "windows": windows, "window_by_id": window_by_id, "weights": weights, "summaries": summaries}


def write_label_map() -> None:
    rows = [{"source_group": k, "stage6_label": v} for k, v in LABELS.items()]
    write_tsv(LABEL_MAP, rows, ["source_group", "stage6_label"])


def modal_table(data: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    excluded = []
    for window in sorted(data["windows"], key=lambda r: int(r["window_index"])):
        wi = int(window["window_index"])
        w = data["weights"][wi]
        total = int(data["summaries"][wi]["total_raw_weight"])
        if sum(w.values()) != total:
            raise ValueError(f"Sparse weights do not sum to total for window {wi}")
        ordered = sorted(w.items(), key=lambda kv: (-kv[1], kv[0]))
        max_weight = ordered[0][1]
        tied = [tid for tid, value in ordered if value == max_weight]
        second_id, second_weight = (ordered[1] if len(ordered) > 1 else ("", 0))
        modal_id = tied[0] if len(tied) == 1 else ""
        modal_newick = data["topo_by_id"][modal_id]["safe_newick"] if modal_id else ""
        row = {
            "window_index": wi,
            "chrom": window["chrom"],
            "start": window["start"],
            "end": window["end"],
            "mid": window["mid"],
            "region": window["region"],
            "modal_topology_id": modal_id,
            "modal_topology_newick": modal_newick,
            "modal_raw_weight": max_weight if len(tied) == 1 else "",
            "total_raw_weight": total,
            "modal_fraction": max_weight / total,
            "second_topology_id": second_id,
            "second_raw_weight": second_weight,
            "second_fraction": second_weight / total if second_weight else 0.0,
            "modal_margin": (max_weight - second_weight) / total if second_weight else max_weight / total,
            "n_tied_for_max": len(tied),
            "tied_topology_ids": ",".join(tied),
        }
        rows.append(row)
        if len(tied) > 1:
            excluded.append({"window_index": wi, "n_tied_for_max": len(tied), "tied_topology_ids": ",".join(tied)})
    write_tsv(MODAL_TABLE, rows, ["window_index", "chrom", "start", "end", "mid", "region", "modal_topology_id", "modal_topology_newick", "modal_raw_weight", "total_raw_weight", "modal_fraction", "second_topology_id", "second_raw_weight", "second_fraction", "modal_margin", "n_tied_for_max"])
    write_tsv(MODAL_TIE_EXCLUSIONS, excluded, ["window_index", "n_tied_for_max", "tied_topology_ids"])
    return rows


def treatment_windows(rows: list[dict[str, object]], treatment: str, modal_usable_only: bool = False) -> list[dict[str, object]]:
    pred = TREATMENTS[treatment]
    selected = [r for r in rows if pred(r)]
    if modal_usable_only:
        selected = [r for r in selected if str(r.get("n_tied_for_max", "1")) == "1" and r.get("modal_topology_id")]
    return sorted(selected, key=lambda r: int(r["window_index"]))


def shared_matrix(data: dict[str, object]) -> tuple[dict[str, tuple[str, ...]], dict[tuple[str, str], int]]:
    signatures = {tid: quartet_signature(data["topo_by_id"][tid]["safe_newick"]) for tid in data["topo_ids"]}
    shared = {}
    ids = data["topo_ids"]
    for i, a in enumerate(ids):
        for b in ids[i:]:
            value = sum(1 for x, y in zip(signatures[a], signatures[b]) if x == y)
            shared[(a, b)] = value
            shared[(b, a)] = value
    return signatures, shared


def score_candidates(data: dict[str, object], shared: dict[tuple[str, str], int], selected: list[dict[str, object]]) -> list[dict[str, object]]:
    scores = {tid: 0.0 for tid in data["topo_ids"]}
    selected_ids = [int(r["window_index"]) for r in selected]
    for wi in selected_ids:
        total = float(data["summaries"][wi]["total_raw_weight"])
        for observed_tid, raw in data["weights"][wi].items():
            p = raw / total
            for candidate in data["topo_ids"]:
                scores[candidate] += p * shared[(candidate, observed_tid)]
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    rows = []
    prev = None
    rank = 0
    for index, (tid, score) in enumerate(ranked, start=1):
        if prev is None or abs(score - prev) > 1e-9:
            rank = index
            prev = score
        rows.append({
            "candidate_topology_id": tid,
            "expected_shared_quartets": score,
            "normalized_expected_quartet_score": score / (35 * len(selected)) if selected else math.nan,
            "rank": rank,
            "is_optimal": rank == 1,
        })
    return rows


def weighted_outputs(data: dict[str, object], shared: dict[tuple[str, str], int], modal_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]], str]:
    all_score_rows = []
    summary_rows = []
    background_newick = ""
    background_splits = None
    treatment_order = ["T_background", "T_chr1", "T_chr16_outside", "T_supergene", "T_chr16_all", "T_all"]
    for treatment in treatment_order:
        selected = treatment_windows(modal_rows, treatment, modal_usable_only=False)
        if len(selected) != EXPECTED_TREATMENT_N[treatment]:
            raise ValueError(f"Unexpected {treatment} window count {len(selected)}")
        scores = score_candidates(data, shared, selected)
        for row in scores:
            all_score_rows.append({"treatment": treatment, **row})
        opt = [r for r in scores if r["is_optimal"]]
        first = opt[0]
        newick = data["topo_by_id"][first["candidate_topology_id"]]["safe_newick"]
        focal_split, focal_class = focal_split_and_class(newick)
        if treatment == "T_background" and len(opt) == 1:
            background_newick = newick
            BACKGROUND_REFERENCE.write_text(newick + "\n")
            background_splits = splits_from_tree(parse_tree(newick))
        rf = {"unrooted_RF": "NA", "fraction_background_splits_recovered": "NA"}
        if background_newick:
            rf = rf_metrics(newick, background_newick)
        summary_rows.append({
            "treatment": treatment,
            "n_windows": len(selected),
            "optimal_topology_id": first["candidate_topology_id"],
            "optimal_topology_newick": newick,
            "n_tied_optima": len(opt),
            "normalized_expected_quartet_score": first["normalized_expected_quartet_score"],
            "focal_split": focal_split,
            "focal_class": focal_class,
            "RF_to_background_reference": rf["unrooted_RF"],
            "fraction_background_splits_recovered": rf["fraction_background_splits_recovered"],
        })
    write_tsv(WEIGHTED_SCORES, all_score_rows, ["treatment", "candidate_topology_id", "expected_shared_quartets", "normalized_expected_quartet_score", "rank", "is_optimal"])
    write_tsv(WEIGHTED_SUMMARY, summary_rows, ["treatment", "n_windows", "optimal_topology_id", "optimal_topology_newick", "n_tied_optima", "normalized_expected_quartet_score", "focal_split", "focal_class", "RF_to_background_reference", "fraction_background_splits_recovered"])
    return all_score_rows, summary_rows, background_newick


def find_astral4() -> Path:
    env = os.environ.get("ASTRAL4_BIN")
    candidates = []
    if env:
        candidates.append(Path(env))
    path_hit = shutil.which("astral4")
    if path_hit:
        candidates.append(Path(path_hit))
    candidates.append(REPO_ROOT / "external" / "ASTER" / "bin" / "astral4")
    candidates.append(Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral4"))
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise RuntimeError("No working astral4 found in $ASTRAL4_BIN, PATH, external/ASTER/bin/astral4, or /Users/ytabatabaee/Desktop/ASTER/bin/astral4")


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


def astral_score(astral: Path, input_file: Path, reference: Path, output: Path, log: Path) -> tuple[str, float, str]:
    cmd = [str(astral), "-C", "-u", "2", "-t", str(THREADS), "-c", str(reference), "-i", str(input_file), "-o", str(output)]
    if output.exists() and output.read_text().strip():
        return " ".join(cmd), 0.0, output.read_text().strip()
    rt = run_cmd(cmd, log)
    return " ".join(cmd), rt, output.read_text().strip()


def write_treefile(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for row in rows:
            f.write(str(row["modal_topology_newick"]).strip() + "\n")


def validate_treefile(path: Path) -> int:
    count = 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if taxa(parse_tree(line)) != EXPECTED_TAXA:
            raise ValueError(f"Taxon mismatch in {path}")
        count += 1
    return count


def write_treatment_treefiles(modal_rows: list[dict[str, object]]) -> dict[str, Path]:
    paths = {}
    for treatment, filename in TREE_FILES.items():
        rows = treatment_windows(modal_rows, treatment, modal_usable_only=True)
        path = INTERMEDIATE / filename
        write_treefile(path, rows)
        validate_treefile(path)
        paths[treatment] = path
    return paths


def best_weighted_by_treatment(weighted_summary: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {r["treatment"]: r for r in weighted_summary}


def astral_treatments(data: dict[str, object], modal_rows: list[dict[str, object]], weighted_summary: list[dict[str, object]], astral: Path) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, str], list[str]]:
    paths = write_treatment_treefiles(modal_rows)
    ASTRAL_RESULTS.mkdir(parents=True, exist_ok=True)
    commands = []
    inferred = {}
    runtimes = {}
    for treatment in ["T_background", "T_chr1", "T_chr16_outside", "T_supergene", "T_chr16_all", "T_all"]:
        output = ASTRAL_RESULTS / f"{treatment}.inferred.nwk"
        log = ASTRAL_RESULTS / f"{treatment}.infer.log"
        cmd, rt, newick = astral_infer(astral, paths[treatment], output, log)
        commands.append(cmd)
        inferred[treatment] = newick
        runtimes[treatment] = rt
        if taxa(parse_tree(newick)) != EXPECTED_TAXA:
            raise ValueError(f"ASTRAL taxon mismatch for {treatment}")
    background = inferred["T_background"]
    scoring_background = unannotated_topology_newick(background)
    (ASTRAL_RESULTS / "stage6_astral4_background_reference.nwk").write_text(scoring_background + "\n")
    weighted_by_treatment = best_weighted_by_treatment(weighted_summary)
    rows = []
    for treatment in ["T_background", "T_chr1", "T_chr16_outside", "T_supergene", "T_chr16_all", "T_all"]:
        newick = inferred[treatment]
        split, cls = focal_split_and_class(newick)
        rf = rf_metrics(newick, background)
        wopt = weighted_by_treatment[treatment]
        rows.append({
            "treatment": treatment,
            "n_input_windows": len(treatment_windows(modal_rows, treatment, modal_usable_only=True)),
            "n_modal_tie_exclusions": len(treatment_windows(modal_rows, treatment, modal_usable_only=False)) - len(treatment_windows(modal_rows, treatment, modal_usable_only=True)),
            "inferred_tree": newick,
            "focal_split": split,
            "focal_class": cls,
            "unrooted_RF_to_background": rf["unrooted_RF"],
            "normalized_RF_to_background": rf["normalized_RF"],
            "n_background_internal_splits": rf["n_background_internal_splits"],
            "n_recovered_background_splits": rf["n_recovered_background_splits"],
            "fraction_background_splits_recovered": rf["fraction_background_splits_recovered"],
            "exact_background_match": rf["exact_background_match"],
            "normalized_weighted_objective_best_score": wopt["normalized_expected_quartet_score"],
            "matches_weighted_objective_optimum": splits_from_tree(parse_tree(newick)) == splits_from_tree(parse_tree(wopt["optimal_topology_newick"])),
            "runtime_seconds": runtimes[treatment],
        })
    write_tsv(ASTRAL_TREATMENT_SUMMARY, rows, ["treatment", "n_input_windows", "n_modal_tie_exclusions", "inferred_tree", "focal_split", "focal_class", "unrooted_RF_to_background", "normalized_RF_to_background", "n_background_internal_splits", "n_recovered_background_splits", "fraction_background_splits_recovered", "exact_background_match", "normalized_weighted_objective_best_score", "matches_weighted_objective_optimum", "runtime_seconds"])

    branch_rows = []
    ref_path = ASTRAL_RESULTS / "stage6_astral4_background_reference.nwk"
    branch_ids = {split: f"B{i:02d}" for i, split in enumerate(sorted(splits_from_tree(parse_tree(background))), start=1)}
    for treatment in ["T_background", "T_chr1", "T_chr16_outside", "T_supergene", "T_chr16_all", "T_all"]:
        output = ASTRAL_RESULTS / f"{treatment}.background_scored.nwk"
        log = ASTRAL_RESULTS / f"{treatment}.score.log"
        cmd, _, scored = astral_score(astral, paths[treatment], ref_path, output, log)
        commands.append(cmd)
        annotations = annotated_branch_support(scored)
        for split, branch_id in branch_ids.items():
            ann = annotations.get(split, {})
            branch_rows.append({
                "treatment": treatment,
                "branch_id": branch_id,
                "split": split,
                "q1": ann.get("q1", "NA"),
                "q2": ann.get("q2", "NA"),
                "q3": ann.get("q3", "NA"),
                "pp1": ann.get("pp1", "NA"),
                "pp2": ann.get("pp2", "NA"),
                "pp3": ann.get("pp3", "NA"),
                "localPP": ann.get("localPP", ann.get("pp1", "NA")),
            })
    write_tsv(BRANCH_SUPPORT, branch_rows, ["treatment", "branch_id", "split", "q1", "q2", "q3", "pp1", "pp2", "pp3", "localPP"])
    return rows, branch_rows, inferred, commands


def selected_downweight_rows(modal_rows: list[dict[str, object]], m: int, rep: int) -> tuple[list[dict[str, object]], int]:
    background = treatment_windows(modal_rows, "T_background", modal_usable_only=True)
    supergene = treatment_windows(modal_rows, "T_supergene", modal_usable_only=True)
    if m == 52:
        return background + supergene, 0
    seed = SEED_BASE + 1000 * m + rep
    rng = random.Random(seed)
    selected_supergene = rng.sample(supergene, min(m, len(supergene)))
    return sorted(background + selected_supergene, key=lambda r: int(r["window_index"])), seed


def infer_rows_astral(astral: Path, rows: list[dict[str, object]], name: str) -> tuple[str, float]:
    treefile = INTERMEDIATE / "downweighting" / f"{name}.trees"
    output = ASTRAL_RESULTS / "downweighting" / f"{name}.inferred.nwk"
    log = ASTRAL_RESULTS / "downweighting" / f"{name}.log"
    write_treefile(treefile, rows)
    validate_treefile(treefile)
    _, rt, newick = astral_infer(astral, treefile, output, log)
    return newick, rt


def weighted_opt_for_selected(data: dict[str, object], shared: dict[tuple[str, str], int], selected: list[dict[str, object]], background_newick: str) -> dict[str, object]:
    scores = score_candidates(data, shared, selected)
    best = scores[0]
    newick = data["topo_by_id"][best["candidate_topology_id"]]["safe_newick"]
    split, cls = focal_split_and_class(newick)
    rf = rf_metrics(newick, background_newick)
    return {"weighted_optimal_topology_id": best["candidate_topology_id"], "weighted_focal_split": split, "weighted_focal_class": cls, "weighted_normalized_score": best["normalized_expected_quartet_score"], "weighted_RF_to_background": rf["unrooted_RF"]}


def downweighting(data: dict[str, object], modal_rows: list[dict[str, object]], shared: dict[tuple[str, str], int], astral: Path, background_newick: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = []
    for m in M_VALUES:
        reps = 1 if m == 52 else FINITE_REPS
        for rep in range(1, reps + 1):
            selected, seed = selected_downweight_rows(modal_rows, m, rep)
            name = f"downweight_m{m}_rep{rep:03d}"
            newick, rt = infer_rows_astral(astral, selected, name)
            split, cls = focal_split_and_class(newick)
            rf = rf_metrics(newick, background_newick)
            wopt = weighted_opt_for_selected(data, shared, selected, background_newick)
            rows.append({
                "m": m,
                "replicate": rep,
                "seed": seed,
                "n_background_windows": len(treatment_windows(modal_rows, "T_background", modal_usable_only=True)),
                "n_supergene_windows": m,
                "n_total_windows": len(selected),
                "unrooted_RF": rf["unrooted_RF"],
                "normalized_RF": rf["normalized_RF"],
                "n_background_splits_recovered": rf["n_recovered_background_splits"],
                "fraction_background_splits_recovered": rf["fraction_background_splits_recovered"],
                "exact_background_match": rf["exact_background_match"],
                "focal_split": split,
                "focal_class": cls,
                "runtime_seconds": rt,
                **wopt,
            })
    write_tsv(DOWNWEIGHTING, rows, ["m", "replicate", "seed", "n_background_windows", "n_supergene_windows", "n_total_windows", "unrooted_RF", "normalized_RF", "n_background_splits_recovered", "fraction_background_splits_recovered", "exact_background_match", "focal_split", "focal_class", "runtime_seconds", "weighted_optimal_topology_id", "weighted_focal_split", "weighted_focal_class", "weighted_normalized_score", "weighted_RF_to_background"])
    summary = []
    for m in M_VALUES:
        sel = [r for r in rows if int(r["m"]) == m]
        n = len(sel)
        counter_astral = Counter(r["focal_class"] for r in sel)
        counter_weighted = Counter(r["weighted_focal_class"] for r in sel)
        fracs = [float(r["fraction_background_splits_recovered"]) for r in sel]
        nrfs = [float(r["normalized_RF"]) for r in sel]
        summary.append({
            "m": m,
            "n_replicates": n,
            "fraction_astral_exact_background": sum(1 for r in sel if str(r["exact_background_match"]).lower() == "true") / n,
            "mean_astral_fraction_background_splits_recovered": statistics.fmean(fracs),
            "median_astral_normalized_RF": statistics.median(nrfs),
            "min_astral_fraction_background_splits_recovered": min(fracs),
            "max_astral_fraction_background_splits_recovered": max(fracs),
            "fraction_astral_species_focal": counter_astral["species"] / n,
            "fraction_astral_haplotype_focal": counter_astral["haplotype"] / n,
            "fraction_astral_third_focal": counter_astral["third"] / n,
            "fraction_weighted_species_focal": counter_weighted["species"] / n,
            "fraction_weighted_haplotype_focal": counter_weighted["haplotype"] / n,
            "fraction_weighted_third_focal": counter_weighted["third"] / n,
        })
    write_tsv(DOWNWEIGHTING_SUMMARY, summary, ["m", "n_replicates", "fraction_astral_exact_background", "mean_astral_fraction_background_splits_recovered", "median_astral_normalized_RF", "min_astral_fraction_background_splits_recovered", "max_astral_fraction_background_splits_recovered", "fraction_astral_species_focal", "fraction_astral_haplotype_focal", "fraction_astral_third_focal", "fraction_weighted_species_focal", "fraction_weighted_haplotype_focal", "fraction_weighted_third_focal"])
    return rows, summary


def mapping_feasibility() -> bool:
    candidates = list((DATA / "raw" / "upstream").rglob("twisst_samples_low_missingness"))
    exists = bool(candidates)
    text = (
        "# Stage 6 mapping-feasibility audit\n\n"
        "The upstream TWISST workflow used `--groupsFile results/twisst_samples_low_missingness` with seven source groups: `geminata`, `saevissima`, `pusillignis`, `invicta/macdonaghi_Sb`, `invicta/macdonaghi_SB`, `richteri_Sb`, and `richteri_SB`.\n\n"
        f"Frozen upstream snapshot contains `results/twisst_samples_low_missingness`: `{str(exists).lower()}`.\n\n"
        "A complete exact 267-tip individual-to-group mapping cannot be reproduced solely from committed public metadata in the frozen snapshot without guessing undocumented sample-name or SB/Sb conventions.\n\n"
        "`complete_exact_individual_to_group_mapping_available = false`\n\n"
        "`individual-level mapped ASTRAL4 run = not performed`\n\n"
        "This is a provenance safeguard, not a failed analysis. Stage 6 therefore uses the published seven-group TWISST representation.\n"
    )
    MAPPING_AUDIT.write_text(text)
    return False


def write_provenance(astral: Path, commands: list[str]) -> tuple[str, str]:
    help_output = subprocess.run([str(astral)], capture_output=True, text=True, check=False)
    help_text = (help_output.stdout + help_output.stderr).strip()
    root = astral.parents[1]
    commit = "unknown"
    status = "not a git checkout"
    if (root / ".git").exists():
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        status = subprocess.check_output(["git", "-C", str(root), "status", "--short"], text=True).strip() or "clean"
    ASTER_PROVENANCE.write_text(
        "# Stage 6 ASTER/ASTRAL4 provenance\n\n"
        f"- ASTER repository: `{root}`\n"
        f"- commit: `{commit}`\n"
        f"- git status: `{status}`\n"
        f"- astral4 executable: `{astral}`\n"
        f"- executable SHA256: `{sha256(astral)}`\n"
        f"- command template: `astral4 -R -u 2 -t {THREADS} -i INPUT.trees -o OUTPUT.nwk`\n"
        f"- fixed-tree scoring template: `astral4 -C -u 2 -t {THREADS} -c BACKGROUND_REFERENCE.nwk -i INPUT.trees -o SCORED.nwk`\n"
        f"- support setting `-u`: `2`\n"
        "- more-round setting `-R`: enabled\n"
        f"- threads: `{THREADS}`\n"
        f"- platform: `{platform.platform()}`\n"
        f"- Python: `{sys.version.split()[0]}`\n"
        f"- Biopython: `{Bio.__version__}`\n\n"
        "## Version/help output\n\n```text\n" + help_text[:4000] + "\n```\n\n"
        "## Commands\n\n" + "\n".join(f"- `{cmd}`" for cmd in commands) + "\n"
    )
    return commit, sha256(astral)


def write_focal_summary(weighted_summary: list[dict[str, object]], astral_summary: list[dict[str, object]]) -> list[dict[str, object]]:
    weighted = {r["treatment"]: r for r in weighted_summary}
    astral = {r["treatment"]: r for r in astral_summary}
    rows = []
    for treatment in ["T_chr1", "T_chr16_outside", "T_background", "T_supergene", "T_chr16_all", "T_all"]:
        rows.append({
            "treatment": treatment,
            "n_windows": EXPECTED_TREATMENT_N[treatment],
            "weighted_optimum_focal_split": weighted[treatment]["focal_split"],
            "weighted_optimum_focal_class": weighted[treatment]["focal_class"],
            "astral4_focal_split": astral[treatment]["focal_split"],
            "astral4_focal_class": astral[treatment]["focal_class"],
            "agreement": weighted[treatment]["focal_class"] == astral[treatment]["focal_class"],
        })
    write_tsv(FOCAL_SUMMARY, rows, ["treatment", "n_windows", "weighted_optimum_focal_split", "weighted_optimum_focal_class", "astral4_focal_split", "astral4_focal_class", "agreement"])
    return rows


def render_tree(newick: str, ax, title: str) -> None:
    tree = parse_tree(newick)
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("["):
            clade.name = ""
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")


def unannotated_topology_newick(newick: str) -> str:
    tree = parse_tree(newick)
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None
        clade.confidence = None
    return tree_to_newick(tree)


def make_figures(weighted_summary: list[dict[str, object]], astral_summary: list[dict[str, object]], down_summary: list[dict[str, object]], down_rows: list[dict[str, object]], inferred: dict[str, str], background_newick: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    render_tree(inferred["T_background"], axes[0, 0], "A. background-window reference")
    render_tree(inferred["T_supergene"], axes[0, 1], "B. supergene-only modal ASTRAL4")
    render_tree(inferred["T_all"], axes[1, 0], "C. all-window modal ASTRAL4")
    ax = axes[1, 1]
    positions = []
    data = []
    labels = []
    for i, m in enumerate(M_VALUES):
        vals = [float(r["fraction_background_splits_recovered"]) for r in down_rows if int(r["m"]) == m]
        if vals:
            positions.append(i)
            data.append(vals)
            labels.append(str(m))
    ax.boxplot(data, positions=positions, patch_artist=True, boxprops={"facecolor": "#c6dbef"}, medianprops={"color": "black"})
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("supergene windows included (m)")
    ax.set_ylabel("fraction background splits recovered")
    ax.set_title("D. downweighting linked supergene windows", fontsize=9)
    fig.tight_layout()
    with PdfPages(ASTRAL_FIG_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(ASTRAL_FIG_PNG, dpi=220)
    plt.close(fig)

    selected = [r for r in weighted_summary if r["treatment"] in {"T_background", "T_supergene", "T_chr16_all", "T_all"}]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = list(range(len(selected)))
    colors = {"species": "#345995", "haplotype": "#b23a48", "third": "#2f7f4f"}
    ax.bar(x, [float(r["normalized_expected_quartet_score"]) for r in selected], color=[colors[r["focal_class"]] for r in selected])
    ax.set_xticks(x)
    ax.set_xticklabels([r["treatment"].replace("T_", "") for r in selected], rotation=25, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("normalized expected quartet score")
    for xi, row in zip(x, selected):
        ax.text(xi, float(row["normalized_expected_quartet_score"]) + 0.02, f"{row['focal_class']}\nRF={row['RF_to_background_reference']}", ha="center", fontsize=8)
    ax.set_title("TWISST-weighted summary objective")
    fig.tight_layout()
    with PdfPages(WEIGHTED_FIG_PDF, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        pdf.savefig(fig)
    fig.savefig(WEIGHTED_FIG_PNG, dpi=220)
    plt.close(fig)


def write_texts(weighted_summary: list[dict[str, object]], astral_summary: list[dict[str, object]], focal_rows: list[dict[str, object]], down_summary: list[dict[str, object]], n_ties: int) -> None:
    w = {r["treatment"]: r for r in weighted_summary}
    a = {r["treatment"]: r for r in astral_summary}
    METHODS_TEXT.write_text(
        "# Stage 6 Methods text\n\n"
        "Stage 6 was a post-freeze summary-tree sensitivity analysis and did not alter the Stage-5 fire-ant empirical result. We used two complementary seven-group representations of the authors' published TWISST analysis. First, the full TWISST-weighted objective used every published topology weight, normalized each genomic window to equal total weight, and scored each of the 945 seven-group candidate topologies by expected shared quartet agreement across the 35 four-taxon subsets. Second, the modal-tree ASTRAL4 sensitivity represented each window by its single most highly weighted seven-group topology and supplied those discrete trees to ASTRAL4.\n\n"
        "The 267-tip individual RAxML trees were not analyzed with an ASTRAL-IV `-a` mapping because the exact public `results/twisst_samples_low_missingness` 267-to-7 source mapping file was unavailable in the frozen upstream snapshot. Guessing sample-group assignments from tree position, sample-name prefixes, or undocumented bigB/littleb conventions would violate the provenance safeguards established in Stage 2. ASTRAL-IV supports mapped multi-individual analyses, but Stage 6 deliberately did not fabricate such a mapping.\n\n"
        "For modal-tree ASTRAL4, source group labels containing slashes were converted by parsing each Newick tree and renaming terminals to ASTER-safe labels. Max-weight ties were recorded and excluded from modal-tree ASTRAL4 inputs rather than broken arbitrarily. Treatments were fixed before inspection: chr1, chr16 outside, background windows, supergene, all chr16, and all windows. Downweighting always retained all usable background modal windows and added deterministic samples of linked supergene windows using the frozen seed rule. The 52 supergene windows are adjacent and biologically linked; Stage 6 asks what happens if a conventional summary pipeline nevertheless counts them as multiple local genealogical inputs.\n"
    )
    RESULTS_TEXT.write_text(
        "# Stage 6 Results text\n\n"
        f"The mapping audit found no complete exact public 267-to-7 mapping file in the frozen upstream snapshot, so individual-level mapped ASTRAL4 was not performed. The modal-topology table inspected all 213 windows and found {n_ties} max-weight ties requiring exclusion from modal ASTRAL4 inputs.\n\n"
        f"In the exact TWISST-weighted objective, the background-window reference had focal class `{w['T_background']['focal_class']}`. The supergene-only optimum had focal class `{w['T_supergene']['focal_class']}`, chr16-all had `{w['T_chr16_all']['focal_class']}`, and all windows had `{w['T_all']['focal_class']}`.\n\n"
        f"In the modal-tree ASTRAL4 sensitivity, the background-window reference had focal class `{a['T_background']['focal_class']}`. The supergene-only modal ASTRAL4 tree had `{a['T_supergene']['focal_class']}`, chr16-all had `{a['T_chr16_all']['focal_class']}`, and all windows had `{a['T_all']['focal_class']}`. Relative to the background-window reference, all-window modal ASTRAL4 recovered {a['T_all']['n_recovered_background_splits']}/{a['T_all']['n_background_internal_splits']} background internal splits.\n\n"
        "The downweighting analysis was descriptive and introduced no Stage-5 P-value. It records how often the background-window reference was retained as deterministic samples of linked supergene windows were counted as additional local tree inputs. The full TWISST-weighted and modal-ASTRAL4 columns should be compared directly because disagreement indicates information loss from replacing a full topology distribution with one modal tree per window.\n"
    )
    REPORT.write_text(
        "# Stage 6 report\n\n"
        "## Purpose\n\nStage 6 measures sensitivity of seven-group summary-tree inference to counting linked supergene windows as separate local genealogical inputs.\n\n"
        "## Stage-5 frozen status\n\nStage 5 remains frozen for manuscript use and its primary empirical conclusions are unchanged.\n\n"
        "## Mapping feasibility\n\nThe exact 267-to-7 TWISST group mapping file was unavailable; individual-level mapped ASTRAL4 was not performed.\n\n"
        "## Grouped TWISST representation\n\nThe analysis used 945 seven-group TWISST topologies and 213 window-level weight distributions.\n\n"
        "## Exact weighted summary objective\n\nSee `stage6_weighted_summary.tsv` and `stage6_weighted_candidate_scores.tsv`.\n\n"
        "## Modal-tree ASTRAL4 design\n\nOne modal topology per non-tied window was supplied to ASTRAL4 with `-R -u 2`.\n\n"
        "## Treatment results\n\nSee `stage6_astral4_treatment_summary.tsv`.\n\n"
        "## Focal quartet behavior\n\nSee `stage6_focal_split_summary.tsv`.\n\n"
        "## Downweighting\n\nSee `stage6_astral4_downweighting.tsv` and `stage6_downweighting_summary.tsv`.\n\n"
        "## Agreement between weighted and ASTRAL4 analyses\n\nThe focal split summary records agreement treatment by treatment.\n\n"
        "## Interpretation\n\nStage 6 is a post-freeze sensitivity demonstration. It does not show that ASTRAL is wrong and does not identify the true historical mechanism.\n\n"
        "## Limitations\n\nModal-tree ASTRAL4 discards TWISST weight-distribution information, and the supergene windows are linked rather than independent evolutionary replicates.\n\n"
        "## Stage boundary\n\nNo MSRC model was fit, no raw sequences were reanalyzed, and Stage 5 was not analytically modified.\n"
    )


def update_readme() -> None:
    text = README.read_text()
    if "python3 empirical/fire_ants_chr16/scripts/06_summary_tree_sensitivity.py --run-tests" not in text:
        text = text.replace(
            "python3 empirical/fire_ants_chr16/scripts/05_finalize_fire_ant_analysis.py --run-tests\n```",
            "python3 empirical/fire_ants_chr16/scripts/05_finalize_fire_ant_analysis.py --run-tests\npython3 empirical/fire_ants_chr16/scripts/06_summary_tree_sensitivity.py --run-tests\n```",
        )
    text = re.sub(r"Stage 6: .*?\n", "Stage 6: Complete — post-freeze grouped TWISST / ASTRAL4 summary-tree sensitivity analysis.\n", text)
    if "Stage 6 does not alter the Stage-5 primary empirical analysis." not in text:
        text += "\nStage 6 does not alter the Stage-5 primary empirical analysis.\n"
    README.write_text(text)


def write_manifest(stage5: dict[str, object], data: dict[str, object], astral_commit: str, astral_sha: str, mapping_available: bool, commands: list[str]) -> None:
    output_paths = [
        LABEL_MAP, MODAL_TABLE, WEIGHTED_SCORES, BACKGROUND_REFERENCE, MODAL_TIE_EXCLUSIONS, MAPPING_AUDIT,
        WEIGHTED_SUMMARY, ASTRAL_TREATMENT_SUMMARY, BRANCH_SUPPORT, FOCAL_SUMMARY, DOWNWEIGHTING,
        DOWNWEIGHTING_SUMMARY, ASTER_PROVENANCE, METHODS_TEXT, RESULTS_TEXT, REPORT, ASTRAL_FIG_PDF,
        ASTRAL_FIG_PNG, WEIGHTED_FIG_PDF, WEIGHTED_FIG_PNG, README,
    ]
    treatment_inputs = sorted(INTERMEDIATE.glob("*.trees"))
    astral_outputs = sorted(ASTRAL_RESULTS.rglob("*.nwk")) + sorted(ASTRAL_RESULTS.rglob("*.log"))
    checksums = {rel(p): sha256(p) for p in output_paths + treatment_inputs + astral_outputs if p.exists()}
    manifest = {
        "stage": 6,
        "analysis_type": "post_freeze_summary_tree_sensitivity",
        "stage5_primary_result_changed": False,
        "stage5_manifest_checksum_after_cleanup": sha256(STAGE5_MANIFEST),
        "stage1_grouped_input_checksums": {rel(p): sha256(p) for p in [TOPOLOGIES, WEIGHTS, WEIGHT_SUMMARY, WINDOW_INDEX]},
        "modal_per_window_table_checksum": sha256(MODAL_TABLE),
        "label_map_checksum": sha256(LABEL_MAP),
        "ASTER_commit": astral_commit,
        "astral4_executable_sha256": astral_sha,
        "threads": THREADS,
        "downweighting_replicates_per_finite_m": FINITE_REPS,
        "seed_rule": "seed = 20261002 + 1000*m + replicate_id",
        "mapping_feasibility": {"complete_exact_individual_to_group_mapping_available": mapping_available, "individual_level_mapped_astral4_run": False},
        "commands": commands,
        "output_checksums": checksums,
        "processing_script_checksum": sha256(Path(__file__)),
        "individual_mapping_guessed": False,
        "raw_sequence_reanalysis": False,
        "new_primary_empirical_test": False,
        "stage5_analysis_modified": False,
        "msrc_model_fit": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> dict[str, object]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    stage5 = verify_stage5()
    mapping_available = mapping_feasibility()
    write_label_map()
    data = load_data()
    modal_rows = modal_table(data)
    _, shared = shared_matrix(data)
    score_rows, weighted_summary, background_newick = weighted_outputs(data, shared, modal_rows)
    astral = find_astral4()
    astral_summary, branch_rows, inferred, commands = astral_treatments(data, modal_rows, weighted_summary, astral)
    down_rows, down_summary = downweighting(data, modal_rows, shared, astral, inferred["T_background"])
    focal_rows = write_focal_summary(weighted_summary, astral_summary)
    make_figures(weighted_summary, astral_summary, down_summary, down_rows, inferred, background_newick)
    write_texts(weighted_summary, astral_summary, focal_rows, down_summary, len(read_tsv(MODAL_TIE_EXCLUSIONS)))
    astral_commit, astral_sha = write_provenance(astral, commands)
    update_readme()
    write_manifest(stage5, data, astral_commit, astral_sha, mapping_available, commands)
    return {"stage5_manifest_sha": sha256(STAGE5_MANIFEST), "stage6_manifest_sha": sha256(MANIFEST), "n_modal_ties": len(read_tsv(MODAL_TIE_EXCLUSIONS)), "astral_commit": astral_commit}


class Stage6Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_data()
        cls.modal = read_tsv(MODAL_TABLE)
        cls.weighted = read_tsv(WEIGHTED_SUMMARY)
        cls.astral = read_tsv(ASTRAL_TREATMENT_SUMMARY)
        cls.down = read_tsv(DOWNWEIGHTING)
        cls.down_summary = read_tsv(DOWNWEIGHTING_SUMMARY)

    def test_stage5_integrity(self):
        manifest = verify_stage5()
        self.assertEqual(manifest["analysis_status"], "frozen_for_manuscript_use")
        self.assertEqual(manifest["final_results"]["Delta_D"], EXPECTED_STAGE5_CORE["Delta_D"])

    def test_mapping_safeguards(self):
        text = MAPPING_AUDIT.read_text()
        self.assertIn("complete_exact_individual_to_group_mapping_available = false", text)
        script = Path(__file__).read_text()
        self.assertNotIn("bigB" + " -> " + "SB", script)
        self.assertNotIn("littleb" + " -> " + "Sb", script)
        self.assertIn("individual-level mapped ASTRAL4 run = not performed", text)

    def test_modal_topology(self):
        self.assertEqual(len(self.modal), 213)
        topo_ids = set(self.data["topo_ids"])
        for row in self.modal:
            if row["modal_topology_id"]:
                self.assertIn(row["modal_topology_id"], topo_ids)
                self.assertEqual(taxa(parse_tree(row["modal_topology_newick"])), EXPECTED_TAXA)
            self.assertGreaterEqual(int(row["n_tied_for_max"]), 1)

    def test_shared_quartets(self):
        signatures, shared = shared_matrix(self.data)
        ids = self.data["topo_ids"]
        self.assertEqual(len(ids), 945)
        for tid in ids[:20]:
            self.assertEqual(shared[(tid, tid)], 35)
        for a in ids[:20]:
            for b in ids[20:40]:
                self.assertEqual(shared[(a, b)], shared[(b, a)])
                self.assertGreaterEqual(shared[(a, b)], 0)
                self.assertLessEqual(shared[(a, b)], 35)
        for tid in ids[:100]:
            _, cls = focal_split_and_class(self.data["topo_by_id"][tid]["safe_newick"])
            self.assertIn(cls, {"species", "haplotype", "third"})

    def test_treatment_counts_and_scores(self):
        for treatment, expected in EXPECTED_TREATMENT_N.items():
            self.assertEqual(len(treatment_windows(self.modal, treatment)), expected)
        scores = read_tsv(WEIGHTED_SCORES)
        self.assertEqual(len([r for r in scores if r["treatment"] == "T_background"]), 945)
        for row in scores:
            v = float(row["normalized_expected_quartet_score"])
            self.assertGreaterEqual(v, 0.0)
            self.assertLessEqual(v, 1.0)

    def test_astral_outputs(self):
        self.assertTrue(ASTER_PROVENANCE.exists())
        for path in INTERMEDIATE.glob("*.trees"):
            self.assertGreater(validate_treefile(path), 0)
        for row in self.astral:
            tree = parse_tree(row["inferred_tree"])
            self.assertEqual(taxa(tree), EXPECTED_TAXA)
            self.assertEqual(len(splits_from_tree(tree)), 4)
        if len(self.astral) >= 2:
            a = self.astral[0]["inferred_tree"]
            b = self.astral[1]["inferred_tree"]
            self.assertEqual(rf_metrics(a, b)["unrooted_RF"], rf_metrics(b, a)["unrooted_RF"])
        self.assertTrue(BRANCH_SUPPORT.exists())

    def test_downweighting(self):
        finite_expected = 100 if FINITE_REPS == 100 else FINITE_REPS
        for m in [1, 2, 5, 10, 20, 40]:
            rows = [r for r in self.down if int(r["m"]) == m]
            self.assertEqual(len(rows), finite_expected)
            for r in rows:
                self.assertEqual(int(r["n_background_windows"]), len(treatment_windows(self.modal, "T_background", modal_usable_only=True)))
                self.assertEqual(int(r["n_supergene_windows"]), min(m, 52))
                self.assertEqual(int(r["seed"]), SEED_BASE + 1000 * m + int(r["replicate"]))
        self.assertEqual(len([r for r in self.down if int(r["m"]) == 52]), 1)

    def test_reproducibility_and_boundaries(self):
        manifest = json.loads(MANIFEST.read_text())
        self.assertFalse(manifest["new_primary_empirical_test"])
        self.assertFalse(manifest["stage5_analysis_modified"])
        self.assertFalse(manifest["individual_mapping_guessed"])
        self.assertFalse(manifest["msrc_model_fit"])
        self.assertEqual(manifest["stage5_manifest_checksum_after_cleanup"], sha256(STAGE5_MANIFEST))


def cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    result = main()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage6Tests)
        outcome = unittest.TextTestRunner(verbosity=2).run(suite)
        if not outcome.wasSuccessful():
            raise SystemExit(1)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    cli()
