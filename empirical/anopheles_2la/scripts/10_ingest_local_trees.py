#!/usr/bin/env python3
"""Stage 1B/2 synthetic-safe local-tree ingestion and quartet analysis.

The module keeps frozen structural predictions separate from observed local-tree
topologies. It contains no MalariaGEN retrieval code and must not be pointed at
real Anopheles topology inputs unless the Stage 1A completion gate exists.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import statistics
import sys
import unittest
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "anopheles_2la"
FIXTURE_ROOT = EMPIRICAL_ROOT / "tests" / "fixtures"
DATA_ROOT = REPO_ROOT / "data" / "anopheles_2la"
REAL_BLOCKED_PATTERNS = ("malariagen", "ag3", "haplot", "snp", "zarr", "vcf", "bcf", "trees", "topolog")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage1b_gate import NO_ELIGIBLE_STATUS, NoEligibleProspectiveDesigns, require_stage1b_designs

CHROM = "2L"
LEFT_BP = 20524058
RIGHT_BP = 42165532
PRIMARY_FLANK_BP = 5_000_000
PRIMARY_WINDOW_BP = 100_000
FINAL_PERMUTATIONS = 10_000
TEST_PERMUTATIONS = 199
SEED = 1729


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_split(a: str, b: str, c: str, d: str) -> str:
    left = ",".join(sorted([a, b]))
    right = ",".join(sorted([c, d]))
    sides = sorted([left, right])
    return f"{sides[0]}|{sides[1]}"


def split_sides(split: str) -> tuple[set[str], set[str]]:
    left, right = split.split("|")
    return set(left.split(",")), set(right.split(","))


def region_for_midpoint(midpoint: int) -> str:
    if LEFT_BP <= midpoint <= RIGHT_BP:
        return "inversion_2La"
    if LEFT_BP - PRIMARY_FLANK_BP <= midpoint < LEFT_BP:
        return "left_flank"
    if RIGHT_BP < midpoint <= RIGHT_BP + PRIMARY_FLANK_BP:
        return "right_flank"
    return "outside_primary_regions"


def make_windows(window_bp: int = PRIMARY_WINDOW_BP) -> list[dict[str, object]]:
    start = LEFT_BP - PRIMARY_FLANK_BP
    end = RIGHT_BP + PRIMARY_FLANK_BP
    rows = []
    pos = start
    idx = 0
    while pos <= end:
        wend = min(pos + window_bp - 1, end)
        midpoint = (pos + wend) // 2
        rows.append(
            {
                "window_id": f"win_{idx:04d}",
                "chromosome": CHROM,
                "start": pos,
                "end": wend,
                "midpoint": midpoint,
                "region_class": region_for_midpoint(midpoint),
                "distance_to_left_breakpoint": abs(midpoint - LEFT_BP),
                "distance_to_right_breakpoint": abs(midpoint - RIGHT_BP),
            }
        )
        pos = wend + 1
        idx += 1
    return rows


def guard_inputs(paths: Iterable[Path], synthetic_validation: bool, real_mode: bool) -> None:
    if synthetic_validation:
        fixture_root = FIXTURE_ROOT.resolve()
        for path in paths:
            resolved = path.resolve()
            if fixture_root not in [resolved, *resolved.parents]:
                raise RuntimeError(f"synthetic validation input must live under {fixture_root}: {path}")
            lowered = str(resolved).lower()
            if "synthetic" not in lowered and any(pattern in lowered for pattern in REAL_BLOCKED_PATTERNS):
                raise RuntimeError(f"synthetic validation rejects likely real-data path: {path}")
    if real_mode:
        require_stage1b_designs(DATA_ROOT)


def filter_quartets_for_enabled_designs(quartets: list[dict[str, str]], enabled_designs: set[str]) -> list[dict[str, str]]:
    prefixes = {"A": "Design_A", "B": "Design_B", "C": "Design_C"}
    enabled_tokens = tuple(prefixes[d] for d in sorted(enabled_designs))
    if not enabled_tokens:
        return []
    if len(enabled_designs) == 1 and quartets and "design_class" not in quartets[0]:
        return quartets
    return [q for q in quartets if q.get("design_class", "").startswith(enabled_tokens)]


@dataclass
class Node:
    name: str | None = None
    children: list["Node"] | None = None

    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(text: str) -> Node:
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    i = 0

    def parse_subtree() -> Node:
        nonlocal i
        if i >= len(text):
            raise ValueError("unexpected end of Newick")
        if text[i] == "(":
            i += 1
            children = [parse_subtree()]
            while i < len(text) and text[i] == ",":
                i += 1
                children.append(parse_subtree())
            if i >= len(text) or text[i] != ")":
                raise ValueError("unbalanced Newick")
            i += 1
            while i < len(text) and text[i] not in ",()":
                i += 1
            return Node(children=children)
        name = []
        while i < len(text) and text[i] not in ",()":
            name.append(text[i])
            i += 1
        label = "".join(name).split(":")[0].strip()
        if not label:
            raise ValueError("empty tip label")
        return Node(name=label)

    root = parse_subtree()
    if i != len(text):
        raise ValueError("trailing Newick content")
    return root


def collect_tips(node: Node) -> list[str]:
    if node.is_leaf():
        return [node.name or ""]
    tips: list[str] = []
    for child in node.children or []:
        tips.extend(collect_tips(child))
    return tips


def build_graph(node: Node) -> tuple[dict[str, set[str]], list[str]]:
    graph: dict[str, set[str]] = defaultdict(set)
    tips: list[str] = []
    counter = 0

    def walk(cur: Node, parent: str | None = None) -> str:
        nonlocal counter
        if cur.is_leaf():
            label = cur.name or ""
            tips.append(label)
            cur_id = label
        else:
            cur_id = f"_n{counter}"
            counter += 1
        if parent is not None:
            graph[cur_id].add(parent)
            graph[parent].add(cur_id)
        for child in cur.children or []:
            walk(child, cur_id)
        return cur_id

    walk(node)
    return graph, tips


def shortest_path_len(graph: dict[str, set[str]], src: str, dst: str) -> int:
    seen = {src}
    queue = [(src, 0)]
    for cur, dist in queue:
        if cur == dst:
            return dist
        for nxt in graph[cur]:
            if nxt not in seen:
                seen.add(nxt)
                queue.append((nxt, dist + 1))
    raise ValueError(f"no path between {src} and {dst}")


def classify_quartet_newick(newick: str, samples: list[str]) -> tuple[str, str, bool, str]:
    try:
        root = parse_newick(newick)
    except ValueError:
        return "", "unusable", False, "malformed_newick"
    graph, tips = build_graph(root)
    counts = Counter(tips)
    missing = [s for s in samples if counts[s] == 0]
    if missing:
        return "", "unusable", False, "missing_tip"
    duplicates = [s for s in samples if counts[s] > 1]
    if duplicates:
        return "", "unusable", False, "duplicate_tip"
    pair_dists: dict[tuple[str, str], int] = {}
    for i, a in enumerate(samples):
        for b in samples[i + 1 :]:
            pair_dists[tuple(sorted([a, b]))] = shortest_path_len(graph, a, b)
    candidate_splits = [
        ((samples[0], samples[1]), (samples[2], samples[3])),
        ((samples[0], samples[2]), (samples[1], samples[3])),
        ((samples[0], samples[3]), (samples[1], samples[2])),
    ]
    scores = []
    for left, right in candidate_splits:
        score = pair_dists[tuple(sorted(left))] + pair_dists[tuple(sorted(right))]
        scores.append((score, canonical_split(left[0], left[1], right[0], right[1])))
    scores.sort()
    if len(scores) > 1 and scores[0][0] == scores[1][0]:
        return "", "unresolved", False, "unresolved_tree"
    return scores[0][1], "resolved", True, ""


def topology_class(observed: str, arrangement: str, species: str) -> str:
    if not observed:
        return "unresolved"
    if observed == arrangement:
        return "arrangement"
    if observed == species:
        return "species"
    return "other"


def ingest_local_trees(windows: list[dict[str, str]], quartets: list[dict[str, str]], tree_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for window in windows:
        for quartet in quartets:
            tree_path = tree_dir / f"{window['window_id']}.{quartet['quartet_id']}.nwk"
            if not tree_path.exists():
                tree_path = tree_dir / f"{window['window_id']}.nwk"
            newick = tree_path.read_text().strip() if tree_path.exists() else ""
            samples = [quartet[f"sample_{i}"] for i in range(1, 5)]
            if not newick:
                observed, status, usable, reason = "", "unusable", False, "missing_tree"
            else:
                observed, status, usable, reason = classify_quartet_newick(newick, samples)
            arrangement = quartet["predicted_arrangement_split"]
            species = quartet["species_tree_split"]
            tclass = topology_class(observed, arrangement, species) if usable else status
            rows.append(
                {
                    "quartet_id": quartet["quartet_id"],
                    "window_id": window["window_id"],
                    "chromosome": window["chromosome"],
                    "start": window["start"],
                    "end": window["end"],
                    "midpoint": window["midpoint"],
                    "region": window["region_class"],
                    "observed_split": observed,
                    "topology_class": tclass,
                    "predicted_arrangement_split": arrangement,
                    "species_tree_split": species,
                    "matches_arrangement": usable and observed == arrangement,
                    "matches_species": usable and observed == species,
                    "usable": usable,
                    "exclusion_reason": reason,
                    "independence_group": quartet.get("independence_group", quartet["quartet_id"]),
                    "design_class": quartet.get("design_class", ""),
                }
            )
    return rows


def collapse_runs(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    runs: list[dict[str, object]] = []
    key_rows = sorted(
        [r for r in rows if r["usable"] is True],
        key=lambda r: (str(r["quartet_id"]), int(r["start"])),
    )
    current: dict[str, object] | None = None
    for row in key_rows:
        key = (row["quartet_id"], row["topology_class"], row["region"])
        if current and current["key"] == key and int(row["start"]) == int(current["last_end"]) + 1:
            current["end"] = row["end"]
            current["last_end"] = row["end"]
            current["n_windows"] = int(current["n_windows"]) + 1
        else:
            current = {
                "key": key,
                "quartet_id": row["quartet_id"],
                "region": row["region"],
                "topology_class": row["topology_class"],
                "start": row["start"],
                "end": row["end"],
                "last_end": row["end"],
                "n_windows": 1,
                "weight": 1.0,
            }
            runs.append(current)
    for run in runs:
        run.pop("key", None)
        run.pop("last_end", None)
    return runs


def summarize_inside_outside(rows: list[dict[str, object]], block_normalized: bool = False) -> list[dict[str, object]]:
    units = collapse_runs(rows) if block_normalized else [r for r in rows if r["usable"] is True]
    by_q: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in units:
        q = str(row["quartet_id"])
        region = str(row["region"])
        if region == "inversion_2La":
            bucket = "inside"
        elif region in {"left_flank", "right_flank"}:
            bucket = "flanks"
        else:
            continue
        weight = float(row.get("weight", 1.0))
        by_q[q][f"{bucket}_total"] += weight
        if row["topology_class"] == "arrangement":
            by_q[q][f"{bucket}_arrangement"] += weight
    summaries = []
    for q, vals in sorted(by_q.items()):
        inside_total = vals["inside_total"]
        flank_total = vals["flanks_total"]
        f_in = vals["inside_arrangement"] / inside_total if inside_total else math.nan
        f_out = vals["flanks_arrangement"] / flank_total if flank_total else math.nan
        summaries.append(
            {
                "quartet_id": q,
                "summary_mode": "block_normalized" if block_normalized else "window_weighted",
                "inside_total": inside_total,
                "inside_arrangement": vals["inside_arrangement"],
                "flank_total": flank_total,
                "flank_arrangement": vals["flanks_arrangement"],
                "f_in": f_in,
                "f_out": f_out,
                "delta_arr": f_in - f_out if not math.isnan(f_in) and not math.isnan(f_out) else math.nan,
                "enrichment": f_in / f_out if f_out and not math.isnan(f_out) else math.nan,
            }
        )
    return summaries


def mean_delta(rows: list[dict[str, object]]) -> float:
    vals = [float(r["delta_arr"]) for r in summarize_inside_outside(rows) if not math.isnan(float(r["delta_arr"]))]
    return statistics.mean(vals) if vals else math.nan


def circular_shift_test(rows: list[dict[str, object]], n_perm: int, seed: int = SEED) -> dict[str, object]:
    usable = [r for r in rows if r["usable"] is True and r["region"] in {"left_flank", "inversion_2La", "right_flank"}]
    observed = mean_delta(usable)
    rng = random.Random(seed)
    by_q: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sorted(usable, key=lambda r: (str(r["quartet_id"]), int(r["start"]))):
        by_q[str(row["quartet_id"])].append(row)
    perm_stats = []
    for _ in range(n_perm):
        shifted = []
        for qrows in by_q.values():
            regions = [r["region"] for r in qrows]
            k = rng.randrange(len(regions)) if regions else 0
            shifted_regions = regions[k:] + regions[:k]
            for row, region in zip(qrows, shifted_regions):
                nr = dict(row)
                nr["region"] = region
                shifted.append(nr)
        perm_stats.append(mean_delta(shifted))
    p = (1 + sum(1 for stat in perm_stats if stat >= observed)) / (n_perm + 1) if not math.isnan(observed) else math.nan
    return {"method": "block_aware_circular_shift", "n_permutations": n_perm, "observed_delta": observed, "p_value": p}


def naive_label_permutation(rows: list[dict[str, object]], n_perm: int, seed: int = SEED) -> dict[str, object]:
    usable = [r for r in rows if r["usable"] is True and r["region"] in {"left_flank", "inversion_2La", "right_flank"}]
    observed = mean_delta(usable)
    rng = random.Random(seed)
    regions = [r["region"] for r in usable]
    perm_stats = []
    for _ in range(n_perm):
        shuffled = regions[:]
        rng.shuffle(shuffled)
        perm_rows = []
        for row, region in zip(usable, shuffled):
            nr = dict(row)
            nr["region"] = region
            perm_rows.append(nr)
        perm_stats.append(mean_delta(perm_rows))
    p = (1 + sum(1 for stat in perm_stats if stat >= observed)) / (n_perm + 1) if not math.isnan(observed) else math.nan
    return {"method": "naive_per_window_label_permutation_diagnostic", "n_permutations": n_perm, "observed_delta": observed, "p_value": p}


def breakpoint_transitions(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    by_q: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["usable"] is True:
            by_q[str(row["quartet_id"])].append(row)
    for q, qrows in by_q.items():
        prev = None
        for row in sorted(qrows, key=lambda r: int(r["start"])):
            if prev and row["topology_class"] != prev["topology_class"]:
                pos = (int(prev["end"]) + int(row["start"])) // 2
                out.append(
                    {
                        "quartet_id": q,
                        "transition_position": pos,
                        "from_topology": prev["topology_class"],
                        "to_topology": row["topology_class"],
                        "distance_to_nearest_2La_breakpoint": min(abs(pos - LEFT_BP), abs(pos - RIGHT_BP)),
                    }
                )
            prev = row
    return out


def breakpoint_localization_test(rows: list[dict[str, object]], n_perm: int = TEST_PERMUTATIONS, seed: int = SEED) -> dict[str, object]:
    transitions = breakpoint_transitions(rows)
    if not transitions:
        return {"method": "transition_distance_circular_position_null", "n_transitions": 0, "observed_mean_distance": math.nan, "p_value": math.nan}
    observed = statistics.mean(float(t["distance_to_nearest_2La_breakpoint"]) for t in transitions)
    rng = random.Random(seed)
    mids = sorted({int(r["midpoint"]) for r in rows if r["usable"] is True})
    distances = []
    for _ in range(n_perm):
        sampled = [rng.choice(mids) for _ in transitions]
        distances.append(statistics.mean(min(abs(x - LEFT_BP), abs(x - RIGHT_BP)) for x in sampled))
    p = (1 + sum(1 for d in distances if d <= observed)) / (n_perm + 1)
    return {
        "method": "transition_distance_circular_position_null",
        "n_transitions": len(transitions),
        "observed_mean_distance": observed,
        "p_value": p,
    }


def design_b_contrasts(rows: list[dict[str, object]], pairs_path: Path) -> list[dict[str, object]]:
    if not pairs_path.exists():
        return []
    pairs = read_tsv(pairs_path)
    summaries = {r["quartet_id"]: r for r in summarize_inside_outside(rows)}
    out = []
    for pair in pairs:
        a = summaries.get(pair["quartet_a"])
        b = summaries.get(pair["quartet_b"])
        if not a or not b:
            continue
        da = float(a["delta_arr"])
        db = float(b["delta_arr"])
        direction = "changed_in_predicted_direction" if (da > 0 and db > 0 and pair["expected_change"] == "opposite_arrangement_prediction") else "no_predicted_change"
        if pair["expected_result"] == "null":
            direction = "null_like" if abs(da - db) <= 0.20 else "changed"
        out.append(
            {
                "contrast_id": pair["contrast_id"],
                "quartet_a": pair["quartet_a"],
                "quartet_b": pair["quartet_b"],
                "delta_a": da,
                "delta_b": db,
                "classification": direction,
                "synthetic_expected_result": pair["expected_result"],
            }
        )
    return out


def geography_summary(quartets: list[dict[str, str]]) -> list[dict[str, object]]:
    rows = []
    for q in quartets:
        countries = [q.get(f"country_{i}", "") for i in range(1, 5)]
        localities = [q.get(f"locality_{i}", "") for i in range(1, 5)]
        rows.append(
            {
                "quartet_id": q["quartet_id"],
                "all_country_present": all(countries),
                "country_matched": len(set(countries)) == 1 if all(countries) else False,
                "all_locality_present": all(localities),
                "locality_matched": len(set(localities)) == 1 if all(localities) else False,
            }
        )
    return rows


def generate_synthetic_fixtures(root: Path = FIXTURE_ROOT) -> None:
    root.mkdir(parents=True, exist_ok=True)
    species = ["arabiensis", "coluzzii", "gambiae", "melas", "merus", "quadriannulatus"]
    states = ["A0_homozygous", "A1_homozygous", "heterokaryotype", "unknown"]
    manifest = []
    prefixes = {
        "arabiensis": "A",
        "coluzzii": "C",
        "gambiae": "G",
        "melas": "ML",
        "merus": "MR",
        "quadriannulatus": "Q",
    }
    for si, sp in enumerate(species):
        for j in range(12):
            state = states[(j + si) % len(states)]
            if j in {0, 1, 4, 5, 8, 9}:
                state = "A0_homozygous" if j % 4 in {0, 1} else "A1_homozygous"
            manifest.append(
                {
                    "sample_id": f"SYN_{prefixes[sp]}_{j + 1:03d}",
                    "species": sp,
                    "taxon": sp,
                    "population": f"SYN_POP_{si % 3}",
                    "cohort": "SYN_COHORT",
                    "country": "SYN_COUNTRY_A" if si < 3 else "SYN_COUNTRY_B",
                    "locality": f"SYN_LOC_{j % 2}",
                    "2La_state": state,
                    "synthetic_fixture": "true",
                }
            )
    write_tsv(root / "synthetic_stage1a_manifest.tsv", manifest, list(manifest[0].keys()))
    quartets = [
        qrow("SYN_Q_A001", "Design_A_four_species_strict_2x2", ["SYN_A_001", "SYN_C_001", "SYN_G_005", "SYN_ML_005"], ["arabiensis", "coluzzii", "gambiae", "melas"], ["A0_homozygous", "A0_homozygous", "A1_homozygous", "A1_homozygous"], "SYN_A_001,SYN_C_001|SYN_G_005,SYN_ML_005", "SYN_A_001,SYN_G_005|SYN_C_001,SYN_ML_005", "grp_A1", "SYN_COUNTRY_A"),
        qrow("SYN_Q_A002", "Design_A_four_species_strict_2x2", ["SYN_A_002", "SYN_MR_009", "SYN_Q_005", "SYN_ML_006"], ["arabiensis", "merus", "quadriannulatus", "melas"], ["A0_homozygous", "A0_homozygous", "A1_homozygous", "A1_homozygous"], "SYN_A_002,SYN_MR_009|SYN_ML_006,SYN_Q_005", "SYN_A_002,SYN_Q_005|SYN_ML_006,SYN_MR_009", "grp_A2", "SYN_COUNTRY_B"),
        qrow("SYN_Q_B_POS_A", "Design_B_arrangement_replacement", ["SYN_G_001", "SYN_C_001", "SYN_A_005", "SYN_ML_005"], ["gambiae", "coluzzii", "arabiensis", "melas"], ["A0_homozygous", "A0_homozygous", "A1_homozygous", "A1_homozygous"], "SYN_C_001,SYN_G_001|SYN_A_005,SYN_ML_005", "SYN_A_005,SYN_G_001|SYN_C_001,SYN_ML_005", "grp_B_pos", "SYN_COUNTRY_A"),
        qrow("SYN_Q_B_POS_B", "Design_B_arrangement_replacement", ["SYN_G_005", "SYN_C_001", "SYN_A_005", "SYN_ML_001"], ["gambiae", "coluzzii", "arabiensis", "melas"], ["A1_homozygous", "A0_homozygous", "A1_homozygous", "A0_homozygous"], "SYN_C_001,SYN_ML_001|SYN_A_005,SYN_G_005", "SYN_A_005,SYN_C_001|SYN_G_005,SYN_ML_001", "grp_B_pos", "SYN_COUNTRY_A"),
        qrow("SYN_Q_B_NULL_A", "Design_B_arrangement_replacement", ["SYN_MR_009", "SYN_A_001", "SYN_C_005", "SYN_Q_005"], ["merus", "arabiensis", "coluzzii", "quadriannulatus"], ["A0_homozygous", "A0_homozygous", "A1_homozygous", "A1_homozygous"], "SYN_A_001,SYN_MR_009|SYN_C_005,SYN_Q_005", "SYN_C_005,SYN_MR_009|SYN_A_001,SYN_Q_005", "grp_B_null", "SYN_COUNTRY_B"),
        qrow("SYN_Q_B_NULL_B", "Design_B_arrangement_replacement", ["SYN_MR_005", "SYN_A_001", "SYN_C_005", "SYN_Q_001"], ["merus", "arabiensis", "coluzzii", "quadriannulatus"], ["A1_homozygous", "A0_homozygous", "A1_homozygous", "A0_homozygous"], "SYN_A_001,SYN_Q_001|SYN_C_005,SYN_MR_005", "SYN_A_001,SYN_C_005|SYN_MR_005,SYN_Q_001", "grp_B_null", "SYN_COUNTRY_B"),
        qrow("SYN_Q_C001", "Design_C_geography_population_matched", ["SYN_A_001", "SYN_C_001", "SYN_G_005", "SYN_ML_005"], ["arabiensis", "coluzzii", "gambiae", "melas"], ["A0_homozygous", "A0_homozygous", "A1_homozygous", "A1_homozygous"], "SYN_A_001,SYN_C_001|SYN_G_005,SYN_ML_005", "SYN_A_001,SYN_G_005|SYN_C_001,SYN_ML_005", "grp_C1", "SYN_COUNTRY_A"),
    ]
    write_tsv(root / "synthetic_frozen_strict_quartets_stage1a.tsv", quartets, list(quartets[0].keys()))
    write_tsv(root / "synthetic_genomic_windows.tsv", make_windows(), list(make_windows()[0].keys()))
    write_tsv(
        root / "synthetic_design_b_pairs.tsv",
        [
            {"contrast_id": "SYN_DESIGN_B_POSITIVE", "quartet_a": "SYN_Q_B_POS_A", "quartet_b": "SYN_Q_B_POS_B", "expected_change": "opposite_arrangement_prediction", "expected_result": "positive"},
            {"contrast_id": "SYN_DESIGN_B_NULL", "quartet_a": "SYN_Q_B_NULL_A", "quartet_b": "SYN_Q_B_NULL_B", "expected_change": "opposite_arrangement_prediction", "expected_result": "null"},
        ],
        ["contrast_id", "quartet_a", "quartet_b", "expected_change", "expected_result"],
    )
    for scenario in ["null_msc_like", "msrc_positive", "broad_alternative"]:
        generate_tree_scenario(root / "synthetic_gene_trees" / scenario, quartets, make_windows(), scenario)
    fasta = root / "synthetic_fasta"
    fasta.mkdir(exist_ok=True)
    (fasta / "win_0000.fasta").write_text(">SYN_A_001\nACGTACGT\n>SYN_C_001\nACGTACGA\n>SYN_G_005\nACGGACGA\n>SYN_ML_005\nACGGACGT\n")


def qrow(qid: str, design: str, samples: list[str], spp: list[str], states: list[str], arr: str, species_split: str, group: str, country: str) -> dict[str, str]:
    arr_left, arr_right = arr.split("|")
    sp_left, sp_right = species_split.split("|")
    row = {
        "quartet_id": qid,
        "design_class": design,
        "predicted_arrangement_split": canonical_split(arr_left.split(",")[0], arr_left.split(",")[1], arr_right.split(",")[0], arr_right.split(",")[1]),
        "species_tree_split": canonical_split(sp_left.split(",")[0], sp_left.split(",")[1], sp_right.split(",")[0], sp_right.split(",")[1]),
        "independence_group": group,
        "synthetic_fixture": "true",
        "state_evidence_type": "synthetic_fixture",
        "geography_matched": "true" if country == "SYN_COUNTRY_A" else "false",
    }
    for i in range(4):
        row[f"sample_{i + 1}"] = samples[i]
        row[f"species_{i + 1}"] = spp[i]
        row[f"arrangement_state_{i + 1}"] = states[i]
        row[f"country_{i + 1}"] = country
        row[f"locality_{i + 1}"] = "SYN_LOC_MATCHED" if country == "SYN_COUNTRY_A" else f"SYN_LOC_{i % 2}"
    return row


def split_to_newick(split: str) -> str:
    left, right = split_sides(split)
    l = sorted(left)
    r = sorted(right)
    return f"(({l[0]},{l[1]}),({r[0]},{r[1]}));"


def alternate_split(samples: list[str], arrangement: str, species: str) -> str:
    candidates = [
        canonical_split(samples[0], samples[1], samples[2], samples[3]),
        canonical_split(samples[0], samples[2], samples[1], samples[3]),
        canonical_split(samples[0], samples[3], samples[1], samples[2]),
    ]
    for split in candidates:
        if split not in {arrangement, species}:
            return split
    return candidates[-1]


def generate_tree_scenario(outdir: Path, quartets: list[dict[str, str]], windows: list[dict[str, object]], scenario: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED + len(scenario))
    for w in windows:
        first_split = None
        for quartet in quartets:
            region = str(w["region_class"])
            x = rng.random()
            samples = [quartet[f"sample_{i}"] for i in range(1, 5)]
            species = quartet["species_tree_split"]
            arrangement = quartet["predicted_arrangement_split"]
            other = alternate_split(samples, arrangement, species)
            if "B_NULL" in quartet["quartet_id"]:
                split = species if x < 0.78 else arrangement
            elif scenario == "msrc_positive" and region == "inversion_2La":
                split = arrangement if x < 0.80 else species
            elif scenario == "msrc_positive":
                split = species if x < 0.78 else arrangement
            elif scenario == "broad_alternative":
                split = arrangement if x < 0.62 else species
            else:
                split = species if x < 0.58 else (arrangement if x < 0.78 else other)
            first_split = first_split or split
            (outdir / f"{w['window_id']}.{quartet['quartet_id']}.nwk").write_text(split_to_newick(split) + "\n")
        if first_split:
            (outdir / f"{w['window_id']}.nwk").write_text(split_to_newick(first_split) + "\n")


def make_figures(track: list[dict[str, object]], summary: list[dict[str, object]], design_b: list[dict[str, object]], outdir: Path) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc_mplconfig")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    paths = []
    fig, ax = plt.subplots(figsize=(9, 2.8))
    ax.axvspan(LEFT_BP, RIGHT_BP, color="#d9d9d9", label="2La")
    ax.hlines(1, LEFT_BP - PRIMARY_FLANK_BP, LEFT_BP - 1, color="#4c78a8", lw=12)
    ax.hlines(1, LEFT_BP, RIGHT_BP, color="#f58518", lw=12)
    ax.hlines(1, RIGHT_BP + 1, RIGHT_BP + PRIMARY_FLANK_BP, color="#54a24b", lw=12)
    ax.text((LEFT_BP + RIGHT_BP) / 2, 1.2, "SYNTHETIC VALIDATION\nprediction frozen before genealogy observation", ha="center")
    ax.set_yticks([])
    ax.set_xlabel("2L coordinate")
    ax.set_title("Figure A: prospective Anopheles design")
    p = outdir / "figure_A_prospective_design.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)

    usable = [r for r in track if r["usable"] is True and r["quartet_id"] == "SYN_Q_A001"]
    ymap = {"species": 0, "arrangement": 1, "other": 2}
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.axvspan(LEFT_BP, RIGHT_BP, color="#eeeeee")
    ax.scatter([int(r["midpoint"]) for r in usable], [ymap.get(str(r["topology_class"]), 3) for r in usable], s=8)
    ax.set_yticks([0, 1, 2], ["species", "arrangement", "other"])
    ax.set_xlabel("2L coordinate")
    ax.set_title("Figure B: spatial topology track - SYNTHETIC VALIDATION")
    p = outdir / "figure_B_spatial_topology_track.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = range(len(summary))
    ax.bar([i - 0.2 for i in x], [float(r["f_out"]) for r in summary], width=0.4, label="flanks")
    ax.bar([i + 0.2 for i in x], [float(r["f_in"]) for r in summary], width=0.4, label="inside")
    ax.set_xticks(list(x), [r["quartet_id"] for r in summary], rotation=45, ha="right")
    ax.set_ylabel("arrangement topology frequency")
    ax.set_title("Figure C: inside versus outside - SYNTHETIC VALIDATION")
    ax.legend()
    p = outdir / "figure_C_inside_outside_enrichment.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)

    if design_b:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.bar([r["contrast_id"] for r in design_b], [float(r["delta_b"]) - float(r["delta_a"]) for r in design_b])
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("paired delta difference")
        ax.set_title("Figure D: arrangement replacement - SYNTHETIC VALIDATION")
        p = outdir / "figure_D_arrangement_replacement.pdf"
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        paths.append(p)
    return paths


TRACK_FIELDS = [
    "quartet_id",
    "window_id",
    "chromosome",
    "start",
    "end",
    "midpoint",
    "region",
    "observed_split",
    "topology_class",
    "predicted_arrangement_split",
    "species_tree_split",
    "matches_arrangement",
    "matches_species",
    "usable",
    "exclusion_reason",
    "independence_group",
    "design_class",
]


class Stage1BTests(unittest.TestCase):
    def test_region_boundaries(self) -> None:
        self.assertEqual(region_for_midpoint(LEFT_BP - 1), "left_flank")
        self.assertEqual(region_for_midpoint(LEFT_BP), "inversion_2La")
        self.assertEqual(region_for_midpoint(RIGHT_BP), "inversion_2La")
        self.assertEqual(region_for_midpoint(RIGHT_BP + 1), "right_flank")

    def test_quartet_classification(self) -> None:
        split, status, usable, reason = classify_quartet_newick("((A,B),(C,D));", ["A", "B", "C", "D"])
        self.assertTrue(usable)
        self.assertEqual(status, "resolved")
        self.assertEqual(reason, "")
        self.assertEqual(split, "A,B|C,D")

    def test_missing_and_polytomy(self) -> None:
        self.assertEqual(classify_quartet_newick("((A,B),(C,E));", ["A", "B", "C", "D"])[3], "missing_tip")
        self.assertEqual(classify_quartet_newick("(A,B,C,D);", ["A", "B", "C", "D"])[3], "unresolved_tree")

    def test_predicted_not_recomputed(self) -> None:
        self.assertEqual(topology_class("A,B|C,D", "A,C|B,D", "A,B|C,D"), "species")

    def test_enabled_design_filter(self) -> None:
        quartets = [
            {"quartet_id": "a", "design_class": "Design_A_four_species_strict_2x2"},
            {"quartet_id": "b", "design_class": "Design_B_arrangement_replacement"},
            {"quartet_id": "c", "design_class": "Design_C_geography_population_matched"},
        ]
        self.assertEqual([q["quartet_id"] for q in filter_quartets_for_enabled_designs(quartets, {"B"})], ["b"])
        self.assertEqual([q["quartet_id"] for q in filter_quartets_for_enabled_designs([{"quartet_id": "b"}], {"B"})], ["b"])

    def test_block_collapse_and_summary(self) -> None:
        rows = [
            {"quartet_id": "q", "start": 1, "end": 2, "region": "inversion_2La", "topology_class": "arrangement", "usable": True},
            {"quartet_id": "q", "start": 3, "end": 4, "region": "inversion_2La", "topology_class": "arrangement", "usable": True},
            {"quartet_id": "q", "start": 5, "end": 6, "region": "left_flank", "topology_class": "species", "usable": True},
        ]
        self.assertEqual(len(collapse_runs(rows)), 2)
        self.assertEqual(summarize_inside_outside(rows)[0]["delta_arr"], 1.0)

    def test_permutation_reproducible_and_breakpoints(self) -> None:
        rows = [
            {"quartet_id": "q", "start": LEFT_BP - 2, "end": LEFT_BP - 1, "midpoint": LEFT_BP - 1, "region": "left_flank", "topology_class": "species", "usable": True},
            {"quartet_id": "q", "start": LEFT_BP, "end": LEFT_BP + 1, "midpoint": LEFT_BP, "region": "inversion_2La", "topology_class": "arrangement", "usable": True},
        ]
        self.assertEqual(circular_shift_test(rows, 10, 1), circular_shift_test(rows, 10, 1))
        self.assertLessEqual(breakpoint_transitions(rows)[0]["distance_to_nearest_2La_breakpoint"], 1)


def run_tests() -> int:
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage1BTests))
    return 0 if result.wasSuccessful() else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", type=Path)
    parser.add_argument("--quartets", type=Path)
    parser.add_argument("--tree-dir", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--synthetic-validation", action="store_true")
    parser.add_argument("--real-mode", action="store_true")
    parser.add_argument("--generate-fixtures", action="store_true")
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args(argv)
    if args.run_tests:
        return run_tests()
    if args.generate_fixtures:
        generate_synthetic_fixtures()
        return 0
    if not (args.windows and args.quartets and args.tree_dir and args.out):
        parser.error("--windows, --quartets, --tree-dir, and --out are required unless generating fixtures/tests")
    try:
        guard_inputs([args.windows, args.quartets, args.tree_dir], args.synthetic_validation, args.real_mode)
        quartets = read_tsv(args.quartets)
        if args.real_mode:
            gate = require_stage1b_designs(DATA_ROOT)
            quartets = filter_quartets_for_enabled_designs(quartets, gate.enabled_designs)
        rows = ingest_local_trees(read_tsv(args.windows), quartets, args.tree_dir)
    except NoEligibleProspectiveDesigns as exc:
        print(exc.result.status)
        print(exc.result.message)
        return 0
    write_tsv(args.out, rows, TRACK_FIELDS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
