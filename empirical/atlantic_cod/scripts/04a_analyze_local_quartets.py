#!/usr/bin/env python3
"""Stage 4A descriptive local-quartet analysis for Atlantic cod.

This is the first stage allowed to inspect the published 250-kb local window
population trees. It verifies frozen Stage-2 and Stage-3 inputs before reading
local trees, then compares local induced quartets to independently frozen
arrangement and baseline splits. It does not run significance tests.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import statistics
import sys
import tempfile
import textwrap
import unittest
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = REPO_ROOT / "data" / "atlantic_cod"
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "atlantic_cod"

STAGE2_CANDIDATES = DATA_ROOT / "processed" / "stage2_structural_quartet_candidates.tsv"
EXPECTED_STAGE2_SHA = "303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c"
STAGE3_QUARTETS = DATA_ROOT / "processed" / "stage3_baseline_quartets.tsv"
EXPECTED_STAGE3_SHA = "5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1"
LOCAL_WINDOWS = DATA_ROOT / "processed" / "cod_window_trees.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"

LOCAL_QUARTETS = DATA_ROOT / "processed" / "stage4a_local_quartets.tsv"
WINDOW_SCORES = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"

INSIDE_OUTSIDE_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4a_inside_outside_summary.tsv"
TOPOLOGY_COMPOSITION = EMPIRICAL_ROOT / "results" / "stage4a_topology_composition.tsv"
QUARTET_EFFECTS = EMPIRICAL_ROOT / "results" / "stage4a_quartet_effects.tsv"
BOUNDARY_DESCRIPTIVE = EMPIRICAL_ROOT / "results" / "stage4a_boundary_descriptive.tsv"
BOUNDARY_SENSITIVITY = EMPIRICAL_ROOT / "results" / "stage4a_boundary_sensitivity.tsv"
LG_EFFECT_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4a_LG_effect_summary.tsv"
BORNHOLM_CHECK = EMPIRICAL_ROOT / "results" / "stage4a_bornholm_check.md"
STAGE4A_MANIFEST = EMPIRICAL_ROOT / "results" / "stage4a_manifest.json"
STAGE4A_REPORT = EMPIRICAL_ROOT / "results" / "stage4a_report.md"

FIG_DIR = EMPIRICAL_ROOT / "figures"
LG_FIGURES = {
    "LG01": FIG_DIR / "stage4a_LG01_topology_track.pdf",
    "LG02": FIG_DIR / "stage4a_LG02_topology_track.pdf",
    "LG07": FIG_DIR / "stage4a_LG07_topology_track.pdf",
    "LG12": FIG_DIR / "stage4a_LG12_topology_track.pdf",
}
FOUR_LG_FIGURE = FIG_DIR / "stage4a_four_LG_summary.pdf"

LGS = ("LG01", "LG02", "LG07", "LG12")
LOCAL_CLASSES = ("arrangement", "baseline", "third", "unresolved")

LOCAL_QUARTET_COLUMNS = [
    "lg", "window_id", "start", "end", "midpoint", "inside_inversion", "overlaps_inversion",
    "quartet_id", "taxon1", "taxon2", "taxon3", "taxon4", "arrangement_split", "baseline_split",
    "informative_for_competing_predictions", "local_split", "local_class", "source_tree", "notes",
]

WINDOW_SCORE_COLUMNS = [
    "lg", "window_id", "start", "end", "midpoint", "inside_inversion", "overlaps_inversion",
    "n_informative_quartets", "n_resolved", "n_arrangement", "n_baseline", "n_third", "n_unresolved",
    "fraction_arrangement", "fraction_baseline", "fraction_third", "D_arrangement_minus_baseline",
]


class UnresolvedQuartetError(ValueError):
    """Raised when a quartet has no unique unrooted split."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum(path: Path, expected: str, label: str) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(f"{label} checksum mismatch: observed {observed}, expected {expected}")
    return observed


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key, "")) for key in fieldnames})


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6f}"
    return str(value)


def bool_from_tsv(value: str) -> bool:
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    raise ValueError(f"Expected boolean string, got {value!r}")


def parse_tree(newick: str):
    tree = Phylo.read(io.StringIO(newick), "newick")
    labels = [terminal.name for terminal in tree.get_terminals()]
    if any(label is None or label == "" for label in labels):
        raise ValueError("Tree contains unnamed terminal.")
    duplicates = [label for label, count in Counter(labels).items() if count > 1]
    if duplicates:
        raise ValueError(f"Tree contains duplicate taxon labels: {duplicates}")
    return tree


def canonical_split(side_a: Iterable[str], side_b: Iterable[str]) -> str:
    left = ",".join(sorted(side_a))
    right = ",".join(sorted(side_b))
    return "|".join(sorted([left, right]))


def canonicalize_existing_split(split: str) -> str:
    if split == "unresolved":
        return split
    parts = split.split("|")
    if len(parts) != 2:
        raise ValueError(f"Invalid split: {split}")
    return canonical_split(parts[0].split(","), parts[1].split(","))


def induced_quartet_split(tree, taxa: list[str], tolerance: float = 1e-10) -> str:
    if len(taxa) != 4:
        raise ValueError("Exactly four taxa are required.")
    if len(set(taxa)) != 4:
        raise ValueError(f"Duplicate taxa in quartet: {taxa}")
    labels = {terminal.name for terminal in tree.get_terminals()}
    missing = sorted(set(taxa) - labels)
    if missing:
        raise ValueError(f"Taxa missing from local tree: {missing}")
    a, b, c, d = taxa
    candidates = [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]
    scores: list[tuple[float, str]] = []
    for left, right in candidates:
        score = float(tree.distance(left[0], left[1])) + float(tree.distance(right[0], right[1]))
        scores.append((score, canonical_split(left, right)))
    min_score = min(score for score, _ in scores)
    best = sorted(split for score, split in scores if abs(score - min_score) <= tolerance)
    if len(best) != 1:
        raise UnresolvedQuartetError(f"Unresolved induced quartet for {taxa}: {scores}")
    return best[0]


def third_possible_split(taxa: list[str], arrangement_split: str, baseline_split: str) -> str:
    a, b, c, d = taxa
    possible = {
        canonical_split((a, b), (c, d)),
        canonical_split((a, c), (b, d)),
        canonical_split((a, d), (b, c)),
    }
    remaining = possible - {canonicalize_existing_split(arrangement_split), canonicalize_existing_split(baseline_split)}
    if len(remaining) != 1:
        raise ValueError(f"Expected exactly one third topology, got {remaining}")
    return next(iter(remaining))


def classify_local_split(local_split: str, arrangement_split: str, baseline_split: str) -> str:
    if local_split == "unresolved":
        return "unresolved"
    arrangement = canonicalize_existing_split(arrangement_split)
    baseline = canonicalize_existing_split(baseline_split)
    if local_split == arrangement:
        return "arrangement"
    if local_split == baseline:
        return "baseline"
    return "third"


def read_stage3_quartets() -> dict[str, list[dict[str, str]]]:
    rows = read_tsv(STAGE3_QUARTETS)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        row = dict(row)
        row["arrangement_split"] = canonicalize_existing_split(row["arrangement_split"])
        row["baseline_split"] = canonicalize_existing_split(row["baseline_split"])
        grouped[row["lg"]].append(row)
    for lg in grouped:
        grouped[lg].sort(key=lambda item: item["quartet_id"])
    return grouped


def read_windows() -> list[dict[str, str]]:
    rows = read_tsv(LOCAL_WINDOWS)
    required = {"lg", "window_id", "start", "end", "midpoint", "inside_inversion", "overlaps_inversion", "tree_newick"}
    if not required.issubset(rows[0].keys() if rows else set()):
        raise ValueError("cod_window_trees.tsv lacks required Stage-4A columns.")
    return sorted([row for row in rows if row["lg"] in LGS], key=lambda r: (r["lg"], int(r["start"]), r["window_id"]))


def read_regions() -> dict[str, dict[str, int | str]]:
    regions: dict[str, dict[str, int | str]] = {}
    for row in read_tsv(REGION_MANIFEST):
        if row["lg"] in LGS:
            regions[row["lg"]] = {
                "region_name": row["region_name"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "coordinate_system": row["coordinate_system"],
            }
    return regions


def build_local_quartets(windows: list[dict[str, str]], quartets_by_lg: dict[str, list[dict[str, str]]]) -> tuple[list[dict[str, object]], list[str]]:
    output: list[dict[str, object]] = []
    parse_issues: list[str] = []
    for window in windows:
        lg = window["lg"]
        try:
            tree = parse_tree(window["tree_newick"])
        except Exception as exc:
            parse_issues.append(f"{window['window_id']}: {exc}")
            continue
        for quartet in quartets_by_lg[lg]:
            taxa = [quartet[f"taxon{i}"] for i in range(1, 5)]
            try:
                local_split = induced_quartet_split(tree, taxa)
                local_class = classify_local_split(local_split, quartet["arrangement_split"], quartet["baseline_split"])
                notes = "Local split induced from published 250-kb population tree."
            except UnresolvedQuartetError:
                local_split = "unresolved"
                local_class = "unresolved"
                notes = "Local quartet unresolved in the published 250-kb tree."
            output.append(
                {
                    "lg": lg,
                    "window_id": window["window_id"],
                    "start": int(window["start"]),
                    "end": int(window["end"]),
                    "midpoint": int(window["midpoint"]),
                    "inside_inversion": window["inside_inversion"],
                    "overlaps_inversion": window["overlaps_inversion"],
                    "quartet_id": quartet["quartet_id"],
                    "taxon1": taxa[0],
                    "taxon2": taxa[1],
                    "taxon3": taxa[2],
                    "taxon4": taxa[3],
                    "arrangement_split": quartet["arrangement_split"],
                    "baseline_split": quartet["baseline_split"],
                    "informative_for_competing_predictions": quartet["informative_for_competing_predictions"],
                    "local_split": local_split,
                    "local_class": local_class,
                    "source_tree": window["source_file"],
                    "notes": notes,
                }
            )
    return output, parse_issues


def compute_window_scores(local_rows: list[dict[str, object]], windows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in local_rows:
        grouped[(str(row["lg"]), str(row["window_id"]))].append(row)
    window_lookup = {(row["lg"], row["window_id"]): row for row in windows}
    output: list[dict[str, object]] = []
    for key in sorted(grouped, key=lambda k: (k[0], int(window_lookup[k]["start"]), k[1])):
        window = window_lookup[key]
        informative = [row for row in grouped[key] if row["informative_for_competing_predictions"] == "true"]
        counts = Counter(row["local_class"] for row in informative)
        n_unresolved = counts["unresolved"]
        n_resolved = len(informative) - n_unresolved
        n_arrangement = counts["arrangement"]
        n_baseline = counts["baseline"]
        n_third = counts["third"]
        if n_resolved:
            f_arr = n_arrangement / n_resolved
            f_base = n_baseline / n_resolved
            f_third = n_third / n_resolved
            d_value = f_arr - f_base
        else:
            f_arr = f_base = f_third = d_value = math.nan
        output.append(
            {
                "lg": window["lg"],
                "window_id": window["window_id"],
                "start": int(window["start"]),
                "end": int(window["end"]),
                "midpoint": int(window["midpoint"]),
                "inside_inversion": window["inside_inversion"],
                "overlaps_inversion": window["overlaps_inversion"],
                "n_informative_quartets": len(informative),
                "n_resolved": n_resolved,
                "n_arrangement": n_arrangement,
                "n_baseline": n_baseline,
                "n_third": n_third,
                "n_unresolved": n_unresolved,
                "fraction_arrangement": f_arr,
                "fraction_baseline": f_base,
                "fraction_third": f_third,
                "D_arrangement_minus_baseline": d_value,
            }
        )
    return output


def region_class(row: dict[str, object], boundary_mode: str = "separate") -> str:
    inside = bool_from_tsv(str(row["inside_inversion"]))
    overlaps = bool_from_tsv(str(row["overlaps_inversion"]))
    if boundary_mode == "inside" and overlaps:
        return "inside"
    if boundary_mode == "outside" and overlaps and not inside:
        return "outside"
    if inside:
        return "inside"
    if overlaps:
        return "boundary-overlap"
    return "outside"


def finite(values: Iterable[float]) -> list[float]:
    return [value for value in values if not math.isnan(value)]


def iqr(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[2] - q[0]


def summarize_values(values: list[float]) -> dict[str, object]:
    vals = finite(values)
    if not vals:
        return {"n_windows": 0, "mean": math.nan, "median": math.nan, "iqr": math.nan, "min": math.nan, "max": math.nan}
    return {
        "n_windows": len(vals),
        "mean": statistics.fmean(vals),
        "median": statistics.median(vals),
        "iqr": iqr(vals),
        "min": min(vals),
        "max": max(vals),
    }


def build_inside_outside_summary(window_scores: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for lg in LGS:
        lg_rows = [row for row in window_scores if row["lg"] == lg]
        for cls in ("inside", "outside", "boundary-overlap"):
            vals = [float(row["D_arrangement_minus_baseline"]) for row in lg_rows if region_class(row) == cls]
            summary = summarize_values(vals)
            rows.append({"lg": lg, "region_class": cls, **{f"D_{k}": v for k, v in summary.items()}})
    return rows


def build_topology_composition(window_scores: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for lg in LGS:
        lg_rows = [row for row in window_scores if row["lg"] == lg]
        for cls in ("inside", "outside", "boundary-overlap"):
            selected = [row for row in lg_rows if region_class(row) == cls]
            rows.append(
                {
                    "lg": lg,
                    "region_class": cls,
                    "n_windows": len(selected),
                    "mean_fraction_arrangement": statistics.fmean([float(row["fraction_arrangement"]) for row in selected]) if selected else math.nan,
                    "mean_fraction_baseline": statistics.fmean([float(row["fraction_baseline"]) for row in selected]) if selected else math.nan,
                    "mean_fraction_third": statistics.fmean([float(row["fraction_third"]) for row in selected]) if selected else math.nan,
                }
            )
    return rows


def build_quartet_effects(local_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    informative = [row for row in local_rows if row["informative_for_competing_predictions"] == "true"]
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in informative:
        grouped[(str(row["lg"]), str(row["quartet_id"]))].append(row)
    output: list[dict[str, object]] = []
    for key in sorted(grouped):
        rows = grouped[key]
        inside_rows = [row for row in rows if region_class(row) == "inside"]
        outside_rows = [row for row in rows if region_class(row) == "outside"]
        def frac(selected: list[dict[str, object]], cls: str) -> float:
            resolved = [row for row in selected if row["local_class"] != "unresolved"]
            if not resolved:
                return math.nan
            return sum(row["local_class"] == cls for row in resolved) / len(resolved)
        def dmean(selected: list[dict[str, object]]) -> float:
            resolved = [row for row in selected if row["local_class"] != "unresolved"]
            if not resolved:
                return math.nan
            return statistics.fmean((1 if row["local_class"] == "arrangement" else 0) - (1 if row["local_class"] == "baseline" else 0) for row in resolved)
        arr_inside = frac(inside_rows, "arrangement")
        arr_outside = frac(outside_rows, "arrangement")
        base_inside = frac(inside_rows, "baseline")
        base_outside = frac(outside_rows, "baseline")
        d_inside = dmean(inside_rows)
        d_outside = dmean(outside_rows)
        first = rows[0]
        output.append(
            {
                "lg": key[0],
                "quartet_id": key[1],
                "arrangement_split": first["arrangement_split"],
                "baseline_split": first["baseline_split"],
                "n_inside": len(inside_rows),
                "n_outside": len(outside_rows),
                "arrangement_fraction_inside": arr_inside,
                "arrangement_fraction_outside": arr_outside,
                "baseline_fraction_inside": base_inside,
                "baseline_fraction_outside": base_outside,
                "delta_arrangement": arr_inside - arr_outside if not math.isnan(arr_inside) and not math.isnan(arr_outside) else math.nan,
                "delta_arrangement_minus_baseline": d_inside - d_outside if not math.isnan(d_inside) and not math.isnan(d_outside) else math.nan,
            }
        )
    return output


def build_lg_effect_summary(window_scores: list[dict[str, object]], quartet_effects: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for lg in LGS:
        inside = [row for row in window_scores if row["lg"] == lg and region_class(row) == "inside"]
        outside = [row for row in window_scores if row["lg"] == lg and region_class(row) == "outside"]
        inside_d = finite(float(row["D_arrangement_minus_baseline"]) for row in inside)
        outside_d = finite(float(row["D_arrangement_minus_baseline"]) for row in outside)
        qfx = [row for row in quartet_effects if row["lg"] == lg and not math.isnan(float(row["delta_arrangement_minus_baseline"]))]
        positive = sum(float(row["delta_arrangement_minus_baseline"]) > 0 for row in qfx)
        output.append(
            {
                "lg": lg,
                "n_informative_quartets": len(qfx),
                "n_inside_windows": len(inside_d),
                "n_outside_windows": len(outside_d),
                "mean_D_inside": statistics.fmean(inside_d) if inside_d else math.nan,
                "mean_D_outside": statistics.fmean(outside_d) if outside_d else math.nan,
                "delta_D": (statistics.fmean(inside_d) - statistics.fmean(outside_d)) if inside_d and outside_d else math.nan,
                "median_D_inside": statistics.median(inside_d) if inside_d else math.nan,
                "median_D_outside": statistics.median(outside_d) if outside_d else math.nan,
                "fraction_quartets_positive_effect": positive / len(qfx) if qfx else math.nan,
            }
        )
    return output


def build_boundary_descriptive(window_scores: list[dict[str, object]], regions: dict[str, dict[str, int | str]], flank_n: int = 3) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for lg in LGS:
        lg_rows = [row for row in window_scores if row["lg"] == lg]
        for boundary_name, pos in (("left", int(regions[lg]["start"])), ("right", int(regions[lg]["end"]))):
            near = sorted(lg_rows, key=lambda row: abs(int(row["midpoint"]) - pos))[:flank_n * 2]
            near = sorted(near, key=lambda row: int(row["start"]))
            for row in near:
                output.append(
                    {
                        "lg": lg,
                        "boundary": boundary_name,
                        "boundary_position": pos,
                        "window_id": row["window_id"],
                        "start": row["start"],
                        "end": row["end"],
                        "midpoint": row["midpoint"],
                        "inside_inversion": row["inside_inversion"],
                        "overlaps_inversion": row["overlaps_inversion"],
                        "distance_midpoint_to_boundary": int(row["midpoint"]) - pos,
                        "D_arrangement_minus_baseline": row["D_arrangement_minus_baseline"],
                    }
                )
    return output


def build_boundary_sensitivity(window_scores: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    modes = {"boundary_excluded": "separate", "boundary_as_inside": "inside", "boundary_as_outside": "outside"}
    for lg in LGS:
        lg_rows = [row for row in window_scores if row["lg"] == lg]
        for label, mode in modes.items():
            inside = [float(row["D_arrangement_minus_baseline"]) for row in lg_rows if region_class(row, mode) == "inside"]
            outside = [float(row["D_arrangement_minus_baseline"]) for row in lg_rows if region_class(row, mode) == "outside"]
            inside = finite(inside)
            outside = finite(outside)
            output.append(
                {
                    "lg": lg,
                    "boundary_treatment": label,
                    "n_inside_windows": len(inside),
                    "n_outside_windows": len(outside),
                    "mean_D_inside": statistics.fmean(inside) if inside else math.nan,
                    "mean_D_outside": statistics.fmean(outside) if outside else math.nan,
                    "delta_D": (statistics.fmean(inside) - statistics.fmean(outside)) if inside and outside else math.nan,
                    "median_D_inside": statistics.median(inside) if inside else math.nan,
                    "median_D_outside": statistics.median(outside) if outside else math.nan,
                }
            )
    return output


def closest_taxon(tree, focal: str, candidates: list[str]) -> str:
    distances = [(float(tree.distance(focal, taxon)), taxon) for taxon in candidates if taxon != focal]
    distances.sort()
    return distances[0][1]


def cherry_partner(tree, focal: str) -> str:
    path = tree.get_path(focal)
    if not path:
        return ""
    parent = None
    for clade in tree.find_clades(order="level"):
        if any(child is path[-1] for child in clade.clades):
            parent = clade
            break
    if parent and len(parent.clades) == 2:
        leaves = [leaf.name for leaf in parent.get_terminals()]
        if len(leaves) == 2 and focal in leaves:
            return next(label for label in leaves if label != focal)
    return ""


def build_bornholm_check(windows: list[dict[str, str]], local_rows: list[dict[str, object]]) -> dict[str, object]:
    lg12 = sorted([row for row in windows if row["lg"] == "LG12"], key=lambda row: int(row["start"]))
    target_id = "LG12_007500001_007750000"
    index = next((i for i, row in enumerate(lg12) if row["window_id"] == target_id), None)
    if index is None:
        return {"recovered": False, "notes": "Target Bornholm window not found."}
    selected_indices = [i for i in (index - 1, index, index + 1) if 0 <= i < len(lg12)]
    stage2_taxa = sorted({str(row[f"taxon{i}"]) for row in local_rows if row["lg"] == "LG12" for i in range(1, 5)})
    summaries = []
    for i in selected_indices:
        window = lg12[i]
        tree = parse_tree(window["tree_newick"])
        bornholm_rows = [row for row in local_rows if row["window_id"] == window["window_id"] and "Gadmor_bor_spc" in [row[f"taxon{j}"] for j in range(1, 5)]]
        class_counts = Counter(row["local_class"] for row in bornholm_rows if row["informative_for_competing_predictions"] == "true")
        summaries.append(
            {
                "window_id": window["window_id"],
                "start": int(window["start"]),
                "end": int(window["end"]),
                "closest_taxon": closest_taxon(tree, "Gadmor_bor_spc", stage2_taxa),
                "cherry_partner": cherry_partner(tree, "Gadmor_bor_spc"),
                "bornholm_informative_arrangement": class_counts["arrangement"],
                "bornholm_informative_baseline": class_counts["baseline"],
                "bornholm_informative_third": class_counts["third"],
                "bornholm_informative_unresolved": class_counts["unresolved"],
            }
        )
    target_summary = next(item for item in summaries if item["window_id"] == target_id)
    neighbor_closest = [item["closest_taxon"] for item in summaries if item["window_id"] != target_id]
    recovered = bool(neighbor_closest) and any(label != target_summary["closest_taxon"] for label in neighbor_closest)
    return {"target_id": target_id, "summaries": summaries, "recovered": recovered, "notes": "Recovered means the target window's closest Bornholm relationship differs from at least one immediate neighboring window."}


def write_bornholm_check(check: dict[str, object]) -> None:
    if "summaries" not in check:
        BORNHOLM_CHECK.write_text(f"# Stage-4A Bornholm LG12 check\n\n{check['notes']}\n")
        return
    lines = [
        "# Stage-4A Bornholm LG12 positive-control check",
        "",
        "This targeted check inspects the Stage-1 target window `LG12_007500001_007750000` and its immediate neighboring published windows. It does not alter any frozen Stage-2 arrangement state.",
        "",
        "Relevant local relationships are summarized by Bornholm's closest taxon in the local tree, its cherry partner when present, and local classifications for informative quartets containing Bornholm.",
        "",
        "| window | start | end | closest taxon | cherry partner | arrangement | baseline | third | unresolved |",
        "|---|---:|---:|---|---|---:|---:|---:|---:|",
    ]
    for item in check["summaries"]:  # type: ignore[index]
        lines.append(
            f"| {item['window_id']} | {item['start']} | {item['end']} | {item['closest_taxon']} | {item['cherry_partner'] or 'none'} | {item['bornholm_informative_arrangement']} | {item['bornholm_informative_baseline']} | {item['bornholm_informative_third']} | {item['bornholm_informative_unresolved']} |"
        )
    lines.extend(
        [
            "",
            f"Reported switch recovered: {'yes' if check['recovered'] else 'no'}. {check['notes']}",
            "",
            "This is a positive-control validation only and is not part of the primary Stage-4A contrast.",
        ]
    )
    BORNHOLM_CHECK.write_text("\n".join(lines) + "\n")


def plot_lg_tracks(window_scores: list[dict[str, object]], regions: dict[str, dict[str, int | str]]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for lg in LGS:
        rows = sorted([row for row in window_scores if row["lg"] == lg], key=lambda row: int(row["midpoint"]))
        x = [int(row["midpoint"]) / 1_000_000 for row in rows]
        d = [float(row["D_arrangement_minus_baseline"]) for row in rows]
        fa = [float(row["fraction_arrangement"]) for row in rows]
        fb = [float(row["fraction_baseline"]) for row in rows]
        f3 = [float(row["fraction_third"]) for row in rows]
        boundary = [bool_from_tsv(str(row["overlaps_inversion"])) and not bool_from_tsv(str(row["inside_inversion"])) for row in rows]
        fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.2), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
        for ax in axes:
            ax.axvspan(int(regions[lg]["start"]) / 1_000_000, int(regions[lg]["end"]) / 1_000_000, color="#d9d9d9", alpha=0.6, label="frozen inversion")
        axes[0].axhline(0, color="black", linewidth=0.8)
        axes[0].plot(x, d, color="#1f77b4", linewidth=1.2, marker="o", markersize=2.6, label="D = fA - fB")
        axes[0].scatter([xx for xx, b in zip(x, boundary) if b], [yy for yy, b in zip(d, boundary) if b], color="#d62728", s=24, zorder=3, label="boundary-overlap")
        axes[0].set_ylabel("D")
        axes[0].set_ylim(-1.05, 1.05)
        axes[0].set_title(f"Atlantic cod {lg}: local quartet topology track")
        axes[0].legend(loc="upper right", fontsize=8)
        axes[1].plot(x, fa, color="#2ca02c", linewidth=1.0, label="arrangement")
        axes[1].plot(x, fb, color="#9467bd", linewidth=1.0, label="baseline")
        axes[1].plot(x, f3, color="#8c564b", linewidth=1.0, label="third")
        axes[1].set_ylabel("fraction")
        axes[1].set_ylim(-0.03, 1.03)
        axes[1].set_xlabel("position (Mb)")
        axes[1].legend(loc="upper right", fontsize=8, ncol=3)
        fig.tight_layout()
        fig.savefig(LG_FIGURES[lg])
        plt.close(fig)


def plot_four_lg_summary(window_scores: list[dict[str, object]], regions: dict[str, dict[str, int | str]]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(8.0, 8.2), sharey=True)
    for ax, lg in zip(axes, LGS):
        rows = sorted([row for row in window_scores if row["lg"] == lg], key=lambda row: int(row["midpoint"]))
        x = [int(row["midpoint"]) / 1_000_000 for row in rows]
        d = [float(row["D_arrangement_minus_baseline"]) for row in rows]
        boundary = [bool_from_tsv(str(row["overlaps_inversion"])) and not bool_from_tsv(str(row["inside_inversion"])) for row in rows]
        ax.axvspan(int(regions[lg]["start"]) / 1_000_000, int(regions[lg]["end"]) / 1_000_000, color="#d9d9d9", alpha=0.6)
        ax.axhline(0, color="black", linewidth=0.7)
        ax.plot(x, d, color="#1f77b4", linewidth=1.1, marker="o", markersize=2.2)
        ax.scatter([xx for xx, b in zip(x, boundary) if b], [yy for yy, b in zip(d, boundary) if b], color="#d62728", s=20, zorder=3)
        ax.set_ylabel(lg)
        ax.set_ylim(-1.05, 1.05)
    axes[-1].set_xlabel("position (Mb)")
    fig.suptitle("Atlantic cod Stage 4A: D = arrangement fraction - baseline fraction", y=0.995)
    fig.tight_layout()
    fig.savefig(FOUR_LG_FIGURE)
    plt.close(fig)


def write_manifest(stage2_sha: str, stage3_sha: str, local_sha: str, output_paths: list[Path], windows: list[dict[str, str]], quartets_by_lg: dict[str, list[dict[str, str]]], local_rows: list[dict[str, object]], regions: dict[str, dict[str, int | str]]) -> None:
    manifest = {
        "stage": "4A",
        "stage2_checksum": stage2_sha,
        "stage3_checksum": stage3_sha,
        "local_window_table": str(LOCAL_WINDOWS.relative_to(REPO_ROOT)),
        "local_window_table_sha256": local_sha,
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "number_of_windows": len(windows),
        "number_informative_quartets_per_lg": {lg: sum(row["informative_for_competing_predictions"] == "true" for row in quartets_by_lg[lg]) for lg in LGS},
        "number_of_local_quartet_observations": len(local_rows),
        "output_checksums": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in output_paths if path.exists()},
        "inversion_coordinates_used": regions,
        "significance_tests_run": False,
        "anti_circularity_note": "Stage 4A used only frozen Stage-2/3 predictions and published local trees after checksum verification; no permutation/spatial-null tests were run.",
    }
    STAGE4A_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_report(stage2_sha: str, stage3_sha: str, local_sha: str, lg_summary: list[dict[str, object]], composition: list[dict[str, object]], quartet_effects: list[dict[str, object]], bornholm: dict[str, object], boundary_sensitivity: list[dict[str, object]], parse_issues: list[str]) -> None:
    summary_lines = []
    for row in lg_summary:
        summary_lines.append(f"| {row['lg']} | {row['mean_D_inside']:.6f} | {row['mean_D_outside']:.6f} | {row['delta_D']:.6f} | {row['median_D_inside']:.6f} | {row['median_D_outside']:.6f} | {row['fraction_quartets_positive_effect']:.6f} |")
    comp_lines = []
    for row in composition:
        if row["region_class"] in {"inside", "outside"}:
            comp_lines.append(f"| {row['lg']} | {row['region_class']} | {row['n_windows']} | {row['mean_fraction_arrangement']:.6f} | {row['mean_fraction_baseline']:.6f} | {row['mean_fraction_third']:.6f} |")
    q_summary = []
    for lg in LGS:
        vals = [float(row["delta_arrangement_minus_baseline"]) for row in quartet_effects if row["lg"] == lg and not math.isnan(float(row["delta_arrangement_minus_baseline"]))]
        pos = sum(value > 0 for value in vals)
        q_summary.append(f"| {lg} | {pos}/{len(vals)} | {(pos / len(vals)) if vals else math.nan:.6f} | {statistics.median(vals) if vals else math.nan:.6f} | {min(vals) if vals else math.nan:.6f} | {max(vals) if vals else math.nan:.6f} |")
    sens_lines = []
    for row in boundary_sensitivity:
        sens_lines.append(f"| {row['lg']} | {row['boundary_treatment']} | {row['delta_D']:.6f} | {row['mean_D_inside']:.6f} | {row['mean_D_outside']:.6f} |")
    visual = [row for row in lg_summary if float(row["delta_D"]) > 0]
    lines = [
        "# Atlantic cod Stage-4A report",
        "",
        "Stage 4A compares published local 250-kb population trees with the independently frozen arrangement and baseline quartet predictions. It is descriptive only; no spatial-null or permutation test was run.",
        "",
        "## Input verification",
        "",
        f"- Stage-2 structural candidate checksum: `{stage2_sha}`.",
        f"- Stage-3 baseline quartet checksum: `{stage3_sha}`.",
        f"- Local-window table checksum: `{local_sha}`.",
        "",
        "## Primary result",
        "",
        "`D = fraction_arrangement - fraction_baseline`, computed per window using informative quartets only.",
        "",
        "| LG | mean D inside | mean D outside | inside - outside | median D inside | median D outside | fraction positive informative quartets |",
        "|---|---:|---:|---:|---:|---:|---:|",
        *summary_lines,
        "",
        "Inside inversion windows show greater alignment with the independently frozen arrangement partition than outside windows for the LGs with positive inside-minus-outside differences. This is descriptive and not a significance claim.",
        "",
        "## Topology composition",
        "",
        "| LG | region | windows | mean arrangement fraction | mean baseline fraction | mean third fraction |",
        "|---|---|---:|---:|---:|---:|",
        *comp_lines,
        "",
        "## Quartet consistency",
        "",
        "| LG | positive informative quartets | fraction positive | median effect | min effect | max effect |",
        "|---|---:|---:|---:|---:|---:|",
        *q_summary,
        "",
        "## Spatial pattern",
        "",
        f"Positive inside-minus-outside D values were observed for {len(visual)} of 4 LGs in the descriptive window-level summary. The topology-track PDFs shade the frozen inversion intervals and mark boundary-overlap windows; visual alignment should be assessed from those tracks before Stage 4B.",
        "",
        "## Bornholm positive control",
        "",
        f"Reported LG12 local-switch behavior recovered: {'yes' if bornholm.get('recovered') else 'no'}. Details are in `results/stage4a_bornholm_check.md`.",
        "",
        "## Boundary sensitivity",
        "",
        "| LG | boundary treatment | delta D | mean D inside | mean D outside |",
        "|---|---|---:|---:|---:|",
        *sens_lines,
        "",
        "The primary comparison excludes partial boundary-overlap windows from both inside and outside classes. The sensitivity table records alternative treatments without selecting a favorable definition.",
        "",
        "## Caveats",
        "",
        "- Adjacent 250-kb windows are spatially linked, so Stage 4A does not use iid tests.",
        "- Multiple quartets from the same local tree are dependent; the primary unit is the genomic window.",
        "- LG12 has only 16 informative quartets and is less powered than LG01, LG02, and LG07.",
        "- LG12 ancestral/derived orientation follows the weaker Stage-2 demographic-inference evidence, unlike the stronger direct outgroup evidence for LG01, LG02, and LG07.",
        "- This is a multi-population Atlantic cod analysis, not a direct multi-species persistence-through-speciation test.",
        "",
        "## Parsing issues",
        "",
        "No unexpected topology/parsing issues were encountered." if not parse_issues else "\n".join(f"- {issue}" for issue in parse_issues),
        "",
        "## Stage boundary",
        "",
        "No circular-shift tests, block permutations, P-values, formal MSRC model fitting, or SNAPP/BEAST reruns were performed.",
        "",
    ]
    STAGE4A_REPORT.write_text("\n".join(lines))


def main() -> None:
    stage2_sha = verify_checksum(STAGE2_CANDIDATES, EXPECTED_STAGE2_SHA, "Stage-2 structural candidates")
    stage3_sha = verify_checksum(STAGE3_QUARTETS, EXPECTED_STAGE3_SHA, "Stage-3 baseline quartets")
    local_sha = sha256(LOCAL_WINDOWS)
    quartets_by_lg = read_stage3_quartets()
    windows = read_windows()
    regions = read_regions()
    local_rows, parse_issues = build_local_quartets(windows, quartets_by_lg)
    window_scores = compute_window_scores(local_rows, windows)
    inside_outside = build_inside_outside_summary(window_scores)
    composition = build_topology_composition(window_scores)
    quartet_effects = build_quartet_effects(local_rows)
    lg_effect = build_lg_effect_summary(window_scores, quartet_effects)
    boundary_desc = build_boundary_descriptive(window_scores, regions)
    boundary_sens = build_boundary_sensitivity(window_scores)
    bornholm = build_bornholm_check(windows, local_rows)

    write_tsv(LOCAL_QUARTETS, local_rows, LOCAL_QUARTET_COLUMNS)
    write_tsv(WINDOW_SCORES, window_scores, WINDOW_SCORE_COLUMNS)
    write_tsv(INSIDE_OUTSIDE_SUMMARY, inside_outside, ["lg", "region_class", "D_n_windows", "D_mean", "D_median", "D_iqr", "D_min", "D_max"])
    write_tsv(TOPOLOGY_COMPOSITION, composition, ["lg", "region_class", "n_windows", "mean_fraction_arrangement", "mean_fraction_baseline", "mean_fraction_third"])
    write_tsv(QUARTET_EFFECTS, quartet_effects, ["lg", "quartet_id", "arrangement_split", "baseline_split", "n_inside", "n_outside", "arrangement_fraction_inside", "arrangement_fraction_outside", "baseline_fraction_inside", "baseline_fraction_outside", "delta_arrangement", "delta_arrangement_minus_baseline"])
    write_tsv(BOUNDARY_DESCRIPTIVE, boundary_desc, ["lg", "boundary", "boundary_position", "window_id", "start", "end", "midpoint", "inside_inversion", "overlaps_inversion", "distance_midpoint_to_boundary", "D_arrangement_minus_baseline"])
    write_tsv(BOUNDARY_SENSITIVITY, boundary_sens, ["lg", "boundary_treatment", "n_inside_windows", "n_outside_windows", "mean_D_inside", "mean_D_outside", "delta_D", "median_D_inside", "median_D_outside"])
    write_tsv(LG_EFFECT_SUMMARY, lg_effect, ["lg", "n_informative_quartets", "n_inside_windows", "n_outside_windows", "mean_D_inside", "mean_D_outside", "delta_D", "median_D_inside", "median_D_outside", "fraction_quartets_positive_effect"])
    write_bornholm_check(bornholm)
    plot_lg_tracks(window_scores, regions)
    plot_four_lg_summary(window_scores, regions)
    write_report(stage2_sha, stage3_sha, local_sha, lg_effect, composition, quartet_effects, bornholm, boundary_sens, parse_issues)
    output_paths = [LOCAL_QUARTETS, WINDOW_SCORES, INSIDE_OUTSIDE_SUMMARY, TOPOLOGY_COMPOSITION, QUARTET_EFFECTS, BOUNDARY_DESCRIPTIVE, BOUNDARY_SENSITIVITY, LG_EFFECT_SUMMARY, BORNHOLM_CHECK, STAGE4A_REPORT, *LG_FIGURES.values(), FOUR_LG_FIGURE]
    write_manifest(stage2_sha, stage3_sha, local_sha, output_paths, windows, quartets_by_lg, local_rows, regions)


class Stage4ATests(unittest.TestCase):
    def parse(self, newick: str):
        return parse_tree(newick)

    def test_local_quartet_extraction(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        self.assertEqual(induced_quartet_split(tree, ["A", "B", "C", "D"]), "A,B|C,D")

    def test_canonical_split_representation(self) -> None:
        self.assertEqual(canonical_split(["B", "A"], ["D", "C"]), "A,B|C,D")
        self.assertEqual(canonicalize_existing_split("D,C|B,A"), "A,B|C,D")

    def test_arrangement_classification(self) -> None:
        self.assertEqual(classify_local_split("A,B|C,D", "A,B|C,D", "A,C|B,D"), "arrangement")

    def test_baseline_classification(self) -> None:
        self.assertEqual(classify_local_split("A,C|B,D", "A,B|C,D", "A,C|B,D"), "baseline")

    def test_third_topology_classification(self) -> None:
        self.assertEqual(classify_local_split("A,D|B,C", "A,B|C,D", "A,C|B,D"), "third")

    def test_unresolved_handling(self) -> None:
        tree = self.parse("(A:1,B:1,C:1,D:1);")
        with self.assertRaises(UnresolvedQuartetError):
            induced_quartet_split(tree, ["A", "B", "C", "D"])
        self.assertEqual(classify_local_split("unresolved", "A,B|C,D", "A,C|B,D"), "unresolved")

    def test_window_level_fraction_calculations_and_d(self) -> None:
        windows = [{"lg": "LGX", "window_id": "W1", "start": "1", "end": "10", "midpoint": "5", "inside_inversion": "true", "overlaps_inversion": "true"}]
        rows = [
            {"lg": "LGX", "window_id": "W1", "informative_for_competing_predictions": "true", "local_class": "arrangement"},
            {"lg": "LGX", "window_id": "W1", "informative_for_competing_predictions": "true", "local_class": "baseline"},
            {"lg": "LGX", "window_id": "W1", "informative_for_competing_predictions": "true", "local_class": "third"},
            {"lg": "LGX", "window_id": "W1", "informative_for_competing_predictions": "true", "local_class": "unresolved"},
            {"lg": "LGX", "window_id": "W1", "informative_for_competing_predictions": "false", "local_class": "arrangement"},
        ]
        out = compute_window_scores(rows, windows)[0]
        self.assertEqual(out["n_resolved"], 3)
        self.assertAlmostEqual(out["fraction_arrangement"], 1 / 3)
        self.assertAlmostEqual(out["fraction_baseline"], 1 / 3)
        self.assertAlmostEqual(out["fraction_third"], 1 / 3)
        self.assertAlmostEqual(out["D_arrangement_minus_baseline"], 0.0)

    def test_stage2_checksum_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "x.tsv"
            path.write_text("abc\n")
            expected = hashlib.sha256(b"abc\n").hexdigest()
            self.assertEqual(verify_checksum(path, expected, "x"), expected)
            with self.assertRaises(ValueError):
                verify_checksum(path, "0" * 64, "x")

    def test_stage3_checksum_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "y.tsv"
            path.write_text("def\n")
            expected = hashlib.sha256(b"def\n").hexdigest()
            self.assertEqual(verify_checksum(path, expected, "y"), expected)
            with self.assertRaises(ValueError):
                verify_checksum(path, "1" * 64, "y")

    def test_inversion_membership_preservation(self) -> None:
        row = {"inside_inversion": "false", "overlaps_inversion": "true"}
        self.assertEqual(region_class(row), "boundary-overlap")
        self.assertEqual(region_class(row, "inside"), "inside")
        self.assertEqual(region_class(row, "outside"), "outside")

    def test_deterministic_output(self) -> None:
        tree = self.parse("((A:1,B:1):1,(C:1,D:1):1);")
        windows = [{"lg": "LGX", "window_id": "W1", "start": "1", "end": "10", "midpoint": "5", "inside_inversion": "false", "overlaps_inversion": "false", "tree_newick": "((A:1,B:1):1,(C:1,D:1):1);", "source_file": "test"}]
        quartets = {"LGX": [{"quartet_id": "Q1", "taxon1": "A", "taxon2": "B", "taxon3": "C", "taxon4": "D", "arrangement_split": "A,B|C,D", "baseline_split": "A,C|B,D", "informative_for_competing_predictions": "true"}]}
        one, _ = build_local_quartets(windows, quartets)
        two, _ = build_local_quartets(windows, quartets)
        self.assertEqual(one, two)
        self.assertEqual(one[0]["local_class"], "arrangement")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="Run self-tests and exit.")
    args = parser.parse_args()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage4ATests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        sys.exit(0 if result.wasSuccessful() else 1)
    main()
