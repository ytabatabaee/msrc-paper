#!/usr/bin/env python3
"""Blind chr4 genealogy change-point proximity to frozen de novo breakpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from statistics import median
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"

LOCUS_PATH = RESULTS_DIR / "locus_table_chr4.tsv"
BREAKPOINT_SETS_PATH = RESULTS_DIR / "denovo_breakpoint_sets.tsv"
CONSENSUS_BREAKPOINTS_PATH = RESULTS_DIR / "denovo_consensus_structural_breakpoints.tsv"
MANIFEST_PATH = RESULTS_DIR / "denovo_breakpoint_manifest.json"

RAW_CP_PATH = RESULTS_DIR / "genealogy_change_points_raw.tsv"
CP_PATH = RESULTS_DIR / "genealogy_change_points.tsv"
SET_SUMMARY_PATH = RESULTS_DIR / "breakpoint_enrichment_by_structural_set.tsv"
METHOD_SUMMARY_PATH = RESULTS_DIR / "breakpoint_enrichment_by_genealogy_method.tsv"
SENSITIVITY_PATH = RESULTS_DIR / "breakpoint_enrichment_sensitivity.tsv"
BOUNDARY_CENTRIC_PATH = RESULTS_DIR / "structural_breakpoint_genealogy_support.tsv"
SHARED_TRANSITIONS_PATH = RESULTS_DIR / "shared_genealogy_transitions.tsv"
SUMMARY_FINAL_PATH = RESULTS_DIR / "breakpoint_enrichment_summary_final.tsv"
BLIND_PDF_PATH = FIGURES_DIR / "chr4_change_points_blind.pdf"
BLIND_PNG_PATH = FIGURES_DIR / "chr4_change_points_blind.png"
STRUCTURE_PDF_PATH = FIGURES_DIR / "chr4_change_points_vs_denovo_structure.pdf"
STRUCTURE_PNG_PATH = FIGURES_DIR / "chr4_change_points_vs_denovo_structure.png"
NULL_PDF_PATH = FIGURES_DIR / "breakpoint_enrichment_null_primary.pdf"
SENSITIVITY_PDF_PATH = FIGURES_DIR / "breakpoint_enrichment_structural_set_sensitivity.pdf"
CLUSTER_AUDIT_PATH = RESULTS_DIR / "genealogy_change_point_cluster_audit.tsv"
SINGLELINK_COMPARISON_PATH = RESULTS_DIR / "stage3a_singlelink_vs_completespan.tsv"
FINAL_MANIFEST_PATH = RESULTS_DIR / "stage3a_final_manifest.json"

CLADES = ["Columbea", "N61", "N62"]
Q_COLS = ("q1", "q2", "q3")
THRESHOLDS = (50_000, 100_000, 250_000, 500_000)
STRUCTURAL_SETS = ("primary", "stringent", "inclusive")
GENEALOGY_METHODS = ("rolling", "binary_segmentation")
PRIMARY_METHOD = "rolling"
PRIMARY_STRUCTURAL_SET = "primary"
OLD_BOUNDARY_INPUTS = {
    "canonical_structural_regions.tsv",
    "canonical_structural_boundaries.tsv",
    "structural_intervals.tsv",
    "outlier-regions-by-maf2synteny.tsv",
    "outlier-regions-by-maf2synteny-summary.tsv",
}
READ_PATHS: list[Path] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    assert_not_old_boundary_input(path)
    READ_PATHS.append(path)
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict[str, object]:
    assert_not_old_boundary_input(path)
    READ_PATHS.append(path)
    with path.open() as handle:
        return json.load(handle)


def assert_not_old_boundary_input(path: Path) -> None:
    if path.name in OLD_BOUNDARY_INPUTS:
        raise AssertionError(f"Stage 3A must not read old structural boundary input: {path}")


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: fmt(row.get(column, "")) for column in columns})


def fmt(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.12g}"
    return str(value)


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def load_loci() -> dict[str, list[dict[str, object]]]:
    rows = read_tsv(LOCUS_PATH)
    if any(row["Chromosome"] != "chr4" for row in rows):
        raise ValueError("Stage 3A requires exact Chromosome == 'chr4' loci only")
    by_clade: dict[str, list[dict[str, object]]] = {clade: [] for clade in CLADES}
    for row in rows:
        if row["clade"] not in by_clade:
            continue
        by_clade[row["clade"]].append(
            {
                "clade": row["clade"],
                "position": float(row["midpoint"]),
                "q": tuple(float(row[q]) for q in Q_COLS),
            }
        )
    for clade in CLADES:
        by_clade[clade].sort(key=lambda row: row["position"])
    return by_clade


def load_frozen_structural_breakpoints() -> dict[str, list[dict[str, object]]]:
    set_rows = read_tsv(BREAKPOINT_SETS_PATH)
    consensus_rows = {row["consensus_breakpoint_id"]: row for row in read_tsv(CONSENSUS_BREAKPOINTS_PATH)}
    manifest = read_json(MANIFEST_PATH)
    out: dict[str, list[dict[str, object]]] = {name: [] for name in STRUCTURAL_SETS}
    for row in set_rows:
        bp_id = row["consensus_breakpoint_id"]
        consensus = consensus_rows[bp_id]
        enriched = {
            "consensus_breakpoint_id": bp_id,
            "position": int(row["reference_position"]),
            "n_species_support": int(row["n_species_support"]),
            "confidence_class": row["confidence_class"],
            "recovery_fraction": float(consensus["parameter_recovery_fraction"]),
        }
        if enriched["confidence_class"] == "strong_consensus":
            assert truthy(row["in_stringent"])
        if enriched["n_species_support"] >= 2:
            assert truthy(row["in_primary"])
        assert truthy(row["in_inclusive"])
        for name in STRUCTURAL_SETS:
            if truthy(row[f"in_{name}"]):
                out[name].append(dict(enriched))
    for name in STRUCTURAL_SETS:
        out[name].sort(key=lambda item: int(item["position"]))
    expected = {
        "primary": manifest.get("primary_breakpoint_count"),
        "stringent": manifest.get("stringent_breakpoint_count"),
        "inclusive": manifest.get("inclusive_breakpoint_count"),
    }
    for name, count in expected.items():
        if count is not None and len(out[name]) != int(count):
            raise AssertionError(f"{name} frozen set has {len(out[name])} rows, manifest says {count}")
    return out


def prefix_vectors(qs: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    prefix = [(0.0, 0.0, 0.0)]
    total = [0.0, 0.0, 0.0]
    for q in qs:
        for i in range(3):
            total[i] += q[i]
        prefix.append(tuple(total))
    return prefix


def mean_between(prefix: list[tuple[float, float, float]], start: int, end: int) -> tuple[float, float, float]:
    n = end - start
    return tuple((prefix[end][i] - prefix[start][i]) / n for i in range(3))


def vec_delta(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(right[i] - left[i] for i in range(3))


def vec_norm(delta: tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in delta))


def dominant(q: tuple[float, float, float]) -> str:
    return max(zip(Q_COLS, q), key=lambda item: (item[1], item[0]))[0]


def rolling_candidates(track: list[dict[str, object]], window: int) -> list[dict[str, object]]:
    qs = [row["q"] for row in track]
    prefix = prefix_vectors(qs)
    candidates = []
    scores = []
    for i in range(window, len(track) - window):
        left = mean_between(prefix, i - window, i)
        right = mean_between(prefix, i, i + window)
        score = vec_norm(vec_delta(left, right))
        scores.append(score)
        candidates.append(candidate_row(track[i]["clade"], "rolling", i, track[i]["position"], score, left, right))
    if not scores:
        return []
    thresh = data_threshold(scores)
    raw = []
    for j, row in enumerate(candidates):
        if row["score"] < thresh:
            continue
        prev_score = candidates[j - 1]["score"] if j > 0 else -1
        next_score = candidates[j + 1]["score"] if j + 1 < len(candidates) else -1
        if row["score"] >= prev_score and row["score"] >= next_score:
            raw.append(row)
    return raw


def data_threshold(values: list[float]) -> float:
    med = median(values)
    mad = median([abs(value - med) for value in values])
    return med + 4.0 * (mad if mad > 0 else 1e-12)


def candidate_row(
    clade: str,
    method: str,
    index: int,
    position: float,
    score: float,
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> dict[str, object]:
    delta = vec_delta(left, right)
    return {
        "clade": clade,
        "method": method,
        "index": index,
        "position": int(round(position)),
        "score": score,
        "left_q1": left[0],
        "left_q2": left[1],
        "left_q3": left[2],
        "right_q1": right[0],
        "right_q2": right[1],
        "right_q3": right[2],
        "delta_q1": delta[0],
        "delta_q2": delta[1],
        "delta_q3": delta[2],
        "dominant_topology_left": dominant(left),
        "dominant_topology_right": dominant(right),
    }


def binary_segmentation_candidates(track: list[dict[str, object]], min_size: int, max_changes: int = 20) -> list[dict[str, object]]:
    qs = [row["q"] for row in track]
    prefix = prefix_vectors(qs)
    prefix_sq = prefix_squares(qs)
    out = []
    queue = [(0, len(track))]
    while queue and len(out) < max_changes:
        start, end = queue.pop(0)
        if end - start < 2 * min_size:
            continue
        parent = segment_sse(prefix, prefix_sq, start, end)
        best = None
        for split in range(start + min_size, end - min_size):
            gain = parent - segment_sse(prefix, prefix_sq, start, split) - segment_sse(prefix, prefix_sq, split, end)
            if best is None or gain > best[0]:
                best = (gain, split)
        if best is None:
            continue
        gain, split = best
        penalty = 3.0 * math.log(end - start)
        if gain <= penalty:
            continue
        left = mean_between(prefix, start, split)
        right = mean_between(prefix, split, end)
        out.append(candidate_row(track[split]["clade"], "binary_segmentation", split, track[split]["position"], gain, left, right))
        queue.extend([(start, split), (split, end)])
    return sorted(out, key=lambda row: row["position"])


def prefix_squares(qs: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    prefix = [(0.0, 0.0, 0.0)]
    total = [0.0, 0.0, 0.0]
    for q in qs:
        for i in range(3):
            total[i] += q[i] * q[i]
        prefix.append(tuple(total))
    return prefix


def segment_sse(prefix: list[tuple[float, float, float]], prefix_sq: list[tuple[float, float, float]], start: int, end: int) -> float:
    n = end - start
    total = 0.0
    for i in range(3):
        s = prefix[end][i] - prefix[start][i]
        ss = prefix_sq[end][i] - prefix_sq[start][i]
        total += ss - s * s / n
    return total


def merge_change_points(rows: list[dict[str, object]], distance: int) -> list[dict[str, object]]:
    merge_distance_bp = int(distance)
    merged = []
    for (clade, method), subset in grouped(rows, ("clade", "method")).items():
        ordered = sorted(subset, key=lambda row: int(row["position"]))
        cluster: list[dict[str, object]] = []
        cluster_min_position = None
        cluster_max_position = None

        def close_cluster() -> None:
            assert cluster
            positions = [int(item["position"]) for item in cluster]
            raw_ids = [str(item.get("change_point_id", "")) for item in cluster]
            span = max(positions) - min(positions)
            assert span <= merge_distance_bp
            representative = dict(max(cluster, key=lambda r: float(r["score"])))
            representative["cluster_min_position"] = min(positions)
            representative["cluster_max_position"] = max(positions)
            representative["cluster_span_bp"] = span
            representative["n_raw_change_points"] = len(cluster)
            representative["contributing_raw_change_point_ids"] = ";".join(raw_ids)
            merged.append(representative)

        for row in ordered:
            candidate_position = int(row["position"])
            if not cluster:
                cluster.append(row)
                cluster_min_position = candidate_position
                cluster_max_position = candidate_position
            elif candidate_position - int(cluster_min_position) <= merge_distance_bp:
                cluster.append(row)
                cluster_max_position = max(int(cluster_max_position), candidate_position)
            else:
                close_cluster()
                cluster = [row]
                cluster_min_position = candidate_position
                cluster_max_position = candidate_position
        if cluster:
            close_cluster()
    for key, subset in grouped(merged, ("clade", "method")).items():
        for i, row in enumerate(sorted(subset, key=lambda r: int(r["position"])), 1):
            row["change_point_id"] = f"{key[0]}_{key[1]}_{i:03d}"
    return sorted(merged, key=lambda row: (row["method"], row["clade"], int(row["position"])))


def grouped(rows: list[dict[str, object]], keys: tuple[str, ...]) -> dict[tuple[object, ...], list[dict[str, object]]]:
    out = defaultdict(list)
    for row in rows:
        out[tuple(row[key] for key in keys)].append(row)
    return out


def nearest_distances(positions: list[int], breakpoints: list[int]) -> list[int]:
    if not breakpoints:
        raise ValueError("nearest-distance calculation requires at least one structural breakpoint")
    return [min(abs(pos - bp) for bp in breakpoints) for pos in positions]


def proximity_stats(positions: list[int], breakpoints: list[int]) -> dict[str, float]:
    if not positions:
        return empty_stats()
    distances = nearest_distances(positions, breakpoints)
    out = {
        "mean_nearest_bp": sum(distances) / len(distances),
        "median_nearest_bp": median(distances),
    }
    for threshold in THRESHOLDS:
        out[f"frac_within_{threshold//1000}kb"] = sum(d <= threshold for d in distances) / len(distances)
    return out


def empty_stats() -> dict[str, float]:
    out = {"mean_nearest_bp": math.nan, "median_nearest_bp": math.nan}
    for threshold in THRESHOLDS:
        out[f"frac_within_{threshold//1000}kb"] = math.nan
    return out


def circular_shift_positions(positions: list[int], shift: int, start: int, end: int) -> list[int]:
    length = end - start
    return [start + ((pos - start + shift) % length) for pos in positions]


def circular_null(
    positions: list[int],
    breakpoints: list[int],
    start: int,
    end: int,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    rng = random.Random(seed)
    observed = proximity_stats(positions, breakpoints)
    null_means = []
    null_fracs = {threshold: [] for threshold in THRESHOLDS}
    for _ in range(permutations):
        shifted = circular_shift_positions(positions, rng.randrange(0, end - start), start, end)
        stats = proximity_stats(shifted, breakpoints)
        null_means.append(stats["mean_nearest_bp"])
        for threshold in THRESHOLDS:
            null_fracs[threshold].append(stats[f"frac_within_{threshold//1000}kb"])
    pvalue = (1 + sum(value <= observed["mean_nearest_bp"] for value in null_means)) / (permutations + 1)
    enrichments = {}
    for threshold in THRESHOLDS:
        obs_frac = observed[f"frac_within_{threshold//1000}kb"]
        expected = sum(null_fracs[threshold]) / len(null_fracs[threshold])
        enrichments[threshold] = obs_frac / expected if expected > 0 else math.nan
    null_mean = sum(null_means) / len(null_means)
    return {
        "observed": observed,
        "null_means": null_means,
        "null_mean_nearest_bp": null_mean,
        "null_median_nearest_bp": median(null_means),
        "observed_to_null_ratio": observed["mean_nearest_bp"] / null_mean if null_mean > 0 else math.nan,
        "pvalue": pvalue,
        "enrichments": enrichments,
    }


def detect_all(by_clade: dict[str, list[dict[str, object]]], window: int) -> list[dict[str, object]]:
    rows = []
    for track in by_clade.values():
        rows.extend(rolling_candidates(track, window))
        rows.extend(binary_segmentation_candidates(track, min_size=window))
    return rows


def add_raw_ids(rows: list[dict[str, object]]) -> None:
    for key, subset in grouped(rows, ("clade", "method")).items():
        for i, row in enumerate(sorted(subset, key=lambda r: int(r["position"])), 1):
            row["change_point_id"] = f"{key[0]}_{key[1]}_raw_{i:03d}"


def usable_range(by_clade: dict[str, list[dict[str, object]]]) -> tuple[int, int]:
    positions = [int(row["position"]) for track in by_clade.values() for row in track]
    return min(positions), max(positions)


def positions_for(merged: list[dict[str, object]], method: str, clade: str) -> list[int]:
    return [
        int(row["position"]) for row in merged
        if row["method"] == method and (clade == "pooled" or row["clade"] == clade)
    ]


def summarize_enrichment(
    merged: list[dict[str, object]],
    method: str,
    structural_set: str,
    breakpoints: list[dict[str, object]],
    by_clade: dict[str, list[dict[str, object]]],
    permutations: int,
    seed: int,
) -> tuple[list[dict[str, object]], dict[str, list[float]]]:
    bp_positions = [int(row["position"]) for row in breakpoints]
    start, end = usable_range(by_clade)
    summaries = []
    nulls_for_plot = {}
    for label in [*CLADES, "pooled"]:
        positions = positions_for(merged, method, label)
        result = circular_null(positions, bp_positions, start, end, permutations, seed + stable_offset(label, method, structural_set))
        row = {
            "structural_set": structural_set,
            "clade": label,
            "genealogy_method": method,
            "n_genealogy_change_points": len(positions),
            "n_structural_breakpoints": len(bp_positions),
            **result["observed"],
            "observed_mean_nearest_bp": result["observed"]["mean_nearest_bp"],
            "null_mean_nearest_bp": result["null_mean_nearest_bp"],
            "null_median_nearest_bp": result["null_median_nearest_bp"],
            "observed_to_null_ratio": result["observed_to_null_ratio"],
            "circular_shift_pvalue": result["pvalue"],
        }
        for threshold in THRESHOLDS:
            row[f"enrichment_{threshold//1000}kb"] = result["enrichments"][threshold]
        summaries.append(row)
        nulls_for_plot[label] = result["null_means"]
    return summaries, nulls_for_plot


def stable_offset(*parts: str) -> int:
    return sum((i + 1) * ord(ch) for part in parts for i, ch in enumerate(part))


def structural_breakpoint_genealogy_support(
    breakpoints: list[dict[str, object]], merged: list[dict[str, object]], method: str
) -> list[dict[str, object]]:
    by_clade_positions = {clade: positions_for(merged, method, clade) for clade in CLADES}
    rows = []
    for bp in breakpoints:
        row = dict(bp)
        distances = {}
        for clade in CLADES:
            distances[clade] = nearest_distances([int(bp["position"])], by_clade_positions[clade])[0] if by_clade_positions[clade] else math.nan
            row[f"nearest_{clade}_bp"] = distances[clade]
        for threshold in (100_000, 250_000, 500_000):
            row[f"n_clades_within_{threshold//1000}kb"] = sum(distances[clade] <= threshold for clade in CLADES)
        rows.append(row)
    return rows


def cluster_shared_genealogy_transitions(
    merged: list[dict[str, object]], primary_breakpoints: list[dict[str, object]], cluster_bp: int
) -> list[dict[str, object]]:
    merge_distance_bp = int(cluster_bp)
    events = [
        {"clade": row["clade"], "position": int(row["position"]), "score": float(row["score"])}
        for row in merged
        if row["method"] == PRIMARY_METHOD and row["clade"] in CLADES
    ]
    events.sort(key=lambda row: row["position"])
    clusters = []
    current = []
    cluster_min_position = None
    cluster_max_position = None

    def close_cluster() -> None:
        assert current
        positions = [int(event["position"]) for event in current]
        span = max(positions) - min(positions)
        assert span <= merge_distance_bp
        clusters.append(list(current))

    for event in events:
        candidate_position = int(event["position"])
        if not current:
            current.append(event)
            cluster_min_position = candidate_position
            cluster_max_position = candidate_position
        elif candidate_position - int(cluster_min_position) <= merge_distance_bp:
            current.append(event)
            cluster_max_position = max(int(cluster_max_position), candidate_position)
        else:
            close_cluster()
            current = [event]
            cluster_min_position = candidate_position
            cluster_max_position = candidate_position
    if current:
        close_cluster()
    bp_positions = [int(row["position"]) for row in primary_breakpoints]
    rows = []
    for i, cluster in enumerate(clusters, 1):
        best_by_clade: dict[str, dict[str, object]] = {}
        for event in cluster:
            clade = str(event["clade"])
            if clade not in best_by_clade or float(event["score"]) > float(best_by_clade[clade]["score"]):
                best_by_clade[clade] = event
        retained = sorted(best_by_clade.values(), key=lambda event: int(event["position"]))
        positions = [int(event["position"]) for event in retained]
        consensus = int(round(median(positions)))
        nearest = min(bp_positions, key=lambda bp: abs(consensus - bp))
        clades = sorted({str(event["clade"]) for event in retained}, key=CLADES.index)
        cluster_min = min(positions)
        cluster_max = max(positions)
        cluster_span = cluster_max - cluster_min
        assert cluster_span <= merge_distance_bp
        rows.append(
            {
                "transition_id": f"GT{cluster_bp//1000}_{i:03d}",
                "consensus_position": consensus,
                "clades_supporting": ",".join(clades),
                "n_clades": len(clades),
                "n_unique_clades": len(clades),
                "cluster_min_position": cluster_min,
                "cluster_max_position": cluster_max,
                "cluster_span_bp": cluster_span,
                "individual_positions": ";".join(f"{event['clade']}:{event['position']}" for event in retained),
                "nearest_primary_structural_breakpoint": nearest,
                "distance_to_structure_bp": abs(consensus - nearest),
            }
        )
    return rows


def sensitivity(
    by_clade: dict[str, list[dict[str, object]]],
    primary_breakpoints: list[dict[str, object]],
    windows: list[int],
    merge_distances: list[int],
    permutations: int,
    seed: int,
) -> list[dict[str, object]]:
    bp_positions = [int(row["position"]) for row in primary_breakpoints]
    start, end = usable_range(by_clade)
    rows = []
    for window in windows:
        raw = []
        for track in by_clade.values():
            raw.extend(rolling_candidates(track, window))
        for distance in merge_distances:
            merged = merge_change_points(raw, distance)
            for label in [*CLADES, "pooled"]:
                positions = positions_for(merged, PRIMARY_METHOD, label)
                result = circular_null(positions, bp_positions, start, end, permutations, seed + stable_offset(str(window), str(distance), label))
                rows.append(
                    {
                        "clade": label,
                        "window_loci": window,
                        "merge_distance_bp": distance,
                        "n_change_points": len(positions),
                        "mean_nearest_bp": result["observed"]["mean_nearest_bp"],
                        "circular_shift_pvalue": result["pvalue"],
                        "enrichment_100kb": result["enrichments"][100_000],
                        "enrichment_250kb": result["enrichments"][250_000],
                    }
                )
    return rows


def make_blind_plot(by_clade: dict[str, list[dict[str, object]]], merged: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"q1": "#222222", "q2": "#1f78b4", "q3": "#e31a1c"}
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, sharey=True)
    for ax, clade in zip(axes, CLADES):
        track = by_clade[clade]
        xs = [row["position"] / 1e6 for row in track]
        for i, q in enumerate(Q_COLS):
            ax.plot(xs, [row["q"][i] for row in track], color=colors[q], linewidth=0.45, alpha=0.55, label=q)
        for cp in [row for row in merged if row["clade"] == clade and row["method"] == PRIMARY_METHOD]:
            ax.axvline(int(cp["position"]) / 1e6, color="#6a3d9a", linewidth=1.0)
        ax.set_title(clade, loc="left", fontsize=10)
        ax.set_ylabel("QQS")
        ax.set_ylim(-0.03, 1.03)
    axes[0].legend(loc="upper right", frameon=False, fontsize=8)
    axes[-1].set_xlabel("Chromosome 4 position (Mb)")
    fig.tight_layout()
    fig.savefig(BLIND_PDF_PATH)
    fig.savefig(BLIND_PNG_PATH, dpi=220)
    plt.close(fig)


def make_structure_plot(
    by_clade: dict[str, list[dict[str, object]]], merged: list[dict[str, object]], breakpoints: list[dict[str, object]]
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"q1": "#222222", "q2": "#1f78b4", "q3": "#e31a1c"}
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, sharey=True)
    for ax, clade in zip(axes, CLADES):
        track = by_clade[clade]
        xs = [row["position"] / 1e6 for row in track]
        for i, q in enumerate(Q_COLS):
            ax.plot(xs, [row["q"][i] for row in track], color=colors[q], linewidth=0.45, alpha=0.55, label=q)
        for bp in breakpoints:
            ax.axvline(int(bp["position"]) / 1e6, color="#7f7f7f", linewidth=0.75, linestyle=":", label="frozen de novo structural breakpoint")
        for cp in [row for row in merged if row["clade"] == clade and row["method"] == PRIMARY_METHOD]:
            ax.axvline(int(cp["position"]) / 1e6, color="#6a3d9a", linewidth=1.1, label="blind genealogy change point")
        ax.set_title(clade, loc="left", fontsize=10)
        ax.set_ylabel("QQS")
        ax.set_ylim(-0.03, 1.03)
    handles, labels = axes[0].get_legend_handles_labels()
    dedup = dict(zip(labels, handles))
    axes[0].legend(dedup.values(), dedup.keys(), loc="upper right", frameon=False, fontsize=8)
    axes[-1].set_xlabel("Chromosome 4 position (Mb)")
    fig.tight_layout()
    fig.savefig(STRUCTURE_PDF_PATH)
    fig.savefig(STRUCTURE_PNG_PATH, dpi=220)
    plt.close(fig)


def make_null_plot(nulls: dict[str, list[float]], summaries: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary_by_clade = {row["clade"]: row for row in summaries}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, label in zip(axes.ravel(), [*CLADES, "pooled"]):
        values = nulls[label]
        ax.hist(values, bins=45, color="#bdbdbd", edgecolor="white")
        obs = float(summary_by_clade[label]["mean_nearest_bp"])
        ax.axvline(obs, color="#6a3d9a", linewidth=1.5)
        ax.set_title(label, loc="left", fontsize=10)
        ax.set_xlabel("Mean nearest frozen breakpoint distance (bp)")
    fig.tight_layout()
    fig.savefig(NULL_PDF_PATH)
    plt.close(fig)


def make_structural_set_sensitivity_plot(rows: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pooled = [row for row in rows if row["clade"] == "pooled" and row["genealogy_method"] == PRIMARY_METHOD]
    fig, ax = plt.subplots(figsize=(6, 4))
    labels = [row["structural_set"] for row in pooled]
    values = [float(row["mean_nearest_bp"]) for row in pooled]
    ax.bar(labels, values, color=["#6a3d9a", "#1f78b4", "#636363"])
    ax.set_ylabel("Pooled mean nearest distance (bp)")
    ax.set_xlabel("Frozen structural set")
    fig.tight_layout()
    fig.savefig(SENSITIVITY_PDF_PATH)
    plt.close(fig)


def output_columns() -> list[str]:
    return [
        "clade", "method", "change_point_id", "position", "score",
        "left_q1", "left_q2", "left_q3", "right_q1", "right_q2", "right_q3",
        "delta_q1", "delta_q2", "delta_q3", "dominant_topology_left", "dominant_topology_right",
        "cluster_min_position", "cluster_max_position", "cluster_span_bp", "n_raw_change_points",
        "contributing_raw_change_point_ids",
    ]


def enrichment_columns(include_set: bool = True, include_method: bool = True) -> list[str]:
    columns = []
    if include_set:
        columns.append("structural_set")
    columns.append("clade")
    if include_method:
        columns.append("genealogy_method")
    columns.extend(
        [
            "n_genealogy_change_points", "n_structural_breakpoints", "mean_nearest_bp", "median_nearest_bp",
            "observed_mean_nearest_bp", "null_mean_nearest_bp", "null_median_nearest_bp", "observed_to_null_ratio",
            "frac_within_50kb", "frac_within_100kb", "frac_within_250kb", "frac_within_500kb",
            "enrichment_50kb", "enrichment_100kb", "enrichment_250kb", "enrichment_500kb",
            "circular_shift_pvalue",
        ]
    )
    return columns


def run_analysis(args: argparse.Namespace) -> None:
    old_cp_rows = read_optional_tsv(CP_PATH)
    old_summary_rows = read_optional_tsv(SUMMARY_FINAL_PATH)
    existing_comparison_rows = read_optional_tsv(SINGLELINK_COMPARISON_PATH)
    by_clade = load_loci()
    structural_sets = load_frozen_structural_breakpoints()
    before_frozen_mtimes = {path: path.stat().st_mtime_ns for path in (BREAKPOINT_SETS_PATH, CONSENSUS_BREAKPOINTS_PATH, MANIFEST_PATH)}

    raw = detect_all(by_clade, args.window_loci)
    add_raw_ids(raw)
    write_tsv(RAW_CP_PATH, raw, output_columns())
    merged = merge_change_points(raw, args.merge_distance_bp)
    write_tsv(CP_PATH, merged, output_columns())
    audit_rows = cluster_audit_rows(merged, args.merge_distance_bp)
    write_tsv(
        CLUSTER_AUDIT_PATH,
        audit_rows,
        [
            "clade", "merged_change_point_id", "representative_position", "cluster_min_position",
            "cluster_max_position", "cluster_span_bp", "merge_distance_bp", "n_raw_change_points",
            "contributing_raw_change_point_ids", "representative_score", "span_valid",
        ],
    )

    set_rows = []
    primary_nulls = {}
    for set_name in STRUCTURAL_SETS:
        rows, nulls = summarize_enrichment(
            merged, PRIMARY_METHOD, set_name, structural_sets[set_name], by_clade, args.permutations, args.seed
        )
        set_rows.extend(rows)
        if set_name == PRIMARY_STRUCTURAL_SET:
            primary_nulls = nulls
            final_rows = rows
    write_tsv(SET_SUMMARY_PATH, set_rows, enrichment_columns())

    method_rows = []
    for method in GENEALOGY_METHODS:
        rows, _ = summarize_enrichment(
            merged, method, PRIMARY_STRUCTURAL_SET, structural_sets[PRIMARY_STRUCTURAL_SET], by_clade, args.permutations, args.seed
        )
        method_rows.extend(rows)
    write_tsv(METHOD_SUMMARY_PATH, method_rows, enrichment_columns())
    write_tsv(SUMMARY_FINAL_PATH, final_rows, enrichment_columns(include_set=False, include_method=False))
    comparison_rows = compare_singlelink_to_completespan(
        old_cp_rows, old_summary_rows, existing_comparison_rows, merged, final_rows
    )
    write_tsv(
        SINGLELINK_COMPARISON_PATH,
        comparison_rows,
        [
            "clade", "old_n_change_points", "new_n_change_points", "old_mean_nearest_bp",
            "new_mean_nearest_bp", "old_circular_shift_pvalue", "new_circular_shift_pvalue",
            "conclusion_changed",
        ],
    )

    support_rows = structural_breakpoint_genealogy_support(structural_sets[PRIMARY_STRUCTURAL_SET], merged, PRIMARY_METHOD)
    write_tsv(
        BOUNDARY_CENTRIC_PATH,
        support_rows,
        [
            "consensus_breakpoint_id", "position", "n_species_support", "confidence_class", "recovery_fraction",
            "nearest_Columbea_bp", "nearest_N61_bp", "nearest_N62_bp",
            "n_clades_within_100kb", "n_clades_within_250kb", "n_clades_within_500kb",
        ],
    )
    shared_rows = []
    for cluster_bp in (100_000, 250_000):
        shared_rows.extend(cluster_shared_genealogy_transitions(merged, structural_sets[PRIMARY_STRUCTURAL_SET], cluster_bp))
    write_tsv(
        SHARED_TRANSITIONS_PATH,
        shared_rows,
        [
            "transition_id", "consensus_position", "clades_supporting", "n_clades", "individual_positions",
            "cluster_min_position", "cluster_max_position", "cluster_span_bp", "n_unique_clades",
            "nearest_primary_structural_breakpoint", "distance_to_structure_bp",
        ],
    )
    sensitivity_rows = sensitivity(
        by_clade, structural_sets[PRIMARY_STRUCTURAL_SET], [25, 50, 100], [100_000, 200_000, 500_000],
        args.sensitivity_permutations, args.seed,
    )
    write_tsv(
        SENSITIVITY_PATH,
        sensitivity_rows,
        [
            "clade", "window_loci", "merge_distance_bp", "n_change_points", "mean_nearest_bp",
            "circular_shift_pvalue", "enrichment_100kb", "enrichment_250kb",
        ],
    )
    make_blind_plot(by_clade, merged)
    make_structure_plot(by_clade, merged, structural_sets[PRIMARY_STRUCTURAL_SET])
    make_null_plot(primary_nulls, final_rows)
    make_structural_set_sensitivity_plot(set_rows)
    write_final_manifest(args, structural_sets, merged, final_rows, sensitivity_rows)

    after_frozen_mtimes = {path: path.stat().st_mtime_ns for path in before_frozen_mtimes}
    if after_frozen_mtimes != before_frozen_mtimes:
        raise AssertionError("Frozen Stage-2.7c structural input files were modified")
    assert_no_old_boundary_access()
    print_report(
        merged, final_rows, set_rows, method_rows, sensitivity_rows, support_rows, shared_rows,
        structural_sets, audit_rows, comparison_rows
    )


def read_optional_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cluster_audit_rows(merged: list[dict[str, object]], merge_distance_bp: int) -> list[dict[str, object]]:
    rows = []
    for row in merged:
        span = int(row["cluster_span_bp"])
        audit_row = {
            "clade": row["clade"],
            "merged_change_point_id": row["change_point_id"],
            "representative_position": row["position"],
            "cluster_min_position": row["cluster_min_position"],
            "cluster_max_position": row["cluster_max_position"],
            "cluster_span_bp": span,
            "merge_distance_bp": merge_distance_bp,
            "n_raw_change_points": row["n_raw_change_points"],
            "contributing_raw_change_point_ids": row["contributing_raw_change_point_ids"],
            "representative_score": row["score"],
            "span_valid": span <= merge_distance_bp,
        }
        if not audit_row["span_valid"]:
            raise AssertionError(f"Invalid genealogy cluster span: {audit_row}")
        rows.append(audit_row)
    return rows


def compare_singlelink_to_completespan(
    old_cp_rows: list[dict[str, str]],
    old_summary_rows: list[dict[str, str]],
    existing_comparison_rows: list[dict[str, str]],
    merged: list[dict[str, object]],
    final_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    old_summary = {row["clade"]: row for row in old_summary_rows}
    existing_comparison = {row["clade"]: row for row in existing_comparison_rows}
    new_summary = {str(row["clade"]): row for row in final_rows}
    rows = []
    for clade in [*CLADES, "pooled"]:
        existing_row = existing_comparison.get(clade, {})
        old_n = existing_row.get("old_n_change_points") or old_summary.get(clade, {}).get("n_genealogy_change_points")
        if old_n in (None, ""):
            if clade == "pooled":
                old_n = sum(1 for row in old_cp_rows if row.get("method") == PRIMARY_METHOD and row.get("clade") in CLADES)
            else:
                old_n = sum(1 for row in old_cp_rows if row.get("method") == PRIMARY_METHOD and row.get("clade") == clade)
        old_mean = existing_row.get("old_mean_nearest_bp") or old_summary.get(clade, {}).get("mean_nearest_bp", "")
        old_p = existing_row.get("old_circular_shift_pvalue") or old_summary.get(clade, {}).get("circular_shift_pvalue", "")
        new_row = new_summary[clade]
        new_n = sum(
            1 for row in merged
            if row["method"] == PRIMARY_METHOD and (clade == "pooled" or row["clade"] == clade)
        )
        rows.append(
            {
                "clade": clade,
                "old_n_change_points": old_n,
                "new_n_change_points": new_n,
                "old_mean_nearest_bp": old_mean,
                "new_mean_nearest_bp": new_row["mean_nearest_bp"],
                "old_circular_shift_pvalue": old_p,
                "new_circular_shift_pvalue": new_row["circular_shift_pvalue"],
                "conclusion_changed": conclusion_changed(old_p, new_row["circular_shift_pvalue"]),
            }
        )
    return rows


def conclusion_changed(old_p: object, new_p: object, alpha: float = 0.05) -> str:
    try:
        return str((float(old_p) < alpha) != (float(new_p) < alpha))
    except (TypeError, ValueError):
        return "unknown"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_final_manifest(
    args: argparse.Namespace,
    structural_sets: dict[str, list[dict[str, object]]],
    merged: list[dict[str, object]],
    final_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
) -> None:
    primary_rows = [row for row in merged if row["method"] == PRIMARY_METHOD]
    final_by_clade = {str(row["clade"]): row for row in final_rows}
    manifest = {
        "merge_method": "complete_span",
        "merge_distance_bp": args.merge_distance_bp,
        "window_loci": args.window_loci,
        "permutations": args.permutations,
        "sensitivity_permutations": args.sensitivity_permutations,
        "random_seed": args.seed,
        "structural_set": PRIMARY_STRUCTURAL_SET,
        "structural_manifest_checksum_sha256": sha256_file(MANIFEST_PATH),
        "locus_table_checksum_sha256": sha256_file(LOCUS_PATH),
        "script_checksum_sha256": sha256_file(Path(__file__).resolve()),
        "final_change_point_coordinates_by_clade": {
            clade: [
                int(row["position"]) for row in primary_rows
                if row["clade"] == clade
            ]
            for clade in CLADES
        },
        "final_change_point_counts_by_clade": {
            clade: sum(1 for row in primary_rows if row["clade"] == clade)
            for clade in CLADES
        },
        "final_per_clade_pvalues": {
            clade: float(final_by_clade[clade]["circular_shift_pvalue"])
            for clade in CLADES
        },
        "pooled_pvalue": float(final_by_clade["pooled"]["circular_shift_pvalue"]),
        "null_mean_distance": float(final_by_clade["pooled"]["null_mean_nearest_bp"]),
        "sensitivity_summary": sensitivity_rows,
        "primary_structural_breakpoint_count": len(structural_sets[PRIMARY_STRUCTURAL_SET]),
    }
    FINAL_MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def assert_no_old_boundary_access() -> None:
    old = [path for path in READ_PATHS if path.name in OLD_BOUNDARY_INPUTS]
    if old:
        raise AssertionError(f"Old structural boundary inputs were accessed: {old}")


def print_report(
    merged: list[dict[str, object]],
    final_rows: list[dict[str, object]],
    set_rows: list[dict[str, object]],
    method_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    support_rows: list[dict[str, object]],
    shared_rows: list[dict[str, object]],
    structural_sets: dict[str, list[dict[str, object]]],
    audit_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
) -> None:
    max_span = max((int(row["cluster_span_bp"]) for row in audit_rows), default=0)
    pooled = next(row for row in final_rows if row["clade"] == "pooled")
    print("Complete-span merge:")
    print(f"  maximum biological cluster span={max_span} bp")
    print(f"  pooled null mean={float(pooled['null_mean_nearest_bp']):.0f}")
    print("Rolling genealogy change points:")
    for clade in CLADES:
        positions = [str(row["position"]) for row in merged if row["method"] == PRIMARY_METHOD and row["clade"] == clade]
        print(f"  {clade}: {len(positions)} ({', '.join(positions)})")
    print(f"Primary frozen structural breakpoints: {len(structural_sets[PRIMARY_STRUCTURAL_SET])}")
    print("Main rolling + primary enrichment:")
    for row in final_rows:
        print(
            f"  {row['clade']}: n={row['n_genealogy_change_points']}, mean={float(row['mean_nearest_bp']):.0f}, "
            f"median={float(row['median_nearest_bp']):.0f}, p={float(row['circular_shift_pvalue']):.4g}, "
            f"E50={float(row['enrichment_50kb']):.3g}, E100={float(row['enrichment_100kb']):.3g}, "
            f"E250={float(row['enrichment_250kb']):.3g}, E500={float(row['enrichment_500kb']):.3g}"
        )
    print("Structural-set sensitivity, pooled:")
    for row in set_rows:
        if row["clade"] == "pooled":
            print(f"  {row['structural_set']}: mean={float(row['mean_nearest_bp']):.0f}, p={float(row['circular_shift_pvalue']):.4g}")
    print("Binary-segmentation sensitivity, pooled:")
    for row in method_rows:
        if row["genealogy_method"] == "binary_segmentation" and row["clade"] == "pooled":
            print(f"  mean={float(row['mean_nearest_bp']):.0f}, p={float(row['circular_shift_pvalue']):.4g}")
    print("Window/merge sensitivity, pooled:")
    for row in sensitivity_rows:
        if row["clade"] == "pooled":
            print(
                f"  window={row['window_loci']} merge={row['merge_distance_bp']}: "
                f"n={row['n_change_points']}, mean={float(row['mean_nearest_bp']):.0f}, "
                f"p={float(row['circular_shift_pvalue']):.4g}, E100={float(row['enrichment_100kb']):.3g}, "
                f"E250={float(row['enrichment_250kb']):.3g}"
            )
    supported = [row for row in support_rows if int(row["n_clades_within_250kb"]) >= 2]
    print("Primary structural breakpoints within 250 kb of >=2 genealogy tracks:")
    print("  " + (", ".join(f"{row['consensus_breakpoint_id']}@{row['position']}" for row in supported) if supported else "none"))
    shared_multi = [row for row in shared_rows if int(row["n_clades"]) >= 2]
    print("Shared cross-clade genealogy transitions:")
    print("  " + (", ".join(f"{row['transition_id']}@{row['consensus_position']}" for row in shared_multi) if shared_multi else "none"))
    far = []
    primary_positions = [int(row["position"]) for row in structural_sets[PRIMARY_STRUCTURAL_SET]]
    for row in merged:
        if row["method"] == PRIMARY_METHOD:
            distance = nearest_distances([int(row["position"])], primary_positions)[0]
            if distance > 500_000:
                far.append(f"{row['clade']}:{row['position']}({distance})")
    print("Rolling genealogy transitions >500 kb from a primary frozen structural breakpoint:")
    print("  " + (", ".join(far) if far else "none"))
    print("Confirmed: no old canonical structural boundaries, structural intervals, or published outlier coordinates were used.")
    print("Old single-link vs corrected complete-span:")
    for row in comparison_rows:
        print(
            f"  {row['clade']}: old_n={row['old_n_change_points']}, new_n={row['new_n_change_points']}, "
            f"old_mean={row['old_mean_nearest_bp']}, new_mean={float(row['new_mean_nearest_bp']):.0f}, "
            f"old_p={row['old_circular_shift_pvalue']}, new_p={float(row['new_circular_shift_pvalue']):.4g}, "
            f"conclusion_changed={row['conclusion_changed']}"
        )
    print("Confirmed: complete-span clustering was used for genealogy and shared-transition clusters.")
    print("Files written:")
    for path in (
        RAW_CP_PATH, CP_PATH, SET_SUMMARY_PATH, METHOD_SUMMARY_PATH, SENSITIVITY_PATH, BOUNDARY_CENTRIC_PATH,
        SHARED_TRANSITIONS_PATH, SUMMARY_FINAL_PATH, CLUSTER_AUDIT_PATH, SINGLELINK_COMPARISON_PATH,
        FINAL_MANIFEST_PATH, BLIND_PDF_PATH, BLIND_PNG_PATH, STRUCTURE_PDF_PATH, STRUCTURE_PNG_PATH,
        NULL_PDF_PATH, SENSITIVITY_PDF_PATH,
    ):
        print(f"  {path.relative_to(REPO_ROOT)}")


class BreakpointTests(unittest.TestCase):
    def synthetic_track(self) -> list[dict[str, object]]:
        rows = []
        for i in range(120):
            q = (0.8, 0.1, 0.1) if i < 60 else (0.1, 0.8, 0.1)
            rows.append({"clade": "test", "position": i * 10_000, "q": q})
        return rows

    def test_rolling_recovers_transition(self) -> None:
        cps = rolling_candidates(self.synthetic_track(), 10)
        self.assertTrue(any(abs(int(cp["position"]) - 600_000) <= 50_000 for cp in cps))

    def test_binary_recovers_transition_but_is_not_primary(self) -> None:
        cps = binary_segmentation_candidates(self.synthetic_track(), 10, max_changes=3)
        self.assertTrue(any(abs(int(cp["position"]) - 600_000) <= 50_000 for cp in cps))
        self.assertEqual(PRIMARY_METHOD, "rolling")

    def test_complete_span_merge_blocks_single_link_chaining(self) -> None:
        rows = [
            {
                "clade": "test",
                "method": "rolling",
                "change_point_id": f"test_rolling_raw_{i:03d}",
                "position": position,
                "score": score,
            }
            for i, (position, score) in enumerate([(0, 1.0), (150_000, 3.0), (300_000, 2.0)], 1)
        ]
        merged = merge_change_points(rows, 200_000)
        self.assertEqual(len(merged), 2)
        observed_clusters = [
            [int(raw_id.rsplit("_", 1)[1]) for raw_id in row["contributing_raw_change_point_ids"].split(";")]
            for row in merged
        ]
        self.assertEqual(observed_clusters, [[1, 2], [3]])
        for row in merged:
            cluster_positions = [
                int(rows[index - 1]["position"])
                for index in [int(raw_id.rsplit("_", 1)[1]) for raw_id in row["contributing_raw_change_point_ids"].split(";")]
            ]
            self.assertLessEqual(max(cluster_positions) - min(cluster_positions), 200_000)

    def test_detection_does_not_read_structural_files(self) -> None:
        def blocked_open(self_path, *args, **kwargs):
            if Path(self_path).name in {BREAKPOINT_SETS_PATH.name, CONSENSUS_BREAKPOINTS_PATH.name, MANIFEST_PATH.name, *OLD_BOUNDARY_INPUTS}:
                raise AssertionError(f"structural file read during genealogy detection: {self_path}")
            return original_open(self_path, *args, **kwargs)

        original_open = Path.open
        with mock.patch.object(Path, "open", blocked_open):
            self.assertTrue(detect_all({"test": self.synthetic_track()}, 10))

    def test_stage3a_old_boundary_guard(self) -> None:
        for filename in OLD_BOUNDARY_INPUTS:
            with self.assertRaises(AssertionError):
                read_tsv(RESULTS_DIR / filename)

    def test_stage3a_structural_inputs_are_frozen_denovo_only(self) -> None:
        READ_PATHS.clear()
        frozen = load_frozen_structural_breakpoints()
        manifest = read_json(MANIFEST_PATH)
        self.assertEqual({path.name for path in READ_PATHS}, {BREAKPOINT_SETS_PATH.name, CONSENSUS_BREAKPOINTS_PATH.name, MANIFEST_PATH.name})
        self.assertEqual(len(frozen["primary"]), int(manifest["primary_breakpoint_count"]))
        self.assertEqual(len(frozen["stringent"]), int(manifest["stringent_breakpoint_count"]))
        self.assertEqual(len(frozen["inclusive"]), int(manifest["inclusive_breakpoint_count"]))

    def test_membership_respected_exactly(self) -> None:
        rows = read_tsv(BREAKPOINT_SETS_PATH)
        frozen = load_frozen_structural_breakpoints()
        for name in STRUCTURAL_SETS:
            expected = {row["consensus_breakpoint_id"] for row in rows if truthy(row[f"in_{name}"])}
            observed = {row["consensus_breakpoint_id"] for row in frozen[name]}
            self.assertEqual(observed, expected)

    def test_circular_shift_spacing(self) -> None:
        shifted = circular_shift_positions([10, 30, 70], 15, 0, 100)
        self.assertEqual(sorted((b - a) % 100 for a, b in zip(shifted, shifted[1:] + shifted[:1])), [20, 40, 40])

    def test_nearest_distances(self) -> None:
        self.assertEqual(nearest_distances([10, 80], [0, 100]), [10, 20])

    def test_pooled_statistics_are_deterministic(self) -> None:
        a = circular_null([10, 50, 90], [0, 100], 0, 100, 100, 7)
        b = circular_null([10, 50, 90], [0, 100], 0, 100, 100, 7)
        self.assertEqual(a["pvalue"], b["pvalue"])
        self.assertEqual(a["observed"], b["observed"])

    def test_main_chr4_only(self) -> None:
        for row in read_tsv(LOCUS_PATH):
            self.assertEqual(row["Chromosome"], "chr4")

    def test_no_published_outlier_coordinates_used(self) -> None:
        self.assertFalse(any(path.name.startswith("outlier-regions-by-maf2synteny") for path in READ_PATHS))


def run_tests() -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(BreakpointTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window-loci", type=int, default=50)
    parser.add_argument("--merge-distance-bp", type=int, default=200_000)
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--sensitivity-permutations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        run_tests()
        return
    run_analysis(args)


if __name__ == "__main__":
    main()
