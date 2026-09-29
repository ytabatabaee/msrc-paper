#!/usr/bin/env python3
"""Stage 4B spatial-null tests for Atlantic cod topology alignment.

Uses frozen Stage-2/3 definitions and Stage-4A window scores. The primary
null is exact circular shifting of the D(w) track within each linkage group.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import random
import statistics
import sys
import tempfile
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
EXPECTED_LOCAL_WINDOWS_SHA = "26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b"
STAGE4A_LOCAL_QUARTETS = DATA_ROOT / "processed" / "stage4a_local_quartets.tsv"
STAGE4A_WINDOW_SCORES = DATA_ROOT / "processed" / "stage4a_window_topology_scores.tsv"
REGION_MANIFEST = DATA_ROOT / "metadata" / "region_manifest.tsv"

WINDOW_TOPOLOGY_IDS = DATA_ROOT / "processed" / "stage4b_window_topology_ids.tsv"

CIRCULAR_NULL = EMPIRICAL_ROOT / "results" / "stage4b_circular_shift_null.tsv"
CIRCULAR_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4b_circular_shift_summary.tsv"
COMPONENT_TESTS = EMPIRICAL_ROOT / "results" / "stage4b_component_tests.tsv"
BLOCK_NULL = EMPIRICAL_ROOT / "results" / "stage4b_block_placement_null.tsv"
BLOCK_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4b_block_placement_summary.tsv"
BOUNDARY_SENSITIVITY = EMPIRICAL_ROOT / "results" / "stage4b_boundary_sensitivity.tsv"
LOPO = EMPIRICAL_ROOT / "results" / "stage4b_leave_one_population_out.tsv"
TOPOLOGY_RUN_SUMMARY = EMPIRICAL_ROOT / "results" / "stage4b_topology_run_summary.tsv"
PRIMARY_PVALUES = EMPIRICAL_ROOT / "results" / "stage4b_primary_pvalues.tsv"
GLOBAL_TEST = EMPIRICAL_ROOT / "results" / "stage4b_global_test.tsv"
MANIFEST = EMPIRICAL_ROOT / "results" / "stage4b_manifest.json"
REPORT = EMPIRICAL_ROOT / "results" / "stage4b_report.md"

FIG_DIR = EMPIRICAL_ROOT / "figures"
NULL_FIGURES = {lg: FIG_DIR / f"stage4b_{lg}_null.pdf" for lg in ("LG01", "LG02", "LG07", "LG12")}
EFFECT_FIGURE = FIG_DIR / "stage4b_effect_summary.pdf"

LGS = ("LG01", "LG02", "LG07", "LG12")
GLOBAL_SEED = 20260929
GLOBAL_N = 100000


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


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    raise ValueError(f"Expected boolean string, got {value!r}")


def read_window_scores() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(STAGE4A_WINDOW_SCORES):
        converted: dict[str, object] = dict(row)
        for key in ("start", "end", "midpoint", "n_informative_quartets", "n_resolved", "n_arrangement", "n_baseline", "n_third", "n_unresolved"):
            converted[key] = int(row[key])
        for key in ("fraction_arrangement", "fraction_baseline", "fraction_third", "D_arrangement_minus_baseline"):
            converted[key] = float(row[key])
        rows.append(converted)
    return sorted(rows, key=lambda r: (str(r["lg"]), int(r["midpoint"])))


def is_partial_boundary(row: dict[str, object]) -> bool:
    return as_bool(row["overlaps_inversion"]) and not as_bool(row["inside_inversion"])


def eligible_rows(rows: list[dict[str, object]], boundary_treatment: str = "exclude") -> list[dict[str, object]]:
    if boundary_treatment == "exclude":
        return [row for row in rows if not is_partial_boundary(row)]
    if boundary_treatment in {"inside", "outside"}:
        return list(rows)
    raise ValueError(f"Unknown boundary treatment: {boundary_treatment}")


def mask_for(rows: list[dict[str, object]], boundary_treatment: str = "exclude") -> list[bool]:
    if boundary_treatment == "inside":
        return [as_bool(row["inside_inversion"]) or is_partial_boundary(row) for row in rows]
    return [as_bool(row["inside_inversion"]) for row in rows]


def delta_for(values: list[float], mask: list[bool]) -> float:
    inside = [value for value, flag in zip(values, mask) if flag]
    outside = [value for value, flag in zip(values, mask) if not flag]
    if not inside or not outside:
        return math.nan
    return statistics.fmean(inside) - statistics.fmean(outside)


def rotate(values: list[float], shift: int) -> list[float]:
    if not values:
        return []
    shift %= len(values)
    return values[shift:] + values[:shift]


def circular_null(rows: list[dict[str, object]], value_key: str = "D_arrangement_minus_baseline", boundary_treatment: str = "exclude") -> list[dict[str, object]]:
    elig = eligible_rows(rows, boundary_treatment)
    values = [float(row[value_key]) for row in elig]
    mask = mask_for(elig, boundary_treatment)
    observed = delta_for(values, mask)
    output: list[dict[str, object]] = []
    for shift in range(len(values)):
        d = delta_for(rotate(values, shift), mask)
        output.append({"shift": shift, "delta": d, "is_observed": shift == 0, "ge_observed": d >= observed})
    return output


def summarize_null(deltas: list[float], observed: float) -> dict[str, object]:
    ge = sum(delta >= observed for delta in deltas)
    greater = sum(delta > observed for delta in deltas)
    n = len(deltas)
    return {
        "n_shifts": n,
        "observed_delta_D": observed,
        "null_mean": statistics.fmean(deltas),
        "null_sd": statistics.pstdev(deltas) if len(deltas) > 1 else 0.0,
        "null_min": min(deltas),
        "null_max": max(deltas),
        "observed_rank": greater + 1,
        "percentile": sum(delta <= observed for delta in deltas) / n,
        "p_one_sided": ge / n,
        "minimum_possible_p": 1 / n,
    }


def bh_adjust(pvalues: dict[str, float]) -> dict[str, float]:
    ordered = sorted(pvalues.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted: dict[str, float] = {}
    running = 1.0
    for rank_from_end, (lg, p) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p * m / rank)
        adjusted[lg] = min(running, 1.0)
    return adjusted


def block_null(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    elig = eligible_rows(rows, "exclude")
    values = [float(row["D_arrangement_minus_baseline"]) for row in elig]
    observed_mask = mask_for(elig, "exclude")
    block_len = sum(observed_mask)
    observed = delta_for(values, observed_mask)
    out: list[dict[str, object]] = []
    for start in range(0, len(values) - block_len + 1):
        mask = [start <= i < start + block_len for i in range(len(values))]
        d = delta_for(values, mask)
        out.append({"block_start_index": start, "block_end_index": start + block_len - 1, "delta_D": d, "is_observed_geometry": mask == observed_mask, "greater_or_equal_observed": d >= observed})
    return out


def read_local_quartets() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(STAGE4A_LOCAL_QUARTETS):
        converted: dict[str, object] = dict(row)
        for key in ("start", "end", "midpoint"):
            converted[key] = int(row[key])
        rows.append(converted)
    return rows


def loo_delta(local_rows: list[dict[str, object]], excluded_population: str, lg: str) -> tuple[int, float]:
    rows = [row for row in local_rows if row["lg"] == lg and row["informative_for_competing_predictions"] == "true" and excluded_population not in [row[f"taxon{i}"] for i in range(1, 5)]]
    remaining_quartets = len({row["quartet_id"] for row in rows})
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["window_id"])].append(row)
    window_ds: list[dict[str, object]] = []
    for window_id, group in grouped.items():
        resolved = [row for row in group if row["local_class"] != "unresolved"]
        if not resolved:
            continue
        n_arr = sum(row["local_class"] == "arrangement" for row in resolved)
        n_base = sum(row["local_class"] == "baseline" for row in resolved)
        first = group[0]
        window_ds.append({"D": n_arr / len(resolved) - n_base / len(resolved), "inside_inversion": first["inside_inversion"], "overlaps_inversion": first["overlaps_inversion"], "midpoint": first["midpoint"]})
    elig = eligible_rows(window_ds, "exclude")
    values = [float(row["D"]) for row in sorted(elig, key=lambda row: int(row["midpoint"]))]
    mask = mask_for(sorted(elig, key=lambda row: int(row["midpoint"])), "exclude")
    return remaining_quartets, delta_for(values, mask)


def parse_tree(newick: str):
    tree = Phylo.read(io.StringIO(newick), "newick")
    labels = [terminal.name for terminal in tree.get_terminals()]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate taxon labels in tree")
    return tree


def canonical_topology(tree) -> str:
    taxa = frozenset(terminal.name for terminal in tree.get_terminals())
    splits: set[str] = set()
    for clade in tree.find_clades():
        side = frozenset(term.name for term in clade.get_terminals())
        if len(side) <= 1 or len(side) >= len(taxa) - 1:
            continue
        other = taxa - side
        left = ",".join(sorted(side))
        right = ",".join(sorted(other))
        splits.add("|".join(sorted([left, right])))
    return ";".join(sorted(splits))


def topology_ids() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    windows = read_tsv(LOCAL_WINDOWS)
    topo_strings: dict[str, str] = {}
    temp_rows: list[dict[str, object]] = []
    for row in sorted(windows, key=lambda r: (r["lg"], int(r["midpoint"]))):
        if row["lg"] not in LGS:
            continue
        topo = canonical_topology(parse_tree(row["tree_newick"]))
        topo_strings[row["window_id"]] = topo
        temp_rows.append({"lg": row["lg"], "window_id": row["window_id"], "midpoint": int(row["midpoint"]), "inside_inversion": row["inside_inversion"], "overlaps_inversion": row["overlaps_inversion"], "topology_string": topo})
    ids = {topo: f"topology_{i:04d}" for i, topo in enumerate(sorted(set(topo_strings.values())), start=1)}
    out_rows = [{"lg": r["lg"], "window_id": r["window_id"], "midpoint": r["midpoint"], "inside_inversion": r["inside_inversion"], "topology_id": ids[str(r["topology_string"])]} for r in temp_rows]
    summary: list[dict[str, object]] = []
    for lg in LGS:
        lg_rows = [r for r in temp_rows if r["lg"] == lg]
        inside = [ids[str(r["topology_string"])] for r in lg_rows if as_bool(r["inside_inversion"])]
        outside = [ids[str(r["topology_string"])] for r in lg_rows if not as_bool(r["inside_inversion"]) and not is_partial_boundary(r)]
        all_ids = [ids[str(r["topology_string"])] for r in lg_rows]
        counts = Counter(all_ids)
        most_common, freq = counts.most_common(1)[0]
        max_run = 0
        transitions = 0
        prev = None
        current = 0
        for topo_id in all_ids:
            if topo_id == prev:
                current += 1
            else:
                if prev is not None:
                    transitions += 1
                current = 1
                prev = topo_id
            max_run = max(max_run, current)
        summary.append({"lg": lg, "n_unique_topologies_inside": len(set(inside)), "n_unique_topologies_outside": len(set(outside)), "most_common_topology_id": most_common, "most_common_topology_frequency": freq, "maximum_consecutive_run_length": max_run, "number_of_topology_transitions": transitions})
    return out_rows, summary


def read_regions() -> dict[str, dict[str, object]]:
    out = {}
    for row in read_tsv(REGION_MANIFEST):
        if row["lg"] in LGS:
            out[row["lg"]] = {"start": int(row["start"]), "end": int(row["end"]), "region_name": row["region_name"], "coordinate_system": row["coordinate_system"]}
    return out


def global_test(primary_null_by_lg: dict[str, list[float]], observed_by_lg: dict[str, float]) -> dict[str, object]:
    observed = statistics.fmean(observed_by_lg[lg] for lg in LGS)
    rng = random.Random(GLOBAL_SEED)
    samples = []
    for _ in range(GLOBAL_N):
        samples.append(statistics.fmean(rng.choice(primary_null_by_lg[lg]) for lg in LGS))
    p = sum(value >= observed for value in samples) / len(samples)
    return {"statistic": "equal_weight_mean_delta_D", "observed_global_delta_D": observed, "n_monte_carlo": GLOBAL_N, "seed": GLOBAL_SEED, "null_mean": statistics.fmean(samples), "null_min": min(samples), "null_max": max(samples), "p_one_sided": p}


def plot_nulls(circular_summary: list[dict[str, object]], null_rows: list[dict[str, object]]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for lg in LGS:
        deltas = [float(row["delta_D"]) for row in null_rows if row["lg"] == lg]
        summary = next(row for row in circular_summary if row["lg"] == lg)
        obs = float(summary["observed_delta_D"])
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        ax.hist(deltas, bins=min(30, max(5, len(set(deltas)))), color="#bdbdbd", edgecolor="white")
        ax.axvline(obs, color="#d62728", linewidth=2, label=f"observed ΔD={obs:.3f}")
        ax.set_title(f"{lg} circular-shift null")
        ax.set_xlabel("ΔD")
        ax.set_ylabel("number of shifts")
        ax.text(0.02, 0.95, f"p={float(summary['p_one_sided']):.4f}\nn={summary['n_shifts']}", transform=ax.transAxes, va="top")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(NULL_FIGURES[lg])
        plt.close(fig)


def plot_effect_summary(circular_summary: list[dict[str, object]]) -> None:
    xs = list(range(len(LGS)))
    obs = [float(next(row for row in circular_summary if row["lg"] == lg)["observed_delta_D"]) for lg in LGS]
    mins = [float(next(row for row in circular_summary if row["lg"] == lg)["null_min"]) for lg in LGS]
    maxs = [float(next(row for row in circular_summary if row["lg"] == lg)["null_max"]) for lg in LGS]
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    for i, lg in enumerate(LGS):
        ax.plot([i, i], [mins[i], maxs[i]], color="#636363", linewidth=3, alpha=0.7)
        ax.scatter([i], [obs[i]], color="#d62728", zorder=3)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(LGS)
    ax.set_ylabel("ΔD")
    ax.set_title("Stage 4B observed effects and circular-shift null ranges")
    fig.tight_layout()
    fig.savefig(EFFECT_FIGURE)
    plt.close(fig)


def main() -> None:
    stage2_sha = verify_checksum(STAGE2_CANDIDATES, EXPECTED_STAGE2_SHA, "Stage-2 structural candidates")
    stage3_sha = verify_checksum(STAGE3_QUARTETS, EXPECTED_STAGE3_SHA, "Stage-3 baseline quartets")
    local_windows_sha = verify_checksum(LOCAL_WINDOWS, EXPECTED_LOCAL_WINDOWS_SHA, "Stage-1 local windows")
    stage4a_local_sha = sha256(STAGE4A_LOCAL_QUARTETS)
    stage4a_scores_sha = sha256(STAGE4A_WINDOW_SCORES)
    scores = read_window_scores()
    scores_by_lg = {lg: [row for row in scores if row["lg"] == lg] for lg in LGS}

    null_rows = []
    circular_summary = []
    primary_null_by_lg = {}
    observed_by_lg = {}
    for lg in LGS:
        null = circular_null(scores_by_lg[lg], "D_arrangement_minus_baseline", "exclude")
        observed = float(null[0]["delta"])
        primary_null_by_lg[lg] = [float(row["delta"]) for row in null]
        observed_by_lg[lg] = observed
        for row in null:
            null_rows.append({"lg": lg, "shift": row["shift"], "delta_D": row["delta"], "is_observed": row["is_observed"], "greater_or_equal_observed": row["ge_observed"]})
        summary = summarize_null(primary_null_by_lg[lg], observed)
        summary.update({"lg": lg, "n_windows": len(eligible_rows(scores_by_lg[lg], "exclude"))})
        circular_summary.append(summary)

    p_raw = {row["lg"]: float(row["p_one_sided"]) for row in circular_summary}
    p_bh = bh_adjust(p_raw)
    primary_rows = [{"lg": lg, "p_raw": p_raw[lg], "p_BH": p_bh[lg]} for lg in LGS]

    component_rows = []
    for lg in LGS:
        for component, key in (("delta_A", "fraction_arrangement"), ("delta_B", "fraction_baseline")):
            null = circular_null(scores_by_lg[lg], key, "exclude")
            observed = float(null[0]["delta"])
            summary = summarize_null([float(row["delta"]) for row in null], observed)
            component_rows.append({"lg": lg, "component": component, **summary})

    block_rows = []
    block_summary_rows = []
    for lg in LGS:
        rows = block_null(scores_by_lg[lg])
        observed = observed_by_lg[lg]
        for row in rows:
            block_rows.append({"lg": lg, **row})
        summary = summarize_null([float(row["delta_D"]) for row in rows], observed)
        block_summary_rows.append({"lg": lg, "null_type": "contiguous_block", **summary})

    boundary_rows = []
    for lg in LGS:
        for treatment in ("exclude", "inside", "outside"):
            null = circular_null(scores_by_lg[lg], "D_arrangement_minus_baseline", treatment)
            observed = float(null[0]["delta"])
            summary = summarize_null([float(row["delta"]) for row in null], observed)
            boundary_rows.append({"lg": lg, "boundary_treatment": treatment, "observed_delta_D": observed, "p_one_sided": summary["p_one_sided"], "rank": summary["observed_rank"], "n_shifts": summary["n_shifts"]})

    local_quartets = read_local_quartets()
    populations = sorted({row[f"taxon{i}"] for row in local_quartets for i in range(1, 5)})
    loo_rows = []
    for lg in LGS:
        full = observed_by_lg[lg]
        for pop in populations:
            n_remaining, d = loo_delta(local_quartets, str(pop), lg)
            loo_rows.append({"lg": lg, "excluded_population": pop, "n_remaining_quartets": n_remaining, "delta_D": d, "difference_from_full": d - full})

    topo_rows, topo_summary = topology_ids()
    global_row = global_test(primary_null_by_lg, observed_by_lg)

    write_tsv(CIRCULAR_NULL, null_rows, ["lg", "shift", "delta_D", "is_observed", "greater_or_equal_observed"])
    write_tsv(CIRCULAR_SUMMARY, circular_summary, ["lg", "n_windows", "n_shifts", "observed_delta_D", "null_mean", "null_sd", "null_min", "null_max", "observed_rank", "percentile", "p_one_sided", "minimum_possible_p"])
    write_tsv(PRIMARY_PVALUES, primary_rows, ["lg", "p_raw", "p_BH"])
    write_tsv(COMPONENT_TESTS, component_rows, ["lg", "component", "n_shifts", "observed_delta_D", "null_mean", "null_sd", "null_min", "null_max", "observed_rank", "percentile", "p_one_sided", "minimum_possible_p"])
    write_tsv(BLOCK_NULL, block_rows, ["lg", "block_start_index", "block_end_index", "delta_D", "is_observed_geometry", "greater_or_equal_observed"])
    write_tsv(BLOCK_SUMMARY, block_summary_rows, ["lg", "null_type", "n_shifts", "observed_delta_D", "null_mean", "null_sd", "null_min", "null_max", "observed_rank", "percentile", "p_one_sided", "minimum_possible_p"])
    write_tsv(BOUNDARY_SENSITIVITY, boundary_rows, ["lg", "boundary_treatment", "observed_delta_D", "p_one_sided", "rank", "n_shifts"])
    write_tsv(LOPO, loo_rows, ["lg", "excluded_population", "n_remaining_quartets", "delta_D", "difference_from_full"])
    write_tsv(WINDOW_TOPOLOGY_IDS, topo_rows, ["lg", "window_id", "midpoint", "inside_inversion", "topology_id"])
    write_tsv(TOPOLOGY_RUN_SUMMARY, topo_summary, ["lg", "n_unique_topologies_inside", "n_unique_topologies_outside", "most_common_topology_id", "most_common_topology_frequency", "maximum_consecutive_run_length", "number_of_topology_transitions"])
    write_tsv(GLOBAL_TEST, [global_row], ["statistic", "observed_global_delta_D", "n_monte_carlo", "seed", "null_mean", "null_min", "null_max", "p_one_sided"])
    plot_nulls(circular_summary, null_rows)
    plot_effect_summary(circular_summary)

    output_paths = [CIRCULAR_NULL, CIRCULAR_SUMMARY, COMPONENT_TESTS, BLOCK_NULL, BLOCK_SUMMARY, BOUNDARY_SENSITIVITY, LOPO, WINDOW_TOPOLOGY_IDS, TOPOLOGY_RUN_SUMMARY, PRIMARY_PVALUES, GLOBAL_TEST, *NULL_FIGURES.values(), EFFECT_FIGURE]
    write_report(circular_summary, component_rows, block_summary_rows, boundary_rows, loo_rows, topo_summary, primary_rows, global_row)
    output_paths.append(REPORT)
    manifest = {
        "stage": "4B",
        "stage2_checksum": stage2_sha,
        "stage3_checksum": stage3_sha,
        "local_windows_checksum": local_windows_sha,
        "stage4a_local_quartets_checksum": stage4a_local_sha,
        "stage4a_window_scores_checksum": stage4a_scores_sha,
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "number_windows_per_lg": {lg: len(scores_by_lg[lg]) for lg in LGS},
        "number_exact_circular_shifts": {row["lg"]: int(row["n_shifts"]) for row in circular_summary},
        "primary_statistic": "Delta_D = mean(D_inside) - mean(D_outside), D = fraction_arrangement - fraction_baseline",
        "primary_null": "Exact circular shifts of ordered D vector with fixed inversion mask, per LG",
        "boundary_rule": "Primary excludes partial boundary-overlap windows",
        "global_test_rng_seed": GLOBAL_SEED,
        "global_test_n_monte_carlo": GLOBAL_N,
        "inversion_coordinates_used": read_regions(),
        "output_checksums": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in output_paths if path.exists()},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_report(circular_summary, component_rows, block_summary_rows, boundary_rows, loo_rows, topo_summary, primary_rows, global_row) -> None:
    primary_lines = []
    for row in circular_summary:
        lg = row["lg"]
        p_bh = next(p["p_BH"] for p in primary_rows if p["lg"] == lg)
        primary_lines.append(f"| {lg} | {row['observed_delta_D']:.6f} | {row['p_one_sided']:.6f} | {p_bh:.6f} | {row['observed_rank']} | {row['null_min']:.6f} to {row['null_max']:.6f} | {row['minimum_possible_p']:.6f} |")
    comp_lines = []
    for lg in LGS:
        da = next(row for row in component_rows if row["lg"] == lg and row["component"] == "delta_A")
        db = next(row for row in component_rows if row["lg"] == lg and row["component"] == "delta_B")
        comp_lines.append(f"| {lg} | {da['observed_delta_D']:.6f} | {da['p_one_sided']:.6f} | {db['observed_delta_D']:.6f} | {db['p_one_sided']:.6f} |")
    block_lines = [f"| {row['lg']} | {row['observed_delta_D']:.6f} | {row['p_one_sided']:.6f} | {row['observed_rank']} | {row['null_min']:.6f} to {row['null_max']:.6f} |" for row in block_summary_rows]
    boundary_lines = [f"| {row['lg']} | {row['boundary_treatment']} | {row['observed_delta_D']:.6f} | {row['p_one_sided']:.6f} | {row['rank']} |" for row in boundary_rows]
    loo_lines = []
    for lg in LGS:
        vals = [float(row["delta_D"]) for row in loo_rows if row["lg"] == lg and not math.isnan(float(row["delta_D"]))]
        diffs = [float(row["difference_from_full"]) for row in loo_rows if row["lg"] == lg and not math.isnan(float(row["difference_from_full"]))]
        undefined = sum(1 for row in loo_rows if row["lg"] == lg and math.isnan(float(row["delta_D"])))
        loo_lines.append(f"| {lg} | {min(vals):.6f} to {max(vals):.6f} | {min(diffs):.6f} to {max(diffs):.6f} | {'yes' if vals and all(v > 0 for v in vals) else 'no'} | {undefined} |")
    topo_lines = [f"| {row['lg']} | {row['n_unique_topologies_inside']} | {row['n_unique_topologies_outside']} | {row['most_common_topology_id']} | {row['most_common_topology_frequency']} | {row['maximum_consecutive_run_length']} | {row['number_of_topology_transitions']} |" for row in topo_summary]
    lines = [
        "# Atlantic cod Stage-4B report",
        "",
        "Stage 4B tests whether the Stage-4A arrangement-concordance signal is unusually aligned with the independently frozen inversion intervals under spatial null models. It does not fit MSRC parameters or alter frozen definitions.",
        "",
        "## Primary spatial test",
        "",
        "The primary statistic is `Delta_D = mean(D_inside) - mean(D_outside)`, using fully inside versus fully outside windows and excluding partial boundary-overlap windows. The primary null exactly circularly shifts the ordered D track while keeping the inversion mask fixed.",
        "",
        "| LG | observed Delta_D | p_raw | p_BH | rank | null range | minimum possible p |",
        "|---|---:|---:|---:|---:|---|---:|",
        *primary_lines,
        "",
        "## Effect decomposition",
        "",
        "| LG | Delta_A | p_A | Delta_B | p_B |",
        "|---|---:|---:|---:|---:|",
        *comp_lines,
        "",
        "## Secondary block-placement null",
        "",
        "| LG | observed Delta_D | p_block | rank | null range |",
        "|---|---:|---:|---:|---|",
        *block_lines,
        "",
        "## Boundary sensitivity",
        "",
        "| LG | boundary treatment | observed Delta_D | p | rank |",
        "|---|---|---:|---:|---:|",
        *boundary_lines,
        "",
        "## Leave-one-population-out robustness",
        "",
        "| LG | finite Delta_D range | finite difference from full range | finite sign always positive | undefined exclusions |",
        "|---|---:|---:|---|---:|",
        *loo_lines,
        "",
        "## Spatial topology structure",
        "",
        "| LG | unique topologies inside | unique topologies outside | most common topology | frequency | max run | transitions |",
        "|---|---:|---:|---|---:|---:|---:|",
        *topo_lines,
        "",
        "## LG07 contrast",
        "",
        "LG07 remains fully reported and shows no corresponding positive spatial association under the primary circular-shift statistic. Stage 4B does not choose among biological explanations for that contrast.",
        "",
        "## Global secondary test",
        "",
        f"The equal-weight four-LG statistic was {global_row['observed_global_delta_D']:.6f}; the seeded Monte Carlo circular-shift summary used seed {global_row['seed']} with {global_row['n_monte_carlo']} draws and gave p = {global_row['p_one_sided']:.6f}. This is secondary to the per-LG tests.",
        "",
        "## Limitations",
        "",
        "- Only four genomic supergene regions are tested.",
        "- Local windows are linked along chromosomes.",
        "- Quartet observations within a tree are highly dependent.",
        "- This is multi-population, not a direct multi-species persistence test.",
        "- Published MCC trees summarize posterior genealogy uncertainty.",
        "- Spatial association alone does not establish a causal inversion effect.",
        "",
        "## Interpretation",
        "",
        "Where supported, local genealogies within the inversion are unusually aligned with the arrangement partition relative to other positions on the same linkage group. This is not a claim that MSRC is proven or that the inversions caused the genealogies.",
        "",
    ]
    REPORT.write_text("\n".join(lines))


class Stage4BTests(unittest.TestCase):
    def test_circular_shift_zero_reproduces_observed(self):
        rows = [{"D_arrangement_minus_baseline": x, "inside_inversion": str(m).lower(), "overlaps_inversion": str(m).lower()} for x, m in [(1, True), (2, False), (3, False)]]
        null = circular_null(rows)
        self.assertEqual(null[0]["delta"], delta_for([1, 2, 3], [True, False, False]))

    def test_circular_shift_preserves_values_and_size(self):
        values = [1, 2, 3, 4]
        self.assertEqual(sorted(rotate(values, 2)), values)
        self.assertEqual(len(rotate(values, 2)), 4)

    def test_exact_p_one_sided(self):
        summary = summarize_null([0, 1, 2, 3], 2)
        self.assertEqual(summary["p_one_sided"], 0.5)
        self.assertEqual(summary["observed_rank"], 2)

    def test_one_sided_direction(self):
        self.assertEqual(summarize_null([0, 1, 2], 2)["p_one_sided"], 1/3)

    def test_contiguous_block_enumeration(self):
        rows = [{"D_arrangement_minus_baseline": x, "inside_inversion": str(i < 2).lower(), "overlaps_inversion": str(i < 2).lower()} for i, x in enumerate([1, 1, 0, 0])]
        self.assertEqual(len(block_null(rows)), 3)

    def test_boundary_exclusion(self):
        rows = [{"inside_inversion": "false", "overlaps_inversion": "true"}, {"inside_inversion": "true", "overlaps_inversion": "true"}]
        self.assertEqual(len(eligible_rows(rows, "exclude")), 1)

    def test_leave_one_population_out_recalculation(self):
        rows = []
        for w in ["W1", "W2"]:
            for q, cls in [("Q1", "arrangement"), ("Q2", "baseline")]:
                rows.append({"lg": "LGX", "window_id": w, "quartet_id": q, "informative_for_competing_predictions": "true", "local_class": cls, "taxon1": "A", "taxon2": "B" if q == "Q1" else "C", "taxon3": "D", "taxon4": "E", "inside_inversion": "true" if w == "W1" else "false", "overlaps_inversion": "true" if w == "W1" else "false", "midpoint": 1 if w == "W1" else 2})
        n, d = loo_delta(rows, "C", "LGX")
        self.assertEqual(n, 1)
        self.assertEqual(d, 0)

    def test_canonical_topology_ignores_branch_lengths(self):
        t1 = canonical_topology(parse_tree("((A:1,B:2):3,(C:4,D:5):6);"))
        t2 = canonical_topology(parse_tree("((A:9,B:8):7,(C:6,D:5):4);"))
        self.assertEqual(t1, t2)

    def test_deterministic_ordering(self):
        self.assertEqual(bh_adjust({"b": 0.02, "a": 0.01})["a"], 0.02)

    def test_checksum_verification(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x"
            path.write_text("abc\n")
            expected = hashlib.sha256(b"abc\n").hexdigest()
            self.assertEqual(verify_checksum(path, expected, "x"), expected)
            with self.assertRaises(ValueError):
                verify_checksum(path, "0" * 64, "x")

    def test_seeded_global_monte_carlo_reproducibility(self):
        nulls = {lg: [0.0, 1.0] for lg in LGS}
        obs = {lg: 0.5 for lg in LGS}
        one = global_test(nulls, obs)
        two = global_test(nulls, obs)
        self.assertEqual(one["p_one_sided"], two["p_one_sided"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="Run self-tests and exit.")
    args = parser.parse_args()
    if args.run_tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Stage4BTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        sys.exit(0 if result.wasSuccessful() else 1)
    main()
