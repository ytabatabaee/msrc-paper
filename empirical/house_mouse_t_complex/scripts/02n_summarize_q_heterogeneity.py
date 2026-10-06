#!/usr/bin/env python3
"""Summarize frozen window heterogeneity and arrangement-mixture contributions."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import RESULTS, read_tsv, write_tsv  # noqa: E402

TOPOLOGIES = ("q_species", "q_t_alt", "q_other")
PATTERNS = ("SSS", "SST", "STS", "STT", "TSS", "TST", "TTS", "TTT")


def all_tip_rows() -> list[dict[str, str]]:
    return [row for row in read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv") if row["treatment"] == "ALL_TIPS"]


def summarize_window_heterogeneity(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    values = {key: np.array([float(row[key]) for row in rows]) for key in TOPOLOGIES}
    delta = values["q_species"] - values["q_t_alt"]
    winner_counts = {key: 0 for key in ("Q_SPECIES", "Q_T_ALT", "Q_OTHER", "TIE_OR_UNRESOLVED")}
    for row in rows:
        q = np.array([float(row[key]) for key in TOPOLOGIES])
        if float(row["q_unresolved"]) > max(q) or np.sum(np.isclose(q, q.max(), atol=1e-12)) != 1:
            winner_counts["TIE_OR_UNRESOLVED"] += 1
        else:
            winner_counts[("Q_SPECIES", "Q_T_ALT", "Q_OTHER")[int(np.argmax(q))]] += 1
    out = []
    for winner, count in winner_counts.items():
        out.append({"section": "dominant_topology", "metric": winner, "statistic": "count", "value": count})
        out.append({"section": "dominant_topology", "metric": winner, "statistic": "fraction", "value": count / len(rows)})
    for key, array in [*values.items(), ("delta_species_alt", delta)]:
        for statistic, value in [("mean", np.mean(array)), ("median", np.median(array)), ("p10", np.percentile(array, 10)), ("p25", np.percentile(array, 25)), ("p75", np.percentile(array, 75)), ("p90", np.percentile(array, 90))]:
            out.append({"section": "distribution", "metric": key, "statistic": statistic, "value": value})
    abs_delta = np.abs(delta)
    for threshold in (0.25, 0.50, 0.75, 0.90):
        out.append({"section": "decisive_delta", "metric": "abs_delta_species_alt", "statistic": f">={threshold:.2f}", "value": np.mean(abs_delta >= threshold)})
    return out


def pattern_contributions() -> list[dict[str, object]]:
    rows = {row["status_pattern"]: row for row in read_tsv(RESULTS / "stage2_arrangement_pattern_summary.tsv") if row["status_pattern"] in PATTERNS}
    total = sum(int(rows[pattern]["n_induced_quartets"]) for pattern in PATTERNS)
    out = []
    for pattern in PATTERNS:
        row = rows[pattern]
        weight = int(row["n_induced_quartets"]) / total
        q_species, q_t_alt, q_other = (float(row[key]) for key in TOPOLOGIES)
        out.append({"status_pattern": pattern, "n_induced_quartets": int(row["n_induced_quartets"]), "pattern_weight": weight, "q_species": q_species, "q_t_alt": q_t_alt, "q_other": q_other, "contribution_q_species": weight * q_species, "contribution_q_t_alt": weight * q_t_alt, "contribution_q_other": weight * q_other, "contribution_delta_species_alt": weight * (q_species - q_t_alt)})
    totals = {key: sum(float(row[key]) for row in out) for key in ("pattern_weight", "contribution_q_species", "contribution_q_t_alt", "contribution_q_other")}
    target = next(row for row in read_tsv(RESULTS / "stage2_fixed_quartet_summary.tsv") if row["treatment"] == "ALL_TIPS")
    assert abs(totals["pattern_weight"] - 1) < 1e-12
    assert abs(totals["contribution_q_species"] - float(target["mean_q_species"])) < 1e-8
    assert abs(totals["contribution_q_t_alt"] - float(target["mean_q_t_alt"])) < 1e-8
    assert abs(totals["contribution_q_other"] - float(target["mean_q_other"])) < 1e-8
    return out


def smoothing_sensitivity(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out = []
    for bin_size in (5_000, 25_000, 50_000, 100_000, 250_000, 500_000):
        groups: dict[int, list[dict[str, str]]] = {}
        for row in rows:
            start = int(row["start_bp"])
            bin_start = 5_000_000 + ((start - 5_000_000) // bin_size) * bin_size
            groups.setdefault(bin_start, []).append(row)
        q = {key: np.array([np.mean([float(row[key]) for row in block]) for _, block in sorted(groups.items())]) for key in TOPOLOGIES}
        delta = q["q_species"] - q["q_t_alt"]
        maximum = np.maximum.reduce([q[key] for key in TOPOLOGIES])
        record = {"bin_size_bp": bin_size, "n_bins": len(groups), "mean_q_species": np.mean(q["q_species"]), "mean_q_t_alt": np.mean(q["q_t_alt"]), "mean_q_other": np.mean(q["q_other"]), "variance_q_species": np.var(q["q_species"]), "variance_q_t_alt": np.var(q["q_t_alt"]), "variance_q_other": np.var(q["q_other"]), "mean_abs_delta_species_alt": np.mean(np.abs(delta))}
        for key, array in q.items():
            record[f"fraction_{key[2:]}_gt_0.5"] = np.mean(array > 0.5)
        for threshold in (0.6, 0.7, 0.8, 0.9):
            record[f"fraction_max_ge_{threshold:.1f}"] = np.mean(maximum >= threshold)
        out.append(record)
    return out


def main() -> int:
    rows = all_tip_rows()
    write_tsv(RESULTS / "stage2_window_topology_heterogeneity.tsv", summarize_window_heterogeneity(rows), ["section", "metric", "statistic", "value"])
    write_tsv(RESULTS / "stage2_arrangement_pattern_contributions.tsv", pattern_contributions(), ["status_pattern", "n_induced_quartets", "pattern_weight", "q_species", "q_t_alt", "q_other", "contribution_q_species", "contribution_q_t_alt", "contribution_q_other", "contribution_delta_species_alt"])
    smoothing = smoothing_sensitivity(rows)
    write_tsv(RESULTS / "stage2_q_smoothing_sensitivity.tsv", smoothing, list(smoothing[0].keys()))
    print("Wrote heterogeneity, arrangement-contribution, and q-smoothing summaries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
