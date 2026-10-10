#!/usr/bin/env python3
"""Low-dimensional, block-held-out validation of the symmetric 2:2 MSRC quartet result."""

from __future__ import annotations

import csv
import math
import random
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "empirical/house_mouse_t_complex/results"
FIGURES = ROOT / "empirical/house_mouse_t_complex/figures"
PATTERNS = ("STT", "TST", "TTS")
TOPOLOGIES = ("q_species", "q_t_alt", "q_other")
ARRANGEMENT_TOPOLOGY = {"STT": "q_species", "TST": "q_other", "TTS": "q_t_alt"}
BLOCK_SIZES = (500_000, 1_000_000, 2_000_000)
PRIMARY_BLOCK_SIZE = 1_000_000
BOOTSTRAP_SEED = 20261009
BOOTSTRAP_REPLICATES = {500_000: 5_000, 1_000_000: 10_000, 2_000_000: 5_000}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: ("NA" if row.get(field) is None else row.get(field)) for field in fields})


def f(value: float) -> str:
    return f"{value:.12g}"


def normalize(q: tuple[float, float, float]) -> tuple[float, float, float]:
    total = sum(q)
    return tuple(x / total for x in q)


def delta_from_q(q_arr: float) -> float:
    return (3 * q_arr - 1) / 2


def predicted_q(delta: float) -> tuple[float, float, float]:
    non = (1 - delta) / 3
    return ((1 + 2 * delta) / 3, non, non)


def vector_metrics(pred: tuple[float, float, float], obs: tuple[float, float, float], arrangement_index: int) -> dict[str, float]:
    errors = [abs(a - b) for a, b in zip(pred, obs)]
    non_indices = [i for i in range(3) if i != arrangement_index]
    return {
        "max_abs_error": max(errors),
        "l1_distance": sum(errors),
        "arrangement_support_error": pred[arrangement_index] - obs[arrangement_index],
        "nonarrangement_symmetry_error": abs(obs[non_indices[0]] - obs[non_indices[1]]),
        "nonarrangement_prediction_error": max(abs(pred[i] - obs[i]) for i in non_indices),
    }


def load_aggregate() -> dict[str, dict[str, float | str]]:
    rows = {r["status_pattern"]: r for r in read_tsv(RESULTS / "stage2_arrangement_pattern_summary.tsv")}
    expected = {"STT", "TST", "TTS"}
    assert expected.issubset(rows)
    out: dict[str, dict[str, float | str]] = {}
    for pattern in PATTERNS:
        raw = tuple(float(rows[pattern][name]) for name in TOPOLOGIES)
        q = normalize(raw)
        unresolved = float(rows[pattern]["q_unresolved"])
        arr_name = ARRANGEMENT_TOPOLOGY[pattern]
        arr_index = TOPOLOGIES.index(arr_name)
        out[pattern] = {
            "pattern": pattern,
            "arrangement_topology": arr_name,
            "q_species": q[0], "q_t_alt": q[1], "q_other": q[2],
            "q_arrangement": q[arr_index],
            "q_nonarr_1": q[(arr_index + 1) % 3],
            "q_nonarr_2": q[(arr_index + 2) % 3],
            "q_unresolved": unresolved,
            "delta_hat": delta_from_q(q[arr_index]),
            "nonarrangement_symmetry_error": abs(q[(arr_index + 1) % 3] - q[(arr_index + 2) % 3]),
        }
    return out


def aggregate_rows(aggregate: dict[str, dict[str, float | str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pattern in PATTERNS:
        a = aggregate[pattern]
        rows.append({
            "analysis": "aggregate",
            "configuration": pattern,
            "training_configurations": "all",
            "arrangement_topology": a["arrangement_topology"],
            "q_species_observed": a["q_species"], "q_t_alt_observed": a["q_t_alt"], "q_other_observed": a["q_other"],
            "q_arrangement_observed": a["q_arrangement"], "q_nonarr_1_observed": a["q_nonarr_1"], "q_nonarr_2_observed": a["q_nonarr_2"],
            "q_unresolved_observed": a["q_unresolved"], "delta_eff": a["delta_hat"],
            "q_species_predicted": None, "q_t_alt_predicted": None, "q_other_predicted": None,
            "max_abs_error": None, "l1_distance": None, "arrangement_support_error": None,
            "nonarrangement_symmetry_error": a["nonarrangement_symmetry_error"], "nonarrangement_prediction_error": None,
        })
    return rows


def prediction_rows(aggregate: dict[str, dict[str, float | str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for heldout in PATTERNS:
        train = [p for p in PATTERNS if p != heldout]
        delta = sum(float(aggregate[p]["delta_hat"]) for p in train) / len(train)
        pred_arr = predicted_q(delta)
        arr_index = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[heldout])
        pred = [0.0, 0.0, 0.0]
        pred[arr_index] = pred_arr[0]
        non = [i for i in range(3) if i != arr_index]
        pred[non[0]], pred[non[1]] = pred_arr[1], pred_arr[2]
        obs = tuple(float(aggregate[heldout][name]) for name in TOPOLOGIES)
        metrics = vector_metrics(tuple(pred), obs, arr_index)
        rows.append({
            "analysis": "leave_one_configuration_out", "configuration": heldout, "training_configurations": ",".join(train),
            "arrangement_topology": ARRANGEMENT_TOPOLOGY[heldout], "q_species_observed": obs[0], "q_t_alt_observed": obs[1], "q_other_observed": obs[2],
            "q_arrangement_observed": obs[arr_index], "q_nonarr_1_observed": aggregate[heldout]["q_nonarr_1"], "q_nonarr_2_observed": aggregate[heldout]["q_nonarr_2"],
            "q_unresolved_observed": aggregate[heldout]["q_unresolved"], "delta_eff": delta,
            "q_species_predicted": pred[0], "q_t_alt_predicted": pred[1], "q_other_predicted": pred[2], **metrics,
        })
    stt_delta = float(aggregate["STT"]["delta_hat"])
    for heldout in ("TST", "TTS"):
        arr_index = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[heldout])
        pred_arr = predicted_q(stt_delta)
        pred = [0.0, 0.0, 0.0]
        pred[arr_index] = pred_arr[0]
        non = [i for i in range(3) if i != arr_index]
        pred[non[0]], pred[non[1]] = pred_arr[1], pred_arr[2]
        obs = tuple(float(aggregate[heldout][name]) for name in TOPOLOGIES)
        metrics = vector_metrics(tuple(pred), obs, arr_index)
        rows.append({
            "analysis": "STT_only_prediction", "configuration": heldout, "training_configurations": "STT",
            "arrangement_topology": ARRANGEMENT_TOPOLOGY[heldout], "q_species_observed": obs[0], "q_t_alt_observed": obs[1], "q_other_observed": obs[2],
            "q_arrangement_observed": obs[arr_index], "q_nonarr_1_observed": aggregate[heldout]["q_nonarr_1"], "q_nonarr_2_observed": aggregate[heldout]["q_nonarr_2"],
            "q_unresolved_observed": aggregate[heldout]["q_unresolved"], "delta_eff": stt_delta,
            "q_species_predicted": pred[0], "q_t_alt_predicted": pred[1], "q_other_predicted": pred[2], **metrics,
        })
    return rows


def load_block_rows(block_size: int) -> dict[int, dict[str, dict[str, object]]]:
    source = read_tsv(RESULTS / "stage2_arrangement_pattern_quartets.tsv")
    blocks: dict[int, dict[str, dict[str, object]]] = defaultdict(dict)
    for r in source:
        pattern = r["status_pattern"]
        if pattern not in PATTERNS:
            continue
        start = int(r["start_bp"])
        block = (start - 5_000_000) // block_size
        bucket = blocks[block].setdefault(pattern, {"n_quartets": 0.0, "weighted": [0.0, 0.0, 0.0], "unresolved": 0.0, "windows": set(), "midpoints": []})
        n = float(r["n_quartets"])
        bucket["n_quartets"] += n
        bucket["weighted"] = [x + n * float(r[q]) for x, q in zip(bucket["weighted"], TOPOLOGIES)]
        bucket["unresolved"] += n * float(r["q_unresolved"])
        bucket["windows"].add(r["locus_id"])
        bucket["midpoints"].append(float(r["midpoint_bp"]))
    for block, config_rows in blocks.items():
        assert set(config_rows) == set(PATTERNS), (block, set(config_rows))
        for pattern, d in config_rows.items():
            q = normalize(tuple(x / d["n_quartets"] for x in d["weighted"]))
            arr_index = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[pattern])
            d.update({"q": q, "delta": delta_from_q(q[arr_index]), "q_unresolved": d["unresolved"] / d["n_quartets"], "n_windows": len(d["windows"]), "midpoint": sum(d["midpoints"]) / len(d["midpoints"])})
    return blocks


def block_cv_rows(blocks: dict[int, dict[str, dict[str, object]]], block_size: int) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    all_blocks = sorted(blocks)
    for test_block in all_blocks:
        train_blocks = [b for b in all_blocks if b != test_block]
        for heldout in PATTERNS:
            train_configs = [p for p in PATTERNS if p != heldout]
            delta_by_config = {p: sum(float(blocks[b][p]["delta"]) for b in train_blocks) / len(train_blocks) for p in train_configs}
            delta = sum(delta_by_config.values()) / len(delta_by_config)
            arr_index = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[heldout])
            pred_arr = predicted_q(delta)
            pred = [0.0, 0.0, 0.0]; pred[arr_index] = pred_arr[0]
            non = [i for i in range(3) if i != arr_index]; pred[non[0]], pred[non[1]] = pred_arr[1], pred_arr[2]
            obs = tuple(blocks[test_block][heldout]["q"])
            metrics = vector_metrics(tuple(pred), obs, arr_index)
            output.append({
                "record_type": "fold", "block_size_bp": block_size, "test_block": test_block,
                "test_midpoint_mb": blocks[test_block][heldout]["midpoint"] / 1e6,
                "heldout_configuration": heldout, "training_configurations": ",".join(train_configs),
                "n_test_windows": blocks[test_block][heldout]["n_windows"], "n_train_blocks": len(train_blocks),
                "delta_eff_train": delta, "q_species_observed": obs[0], "q_t_alt_observed": obs[1], "q_other_observed": obs[2],
                "q_species_predicted": pred[0], "q_t_alt_predicted": pred[1], "q_other_predicted": pred[2],
                "q_unresolved_observed": blocks[test_block][heldout]["q_unresolved"], **metrics,
            })
    for heldout in PATTERNS:
        subset = [r for r in output if r["heldout_configuration"] == heldout]
        summary = {"record_type": "summary", "block_size_bp": block_size, "test_block": "ALL", "test_midpoint_mb": None,
                   "heldout_configuration": heldout, "training_configurations": "other two", "n_test_windows": sum(int(r["n_test_windows"]) for r in subset),
                   "n_train_blocks": None, "n_folds": len(subset), "delta_eff_train": sum(float(r["delta_eff_train"]) for r in subset) / len(subset),
                   "q_species_observed": None, "q_t_alt_observed": None, "q_other_observed": None, "q_species_predicted": None, "q_t_alt_predicted": None, "q_other_predicted": None,
                   "q_unresolved_observed": None, "max_abs_error": sum(float(r["max_abs_error"]) for r in subset) / len(subset),
                   "l1_distance": sum(float(r["l1_distance"]) for r in subset) / len(subset),
                   "arrangement_support_error": sum(float(r["arrangement_support_error"]) for r in subset) / len(subset),
                   "nonarrangement_symmetry_error": sum(float(r["nonarrangement_symmetry_error"]) for r in subset) / len(subset),
                   "nonarrangement_prediction_error": sum(float(r["nonarrangement_prediction_error"]) for r in subset) / len(subset)}
        output.append(summary)
    return output


def bootstrap_one(blocks: dict[int, dict[str, dict[str, object]]], sampled_blocks: list[int], replicate: int, block_size: int) -> list[dict[str, object]]:
    aggregate: dict[str, tuple[float, list[float], float]] = {}
    for pattern in PATTERNS:
        n_quartets = sum(float(blocks[b][pattern]["n_quartets"]) for b in sampled_blocks)
        weighted = [sum(float(blocks[b][pattern]["weighted"][i]) for b in sampled_blocks) for i in range(3)]
        q = normalize(tuple(x / n_quartets for x in weighted))
        unresolved = sum(float(blocks[b][pattern]["unresolved"]) for b in sampled_blocks) / n_quartets
        aggregate[pattern] = (n_quartets, list(q), unresolved)
    deltas = {p: delta_from_q(aggregate[p][1][TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[p])]) for p in PATTERNS}
    rows: list[dict[str, object]] = []
    for heldout in PATTERNS:
        train = [p for p in PATTERNS if p != heldout]
        delta_train = sum(deltas[p] for p in train) / 2
        arr_index = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[heldout])
        p_arr = predicted_q(delta_train)
        pred = [0.0, 0.0, 0.0]
        pred[arr_index] = p_arr[0]
        non = [i for i in range(3) if i != arr_index]
        pred[non[0]], pred[non[1]] = p_arr[1], p_arr[2]
        obs = tuple(aggregate[heldout][1])
        metrics = vector_metrics(tuple(pred), obs, arr_index)
        rows.append({
            "block_size_bp": block_size, "bootstrap_replicate": replicate, "heldout_configuration": heldout,
            "training_configurations": ",".join(train), "n_sampled_blocks": len(sampled_blocks),
            "delta_STT": deltas["STT"], "delta_TST": deltas["TST"], "delta_TTS": deltas["TTS"], "delta_train": delta_train,
            "q_species_observed": obs[0], "q_t_alt_observed": obs[1], "q_other_observed": obs[2],
            "q_species_predicted": pred[0], "q_t_alt_predicted": pred[1], "q_other_predicted": pred[2],
            "max_abs_error": metrics["max_abs_error"], "l1_distance": metrics["l1_distance"],
            "arrangement_support_error": metrics["arrangement_support_error"],
            "nonarrangement_symmetry_error": metrics["nonarrangement_symmetry_error"],
        })
    return rows


def bootstrap_rows(blocks: dict[int, dict[str, dict[str, object]]], block_size: int, n_bootstrap: int, seed: int) -> list[dict[str, object]]:
    block_ids = sorted(blocks)
    rng = random.Random(seed + block_size)
    rows: list[dict[str, object]] = []
    for replicate in range(n_bootstrap):
        sampled = [block_ids[rng.randrange(len(block_ids))] for _ in block_ids]
        rows.extend(bootstrap_one(blocks, sampled, replicate, block_size))
    return rows


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def bootstrap_summary(rows: list[dict[str, object]], block_size: int, n_bootstrap: int) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for pattern in PATTERNS:
        subset = [r for r in rows if r["heldout_configuration"] == pattern]
        def stats(field: str) -> tuple[float, float, float]:
            values = [float(r[field]) for r in subset]
            return percentile(values, 0.5), percentile(values, 0.025), percentile(values, 0.975)
        d = stats("delta_train")
        e = stats("max_abs_error")
        l1 = stats("l1_distance")
        a = stats("arrangement_support_error")
        s = stats("nonarrangement_symmetry_error")
        output.append({
            "block_size_bp": block_size, "record_type": "heldout_configuration", "configuration": pattern, "n_bootstrap": n_bootstrap,
            "delta_train_median": d[0], "delta_train_ci025": d[1], "delta_train_ci975": d[2],
            "max_abs_error_median": e[0], "max_abs_error_ci025": e[1], "max_abs_error_ci975": e[2],
            "l1_median": l1[0], "l1_ci025": l1[1], "l1_ci975": l1[2],
            "arrangement_support_error_median": a[0], "arrangement_support_error_ci025": a[1], "arrangement_support_error_ci975": a[2],
            "nonarrangement_symmetry_error_median": s[0], "nonarrangement_symmetry_error_ci025": s[1], "nonarrangement_symmetry_error_ci975": s[2],
        })
    for pattern in PATTERNS:
        values = [float(r[f"delta_{pattern}"]) for r in rows if r["heldout_configuration"] == "STT"]
        output.append({
            "block_size_bp": block_size, "record_type": "configuration_delta", "configuration": pattern, "n_bootstrap": n_bootstrap,
            "delta_train_median": percentile(values, 0.5), "delta_train_ci025": percentile(values, 0.025), "delta_train_ci975": percentile(values, 0.975),
            "max_abs_error_median": None, "max_abs_error_ci025": None, "max_abs_error_ci975": None,
            "l1_median": None, "l1_ci025": None, "l1_ci975": None,
            "arrangement_support_error_median": None, "arrangement_support_error_ci025": None, "arrangement_support_error_ci975": None,
            "nonarrangement_symmetry_error_median": None, "nonarrangement_symmetry_error_ci025": None, "nonarrangement_symmetry_error_ci975": None,
        })
    return output


def rank_values(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2
        for k in range(i, j):
            ranks[order[k]] = rank
        i = j
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    rx, ry = rank_values(x), rank_values(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def local_heterogeneity_rows(blocks: dict[int, dict[str, dict[str, object]]], block_size: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for block in sorted(blocks):
        for pattern in PATTERNS:
            d = blocks[block][pattern]
            q = d["q"]
            rows.append({"block_size_bp": block_size, "block_id": block, "block_midpoint_mb": d["midpoint"] / 1e6, "configuration": pattern,
                         "n_windows": d["n_windows"], "q_species": q[0], "q_t_alt": q[1], "q_other": q[2],
                         "q_unresolved": d["q_unresolved"], "delta_block": d["delta"]})
    return rows


def local_heterogeneity_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    by_pattern = {p: [r for r in rows if r["configuration"] == p] for p in PATTERNS}
    stats: dict[str, object] = {}
    for p in PATTERNS:
        values = [float(r["delta_block"]) for r in by_pattern[p]]
        q1, med, q3 = percentile(values, 0.25), percentile(values, 0.5), percentile(values, 0.75)
        stats[p] = {"range": (min(values), max(values)), "median": med, "iqr": (q1, q3), "fraction_negative": sum(v < 0 for v in values) / len(values), "fraction_gt1": sum(v > 1 for v in values) / len(values)}
    common = {int(r["block_id"]): r for r in by_pattern["STT"]}
    for p in ("TST", "TTS"):
        other = {int(r["block_id"]): r for r in by_pattern[p]}
        ids = sorted(set(common) & set(other))
        stats[f"spearman_STT_{p}"] = spearman([float(common[i]["delta_block"]) for i in ids], [float(other[i]["delta_block"]) for i in ids])
    return stats


def build_figure(aggregate: dict[str, dict[str, float | str]], loo_rows: list[dict[str, object]], bootstrap_summary_rows: list[dict[str, object]]) -> tuple[Path, Path]:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = {"STT": "#0072B2", "TST": "#D55E00", "TTS": "#009E73"}
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.25), constrained_layout=True)

    ax = axes[0]
    x = list(range(3)); width = 0.24
    topo_labels = ["Q_SPECIES", "Q_T_ALT", "Q_OTHER"]
    for j, topo in enumerate(TOPOLOGIES):
        vals = [float(aggregate[p][topo]) for p in PATTERNS]
        bars = ax.bar([i + (j - 1) * width for i in x], vals, width, label=topo, color=["#0072B2", "#D55E00", "#666666"][j], alpha=0.85)
    for i, p in enumerate(PATTERNS):
        arr_idx = TOPOLOGIES.index(ARRANGEMENT_TOPOLOGY[p])
        ax.bar(i + (arr_idx - 1) * width, float(aggregate[p][TOPOLOGIES[arr_idx]]), width, fill=False, edgecolor=colors[p], linewidth=2.1)
    ax.set_xticks(x, ["STT\nQ_SPECIES", "TST\nQ_OTHER", "TTS\nQ_T_ALT"]); ax.set_ylim(0, 0.78)
    ax.set_ylabel("Resolved quartet frequency"); ax.set_title("A. Rotating arrangement-concordant topology")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02), handlelength=1)
    ax.grid(axis="y", alpha=0.18, lw=0.45)

    ax = axes[1]
    for row in [r for r in loo_rows if r["analysis"] == "leave_one_configuration_out"]:
        for topo in TOPOLOGIES:
            ax.scatter(float(row[f"q_{topo[2:]}_predicted"]), float(row[f"q_{topo[2:]}_observed"]), s=35, marker={"STT": "o", "TST": "s", "TTS": "D"}[row["configuration"]], color=colors[row["configuration"]], edgecolor="white", linewidth=0.4)
    ax.plot([0.1, 0.75], [0.1, 0.75], "--", color="0.45", lw=0.9)
    ax.set(xlim=(0.1, 0.75), ylim=(0.1, 0.75), xlabel="Predicted quartet frequency", ylabel="Observed quartet frequency", title="B. Leave-one-configuration-out prediction")
    handles = [plt.Line2D([], [], marker=m, color="none", markerfacecolor=colors[p], markeredgecolor="white", markersize=5, label=p) for p, m in (("STT", "o"), ("TST", "s"), ("TTS", "D"))]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right", ncol=3, handletextpad=0.2, columnspacing=0.5); ax.grid(alpha=0.18, lw=0.45)

    ax = axes[2]
    y = list(range(len(PATTERNS)))[::-1]
    for ypos, p in zip(y, PATTERNS):
        summary = next(r for r in bootstrap_summary_rows if r["record_type"] == "heldout_configuration" and r["configuration"] == p)
        median = float(summary["max_abs_error_median"])
        low = float(summary["max_abs_error_ci025"])
        high = float(summary["max_abs_error_ci975"])
        observed = next(float(r["max_abs_error"]) for r in loo_rows if r["analysis"] == "leave_one_configuration_out" and r["configuration"] == p)
        ax.plot([low, high], [ypos, ypos], color=colors[p], lw=3, alpha=0.65, solid_capstyle="round")
        ax.scatter(median, ypos, s=42, color=colors[p], marker="o", edgecolor="white", linewidth=0.5, zorder=3)
        ax.scatter(observed, ypos, s=48, color="black", marker="D", edgecolor="white", linewidth=0.5, zorder=4)
    ax.set(yticks=y, yticklabels=PATTERNS, xlabel="Maximum absolute prediction error", title="C. Aggregate error under block bootstrap")
    ax.set_xlim(left=0); ax.grid(axis="x", alpha=0.18, lw=0.45)
    ax.legend([plt.Line2D([], [], color="0.35", lw=3), plt.Line2D([], [], marker="o", color="none", markerfacecolor="0.35", markersize=6), plt.Line2D([], [], marker="D", color="none", markerfacecolor="black", markersize=6)], ["95% block-bootstrap CI", "bootstrap median", "aggregate observed"], frameon=False, fontsize=6.7, loc="lower right")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"stage2c_effective_msrc_validation.{suffix}", dpi=300 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return FIGURES / "stage2c_effective_msrc_validation.png", FIGURES / "stage2c_effective_msrc_validation.pdf"


def write_report(aggregate: dict[str, dict[str, float | str]], loo_rows: list[dict[str, object]], bootstrap_rows_all: list[dict[str, object]], bootstrap_summary_rows: list[dict[str, object]], local_stats: dict[str, object], mouse_tip_check: str) -> Path:
    report = RESULTS / "stage2c_effective_msrc_report.md"
    lines = [
        "# Stage 2C effective MSRC validation",
        "",
        "The symmetric 2:2 MSRC quartet result has one identifiable effective parameter, `Delta`, with `q_arr = (1 + 2 Delta)/3` and each non-arrangement topology `(1 - Delta)/3`. The empirical validation estimates this composite parameter only; it does not estimate biological `t` and `m` separately.",
        "",
        "The three clean configurations rotate which biological topology is arrangement-concordant: STT maps to Q_SPECIES, TST maps to Q_OTHER, and TTS maps to Q_T_ALT. Mus spretus is treated as `outgroup_not_t_haplotype` according to the frozen tip metadata.",
        "",
        "## Aggregate configuration estimates",
        "",
        "| configuration | arrangement topology | q_arr | non-arr 1 | non-arr 2 | unresolved | Delta |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for p in PATTERNS:
        a = aggregate[p]
        lines.append(f"| {p} | {a['arrangement_topology']} | {float(a['q_arrangement']):.6f} | {float(a['q_nonarr_1']):.6f} | {float(a['q_nonarr_2']):.6f} | {float(a['q_unresolved']):.6f} | {float(a['delta_hat']):.6f} |")
    lines += ["", "## Held-out validation", "", "The primary leave-one-configuration-out fit gives the two training configurations equal weight. STT-only predictions are reported separately as a simple demonstration. The two non-arrangement frequencies are expected to be equal under the symmetric model; their observed difference is retained as a model-departure diagnostic.", "", "| analysis | held-out | training | Delta | max abs error | L1 | arrangement-support error | non-arrangement symmetry error |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for r in loo_rows:
        if r["analysis"] in {"leave_one_configuration_out", "STT_only_prediction"}:
            lines.append(f"| {r['analysis']} | {r['configuration']} | {r['training_configurations']} | {float(r['delta_eff']):.6f} | {float(r['max_abs_error']):.6f} | {float(r['l1_distance']):.6f} | {float(r['arrangement_support_error']):+.6f} | {float(r['nonarrangement_symmetry_error']):.6f} |")
    lines += ["", "## Linkage-aware block bootstrap", "", "Physical blocks were used as resampling units to quantify uncertainty in the aggregate held-out predictions while preserving local linkage and spatial heterogeneity. Blocks are nonoverlapping intervals aligned to the first retained 5-kb window at 5 Mb. The primary analysis uses 1-Mb blocks and 10,000 bootstrap replicates; fixed 500-kb and 2-Mb analyses use 5,000 replicates each. Each replicate samples complete blocks with replacement until the original number of blocks is restored, aggregates quartet counts, and gives the two training configurations equal weight.", "", "| block size | held-out | n_bootstrap | Delta train median [95% CI] | max error median [95% CI] | L1 median [95% CI] |", "|---:|---|---:|---|---|---|"]
    for r in bootstrap_summary_rows:
        if r["record_type"] == "heldout_configuration":
            lines.append(f"| {r['block_size_bp']} | {r['configuration']} | {r['n_bootstrap']} | {float(r['delta_train_median']):.6f} [{float(r['delta_train_ci025']):.6f}, {float(r['delta_train_ci975']):.6f}] | {float(r['max_abs_error_median']):.6f} [{float(r['max_abs_error_ci025']):.6f}, {float(r['max_abs_error_ci975']):.6f}] | {float(r['l1_median']):.6f} [{float(r['l1_ci025']):.6f}, {float(r['l1_ci975']):.6f}] |")
    lines += ["", "## Local spatial heterogeneity diagnostic", "", "A constant Delta is not intended to predict every linked physical block. For each 1-Mb block, `Delta_block = (3 q_arr - 1)/2` is retained without truncation; negative values therefore indicate local model departure rather than an invalid estimate. The bootstrap is the aggregate robustness analysis, while these block estimates describe local heterogeneity.", "", "| configuration | Delta range | median | IQR | fraction < 0 | fraction > 1 |", "|---|---|---:|---|---:|---:|"]
    for p in PATTERNS:
        s = local_stats[p]
        lines.append(f"| {p} | [{s['range'][0]:.6f}, {s['range'][1]:.6f}] | {s['median']:.6f} | [{s['iqr'][0]:.6f}, {s['iqr'][1]:.6f}] | {s['fraction_negative']:.3f} | {s['fraction_gt1']:.3f} |")
    lines += ["", f"Spearman correlations of block Delta estimates: STT/TST = {local_stats['spearman_STT_TST']:.4f}; STT/TTS = {local_stats['spearman_STT_TTS']:.4f}.", "", "No conventional multinomial p-values are reported: induced quartets reuse individuals, windows, and linked genomic segments, so they are not independent replicates. The old leave-one-single-block-out table is retained for provenance, but it is not interpreted as the preferred robustness analysis.", "", "The held-out validation is a low-dimensional effective MSRC demonstration. TTS has a modest non-arrangement asymmetry (~0.031), indicating departure from the simplest symmetric model. The result does not establish that the entire mouse genealogy is generated solely by MSRC, and it does not separate the biological switching time and migration parameters.", "", mouse_tip_check]
    report.write_text("\n".join(lines) + "\n")
    return report


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    # Explicitly verify the frozen Mus spretus metadata before using the 2:2 mapping.
    tip_rows = read_tsv(ROOT / "data/house_mouse_t_complex/metadata/tip_mapping.tsv")
    spretus = [r for r in tip_rows if r["subspecies"] == "Mus spretus"]
    assert spretus and all(r["t_status"] == "outgroup_not_t_haplotype" for r in spretus)
    mouse_tip_check = f"Frozen tip metadata check: {len(spretus)} Mus spretus tips are `outgroup_not_t_haplotype`."
    aggregate = load_aggregate()
    validation_rows = aggregate_rows(aggregate) + prediction_rows(aggregate)
    validation_fields = ["analysis", "configuration", "training_configurations", "arrangement_topology", "q_species_observed", "q_t_alt_observed", "q_other_observed", "q_arrangement_observed", "q_nonarr_1_observed", "q_nonarr_2_observed", "q_unresolved_observed", "delta_eff", "q_species_predicted", "q_t_alt_predicted", "q_other_predicted", "max_abs_error", "l1_distance", "arrangement_support_error", "nonarrangement_symmetry_error", "nonarrangement_prediction_error"]
    write_tsv(RESULTS / "stage2c_effective_msrc_validation.tsv", validation_rows, validation_fields)
    block_all: list[dict[str, object]] = []
    block_maps: dict[int, dict[int, dict[str, dict[str, object]]]] = {}
    for block_size in BLOCK_SIZES:
        block_map = load_block_rows(block_size); block_maps[block_size] = block_map; block_all.extend(block_cv_rows(block_map, block_size))
    block_fields = ["record_type", "block_size_bp", "test_block", "test_midpoint_mb", "heldout_configuration", "training_configurations", "n_test_windows", "n_train_blocks", "n_folds", "delta_eff_train", "q_species_observed", "q_t_alt_observed", "q_other_observed", "q_species_predicted", "q_t_alt_predicted", "q_other_predicted", "q_unresolved_observed", "max_abs_error", "l1_distance", "arrangement_support_error", "nonarrangement_symmetry_error", "nonarrangement_prediction_error"]
    write_tsv(RESULTS / "stage2c_effective_msrc_block_cv.tsv", block_all, block_fields)
    bootstrap_all: list[dict[str, object]] = []
    bootstrap_summaries: list[dict[str, object]] = []
    for block_size in BLOCK_SIZES:
        n_bootstrap = BOOTSTRAP_REPLICATES[block_size]
        rows = bootstrap_rows(block_maps[block_size], block_size, n_bootstrap, BOOTSTRAP_SEED)
        bootstrap_all.extend(rows)
        bootstrap_summaries.extend(bootstrap_summary(rows, block_size, n_bootstrap))
    bootstrap_fields = ["block_size_bp", "bootstrap_replicate", "heldout_configuration", "training_configurations", "n_sampled_blocks", "delta_STT", "delta_TST", "delta_TTS", "delta_train", "q_species_observed", "q_t_alt_observed", "q_other_observed", "q_species_predicted", "q_t_alt_predicted", "q_other_predicted", "max_abs_error", "l1_distance", "arrangement_support_error", "nonarrangement_symmetry_error"]
    write_tsv(RESULTS / "stage2c_effective_msrc_block_bootstrap.tsv", bootstrap_all, bootstrap_fields)
    summary_fields = ["block_size_bp", "record_type", "configuration", "n_bootstrap", "delta_train_median", "delta_train_ci025", "delta_train_ci975", "max_abs_error_median", "max_abs_error_ci025", "max_abs_error_ci975", "l1_median", "l1_ci025", "l1_ci975", "arrangement_support_error_median", "arrangement_support_error_ci025", "arrangement_support_error_ci975", "nonarrangement_symmetry_error_median", "nonarrangement_symmetry_error_ci025", "nonarrangement_symmetry_error_ci975"]
    write_tsv(RESULTS / "stage2c_effective_msrc_bootstrap_summary.tsv", bootstrap_summaries, summary_fields)
    local_rows = local_heterogeneity_rows(block_maps[PRIMARY_BLOCK_SIZE], PRIMARY_BLOCK_SIZE)
    local_fields = ["block_size_bp", "block_id", "block_midpoint_mb", "configuration", "n_windows", "q_species", "q_t_alt", "q_other", "q_unresolved", "delta_block"]
    write_tsv(RESULTS / "stage2c_effective_msrc_local_heterogeneity.tsv", local_rows, local_fields)
    local_stats = local_heterogeneity_summary(local_rows)
    loo_rows = prediction_rows(aggregate)
    png, pdf = build_figure(aggregate, loo_rows, [r for r in bootstrap_summaries if r["block_size_bp"] == PRIMARY_BLOCK_SIZE])
    report = write_report(aggregate, loo_rows, bootstrap_all, bootstrap_summaries, local_stats, mouse_tip_check)
    print("Aggregate Delta estimates:")
    for p in PATTERNS: print(p, f(float(aggregate[p]["delta_hat"])))
    print("\nLeave-one-configuration-out and STT-only predictions:")
    for r in loo_rows: print(r["analysis"], r["configuration"], f(float(r["delta_eff"])), f(float(r["max_abs_error"])), f(float(r["l1_distance"])))
    print("\nBlock-CV summary:")
    for r in block_all:
        if r["record_type"] == "summary": print(r["block_size_bp"], r["heldout_configuration"], f(float(r["max_abs_error"])), f(float(r["l1_distance"])))
    print("\nLocal heterogeneity summary:")
    for p in PATTERNS: print(p, local_stats[p])
    print("Spearman STT/TST", f(local_stats["spearman_STT_TST"]), "STT/TTS", f(local_stats["spearman_STT_TTS"]))
    print("\nOutputs:", RESULTS / "stage2c_effective_msrc_validation.tsv", RESULTS / "stage2c_effective_msrc_block_cv.tsv", RESULTS / "stage2c_effective_msrc_block_bootstrap.tsv", RESULTS / "stage2c_effective_msrc_bootstrap_summary.tsv", RESULTS / "stage2c_effective_msrc_local_heterogeneity.tsv", report, png, pdf)


if __name__ == "__main__":
    main()
