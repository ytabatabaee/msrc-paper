#!/usr/bin/env python3
"""Low-dimensional, block-held-out validation of the symmetric 2:2 MSRC quartet result."""

from __future__ import annotations

import csv
import math
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


def build_figure(aggregate: dict[str, dict[str, float | str]], loo_rows: list[dict[str, object]], blocks_1mb: dict[int, dict[str, dict[str, object]]]) -> tuple[Path, Path]:
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
    for p in PATTERNS:
        ys = [float(blocks_1mb[b][p]["delta"]) for b in sorted(blocks_1mb)]
        xs = [blocks_1mb[b][p]["midpoint"] / 1e6 for b in sorted(blocks_1mb)]
        ax.scatter(xs, ys, s=18, color=colors[p], alpha=0.78, label=p)
        ax.axhline(float(aggregate[p]["delta_hat"]), color=colors[p], lw=1, ls="--", alpha=0.7)
    ax.set(xlabel="chr17 block midpoint (Mb)", ylabel="Effective Δ", title="C. Effective Δ across 1-Mb blocks")
    ax.set_ylim(0.35, 0.75); ax.grid(alpha=0.18, lw=0.45); ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper right")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"stage2c_effective_msrc_validation.{suffix}", dpi=300 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return FIGURES / "stage2c_effective_msrc_validation.png", FIGURES / "stage2c_effective_msrc_validation.pdf"


def write_report(aggregate: dict[str, dict[str, float | str]], loo_rows: list[dict[str, object]], block_rows: list[dict[str, object]], mouse_tip_check: str) -> Path:
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
    lines += ["", "## Spatially blocked validation", "", "Blocks are nonoverlapping physical intervals aligned to the first retained 5-kb window at 5 Mb. The primary block size is 1 Mb because it is large relative to the 5-kb window spacing while retaining multiple spatial folds; 500-kb and 2-Mb runs are fixed sensitivity analyses, not choices optimized for accuracy. Complete blocks, never individual windows, are assigned to the held-out fold. For each held-out block and configuration, Delta is fit from the other two configurations in all remaining blocks with equal configuration weight.", "", "| block size | configuration | folds | mean max abs error | mean L1 | mean arrangement-support error | mean non-arr symmetry error |", "|---:|---|---:|---:|---:|---:|---:|"]
    for r in block_rows:
        if r["record_type"] == "summary":
            lines.append(f"| {r['block_size_bp']} | {r['heldout_configuration']} | {r['n_folds']} | {float(r['max_abs_error']):.6f} | {float(r['l1_distance']):.6f} | {float(r['arrangement_support_error']):+.6f} | {float(r['nonarrangement_symmetry_error']):.6f} |")
    lines += ["", "No conventional multinomial p-values are reported: induced quartets reuse individuals, windows, and linked genomic segments, so they are not independent replicates. The block-CV results are the preferred robustness analysis.", "", "The held-out validation is a low-dimensional effective MSRC demonstration. It does not establish that the mouse affected genealogy is generated solely by MSRC, and it does not separate the biological switching time and migration parameters.", "", mouse_tip_check]
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
    loo_rows = prediction_rows(aggregate)
    png, pdf = build_figure(aggregate, loo_rows, block_maps[PRIMARY_BLOCK_SIZE])
    report = write_report(aggregate, loo_rows, block_all, mouse_tip_check)
    print("Aggregate Delta estimates:")
    for p in PATTERNS: print(p, f(float(aggregate[p]["delta_hat"])))
    print("\nLeave-one-configuration-out and STT-only predictions:")
    for r in loo_rows: print(r["analysis"], r["configuration"], f(float(r["delta_eff"])), f(float(r["max_abs_error"])), f(float(r["l1_distance"])))
    print("\nBlock-CV summary:")
    for r in block_all:
        if r["record_type"] == "summary": print(r["block_size_bp"], r["heldout_configuration"], f(float(r["max_abs_error"])), f(float(r["l1_distance"])))
    print("\nOutputs:", RESULTS / "stage2c_effective_msrc_validation.tsv", RESULTS / "stage2c_effective_msrc_block_cv.tsv", report, png, pdf)


if __name__ == "__main__":
    main()
