#!/usr/bin/env python3
"""Test frozen Stage 3B structural predictions against chr4 QQS tracks."""

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
from collections import Counter, defaultdict
from pathlib import Path


RESULTS = Path(__file__).resolve().parents[1] / "results"
FIGURES = Path(__file__).resolve().parents[1] / "figures"
PRED_MANIFEST = RESULTS / "stage3b_structural_prediction_manifest.json"
PRED_PRIMARY = RESULTS / "structural_quartet_predictions.tsv"
PRED_ALL = RESULTS / "structural_quartet_predictions_all.tsv"
LOCI = RESULTS / "locus_table_chr4.tsv"
STAGE27C_MANIFEST = RESULTS / "denovo_breakpoint_manifest.json"
STAGE3A_MANIFEST = RESULTS / "stage3a_final_manifest.json"
Q = ["q1", "q2", "q3"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


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


def mean(vals: list[float]) -> float:
    return statistics.mean(vals) if vals else float("nan")


def margin(row: dict[str, object], pred: str) -> float:
    if pred not in Q:
        return float("nan")
    others = [q for q in Q if q != pred]
    return float(row[f"mean_{pred}"]) - (float(row[f"mean_{others[0]}"]) + float(row[f"mean_{others[1]}"])) / 2


def support_for_predictions(preds: list[dict[str, str]], loci: list[dict[str, str]]) -> list[dict[str, object]]:
    by_clade = defaultdict(list)
    for row in loci:
        by_clade[row["clade"]].append(row)
    out = []
    seen = set()
    for pred in preds:
        key = (pred["clade"], pred["structural_set"], pred["anchor_size"], pred["segment_id"])
        if key in seen:
            continue
        seen.add(key)
        start, end = int(pred["start"]), int(pred["end"])
        rows = [r for r in by_clade[pred["clade"]] if start <= float(r["midpoint"]) <= end]
        means = {q: mean([float(r[q]) for r in rows]) for q in Q}
        finite = [(q, means[q]) for q in Q if not math.isnan(means[q])]
        dominant = max(finite, key=lambda x: x[1])[0] if finite else ""
        out.append(
            {
                "clade": pred["clade"],
                "structural_set": pred["structural_set"],
                "anchor_size": pred["anchor_size"],
                "segment_id": pred["segment_id"],
                "start": start,
                "end": end,
                "n_loci": len(rows),
                "mean_q1": means["q1"],
                "mean_q2": means["q2"],
                "mean_q3": means["q3"],
                "dominant_observed_topology": dominant,
            }
        )
    return out


def empirical_p(obs: float, null: list[float]) -> float:
    if math.isnan(obs) or not null:
        return float("nan")
    return (1 + sum(x >= obs for x in null)) / (len(null) + 1)


def summarize(rows: list[dict[str, object]], pred_col: str, margin_col: str, analysis_type: str, permutations: int, seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary = []
    perm_rows = []
    rng = random.Random(seed)
    for (struct_set, anchor, clade), group in sorted(defaultdict(list, {k: [] for k in []}).items()):
        pass
    by_group = defaultdict(list)
    for row in rows:
        if row.get(pred_col) in Q:
            by_group[(row["structural_set"], row["anchor_size"], row["clade"])].append(row)
    for (struct_set, anchor, clade), group in sorted(by_group.items()):
        obs_m = [float(r[margin_col]) for r in group if not math.isnan(float(r[margin_col]))]
        obs = mean(obs_m)
        labels = [r[pred_col] for r in group]
        null = []
        for i in range(permutations):
            shuf = labels[:]
            rng.shuffle(shuf)
            vals = [margin(r, lab) for r, lab in zip(group, shuf)]
            stat = mean([v for v in vals if not math.isnan(v)])
            null.append(stat)
            if struct_set == "primary" and int(anchor) == 100000:
                perm_rows.append({"analysis_type": analysis_type, "clade": clade, "permutation": i + 1, "statistic": stat})
        summary.append(make_summary_row(clade, analysis_type, group, pred_col, margin_col, obs, empirical_p(obs, null)))
    return summary, perm_rows


def make_summary_row(clade: str, analysis_type: str, group: list[dict[str, object]], pred_col: str, margin_col: str, obs: float, p: float, circ_p: float = float("nan")) -> dict[str, object]:
    margins = [float(r[margin_col]) for r in group if not math.isnan(float(r[margin_col]))]
    pred_q = [float(r[f"mean_{r[pred_col]}"]) for r in group if r[pred_col] in Q]
    alt_q = []
    for r in group:
        if r[pred_col] in Q:
            alt_q.extend(float(r[f"mean_{q}"]) for q in Q if q != r[pred_col])
    return {
        "structural_set": group[0]["structural_set"] if group else "",
        "anchor_size": group[0]["anchor_size"] if group else "",
        "clade": clade,
        "analysis_type": analysis_type,
        "n_segments": len(group),
        "mean_structural_margin": mean([float(r.get("generalized_margin", "nan") or "nan") for r in group]),
        "median_structural_margin": statistics.median([float(r.get("generalized_margin", "nan") or "nan") for r in group]) if group else float("nan"),
        "fraction_positive_margin": sum(m > 0 for m in margins) / len(margins) if margins else float("nan"),
        "topology_prediction_accuracy": sum(r[pred_col] == r["dominant_observed_topology"] for r in group) / len(group) if group else float("nan"),
        "mean_predicted_q": mean(pred_q),
        "mean_alternative_q": mean(alt_q),
        "mean_segment_margin": obs,
        "median_margin": statistics.median(margins) if margins else float("nan"),
        "permutation_pvalue": p,
        "circular_shift_pvalue": circ_p,
    }


def circular_p(group: list[dict[str, object]], pred_col: str, margin_col: str) -> float:
    if len(group) < 3:
        return float("nan")
    group = sorted(group, key=lambda r: int(r["start"]))
    obs = mean([float(r[margin_col]) for r in group])
    shifted = []
    for shift in range(1, len(group)):
        vals = []
        for i, row in enumerate(group):
            qrow = group[(i + shift) % len(group)]
            fake = dict(row)
            for q in Q:
                fake[f"mean_{q}"] = qrow[f"mean_{q}"]
            vals.append(margin(fake, row[pred_col]))
        shifted.append(mean(vals))
    return empirical_p(obs, shifted)


def combined_summary(rows: list[dict[str, object]], pred_col: str, margin_col: str, analysis_type: str, permutations: int, seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    out = []
    perm_out = []
    rng = random.Random(seed + 101)
    by_set_anchor = defaultdict(list)
    for row in rows:
        if row.get(pred_col) in Q:
            by_set_anchor[(row["structural_set"], row["anchor_size"])].append(row)
    for (struct_set, anchor), group in sorted(by_set_anchor.items()):
        obs = mean([float(r[margin_col]) for r in group])
        by_seg = defaultdict(dict)
        rows_by_seg_clade = {}
        for r in group:
            by_seg[r["segment_id"]][r["clade"]] = r[pred_col]
            rows_by_seg_clade[(r["segment_id"], r["clade"])] = r
        seg_ids = sorted(by_seg)
        vectors = [by_seg[s] for s in seg_ids]
        null = []
        for i in range(permutations):
            shuf = vectors[:]
            rng.shuffle(shuf)
            vals = []
            for seg_id, vec in zip(seg_ids, shuf):
                for clade, label in vec.items():
                    r = rows_by_seg_clade.get((seg_id, clade))
                    if r is not None:
                        vals.append(margin(r, label))
            stat = mean([v for v in vals if not math.isnan(v)])
            null.append(stat)
            if struct_set == "primary" and int(anchor) == 100000:
                perm_out.append({"analysis_type": analysis_type, "clade": "combined", "permutation": i + 1, "statistic": stat})
        out.append(make_summary_row("combined", analysis_type, group, pred_col, margin_col, obs, empirical_p(obs, null)))
    return out, perm_out


def continuous(rows: list[dict[str, object]], permutations: int, seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rng = random.Random(seed + 202)
    summary = []
    perm_rows = []
    by_group = defaultdict(list)
    for r in rows:
        if all(r.get(g, "") not in {"", None} for g in ["G1", "G2", "G3"]) and all(r.get(f"mean_{q}", "") not in {"", None} for q in Q):
            by_group[(r["structural_set"], r["anchor_size"], r["clade"])].append(r)
    for (struct_set, anchor, clade), group in sorted(by_group.items()):
        def stat(gvecs):
            vals = []
            for r, gs in zip(group, gvecs):
                gmean = sum(gs) / 3
                vals.append(sum((gs[i] - gmean) * (float(r[f"mean_q{i+1}"]) - 1 / 3) for i in range(3)))
            return mean(vals)
        gvecs = [[float(r["G1"]), float(r["G2"]), float(r["G3"])] for r in group]
        obs = stat(gvecs)
        null = []
        for i in range(permutations):
            shuf = gvecs[:]
            rng.shuffle(shuf)
            val = stat(shuf)
            null.append(val)
            if struct_set == "primary" and int(anchor) == 100000:
                perm_rows.append({"analysis_type": "continuous_generalized", "clade": clade, "permutation": i + 1, "statistic": val})
        summary.append(
            {
                "structural_set": struct_set,
                "anchor_size": anchor,
                "clade": clade,
                "analysis_type": "continuous_generalized",
                "n_segments": len(group),
                "mean_structural_margin": "",
                "median_structural_margin": "",
                "fraction_positive_margin": "",
                "topology_prediction_accuracy": "",
                "mean_predicted_q": "",
                "mean_alternative_q": "",
                "mean_segment_margin": obs,
                "median_margin": "",
                "permutation_pvalue": empirical_p(obs, null),
                "circular_shift_pvalue": "",
            }
        )
    return summary, perm_rows


def join_predictions_support(preds: list[dict[str, str]], supports: list[dict[str, object]]) -> list[dict[str, object]]:
    support_key = {(r["clade"], r["structural_set"], str(r["anchor_size"]), r["segment_id"]): r for r in supports}
    out = []
    for p in preds:
        s = support_key.get((p["clade"], p["structural_set"], str(p["anchor_size"]), p["segment_id"]))
        if not s:
            continue
        row = dict(p)
        row.update(s)
        row["strict_margin"] = margin(row, row.get("predicted_discrete_topology", ""))
        row["generalized_margin_q"] = margin(row, row.get("predicted_generalized_topology", ""))
        row["prediction_matches_observed"] = row.get("predicted_generalized_topology", "") == row.get("dominant_observed_topology", "")
        out.append(row)
    return out


def sensitivity(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    by_set = defaultdict(list)
    for r in rows:
        by_set[(r["structural_set"], r["clade"])].append(r)
    for (struct_set, clade), group in sorted(by_set.items()):
        by_seg = defaultdict(list)
        for r in group:
            by_seg[r["segment_id"]].append(r)
        stable_disc = stable_gen = comparable = 0
        for sg in by_seg.values():
            if {int(r["anchor_size"]) for r in sg} >= {50000, 100000, 200000}:
                comparable += 1
                stable_disc += len({r.get("predicted_discrete_topology", "") for r in sg}) == 1
                stable_gen += len({r.get("predicted_generalized_topology", "") for r in sg}) == 1
        out.append(
            {
                "structural_set": struct_set,
                "clade": clade,
                "comparable_segments": comparable,
                "fraction_discrete_predictions_unchanged": stable_disc / comparable if comparable else float("nan"),
                "fraction_generalized_predictions_unchanged": stable_gen / comparable if comparable else float("nan"),
                "strict_predictions_stable_all_anchor_sizes": stable_disc,
                "generalized_predictions_stable_all_anchor_sizes": stable_gen,
            }
        )
    return out


def make_figures(rows: list[dict[str, object]], perm_rows: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))
    import matplotlib.pyplot as plt

    primary = [r for r in rows if r["structural_set"] == "primary" and int(r["anchor_size"]) == 100000]
    fig, axes = plt.subplots(4, 1, figsize=(11, 7), sharex=True)
    clades = ["Columbea", "N61", "N62"]
    colors = {"q1": "#2b8cbe", "q2": "#f03b20", "q3": "#31a354", "": "#d9d9d9"}
    for r in primary:
        y = clades.index(r["clade"]) if r["clade"] in clades else 0
        axes[0].plot([int(r["start"]), int(r["end"])], [y, y], lw=6, color="#636363")
        axes[1].plot([int(r["start"]), int(r["end"])], [y, y], lw=8, color=colors.get(r.get("predicted_generalized_topology", ""), "#d9d9d9"))
        axes[3].scatter((int(r["start"]) + int(r["end"])) / 2, float(r["generalized_margin_q"]) if r.get("generalized_margin_q") not in {"", None} else float("nan"), c=colors.get(r.get("predicted_generalized_topology", ""), "#d9d9d9"), s=16)
    for clade in clades:
        cg = sorted([r for r in primary if r["clade"] == clade], key=lambda x: int(x["start"]))
        mids = [(int(r["start"]) + int(r["end"])) / 2 for r in cg]
        for q in Q:
            axes[2].plot(mids, [float(r[f"mean_{q}"]) for r in cg], label=f"{clade} {q}", lw=1)
    axes[0].set_title("A. Primary structural segments")
    axes[1].set_title("B. Frozen generalized structural topology prediction")
    axes[2].set_title("C. Observed segment mean quartet support")
    axes[3].set_title("D. Predicted-topology support margin")
    axes[0].set_yticks(range(len(clades)), clades)
    axes[1].set_yticks(range(len(clades)), clades)
    axes[2].set_ylabel("mean q")
    axes[3].set_ylabel("M")
    axes[3].set_xlabel("GalGal6 chr4 coordinate")
    axes[2].legend(ncol=3, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "chr4_structural_state_predicts_topology.pdf")
    fig.savefig(FIGURES / "chr4_structural_state_predicts_topology.png", dpi=200)
    plt.close(fig)

    gen_perm = [float(r["statistic"]) for r in perm_rows if r["analysis_type"] == "generalized" and r["clade"] in {"N61", "N62"}]
    obs = [float(r["generalized_margin_q"]) for r in primary if r.get("predicted_generalized_topology") in Q]
    fig, ax = plt.subplots(figsize=(6, 4))
    if gen_perm:
        ax.hist(gen_perm, bins=40, color="#bdbdbd", edgecolor="white")
    if obs:
        ax.axvline(mean(obs), color="#cb181d", lw=2, label="observed")
    ax.set_title("Topology-label permutation null")
    ax.set_xlabel("mean segment margin")
    ax.set_ylabel("permutations")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "structural_topology_prediction_null.pdf")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=8675309)
    args = parser.parse_args()

    before = {p: sha256(p) for p in [PRED_MANIFEST, PRED_PRIMARY, PRED_ALL]}
    manifest = json.loads(PRED_MANIFEST.read_text())
    preds = read_tsv(PRED_ALL)
    loci = read_tsv(LOCI)
    support = support_for_predictions(preds, loci)
    primary_support = [r for r in support if r["structural_set"] == "primary" and int(r["anchor_size"]) == 100000]
    write_tsv(RESULTS / "segment_quartet_support.tsv", primary_support, ["clade", "structural_set", "anchor_size", "segment_id", "start", "end", "n_loci", "mean_q1", "mean_q2", "mean_q3", "dominant_observed_topology"])

    rows = join_predictions_support(preds, support)
    fields = ["clade", "structural_set", "anchor_size", "segment_id", "start", "end", "n_loci", "predicted_discrete_topology", "discrete_prediction_status", "predicted_generalized_topology", "prediction_stable_50_100_200kb", "G1", "G2", "G3", "generalized_margin", "mean_q1", "mean_q2", "mean_q3", "dominant_observed_topology", "strict_margin", "generalized_margin_q", "prediction_matches_observed", "structural_confidence", "role_status"]
    write_tsv(RESULTS / "structural_topology_prediction_segment_results.tsv", rows, fields)

    strict_rows = [r for r in rows if r.get("discrete_prediction_status") == "unique_2:2" and r.get("predicted_discrete_topology") in Q]
    gen_rows = [r for r in rows if r.get("predicted_generalized_topology") in Q]
    stable_gen_rows = [r for r in gen_rows if r.get("prediction_stable_50_100_200kb") == "true" or r.get("prediction_stable_50_100_200kb") is True]
    summary = []
    perm_rows = []
    for label, subset, pred_col, margin_col in [
        ("strict_discrete", strict_rows, "predicted_discrete_topology", "strict_margin"),
        ("generalized", gen_rows, "predicted_generalized_topology", "generalized_margin_q"),
        ("generalized_stable_anchors", stable_gen_rows, "predicted_generalized_topology", "generalized_margin_q"),
    ]:
        s, p = summarize(subset, pred_col, margin_col, label, args.permutations, args.seed)
        for row in s:
            group = [r for r in subset if r["clade"] == row["clade"] and row["analysis_type"] == label and r.get(pred_col) in Q and r["structural_set"] == "primary" and int(r["anchor_size"]) == 100000]
            row["circular_shift_pvalue"] = circular_p(group, pred_col, margin_col)
        summary.extend(s)
        perm_rows.extend(p)
        cs, cp = combined_summary(subset, pred_col, margin_col, label, args.permutations, args.seed)
        summary.extend(cs)
        perm_rows.extend(cp)
    cs, cp = continuous(rows, args.permutations, args.seed)
    summary.extend(cs)
    perm_rows.extend(cp)
    write_tsv(RESULTS / "structural_topology_prediction_summary.tsv", summary, ["structural_set", "anchor_size", "clade", "analysis_type", "n_segments", "mean_structural_margin", "median_structural_margin", "fraction_positive_margin", "topology_prediction_accuracy", "mean_predicted_q", "mean_alternative_q", "mean_segment_margin", "median_margin", "permutation_pvalue", "circular_shift_pvalue"])
    write_tsv(RESULTS / "structural_topology_prediction_permutation.tsv", perm_rows, ["analysis_type", "clade", "permutation", "statistic"])
    write_tsv(RESULTS / "structural_topology_prediction_sensitivity.tsv", sensitivity(rows), ["structural_set", "clade", "comparable_segments", "fraction_discrete_predictions_unchanged", "fraction_generalized_predictions_unchanged", "strict_predictions_stable_all_anchor_sizes", "generalized_predictions_stable_all_anchor_sizes"])
    make_figures(rows, perm_rows)

    after = {p: sha256(p) for p in [PRED_MANIFEST, PRED_PRIMARY, PRED_ALL]}
    if before != after:
        raise RuntimeError("Phase B modified frozen structural prediction files")

    primary_rows = [r for r in rows if r["structural_set"] == "primary" and int(r["anchor_size"]) == 100000]
    final = {
        "phase_a_script_checksum": manifest["phase_a_script_sha256"],
        "phase_b_script_checksum": sha256(Path(__file__)),
        "structural_prediction_manifest_checksum": sha256(PRED_MANIFEST),
        "stage27c_manifest_checksum": sha256(STAGE27C_MANIFEST),
        "stage3a_manifest_checksum": sha256(STAGE3A_MANIFEST) if STAGE3A_MANIFEST.exists() else None,
        "exact_role_mapping": manifest["role_mapping"],
        "structural_set_definitions": ["primary", "stringent", "inclusive"],
        "anchor_size_definitions": [50000, 100000, 200000],
        "strict_prediction_rule": "C1=C2!=S=O -> q1; C1=S!=C2=O -> q2; C2=S!=C1=O -> q3",
        "generalized_score_formula": manifest["generalized_score_formula"],
        "number_of_eligible_segments": {
            "primary_100kb_strict": len([r for r in primary_rows if r.get("predicted_discrete_topology") in Q]),
            "primary_100kb_generalized": len([r for r in primary_rows if r.get("predicted_generalized_topology") in Q]),
        },
        "permutation_count": args.permutations,
        "random_seed": args.seed,
        "final_effect_sizes_and_pvalues": summary,
        "sensitivity_summary_sha256": sha256(RESULTS / "structural_topology_prediction_sensitivity.tsv"),
        "anti_circularity": {
            "structural_predictions_checksummed_before_q_read": True,
            "phase_b_modified_prediction_files": False,
            "topology_label_permutations_preserve_boundaries_and_counts": True,
            "combined_track_permutations_preserve_cross_clade_prediction_vector": True,
        },
    }
    with (RESULTS / "stage3b_final_manifest.json").open("w") as handle:
        json.dump(final, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
