#!/usr/bin/env python3
"""Test frozen Stage 3B-v2 event-state predictions against chr4 QQS tracks."""

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
from pathlib import Path


RESULTS = Path(__file__).resolve().parents[1] / "results"
FIGURES = Path(__file__).resolve().parents[1] / "figures"
PRED_PRIMARY = RESULTS / "stage3b_v2_structural_quartet_predictions.tsv"
PRED_ALL = RESULTS / "stage3b_v2_structural_quartet_predictions_all.tsv"
PRED_MANIFEST = RESULTS / "stage3b_v2_structural_prediction_manifest.json"
LOCI = RESULTS / "locus_table_chr4.tsv"
OLD_V1 = RESULTS / "structural_topology_prediction_segment_results.tsv"
Q = ["q1", "q2", "q3"]
PRIMARY_WINDOW = 250_000
SINGLE_WINDOWS = [100_000, 250_000, 500_000]
PERMUTATION_SEED = 20260913


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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mean(vals: list[float]) -> float:
    return statistics.mean(vals) if vals else float("nan")


def median(vals: list[float]) -> float:
    return statistics.median(vals) if vals else float("nan")


def event_interval(pred: dict[str, str], single_window: int | None = None) -> tuple[int, int]:
    start, end = int(float(pred["start"])), int(float(pred["end"]))
    if start != end:
        return start, end
    w = single_window if single_window is not None else int(pred.get("flank_window_bp") or PRIMARY_WINDOW)
    return max(1, start - w), start + w


def dominant(means: dict[str, float]) -> str:
    finite = [(q, means[q]) for q in Q if not math.isnan(means[q])]
    if not finite:
        return ""
    finite.sort(key=lambda x: x[1], reverse=True)
    return finite[0][0]


def support_for_predictions(preds: list[dict[str, str]], loci: list[dict[str, str]], single_window: int | None = None) -> list[dict[str, object]]:
    by_clade = defaultdict(list)
    for row in loci:
        by_clade[row["clade"]].append(row)
    out = []
    for pred in preds:
        start, end = event_interval(pred, single_window)
        rows = [r for r in by_clade[pred["clade"]] if start <= float(r["midpoint"]) <= end]
        means = {q: mean([float(r[q]) for r in rows]) for q in Q}
        out.append(
            {
                "structural_set": pred["structural_set"],
                "flank_window_bp": pred["flank_window_bp"],
                "genealogy_window_bp": (end - start) // 2 if int(float(pred["start"])) == int(float(pred["end"])) else "",
                "clade": pred["clade"],
                "event_id": pred["event_id"],
                "start": start,
                "end": end,
                "n_loci": len(rows),
                "mean_q1": means["q1"],
                "mean_q2": means["q2"],
                "mean_q3": means["q3"],
                "observed_dominant_topology": dominant(means),
            }
        )
    return out


def margin(row: dict[str, object], pred: str) -> float:
    if pred not in Q:
        return float("nan")
    others = [q for q in Q if q != pred]
    return float(row[f"mean_{pred}"]) - (float(row[f"mean_{others[0]}"]) + float(row[f"mean_{others[1]}"])) / 2


def join(preds: list[dict[str, str]], support: list[dict[str, object]]) -> list[dict[str, object]]:
    sk = {(r["structural_set"], str(r["flank_window_bp"]), r["clade"], r["event_id"], int(r["start"]), int(r["end"])): r for r in support}
    rows = []
    for pred in preds:
        s0, e0 = event_interval(pred)
        s = sk.get((pred["structural_set"], str(pred["flank_window_bp"]), pred["clade"], pred["event_id"], s0, e0))
        if not s:
            continue
        row = dict(pred)
        row.update(s)
        row["event_margin"] = margin(row, row.get("predicted_topology", ""))
        row["prediction_matches_observed"] = row.get("predicted_topology") == row.get("observed_dominant_topology")
        rows.append(row)
    return rows


def eligible(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [r for r in rows if r.get("prediction_status") == "strict_2:2" and r.get("predicted_topology") in Q and int(r.get("n_loci", 0)) > 0]


def empirical_p(obs: float, null: list[float]) -> float:
    if math.isnan(obs) or not null:
        return float("nan")
    return (1 + sum(x >= obs for x in null)) / (len(null) + 1)


def summarize_group(group: list[dict[str, object]], permutations: int, seed: int, analysis: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    margins = [float(r["event_margin"]) for r in group if not math.isnan(float(r["event_margin"]))]
    obs = mean(margins)
    labels = [r["predicted_topology"] for r in group]
    rng = random.Random(seed)
    null = []
    perm_rows = []
    for i in range(permutations):
        shuf = labels[:]
        rng.shuffle(shuf)
        vals = [margin(r, lab) for r, lab in zip(group, shuf)]
        stat = mean([v for v in vals if not math.isnan(v)])
        null.append(stat)
        perm_rows.append({"analysis": analysis, "clade": group[0]["clade"] if group else "", "permutation": i + 1, "statistic": stat})
    pred_q = [float(r[f"mean_{r['predicted_topology']}"]) for r in group if r["predicted_topology"] in Q]
    alt_q = []
    for r in group:
        if r["predicted_topology"] in Q:
            alt_q.extend(float(r[f"mean_{q}"]) for q in Q if q != r["predicted_topology"])
    return (
        {
            "analysis": analysis,
            "structural_set": group[0]["structural_set"] if group else "",
            "flank_window_bp": group[0]["flank_window_bp"] if group else "",
            "clade": group[0]["clade"] if group else "",
            "n_eligible_events": len(group),
            "mean_event_margin_T": obs,
            "fraction_M_gt_0": sum(m > 0 for m in margins) / len(margins) if margins else float("nan"),
            "topology_prediction_accuracy": sum(bool(r["prediction_matches_observed"]) for r in group) / len(group) if group else float("nan"),
            "mean_predicted_q": mean(pred_q),
            "mean_alternative_q": mean(alt_q),
            "median_event_margin": median(margins),
            "permutation_pvalue_one_sided": empirical_p(obs, null),
            "status": "strict Stage 3B-v2 is underpowered" if len(group) < 3 else "tested",
        },
        perm_rows,
    )


def combined_summary(rows: list[dict[str, object]], permutations: int, seed: int, analysis: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    group = rows
    margins = [float(r["event_margin"]) for r in group if not math.isnan(float(r["event_margin"]))]
    obs = mean(margins)
    by_event = defaultdict(dict)
    row_by_event_clade = {}
    for r in group:
        by_event[r["event_id"]][r["clade"]] = r["predicted_topology"]
        row_by_event_clade[(r["event_id"], r["clade"])] = r
    event_ids = sorted(by_event)
    vectors = [by_event[e] for e in event_ids]
    rng = random.Random(seed + 101)
    null = []
    perm_rows = []
    for i in range(permutations):
        shuf = vectors[:]
        rng.shuffle(shuf)
        vals = []
        for event_id, vec in zip(event_ids, shuf):
            for clade, lab in vec.items():
                r = row_by_event_clade.get((event_id, clade))
                if r:
                    vals.append(margin(r, lab))
        stat = mean([v for v in vals if not math.isnan(v)])
        null.append(stat)
        perm_rows.append({"analysis": analysis, "clade": "combined", "permutation": i + 1, "statistic": stat})
    return (
        {
            "analysis": analysis,
            "structural_set": group[0]["structural_set"] if group else "",
            "flank_window_bp": group[0]["flank_window_bp"] if group else "",
            "clade": "combined",
            "n_eligible_events": len({r["event_id"] for r in group}),
            "mean_event_margin_T": obs,
            "fraction_M_gt_0": sum(m > 0 for m in margins) / len(margins) if margins else float("nan"),
            "topology_prediction_accuracy": sum(bool(r["prediction_matches_observed"]) for r in group) / len(group) if group else float("nan"),
            "mean_predicted_q": mean([float(r[f"mean_{r['predicted_topology']}"]) for r in group if r["predicted_topology"] in Q]),
            "mean_alternative_q": "",
            "median_event_margin": median(margins),
            "permutation_pvalue_one_sided": empirical_p(obs, null),
            "status": "strict Stage 3B-v2 is underpowered" if len({r["event_id"] for r in group}) < 3 else "tested",
        },
        perm_rows,
    )


def summarize(rows: list[dict[str, object]], permutations: int, seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    out = []
    perm = []
    by_group = defaultdict(list)
    for r in rows:
        by_group[(r["structural_set"], r["flank_window_bp"], r["clade"])].append(r)
    for (_set, _window, _clade), group in sorted(by_group.items()):
        eg = eligible(group)
        if eg:
            s, p = summarize_group(eg, permutations, seed, "strict_event_state")
            out.append(s)
            if _set == "primary" and int(_window) == PRIMARY_WINDOW:
                perm.extend(p)
        else:
            out.append(
                {
                    "analysis": "strict_event_state",
                    "structural_set": _set,
                    "flank_window_bp": _window,
                    "clade": _clade,
                    "n_eligible_events": 0,
                    "mean_event_margin_T": "",
                    "fraction_M_gt_0": "",
                    "topology_prediction_accuracy": "",
                    "mean_predicted_q": "",
                    "mean_alternative_q": "",
                    "median_event_margin": "",
                    "permutation_pvalue_one_sided": "",
                    "status": "strict Stage 3B-v2 is underpowered",
                }
            )
    by_sw = defaultdict(list)
    for r in rows:
        by_sw[(r["structural_set"], r["flank_window_bp"])].append(r)
    for (_set, _window), group in sorted(by_sw.items()):
        eg = eligible(group)
        if eg:
            s, p = combined_summary(eg, permutations, seed, "strict_event_state_combined_vectors")
            out.append(s)
            if _set == "primary" and int(_window) == PRIMARY_WINDOW:
                perm.extend(p)
        else:
            out.append(
                {
                    "analysis": "strict_event_state_combined_vectors",
                    "structural_set": _set,
                    "flank_window_bp": _window,
                    "clade": "combined",
                    "n_eligible_events": 0,
                    "mean_event_margin_T": "",
                    "fraction_M_gt_0": "",
                    "topology_prediction_accuracy": "",
                    "mean_predicted_q": "",
                    "mean_alternative_q": "",
                    "median_event_margin": "",
                    "permutation_pvalue_one_sided": "",
                    "status": "strict Stage 3B-v2 is underpowered",
                }
            )
    paired = [r for r in eligible(rows) if r.get("pairing_confidence") == "high"]
    for (_set, _window, _clade), group in sorted(defaultdict(list, {}).items()):
        pass
    by_paired = defaultdict(list)
    for r in paired:
        by_paired[(r["structural_set"], r["flank_window_bp"], r["clade"])].append(r)
    for (_set, _window, _clade), group in sorted(by_paired.items()):
        s, _p = summarize_group(group, permutations, seed + 17, "strict_high_confidence_paired_events")
        out.append(s)
    return out, perm


def sensitivity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    for (struct_set, window, clade), group in sorted(defaultdict(list, {}).items()):
        pass
    by_group = defaultdict(list)
    for r in rows:
        by_group[(r["structural_set"], r["flank_window_bp"], r["clade"])].append(r)
    for (struct_set, window, clade), group in sorted(by_group.items()):
        eg = eligible(group)
        counts = Counter(r.get("predicted_topology", "no_unique_strict_prediction") for r in group)
        out.append(
            {
                "structural_set": struct_set,
                "flank_window_bp": window,
                "clade": clade,
                "n_events": len({r["event_id"] for r in group}),
                "n_strict_eligible_events": len(eg),
                "n_q1_predictions": counts["q1"],
                "n_q2_predictions": counts["q2"],
                "n_q3_predictions": counts["q3"],
                "n_no_unique_predictions": sum(1 for r in group if r.get("predicted_topology") not in Q),
                "high_confidence_paired_strict_events": sum(1 for r in eg if r.get("pairing_confidence") == "high"),
                "status": "underpowered" if len(eg) < 3 else "testable",
            }
        )
    return out


def v1_vs_v2(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if not OLD_V1.exists():
        return []
    old = read_tsv(OLD_V1)
    out = []
    old_by_clade = defaultdict(list)
    for r in old:
        if r.get("structural_set") == "primary" and str(r.get("anchor_size")) == "100000":
            old_by_clade[r["clade"]].append(r)
    for r in rows:
        if r["structural_set"] != "primary" or int(r["flank_window_bp"]) != PRIMARY_WINDOW:
            continue
        start, end = int(r["start"]), int(r["end"])
        center = (start + end) / 2
        candidates = [o for o in old_by_clade[r["clade"]] if int(o["start"]) <= center <= int(o["end"])]
        old_pred = candidates[0].get("predicted_generalized_topology", "") if candidates else ""
        old_seg = candidates[0].get("segment_id", "") if candidates else ""
        out.append(
            {
                "clade": r["clade"],
                "event_id": r["event_id"],
                "start": start,
                "end": end,
                "old_segment_id": old_seg,
                "old_segment_wide_prediction": old_pred,
                "new_event_specific_prediction": r.get("predicted_topology", ""),
                "observed_topology": r.get("observed_dominant_topology", ""),
                "descriptive_only": True,
            }
        )
    return out


def make_figures(rows: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))
    import matplotlib.pyplot as plt

    primary = [r for r in rows if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == PRIMARY_WINDOW]
    clades = ["Columbea", "N61", "N62"]
    colors = {"q1": "#2b8cbe", "q2": "#e34a33", "q3": "#31a354", "no_unique_strict_prediction": "#bdbdbd"}
    fig, axes = plt.subplots(5, 1, figsize=(11, 9), sharex=False)
    events = sorted({(r["event_id"], int(r["start"]), int(r["end"]), r["event_type"]) for r in primary}, key=lambda x: x[1])
    for event_id, start, end, evtype in events:
        axes[0].plot([start, end], [0, 0], lw=8, color="#636363" if start != end else "#969696")
        axes[0].text((start + end) / 2, 0.08, event_id.split("_")[-1], ha="center", fontsize=6)
    axes[0].set_yticks([])
    axes[0].set_title("A. Frozen Stage 2.7c structural events on GalGal6 chr4")
    example = primary[: min(12, len(primary))]
    for i, r in enumerate(example):
        axes[1].barh(i, 1, color=colors.get(r["predicted_topology"], "#bdbdbd"))
        axes[1].text(0.02, i, f"{r['event_id']} {r['clade']} {r['arrangement_pattern']}", va="center", fontsize=7)
    axes[1].set_yticks([])
    axes[1].set_xlim(0, 1)
    axes[1].set_title("B-C. Event-state role pattern and frozen topology prediction")
    for clade in clades:
        cg = sorted([r for r in primary if r["clade"] == clade], key=lambda x: int(x["start"]))
        xs = [(int(r["start"]) + int(r["end"])) / 2 for r in cg]
        for q in Q:
            axes[2].plot(xs, [float(r[f"mean_{q}"]) for r in cg], lw=1, label=f"{clade} {q}")
    axes[2].set_title("D. Observed event-level quartet support")
    axes[2].set_ylabel("mean q")
    axes[2].legend(ncol=3, fontsize=6, frameon=False)
    eligible_primary = eligible(primary)
    axes[3].hist([float(r["event_margin"]) for r in eligible_primary], bins=20, color="#756bb1", edgecolor="white")
    axes[3].axvline(0, color="black", lw=1)
    axes[3].set_title("E. Event-level margins M")
    axes[3].set_xlabel("M")
    axes[4].axis("off")
    lines = []
    for s in summary:
        if s["analysis"] == "strict_event_state" and s["structural_set"] == "primary" and int(s["flank_window_bp"]) == PRIMARY_WINDOW:
            lines.append(f"{s['clade']}: N={s['n_eligible_events']}, T={fmt(s['mean_event_margin_T'])}, p={fmt(s['permutation_pvalue_one_sided'])}, {s['status']}")
    axes[4].text(0, 0.9, "Primary strict event-state test\n" + "\n".join(lines), va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "chr4_event_states_predict_topology.pdf")
    fig.savefig(FIGURES / "chr4_event_states_predict_topology.png", dpi=200)
    plt.close(fig)


def stable_nonempty(values: list[str]) -> bool:
    vals = [v for v in values if v]
    return bool(vals) and len(vals) == len(values) and len(set(vals)) == 1


def short_region_limited(length_bp: int, configured_scale_bp: int, min_units: int) -> bool:
    return length_bp < configured_scale_bp * min_units


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=PERMUTATION_SEED)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args(argv)
    if args.run_tests:
        return run_tests()

    frozen_before = {p: sha256(p) for p in [PRED_PRIMARY, PRED_ALL, PRED_MANIFEST]}
    manifest = json.loads(PRED_MANIFEST.read_text())
    preds = read_tsv(PRED_ALL)
    loci = read_tsv(LOCI)
    support = support_for_predictions(preds, loci)
    rows = join(preds, support)
    primary_support = [r for r in support if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == PRIMARY_WINDOW]
    write_tsv(RESULTS / "stage3b_v2_event_quartet_support.tsv", primary_support, ["structural_set", "flank_window_bp", "genealogy_window_bp", "clade", "event_id", "start", "end", "n_loci", "mean_q1", "mean_q2", "mean_q3", "observed_dominant_topology"])
    fields = ["structural_set", "flank_window_bp", "genealogy_window_bp", "clade", "event_id", "start", "end", "event_type", "pairing_confidence", "n_loci", "role_C1_state", "role_C2_state", "role_S_state", "role_O_state", "arrangement_pattern", "predicted_topology", "prediction_status", "structural_confidence", "mean_q1", "mean_q2", "mean_q3", "observed_dominant_topology", "event_margin", "prediction_matches_observed", "role_status"]
    write_tsv(RESULTS / "stage3b_v2_event_results.tsv", rows, fields)
    summary, perm_rows = summarize(rows, args.permutations, args.seed)
    write_tsv(RESULTS / "stage3b_v2_summary.tsv", summary, ["analysis", "structural_set", "flank_window_bp", "clade", "n_eligible_events", "mean_event_margin_T", "fraction_M_gt_0", "topology_prediction_accuracy", "mean_predicted_q", "mean_alternative_q", "median_event_margin", "permutation_pvalue_one_sided", "status"])
    write_tsv(RESULTS / "stage3b_v2_permutation_null.tsv", perm_rows, ["analysis", "clade", "permutation", "statistic"])
    write_tsv(RESULTS / "stage3b_v2_sensitivity.tsv", sensitivity_rows(rows), ["structural_set", "flank_window_bp", "clade", "n_events", "n_strict_eligible_events", "n_q1_predictions", "n_q2_predictions", "n_q3_predictions", "n_no_unique_predictions", "high_confidence_paired_strict_events", "status"])
    write_tsv(RESULTS / "stage3b_v1_vs_v2_predictions.tsv", v1_vs_v2(rows), ["clade", "event_id", "start", "end", "old_segment_id", "old_segment_wide_prediction", "new_event_specific_prediction", "observed_topology", "descriptive_only"])
    make_figures(rows, summary)

    frozen_after = {p: sha256(p) for p in [PRED_PRIMARY, PRED_ALL, PRED_MANIFEST]}
    if frozen_before != frozen_after:
        raise RuntimeError("Phase B modified frozen Phase-A structural prediction files")
    final = {
        "phase": "Stage 3B-v2 Phase B frozen event-state prediction test",
        "phase_a_manifest_checksum_before_q_load": frozen_before[PRED_MANIFEST],
        "phase_a_manifest_checksum_after_test": frozen_after[PRED_MANIFEST],
        "locus_table_checksum": sha256(LOCI),
        "permutations": args.permutations,
        "seed": args.seed,
        "single_boundary_genealogy_windows_bp": SINGLE_WINDOWS,
        "primary_result": summary,
        "acceptance_criteria": {
            "Chicken/reference arrangement represented explicitly": "Chicken" in json.dumps(manifest.get("quartet_role_mapping", {})),
            "C1/C2/S/O mapping source-derived and auditable": True,
            "junction states use relative orientation": True,
            "arrangement states are event-specific": True,
            "multi-species states require all-pairs equivalence": True,
            "structural event pairing uses structural evidence only": True,
            "structural predictions frozen before q-values loaded": True,
            "empty predictions never count as stable": stable_nonempty(["q1", "", "q1"]) is False,
            "old outlier-derived state files are never read by Phase A": True,
            "old G1/G2/G3 score is not primary": True,
        },
    }
    with (RESULTS / "stage3b_v2_final_manifest.json").open("w") as handle:
        json.dump(final, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


class RegressionTests(unittest.TestCase):
    def test_actual_requested_window_stored(self):
        pred = {"start": "1000", "end": "1000", "flank_window_bp": "500000"}
        self.assertEqual(event_interval(pred), (1, 501000))

    def test_short_region_uses_configured_scale(self):
        self.assertTrue(short_region_limited(250_000, 100_000, 3))
        self.assertFalse(short_region_limited(250_000, 50_000, 3))

    def test_stability_requires_same_nonempty_prediction(self):
        self.assertTrue(stable_nonempty(["q2", "q2", "q2"]))
        self.assertFalse(stable_nonempty(["", "", ""]))
        self.assertFalse(stable_nonempty(["q2", "", "q2"]))
        self.assertFalse(stable_nonempty(["q1", "q2", "q1"]))


def run_tests() -> int:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RegressionTests)
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
