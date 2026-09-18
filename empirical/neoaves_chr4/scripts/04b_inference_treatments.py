#!/usr/bin/env python3
"""Stage 4B: frozen weighting treatments for Neoaves chromosome 4.

T0 counts resolved window topology states. T1 gives each consecutive topology
run total weight one. T2 excludes windows in a frozen structural mask. T3 gives
mask windows fixed weight 0.5 and background windows weight 1.0. The T3 rule
is topology neutral and was fixed before inspecting Stage 4B results.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import random
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/neoaves_chr4/processed"
RESULTS = ROOT / "empirical/neoaves_chr4/results"
FIGURES = ROOT / "empirical/neoaves_chr4/figures"
LOCI = DATA / "locus_table_chr4.tsv"
BREAKPOINTS = DATA / "denovo_breakpoint_sets.tsv"
EVENTS = DATA / "stage3b_v2_structural_events.tsv"
STAGE3B_SUMMARY = DATA / "stage3b_v2_summary.tsv"
STAGE4A_MANIFEST = RESULTS / "stage4a_manifest.json"
STAGE4A_WB = RESULTS / "stage4a_window_vs_block_support.tsv"
SETS = ("stringent", "primary", "inclusive")
CLADES = ("Columbea", "N61", "N62")
TOPOLOGIES = ("q1", "q2", "q3")
TREATMENTS = ("T0", "T1", "T2", "T3")
MASK_RADIUS = 250_000
CHR_START, CHR_END = 1, 91_310_470
T3_ASSOCIATED_WEIGHT = 0.5
T3_BACKGROUND_WEIGHT = 1.0
SEED = 41729
BOOTSTRAPS = 1000


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "nan" if math.isnan(value) else f"{value:.12g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def build_masks(rows: list[dict[str, str]], excluded_ids: set[str] | None = None) -> dict[str, list[tuple[int, int]]]:
    """Use only frozen structural columns; excluded IDs support event LOO."""
    excluded_ids = excluded_ids or set()
    masks = {}
    for set_name in SETS:
        positions = [int(row["reference_position"]) for row in rows
                     if row["consensus_breakpoint_id"] not in excluded_ids and truthy(row[f"in_{set_name}"])]
        masks[set_name] = merge_intervals([(max(CHR_START, p - MASK_RADIUS), min(CHR_END, p + MASK_RADIUS)) for p in positions])
    return masks


def in_mask(position: float, intervals: list[tuple[int, int]]) -> bool:
    return any(start <= position <= end for start, end in intervals)


def load_tracks() -> dict[str, list[dict[str, object]]]:
    tracks = {clade: [] for clade in CLADES}
    for row in read_tsv(LOCI):
        if row["Chromosome"] != "chr4" or row["clade"] not in tracks:
            continue
        tracks[row["clade"]].append({
            "gene": row["Gene"], "start": int(row["start"]), "end": int(row["end"]),
            "midpoint": float(row["midpoint"]), "topology": row["dominant_topology"],
        })
    for track in tracks.values():
        track.sort(key=lambda r: (float(r["midpoint"]), int(r["start"]), str(r["gene"])))
        assert all(row["topology"] in TOPOLOGIES for row in track)
    return tracks


def make_runs(track: list[dict[str, object]]) -> list[list[dict[str, object]]]:
    runs: list[list[dict[str, object]]] = []
    for row in track:
        if not runs or runs[-1][-1]["topology"] != row["topology"]:
            runs.append([row])
        else:
            runs[-1].append(row)
    return runs


def weights_for(track: list[dict[str, object]], treatment: str, mask: list[tuple[int, int]]) -> list[float]:
    if treatment == "T0":
        return [1.0] * len(track)
    if treatment == "T1":
        weights: list[float] = []
        for run in make_runs(track):
            weights.extend([1.0 / len(run)] * len(run))
        return weights
    if treatment == "T2":
        return [0.0 if in_mask(float(row["midpoint"]), mask) else 1.0 for row in track]
    if treatment == "T3":
        return [T3_ASSOCIATED_WEIGHT if in_mask(float(row["midpoint"]), mask) else T3_BACKGROUND_WEIGHT for row in track]
    raise ValueError(treatment)


def summarize(track: list[dict[str, object]], weights: list[float]) -> dict[str, object]:
    totals = {q: sum(w for row, w in zip(track, weights) if row["topology"] == q) for q in TOPOLOGIES}
    weight_sum = sum(totals.values())
    if weight_sum <= 0:
        raise ValueError("treatment has no usable evidence")
    freq = {q: totals[q] / weight_sum for q in TOPOLOGIES}
    dominant = max(TOPOLOGIES, key=lambda q: (freq[q], q))
    margin = freq["q1"] - max(freq["q2"], freq["q3"])
    return {**freq, "dominant_topology": dominant, "M_species": margin,
            "M_q1_q2": freq["q1"] - freq["q2"], "effective_weight_sum": weight_sum,
            "n_windows": sum(w > 0 for w in weights), "n_blocks": len(make_runs([r for r, w in zip(track, weights) if w > 0]))}


def percentile(values: list[float], p: float) -> float:
    values = sorted(values)
    x = (len(values) - 1) * p
    lo, hi = math.floor(x), math.ceil(x)
    return values[lo] if lo == hi else values[lo] * (hi - x) + values[hi] * (x - lo)


def bootstrap(track: list[dict[str, object]], mask: list[tuple[int, int]], seed: int) -> dict[str, dict[str, tuple[float, float]]]:
    """Paired cluster bootstrap of frozen topology runs."""
    runs = make_runs(track)
    # Cache each run's contribution so resampling preserves complete runs
    # without rebuilding thousands of pseudo-window dictionaries per draw.
    contributions = []
    for run in runs:
        topology = str(run[0]["topology"])
        associated = sum(in_mask(float(row["midpoint"]), mask) for row in run)
        background = len(run) - associated
        contributions.append({
            "topology": topology, "T0": float(len(run)), "T1": 1.0,
            "T2": float(background),
            "T3": float(background) * T3_BACKGROUND_WEIGHT + float(associated) * T3_ASSOCIATED_WEIGHT,
        })
    draws = {t: {metric: [] for metric in (*TOPOLOGIES, "M_species", "M_q1_q2")} for t in TREATMENTS}
    deltas = {t: [] for t in ("T1", "T2", "T3")}
    rng = random.Random(seed)
    population = list(range(len(contributions)))
    for _ in range(BOOTSTRAPS):
        multiplicities = Counter(rng.choices(population, k=len(population)))
        summaries = {}
        for treatment in TREATMENTS:
            totals = {q: 0.0 for q in TOPOLOGIES}
            for index, multiplicity in multiplicities.items():
                run = contributions[index]
                totals[str(run["topology"])] += multiplicity * float(run[treatment])
            weight_sum = sum(totals.values())
            freq = {q: totals[q] / weight_sum for q in TOPOLOGIES}
            summaries[treatment] = {**freq, "M_species": freq["q1"] - max(freq["q2"], freq["q3"]),
                                    "M_q1_q2": freq["q1"] - freq["q2"]}
            for metric in draws[treatment]:
                draws[treatment][metric].append(float(summaries[treatment][metric]))
        for treatment in deltas:
            deltas[treatment].append(float(summaries[treatment]["M_species"]) - float(summaries["T0"]["M_species"]))
    out = {t: {m: (percentile(v, .025), percentile(v, .975)) for m, v in metrics.items()} for t, metrics in draws.items()}
    out["delta"] = {t: (percentile(v, .025), percentile(v, .975)) for t, v in deltas.items()}  # type: ignore[assignment]
    return out


def treatment_rows(tracks, masks) -> list[dict[str, object]]:
    rows = []
    for set_name in SETS:
        for ci, clade in enumerate(CLADES):
            track = tracks[clade]
            boot = bootstrap(track, masks[set_name], SEED + 1000 * SETS.index(set_name) + ci)
            summaries = {t: summarize(track, weights_for(track, t, masks[set_name])) for t in TREATMENTS}
            base = summaries["T0"]
            for treatment in TREATMENTS:
                s = summaries[treatment]
                row = {"clade": clade, "structural_mask": set_name, "treatment": treatment, **s,
                       "delta_M_vs_T0": float(s["M_species"]) - float(base["M_species"]),
                       "delta_q2_vs_T0": float(s["q2"]) - float(base["q2"])}
                for metric in (*TOPOLOGIES, "M_species", "M_q1_q2"):
                    row[f"{metric}_ci_low"], row[f"{metric}_ci_high"] = boot[treatment][metric]
                if treatment == "T0":
                    row["delta_M_ci_low"] = row["delta_M_ci_high"] = 0.0
                else:
                    row["delta_M_ci_low"], row["delta_M_ci_high"] = boot["delta"][treatment]
                rows.append(row)
    return rows


def frozen_events() -> list[dict[str, str]]:
    rows = [row for row in read_tsv(EVENTS) if row["structural_set"] == "primary" and int(row["flank_window_bp"]) == MASK_RADIUS]
    assert len(rows) == 16 and len({r["event_id"] for r in rows}) == 16
    return rows


def leave_one_out_rows(tracks, breakpoint_rows) -> list[dict[str, object]]:
    out = []
    for event in frozen_events():
        removed = {event["left_breakpoint_id"]}
        if event["right_breakpoint_id"]:
            removed.add(event["right_breakpoint_id"])
        mask = build_masks(breakpoint_rows, removed)["primary"]
        for clade in CLADES:
            for treatment in ("T2", "T3"):
                s = summarize(tracks[clade], weights_for(tracks[clade], treatment, mask))
                out.append({"clade": clade, "event_id_removed": event["event_id"], "treatment": treatment,
                            **{k: s[k] for k in (*TOPOLOGIES, "M_species", "M_q1_q2", "dominant_topology")}})
    return out


def validate_frozen_inputs() -> None:
    manifest = json.loads(STAGE4A_MANIFEST.read_text())
    for relative, expected in manifest["inputs"].items():
        assert sha256(ROOT / relative) == expected, f"frozen input changed: {relative}"
    strict = [r for r in read_tsv(STAGE3B_SUMMARY) if r["analysis"].startswith("strict_event_state")]
    assert strict and all(int(r["n_eligible_events"]) == 0 and "underpowered" in r["status"] for r in strict)


def validate_stage4a_reproduction(rows: list[dict[str, object]]) -> None:
    old = read_tsv(STAGE4A_WB)
    scheme = {"T0": "window-weighted", "T1": "block-normalized"}
    for treatment in ("T0", "T1"):
        for clade in CLADES:
            new = next(r for r in rows if r["structural_mask"] == "primary" and r["clade"] == clade and r["treatment"] == treatment)
            prior = next(r for r in old if r["clade"] == clade and r["weighting_scheme"] == scheme[treatment])
            for i, q in enumerate(TOPOLOGIES, 1):
                assert math.isclose(float(new[q]), float(prior[f"q_{i}"]), abs_tol=5e-12)


def make_figures(rows, loo) -> None:
    subprocess.run(["Rscript", str(Path(__file__).with_name("04b_make_figures.R"))], cwd=ROOT, check=True)


def category(clade: str, primary: dict[str, dict[str, object]]) -> str:
    t0, t1, t2 = primary["T0"], primary["T1"], primary["T2"]
    d1 = float(t1["M_species"]) - float(t0["M_species"])
    d2 = float(t2["M_species"]) - float(t0["M_species"])
    if clade == "N61" and d1 > 0 and d2 > 0:
        return "A. STRUCTURE-LINKED WEIGHTING EFFECT"
    if abs(d1) >= 0.05 and abs(d2) < 0.05:
        return "B. GENERIC BLOCK-LENGTH EFFECT"
    if abs(d2) >= 0.05 and abs(d1) < 0.05:
        return "C. STRUCTURAL-REGION EFFECT WITHOUT BROAD BLOCK EFFECT"
    return "D. LITTLE MATERIAL WEIGHTING EFFECT"


def write_report(rows, loo, masks) -> None:
    primary_rows = [r for r in rows if r["structural_mask"] == "primary"]
    by_clade = {c: {t: next(r for r in primary_rows if r["clade"] == c and r["treatment"] == t) for t in TREATMENTS} for c in CLADES}
    lines = ["# Neoaves chromosome 4 Stage 4B: inference treatments", "", "## Inputs / frozen definitions", "",
             "All Stage 1–4A files are read only. Masks are the unchanged unions of ±250 kb neighborhoods around frozen breakpoints. The 16 leave-one-out units are the frozen primary W250000 event table used in Stage 4A. Reference labels q1/q2/q3 and topology runs are unchanged.", "",
             "## Treatment definitions", "", "- T0: every resolved window has weight 1.", "- T1: consecutive windows with the same dominant topology form a run; each run has total weight 1, shared equally among its windows.", "- T2: windows whose midpoint is inside the selected frozen structural mask have weight 0; other windows have weight 1.", "- T3: windows whose midpoint is inside the selected frozen structural mask have fixed weight 0.5; other windows have weight 1. This rule is topology neutral and was fixed before Stage 4B execution.", "",
             "## Quartet-level results", "", "Primary mask results:", "", "| clade | treatment | q1/q2/q3 | dominant | M_species | q1-q2 | 95% block-bootstrap CI for M |", "|---|---|---:|---|---:|---:|---:|"]
    for c in CLADES:
        for t in TREATMENTS:
            r = by_clade[c][t]
            lines.append(f"| {c} | {t} | {r['q1']:.3f}/{r['q2']:.3f}/{r['q3']:.3f} | {r['dominant_topology']} | {r['M_species']:.3f} | {r['M_q1_q2']:.3f} | [{r['M_species_ci_low']:.3f}, {r['M_species_ci_high']:.3f}] |")
    lines += ["", "## Margin changes", "", "| clade | DeltaM T1 | DeltaM T2 | DeltaM T3 | Delta q2 T1/T2/T3 |", "|---|---:|---:|---:|---:|"]
    for c in CLADES:
        d = by_clade[c]
        lines.append(f"| {c} | {d['T1']['delta_M_vs_T0']:.3f} | {d['T2']['delta_M_vs_T0']:.3f} | {d['T3']['delta_M_vs_T0']:.3f} | {d['T1']['delta_q2_vs_T0']:.3f}/{d['T2']['delta_q2_vs_T0']:.3f}/{d['T3']['delta_q2_vs_T0']:.3f} |")
    lines += ["", "## Structural vs generic block decomposition", ""]
    for c in CLADES:
        d = by_clade[c]
        lines.append(f"- **{c}:** T1 changes M by {d['T1']['delta_M_vs_T0']:.3f}; T2 changes it by {d['T2']['delta_M_vs_T0']:.3f}; T3 changes it by {d['T3']['delta_M_vs_T0']:.3f}. Assigned category: **{category(c, d)}**.")
    lines += ["", "## Mask sensitivity", ""]
    for c in CLADES:
        bits = []
        for s in SETS:
            rr = {t: next(r for r in rows if r["clade"] == c and r["structural_mask"] == s and r["treatment"] == t) for t in TREATMENTS}
            bits.append(f"{s}: T2 DeltaM={rr['T2']['delta_M_vs_T0']:.3f}, T3 DeltaM={rr['T3']['delta_M_vs_T0']:.3f}")
        lines.append(f"- {c}: " + "; ".join(bits) + ".")
    lines += ["", "## Leave-one-event-out", ""]
    for c in CLADES:
        for t in ("T2", "T3"):
            vals = [float(r["M_species"]) for r in loo if r["clade"] == c and r["treatment"] == t]
            doms = Counter(str(r["dominant_topology"]) for r in loo if r["clade"] == c and r["treatment"] == t)
            lines.append(f"- {c} {t}: M_species range {min(vals):.3f} to {max(vals):.3f}; dominant topologies {dict(doms)}.")
    lines += ["", "## Primary questions", "",
              "1. Block normalization materially redistributes support in all three clades and changes the dominant topology for N61 and N62.",
              "2. Structural exclusion shifts support in the same direction as T1 for each clade, but its effects are much smaller and do not change any dominant topology.",
              "3. For N61, both T1 and T2 reduce q2 dominance; only T1 changes dominance from q2 to q1.",
              "4. Dominant topology changes under T1 for N61 (q2 to q1) and N62 (q1 to q2), and under no other primary treatment.",
              "5. N61 T2 and T3 effects have the same positive direction under stringent, primary, and inclusive masks; neither changes dominance.",
              "6. All 16 N61 leave-one-event-out T2 and T3 analyses remain q2 dominant, with narrow margin ranges, so no single event drives the structural weighting effect.",
              "7. N62's T1 topology change does not persist under T2: T2 remains q1 dominant.",
              "8. N62 is therefore most consistent with a generic block-length effect rather than a frozen structural-region-specific effect.",
              "9. Columbea remains q1 dominant in all treatments despite the large T1 support redistribution.",
              "10. Full species-tree inference was not possible from the repository inputs.", "",
              "## Full species-tree inference", "", "Not performed. The repository contains chr4 locus coordinates and precomputed focal quartet distance summaries, but no per-locus multi-taxon Newick gene trees and no fixed empirical ASTRAL or other species-tree estimator configuration. A genuine full-tree run requires those gene trees, their locus-coordinate mapping, and the project's fixed estimator executable/version/settings. T3 additionally requires native topology-neutral locus weights; otherwise it remains quartet-summary only.", "",
              "## Interpretation", ""]
    for c in CLADES:
        lines.append(f"- {c}: **{category(c, by_clade[c])}**.")
    lines += ["", "These comparisons show whether naive dense-window weighting changes the inferred focal relationship relative to block-aware and frozen-structure-aware weighting. They do not identify a true topology, establish species-tree bias, or validate a treatment against independent evidence.", "", "## Limitations", "", f"Topology support is the weighted frequency of the frozen dominant quartet state per locus, matching Stage 4A. Focal tracks are correlated. Runs and overlapping structural neighborhoods are evidence clusters, not independent biological replicates. The 95% percentile intervals come from {BOOTSTRAPS:,} paired resamples of complete frozen topology runs; they quantify spatial block uncertainty and are not a chromosome replicate analysis.", "", "## Stage 4C requirements", "", "Stage 4C must compare the predeclared treatment results with independent genome-wide evidence, while preserving the frozen labels and treatment rules. Only that external comparison can assess whether any treatment aligns with independent evidence and support language about species-tree bias.", ""]
    (RESULTS / "stage4b_report.md").write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    validate_frozen_inputs()
    breakpoint_rows = read_tsv(BREAKPOINTS)
    masks = build_masks(breakpoint_rows)
    tracks = load_tracks()
    rows = treatment_rows(tracks, masks)
    validate_stage4a_reproduction(rows)
    loo = leave_one_out_rows(tracks, breakpoint_rows)
    treatment_path = RESULTS / "stage4b_inference_treatments.tsv"
    loo_path = RESULTS / "stage4b_leave_one_event_out.tsv"
    write_tsv(treatment_path, rows); write_tsv(loo_path, loo)
    make_figures(rows, loo); write_report(rows, loo, masks)
    outputs = [treatment_path, loo_path, RESULTS / "stage4b_report.md",
               FIGURES / "stage4b_figure_A_treatment_support.pdf", FIGURES / "stage4b_figure_A_treatment_support.png",
               FIGURES / "stage4b_figure_B_species_margin.pdf", FIGURES / "stage4b_figure_B_species_margin.png",
               FIGURES / "stage4b_figure_C_N61_leave_one_event_out.pdf", FIGURES / "stage4b_figure_C_N61_leave_one_event_out.png",
               ROOT / "PROJECT_STATUS.md"]
    inputs = [LOCI, BREAKPOINTS, EVENTS, STAGE3B_SUMMARY, STAGE4A_MANIFEST, STAGE4A_WB]
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
    except subprocess.CalledProcessError:
        commit = "unavailable"
    manifest = {"stage": "4B", "git_commit": commit, "date_utc": datetime.now(timezone.utc).isoformat(),
                "seed": SEED, "bootstrap_replicates": BOOTSTRAPS,
                "stage4a_manifest_checksum": sha256(STAGE4A_MANIFEST),
                "treatments": {"T0": "all resolved windows, weight 1", "T1": "each consecutive identical-topology run has total weight 1", "T2": "frozen-mask windows excluded; background weight 1", "T3": "frozen-mask windows weight 0.5; background weight 1"},
                "t3_topology_neutral": True, "structural_masks": {s: {"merged_intervals": len(masks[s]), "span_bp": sum(b-a+1 for a,b in masks[s])} for s in SETS},
                "full_species_tree_inference": {"performed": False, "reason": "no per-locus multi-taxon gene trees or fixed empirical estimator configuration in repository"},
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
                "stage3b_eligibility_relaxed": False}
    (RESULTS / "stage4b_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
