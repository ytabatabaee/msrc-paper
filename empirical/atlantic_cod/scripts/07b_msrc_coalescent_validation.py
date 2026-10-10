#!/usr/bin/env python3
"""Post-freeze theorem-linked MSRC coalescence-time validation."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import statistics
import sys
import unittest
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/atlantic_cod"
RES = EMP / "results"
FIG = EMP / "figures"
S7_PATH = EMP / "scripts/07_divergence_time_sensitivity.py"
MRCA = ROOT / "data/atlantic_cod/processed/stage7_pairwise_mrca_times.tsv"
LGS = ("LG01", "LG02", "LG07", "LG12")

OUT = {name: RES / f"stage7b_msrc_coalescent_{name}" for name in (
    "window_signal.tsv", "pair_summary.tsv", "lg_summary.tsv", "population_jackknife.tsv",
    "spatial_null.tsv", "topology_comparison.tsv", "parameter_grid.tsv", "report.md",
    "methods_text.md", "results_text.md", "figure_caption.md")}
FIG_PNG = FIG / "atlantic_cod_stage7b_msrc_coalescent.png"
FIG_PDF = FIG / "atlantic_cod_stage7b_msrc_coalescent.pdf"


def load_stage7():
    spec = importlib.util.spec_from_file_location("stage7_divergence_time_sensitivity", S7_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


S7 = load_stage7()


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: "NA" if isinstance(row.get(k), float) and math.isnan(row[k]) else row.get(k, "") for k in fields})


def fmt(x):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.9g}" if isinstance(x, float) else str(x)


def pair_key(r):
    return (str(r["lg"]), str(r["population1"]), str(r["population2"]))


def frozen_mrca():
    return [{**r, "start": int(r["start"]), "end": int(r["end"]), "mrca_time": float(r["mrca_time"])} for r in read_tsv(MRCA)]


def assert_mrca_unchanged():
    expected, _ = S7.build_pairwise_times(S7.read_windows(), S7.read_arrangements())
    actual = frozen_mrca()
    def compact(rows):
        return [(r["lg"], r["window_id"], r["population1"], r["population2"], float(r["mrca_time"])) for r in rows]
    if len(expected) != len(actual):
        raise AssertionError("Frozen Stage-7 MRCA row count changed")
    for a, b in zip(compact(expected), compact(actual)):
        if a[:4] != b[:4] or not math.isclose(a[4], b[4], rel_tol=0, abs_tol=5e-9):
            raise AssertionError(f"Frozen MRCA mismatch for {a[:4]}")
    return actual


def residual_rows(rows, rule="median"):
    outside = defaultdict(list)
    for r in rows:
        if r["region_class"] == "outside":
            outside[pair_key(r)].append(float(r["mrca_time"]))
    baseline = {k: statistics.median(v) if rule == "median" else statistics.fmean(v) for k, v in outside.items()}
    if any(len(v) < 3 for v in outside.values()):
        raise ValueError("A pair has too few fully outside windows")
    return [{**r, "pair_baseline": baseline[pair_key(r)], "residual": float(r["mrca_time"]) - baseline[pair_key(r)], "baseline_rule": rule} for r in rows]


def eligible(rows):
    return [r for r in rows if r["region_class"] != "boundary"]


def contrast(values, mask):
    inside = [v for v, m in zip(values, mask) if m]
    outside = [v for v, m in zip(values, mask) if not m]
    return statistics.fmean(inside) - statistics.fmean(outside) if inside and outside else math.nan


def window_signal(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["lg"], r["window_id"])].append(r)
    out = []
    for (lg, wid), g in sorted(groups.items(), key=lambda x: (x[0][0], x[1][0]["start"])):
        same = [float(r["residual"]) for r in g if r["pair_class"] == "same_arrangement"]
        opp = [float(r["residual"]) for r in g if r["pair_class"] == "opposite_arrangement"]
        f = g[0]
        out.append({"lg": lg, "window_id": wid, "start": f["start"], "end": f["end"], "midpoint": (f["start"] + f["end"]) // 2,
                    "region_class": f["region_class"], "n_same_pairs": len(same), "n_opposite_pairs": len(opp),
                    "mean_residual_same": statistics.fmean(same), "mean_residual_opposite": statistics.fmean(opp),
                    "A_res": statistics.fmean(opp) - statistics.fmean(same)})
    return out


def circular(rows):
    rows = eligible(rows)
    values = [float(r["A_res"]) for r in rows]
    mask = [r["region_class"] == "inside" for r in rows]
    obs = contrast(values, mask)
    null = [contrast(S7.rotate(values, s), mask) for s in range(len(values))]
    return obs, sum(v >= obs for v in null) / len(null), len(null)


def pair_summary(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[pair_key(r)].append(r)
    out = []
    for k, g in sorted(groups.items()):
        inside = [float(r["residual"]) for r in g if r["region_class"] == "inside"]
        outside = [float(r["residual"]) for r in g if r["region_class"] == "outside"]
        if not inside or not outside:
            continue
        q = statistics.quantiles(inside, n=4, method="inclusive") if len(inside) > 1 else [inside[0]] * 3
        f = g[0]
        out.append({"lg": k[0], "population1": k[1], "population2": k[2], "arrangement1": f["arrangement1"], "arrangement2": f["arrangement2"], "pair_class": f["pair_class"], "baseline": f["pair_baseline"],
                    "residual_inside_mean": statistics.fmean(inside), "residual_inside_median": statistics.median(inside), "residual_inside_iqr": q[2] - q[0],
                    "residual_outside_mean": statistics.fmean(outside), "residual_outside_median": statistics.median(outside), "fraction_positive_inside": sum(x > 0 for x in inside) / len(inside), "n_inside": len(inside), "n_outside": len(outside)})
    return out


def class_summary(pairs):
    out = []
    for (lg, cls), g in sorted(((k, list(v)) for k, v in __import__("itertools").groupby(sorted(pairs, key=lambda r: (r["lg"], r["pair_class"])), key=lambda r: (r["lg"], r["pair_class"]))), key=lambda x: x[0]):
        vals = [float(r["residual_inside_mean"]) for r in g]
        q = statistics.quantiles(vals, n=4, method="inclusive") if len(vals) > 1 else [vals[0]] * 3
        out.append({"lg": lg, "pair_class": cls, "n_pairs": len(vals), "mean_inside_residual": statistics.fmean(vals), "median_inside_residual": statistics.median(vals), "iqr_inside_residual": q[2] - q[0], "fraction_pairs_positive_inside": sum(v > 0 for v in vals) / len(vals)})
    return out


def lg_summary(windows, pairs, rule):
    out = []
    raw = {}
    for lg in LGS:
        w = [r for r in windows if r["lg"] == lg]
        e = eligible(w)
        obs, p, nshift = circular(w)
        inside = [float(r["A_res"]) for r in e if r["region_class"] == "inside"]
        outside = [float(r["A_res"]) for r in e if r["region_class"] == "outside"]
        g = [r for r in pairs if r["lg"] == lg]
        def vals(cls): return [float(r["residual_inside_mean"]) for r in g if r["pair_class"] == cls]
        same, opp = vals("same_arrangement"), vals("opposite_arrangement")
        out.append({"lg": lg, "baseline_rule": rule, "n_inside_windows": len(inside), "n_outside_windows": len(outside), "mean_A_res_inside": statistics.fmean(inside), "mean_A_res_outside": statistics.fmean(outside), "delta_coal": obs, "circular_p": p, "n_circular_shifts": nshift, "m_eff": 1 / (2 * obs) if obs > 0 else math.nan, "same_inside_residual_mean": statistics.fmean(same), "same_inside_residual_median": statistics.median(same), "opposite_inside_residual_mean": statistics.fmean(opp), "opposite_inside_residual_median": statistics.median(opp)})
        raw[lg] = p
    bh = S7.bh_adjust(raw)
    for r in out: r["circular_BH_p"] = bh[r["lg"]]
    return out


def delta_T(m, lam, tau):
    if m <= 0 or lam <= 0 or tau < 0:
        raise ValueError("m and lambda must be positive; tau must be nonnegative")
    r = math.sqrt(lam * lam + 16 * m * m)
    a, b = (lam + 4 * m) / 2, r / 2
    bracket = .5 * (1 + lam / r) * math.exp(-(a - b) * tau) + .5 * (1 - lam / r) * math.exp(-(a + b) * tau)
    return (1 / (2 * m)) * (1 - bracket)


def delta_T_dimensionless(mu, x):
    return delta_T(mu, 1.0, x)


def jackknife(rows, populations):
    all_windows = window_signal(rows)
    all_pairs = pair_summary(rows)
    full = {r["lg"]: r["delta_coal"] for r in lg_summary(all_windows, all_pairs, "median")}
    output = []
    for lg in LGS:
        estimates = []
        for pop in populations:
            kept = [r for r in rows if r["population1"] != pop and r["population2"] != pop]
            status = "defined"
            try:
                estimate = next(r["delta_coal"] for r in lg_summary(window_signal(kept), pair_summary(kept), "median") if r["lg"] == lg)
                estimates.append(estimate)
            except (StopIteration, statistics.StatisticsError, ValueError, ZeroDivisionError):
                estimate, status = math.nan, "undefined_after_population_removal"
            output.append({"lg": lg, "removed_population": pop, "delta_coal_leave_one_out": estimate, "status": status})
        defined = [r["delta_coal_leave_one_out"] for r in output if r["lg"] == lg and r["status"] == "defined"]
        mean = statistics.fmean(defined) if defined else math.nan
        se = math.sqrt((len(populations) - 1) / len(populations) * sum((x - mean) ** 2 for x in defined)) if len(defined) == len(populations) else math.nan
        for r in output:
            if r["lg"] == lg:
                r.update({"full_delta_coal": full[lg], "jackknife_n": len(populations), "jackknife_mean": mean, "jackknife_min": min(defined) if defined else math.nan, "jackknife_max": max(defined) if defined else math.nan, "jackknife_SE": se})
    return output


def physical_null(windows):
    regions = S7.read_regions()
    out = []
    for lg in LGS:
        ordered = eligible([r for r in windows if r["lg"] == lg])
        values = [float(r["A_res"]) for r in ordered]
        mask = [r["region_class"] == "inside" for r in ordered]
        observed = contrast(values, mask)
        for shift in range(len(values)):
            value = contrast(S7.rotate(values, shift), mask)
            out.append({"lg": lg, "null_type": "circular_shift", "shift_or_candidate": shift, "delta_coal_null": value,
                        "greater_or_equal_observed": value >= observed, "physical_p": math.nan})
        compatible = [{**r, "A_opposite_minus_same": r["A_res"]} for r in windows if r["lg"] == lg]
        p, rows = S7.physical_test_for_lg(compatible, regions[lg])
        out.extend({"lg": lg, "null_type": "physical_coordinate", "shift_or_candidate": r["candidate_index"], "delta_coal_null": r["delta_A"],
                    "greater_or_equal_observed": r["greater_or_equal_observed"], "physical_p": p} for r in rows)
    return out


def topology_comparison(lg_rows):
    old = {r["lg"]: r for r in read_tsv(RES / "stage4b_circular_shift_summary.tsv")}
    out = []
    for r in lg_rows:
        lg = r["lg"]; dd = float(old[lg]["observed_delta_D"]); tp = float(old[lg]["p_one_sided"]); dc = float(r["delta_coal"]); cp = float(r["circular_p"])
        cls = "topology_and_time_shift" if dd > 0 and tp <= .05 and dc > 0 and cp <= .05 else "time_shift_without_topology_shift" if dc > 0 and cp <= .05 and (dd <= 0 or tp > .05) else "weak_or_no_shift"
        out.append({"lg": lg, "delta_D": dd, "topology_p": tp, "delta_coal": dc, "coalescent_circular_p": cp, "coalescent_BH_p": r["circular_BH_p"], "qualitative_class": cls})
    return out


def grid():
    return [{"mu": mu, "x": x, "scaled_delta_lambda_T": delta_T_dimensionless(mu, x), "interpretation": "lambda_times_delta_T; identifiability_grid_only"}
            for mu in (.01, .03, .1, .3, 1., 3., 10.) for x in (.1, .3, 1., 3., 10., 30.)]


def plot(pairs, jack, topo):
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), gridspec_kw={"width_ratios": [1.25, 1.15, 1]})
    colors = {"same_arrangement": "#4c78a8", "opposite_arrangement": "#e45756"}
    for i, lg in enumerate(LGS):
        for cls, color in colors.items():
            v = [float(r["residual_inside_mean"]) for r in pairs if r["lg"] == lg and r["pair_class"] == cls]
            x = [i + (-.16 if cls.startswith("same") else .16) + (j % 5 - 2) * .018 for j in range(len(v))]
            ax[0].scatter(x, v, s=16, alpha=.7, color=color, label=cls.replace("_", " ") if i == 0 else None)
    ax[0].axhline(0, color="black", lw=.7); ax[0].set_xticks(range(4), LGS); ax[0].set_ylabel("inside residual MRCA time"); ax[0].set_title("A  Pair-level shifts"); ax[0].legend(frameon=False, fontsize=8)
    for i, lg in enumerate(LGS):
        v = [float(r["delta_coal_leave_one_out"]) for r in jack if r["lg"] == lg and r["status"] == "defined"]
        ax[1].scatter([i] * len(v), v, s=13, color="#777777", alpha=.65)
        ax[1].scatter([i], [next(r["full_delta_coal"] for r in jack if r["lg"] == lg)], s=70, color="#111111", zorder=3)
    ax[1].axhline(0, color="black", lw=.7); ax[1].set_xticks(range(4), LGS); ax[1].set_ylabel("Delta_coal"); ax[1].set_title("B  Population jackknife")
    for r in topo:
        ax[2].scatter(float(r["delta_D"]), float(r["delta_coal"]), s=45, color="#e45756" if r["lg"] == "LG07" else "#4c78a8"); ax[2].annotate(r["lg"], (float(r["delta_D"]), float(r["delta_coal"])), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax[2].axhline(0, color="black", lw=.7); ax[2].axvline(0, color="black", lw=.7); ax[2].set_xlabel("Delta_D"); ax[2].set_ylabel("Delta_coal"); ax[2].set_title("C  Topology vs time")
    fig.tight_layout(); fig.savefig(FIG_PDF); fig.savefig(FIG_PNG, dpi=300); plt.close(fig)


def documents(lg, topo, classes, jack):
    OUT["methods_text.md"].write_text("""This post-freeze theorem-linked extension reused the frozen 250-kb SNAPP MCC window trees, inversion coordinates, population arrangement states, and Stage-7 MRCA pipeline. The Stage-7 MRCA values were verified unchanged. For each linkage group and population pair, the primary baseline was the median MRCA time across fully outside windows; residuals were MRCA time minus that baseline. Boundary-overlap windows were excluded. A_res(w) is mean residual for opposite-arrangement pairs minus mean residual for same-arrangement pairs, and Delta_coal is its inside-minus-outside mean. The exact chromosome-aware circular-shift null preserved track order, kept the inversion mask fixed, included the observed alignment, and used the one-sided predeclared direction. BH correction was across LG01, LG02, LG07, and LG12. Mean-outside baselines, delete-one-population jackknife, and physical-coordinate placement were sensitivities. Published MCC point estimates were used because per-window posterior samples were unavailable. No SNAPP, BEAST, ASTRAL, or local gene-tree inference was run.

The finite-duration MSRC formula was implemented as delta_T(m, lambda, tau), with dimensionless mu=m/lambda and x=lambda*tau. The parameter grid is descriptive only: one observed contrast cannot separately identify m, lambda, and tau.
""")
    lines = ["# Stage7b MSRC coalescent validation", "", "The predeclared theorem-linked direction is `Delta_coal > 0`: opposite-arrangement pairs have extra residual coalescence depth inside the supergene after pair-specific collinear baselines are removed.", "", "| LG | Delta_coal | circular p | BH p | m_eff |", "|---|---:|---:|---:|---:|"]
    for r in lg: lines.append(f"| {r['lg']} | {float(r['delta_coal']):.6g} | {float(r['circular_p']):.6g} | {float(r['circular_BH_p']):.6g} | {fmt(r['m_eff'])} |")
    lines += ["", "m_eff is only the effective long-duration approximation 1/(2 Delta_coal), in inverse published SNAPP time units. It is not a recombination rate and does not identify biological m.", "", "The circular null is chromosome-aware; the jackknife uses populations as dependence units and does not treat population pairs as iid. LG07 is the key example of near-zero topology shift with positive coalescence-time shift. This is a theory-consistent empirical signature, not proof that MSRC is the sole historical process.", "", "## Arrangement-class pair summaries", "", "| LG | class | n pairs | mean inside residual | median inside residual | IQR | fraction positive |", "|---|---|---:|---:|---:|---:|---:|"]
    for r in classes: lines.append(f"| {r['lg']} | {r['pair_class']} | {r['n_pairs']} | {float(r['mean_inside_residual']):.6g} | {float(r['median_inside_residual']):.6g} | {float(r['iqr_inside_residual']):.6g} | {float(r['fraction_pairs_positive_inside']):.6g} |")
    lines += ["", "The mean-outside baseline sensitivity preserved the sign and qualitative LG ordering of Delta_coal; because the same pair set contributes to every eligible window, subtracting a different constant baseline per pair changes absolute residual levels but leaves this inside-minus-outside contrast unchanged up to rounding.", "", "| LG | Delta_D | topology p | Delta_coal | coalescent p | BH p | class |", "|---|---:|---:|---:|---:|---:|---|"]
    for r in topo: lines.append(f"| {r['lg']} | {float(r['delta_D']):.6g} | {float(r['topology_p']):.6g} | {float(r['delta_coal']):.6g} | {float(r['coalescent_circular_p']):.6g} | {float(r['coalescent_BH_p']):.6g} | {r['qualitative_class']} |")
    lines += ["", "The data do not separately identify m, lambda, and tau under the finite-duration model. The analysis uses published per-window SNAPP MCC point estimates because posterior samples are unavailable."]
    OUT["report.md"].write_text("\n".join(lines) + "\n")
    result = ["After pair-specific median-outside correction, opposite-arrangement residual depth exceeded same-arrangement residual depth in the primary inside-versus-outside contrast for: " + ", ".join(r["lg"] for r in lg if float(r["delta_coal"]) > 0) + ".", "", "| LG | same mean inside residual | opposite mean inside residual | jackknife range |", "|---|---:|---:|---:|"]
    for r in lg:
        j = [x for x in jack if x["lg"] == r["lg"] and x["status"] == "defined"]
        result.append(f"| {r['lg']} | {float(r['same_inside_residual_mean']):.6g} | {float(r['opposite_inside_residual_mean']):.6g} | {min(float(x['delta_coal_leave_one_out']) for x in j):.6g} to {max(float(x['delta_coal_leave_one_out']) for x in j):.6g} |")
    OUT["results_text.md"].write_text("\n".join(result) + "\n")
    OUT["figure_caption.md"].write_text("Figure. Theorem-linked Atlantic cod coalescence-time validation. Panel A shows population-pair inside residual MRCA times after subtracting each pair's fully outside baseline. Panel B shows the full Delta_coal and delete-one-population estimates. Panel C combines the existing topology contrast Delta_D with Delta_coal; LG07 illustrates positive time distortion with little topology distortion. Published SNAPP MCC point estimates and a chromosome-aware circular-shift null were used.\n")


def main():
    frozen = assert_mrca_unchanged()
    populations = sorted({r["population1"] for r in frozen} | {r["population2"] for r in frozen})
    primary = residual_rows(frozen, "median")
    mean_rows = residual_rows(frozen, "mean")
    primary_windows, mean_windows = window_signal(primary), window_signal(mean_rows)
    primary_pairs, mean_pairs = pair_summary(primary), pair_summary(mean_rows)
    lg = lg_summary(primary_windows, primary_pairs, "median")
    mean_lg = lg_summary(mean_windows, mean_pairs, "mean")
    jack = jackknife(primary, populations)
    spatial = physical_null(primary_windows)
    topo = topology_comparison(lg)
    classes = class_summary(primary_pairs)
    write_tsv(OUT["window_signal.tsv"], primary_windows, ["lg", "window_id", "start", "end", "midpoint", "region_class", "n_same_pairs", "n_opposite_pairs", "mean_residual_same", "mean_residual_opposite", "A_res"])
    write_tsv(OUT["pair_summary.tsv"], primary_pairs, ["lg", "population1", "population2", "arrangement1", "arrangement2", "pair_class", "baseline", "residual_inside_mean", "residual_inside_median", "residual_inside_iqr", "residual_outside_mean", "residual_outside_median", "fraction_positive_inside", "n_inside", "n_outside"])
    write_tsv(OUT["lg_summary.tsv"], lg + [{**r, "baseline_rule": "mean"} for r in mean_lg], ["lg", "baseline_rule", "n_inside_windows", "n_outside_windows", "mean_A_res_inside", "mean_A_res_outside", "delta_coal", "circular_p", "circular_BH_p", "n_circular_shifts", "m_eff", "same_inside_residual_mean", "same_inside_residual_median", "opposite_inside_residual_mean", "opposite_inside_residual_median"])
    write_tsv(OUT["population_jackknife.tsv"], jack, ["lg", "removed_population", "delta_coal_leave_one_out", "status", "full_delta_coal", "jackknife_n", "jackknife_mean", "jackknife_min", "jackknife_max", "jackknife_SE"])
    write_tsv(OUT["spatial_null.tsv"], spatial, ["lg", "null_type", "shift_or_candidate", "delta_coal_null", "greater_or_equal_observed", "physical_p"])
    write_tsv(OUT["topology_comparison.tsv"], topo, ["lg", "delta_D", "topology_p", "delta_coal", "coalescent_circular_p", "coalescent_BH_p", "qualitative_class"])
    write_tsv(OUT["parameter_grid.tsv"], grid(), ["mu", "x", "scaled_delta_lambda_T", "interpretation"])
    documents(lg, topo, classes, jack)
    plot(primary_pairs, jack, topo)
    print("Delta_coal for LG01, LG02, LG07, LG12:")
    for r, mr in zip(lg, mean_lg):
        j = [x for x in jack if x["lg"] == r["lg"] and x["status"] == "defined"]
        print(f"{r['lg']} delta_coal={r['delta_coal']:.9g} circular_p={r['circular_p']:.9g} BH_p={r['circular_BH_p']:.9g}; primary median-baseline same inside median/mean={r['same_inside_residual_median']:.9g}/{r['same_inside_residual_mean']:.9g}, opposite={r['opposite_inside_residual_median']:.9g}/{r['opposite_inside_residual_mean']:.9g}; mean-baseline same inside median/mean={mr['same_inside_residual_median']:.9g}/{mr['same_inside_residual_mean']:.9g}, opposite={mr['opposite_inside_residual_median']:.9g}/{mr['opposite_inside_residual_mean']:.9g}; jackknife range={min(x['delta_coal_leave_one_out'] for x in j):.9g}..{max(x['delta_coal_leave_one_out'] for x in j):.9g}; mean-baseline delta={mr['delta_coal']:.9g}; m_eff={fmt(r['m_eff'])}")
    print("Delta_D versus Delta_coal table:")
    for r in topo: print(f"{r['lg']} delta_D={r['delta_D']:.9g} delta_coal={r['delta_coal']:.9g} topology_p={r['topology_p']:.9g} coalescent_p={r['coalescent_circular_p']:.9g} BH={r['coalescent_BH_p']:.9g} class={r['qualitative_class']}")
    print("Outputs:")
    for p in (*OUT.values(), FIG_PNG, FIG_PDF): print(p)
    print("Existing Stage-7 files were unchanged; no SNAPP/BEAST/ASTRAL/gene-tree inference was rerun.")


class Stage7bTests(unittest.TestCase):
    def test_mrca_values_unchanged(self): assert_mrca_unchanged()

    def test_primary_baseline_is_median_outside(self):
        rows = frozen_mrca(); keyed = pair_key(rows[0]); values = [r["mrca_time"] for r in rows if pair_key(r) == keyed and r["region_class"] == "outside"]
        primary = residual_rows(rows, "median")
        self.assertEqual(primary[0]["pair_baseline"], statistics.median(values))

    def test_boundary_and_residual_contrast(self):
        self.assertEqual(len(eligible([{"region_class": "inside"}, {"region_class": "boundary"}, {"region_class": "outside"}])), 2)
        rows = residual_rows(frozen_mrca()); windows = window_signal(rows); row = next(r for r in rows if r["region_class"] == "outside")
        self.assertAlmostEqual(row["residual"], row["mrca_time"] - row["pair_baseline"])
        w = next(r for r in windows if r["region_class"] != "boundary")
        self.assertAlmostEqual(w["A_res"], w["mean_residual_opposite"] - w["mean_residual_same"])

    def test_delta_inside_minus_outside(self):
        rows = residual_rows(frozen_mrca()); windows = window_signal(rows); pairs = pair_summary(rows); got = next(r for r in lg_summary(windows, pairs, "median") if r["lg"] == "LG01")
        e = eligible([r for r in windows if r["lg"] == "LG01"])
        self.assertAlmostEqual(got["delta_coal"], contrast([r["A_res"] for r in e], [r["region_class"] == "inside" for r in e]))

    def test_jackknife_is_population_level(self):
        rows = residual_rows(frozen_mrca()); pops = sorted({r["population1"] for r in rows} | {r["population2"] for r in rows}); result = jackknife(rows, pops)
        self.assertEqual(len([r for r in result if r["lg"] == "LG01"]), len(pops))
        self.assertEqual(next(r for r in result if r["lg"] == "LG01")["jackknife_n"], len(pops))

    def test_circular_shift_preserves_track_length_and_mask_count(self):
        rows = [r for r in window_signal(residual_rows(frozen_mrca())) if r["lg"] == "LG01"]; e = eligible(rows); values = [r["A_res"] for r in e]; mask = [r["region_class"] == "inside" for r in e]
        self.assertEqual(len(S7.rotate(values, 7)), len(values)); self.assertEqual(sum(mask), len([r for r in e if r["region_class"] == "inside"]))

    def test_bh_across_four_lgs(self): self.assertEqual(set(S7.bh_adjust({lg: .01 for lg in LGS})), set(LGS))

    def test_finite_duration_curve(self):
        self.assertGreater(delta_T(.1, 1., 1.), 0); self.assertAlmostEqual(delta_T(.1, 1., 0), 0, places=12); self.assertAlmostEqual(delta_T(.1, 1., 1e4), 5., places=8)

    def test_no_upstream_inference_invocation(self):
        text = Path(__file__).read_text().lower(); self.assertNotIn("sub" + "process", text); self.assertNotIn("os." + "system(", text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--run-tests", action="store_true"); args = parser.parse_args()
    if args.run_tests:
        raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage7bTests)).wasSuccessful() else 1)
    main()
