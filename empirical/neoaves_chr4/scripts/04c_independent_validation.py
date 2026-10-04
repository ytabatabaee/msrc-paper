#!/usr/bin/env python3
"""Stage 4C: validate frozen chr4 treatments against non-chr4 evidence."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import lzma
import math
import random
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/neoaves_chr4/raw/genetreesupport"
RESULTS = ROOT / "empirical/neoaves_chr4/results"
FIGURES = ROOT / "empirical/neoaves_chr4/figures"
META = RAW / "63K_trees.names_header.txt.xz"
REC = RAW / "clade-rec.stat.xz"
RAW_README = RAW / "README.md"
CLADE_README = RAW / "clade-analysis/README.md"
DRAW_SCRIPT = RAW / "draw-movingaverage.r"
STAGE4B = RESULTS / "stage4b_inference_treatments.tsv"
STAGE4B_MANIFEST = RESULTS / "stage4b_manifest.json"
STAGE3B_ROLES = RESULTS / "stage3b_v2_quartet_role_mapping.tsv"
STAGE4B_SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/04b_inference_treatments.py"

CLADES = ("Columbea", "N61", "N62")
TOPOLOGIES = ("q1", "q2", "q3")
TREATMENTS = ("T0", "T1", "T2", "T3")
FOCAL_MAP = {
    ("Columbimorphae", "Phoenicopteriformes"): "Columbea",
    ("Otidimorphae", "Columbimorphae"): "N61",
    ("Columbiformes", "OtherColumbimorphae"): "N62",
}
REFERENCE_ID = "nonchr4_genomewide_QQS_v1"
SEED = 431729
BOOTSTRAPS = 10_000


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


def chromosome_group(label: str) -> str:
    stem = label.split("_", 1)[0]
    return stem[3:] if stem.startswith("chr") else stem


def is_chr4(label: str) -> bool:
    return chromosome_group(label) == "4"


def is_main_chromosome(label: str) -> bool:
    return re.fullmatch(r"chr(?:[0-9]+|W|Z)", label) is not None


def compute_qqs(x: dict[str, float]) -> dict[str, float] | None:
    if set(x) != {"1", "2", "3", "4"}:
        return None
    d1, d2, d3 = x["1"] - x["4"], x["2"] - x["4"], x["3"] - x["4"]
    denominator = d1 + d2 + d3
    if denominator == 0 or not math.isfinite(denominator):
        return None
    return {"q1": (d2 + d3 - d1) / denominator,
            "q2": (d1 + d3 - d2) / denominator,
            "q3": (d1 + d2 - d3) / denominator}


def load_genomewide_rows() -> list[dict[str, object]]:
    metadata: dict[str, dict[str, object]] = {}
    with lzma.open(META, "rt") as handle:
        for row in csv.DictReader(handle, delimiter=" ", skipinitialspace=True):
            metadata[row["Gene"]] = {"chromosome": row["Chromosome"], "start": int(row["ws"]),
                                      "end": int(row["we"]), "gene": row["Gene"]}
    values: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    with lzma.open(REC, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue
            c1, c2, topology_id, gene, value = line.split()
            clade = FOCAL_MAP.get((c1, c2))
            if clade and gene in metadata:
                values[(gene, clade)][topology_id] = float(value)
    rows = []
    for (gene, clade), x in values.items():
        q = compute_qqs(x)
        if q is None:
            continue
        dominant = max(TOPOLOGIES, key=lambda topology: (q[topology], topology))
        meta = metadata[gene]
        rows.append({**meta, "clade": clade, **q, "topology": dominant,
                     "chromosome_group": chromosome_group(str(meta["chromosome"])),
                     "main_chromosome": is_main_chromosome(str(meta["chromosome"]))})
    rows.sort(key=lambda r: (CLADES.index(str(r["clade"])), str(r["chromosome"]), int(r["start"]), int(r["end"]), int(r["gene"])))
    return rows


def state_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    counts = Counter(str(row["topology"]) for row in rows)
    n = sum(counts.values())
    q = {topology: counts[topology] / n for topology in TOPOLOGIES}
    dominant = max(TOPOLOGIES, key=lambda topology: (q[topology], topology))
    margin = q[dominant] - max(q[t] for t in TOPOLOGIES if t != dominant)
    return {**q, "dominant": dominant, "margin": margin, "n_loci": n}


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low, high = math.floor(position), math.ceil(position)
    return ordered[low] if low == high else ordered[low] * (high - position) + ordered[high] * (position - low)


def chromosome_bootstrap(rows: list[dict[str, object]], seed: int) -> dict[str, tuple[float, float]]:
    groups = defaultdict(Counter)
    for row in rows:
        groups[str(row["chromosome_group"])][str(row["topology"])] += 1
    names = sorted(groups)
    rng = random.Random(seed)
    draws = {metric: [] for metric in (*TOPOLOGIES, "margin")}
    reference_topology = str(state_summary(rows)["dominant"])
    for _ in range(BOOTSTRAPS):
        sampled_names = Counter(rng.choices(names, k=len(names)))
        counts = Counter()
        for name, multiplicity in sampled_names.items():
            for topology in TOPOLOGIES:
                counts[topology] += multiplicity * groups[name][topology]
        total = sum(counts.values())
        q = {topology: counts[topology] / total for topology in TOPOLOGIES}
        for topology in TOPOLOGIES:
            draws[topology].append(q[topology])
        draws["margin"].append(q[reference_topology] - max(q[t] for t in TOPOLOGIES if t != reference_topology))
    return {metric: (percentile(values, .025), percentile(values, .975)) for metric, values in draws.items()}


def reference_results(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    out = {}
    for index, clade in enumerate(CLADES):
        selected = [row for row in rows if row["clade"] == clade and not is_chr4(str(row["chromosome"]))]
        summary = state_summary(selected)
        ci = chromosome_bootstrap(selected, SEED + index)
        out[clade] = {**summary,
                      **{f"{metric}_ci_low": bounds[0] for metric, bounds in ci.items()},
                      **{f"{metric}_ci_high": bounds[1] for metric, bounds in ci.items()},
                      "n_chromosome_units": len({str(row["chromosome_group"]) for row in selected})}
    return out


def reference_margin(q: dict[str, float], topology: str) -> float:
    return q[topology] - max(q[t] for t in TOPOLOGIES if t != topology)


def treatment_reference_bootstrap(clade: str, reference_topology: str) -> dict[str, dict[str, tuple[float, float]]]:
    """Repeat the frozen Stage 4B run bootstrap for an arbitrary reference."""
    spec = importlib.util.spec_from_file_location("stage4b_for_stage4c", STAGE4B_SCRIPT)
    stage4b = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(stage4b)
    track = stage4b.load_tracks()[clade]
    mask = stage4b.build_masks(stage4b.read_tsv(stage4b.BREAKPOINTS))["primary"]
    runs = stage4b.make_runs(track)
    contributions = []
    for run in runs:
        associated = sum(stage4b.in_mask(float(row["midpoint"]), mask) for row in run)
        background = len(run) - associated
        contributions.append({"topology": str(run[0]["topology"]), "T0": float(len(run)), "T1": 1.0,
                              "T2": float(background), "T3": float(background) + .5 * float(associated)})
    seed = stage4b.SEED + 1000 + CLADES.index(clade)
    rng = random.Random(seed)
    population = list(range(len(contributions)))
    margins = {t: [] for t in TREATMENTS}
    deltas = {t: [] for t in TREATMENTS}
    for _ in range(stage4b.BOOTSTRAPS):
        multiplicities = Counter(rng.choices(population, k=len(population)))
        replicate = {}
        for treatment in TREATMENTS:
            counts = {q: 0.0 for q in TOPOLOGIES}
            for index, multiplicity in multiplicities.items():
                run = contributions[index]
                counts[str(run["topology"])] += multiplicity * float(run[treatment])
            total = sum(counts.values())
            q = {topology: counts[topology] / total for topology in TOPOLOGIES}
            replicate[treatment] = reference_margin(q, reference_topology)
            margins[treatment].append(replicate[treatment])
        for treatment in TREATMENTS:
            deltas[treatment].append(replicate[treatment] - replicate["T0"])
    return {t: {"margin": (percentile(margins[t], .025), percentile(margins[t], .975)),
                "delta": (percentile(deltas[t], .025), percentile(deltas[t], .975))} for t in TREATMENTS}


def primary_rows(references: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    stage4b = [row for row in read_tsv(STAGE4B) if row["structural_mask"] == "primary"]
    out = []
    categories = {}
    for clade in CLADES:
        ref = references[clade]
        ref_top = str(ref["dominant"])
        treatments = {t: next(row for row in stage4b if row["clade"] == clade and row["treatment"] == t) for t in TREATMENTS}
        margins = {t: reference_margin({q: float(treatments[t][q]) for q in TOPOLOGIES}, ref_top) for t in TREATMENTS}
        boot = treatment_reference_bootstrap(clade, ref_top)
        if treatments["T0"]["dominant_topology"] != ref_top and any(margins[t] > margins["T0"] for t in ("T1", "T2", "T3")) and float(ref["margin_ci_low"]) > 0:
            category = "A — EMPIRICAL BIAS / PARTIAL CORRECTION"
        elif treatments["T0"]["dominant_topology"] == ref_top and any(treatments[t]["dominant_topology"] != ref_top for t in ("T1", "T2", "T3")):
            category = "D — CORRECTION MOVES AWAY FROM REFERENCE"
        elif treatments["T0"]["dominant_topology"] == ref_top:
            category = "C — NAIVE CHR4 AGREES WITH REFERENCE"
        else:
            category = "B — WEIGHTING SENSITIVITY WITHOUT VALIDATED BIAS"
        categories[clade] = category
        for treatment in TREATMENTS:
            row = treatments[treatment]
            q = {topology: float(row[topology]) for topology in TOPOLOGIES}
            m_low, m_high = boot[treatment]["margin"]
            d_low, d_high = boot[treatment]["delta"]
            out.append({
                "clade": clade, "reference_id": REFERENCE_ID,
                "reference_source": "63K locus metadata + clade-rec QQS; all chr4-labelled sequences excluded",
                "reference_topology": ref_top,
                **{f"reference_{q}": ref[q] for q in TOPOLOGIES}, "reference_margin": ref["margin"],
                **{f"reference_{q}_ci_low": ref[f"{q}_ci_low"] for q in TOPOLOGIES},
                **{f"reference_{q}_ci_high": ref[f"{q}_ci_high"] for q in TOPOLOGIES},
                "reference_margin_ci_low": ref["margin_ci_low"], "reference_margin_ci_high": ref["margin_ci_high"],
                "reference_n_loci": ref["n_loci"], "reference_n_chromosome_units": ref["n_chromosome_units"],
                "treatment": treatment, **{f"treatment_{topology}": q[topology] for topology in TOPOLOGIES},
                "treatment_dominant": row["dominant_topology"], "M_reference": margins[treatment],
                "M_reference_ci_low": m_low, "M_reference_ci_high": m_high,
                "delta_M_reference_vs_T0": margins[treatment] - margins["T0"],
                "delta_M_reference_ci_low": d_low, "delta_M_reference_ci_high": d_high,
                "matches_reference": row["dominant_topology"] == ref_top,
                "interpretation_category": category,
            })
    return out


def make_runs(rows: list[dict[str, object]]) -> list[list[dict[str, object]]]:
    runs: list[list[dict[str, object]]] = []
    for row in sorted(rows, key=lambda r: (int(r["start"]), int(r["end"]), int(r["gene"]))):
        if not runs or runs[-1][-1]["topology"] != row["topology"]:
            runs.append([row])
        else:
            runs[-1].append(row)
    return runs


def control_rows(rows: list[dict[str, object]], references: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    chromosomes = sorted({str(row["chromosome"]) for row in rows if bool(row["main_chromosome"])}, key=lambda x: (chromosome_group(x) not in {str(i) for i in range(1, 100)}, int(chromosome_group(x)) if chromosome_group(x).isdigit() else 999, x))
    out = []
    for chromosome in chromosomes:
        for clade in CLADES:
            selected = [row for row in rows if row["chromosome"] == chromosome and row["clade"] == clade]
            if not selected:
                continue
            window = state_summary(selected)
            runs = make_runs(selected)
            blocks = state_summary([{"topology": run[0]["topology"]} for run in runs])
            ref_top = str(references[clade]["dominant"])
            window_q = {q: float(window[q]) for q in TOPOLOGIES}
            block_q = {q: float(blocks[q]) for q in TOPOLOGIES}
            out.append({"chromosome": chromosome, "focal_clade": clade,
                        **{f"window_{q}": window_q[q] for q in TOPOLOGIES},
                        **{f"block_{q}": block_q[q] for q in TOPOLOGIES},
                        "window_dominant": window["dominant"], "block_dominant": blocks["dominant"],
                        "M_reference_window": reference_margin(window_q, ref_top),
                        "M_reference_block": reference_margin(block_q, ref_top),
                        "delta_M_reference_block_vs_window": reference_margin(block_q, ref_top) - reference_margin(window_q, ref_top),
                        "topology_changed": window["dominant"] != blocks["dominant"],
                        "n_loci": len(selected), "n_blocks": len(runs),
                        "reference_topology": ref_top})
    return out


def jackknife_rows(rows: list[dict[str, object]], references: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    out = []
    for clade in CLADES:
        base = [row for row in rows if row["clade"] == clade and not is_chr4(str(row["chromosome"]))]
        groups = sorted({str(row["chromosome_group"]) for row in base})
        for group in groups:
            summary = state_summary([row for row in base if row["chromosome_group"] != group])
            ref_top = str(references[clade]["dominant"])
            q = {topology: float(summary[topology]) for topology in TOPOLOGIES}
            out.append({"clade": clade, "chromosome_unit_removed": group,
                        **q, "dominant_topology": summary["dominant"],
                        "M_reference": reference_margin(q, ref_top),
                        "reference_topology": ref_top, "n_loci": summary["n_loci"]})
    return out


def candidate_rows() -> list[dict[str, object]]:
    return [
        {"reference_id": REFERENCE_ID, "path": f"{META.relative_to(ROOT)};{REC.relative_to(ROOT)}", "taxa": "48-taxon source trees summarized into frozen focal quadripartitions", "genomic_coverage": "63,430 loci across genome labels", "chr4_in_source": True, "coordinates_available": True, "inference_method": "published QQS from quartet distances to four quadripartitions", "predates_stage4c": True, "independent_of_chr4": "yes after prospective exclusion of every chr4-labelled sequence", "available": True, "hierarchy_rank": "F (genome-wide quartet summaries)"},
        {"reference_id": "published_q1_branch_definitions", "path": f"{DRAW_SCRIPT.relative_to(ROOT)};{CLADE_README.relative_to(ROOT)};{STAGE3B_ROLES.relative_to(ROOT)}", "taxa": "focal clade definitions and source-tree taxon mapping", "genomic_coverage": "published branch definitions; no standalone reference tree", "chr4_in_source": "not separable", "coordinates_available": False, "inference_method": "quadripartition labels, not an extractable independent tree", "predates_stage4c": True, "independent_of_chr4": "cannot establish independence or extract a reference topology", "available": False, "hierarchy_rank": "E (not usable as a reference tree)"},
        {"reference_id": "resolved_genetrees_absent", "path": "resolved-genetrees.tre.gz mentioned in clade-analysis/README.md but absent", "taxa": "multi-taxon", "genomic_coverage": "would be genome-wide", "chr4_in_source": "unknown", "coordinates_available": "mapping exists separately", "inference_method": "source gene trees", "predates_stage4c": True, "independent_of_chr4": "potentially, but unavailable", "available": False, "hierarchy_rank": "A/C (unavailable)"},
    ]


def validate_frozen() -> None:
    manifest = json.loads(STAGE4B_MANIFEST.read_text())
    for relative, expected in manifest["inputs"].items():
        assert sha256(ROOT / relative) == expected
    for relative, expected in manifest["outputs"].items():
        # PROJECT_STATUS.md is a cumulative project log explicitly updated by
        # each later stage; Stage 4B scientific artifacts remain immutable.
        if relative == "PROJECT_STATUS.md":
            continue
        assert sha256(ROOT / relative) == expected
    assert manifest["treatments"]["T0"] == "all resolved windows, weight 1"
    assert manifest["treatments"]["T1"] == "each consecutive identical-topology run has total weight 1"
    assert manifest["treatments"]["T2"] == "frozen-mask windows excluded; background weight 1"
    assert manifest["treatments"]["T3"] == "frozen-mask windows weight 0.5; background weight 1"


def write_report(validation, references, controls, jackknife, candidates) -> None:
    by_clade = {clade: [r for r in validation if r["clade"] == clade] for clade in CLADES}
    lines = ["# Neoaves chromosome 4 Stage 4C: independent validation", "",
             "## Independent reference search", "", "The repository search identified the following candidates before treatment agreement was evaluated:", "",
             "| candidate | available | evidence | independence |", "|---|---|---|---|"]
    for row in candidates:
        lines.append(f"| {row['reference_id']} | {row['available']} | {row['hierarchy_rank']}: `{row['path']}` | {row['independent_of_chr4']} |")
    lines += ["", "No per-locus Newick gene trees, non-chr4 species tree, frozen ASTRAL tree, or fixed empirical species-tree estimator configuration is present. The README mentions `resolved-genetrees.tre.gz`, but that file is absent.", "",
              "## Reference selection rule", "", "The rule was fixed before comparison: choose the first available item in A→F order, preferring an analysis that can exclude chr4. No usable A–E item exists: the repository has no gene trees, species tree, or extractable published reference tree. The q1 quadripartition definitions label tested branches but do not provide a standalone independent topology. The primary reference is therefore the available F item, `nonchr4_genomewide_QQS_v1`, after excluding every metadata label whose chromosome group is 4. Selection did not use agreement with T0–T3.", "",
              "## Reference topology/support", "", "Support is the frequency of the dominant q1/q2/q3 state per non-chr4 locus, matching the Stage 4B estimand. Intervals are 95% chromosome-cluster bootstrap intervals.", "", "| clade | non-chr4 q1/q2/q3 | topology | margin | 95% margin CI | loci | chromosome units |", "|---|---:|---|---:|---:|---:|---:|"]
    for clade in CLADES:
        r = references[clade]
        lines.append(f"| {clade} | {r['q1']:.3f}/{r['q2']:.3f}/{r['q3']:.3f} | {r['dominant']} | {r['margin']:.3f} | [{r['margin_ci_low']:.3f}, {r['margin_ci_high']:.3f}] | {r['n_loci']} | {r['n_chromosome_units']} |")
    lines += ["", "## T0-T3 comparison", "", "| clade | treatment | dominant | matches reference | M_reference | Delta vs T0 | 95% M_reference CI |", "|---|---|---|---|---:|---:|---:|"]
    for clade in CLADES:
        for row in by_clade[clade]:
            lines.append(f"| {clade} | {row['treatment']} | {row['treatment_dominant']} | {row['matches_reference']} | {row['M_reference']:.3f} | {row['delta_M_reference_vs_T0']:.3f} | [{row['M_reference_ci_low']:.3f}, {row['M_reference_ci_high']:.3f}] |")
    lines += ["", "## N61 focal interpretation", ""]
    n61 = {r["treatment"]: r for r in by_clade["N61"]}
    lines += [f"The independent non-chr4 reference is {n61['T0']['reference_topology']}. T0 is q2 dominant and disagrees. T1 increases M_reference by {n61['T1']['delta_M_reference_vs_T0']:.3f} and becomes slightly q1 dominant. T2 and T3 increase M_reference by {n61['T2']['delta_M_reference_vs_T0']:.3f} and {n61['T3']['delta_M_reference_vs_T0']:.3f}, respectively, while remaining q2 dominant. All three changes are in the direction predicted by the earlier structural q2 enrichment result; the structural-only effects are modest.", "",
              "## Columbea", "", "The non-chr4 reference and T0–T3 are all q1 dominant. T1 substantially reduces the q1 margin, while T2/T3 have smaller negative effects. T0 already agrees with the independent reference.", "",
              "## N62", "", "The non-chr4 reference is q2, whereas T0 is q1 dominant. T1 moves strongly toward the reference and changes dominance to q2. T2 and T3 move modestly toward q2 but remain q1 dominant. Thus N62 also shows empirical weighting bias relative to non-chr4 evidence, but its improvement is specifically block-aware: Stage 4A found structural-region q2 depletion, and structural exclusion does not reproduce the topology change.", "",
              "## Genome-wide/full-tree analysis if available", "", "A full species-tree inference was not run. Per-locus multi-taxon Newick trees and a fixed empirical estimator/configuration are absent. This is topology-reference validation using independent non-chr4 focal quartet summaries.", "",
              "## Control chromosomes if available", ""]
    for clade in CLADES:
        subset = [r for r in controls if r["focal_clade"] == clade and r["chromosome"] != "chr4"]
        changed = sum(str(r["topology_changed"]).lower() == "true" if isinstance(r["topology_changed"], str) else bool(r["topology_changed"]) for r in subset)
        chr4 = next(r for r in controls if r["focal_clade"] == clade and r["chromosome"] == "chr4")
        effects = sorted(abs(float(r["delta_M_reference_block_vs_window"])) for r in subset)
        percentile_rank = sum(x <= abs(float(chr4["delta_M_reference_block_vs_window"])) for x in effects) / len(effects)
        lines.append(f"- {clade}: {changed}/{len(subset)} non-chr4 main chromosomes change dominant topology after run normalization. Chr4 |Delta M_reference|={abs(float(chr4['delta_M_reference_block_vs_window'])):.3f}, empirical percentile={percentile_rank:.1%} among non-chr4 controls.")
    lines += ["", "## Uncertainty", "", f"Reference intervals use {BOOTSTRAPS:,} chromosome-cluster bootstrap replicates. Dense loci are never resampled independently. Treatment intervals and paired deltas are the frozen Stage 4B complete-run bootstrap results. Chromosome jackknifing removes one non-chr4 chromosome unit at a time.", ""]
    for clade in CLADES:
        subset = [r for r in jackknife if r["clade"] == clade]
        dom = Counter(str(r["dominant_topology"]) for r in subset)
        margins = [float(r["M_reference"]) for r in subset]
        lines.append(f"- {clade}: jackknife topology counts {dict(dom)}; M_reference range {min(margins):.3f} to {max(margins):.3f}.")
    lines += ["", "## Interpretation category", ""]
    for clade in CLADES:
        lines.append(f"- {clade}: **{by_clade[clade][0]['interpretation_category']}**")
    lines += ["", "## What claim is now justified", "", "For N61, naive dense chr4 locus weighting produces a focal relationship that disagrees with stable non-chr4 genome-wide quartet evidence. Block normalization reverses the dominant relationship toward that reference; independently frozen structural exclusion/downweighting moves support in the same direction but does not reverse dominance. This supports focal empirical species-tree bias from dense chr4 weighting and strong block-aware, partial structural-aware correction relative to the independent reference.", "", "N62 independently shows naive chr4 disagreement and block-aware correction toward its q2 non-chr4 reference, while its structural-aware effects are modest and do not reverse dominance. Columbea already agrees with its q1 reference, and all weighting treatments reduce rather than improve its reference margin.", "",
              "## Remaining limitations", "", "The independent evidence consists of published per-locus focal quartet summaries, not recoverable multi-taxon gene trees or a re-estimated full species tree. Locus states on the same chromosome remain linked, which is why uncertainty uses chromosomes. The reference and chr4 tracks derive from the same original genome-wide gene-tree study, although the primary baseline excludes chr4 completely. The result validates focal relationships and weighting behavior; it does not establish a causal rearrangement mechanism or resolve the underpowered Stage 3B 2:2 test.", ""]
    (RESULTS / "stage4c_report.md").write_text("\n".join(lines))


def main() -> None:
    validate_frozen()
    rows = load_genomewide_rows()
    references = reference_results(rows)
    validation = primary_rows(references)
    controls = control_rows(rows, references)
    jackknife = jackknife_rows(rows, references)
    candidates = candidate_rows()
    paths = {
        "validation": RESULTS / "stage4c_independent_validation.tsv",
        "controls": RESULTS / "stage4c_control_chromosomes.tsv",
        "jackknife": RESULTS / "stage4c_chromosome_jackknife.tsv",
        "candidates": RESULTS / "stage4c_reference_candidates.tsv",
    }
    write_tsv(paths["validation"], validation); write_tsv(paths["controls"], controls)
    write_tsv(paths["jackknife"], jackknife); write_tsv(paths["candidates"], candidates)
    write_report(validation, references, controls, jackknife, candidates)
    subprocess.run(["Rscript", str(Path(__file__).with_name("04c_make_figures.R"))], cwd=ROOT, check=True)
    outputs = [*paths.values(), RESULTS / "stage4c_report.md",
               FIGURES / "stage4c_figure_A_reference_comparison.pdf", FIGURES / "stage4c_figure_A_reference_comparison.png",
               FIGURES / "stage4c_figure_B_N61_focal.pdf", FIGURES / "stage4c_figure_B_N61_focal.png",
               FIGURES / "stage4c_figure_C_control_context.pdf", FIGURES / "stage4c_figure_C_control_context.png",
               ROOT / "PROJECT_STATUS.md"]
    inputs = [META, REC, RAW_README, CLADE_README, DRAW_SCRIPT, STAGE4B, STAGE4B_MANIFEST, STAGE3B_ROLES]
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
    manifest = {"stage": "4C", "git_commit": commit, "date_utc": datetime.now(timezone.utc).isoformat(),
                "stage4b_manifest_checksum": sha256(STAGE4B_MANIFEST), "seed": SEED,
                "reference_bootstrap_replicates": BOOTSTRAPS,
                "reference_selection_rule": "first available A-to-F hierarchy item, prioritizing chr4 exclusion; selected non-chr4 genome-wide focal QQS summaries without inspecting T0-T3 agreement",
                "primary_reference_id": REFERENCE_ID, "chr4_excluded": True,
                "reference_estimand": "frequency of dominant q1/q2/q3 state per locus",
                "uncertainty_unit": "chromosome group; dense loci are not independently bootstrapped",
                "full_species_tree_inference": {"performed": False, "reason": "per-locus multi-taxon Newick trees and fixed empirical estimator configuration absent"},
                "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
                "outputs": {str(path.relative_to(ROOT)): sha256(path) for path in outputs},
                "prior_treatments_redefined": False, "stage3b_eligibility_relaxed": False}
    (RESULTS / "stage4c_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
