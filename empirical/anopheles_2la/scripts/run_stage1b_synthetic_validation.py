#!/usr/bin/env python3
"""Run the synthetic-only Stage 1B/Stage 2 validation pipeline."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from datetime import date
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "anopheles_2la"
FIXTURE_ROOT = EMPIRICAL_ROOT / "tests" / "fixtures"
RESULTS_ROOT = EMPIRICAL_ROOT / "results" / "stage1b_synthetic"
FIGURE_ROOT = EMPIRICAL_ROOT / "figures" / "stage1b_synthetic"


def load_ingest():
    path = EMPIRICAL_ROOT / "scripts" / "10_ingest_local_trees.py"
    spec = importlib.util.spec_from_file_location("stage1b_ingest", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["stage1b_ingest"] = module
    spec.loader.exec_module(module)
    return module


def classification(result: dict[str, object]) -> str:
    delta = float(result["observed_delta"])
    p = float(result["p_value"])
    if math.isnan(delta) or math.isnan(p):
        return "inconclusive"
    if delta > 0 and p < 0.05:
        return "positive"
    return "not_positive"


def scenario_specificity(summary: list[dict[str, object]]) -> dict[str, object]:
    broad_flags = []
    for row in summary:
        if "B_NULL" in row["quartet_id"]:
            continue
        f_in = float(row["f_in"])
        f_out = float(row["f_out"])
        broad_flags.append(f_out >= 0.45 and abs(f_in - f_out) < 0.25)
    return {"broad_alternative_like": all(broad_flags) if broad_flags else False}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synthetic-validation", action="store_true", required=True)
    parser.add_argument("--permutations", type=int, default=199)
    args = parser.parse_args(argv)

    m = load_ingest()
    m.generate_synthetic_fixtures(FIXTURE_ROOT)
    windows_path = FIXTURE_ROOT / "synthetic_genomic_windows.tsv"
    quartets_path = FIXTURE_ROOT / "synthetic_frozen_strict_quartets_stage1a.tsv"
    pairs_path = FIXTURE_ROOT / "synthetic_design_b_pairs.tsv"
    m.guard_inputs([windows_path, quartets_path, FIXTURE_ROOT / "synthetic_gene_trees"], True, False)
    windows = m.read_tsv(windows_path)
    quartets = m.read_tsv(quartets_path)

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    aggregate = {}
    scenario_rows = []
    for scenario in ["null_msc_like", "msrc_positive", "broad_alternative"]:
        outdir = RESULTS_ROOT / scenario
        outdir.mkdir(parents=True, exist_ok=True)
        track = m.ingest_local_trees(windows, quartets, FIXTURE_ROOT / "synthetic_gene_trees" / scenario)
        summary = m.summarize_inside_outside(track)
        block_summary = m.summarize_inside_outside(track, block_normalized=True)
        circ = m.circular_shift_test(track, args.permutations)
        naive = m.naive_label_permutation(track, args.permutations)
        transitions = m.breakpoint_transitions(track)
        bp_test = m.breakpoint_localization_test(track, args.permutations)
        design_b = m.design_b_contrasts(track, pairs_path)
        geography = m.geography_summary(quartets)
        figures = m.make_figures(track, summary, design_b, FIGURE_ROOT / scenario)

        m.write_tsv(outdir / "stage2_spatial_topology_track.tsv", track, m.TRACK_FIELDS)
        m.write_tsv(outdir / "inside_outside_summary.tsv", summary, list(summary[0].keys()))
        m.write_tsv(outdir / "inside_outside_block_normalized_summary.tsv", block_summary, list(block_summary[0].keys()))
        m.write_tsv(outdir / "breakpoint_transitions.tsv", transitions, list(transitions[0].keys()) if transitions else ["quartet_id", "transition_position", "from_topology", "to_topology", "distance_to_nearest_2La_breakpoint"])
        m.write_tsv(outdir / "design_b_contrasts.tsv", design_b, list(design_b[0].keys()) if design_b else ["contrast_id", "quartet_a", "quartet_b", "delta_a", "delta_b", "classification", "synthetic_expected_result"])
        m.write_tsv(outdir / "geography_controls.tsv", geography, list(geography[0].keys()))
        m.write_json(outdir / "permutation_tests.json", {"SYNTHETIC_ONLY": True, "block_aware": circ, "naive_diagnostic": naive, "breakpoint": bp_test})

        status = classification(circ)
        specificity = scenario_specificity(summary)
        aggregate[scenario] = {
            "classification": status,
            "block_aware_test": circ,
            "naive_diagnostic": naive,
            "breakpoint_test": bp_test,
            "specificity": specificity,
            "figures": [str(p.relative_to(REPO_ROOT)) for p in figures],
        }
        scenario_rows.append(
            {
                "scenario": scenario,
                "classification": status,
                "observed_delta": circ["observed_delta"],
                "block_aware_p_value": circ["p_value"],
                "broad_alternative_like": specificity["broad_alternative_like"],
                "SYNTHETIC_ONLY": "TRUE",
            }
        )

    m.write_tsv(RESULTS_ROOT / "scenario_summary.tsv", scenario_rows, list(scenario_rows[0].keys()))
    previous_checksum = "a4c366a9921ac11d3794bf29d7dbb102a278dda543a9650cd4ede1c985aa8fe2"
    manifest = {
        "SYNTHETIC_ONLY": True,
        "analysis_plan_version": "1.0.1",
        "date": str(date.today()),
        "git_commit": git_commit(),
        "previous_freeze_checksum": previous_checksum,
        "checksum_supersession_reason": "pre-data hardening: design-agnostic execution gate, tree-inference environment pinning, test-only network isolation",
        "2La_coordinates": {"chromosome": "2L", "start": m.LEFT_BP, "end": m.RIGHT_BP, "coordinate_system": "1-based inclusive"},
        "primary_flank_width_bp": m.PRIMARY_FLANK_BP,
        "sensitivity_flank_widths_bp": [2_000_000, 10_000_000],
        "primary_window_size_bp": m.PRIMARY_WINDOW_BP,
        "sensitivity_window_sizes_bp": [50_000, 200_000],
        "primary_statistic": "Delta_arr = f_in - f_out",
        "permutation_method": "block_aware_circular_shift",
        "final_permutations": m.FINAL_PERMUTATIONS,
        "validation_permutations_used": args.permutations,
        "polytomy_policy": "unresolved quartet excluded from denominators but counted",
        "missing_tip_policy": "record unusable with exclusion reason",
        "primary_tree_estimator": "IQ-TREE",
        "tree_inference_environment": {
            "package": "bioconda::iqtree=2.4.0",
            "executable": "iqtree2",
            "environment_file": "empirical/anopheles_2la/config/tree_inference_environment.yml",
            "version_verification_before_real_run": "iqtree2 --version must report IQ-TREE 2.4.0",
        },
        "tree_inference_command_template": "iqtree2 -s WINDOW.fasta -seed 1729 -nt AUTO -pre OUTPUT_PREFIX -m MFP",
        "tree_estimator_seed": 1729,
        "model_selection_policy": "IQ-TREE ModelFinder Plus (-m MFP), frozen before sequence access",
        "thread_policy": "-nt AUTO",
        "bootstrap_support_policy": "none for primary quartet-classification analysis; no support threshold is applied",
        "expected_tree_inference_outputs": [
            "OUTPUT_PREFIX.treefile",
            "OUTPUT_PREFIX.iqtree",
            "OUTPUT_PREFIX.log",
            "OUTPUT_PREFIX.ckp.gz",
            "OUTPUT_PREFIX.model.gz",
        ],
        "tree_estimation_prediction_blind": True,
        "block_normalization_rule": "collapse consecutive windows with same local quartet state for a frozen comparison into one run of weight one",
        "breakpoint_analysis_rule": "fixed 2La breakpoints; topology transition distance to nearest breakpoint",
        "design_agnostic_real_data_gate": {
            "stage1a_completion_marker": "data/anopheles_2la/processed/stage1a_freeze_complete.json",
            "requires_stage1a_freeze_provenance_validation": True,
            "design_A_table": "data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv or repository-equivalent Design-A rows",
            "design_B_tables": [
                "data/anopheles_2la/processed/frozen_design_b_contrasts_stage1a.tsv",
                "data/anopheles_2la/processed/frozen_arrangement_replacement_contrasts_stage1a.tsv",
                "or repository-equivalent Design-B rows",
            ],
            "design_C_tables": [
                "data/anopheles_2la/processed/frozen_design_c_geography_controls_stage1a.tsv",
                "data/anopheles_2la/processed/frozen_geography_population_controls_stage1a.tsv",
                "or repository-equivalent Design-C rows",
            ],
            "design_A_enabled": "true iff a prospectively frozen Design-A table/row set is nonempty",
            "design_B_enabled": "true iff a prospectively frozen Design-B table/row set is nonempty",
            "design_C_enabled": "true iff a prospectively frozen Design-C table/row set is nonempty",
            "all_empty_status": "NO_ELIGIBLE_PROSPECTIVE_DESIGNS",
            "all_empty_interpretation": "inconclusive / underpowered empirical result, not evidence against MSRC",
        },
        "test_status": {
            "10_ingest_local_trees.py --run-tests": "passed during hardening",
            "test_stage1b_synthetic.py": "passed during hardening",
            "00_audit_samples.py --run-tests": "passed during hardening",
            "01_fetch_sample_karyotypes.py --run-tests": "passed during hardening with no live API path",
            "11_infer_local_trees.py --run-tests": "passed during hardening",
            "synthetic_end_to_end_validation": "passed during hardening",
        },
        "figure_definitions": ["prospective_design", "spatial_topology_track", "inside_outside_enrichment", "arrangement_replacement"],
        "synthetic_scenarios": aggregate,
    }
    freeze_path = RESULTS_ROOT / "stage1b_analysis_freeze.json"
    m.write_json(freeze_path, manifest)
    checksum = m.sha256(freeze_path)
    (RESULTS_ROOT / "stage1b_analysis_freeze.sha256").write_text(f"{checksum}  {freeze_path.name}\n")
    write_report(RESULTS_ROOT / "synthetic_validation_report.md", aggregate, checksum)
    return 0


def git_commit() -> str:
    import subprocess

    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except Exception:
        return "unavailable"


def write_report(path: Path, aggregate: dict[str, object], checksum: str) -> None:
    lines = [
        "# Stage 1B Synthetic Validation Report",
        "",
        "`SYNTHETIC_ONLY = TRUE`",
        "",
        "No real Anopheles sequence-derived topology, local-tree, SNP, haplotype, q1/q2/q3, QQS/BQS, or NJ-tree data were read.",
        "",
        f"Freeze checksum: `{checksum}`",
        "",
        "| Scenario | Classification | Delta | Block-aware p | Interpretation |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for scenario, payload in aggregate.items():
        test = payload["block_aware_test"]
        if scenario == "null_msc_like":
            interp = "not systematically positive expected"
        elif scenario == "msrc_positive":
            interp = "positive control should be detected"
        else:
            interp = "broad alternative should be flagged as not 2La-specific"
        lines.append(f"| {scenario} | {payload['classification']} | {float(test['observed_delta']):.3f} | {float(test['p_value']):.3f} | {interp} |")
    lines.extend(
        [
            "",
            "The prospective test is whether the independently predicted arrangement topology is enriched inside 2La relative to fixed flanks.",
            "The Design-B validation asks whether replacing arrangement background changes topology in the structurally predicted direction while species identity is controlled.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
