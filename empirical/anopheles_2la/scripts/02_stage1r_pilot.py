#!/usr/bin/env python3
"""Stage 1R MalariaGEN regional SNP-access and homozygote pilot for Anopheles 2La.

This script intentionally performs only a small predefined regional pilot. It does
not run a full 2L scan, does not process raw reads, and does not compute final
manuscript-wide statistics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import itertools
import json
import math
import os
import platform
import sys
import time
import unittest
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "data" / "anopheles_2la"
META = DATA / "metadata"
PROCESSED = DATA / "processed"
RESULTS = REPO / "empirical" / "anopheles_2la" / "results"
FIGURES = REPO / "empirical" / "anopheles_2la" / "figures"
README = REPO / "empirical" / "anopheles_2la" / "README.md"
PROJECT_STATUS = REPO / "PROJECT_STATUS.md"

SAMPLE_MANIFEST = META / "sample_manifest.tsv"
RAW_KARYOTYPES = META / "stage1a_raw_2la_karyotypes.tsv"
STAGE1A_RECON = RESULTS / "stage1a_sample_reconciliation.tsv"
STAGE1A_REPORT = RESULTS / "stage1a_metadata_karyotype_freeze.md"
STAGE1A_PROV = META / "stage1a_malariagen_api_provenance.json"

PRIMARY_SAMPLES = PROCESSED / "stage1r_primary_homozygote_samples.tsv"
PILOT_WINDOWS = PROCESSED / "stage1r_pilot_windows.tsv"
PAIRWISE = PROCESSED / "stage1r_pairwise_distances.tsv"
TREES = PROCESSED / "stage1r_pilot_trees.nwk"
ENV_PROV = RESULTS / "stage1r_environment_provenance.json"
WINDOW_QC = RESULTS / "stage1r_window_qc.tsv"
FILTER_PLAN = RESULTS / "stage1r_filtering_plan.md"
CROSSED = RESULTS / "stage1r_crossed_distance_test.tsv"
CLASS_SUMMARY = RESULTS / "stage1r_distance_class_summary.tsv"
QUARTETS = RESULTS / "stage1r_pilot_quartet_support.tsv"
HAP_MEMBERSHIP = RESULTS / "stage1r_haplotype_membership.tsv"
REPORT = RESULTS / "stage1r_report.md"
MANIFEST = RESULTS / "stage1r_manifest.json"
STAGE2_PLAN = RESULTS / "stage2r_preanalysis_plan.md"
FIG_PDF = FIGURES / "stage1r_pilot_signal.pdf"
FIG_PNG = FIGURES / "stage1r_pilot_signal.png"

CHROM = "2L"
INV_START = 20_524_058
INV_END = 42_165_532
WIN = 50_000
RELEASE = "3.10"
SAMPLE_SET = "fontaine-2015-rebuild"
BUCKET = "gs://vo_agam_release_master_us_central1"
EXPECTED_COUNTS = {
    "arabiensis": {"n": 12, "2L+a/2L+a": 0, "heterokaryotype": 0, "2La/2La": 12},
    "coluzzii": {"n": 11, "2L+a/2L+a": 8, "heterokaryotype": 0, "2La/2La": 3},
    "gambiae": {"n": 26, "2L+a/2L+a": 8, "heterokaryotype": 6, "2La/2La": 12},
    "melas": {"n": 4, "2L+a/2L+a": 0, "heterokaryotype": 4, "2La/2La": 0},
    "merus": {"n": 9, "2L+a/2L+a": 0, "heterokaryotype": 9, "2La/2La": 0},
    "quadriannulatus": {"n": 10, "2L+a/2L+a": 10, "heterokaryotype": 0, "2La/2La": 0},
}
ARR = {"0": "2L+a/2L+a", "1": "heterokaryotype", "2": "2La/2La"}
STATE_FROM_ARR = {"2L+a/2L+a": "A0_homozygous", "2La/2La": "A1_homozygous", "heterokaryotype": "heterokaryotype"}
PRIMARY_ARR = {"2L+a/2L+a", "2La/2La"}
MAX_SITE_MISSING = 0.25
MAX_SAMPLE_MISSING = 0.50
MIN_MINOR_AC = 2
NO_FULL_SCAN_MAX_TOTAL_BP = 10 * WIN


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: format_value(row.get(k, "")) for k in fields})


def format_value(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        return f"{v:.12g}"
    if isinstance(v, (np.floating,)):
        if np.isnan(v):
            return "nan"
        return f"{float(v):.12g}"
    if isinstance(v, (np.integer,)):
        return str(int(v))
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fixed_windows() -> list[dict[str, Any]]:
    length = INV_END - INV_START + 1
    def inside_centered(window_id: str, fraction: float) -> dict[str, Any]:
        center = INV_START + round(length * fraction)
        start = center - WIN // 2
        end = start + WIN - 1
        return {
            "window_id": window_id,
            "chrom": CHROM,
            "start": start,
            "end": end,
            "position_class": "inside_deep",
            "distance_to_nearest_breakpoint": min(start - INV_START, INV_END - end),
            "selection_rule": f"50-kb window centered at {fraction:.2f} of frozen 2La interval; selected before SNP/topology inspection",
        }
    rows = [
        {
            "window_id": "outside_left_distal_10Mb",
            "chrom": CHROM,
            "start": 10_000_000,
            "end": 10_049_999,
            "position_class": "outside_left",
            "distance_to_nearest_breakpoint": INV_START - 10_049_999,
            "selection_rule": "fixed distal left 2L control at 10.00 Mb; selected before SNP/topology inspection",
        },
        {
            "window_id": "outside_left_flank_50kb",
            "chrom": CHROM,
            "start": INV_START - WIN,
            "end": INV_START - 1,
            "position_class": "outside_left",
            "distance_to_nearest_breakpoint": 1,
            "selection_rule": "50-kb window immediately left of frozen 2La interval without boundary overlap",
        },
        {
            "window_id": "inside_left_boundary_50kb",
            "chrom": CHROM,
            "start": INV_START + 1,
            "end": INV_START + WIN,
            "position_class": "inside_near_boundary",
            "distance_to_nearest_breakpoint": 1,
            "selection_rule": "50-kb window immediately inside left 2La boundary without crossing breakpoint",
        },
        inside_centered("inside_deep_q1", 0.25),
        inside_centered("inside_deep_center", 0.50),
        inside_centered("inside_deep_q3", 0.75),
        {
            "window_id": "inside_right_boundary_50kb",
            "chrom": CHROM,
            "start": INV_END - WIN,
            "end": INV_END - 1,
            "position_class": "inside_near_boundary",
            "distance_to_nearest_breakpoint": 1,
            "selection_rule": "50-kb window immediately inside right 2La boundary without crossing breakpoint",
        },
        {
            "window_id": "outside_right_flank_50kb",
            "chrom": CHROM,
            "start": INV_END + 1,
            "end": INV_END + WIN,
            "position_class": "outside_right",
            "distance_to_nearest_breakpoint": 1,
            "selection_rule": "50-kb window immediately right of frozen 2La interval without boundary overlap",
        },
    ]
    rows.sort(key=lambda r: r["start"])
    return rows


def region_string(w: dict[str, Any]) -> str:
    return f"{w['chrom']}:{w['start']}-{w['end']}"


def load_primary_rows() -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    rows = read_tsv(SAMPLE_MANIFEST)
    raw = read_tsv(RAW_KARYOTYPES)
    raw_by_id = {r["sample_id"]: r for r in raw}
    out = []
    for r in rows:
        num = r["2La_karyotype_raw"]
        arr = ARR.get(num, "unknown")
        include = arr in PRIMARY_ARR
        out.append({
            "sample_id": r["sample_id"],
            "species": r["species"],
            "karyotype_numeric": num,
            "arrangement": arr,
            "sample_set": r["sample_set"],
            "study_id": SAMPLE_SET,
            "country": r["country"],
            "population": r["population"],
            "include_primary": include,
            "exclusion_reason": "" if include else "heterokaryotype excluded from primary unphased diploid Stage-1R analysis",
            "karyotype_2La_mean": raw_by_id.get(r["sample_id"], {}).get("karyotype_2La_mean", ""),
        })
    out.sort(key=lambda r: (not r["include_primary"], r["species"], r["arrangement"], r["sample_id"]))
    return out, rows


def validate_primary(primary: list[dict[str, Any]], source: list[dict[str, str]]) -> None:
    assert len(source) == 72, f"expected 72 Fontaine samples, got {len(source)}"
    ids = [r["sample_id"] for r in primary]
    assert len(ids) == len(set(ids)), "duplicated sample IDs in primary table"
    included = [r for r in primary if r["include_primary"]]
    excluded = [r for r in primary if not r["include_primary"]]
    assert len(included) == 53, f"expected 53 primary homozygotes, got {len(included)}"
    assert len(excluded) == 19, f"expected 19 excluded heterokaryotypes, got {len(excluded)}"
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for r in primary:
        counts[r["species"]][r["arrangement"]] += 1
        counts[r["species"]]["n"] += 1
    for sp, exp in EXPECTED_COUNTS.items():
        got = counts[sp]
        for k, v in exp.items():
            assert got[k] == v, f"{sp} {k}: expected {v}, got {got[k]}"


def validate_windows(windows: list[dict[str, Any]]) -> None:
    assert INV_START == 20_524_058 and INV_END == 42_165_532, "frozen 2La coordinates changed"
    assert len(windows) == 8, "Stage 1R pilot must use exactly 8 predefined windows"
    total_bp = sum(int(w["end"]) - int(w["start"]) + 1 for w in windows)
    assert total_bp <= NO_FULL_SCAN_MAX_TOTAL_BP, "pilot exceeds no-full-scan bp limit"
    for w in windows:
        assert int(w["end"]) - int(w["start"]) + 1 == WIN
        inside = int(w["start"]) >= INV_START and int(w["end"]) <= INV_END
        outside = int(w["end"]) < INV_START or int(w["start"]) > INV_END
        assert inside or outside, f"boundary-overlap window classified as clean: {w['window_id']}"
        if w["position_class"].startswith("inside"):
            assert inside, f"inside class outside interval: {w['window_id']}"
        else:
            assert outside, f"outside class overlaps interval: {w['window_id']}"


def env_provenance(ag3: Any, migration_note: str) -> dict[str, Any]:
    import malariagen_data
    try:
        import numpy
        npv = numpy.__version__
    except Exception:
        npv = "not-importable"
    try:
        import pandas
        pdv = pandas.__version__
    except Exception:
        pdv = "not-importable"
    try:
        import allel
        allev = allel.__version__
    except Exception:
        allev = "not-used/not-importable"
    try:
        import Bio
        biov = Bio.__version__
    except Exception:
        biov = "not-used/not-importable"
    return {
        "python_version": sys.version.replace("\n", " "),
        "python_executable": sys.executable,
        "python_3_11_migration": migration_note,
        "malariagen_data_version": malariagen_data.__version__,
        "numpy_version": npv,
        "pandas_version": pdv,
        "scikit_allel_version_if_used": allev,
        "biopython_version_if_used": biov,
        "gcs_auth_method": "Google Application Default Credentials; credentials/tokens not stored",
        "Ag3_release": RELEASE,
        "sample_set": SAMPLE_SET,
        "bucket": BUCKET,
        "api_methods_used_stage1r": ["sample_metadata", "snp_calls", "biallelic_snp_calls", "haplotypes"],
        "api_method_signatures": {name: str(inspect.signature(getattr(ag3, name))) for name in ["sample_metadata", "snp_calls", "biallelic_snp_calls", "haplotypes"] if hasattr(ag3, name)},
    }


def dosage_and_filters(gt: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # gt: variants x samples x ploidy, allele codes 0/1 for biallelic, -1 missing.
    called = np.all(gt >= 0, axis=2)
    dosage = np.where(called, gt.sum(axis=2), -1).astype(np.int16)
    n_called = called.sum(axis=1)
    alt_ac = np.where(called, dosage, 0).sum(axis=1)
    an = 2 * n_called
    minor_ac = np.minimum(alt_ac, an - alt_ac)
    segregating = (minor_ac > 0) & (n_called > 0)
    keep = segregating & ((1 - n_called / gt.shape[1]) <= MAX_SITE_MISSING) & (minor_ac >= MIN_MINOR_AC)
    return dosage, called, alt_ac, minor_ac, keep


def distance_matrix(dosage: np.ndarray, called: np.ndarray, keep: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    idx = np.where(keep)[0]
    n = dosage.shape[1]
    dist = np.full((n, n), np.nan, dtype=float)
    compared = np.zeros((n, n), dtype=int)
    for i in range(n):
        dist[i, i] = 0.0
        compared[i, i] = int(len(idx))
        for j in range(i + 1, n):
            both = called[idx, i] & called[idx, j]
            m = int(both.sum())
            compared[i, j] = compared[j, i] = m
            if m:
                d = np.abs(dosage[idx[both], i] - dosage[idx[both], j]).mean() / 2.0
                dist[i, j] = dist[j, i] = float(d)
    return dist, compared


def pairwise_rows(window_id: str, samples: list[dict[str, Any]], dist: np.ndarray, compared: np.ndarray) -> list[dict[str, Any]]:
    rows = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            rows.append({
                "window_id": window_id,
                "sample1": samples[i]["sample_id"],
                "sample2": samples[j]["sample_id"],
                "species1": samples[i]["species"],
                "species2": samples[j]["species"],
                "arrangement1": samples[i]["arrangement"],
                "arrangement2": samples[j]["arrangement"],
                "distance": dist[i, j],
                "n_sites_compared": compared[i, j],
            })
    return rows


def class_summaries(window: dict[str, Any], pairs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    buckets: dict[str, list[float]] = defaultdict(list)
    sites: dict[str, list[int]] = defaultdict(list)
    for p in pairs:
        if p["species1"] == p["species2"] and p["arrangement1"] == p["arrangement2"]:
            cls = "same_species_same_arrangement"
        elif p["species1"] == p["species2"]:
            cls = "same_species_opposite_arrangement"
        elif p["arrangement1"] == p["arrangement2"]:
            cls = "different_species_same_arrangement"
        else:
            cls = "different_species_opposite_arrangement"
        d = float(p["distance"])
        if not math.isnan(d):
            buckets[cls].append(d)
            sites[cls].append(int(p["n_sites_compared"]))
    out = []
    for cls in ["same_species_same_arrangement", "same_species_opposite_arrangement", "different_species_same_arrangement", "different_species_opposite_arrangement"]:
        vals = buckets.get(cls, [])
        out.append({
            "window_id": window["window_id"],
            "position_class": window["position_class"],
            "distance_class": cls,
            "n_pairs": len(vals),
            "mean_distance": float(np.mean(vals)) if vals else np.nan,
            "median_distance": float(np.median(vals)) if vals else np.nan,
            "min_distance": float(np.min(vals)) if vals else np.nan,
            "max_distance": float(np.max(vals)) if vals else np.nan,
            "mean_sites_compared": float(np.mean(sites[cls])) if vals else np.nan,
        })
    return out


def crossed_contrast(window: dict[str, Any], pairs: list[dict[str, Any]]) -> dict[str, Any]:
    within = []
    cross = []
    for p in pairs:
        sp = {p["species1"], p["species2"]}
        arr = {p["arrangement1"], p["arrangement2"]}
        if not sp <= {"gambiae", "coluzzii"}:
            continue
        d = float(p["distance"])
        if math.isnan(d):
            continue
        if p["species1"] == p["species2"] and p["arrangement1"] != p["arrangement2"]:
            within.append(d)
        if p["species1"] != p["species2"] and p["arrangement1"] == p["arrangement2"]:
            cross.append(d)
    mw = float(np.mean(within)) if within else np.nan
    mc = float(np.mean(cross)) if cross else np.nan
    return {
        "window_id": window["window_id"],
        "position_class": window["position_class"],
        "n_within_species_opposite_pairs": len(within),
        "mean_within_species_opposite_distance": mw,
        "n_cross_species_same_pairs": len(cross),
        "mean_cross_species_same_distance": mc,
        "C": mw - mc if within and cross else np.nan,
    }


def quartet_support(window: dict[str, Any], samples: list[dict[str, Any]], dist: np.ndarray) -> dict[str, Any]:
    groups = {
        "G_s": [i for i, s in enumerate(samples) if s["species"] == "gambiae" and s["arrangement"] == "2L+a/2L+a"],
        "G_i": [i for i, s in enumerate(samples) if s["species"] == "gambiae" and s["arrangement"] == "2La/2La"],
        "C_s": [i for i, s in enumerate(samples) if s["species"] == "coluzzii" and s["arrangement"] == "2L+a/2L+a"],
        "C_i": [i for i, s in enumerate(samples) if s["species"] == "coluzzii" and s["arrangement"] == "2La/2La"],
    }
    counts = Counter()
    n = 0
    for gs, gi, cs, ci in itertools.product(groups["G_s"], groups["G_i"], groups["C_s"], groups["C_i"]):
        s_species = dist[gs, gi] + dist[cs, ci]
        s_arr = dist[gs, cs] + dist[gi, ci]
        s_third = dist[gs, ci] + dist[gi, cs]
        vals = {"species": s_species, "arrangement": s_arr, "third": s_third}
        if any(math.isnan(v) for v in vals.values()):
            continue
        n += 1
        best = min(vals.values())
        winners = [k for k, v in vals.items() if abs(v - best) < 1e-12]
        if len(winners) == 1:
            counts[winners[0]] += 1
        else:
            counts["tie"] += 1
    q_species = counts["species"] / n if n else np.nan
    q_arr = counts["arrangement"] / n if n else np.nan
    q_third = (counts["third"] + counts["tie"]) / n if n else np.nan
    return {
        "window_id": window["window_id"],
        "position_class": window["position_class"],
        "n_informative_quartets": n,
        "q_species": q_species,
        "q_arrangement": q_arr,
        "q_third": q_third,
        "D_arrangement_minus_species": q_arr - q_species if n else np.nan,
    }


def nj_newick(labels: list[str], D: np.ndarray) -> str:
    # Simple deterministic neighbor joining. Replaces NaN with max finite distance.
    finite = D[np.isfinite(D)]
    Dm = D.copy().astype(float)
    fill = float(finite.max()) if finite.size else 1.0
    Dm[~np.isfinite(Dm)] = fill
    np.fill_diagonal(Dm, 0.0)
    clusters = {i: labels[i] for i in range(len(labels))}
    active = list(range(len(labels)))
    next_id = len(labels)
    lengths = {}
    while len(active) > 2:
        n = len(active)
        total = {i: sum(Dm[i, j] for j in active if j != i) for i in active}
        best = None
        for a_i, i in enumerate(active):
            for j in active[a_i + 1:]:
                q = (n - 2) * Dm[i, j] - total[i] - total[j]
                key = (q, min(i, j), max(i, j))
                if best is None or key < best[0]:
                    best = (key, i, j)
        _, i, j = best
        li = 0.5 * Dm[i, j] + (total[i] - total[j]) / (2 * (n - 2))
        lj = Dm[i, j] - li
        li = max(0.0, float(li)); lj = max(0.0, float(lj))
        new_label = f"({clusters[i]}:{li:.6g},{clusters[j]}:{lj:.6g})"
        new_id = next_id; next_id += 1
        if new_id >= Dm.shape[0]:
            old = Dm
            Dm = np.zeros((old.shape[0] + 1, old.shape[1] + 1), dtype=float)
            Dm[:old.shape[0], :old.shape[1]] = old
        for k in active:
            if k not in (i, j):
                Dm[new_id, k] = Dm[k, new_id] = 0.5 * (Dm[i, k] + Dm[j, k] - Dm[i, j])
        clusters[new_id] = new_label
        active = [k for k in active if k not in (i, j)] + [new_id]
    i, j = active
    bl = max(0.0, Dm[i, j] / 2)
    return f"({clusters[i]}:{bl:.6g},{clusters[j]}:{bl:.6g});"


def sanitize_label(sample: dict[str, Any]) -> str:
    sp = sample["species"][:4]
    ar = "std" if sample["arrangement"] == "2L+a/2L+a" else "inv"
    return f"{sample['sample_id']}|{sp}|{ar}".replace(" ", "_")


def query_window(ag3: Any, window: dict[str, Any], sample_indices: list[int], primary: list[dict[str, Any]]) -> dict[str, Any]:
    t0 = time.perf_counter()
    reg = region_string(window)
    raw = ag3.snp_calls(region=reg, sample_sets=SAMPLE_SET, sample_indices=sample_indices, site_mask=None, inline_array=True, chunks="native")
    bial = ag3.biallelic_snp_calls(region=reg, sample_sets=SAMPLE_SET, sample_indices=sample_indices, site_mask=None, inline_array=True, chunks="native")
    raw_samples = [str(x) for x in raw["sample_id"].values.tolist()]
    bial_samples = [str(x) for x in bial["sample_id"].values.tolist()]
    expected_ids = [s["sample_id"] for s in primary]
    if raw_samples != expected_ids or bial_samples != expected_ids:
        raise RuntimeError(f"sample IDs returned out of order or incomplete for {window['window_id']}")
    gt = bial["call_genotype"].values
    dosage, called, alt_ac, minor_ac, keep = dosage_and_filters(gt)
    dist, compared = distance_matrix(dosage, called, keep)
    n_sites_raw = int(raw.sizes.get("variants", 0))
    n_bial = int(bial.sizes.get("variants", 0))
    n_called = called.sum(axis=1) if called.size else np.array([])
    sample_missing = 1 - called.mean(axis=0) if called.size else np.ones(len(primary))
    het = ((dosage == 1) & called).sum() / called.sum() if called.sum() else np.nan
    filt_sample_missing = 1 - called[keep].mean(axis=0) if keep.any() else np.ones(len(primary))
    qc = {
        "window_id": window["window_id"],
        "position_class": window["position_class"],
        "n_samples_requested": len(primary),
        "n_samples_returned": len(bial_samples),
        "n_sites_raw": n_sites_raw,
        "n_biallelic_snps": n_bial,
        "n_segregating_biallelic_snps": int(((minor_ac > 0) & (n_called > 0)).sum()) if n_bial else 0,
        "n_post_filter_snps": int(keep.sum()),
        "mean_missingness": float(1 - called.mean()) if called.size else 1.0,
        "max_sample_missingness": float(sample_missing.max()) if sample_missing.size else 1.0,
        "mean_heterozygosity": float(het) if not math.isnan(het) else np.nan,
        "callable": bool(keep.sum() >= 25 and (float(sample_missing.max()) if sample_missing.size else 1.0) <= MAX_SAMPLE_MISSING),
        "status": "ok" if keep.sum() >= 25 and (float(sample_missing.max()) if sample_missing.size else 1.0) <= MAX_SAMPLE_MISSING else "low_signal_or_missing",
        "notes": "direct MalariaGEN regional biallelic_snp_calls; no raw reads; filter applied after fixed-window query",
        "runtime_seconds": time.perf_counter() - t0,
        "approx_dataset_bytes": int(getattr(raw, "nbytes", 0)) + int(getattr(bial, "nbytes", 0)),
        "mean_post_filter_sample_missingness": float(filt_sample_missing.mean()) if filt_sample_missing.size else 1.0,
        "alt_ac_mean": float(np.mean(alt_ac[keep])) if keep.any() else np.nan,
        "minor_ac_mean": float(np.mean(minor_ac[keep])) if keep.any() else np.nan,
    }
    return {"qc": qc, "dist": dist, "compared": compared, "pairs": pairwise_rows(window["window_id"], primary, dist, compared), "newick": nj_newick([sanitize_label(s) for s in primary], dist)}


def haplotype_membership(ag3: Any, all_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    region = f"{CHROM}:{(INV_START + INV_END)//2}-{(INV_START + INV_END)//2 + 999}"
    ids = [r["sample_id"] for r in all_rows]
    by_id = {r["sample_id"]: r for r in all_rows}
    out = {sid: {"sample_id": sid, "species": by_id[sid]["species"], "arrangement": ARR.get(by_id[sid]["2La_karyotype_raw"], "unknown"), "haplotype_panel": "Ag3 default", "phased_available": False, "notes": "not returned by regional haplotypes query"} for sid in ids}
    try:
        ds = ag3.haplotypes(region=region, sample_sets=SAMPLE_SET, inline_array=True, chunks="native")
        returned = [str(x) for x in ds["sample_id"].values.tolist()]
        for sid in returned:
            if sid in out:
                out[sid]["phased_available"] = True
                out[sid]["notes"] = f"returned by regional haplotypes query {region}; both homologues represented in call_genotype_haplotypes"
        panel = f"Ag3 default haplotypes; region {region}; samples_returned={len(returned)}"
        for row in out.values():
            row["haplotype_panel"] = panel
    except Exception as e:
        msg = f"regional haplotypes query failed: {type(e).__name__}: {e}"
        for row in out.values():
            row["haplotype_panel"] = "Ag3 default haplotypes"
            row["notes"] = msg
    return [out[sid] for sid in ids]


def write_filter_plan(qc_rows: list[dict[str, Any]]) -> None:
    text = f"""# Stage 1R pilot SNP filtering plan

This filter was frozen for the Stage 1R pilot before interpreting distance or tree topology results.

Input data are direct MalariaGEN Ag3 regional calls from `biallelic_snp_calls`, restricted to the 53 frozen homozygous Fontaine-associated samples. Raw-read processing is not used.

Filter:

- use biallelic SNPs only;
- exclude invariant sites by requiring minor allele count > 0;
- require minor allele count >= {MIN_MINOR_AC};
- require site missingness <= {MAX_SITE_MISSING:.2f};
- record sample missingness and flag windows if any sample has missingness > {MAX_SAMPLE_MISSING:.2f};
- do not LD-prune for this initial within-window genealogy/distance pilot;
- do not choose, remove, or reorder windows based on topology.

Pairwise distance formula:

For each retained biallelic SNP, convert an unphased diploid genotype to alternate-allele dosage 0, 1, or 2. For samples `i` and `j`, compare only sites callable in both samples and compute

```text
d(i,j) = mean_s |dosage_i(s) - dosage_j(s)| / 2
```

The denominator 2 scales the per-site distance to [0,1].

Pre-filter and post-filter counts by window are recorded in `stage1r_window_qc.tsv`.
"""
    FILTER_PLAN.write_text(text)


def summarize_signal(cross_rows: list[dict[str, Any]], quartet_rows: list[dict[str, Any]]) -> dict[str, Any]:
    def inside(cls: str) -> bool:
        return cls.startswith("inside")
    def mean_vals(rows, key, pred):
        vals = [float(r[key]) for r in rows if pred(r["position_class"]) and not math.isnan(float(r[key]))]
        return float(np.mean(vals)) if vals else np.nan
    return {
        "mean_C_inside": mean_vals(cross_rows, "C", inside),
        "mean_C_outside": mean_vals(cross_rows, "C", lambda c: not inside(c)),
        "mean_D_inside": mean_vals(quartet_rows, "D_arrangement_minus_species", inside),
        "mean_D_outside": mean_vals(quartet_rows, "D_arrangement_minus_species", lambda c: not inside(c)),
    }


def make_figure(windows, cross_rows, quartet_rows, tree_by_window):
    FIGURES.mkdir(parents=True, exist_ok=True)
    x = [(int(w["start"]) + int(w["end"])) / 2 for w in windows]
    wid = [w["window_id"] for w in windows]
    c_by = {r["window_id"]: float(r["C"]) for r in cross_rows}
    d_by = {r["window_id"]: float(r["D_arrangement_minus_species"]) for r in quartet_rows}
    colors = ["#d95f02" if w["position_class"].startswith("inside") else "#1b9e77" for w in windows]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 8))
        gs = fig.add_gridspec(3, 2, height_ratios=[0.8, 1, 1.5])
        ax = fig.add_subplot(gs[0, :])
        ax.hlines(0, 0, 50_000_000, color="0.8", lw=8)
        ax.hlines(0, INV_START, INV_END, color="#7570b3", lw=12, label="2La")
        for w, xx, col in zip(windows, x, colors):
            ax.vlines(xx, -0.25, 0.25, color=col, lw=3)
        ax.set_xlim(8_000_000, 44_000_000)
        ax.set_yticks([])
        ax.set_xlabel("2L coordinate (AgamP4)")
        ax.set_title("A. fixed Stage-1R pilot windows")
        ax.legend(loc="upper right")
        ax = fig.add_subplot(gs[1, 0])
        ax.axhline(0, color="0.5", lw=1)
        ax.scatter(x, [c_by[i] for i in wid], c=colors)
        ax.set_title("B. crossed distance contrast C(w)")
        ax.set_xlabel("2L coordinate")
        ax.set_ylabel("C")
        ax = fig.add_subplot(gs[1, 1])
        ax.axhline(0, color="0.5", lw=1)
        ax.scatter(x, [d_by[i] for i in wid], c=colors)
        ax.set_title("C. crossed quartet D(w)")
        ax.set_xlabel("2L coordinate")
        ax.set_ylabel("D = q_arrangement - q_species")
        ax = fig.add_subplot(gs[2, :])
        inside_id = "inside_deep_center"
        outside_id = "outside_left_flank_50kb"
        txt = f"D. predetermined NJ trees (Newick excerpt)\n\nOutside control {outside_id}:\n{tree_by_window.get(outside_id, '')[:650]}...\n\nCentral inside {inside_id}:\n{tree_by_window.get(inside_id, '')[:650]}..."
        ax.text(0.01, 0.99, txt, va="top", ha="left", family="monospace", fontsize=7, wrap=True)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(FIG_PDF)
        fig.savefig(FIG_PNG, dpi=200)
        plt.close(fig)
        return
    except ModuleNotFoundError:
        pass

    from PIL import Image, ImageDraw, ImageFont
    W, H = 1800, 1200
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 24)
        small = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 18)
        mono = ImageFont.truetype("/System/Library/Fonts/Menlo.ttc", 14)
    except Exception:
        font = small = mono = ImageFont.load_default()
    def sx(pos):
        return int(120 + (pos - 8_000_000) / (44_000_000 - 8_000_000) * (W - 240))
    draw.text((40, 25), "Stage 1R pilot signal", fill="black", font=font)
    y0 = 130
    draw.text((40, y0 - 55), "A. fixed Stage-1R pilot windows", fill="black", font=font)
    draw.line((sx(8_000_000), y0, sx(44_000_000), y0), fill=(210,210,210), width=12)
    draw.line((sx(INV_START), y0, sx(INV_END), y0), fill=(117,112,179), width=20)
    for w, xx, col in zip(windows, x, colors):
        rgb = (217,95,2) if col == "#d95f02" else (27,158,119)
        draw.line((sx(xx), y0-40, sx(xx), y0+40), fill=rgb, width=5)
    draw.text((sx(INV_START), y0+35), "2La", fill=(117,112,179), font=small)
    def panel(points, title, ytop, key):
        draw.text((40, ytop), title, fill="black", font=font)
        x0, x1 = 120, W-120
        y1, y2 = ytop+70, ytop+320
        vals = [points[i] for i in wid]
        mn, mx = min(vals), max(vals)
        if abs(mx-mn) < 1e-12:
            mn -= 0.01; mx += 0.01
        if mn < 0 < mx:
            yz = int(y2 - (0-mn)/(mx-mn)*(y2-y1))
            draw.line((x0, yz, x1, yz), fill=(150,150,150), width=2)
        draw.rectangle((x0,y1,x1,y2), outline=(0,0,0), width=1)
        for w, xx, col in zip(windows, x, colors):
            rgb = (217,95,2) if col == "#d95f02" else (27,158,119)
            yy = int(y2 - (points[w['window_id']]-mn)/(mx-mn)*(y2-y1))
            draw.ellipse((sx(xx)-7, yy-7, sx(xx)+7, yy+7), fill=rgb)
        draw.text((x0, y2+8), "2L coordinate", fill="black", font=small)
        draw.text((x1-230, y1-28), f"range {mn:.3g} to {mx:.3g}", fill="black", font=small)
    panel(c_by, "B. crossed distance contrast C(w)", 230, "C")
    panel(d_by, "C. crossed quartet D(w)", 580, "D")
    y = 930
    outside_id = "outside_left_flank_50kb"
    inside_id = "inside_deep_center"
    draw.text((40, y), "D. predetermined NJ trees (Newick excerpts)", fill="black", font=font)
    text = f"Outside control {outside_id}:\n{tree_by_window.get(outside_id, '')[:520]}...\n\nCentral inside {inside_id}:\n{tree_by_window.get(inside_id, '')[:520]}..."
    yy = y + 45
    for line in text.splitlines():
        draw.text((50, yy), line, fill="black", font=mono)
        yy += 18
    img.save(FIG_PNG)
    img.save(FIG_PDF, "PDF", resolution=200.0)

def update_readme_and_status(decision: str, signal: dict[str, Any]) -> None:
    block = f"""

## Stage 1R — MalariaGEN regional SNP-access and homozygote pilot

Stage 1R uses authenticated MalariaGEN Ag3 release 3.10 data from `fontaine-2015-rebuild` via `gs://vo_agam_release_master_us_central1`. The primary pilot cohort is the frozen set of 53 2La homozygotes; 19 heterokaryotypes remain in metadata but are excluded from the unphased diploid pilot.

The fixed pilot queried eight predefined 50-kb windows on 2L, used direct regional `snp_calls` and `biallelic_snp_calls`, built IBS-style diploid distances, generated NJ trees for visualization, and scored crossed gambiae/coluzzii quartet support. No raw reads were downloaded and no full 2L scan was performed.

Stage 1R decision: **{decision}**. Mean pilot C inside/outside: {signal['mean_C_inside']:.6g} / {signal['mean_C_outside']:.6g}. Mean pilot D inside/outside: {signal['mean_D_inside']:.6g} / {signal['mean_D_outside']:.6g}.
"""
    for path in [README, PROJECT_STATUS]:
        old = path.read_text() if path.exists() else ""
        marker = "## Stage 1R — MalariaGEN regional SNP-access and homozygote pilot"
        if marker in old:
            old = old.split(marker)[0].rstrip() + "\n"
        path.write_text(old.rstrip() + block)


def write_report(decision: str, signal: dict[str, Any], qc_rows, cross_rows, quartet_rows, hap_rows, env):
    hap_available = sum(1 for r in hap_rows if r["phased_available"])
    qc_table = "\n".join(f"- {r['window_id']}: {r['n_post_filter_snps']} post-filter SNPs, mean missingness {float(r['mean_missingness']):.4f}, status {r['status']}" for r in qc_rows)
    cross_table = "\n".join(f"- {r['window_id']}: C={float(r['C']):.6g}" for r in cross_rows)
    quart_table = "\n".join(f"- {r['window_id']}: D={float(r['D_arrangement_minus_species']):.6g}, q_arr={float(r['q_arrangement']):.3f}, q_species={float(r['q_species']):.3f}" for r in quartet_rows)
    REPORT.write_text(f"""# Stage 1R — MalariaGEN regional SNP-access and homozygote pilot

Stage 1R tested whether the authenticated MalariaGEN sample-level data are technically usable for a 2La arrangement-vs-species analysis before any full 2L scan. It used the frozen Stage 1A sample and karyotype data, the fixed AgamP4 interval `{CHROM}:{INV_START}-{INV_END}`, and eight coordinate-predefined 50-kb pilot windows. No raw reads were downloaded, no full-chromosome scan was run, and no final manuscript-wide statistic was calculated.

## Environment

- Python: `{env['python_version']}`
- executable: `{env['python_executable']}`
- Python 3.11 migration: {env['python_3_11_migration']}
- malariagen-data: `{env['malariagen_data_version']}`
- Ag3 release: `{RELEASE}`
- sample set: `{SAMPLE_SET}`
- bucket: `{BUCKET}`

## Primary cohort

The primary Stage 1R cohort contains exactly 53 homozygotes. The 19 heterokaryotypes remain in metadata and are excluded from the primary unphased diploid analysis.

## SNP QC

{qc_table}

## Crossed gambiae-coluzzii pilot C(w)

`C(w) = mean same-species/opposite-arrangement distance - mean different-species/same-arrangement distance`. Positive values indicate that arrangement similarity can overcome species identity in this pilot contrast.

{cross_table}

## Crossed quartet pilot D(w)

`D(w)=q_arrangement-q_species` from all gambiae-standard, gambiae-inverted, coluzzii-standard, coluzzii-inverted individual quartets.

{quart_table}

Mean C inside/outside: `{signal['mean_C_inside']:.6g}` / `{signal['mean_C_outside']:.6g}`.
Mean D inside/outside: `{signal['mean_D_inside']:.6g}` / `{signal['mean_D_outside']:.6g}`.

## Haplotype membership

Regional haplotype availability was checked but not used as the primary Stage 1R analysis. `{hap_available}` of 72 Fontaine-associated samples were returned by the regional haplotypes query.

## Interpretation

The primary biological inference for Stage 1R is restricted to the crossed gambiae + coluzzii design, because those two species contain both homozygous arrangements. Arabiensis and quadriannulatus are useful anchors in all-53 visualizations, but arrangement and species are fully confounded within each, so they are not allowed to drive the claim that arrangement overrides species.

The pilot does not claim formal significance. It asks whether direct regional SNP access works, whether the fixed homozygote cohort produces usable distances and NJ visualizations, and whether inside windows show a qualitatively stronger arrangement-associated signal than outside controls.

Decision: **{decision}**.

Stage 2R has {'been frozen in `stage2r_preanalysis_plan.md`' if decision in {'PASS', 'TECHNICAL PASS / BIOLOGICAL UNCERTAIN'} else 'not been frozen because the pilot failed'}.
""")


def write_stage2_plan(decision: str) -> None:
    if decision not in {"PASS", "TECHNICAL PASS / BIOLOGICAL UNCERTAIN"}:
        return
    STAGE2_PLAN.write_text(f"""# Stage 2R preanalysis plan — full 2L 50-kb scan, not executed

Stage 2R will use the same authenticated MalariaGEN Ag3 release `{RELEASE}` data resource and the same frozen 53 homozygous Fontaine-associated samples from Stage 1R.

Predeclared design:

- construct non-overlapping 50-kb windows across chromosome arm 2L in AgamP4 coordinates;
- classify windows fully inside frozen 2La interval `{CHROM}:{INV_START}-{INV_END}` as inside;
- classify windows fully outside the interval as collinear outside controls;
- exclude any boundary-overlap window from inside/outside tests;
- use the same biallelic SNP filter as Stage 1R unless a change is documented before analysis: minor allele count >= {MIN_MINOR_AC}, site missingness <= {MAX_SITE_MISSING:.2f}, no LD pruning for the primary distance/tree scan;
- mark low-SNP windows as missing rather than replacing them after inspecting topology;
- primary statistic: crossed gambiae/coluzzii `C(w)` and quartet `D(w)` using only samples from the four crossed classes;
- secondary analysis: all-53 descriptive distances and NJ visualization summaries;
- spatial null/control: compare inside 2La windows against all eligible collinear 2L windows, with sensitivity to matched controls by SNP density/callability only if the matching rule is fixed before topology inspection;
- fixed seed `20261003` for any deterministic quartet subsampling, although exhaustive crossed quartets are preferred;
- no raw-read processing;
- no sample removal based on local topology.

Stage 2R is not executed in Stage 1R.
""")


def run(args):
    import malariagen_data
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    primary_all, source_rows = load_primary_rows()
    validate_primary(primary_all, source_rows)
    primary = [r for r in primary_all if r["include_primary"]]
    write_tsv(PRIMARY_SAMPLES, primary_all, ["sample_id", "species", "karyotype_numeric", "arrangement", "sample_set", "study_id", "country", "population", "include_primary", "exclusion_reason"])
    windows = fixed_windows()
    validate_windows(windows)
    write_tsv(PILOT_WINDOWS, windows, ["window_id", "chrom", "start", "end", "position_class", "distance_to_nearest_breakpoint", "selection_rule"])
    ag3 = malariagen_data.Ag3(url=args.url, check_location=False, show_progress=False)
    migration_note = "not attempted beyond PATH probe: no python3.11 executable available; Stage 1R used validated Python 3.10 environment"
    env = env_provenance(ag3, migration_note)
    ENV_PROV.write_text(json.dumps(env, indent=2, sort_keys=True) + "\n")
    md = ag3.sample_metadata(sample_sets=SAMPLE_SET)
    md_ids = [str(x) for x in md["sample_id"].tolist()]
    id_to_idx = {sid: i for i, sid in enumerate(md_ids)}
    sample_indices = [id_to_idx[s["sample_id"]] for s in primary]
    assert len(sample_indices) == 53 and len(set(sample_indices)) == 53
    qc_rows = []
    all_pair_rows = []
    class_rows = []
    cross_rows = []
    quartet_rows = []
    tree_lines = []
    tree_by_window = {}
    for w in windows:
        result = query_window(ag3, w, sample_indices, primary)
        qc_rows.append(result["qc"])
        all_pair_rows.extend(result["pairs"])
        class_rows.extend(class_summaries(w, result["pairs"]))
        cross_rows.append(crossed_contrast(w, result["pairs"]))
        quartet_rows.append(quartet_support(w, primary, result["dist"]))
        tree_by_window[w["window_id"]] = result["newick"]
        tree_lines.append(f"{w['window_id']}\t{result['newick']}")
    write_tsv(WINDOW_QC, qc_rows, ["window_id", "position_class", "n_samples_requested", "n_samples_returned", "n_sites_raw", "n_biallelic_snps", "n_segregating_biallelic_snps", "n_post_filter_snps", "mean_missingness", "max_sample_missingness", "mean_heterozygosity", "callable", "status", "notes", "runtime_seconds", "approx_dataset_bytes", "mean_post_filter_sample_missingness", "alt_ac_mean", "minor_ac_mean"])
    write_filter_plan(qc_rows)
    write_tsv(PAIRWISE, all_pair_rows, ["window_id", "sample1", "sample2", "species1", "species2", "arrangement1", "arrangement2", "distance", "n_sites_compared"])
    TREES.write_text("\n".join(tree_lines) + "\n")
    write_tsv(CROSSED, cross_rows, ["window_id", "position_class", "n_within_species_opposite_pairs", "mean_within_species_opposite_distance", "n_cross_species_same_pairs", "mean_cross_species_same_distance", "C"])
    write_tsv(CLASS_SUMMARY, class_rows, ["window_id", "position_class", "distance_class", "n_pairs", "mean_distance", "median_distance", "min_distance", "max_distance", "mean_sites_compared"])
    write_tsv(QUARTETS, quartet_rows, ["window_id", "position_class", "n_informative_quartets", "q_species", "q_arrangement", "q_third", "D_arrangement_minus_species"])
    hap_rows = haplotype_membership(ag3, source_rows)
    write_tsv(HAP_MEMBERSHIP, hap_rows, ["sample_id", "species", "arrangement", "haplotype_panel", "phased_available", "notes"])
    validate_outputs(primary_all, windows, qc_rows, all_pair_rows, quartet_rows)
    signal = summarize_signal(cross_rows, quartet_rows)
    technical_ok = all(r["status"] == "ok" for r in qc_rows)
    inside_stronger = signal["mean_C_inside"] > signal["mean_C_outside"] and signal["mean_D_inside"] > signal["mean_D_outside"]
    if technical_ok and inside_stronger:
        decision = "PASS"
    elif technical_ok:
        decision = "TECHNICAL PASS / BIOLOGICAL UNCERTAIN"
    else:
        decision = "FAIL"
    make_figure(windows, cross_rows, quartet_rows, tree_by_window)
    write_report(decision, signal, qc_rows, cross_rows, quartet_rows, hap_rows, env)
    write_stage2_plan(decision)
    update_readme_and_status(decision, signal)
    outputs = [PRIMARY_SAMPLES, PILOT_WINDOWS, PAIRWISE, TREES, ENV_PROV, WINDOW_QC, FILTER_PLAN, CROSSED, CLASS_SUMMARY, QUARTETS, HAP_MEMBERSHIP, REPORT, FIG_PDF, FIG_PNG]
    if STAGE2_PLAN.exists():
        outputs.append(STAGE2_PLAN)
    manifest = {
        "stage": "Stage 1R",
        "decision": decision,
        "no_full_scan": True,
        "no_raw_read_processing": True,
        "stage2r_executed": False,
        "total_pilot_bp": len(windows) * WIN,
        "outputs": {str(p.relative_to(REPO)): sha256(p) for p in outputs if p.exists()},
        "signal_summary": signal,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return 0


def validate_outputs(primary_all, windows, qc_rows, pair_rows, quartet_rows):
    validate_primary(primary_all, read_tsv(SAMPLE_MANIFEST))
    validate_windows(windows)
    assert all(r["n_samples_returned"] == 53 for r in qc_rows), "not all requested samples returned"
    assert len({r["window_id"] for r in qc_rows}) == 8
    d: dict[tuple[str, str], float] = {}
    for r in pair_rows:
        val = float(r["distance"])
        assert val >= 0 or math.isnan(val), "negative distance"
        assert r["sample1"] != r["sample2"], "self-pair written"
    for r in quartet_rows:
        qsum = float(r["q_species"]) + float(r["q_arrangement"]) + float(r["q_third"])
        assert abs(qsum - 1.0) < 1e-9, f"quartet fractions do not sum to 1 for {r['window_id']}"
    assert len(windows) * WIN <= NO_FULL_SCAN_MAX_TOTAL_BP, "full-chromosome scan suspected"
    assert not any(str(p).endswith((".fastq", ".fq", ".bam", ".cram")) for p in RESULTS.glob("**/*")), "raw-read processing output detected in Stage 1R results"


class Stage1RTests(unittest.TestCase):
    def test_fixed_windows_do_not_overlap_boundaries(self):
        validate_windows(fixed_windows())
    def test_primary_counts(self):
        primary, source = load_primary_rows()
        validate_primary(primary, source)
    def test_distance_formula(self):
        gt = np.array([[[0,0],[1,1]], [[0,1],[0,1]], [[-1,-1],[1,1]]], dtype=np.int8)
        dosage, called, alt_ac, minor_ac, keep = dosage_and_filters(gt)
        keep[:] = True
        dist, comp = distance_matrix(dosage, called, keep)
        self.assertEqual(comp[0,1], 2)
        self.assertAlmostEqual(dist[0,1], 0.5)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--url", default=BUCKET)
    p.add_argument("--run-tests", action="store_true")
    args = p.parse_args(argv)
    if args.run_tests:
        res = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Stage1RTests))
        return 0 if res.wasSuccessful() else 1
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
