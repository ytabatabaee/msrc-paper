#!/usr/bin/env python3
"""Build chromosome-4 local QQS tracks for focal Neoaves clades."""

from __future__ import annotations

import argparse
import csv
import lzma
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from statistics import mean
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
META_PATH = REPO_ROOT / "genetreesupport" / "63K_trees.names_header.txt.xz"
REC_PATH = REPO_ROOT / "genetreesupport" / "clade-rec.stat.xz"
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"
TABLE_PATH = RESULTS_DIR / "locus_table_chr4.tsv"
VALIDATION_PATH = RESULTS_DIR / "validation_summary.tsv"
MEANS_PATH = RESULTS_DIR / "chr4_clade_means.tsv"
PDF_PATH = FIGURES_DIR / "chr4_raw_quartet_tracks.pdf"
PNG_PATH = FIGURES_DIR / "chr4_raw_quartet_tracks.png"

FOCAL_CLADE_MAP = {
    ("Columbimorphae", "Phoenicopteriformes"): "Columbea",
    ("Otidimorphae", "Columbimorphae"): "N61",
    ("Columbiformes", "OtherColumbimorphae"): "N62",
}
FOCAL_ORDER = ["Columbea", "N61", "N62"]

TOPOLOGY_LABELS = {
    "q1": "(C1,C2)|(S,O)",
    "q2": "(C1,S)|(C2,O)",
    "q3": "(C2,S)|(C1,O)",
}

OUT_COLUMNS = [
    "Gene",
    "clade",
    "C1",
    "C2",
    "Chromosome",
    "start",
    "end",
    "midpoint",
    "q1",
    "q2",
    "q3",
    "dominant_topology",
    "dominant_support",
    "x1",
    "x2",
    "x3",
    "x4",
    "d1",
    "d2",
    "d3",
]


def parse_chr(chromosome: str) -> str:
    """Reproduce draw-movingaverage.r chromosome parsing for factor labels."""
    stem = chromosome.split("_", 1)[0].replace("chr", "", 1)
    try:
        return str(int(float(stem)))
    except ValueError:
        return stem


def compute_qqs(x1: float, x2: float, x3: float, x4: float) -> dict[str, float] | None:
    d1 = x1 - x4
    d2 = x2 - x4
    d3 = x3 - x4
    denom = d1 + d2 + d3
    if denom == 0 or not math.isfinite(denom):
        return None
    return {
        "d1": d1,
        "d2": d2,
        "d3": d3,
        "q1": (d2 + d3 - d1) / denom,
        "q2": (d1 + d3 - d2) / denom,
        "q3": (d1 + d2 - d3) / denom,
    }


def fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def load_chr4_metadata(meta_path: Path) -> dict[str, dict[str, object]]:
    metadata = {}
    with lzma.open(meta_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        required = {"Gene", "Chromosome", "ws", "we"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{meta_path} is missing required columns: {sorted(missing)}")
        for row in reader:
            if row["Chromosome"] != "chr4":
                continue
            start = int(row["ws"])
            end = int(row["we"])
            metadata[row["Gene"]] = {
                "Gene": row["Gene"],
                "Chromosome": row["Chromosome"],
                "start": start,
                "end": end,
                "midpoint": (start + end) / 2,
            }
    return metadata


def load_focal_chr4_xvalues(
    rec_path: Path, chr4_genes: set[str]
) -> dict[tuple[str, str], dict[str, object]]:
    values = defaultdict(lambda: {"x": {}})
    with lzma.open(rec_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) != 5:
                raise ValueError(f"Unexpected field count at {rec_path}:{line_number}: {line!r}")
            c1, c2, topology_id, gene, x = parts
            clade = FOCAL_CLADE_MAP.get((c1, c2))
            if clade is None or gene not in chr4_genes:
                continue
            if topology_id not in {"1", "2", "3", "4"}:
                raise ValueError(f"Unexpected topology_id {topology_id!r} at line {line_number}")
            key = (gene, clade)
            values[key]["Gene"] = gene
            values[key]["clade"] = clade
            values[key]["C1"] = c1
            values[key]["C2"] = c2
            values[key]["x"][topology_id] = float(x)
    return values


def build_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    metadata = load_chr4_metadata(META_PATH)
    xvalues = load_focal_chr4_xvalues(REC_PATH, set(metadata))

    rows = []
    undefined = 0
    missing_topologies = 0
    for (gene, _clade), record in xvalues.items():
        xs = record["x"]
        if set(xs) != {"1", "2", "3", "4"}:
            missing_topologies += 1
            continue
        q = compute_qqs(xs["1"], xs["2"], xs["3"], xs["4"])
        if q is None:
            undefined += 1
            continue
        supports = {"q1": q["q1"], "q2": q["q2"], "q3": q["q3"]}
        dominant_topology, dominant_support = max(
            supports.items(), key=lambda item: (item[1], item[0])
        )
        rows.append(
            {
                **{k: record[k] for k in ("Gene", "clade", "C1", "C2")},
                **metadata[gene],
                "q1": q["q1"],
                "q2": q["q2"],
                "q3": q["q3"],
                "dominant_topology": dominant_topology,
                "dominant_support": dominant_support,
                "x1": xs["1"],
                "x2": xs["2"],
                "x3": xs["3"],
                "x4": xs["4"],
                "d1": q["d1"],
                "d2": q["d2"],
                "d3": q["d3"],
            }
        )

    rows.sort(key=lambda r: (FOCAL_ORDER.index(r["clade"]), r["start"], r["end"], int(r["Gene"])))
    stats = validate_rows(rows, undefined, missing_topologies)
    return rows, stats


def validate_rows(
    rows: list[dict[str, object]], undefined: int, missing_topologies: int
) -> dict[str, object]:
    counts = {clade: 0 for clade in FOCAL_ORDER}
    for row in rows:
        counts[row["clade"]] += 1

    sum_deviations = [
        abs(float(row["q1"]) + float(row["q2"]) + float(row["q3"]) - 1.0) for row in rows
    ]
    max_sum_deviation = max(sum_deviations, default=0.0)
    outside_values = sum(
        q < -1e-12 or q > 1 + 1e-12
        for row in rows
        for q in (float(row["q1"]), float(row["q2"]), float(row["q3"]))
    )
    duplicated = len(rows) - len({(row["Gene"], row["clade"]) for row in rows})
    coord_min = min(int(row["start"]) for row in rows)
    coord_max = max(int(row["end"]) for row in rows)

    if max_sum_deviation > 1e-10:
        raise ValueError(
            f"QQS probabilities fail sum-to-one tolerance: max deviation {max_sum_deviation}"
        )

    return {
        "chr4_locus_counts": counts,
        "undefined_qqs_dropped": undefined,
        "missing_topology_rows_dropped": missing_topologies,
        "max_abs_q_sum_minus_1": max_sum_deviation,
        "outside_0_1_values": outside_values,
        "duplicate_gene_clade_rows": duplicated,
        "coordinate_min_start": coord_min,
        "coordinate_max_end": coord_max,
    }


def write_table(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: fmt(row[column]) for column in OUT_COLUMNS})


def summarize_means(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summaries = []
    for clade in FOCAL_ORDER:
        subset = [row for row in rows if row["clade"] == clade]
        qmeans = {q: mean(float(row[q]) for row in subset) for q in ("q1", "q2", "q3")}
        dominant = max(qmeans.items(), key=lambda item: (item[1], item[0]))[0]
        summaries.append(
            {
                "clade": clade,
                "n_loci": len(subset),
                "mean_q1": qmeans["q1"],
                "mean_q2": qmeans["q2"],
                "mean_q3": qmeans["q3"],
                "chr4_mean_dominant_topology": dominant,
                "chr4_mean_dominant_support": qmeans[dominant],
            }
        )
    return summaries


def write_summaries(stats: dict[str, object], means: list[dict[str, object]]) -> None:
    with VALIDATION_PATH.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for clade, count in stats["chr4_locus_counts"].items():
            writer.writerow([f"chr4_loci_{clade}", count])
        for key in (
            "undefined_qqs_dropped",
            "missing_topology_rows_dropped",
            "max_abs_q_sum_minus_1",
            "outside_0_1_values",
            "duplicate_gene_clade_rows",
            "coordinate_min_start",
            "coordinate_max_end",
        ):
            writer.writerow([key, fmt(stats[key])])

    with MEANS_PATH.open("w", newline="") as handle:
        columns = [
            "clade",
            "n_loci",
            "mean_q1",
            "mean_q2",
            "mean_q3",
            "chr4_mean_dominant_topology",
            "chr4_mean_dominant_support",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in means:
            writer.writerow({column: fmt(row[column]) for column in columns})


def rolling_mean(points: list[tuple[float, float]], window: int = 50) -> tuple[list[float], list[float]]:
    if not points:
        return [], []
    xs, ys = zip(*points)
    half = window // 2
    smoothed = []
    for i in range(len(ys)):
        lo = max(0, i - half)
        hi = min(len(ys), i + half + 1)
        smoothed.append(sum(ys[lo:hi]) / (hi - lo))
    return list(xs), smoothed


def make_plot(rows: list[dict[str, object]]) -> None:
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"q1": "#222222", "q2": "#1f78b4", "q3": "#e31a1c"}
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, sharey=True)
    for axis, clade in zip(axes, FOCAL_ORDER):
        subset = [row for row in rows if row["clade"] == clade]
        for q in ("q1", "q2", "q3"):
            points = [(float(row["midpoint"]), float(row[q])) for row in subset]
            xs = [point[0] for point in points]
            ys = [point[1] for point in points]
            axis.scatter(xs, ys, s=6, alpha=0.18, color=colors[q], linewidths=0)
            smooth_x, smooth_y = rolling_mean(points, window=50)
            axis.plot(
                smooth_x,
                smooth_y,
                color=colors[q],
                linewidth=1.2,
                label=f"{q}: {TOPOLOGY_LABELS[q]} rolling mean (n=50 loci)",
            )
        axis.axhline(1 / 3, color="#33a02c", linestyle=":", linewidth=1.0)
        axis.set_title(clade, loc="left", fontsize=11)
        axis.set_ylabel("QQS")
        axis.set_ylim(-0.03, 1.03)
        axis.grid(axis="y", color="#dddddd", linewidth=0.5)
    axes[-1].set_xlabel("Chromosome 4 position (midpoint bp; ws/we from metadata)")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=220)
    plt.close(fig)


def print_report(stats: dict[str, object], means: list[dict[str, object]]) -> None:
    print("Chr4 focal locus counts:")
    for clade, count in stats["chr4_locus_counts"].items():
        print(f"  {clade}: {count}")
    print(
        "Coordinate range: "
        f"{stats['coordinate_min_start']}..{stats['coordinate_max_end']} bp"
    )
    print("QQS validation:")
    for key in (
        "undefined_qqs_dropped",
        "missing_topology_rows_dropped",
        "max_abs_q_sum_minus_1",
        "outside_0_1_values",
        "duplicate_gene_clade_rows",
    ):
        print(f"  {key}: {fmt(stats[key])}")
    print("Mean QQS and chr4 mean dominant topology:")
    for row in means:
        print(
            f"  {row['clade']}: "
            f"mean q1={row['mean_q1']:.6f}, q2={row['mean_q2']:.6f}, "
            f"q3={row['mean_q3']:.6f}; dominant={row['chr4_mean_dominant_topology']} "
            f"({row['chr4_mean_dominant_support']:.6f})"
        )
    print("Files written:")
    for path in (TABLE_PATH, VALIDATION_PATH, MEANS_PATH, PDF_PATH, PNG_PATH):
        print(f"  {path.relative_to(REPO_ROOT)}")


class QQSTests(unittest.TestCase):
    def test_synthetic_sum_to_one(self) -> None:
        result = compute_qqs(x1=9.0, x2=8.0, x3=7.0, x4=3.0)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["q1"] + result["q2"] + result["q3"], 1.0)

    def test_expected_supports_from_distances(self) -> None:
        # a1=5, a2=3, a3=2 implies d1=5, d2=7, d3=8 and QQS=(0.5,0.3,0.2).
        result = compute_qqs(x1=6.0, x2=8.0, x3=9.0, x4=1.0)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["q1"], 0.5)
        self.assertAlmostEqual(result["q2"], 0.3)
        self.assertAlmostEqual(result["q3"], 0.2)

    def test_undefined_zero_denominator(self) -> None:
        self.assertIsNone(compute_qqs(x1=4.0, x2=4.0, x3=4.0, x4=4.0))


def run_tests() -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(QQSTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run lightweight QQS tests first")
    args = parser.parse_args()

    if args.run_tests:
        run_tests()

    rows, stats = build_rows()
    means = summarize_means(rows)
    write_table(rows, TABLE_PATH)
    write_summaries(stats, means)
    make_plot(rows)
    print_report(stats, means)


if __name__ == "__main__":
    main()
