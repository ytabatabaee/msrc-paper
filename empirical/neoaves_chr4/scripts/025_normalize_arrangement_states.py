#!/usr/bin/env python3
"""Normalize raw chr4 synteny signatures into comparable arrangement states."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"

BLOCKS_PATH = RESULTS_DIR / "structural_blocks.tsv"
INTERVALS_PATH = RESULTS_DIR / "structural_intervals.tsv"
RAW_STATES_PATH = RESULTS_DIR / "arrangement_states.tsv"

NORMALIZED_STATES_PATH = RESULTS_DIR / "normalized_arrangement_states.tsv"
NORMALIZED_DEFS_PATH = RESULTS_DIR / "normalized_arrangement_state_definitions.tsv"
PAIRWISE_PATH = RESULTS_DIR / "arrangement_state_pairwise_similarity.tsv"
PARTITIONS_PATH = RESULTS_DIR / "structural_partitions.tsv"
VALIDATION_PATH = RESULTS_DIR / "normalized_arrangement_validation.tsv"
PDF_PATH = FIGURES_DIR / "chr4_normalized_arrangement_states.pdf"
PNG_PATH = FIGURES_DIR / "chr4_normalized_arrangement_states.png"

FOCAL_SPECIES = ["Stork", "Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo"]
PLOT_ORDER = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]
BOUNDARY_TOLERANCE_BP = 50_000
MIN_ANCHOR_SIZE_BP = 100_000
MIN_ANCHOR_OVERLAP_BP = 1
LOW_SIMILARITY_THRESHOLD = 0.50


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: fmt(row.get(column, "")) for column in columns})


def fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def as_int(row: dict[str, object], key: str) -> int:
    return int(float(str(row[key])))


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def build_anchors(
    interval: dict[str, str], blocks_by_species: dict[str, list[dict[str, str]]]
) -> list[dict[str, object]]:
    start = as_int(interval, "start")
    end = as_int(interval, "end")
    boundaries = [start, end]
    for species in FOCAL_SPECIES:
        for block in blocks_by_species.get(species, []):
            b_start = as_int(block, "reference_start")
            b_end = as_int(block, "reference_end")
            if overlap(start, end, b_start, b_end) <= 0:
                continue
            boundaries.append(max(start, b_start))
            boundaries.append(min(end, b_end))
    merged = merge_boundaries(boundaries, start, end)
    anchors = []
    for i, (anchor_start, anchor_end) in enumerate(zip(merged, merged[1:]), 1):
        if anchor_end - anchor_start < MIN_ANCHOR_SIZE_BP:
            continue
        anchors.append(
            {
                "anchor_id": f"{interval['structural_interval_id']}_ANCHOR_{len(anchors) + 1:03d}",
                "anchor_index": len(anchors) + 1,
                "start": anchor_start,
                "end": anchor_end,
                "length_bp": anchor_end - anchor_start,
            }
        )
    if not anchors and end > start:
        anchors.append(
            {
                "anchor_id": f"{interval['structural_interval_id']}_ANCHOR_001",
                "anchor_index": 1,
                "start": start,
                "end": end,
                "length_bp": end - start,
            }
        )
    return anchors


def merge_boundaries(boundaries: list[int], start: int, end: int) -> list[int]:
    values = sorted(set(boundaries))
    clusters: list[list[int]] = []
    for value in values:
        if not clusters or value - clusters[-1][0] > BOUNDARY_TOLERANCE_BP:
            clusters.append([value])
        else:
            clusters[-1].append(value)
    merged = []
    for cluster in clusters:
        if start in cluster:
            merged.append(start)
        elif end in cluster:
            merged.append(end)
        else:
            merged.append(round(sum(cluster) / len(cluster)))
    return sorted(set([start, end, *merged]))


def project_species(
    species: str, anchors: list[dict[str, object]], blocks: list[dict[str, str]]
) -> dict[str, object]:
    observations = []
    for anchor in anchors:
        anchor_start = as_int(anchor, "start")
        anchor_end = as_int(anchor, "end")
        hits = []
        for block in blocks:
            hit = overlap(
                anchor_start,
                anchor_end,
                as_int(block, "reference_start"),
                as_int(block, "reference_end"),
            )
            if hit >= MIN_ANCHOR_OVERLAP_BP:
                hits.append((hit, block))
        if not hits:
            observations.append({"status": "missing", "anchor": anchor})
            continue
        chroms = {hit[1]["query_chrom"] for hit in hits}
        orientations = {hit[1]["orientation"] for hit in hits}
        if len(chroms) > 1 or len(orientations) > 1:
            best = max(hits, key=lambda item: item[0])[1]
            observations.append(
                {
                    "status": "ambiguous",
                    "anchor": anchor,
                    "query_chrom": best["query_chrom"],
                    "query_midpoint": float(best["query_midpoint"]),
                    "orientation": best["orientation"],
                }
            )
            continue
        best = max(hits, key=lambda item: item[0])[1]
        observations.append(
            {
                "status": "observed",
                "anchor": anchor,
                "query_chrom": best["query_chrom"],
                "query_midpoint": float(best["query_midpoint"]),
                "orientation": best["orientation"],
            }
        )

    return build_signatures(species, observations)


def build_signatures(species: str, observations: list[dict[str, object]]) -> dict[str, object]:
    observed = [obs for obs in observations if obs["status"] == "observed"]
    chrom_alias: dict[str, str] = {}
    for obs in observed:
        chrom = str(obs["query_chrom"])
        if chrom not in chrom_alias:
            chrom_alias[chrom] = f"Q{len(chrom_alias) + 1}"

    rank_by_anchor = {}
    by_chrom = defaultdict(list)
    for obs in observed:
        by_chrom[str(obs["query_chrom"])].append(obs)
    for chrom, chrom_obs in by_chrom.items():
        ordered = sorted(chrom_obs, key=lambda obs: float(obs["query_midpoint"]))
        for rank, obs in enumerate(ordered, 1):
            rank_by_anchor[int(obs["anchor"]["anchor_index"])] = rank

    order_parts = []
    orientation_parts = []
    compatibility_tokens = {}
    for obs in observations:
        anchor_index = int(obs["anchor"]["anchor_index"])
        if obs["status"] == "missing":
            order_token = "M"
            orientation_token = "?"
        elif obs["status"] == "ambiguous":
            order_token = "X"
            orientation_token = "!"
        else:
            alias = chrom_alias[str(obs["query_chrom"])]
            order_token = f"{alias}:{rank_by_anchor[anchor_index]}"
            orientation_token = str(obs["orientation"])
            compatibility_tokens[anchor_index] = (order_token, orientation_token)
        order_parts.append(f"{anchor_index}:{order_token}")
        orientation_parts.append(f"{anchor_index}:{orientation_token}")

    n = len(observations) or 1
    missing = sum(obs["status"] == "missing" for obs in observations)
    ambiguous = sum(obs["status"] == "ambiguous" for obs in observations)
    return {
        "species": species,
        "anchor_order_signature": "|".join(order_parts),
        "anchor_orientation_signature": "|".join(orientation_parts),
        "n_shared_anchors": len(observations),
        "missing_anchor_fraction": missing / n,
        "ambiguous_anchor_fraction": ambiguous / n,
        "compatibility_tokens": compatibility_tokens,
        "anchor_observations": observations,
    }


def compatible(a: dict[str, object], b: dict[str, object]) -> bool:
    a_tokens = a["compatibility_tokens"]
    b_tokens = b["compatibility_tokens"]
    shared = set(a_tokens) & set(b_tokens)
    if not shared:
        return False
    return all(a_tokens[index] == b_tokens[index] for index in shared)


def assign_states(records: list[dict[str, object]]) -> None:
    state_members: list[list[dict[str, object]]] = []
    for record in sorted(records, key=lambda row: FOCAL_SPECIES.index(str(row["species"]))):
        placed = False
        for members in state_members:
            if all(compatible(record, member) for member in members):
                members.append(record)
                placed = True
                break
        if not placed:
            state_members.append([record])
    for i, members in enumerate(state_members):
        state = f"A{i}"
        for member in members:
            member["normalized_arrangement_state"] = state
            member["equivalence_reason"] = (
                "exact observed shared-anchor order/orientation compatibility; "
                "missing anchors ignored; ambiguous anchors not used for grouping"
            )


def pairwise_similarity(interval_id: str, records: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    by_species = {str(row["species"]): row for row in records}
    for i, species1 in enumerate(FOCAL_SPECIES):
        for species2 in FOCAL_SPECIES[i + 1 :]:
            first = by_species[species1]
            second = by_species[species2]
            first_tokens = first["compatibility_tokens"]
            second_tokens = second["compatibility_tokens"]
            shared = sorted(set(first_tokens) & set(second_tokens))
            shared_count = len(shared)
            if shared_count:
                order_match = sum(first_tokens[idx][0] == second_tokens[idx][0] for idx in shared)
                orient_match = sum(first_tokens[idx][1] == second_tokens[idx][1] for idx in shared)
                order_concordance = order_match / shared_count
                orientation_concordance = orient_match / shared_count
                overall = (order_concordance + orientation_concordance) / 2
            else:
                order_concordance = orientation_concordance = overall = 0.0
            rows.append(
                {
                    "structural_interval_id": interval_id,
                    "species1": species1,
                    "species2": species2,
                    "shared_anchor_count": shared_count,
                    "order_concordance": order_concordance,
                    "orientation_concordance": orientation_concordance,
                    "overall_similarity": overall,
                    "same_normalized_state": str(
                        first["normalized_arrangement_state"]
                        == second["normalized_arrangement_state"]
                    ).lower(),
                }
            )
    return rows


def normalize_states() -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    blocks = read_tsv(BLOCKS_PATH)
    intervals = read_tsv(INTERVALS_PATH)
    raw_states = read_tsv(RAW_STATES_PATH)
    blocks_by_species = defaultdict(list)
    for block in blocks:
        if block["species"] in FOCAL_SPECIES:
            blocks_by_species[block["species"]].append(block)
    raw_by_species_interval = {
        (row["species"], row["structural_interval_id"]): row["arrangement_state"] for row in raw_states
    }

    normalized_rows = []
    definition_rows = []
    pairwise_rows = []
    partition_rows = []
    validation_rows = []

    for interval in intervals:
        anchors = build_anchors(interval, blocks_by_species)
        records = []
        for species in FOCAL_SPECIES:
            projected = project_species(species, anchors, blocks_by_species[species])
            raw_label = raw_by_species_interval.get((species, interval["structural_interval_id"]), "")
            if not raw_label:
                raw_label = overlapping_raw_labels(species, interval, raw_states, intervals)
            projected.update(
                {
                    "structural_interval_id": interval["structural_interval_id"],
                    "raw_state_label": raw_label,
                }
            )
            records.append(projected)
        assign_states(records)

        state_groups = defaultdict(list)
        for record in records:
            state_groups[record["normalized_arrangement_state"]].append(record)
            normalized_rows.append(
                {
                    "species": record["species"],
                    "structural_interval_id": interval["structural_interval_id"],
                    "normalized_arrangement_state": record["normalized_arrangement_state"],
                    "anchor_order_signature": record["anchor_order_signature"],
                    "anchor_orientation_signature": record["anchor_orientation_signature"],
                    "n_shared_anchors": record["n_shared_anchors"],
                    "missing_anchor_fraction": record["missing_anchor_fraction"],
                    "ambiguous_anchor_fraction": record["ambiguous_anchor_fraction"],
                    "raw_state_label": record["raw_state_label"],
                    "equivalence_reason": record["equivalence_reason"],
                }
            )
        for state, members in sorted(state_groups.items()):
            species_members = [str(member["species"]) for member in members]
            partition_rows.append(
                {
                    "structural_interval_id": interval["structural_interval_id"],
                    "state": state,
                    "species_members": ",".join(species_members),
                    "n_species": len(species_members),
                }
            )
            example = members[0]
            definition_rows.append(
                {
                    "structural_interval_id": interval["structural_interval_id"],
                    "normalized_arrangement_state": state,
                    "species_members": ",".join(species_members),
                    "n_species": len(species_members),
                    "n_shared_anchors": example["n_shared_anchors"],
                    "definition": (
                        f"State {state} over {interval['structural_interval_id']} has "
                        f"{len(species_members)} species with compatible observed anchor "
                        f"order/orientation. Missing anchors are explicit and ignored for grouping."
                    ),
                    "anchor_order_signature": example["anchor_order_signature"],
                    "anchor_orientation_signature": example["anchor_orientation_signature"],
                }
            )
        pairwise_rows.extend(pairwise_similarity(interval["structural_interval_id"], records))
        raw_count = len({row["raw_state_label"] for row in records if row["raw_state_label"]})
        norm_count = len(state_groups)
        validation_rows.append(
            {
                "structural_interval_id": interval["structural_interval_id"],
                "n_anchors": len(anchors),
                "raw_state_count": raw_count,
                "normalized_state_count": norm_count,
                "n_species_in_shared_normalized_states": sum(
                    len(members) for members in state_groups.values() if len(members) > 1
                ),
                "all_species_unique": str(norm_count == len(FOCAL_SPECIES)).lower(),
                "max_missing_anchor_fraction": max(float(row["missing_anchor_fraction"]) for row in records),
                "max_ambiguous_anchor_fraction": max(float(row["ambiguous_anchor_fraction"]) for row in records),
                "normalization_changed_raw_labels": str(norm_count < raw_count).lower(),
                "grouping_depends_on_tolerance": "false",
            }
        )

    return normalized_rows, definition_rows, pairwise_rows, partition_rows, validation_rows


def overlapping_raw_labels(
    species: str,
    interval: dict[str, str],
    raw_states: list[dict[str, str]],
    intervals: list[dict[str, str]],
) -> str:
    interval_by_id = {row["structural_interval_id"]: row for row in intervals}
    labels = []
    start = as_int(interval, "start")
    end = as_int(interval, "end")
    for state in raw_states:
        if state["species"] != species:
            continue
        raw_interval = interval_by_id.get(state["structural_interval_id"])
        if raw_interval and overlap(start, end, as_int(raw_interval, "start"), as_int(raw_interval, "end")) > 0:
            labels.append(f"{state['structural_interval_id']}:{state['arrangement_state']}")
    return ";".join(labels)


def write_outputs(
    normalized_rows: list[dict[str, object]],
    definition_rows: list[dict[str, object]],
    pairwise_rows: list[dict[str, object]],
    partition_rows: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
) -> None:
    write_tsv(
        NORMALIZED_STATES_PATH,
        normalized_rows,
        [
            "species",
            "structural_interval_id",
            "normalized_arrangement_state",
            "anchor_order_signature",
            "anchor_orientation_signature",
            "n_shared_anchors",
            "missing_anchor_fraction",
            "ambiguous_anchor_fraction",
            "raw_state_label",
            "equivalence_reason",
        ],
    )
    write_tsv(
        NORMALIZED_DEFS_PATH,
        definition_rows,
        [
            "structural_interval_id",
            "normalized_arrangement_state",
            "species_members",
            "n_species",
            "n_shared_anchors",
            "definition",
            "anchor_order_signature",
            "anchor_orientation_signature",
        ],
    )
    write_tsv(
        PAIRWISE_PATH,
        pairwise_rows,
        [
            "structural_interval_id",
            "species1",
            "species2",
            "shared_anchor_count",
            "order_concordance",
            "orientation_concordance",
            "overall_similarity",
            "same_normalized_state",
        ],
    )
    write_tsv(PARTITIONS_PATH, partition_rows, ["structural_interval_id", "state", "species_members", "n_species"])
    write_tsv(
        VALIDATION_PATH,
        validation_rows,
        [
            "structural_interval_id",
            "n_anchors",
            "raw_state_count",
            "normalized_state_count",
            "n_species_in_shared_normalized_states",
            "all_species_unique",
            "max_missing_anchor_fraction",
            "max_ambiguous_anchor_fraction",
            "normalization_changed_raw_labels",
            "grouping_depends_on_tolerance",
        ],
    )


def make_plot(normalized_rows: list[dict[str, object]], intervals: list[dict[str, str]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    state_palette = [
        "#4c78a8",
        "#f58518",
        "#54a24b",
        "#e45756",
        "#72b7b2",
        "#b279a2",
        "#ff9da6",
        "#9d755d",
        "#bab0ac",
        "#59a14f",
        "#edc948",
        "#af7aa1",
    ]
    interval_by_id = {row["structural_interval_id"]: row for row in intervals}
    state_colors = {}
    y_positions = {species: i for i, species in enumerate(PLOT_ORDER)}
    fig, ax = plt.subplots(figsize=(12, 5.5))
    for row in normalized_rows:
        species = str(row["species"])
        interval = interval_by_id[str(row["structural_interval_id"])]
        state_key = f"{row['structural_interval_id']}:{row['normalized_arrangement_state']}"
        if state_key not in state_colors:
            state_colors[state_key] = state_palette[len(state_colors) % len(state_palette)]
        y = y_positions[species]
        x0 = as_int(interval, "start") / 1e6
        x1 = as_int(interval, "end") / 1e6
        ax.barh(y, x1 - x0, left=x0, height=0.56, color=state_colors[state_key], edgecolor="white", linewidth=0.3)
        if x1 - x0 >= 0.25:
            ax.text((x0 + x1) / 2, y, str(row["normalized_arrangement_state"]), ha="center", va="center", fontsize=6)
    ax.set_yticks([y_positions[species] for species in PLOT_ORDER], PLOT_ORDER)
    ax.set_xlabel("Chicken GalGal6 chromosome 4 coordinate (Mb)")
    ax.set_ylabel("")
    ax.set_title("Normalized structural arrangement states", loc="left")
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.6)
    ax.set_ylim(-0.7, len(PLOT_ORDER) - 0.3)
    fig.tight_layout()
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=220)
    plt.close(fig)


def validate_no_stage1_columns() -> None:
    forbidden = {"q1", "q2", "q3", "dominant_topology", "dominant_support"}
    for path in (BLOCKS_PATH, INTERVALS_PATH, RAW_STATES_PATH):
        with path.open() as handle:
            header = next(csv.reader(handle, delimiter="\t"))
        overlap_cols = forbidden & set(header)
        if overlap_cols:
            raise ValueError(f"Stage-2.5 structural input {path} contains forbidden Stage-1 columns: {sorted(overlap_cols)}")


def print_report(
    validation_rows: list[dict[str, object]],
    partition_rows: list[dict[str, object]],
    pairwise_rows: list[dict[str, object]],
) -> None:
    print("Raw vs normalized state counts:")
    for row in validation_rows:
        print(
            f"  {row['structural_interval_id']}: raw={row['raw_state_count']}, "
            f"normalized={row['normalized_state_count']}, anchors={row['n_anchors']}"
        )
    print("Shared normalized states:")
    for row in partition_rows:
        if int(row["n_species"]) > 1:
            print(f"  {row['structural_interval_id']} {row['state']}: {row['species_members']}")
    low = [row for row in pairwise_rows if float(row["overall_similarity"]) < LOW_SIMILARITY_THRESHOLD]
    print(f"Low pairwise structural similarities (<{LOW_SIMILARITY_THRESHOLD}): {len(low)}")
    print("Files written:")
    for path in (NORMALIZED_STATES_PATH, NORMALIZED_DEFS_PATH, PAIRWISE_PATH, PARTITIONS_PATH, VALIDATION_PATH, PDF_PATH, PNG_PATH):
        print(f"  {path.relative_to(REPO_ROOT)}")


class NormalizeTests(unittest.TestCase):
    def test_fragmentation_invariance(self) -> None:
        anchors = [{"anchor_index": 1, "start": 0, "end": 1000}]
        whole = [
            {"reference_start": "0", "reference_end": "1000", "query_chrom": "chrA", "query_midpoint": "500", "orientation": "+"}
        ]
        fragmented = [
            {"reference_start": "0", "reference_end": "600", "query_chrom": "chrA", "query_midpoint": "300", "orientation": "+"},
            {"reference_start": "600", "reference_end": "1000", "query_chrom": "chrA", "query_midpoint": "800", "orientation": "+"},
        ]
        self.assertEqual(
            project_species("sp1", anchors, whole)["compatibility_tokens"],
            project_species("sp2", anchors, fragmented)["compatibility_tokens"],
        )

    def test_reversed_orientation_not_merged(self) -> None:
        first = {"compatibility_tokens": {1: ("Q1:1", "+")}}
        second = {"compatibility_tokens": {1: ("Q1:1", "-")}}
        self.assertFalse(compatible(first, second))

    def test_missing_anchor_compatible(self) -> None:
        first = {"compatibility_tokens": {1: ("Q1:1", "+"), 2: ("Q1:2", "+")}}
        second = {"compatibility_tokens": {1: ("Q1:1", "+")}}
        self.assertTrue(compatible(first, second))


def run_tests() -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(NormalizeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run normalization tests first")
    args = parser.parse_args()
    if args.run_tests:
        run_tests()
    validate_no_stage1_columns()
    normalized_rows, definition_rows, pairwise_rows, partition_rows, validation_rows = normalize_states()
    write_outputs(normalized_rows, definition_rows, pairwise_rows, partition_rows, validation_rows)
    make_plot(normalized_rows, read_tsv(INTERVALS_PATH))
    print_report(validation_rows, partition_rows, pairwise_rows)


if __name__ == "__main__":
    main()
