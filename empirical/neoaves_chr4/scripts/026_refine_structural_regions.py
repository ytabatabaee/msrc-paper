#!/usr/bin/env python3
"""Refine chr4 structural intervals into canonical multi-anchor regions."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
import unittest
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"

INTERVALS_PATH = RESULTS_DIR / "structural_intervals.tsv"
BLOCKS_PATH = RESULTS_DIR / "structural_blocks.tsv"
BREAKPOINTS_PATH = RESULTS_DIR / "structural_breakpoints.tsv"

REGIONS_PATH = RESULTS_DIR / "canonical_structural_regions.tsv"
ANCHORS_PATH = RESULTS_DIR / "canonical_anchor_states.tsv"
SIMILARITY_PATH = RESULTS_DIR / "canonical_region_pairwise_similarity.tsv"
STATES_PATH = RESULTS_DIR / "canonical_arrangement_states.tsv"
STATE_DEFS_PATH = RESULTS_DIR / "canonical_arrangement_state_definitions.tsv"
SUMMARY_PATH = RESULTS_DIR / "canonical_region_summary.tsv"
SENSITIVITY_PATH = RESULTS_DIR / "anchor_size_sensitivity.tsv"
VALIDATION_PATH = RESULTS_DIR / "canonical_region_validation.tsv"
PDF_PATH = FIGURES_DIR / "chr4_canonical_structural_regions.pdf"
PNG_PATH = FIGURES_DIR / "chr4_canonical_structural_regions.png"

FOCAL_SPECIES = ["Stork", "Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo"]
PLOT_ORDER = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]
BOUNDARY_TOLERANCE_BP = 200_000
SPAN_RECIPROCAL_OVERLAP = 0.90
SPAN_LENGTH_RATIO = 0.80
MIN_COVERAGE_FRACTION = 0.20
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


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def intervals_cluster(a: dict[str, str], b: dict[str, str], tolerance: int = BOUNDARY_TOLERANCE_BP) -> bool:
    a_start, a_end = as_int(a, "start"), as_int(a, "end")
    b_start, b_end = as_int(b, "start"), as_int(b, "end")
    if abs(a_start - b_start) <= tolerance and abs(a_end - b_end) <= tolerance:
        return True
    ov = overlap_len(a_start, a_end, b_start, b_end)
    if ov == 0:
        return False
    a_len = a_end - a_start
    b_len = b_end - b_start
    reciprocal = min(ov / a_len, ov / b_len)
    length_ratio = min(a_len, b_len) / max(a_len, b_len)
    return reciprocal >= SPAN_RECIPROCAL_OVERLAP and length_ratio >= SPAN_LENGTH_RATIO


def canonicalize_regions(intervals: list[dict[str, str]]) -> list[dict[str, object]]:
    adjacency = defaultdict(set)
    for i, j in combinations(range(len(intervals)), 2):
        if intervals_cluster(intervals[i], intervals[j]):
            adjacency[i].add(j)
            adjacency[j].add(i)
    seen = set()
    components = []
    for i in range(len(intervals)):
        if i in seen:
            continue
        queue = deque([i])
        seen.add(i)
        component = []
        while queue:
            idx = queue.popleft()
            component.append(intervals[idx])
            for nxt in adjacency[idx]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        components.append(component)
    components.sort(key=lambda rows: (min(as_int(r, "start") for r in rows), max(as_int(r, "end") for r in rows), len(rows)))

    regions = []
    for idx, rows in enumerate(components, 1):
        start = min(as_int(row, "start") for row in rows)
        end = max(as_int(row, "end") for row in rows)
        species = sorted({row["species"] for row in rows}, key=FOCAL_SPECIES.index)
        event_types = Counter(row["event_type"] for row in rows)
        region_type = "mixed" if len(event_types) > 1 else event_types.most_common(1)[0][0]
        regions.append(
            {
                "canonical_region_id": f"CR_{idx:02d}",
                "reference_chrom": rows[0]["reference_chrom"],
                "canonical_start": start,
                "canonical_end": end,
                "length_bp": end - start,
                "n_source_species": len(species),
                "species_members": ",".join(species),
                "source_interval_ids": ";".join(row["structural_interval_id"] for row in sorted(rows, key=lambda r: (r["species"], r["start"]))),
                "region_type": region_type,
                "notes": (
                    f"Clustered by boundary proximity <= {BOUNDARY_TOLERANCE_BP} bp or "
                    f"near-identical reciprocal span overlap >= {SPAN_RECIPROCAL_OVERLAP}; "
                    "nested/subregion calls are preserved unless spans are near-identical."
                ),
            }
        )
    return regions


def make_anchor_bins(region: dict[str, object], anchor_size: int) -> list[dict[str, object]]:
    start = int(region["canonical_start"])
    end = int(region["canonical_end"])
    anchors = []
    pos = start
    while pos < end:
        nxt = min(pos + anchor_size, end)
        if anchors and nxt - pos < anchor_size * 0.25:
            anchors[-1]["reference_end"] = nxt
            anchors[-1]["length_bp"] = int(anchors[-1]["reference_end"]) - int(anchors[-1]["reference_start"])
        else:
            anchors.append(
                {
                    "canonical_region_id": region["canonical_region_id"],
                    "anchor_id": f"{region['canonical_region_id']}_A{len(anchors) + 1:03d}",
                    "anchor_index": len(anchors) + 1,
                    "reference_start": pos,
                    "reference_end": nxt,
                    "length_bp": nxt - pos,
                }
            )
        pos = nxt
    return anchors


def dominant_anchor_state(species: str, anchor: dict[str, object], blocks: list[dict[str, str]]) -> dict[str, object]:
    start = int(anchor["reference_start"])
    end = int(anchor["reference_end"])
    length = end - start
    hits = []
    for block in blocks:
        hit = overlap_len(start, end, as_int(block, "reference_start"), as_int(block, "reference_end"))
        if hit > 0:
            hits.append((hit, block))
    total = sum(hit for hit, _ in hits)
    coverage = min(1.0, total / length) if length else 0.0
    if not hits or coverage < MIN_COVERAGE_FRACTION:
        return anchor_row(species, anchor, "", "", "", "", coverage, len(hits), False, True)

    by_chrom = defaultdict(int)
    by_orientation = defaultdict(int)
    for hit, block in hits:
        by_chrom[block["query_chrom"]] += hit
        by_orientation[block["orientation"]] += hit
    top_chrom, top_chrom_bp = sorted(by_chrom.items(), key=lambda item: (-item[1], item[0]))[0]
    top_orient, top_orient_bp = sorted(by_orientation.items(), key=lambda item: (-item[1], item[0]))[0]
    chrom_tie = list(by_chrom.values()).count(top_chrom_bp) > 1
    orient_tie = list(by_orientation.values()).count(top_orient_bp) > 1
    dominant_hits = [block for hit, block in hits if block["query_chrom"] == top_chrom]
    q_start = min(as_int(block, "query_start") for block in dominant_hits)
    q_end = max(as_int(block, "query_end") for block in dominant_hits)
    ambiguous = chrom_tie or orient_tie or (top_chrom_bp / total < 0.60) or (top_orient_bp / total < 0.60)
    return anchor_row(species, anchor, top_chrom, q_start, q_end, top_orient, coverage, len(hits), ambiguous, False)


def anchor_row(
    species: str,
    anchor: dict[str, object],
    chrom: object,
    q_start: object,
    q_end: object,
    orientation: object,
    coverage: float,
    n_blocks: int,
    ambiguous: bool,
    missing: bool,
) -> dict[str, object]:
    return {
        "species": species,
        "canonical_region_id": anchor["canonical_region_id"],
        "anchor_id": anchor["anchor_id"],
        "reference_start": anchor["reference_start"],
        "reference_end": anchor["reference_end"],
        "dominant_query_chrom": chrom,
        "dominant_query_start": q_start,
        "dominant_query_end": q_end,
        "orientation": orientation,
        "coverage_fraction": coverage,
        "n_raw_blocks": n_blocks,
        "ambiguous": str(ambiguous).lower(),
        "missing": str(missing).lower(),
    }


def build_anchor_states(regions: list[dict[str, object]], blocks: list[dict[str, str]], anchor_size: int) -> list[dict[str, object]]:
    blocks_by_species = defaultdict(list)
    for block in blocks:
        blocks_by_species[block["species"]].append(block)
    rows = []
    for region in regions:
        for anchor in make_anchor_bins(region, anchor_size):
            for species in FOCAL_SPECIES:
                rows.append(dominant_anchor_state(species, anchor, blocks_by_species[species]))
    return rows


def species_region_signature(species: str, region_id: str, rows: list[dict[str, object]]) -> dict[str, object]:
    subset = sorted(
        [row for row in rows if row["species"] == species and row["canonical_region_id"] == region_id],
        key=lambda row: row["reference_start"],
    )
    informative = [row for row in subset if row["missing"] == "false" and row["ambiguous"] == "false"]
    by_chrom = defaultdict(list)
    for row in informative:
        by_chrom[str(row["dominant_query_chrom"])].append(row)
    ranks = {}
    for chrom, chrom_rows in by_chrom.items():
        ordered = sorted(chrom_rows, key=lambda row: (float(row["dominant_query_start"]) + float(row["dominant_query_end"])) / 2)
        for rank, row in enumerate(ordered, 1):
            ranks[row["anchor_id"]] = rank
    order_tokens = []
    orient_tokens = []
    position_tokens = []
    for row in subset:
        anchor_id = row["anchor_id"]
        if row["missing"] == "true":
            order_tokens.append(f"{anchor_id}:M")
            orient_tokens.append(f"{anchor_id}:?")
            position_tokens.append(f"{anchor_id}:M")
        elif row["ambiguous"] == "true":
            order_tokens.append(f"{anchor_id}:X")
            orient_tokens.append(f"{anchor_id}:!")
            position_tokens.append(f"{anchor_id}:X")
        else:
            token = f"{anchor_id}:{ranks[anchor_id]}"
            order_tokens.append(token)
            orient_tokens.append(f"{anchor_id}:{row['orientation']}")
            position_tokens.append(f"{anchor_id}:{row['dominant_query_chrom']}:{row['dominant_query_start']}-{row['dominant_query_end']}")
    n = len(subset) or 1
    return {
        "canonical_region_id": region_id,
        "species": species,
        "informative_anchor_count": len(informative),
        "anchor_order_signature": "|".join(order_tokens),
        "orientation_signature": "|".join(orient_tokens),
        "mapped_query_position_signature": "|".join(position_tokens),
        "missing_anchor_fraction": sum(row["missing"] == "true" for row in subset) / n,
        "ambiguous_anchor_fraction": sum(row["ambiguous"] == "true" for row in subset) / n,
        "rank_by_anchor": {row["anchor_id"]: ranks[row["anchor_id"]] for row in informative},
        "orientation_by_anchor": {row["anchor_id"]: row["orientation"] for row in informative},
    }


def compare_signatures(a: dict[str, object], b: dict[str, object], min_shared: int) -> dict[str, object]:
    shared = sorted(set(a["rank_by_anchor"]) & set(b["rank_by_anchor"]))
    n = len(shared)
    missing_fraction = max(float(a["missing_anchor_fraction"]), float(b["missing_anchor_fraction"]))
    if n < min_shared:
        return {
            "informative_shared_anchors": n,
            "order_concordance": "",
            "orientation_concordance": "",
            "missing_fraction": missing_fraction,
            "overall_similarity": "",
            "equivalent_arrangement": "false",
            "equivalence_reason": "insufficient_information",
        }
    orientation_matches = sum(a["orientation_by_anchor"][anchor] == b["orientation_by_anchor"][anchor] for anchor in shared)
    orientation_concordance = orientation_matches / n
    pair_total = 0
    pair_matches = 0
    for left, right in combinations(shared, 2):
        a_delta = sign(int(a["rank_by_anchor"][right]) - int(a["rank_by_anchor"][left]))
        b_delta = sign(int(b["rank_by_anchor"][right]) - int(b["rank_by_anchor"][left]))
        pair_total += 1
        pair_matches += a_delta == b_delta
    order_concordance = pair_matches / pair_total if pair_total else 1.0
    overall = (order_concordance + orientation_concordance) / 2
    equivalent = order_concordance == 1.0 and orientation_concordance == 1.0
    return {
        "informative_shared_anchors": n,
        "order_concordance": order_concordance,
        "orientation_concordance": orientation_concordance,
        "missing_fraction": missing_fraction,
        "overall_similarity": overall,
        "equivalent_arrangement": str(equivalent).lower(),
        "equivalence_reason": "exact_multi_anchor_order_orientation_match" if equivalent else "discordant_order_or_orientation",
    }


def sign(value: int) -> int:
    return (value > 0) - (value < 0)


def assign_arrangement_states(signatures: list[dict[str, object]], pairwise: list[dict[str, object]]) -> None:
    by_region = defaultdict(list)
    equivalent_pairs = defaultdict(set)
    for row in pairwise:
        if row["equivalent_arrangement"] == "true":
            key = row["canonical_region_id"]
            equivalent_pairs[key].add(frozenset((row["species1"], row["species2"])))
    for sig in signatures:
        by_region[sig["canonical_region_id"]].append(sig)
    for region_id, rows in by_region.items():
        components: list[list[str]] = []
        for species in FOCAL_SPECIES:
            placed = False
            for component in components:
                if all(frozenset((species, member)) in equivalent_pairs[region_id] for member in component):
                    component.append(species)
                    placed = True
                    break
            if not placed:
                components.append([species])
        state_for_species = {}
        for i, component in enumerate(components):
            for species in component:
                state_for_species[species] = f"A{i}"
        for sig in rows:
            sig["normalized_arrangement_state"] = state_for_species[sig["species"]]
            sig["confidence"] = "confident" if len([s for s, state in state_for_species.items() if state == state_for_species[sig["species"]]]) > 1 else "singleton_or_unresolved"


def run_refinement(anchor_size: int, min_shared_anchors: int) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    validate_inputs()
    intervals = read_tsv(INTERVALS_PATH)
    blocks = read_tsv(BLOCKS_PATH)
    regions = canonicalize_regions(intervals)
    anchor_rows = build_anchor_states(regions, blocks, anchor_size)
    signatures = []
    for region in regions:
        for species in FOCAL_SPECIES:
            signatures.append(species_region_signature(species, str(region["canonical_region_id"]), anchor_rows))
    pairwise = []
    sig_by_key = {(sig["canonical_region_id"], sig["species"]): sig for sig in signatures}
    for region in regions:
        region_id = str(region["canonical_region_id"])
        for species1, species2 in combinations(FOCAL_SPECIES, 2):
            comparison = compare_signatures(sig_by_key[(region_id, species1)], sig_by_key[(region_id, species2)], min_shared_anchors)
            pairwise.append({"canonical_region_id": region_id, "species1": species1, "species2": species2, **comparison})
    assign_arrangement_states(signatures, pairwise)
    defs, summary = summarize_regions(regions, signatures, anchor_rows, min_shared_anchors)
    return regions, anchor_rows, pairwise, signatures, defs, summary


def summarize_regions(
    regions: list[dict[str, object]],
    signatures: list[dict[str, object]],
    anchor_rows: list[dict[str, object]],
    min_shared_anchors: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    defs = []
    summary = []
    sig_by_region = defaultdict(list)
    anchors_by_region = defaultdict(set)
    for sig in signatures:
        sig_by_region[sig["canonical_region_id"]].append(sig)
    for row in anchor_rows:
        anchors_by_region[row["canonical_region_id"]].add(row["anchor_id"])
    for region in regions:
        region_id = str(region["canonical_region_id"])
        rows = sig_by_region[region_id]
        by_state = defaultdict(list)
        for row in rows:
            by_state[row["normalized_arrangement_state"]].append(row)
        state_parts = []
        for state, members in sorted(by_state.items()):
            species = [m["species"] for m in members]
            state_parts.append(f"{state}:{','.join(species)}")
            example = members[0]
            defs.append(
                {
                    "canonical_region_id": region_id,
                    "normalized_arrangement_state": state,
                    "species_members": ",".join(species),
                "n_species": len(species),
                    "representative_informative_anchor_count": example["informative_anchor_count"],
                    "representative_anchor_order_signature": example["anchor_order_signature"],
                    "representative_orientation_signature": example["orientation_signature"],
                    "description": f"{len(species)} species with equivalent multi-anchor order/orientation in {region_id}.",
                }
            )
        informative_counts = [int(row["informative_anchor_count"]) for row in rows]
        missingness = sum(float(row["missing_anchor_fraction"]) for row in rows) / len(rows)
        ambiguity = sum(float(row["ambiguous_anchor_fraction"]) for row in rows) / len(rows)
        n_states = len(by_state)
        multi_species_states = [members for members in by_state.values() if len(members) > 1]
        if min(informative_counts) < min_shared_anchors:
            complexity = "insufficient_information"
        elif n_states == 2 and len(multi_species_states) == 2 and missingness < 0.25 and ambiguity < 0.25:
            complexity = "clean_two_state"
        elif n_states > 2 and multi_species_states:
            complexity = "multistate_with_shared_state"
        else:
            complexity = "structurally_complex"
        summary.append(
            {
                "canonical_region_id": region_id,
                "start": region["canonical_start"],
                "end": region["canonical_end"],
                "length_bp": region["length_bp"],
                "n_source_species": region["n_source_species"],
                "n_species_evaluated": len(FOCAL_SPECIES),
                "n_normalized_states": n_states,
                "state_partition": ";".join(state_parts),
                "mean_informative_anchor_count": sum(informative_counts) / len(informative_counts),
                "min_informative_anchor_count": min(informative_counts),
                "missingness": missingness,
                "ambiguity": ambiguity,
                "n_anchors": len(anchors_by_region[region_id]),
                "complexity_class": complexity,
                "usable_for_simple_two_state_test": str(complexity == "clean_two_state").lower(),
            }
        )
    return defs, summary


def validate_inputs() -> None:
    forbidden = {"q1", "q2", "q3", "dominant_topology", "dominant_support", "topology"}
    for path in (INTERVALS_PATH, BLOCKS_PATH, BREAKPOINTS_PATH):
        with path.open() as handle:
            header = next(csv.reader(handle, delimiter="\t"))
        bad = forbidden & set(header)
        if bad:
            raise ValueError(f"Forbidden genealogy columns in structural input {path}: {sorted(bad)}")


def write_primary_outputs(
    regions: list[dict[str, object]],
    anchors: list[dict[str, object]],
    pairwise: list[dict[str, object]],
    signatures: list[dict[str, object]],
    defs: list[dict[str, object]],
    summary: list[dict[str, object]],
    validation: list[dict[str, object]],
) -> None:
    write_tsv(REGIONS_PATH, regions, ["canonical_region_id", "reference_chrom", "canonical_start", "canonical_end", "length_bp", "n_source_species", "species_members", "source_interval_ids", "region_type", "notes"])
    write_tsv(ANCHORS_PATH, anchors, ["species", "canonical_region_id", "anchor_id", "reference_start", "reference_end", "dominant_query_chrom", "dominant_query_start", "dominant_query_end", "orientation", "coverage_fraction", "n_raw_blocks", "ambiguous", "missing"])
    write_tsv(SIMILARITY_PATH, pairwise, ["canonical_region_id", "species1", "species2", "informative_shared_anchors", "order_concordance", "orientation_concordance", "missing_fraction", "overall_similarity", "equivalent_arrangement", "equivalence_reason"])
    write_tsv(STATES_PATH, signatures, ["canonical_region_id", "species", "normalized_arrangement_state", "informative_anchor_count", "anchor_order_signature", "orientation_signature", "mapped_query_position_signature", "missing_anchor_fraction", "ambiguous_anchor_fraction", "confidence"])
    write_tsv(STATE_DEFS_PATH, defs, ["canonical_region_id", "normalized_arrangement_state", "species_members", "n_species", "representative_informative_anchor_count", "description", "representative_anchor_order_signature", "representative_orientation_signature"])
    write_tsv(SUMMARY_PATH, summary, ["canonical_region_id", "start", "end", "length_bp", "n_source_species", "n_species_evaluated", "n_normalized_states", "state_partition", "mean_informative_anchor_count", "min_informative_anchor_count", "missingness", "ambiguity", "n_anchors", "complexity_class", "usable_for_simple_two_state_test"])
    write_tsv(VALIDATION_PATH, validation, ["metric", "value"])


def sensitivity(anchor_sizes: list[int], min_shared_anchors: int) -> list[dict[str, object]]:
    rows_by_size = {}
    for size in anchor_sizes:
        regions, anchors, pairwise, signatures, defs, summary = run_refinement(size, min_shared_anchors)
        partition = {}
        for row in summary:
            partition[row["canonical_region_id"]] = row["state_partition"]
        rows_by_size[size] = (summary, partition)
    all_regions = sorted({region for _, partition in rows_by_size.values() for region in partition})
    out = []
    baseline = rows_by_size[100000][1] if 100000 in rows_by_size else rows_by_size[anchor_sizes[0]][1]
    for region_id in all_regions:
        partitions = {size: rows_by_size[size][1].get(region_id, "") for size in anchor_sizes}
        stable = len(set(partitions.values())) == 1
        out.append(
            {
                "canonical_region_id": region_id,
                "partition_50000": partitions.get(50000, ""),
                "partition_100000": partitions.get(100000, ""),
                "partition_200000": partitions.get(200000, ""),
                "stable_across_anchor_sizes": str(stable).lower(),
                "matches_100kb_baseline": str(all(value == baseline.get(region_id, "") for value in partitions.values())).lower(),
            }
        )
    return out


def validation_rows(intervals: list[dict[str, str]], regions: list[dict[str, object]], anchors: list[dict[str, object]], summary: list[dict[str, object]], sensitivity_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    anchors_by_region = Counter(row["canonical_region_id"] for row in anchors if row["species"] == FOCAL_SPECIES[0])
    return [
        {"metric": "species_specific_intervals_before_canonicalization", "value": len(intervals)},
        {"metric": "canonical_regions_after_clustering", "value": len(regions)},
        {"metric": "regions_with_one_anchor", "value": sum(count == 1 for count in anchors_by_region.values())},
        {"metric": "regions_with_one_min_informative_anchor", "value": sum(float(row["min_informative_anchor_count"]) <= 1 for row in summary)},
        {"metric": "clean_two_state_candidate_regions", "value": ",".join(row["canonical_region_id"] for row in summary if row["usable_for_simple_two_state_test"] == "true")},
        {"metric": "insufficient_information_regions", "value": ",".join(row["canonical_region_id"] for row in summary if row["complexity_class"] == "insufficient_information")},
        {"metric": "stable_regions_50_100_200kb", "value": ",".join(row["canonical_region_id"] for row in sensitivity_rows if row["stable_across_anchor_sizes"] == "true")},
        {"metric": "low_pairwise_similarity_threshold", "value": LOW_SIMILARITY_THRESHOLD},
    ]


def make_plot(regions: list[dict[str, object]], signatures: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = ["#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2", "#b279a2", "#ff9da6", "#9d755d"]
    colors = {}
    region_by_id = {row["canonical_region_id"]: row for row in regions}
    y_positions = {species: i for i, species in enumerate(PLOT_ORDER)}
    fig, ax = plt.subplots(figsize=(12, 5.5))
    for sig in signatures:
        region = region_by_id[sig["canonical_region_id"]]
        key = f"{sig['canonical_region_id']}:{sig['normalized_arrangement_state']}"
        if key not in colors:
            colors[key] = palette[len(colors) % len(palette)]
        y = y_positions[sig["species"]]
        x0 = int(region["canonical_start"]) / 1e6
        x1 = int(region["canonical_end"]) / 1e6
        ax.barh(y, x1 - x0, left=x0, height=0.48, color=colors[key], edgecolor="white", linewidth=0.3, alpha=0.9)
        if x1 - x0 > 0.3:
            ax.text((x0 + x1) / 2, y, sig["normalized_arrangement_state"], ha="center", va="center", fontsize=6)
    for region in regions:
        ax.axvline(int(region["canonical_start"]) / 1e6, color="#333333", linewidth=0.4, alpha=0.45)
        ax.axvline(int(region["canonical_end"]) / 1e6, color="#333333", linewidth=0.4, alpha=0.45)
    ax.set_yticks([y_positions[species] for species in PLOT_ORDER], PLOT_ORDER)
    ax.set_xlabel("Chicken GalGal6 chromosome 4 coordinate (Mb)")
    ax.set_title("Canonical structural regions and normalized arrangement states", loc="left")
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.6)
    ax.set_ylim(-0.7, len(PLOT_ORDER) - 0.3)
    fig.tight_layout()
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=220)
    plt.close(fig)


def print_report(regions, summary, sensitivity_rows, validation) -> None:
    print(f"Species-specific intervals -> canonical regions: {validation[0]['value']} -> {validation[1]['value']}")
    print("Anchors and partitions:")
    for row in summary:
        print(f"  {row['canonical_region_id']} {row['start']}-{row['end']}: anchors={row['n_anchors']}, states={row['n_normalized_states']}, {row['complexity_class']}, {row['state_partition']}")
    stable = [row["canonical_region_id"] for row in sensitivity_rows if row["stable_across_anchor_sizes"] == "true"]
    print(f"Stable across 50/100/200 kb: {','.join(stable) if stable else 'none'}")
    clean = [row["canonical_region_id"] for row in summary if row["usable_for_simple_two_state_test"] == "true"]
    print(f"Clean two-state candidates: {','.join(clean) if clean else 'none'}")
    print("Files written:")
    for path in (REGIONS_PATH, ANCHORS_PATH, SIMILARITY_PATH, STATES_PATH, STATE_DEFS_PATH, SUMMARY_PATH, SENSITIVITY_PATH, VALIDATION_PATH, PDF_PATH, PNG_PATH):
        print(f"  {path.relative_to(REPO_ROOT)}")


class RefinementTests(unittest.TestCase):
    def test_long_region_many_anchors(self) -> None:
        anchors = make_anchor_bins({"canonical_region_id": "CR_X", "canonical_start": 0, "canonical_end": 7_500_000}, 100_000)
        self.assertGreaterEqual(len(anchors), 70)

    def test_fragmentation_same_state(self) -> None:
        anchor = {"canonical_region_id": "CR_X", "anchor_id": "A1", "reference_start": 0, "reference_end": 100_000}
        whole = [{"species": "A", "reference_start": "0", "reference_end": "100000", "query_chrom": "c", "query_start": "0", "query_end": "100000", "orientation": "+"}]
        frag = [
            {"species": "B", "reference_start": "0", "reference_end": "50000", "query_chrom": "c", "query_start": "0", "query_end": "50000", "orientation": "+"},
            {"species": "B", "reference_start": "50000", "reference_end": "100000", "query_chrom": "c", "query_start": "50000", "query_end": "100000", "orientation": "+"},
        ]
        self.assertEqual(dominant_anchor_state("A", anchor, whole)["orientation"], dominant_anchor_state("B", anchor, frag)["orientation"])

    def test_inversion_differs(self) -> None:
        a = {"rank_by_anchor": {"a": 1, "b": 2, "c": 3, "d": 4, "e": 5}, "orientation_by_anchor": {x: "+" for x in "abcde"}, "missing_anchor_fraction": 0}
        b = {"rank_by_anchor": {"a": 5, "b": 4, "c": 3, "d": 2, "e": 1}, "orientation_by_anchor": {x: "-" for x in "abcde"}, "missing_anchor_fraction": 0}
        self.assertEqual(compare_signatures(a, b, 5)["equivalent_arrangement"], "false")

    def test_missing_not_automatically_distinct(self) -> None:
        a = {"rank_by_anchor": {"a": 1, "b": 2, "c": 3, "d": 4, "e": 5}, "orientation_by_anchor": {x: "+" for x in "abcde"}, "missing_anchor_fraction": 0}
        b = {"rank_by_anchor": {"a": 1, "b": 2, "c": 3, "d": 4, "e": 5}, "orientation_by_anchor": {x: "+" for x in "abcde"}, "missing_anchor_fraction": 0.5}
        self.assertEqual(compare_signatures(a, b, 5)["equivalent_arrangement"], "true")

    def test_canonical_groups_near_identical(self) -> None:
        a = {"start": "1000", "end": "5000"}
        b = {"start": "1100", "end": "5100"}
        self.assertTrue(intervals_cluster(a, b, tolerance=200))

    def test_nested_not_collapsed(self) -> None:
        a = {"start": "0", "end": "1000000"}
        b = {"start": "400000", "end": "600000"}
        self.assertFalse(intervals_cluster(a, b, tolerance=1000))

    def test_validate_no_genealogy(self) -> None:
        validate_inputs()

    def test_state_groups_require_all_pairs(self) -> None:
        signatures = [
            {"canonical_region_id": "CR_X", "species": "Stork"},
            {"canonical_region_id": "CR_X", "species": "Flamingo"},
            {"canonical_region_id": "CR_X", "species": "Dove"},
            {"canonical_region_id": "CR_X", "species": "Sandgrouse"},
            {"canonical_region_id": "CR_X", "species": "Turaco"},
            {"canonical_region_id": "CR_X", "species": "Cuckoo"},
        ]
        pairwise = [
            {
                "canonical_region_id": "CR_X",
                "species1": "Stork",
                "species2": "Flamingo",
                "equivalent_arrangement": "true",
            },
            {
                "canonical_region_id": "CR_X",
                "species1": "Flamingo",
                "species2": "Dove",
                "equivalent_arrangement": "true",
            },
        ]
        assign_arrangement_states(signatures, pairwise)
        state_by_species = {row["species"]: row["normalized_arrangement_state"] for row in signatures}
        self.assertEqual(state_by_species["Stork"], state_by_species["Flamingo"])
        self.assertNotEqual(state_by_species["Stork"], state_by_species["Dove"])


def run_tests() -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RefinementTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-size", type=int, default=100_000)
    parser.add_argument("--min-shared-anchors", type=int, default=5)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        run_tests()

    regions, anchors, pairwise, signatures, defs, summary = run_refinement(args.anchor_size, args.min_shared_anchors)
    sensitivity_rows = sensitivity([50_000, 100_000, 200_000], args.min_shared_anchors)
    validation = validation_rows(read_tsv(INTERVALS_PATH), regions, anchors, summary, sensitivity_rows)
    write_primary_outputs(regions, anchors, pairwise, signatures, defs, summary, validation)
    write_tsv(SENSITIVITY_PATH, sensitivity_rows, ["canonical_region_id", "partition_50000", "partition_100000", "partition_200000", "stable_across_anchor_sizes", "matches_100kb_baseline"])
    make_plot(regions, signatures)
    print_report(regions, summary, sensitivity_rows, validation)


if __name__ == "__main__":
    main()
