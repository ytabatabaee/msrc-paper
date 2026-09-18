#!/usr/bin/env python3
"""Build structural-only chromosome-4 synteny states for MSRC analysis."""

from __future__ import annotations

import argparse
import csv
import lzma
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
REARRANGEMENTS_DIR = REPO_ROOT / "rearrangements"
RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"
FIGURES_DIR = Path(__file__).resolve().parents[1] / "figures"

TOTAB_PATH = REARRANGEMENTS_DIR / "totab_chr4.csv.xz"
PUBLISHED_BLOCKS_PATH = REARRANGEMENTS_DIR / "outlier-regions-by-maf2synteny.tsv"
PUBLISHED_SUMMARY_PATH = REARRANGEMENTS_DIR / "outlier-regions-by-maf2synteny-summary.tsv"

BLOCKS_PATH = RESULTS_DIR / "structural_blocks.tsv"
INTERVALS_PATH = RESULTS_DIR / "structural_intervals.tsv"
BREAKPOINTS_PATH = RESULTS_DIR / "structural_breakpoints.tsv"
STATES_PATH = RESULTS_DIR / "arrangement_states.tsv"
STATE_DEFS_PATH = RESULTS_DIR / "arrangement_state_definitions.tsv"
VALIDATION_PATH = RESULTS_DIR / "structural_interval_validation.tsv"
SUMMARY_PATH = RESULTS_DIR / "structural_validation_summary.tsv"
METADATA_PATH = RESULTS_DIR / "structural_source_metadata.tsv"
PDF_PATH = FIGURES_DIR / "chr4_structural_track.pdf"
PNG_PATH = FIGURES_DIR / "chr4_structural_track.png"

REFERENCE_DESCRIPTION = "GCF_000002315.5_GalGal6.chr4"
REFERENCE_SPECIES = "Chicken"
REFERENCE_ASSEMBLY = "GCF_000002315.5_GalGal6"
STRUCTURAL_ANALYSIS_REFERENCE = "GCA_017639555.1_bCicMag1.pri CM030196.1"
REFERENCE_CHROM = "chr4"

SPECIES_LABELS = {
    "GCA_017639555": "Stork",
    "GCA_009819775": "Flamingo",
    "GCA_901699155": "Dove",
    "GCA_009769525": "Sandgrouse",
    "GCA_009769465": "Turaco",
    "GCA_017976375": "Cuckoo",
    "GCF_003957565": "Finch",
}
FOCAL_SPECIES = ["Stork", "Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo"]
PLOT_ORDER = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]
MERGE_GAP_BP = 1_000_000
APPROX_MATCH_BP = 1_000_000


def parse_description(description: str) -> dict[str, str]:
    parts = description.split(".", 2)
    return {
        "assembly": parts[0] if parts else "",
        "genome": parts[1] if len(parts) > 1 else "",
        "chrom": parts[2] if len(parts) > 2 else "",
    }


def clean_int(value: str) -> int:
    return int(value.replace(",", "").replace('"', "").strip())


def normalize_interval(start: int, end: int) -> tuple[int, int]:
    return (start, end) if start <= end else (end, start)


def load_paired_blocks() -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with lzma.open(TOTAB_PATH, "rt") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed = parse_description(row["Description"])
            row.update(
                {
                    "assembly": parsed["assembly"],
                    "genome": parsed["genome"],
                    "chrom": parsed["chrom"],
                }
            )
            by_pair[(row["sp"], row["Block"])].append(row)

    blocks = []
    for (sp, block_id), rows in by_pair.items():
        ref_rows = [
            row
            for row in rows
            if row["Description"] == REFERENCE_DESCRIPTION and row["chrom"] == REFERENCE_CHROM
        ]
        query_rows = [row for row in rows if row["Description"] != REFERENCE_DESCRIPTION]
        if len(ref_rows) != 1 or len(query_rows) != 1:
            continue
        if sp not in SPECIES_LABELS or SPECIES_LABELS[sp] not in FOCAL_SPECIES:
            continue

        ref = ref_rows[0]
        query = query_rows[0]
        ref_start, ref_end = normalize_interval(clean_int(ref["Start"]), clean_int(ref["End"]))
        query_start, query_end = normalize_interval(clean_int(query["Start"]), clean_int(query["End"]))
        if ref["chrom"] != REFERENCE_CHROM:
            continue

        blocks.append(
            {
                "species": SPECIES_LABELS[sp],
                "genome_accession": sp,
                "reference_species": REFERENCE_SPECIES,
                "reference_assembly": REFERENCE_ASSEMBLY,
                "reference_chrom": ref["chrom"],
                "reference_start": ref_start,
                "reference_end": ref_end,
                "reference_midpoint": (ref_start + ref_end) / 2,
                "query_assembly": query["assembly"],
                "query_genome": query["genome"],
                "query_chrom": query["chrom"],
                "query_start": query_start,
                "query_end": query_end,
                "query_midpoint": (query_start + query_end) / 2,
                "block_id": block_id,
                "orientation": query["Strand"],
                "reference_orientation": ref["Strand"],
                "query_original_start": clean_int(query["Start"]),
                "query_original_end": clean_int(query["End"]),
                "reference_original_start": clean_int(ref["Start"]),
                "reference_original_end": clean_int(ref["End"]),
                "block_length": clean_int(query["Length"]),
                "reference_block_length": clean_int(ref["Length"]),
                "query_seq_id": query["Seq_id"],
                "reference_seq_id": ref["Seq_id"],
                "query_seq_size": clean_int(query["Size"]),
                "source_file": str(TOTAB_PATH.relative_to(REPO_ROOT)),
                "structural_interval_id": "",
            }
        )

    blocks.sort(key=lambda row: (row["species"], row["reference_start"], row["reference_end"]))
    for species, subset in group_by(blocks, "species").items():
        for i, row in enumerate(sorted(subset, key=lambda r: (r["reference_start"], r["reference_end"])), 1):
            row["order_index"] = i
        chrom_groups = defaultdict(list)
        for row in subset:
            chrom_groups[row["query_chrom"]].append(row)
        for chrom_subset in chrom_groups.values():
            for i, row in enumerate(sorted(chrom_subset, key=lambda r: r["query_midpoint"]), 1):
                row["query_order_index"] = i
    return blocks


def group_by(rows: list[dict[str, object]], key: str) -> dict[object, list[dict[str, object]]]:
    grouped: dict[object, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    return grouped


def read_published_outlier_blocks() -> list[dict[str, object]]:
    rows = []
    with PUBLISHED_BLOCKS_PATH.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        has_row_names = len(header) == 15
        for fields in reader:
            if not fields:
                continue
            if has_row_names and len(fields) == 16:
                fields = fields[1:]
            if len(fields) != 15:
                raise ValueError(f"Unexpected field count in {PUBLISHED_BLOCKS_PATH}: {fields!r}")
            species, ident, chrom, start, end, length, block, strand, description, sp, seq_id, size, tag, gal_start, gal_end = fields
            if species == "Chicken" or species.startswith("Finch"):
                continue
            species = species.split()[0]
            if species not in FOCAL_SPECIES:
                continue
            ref_start, ref_end = normalize_interval(clean_int(gal_start), clean_int(gal_end))
            query_start, query_end = normalize_interval(clean_int(start), clean_int(end))
            rows.append(
                {
                    "species": species,
                    "genome_accession": sp,
                    "query_assembly": ident,
                    "query_chrom": chrom,
                    "query_start": query_start,
                    "query_end": query_end,
                    "block_id": block,
                    "orientation": strand,
                    "reference_chrom": tag,
                    "reference_start": ref_start,
                    "reference_end": ref_end,
                    "source_file": str(PUBLISHED_BLOCKS_PATH.relative_to(REPO_ROOT)),
                    "description": description,
                }
            )
    return rows


def infer_intervals(published_blocks: list[dict[str, object]]) -> list[dict[str, object]]:
    intervals = []
    for species in FOCAL_SPECIES:
        groups = []
        species_rows = sorted(
            [row for row in published_blocks if row["species"] == species],
            key=lambda row: (row["reference_start"], row["reference_end"]),
        )
        reference_clusters = []
        if species_rows:
            current_ref_cluster = [species_rows[0]]
            current_ref_end = int(species_rows[0]["reference_end"])
            for row in species_rows[1:]:
                if int(row["reference_start"]) - current_ref_end <= MERGE_GAP_BP:
                    current_ref_cluster.append(row)
                    current_ref_end = max(current_ref_end, int(row["reference_end"]))
                else:
                    reference_clusters.append(current_ref_cluster)
                    current_ref_cluster = [row]
                    current_ref_end = int(row["reference_end"])
            reference_clusters.append(current_ref_cluster)

        for reference_cluster in reference_clusters:
            by_query_chrom = defaultdict(list)
            for row in reference_cluster:
                by_query_chrom[row["query_chrom"]].append(row)
            for query_chrom, subset in by_query_chrom.items():
                ordered = sorted(subset, key=lambda row: (row["query_start"], row["query_end"]))
                if not ordered:
                    continue
                current = [ordered[0]]
                current_end = int(ordered[0]["query_end"])
                for row in ordered[1:]:
                    if int(row["query_start"]) - current_end <= MERGE_GAP_BP:
                        current.append(row)
                        current_end = max(current_end, int(row["query_end"]))
                    else:
                        groups.append(current)
                        current = [row]
                        current_end = int(row["query_end"])
                groups.append(current)

        groups.sort(
            key=lambda group: (
                min(int(row["reference_start"]) for row in group),
                max(int(row["reference_end"]) for row in group),
                str(group[0]["query_chrom"]),
            )
        )
        for index, group in enumerate(groups, 1):
            ref_start = min(int(row["reference_start"]) for row in group)
            ref_end = max(int(row["reference_end"]) for row in group)
            chrom_counts = Counter(str(row["query_chrom"]) for row in group)
            orient_counts = Counter(str(row["orientation"]) for row in group)
            query_chrom = str(group[0]["query_chrom"])
            q_starts = [int(row["query_start"]) for row in group]
            q_ends = [int(row["query_end"]) for row in group]
            query_start = min(q_starts)
            query_end = max(q_ends)
            interval_id = f"SI_{species}_{index:02d}"
            event_type = classify_event(group)
            intervals.append(
                {
                    "species": species,
                    "structural_interval_id": interval_id,
                    "reference_chrom": REFERENCE_CHROM,
                    "start": ref_start,
                    "end": ref_end,
                    "length_bp": ref_end - ref_start,
                    "event_type": event_type,
                    "n_blocks": len({row["block_id"] for row in group}),
                    "evidence": (
                        f"published maf2synteny outlier blocks; "
                        f"query_chroms={format_counter(chrom_counts)}; "
                        f"orientations={format_counter(orient_counts)}"
                    ),
                    "source": str(PUBLISHED_BLOCKS_PATH.relative_to(REPO_ROOT)),
                    "confidence_or_status": "published_structural_region_boundary_not_nucleotide_resolved",
                    "query_chrom": query_chrom,
                    "query_start": query_start,
                    "query_end": query_end,
                    "query_chrom_count": len(chrom_counts),
                    "structural_analysis_reference": STRUCTURAL_ANALYSIS_REFERENCE,
                    "main_coordinate_reference": REFERENCE_DESCRIPTION,
                }
            )
    return intervals


def classify_event(group: list[dict[str, object]]) -> str:
    chrom_count = len({row["query_chrom"] for row in group})
    orientations = {row["orientation"] for row in group}
    if chrom_count > 1:
        return "complex"
    ordered = sorted(group, key=lambda row: (row["reference_start"], row["reference_end"]))
    query_midpoints = [(int(row["query_start"]) + int(row["query_end"])) / 2 for row in ordered]
    if len(query_midpoints) < 2:
        return "unknown"
    increasing = sum(b >= a for a, b in zip(query_midpoints, query_midpoints[1:]))
    decreasing = sum(b <= a for a, b in zip(query_midpoints, query_midpoints[1:]))
    total = len(query_midpoints) - 1
    if orientations == {"-"} or decreasing / total >= 0.8:
        return "inversion"
    if increasing / total < 0.8:
        return "reordered_block"
    return "unknown"


def format_counter(counter: Counter) -> str:
    return ",".join(f"{key}:{value}" for key, value in sorted(counter.items()))


def assign_interval_ids_to_blocks(
    blocks: list[dict[str, object]], intervals: list[dict[str, object]]
) -> None:
    intervals_by_species = group_by(intervals, "species")
    for row in blocks:
        matches = [
            interval["structural_interval_id"]
            for interval in intervals_by_species.get(row["species"], [])
            if ranges_overlap(
                int(row["reference_start"]),
                int(row["reference_end"]),
                int(interval["start"]),
                int(interval["end"]),
            )
        ]
        row["structural_interval_id"] = ";".join(matches)


def ranges_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return max(a_start, b_start) <= min(a_end, b_end)


def build_breakpoints(intervals: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for interval in intervals:
        for side, position in (("left_boundary", interval["start"]), ("right_boundary", interval["end"])):
            rows.append(
                {
                    "species": interval["species"],
                    "structural_interval_id": interval["structural_interval_id"],
                    "breakpoint_id": f"{interval['structural_interval_id']}_{side}",
                    "reference_chrom": interval["reference_chrom"],
                    "position": position,
                    "breakpoint_type": "published_or_inferred_structural_boundary_not_nucleotide_resolved",
                    "source": interval["source"],
                }
            )
    return rows


def build_arrangement_states(
    intervals: list[dict[str, object]], blocks: list[dict[str, object]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    blocks_by_species = group_by(blocks, "species")
    signature_to_state: dict[tuple[str, str], str] = {}
    rows = []
    for interval in sorted(intervals, key=lambda r: (r["start"], r["end"], r["species"])):
        subset = [
            row
            for row in blocks_by_species.get(interval["species"], [])
            if ranges_overlap(
                int(row["reference_start"]),
                int(row["reference_end"]),
                int(interval["start"]),
                int(interval["end"]),
            )
        ]
        block_order_signature, orientation_signature = normalized_signature(subset)
        key = (block_order_signature, orientation_signature)
        if key not in signature_to_state:
            signature_to_state[key] = f"A{len(signature_to_state)}"
        rows.append(
            {
                "species": interval["species"],
                "structural_interval_id": interval["structural_interval_id"],
                "arrangement_state": signature_to_state[key],
                "block_order_signature": block_order_signature,
                "orientation_signature": orientation_signature,
            }
        )

    definitions = [
        {
            "arrangement_state": state,
            **describe_signature(key[0], key[1]),
            "block_order_signature": key[0],
            "orientation_signature": key[1],
            "example_intervals": ";".join(
                f"{row['species']}:{row['structural_interval_id']}"
                for row in rows
                if row["arrangement_state"] == state
            ),
            "definition": describe_state(key[0], key[1]),
        }
        for key, state in sorted(signature_to_state.items(), key=lambda item: item[1])
    ]
    return rows, definitions


def normalized_signature(blocks: list[dict[str, object]]) -> tuple[str, str]:
    ordered = sorted(blocks, key=lambda row: (row["reference_start"], row["reference_end"]))
    chrom_alias: dict[str, str] = {}
    chrom_counters = defaultdict(int)
    query_sorted = defaultdict(list)
    for row in ordered:
        query_sorted[row["query_chrom"]].append(row)
    query_rank = {}
    for chrom, rows in query_sorted.items():
        for rank, row in enumerate(sorted(rows, key=lambda r: r["query_midpoint"]), 1):
            query_rank[(row["block_id"], chrom)] = rank
    order_parts = []
    orient_parts = []
    for row in ordered:
        chrom = str(row["query_chrom"])
        if chrom not in chrom_alias:
            chrom_alias[chrom] = f"Q{len(chrom_alias) + 1}"
        alias = chrom_alias[chrom]
        chrom_counters[alias] += 1
        rank = query_rank[(row["block_id"], chrom)]
        order_parts.append(f"{alias}:{rank}")
        orient_parts.append(str(row["orientation"]))
    return "|".join(order_parts), "".join(orient_parts)


def describe_signature(block_order_signature: str, orientation_signature: str) -> dict[str, object]:
    parts = block_order_signature.split("|") if block_order_signature else []
    aliases = sorted({part.split(":", 1)[0] for part in parts})
    ranks = [int(part.split(":", 1)[1]) for part in parts]
    increasing = all(b >= a for a, b in zip(ranks, ranks[1:]))
    decreasing = all(b <= a for a, b in zip(ranks, ranks[1:]))
    if increasing:
        order_pattern = "increasing"
    elif decreasing:
        order_pattern = "decreasing"
    else:
        order_pattern = "mixed"
    orient_counts = Counter(orientation_signature)
    return {
        "n_blocks": len(parts),
        "query_chrom_alias_count": len(aliases),
        "order_pattern": order_pattern,
        "orientation_counts": format_counter(orient_counts),
    }


def describe_state(block_order_signature: str, orientation_signature: str) -> str:
    details = describe_signature(block_order_signature, orientation_signature)
    return (
        f"{details['n_blocks']} blocks across {details['query_chrom_alias_count']} normalized "
        f"query chromosome/scaffold aliases; reference-ordered blocks have "
        f"{details['order_pattern']} target-coordinate ranks and orientation counts "
        f"{details['orientation_counts']}."
    )


def read_published_summary() -> list[dict[str, object]]:
    rows = []
    with PUBLISHED_SUMMARY_PATH.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "species": row["Species"],
                    "genome": row["Genome"],
                    "chromosome": row["Chromosome"],
                    "start": clean_int(row["start"]),
                    "end": clean_int(row["end"]),
                    "description": row.get("Description", ""),
                    "source": str(PUBLISHED_SUMMARY_PATH.relative_to(REPO_ROOT)),
                }
            )
    return rows


def validate_against_published(
    intervals: list[dict[str, object]], published: list[dict[str, object]]
) -> list[dict[str, object]]:
    by_species = group_by(intervals, "species")
    validation = []
    for pub in published:
        species = str(pub["species"])
        candidates = by_species.get(species, [])
        if species == REFERENCE_SPECIES:
            candidates = [
                interval
                for subset in by_species.values()
                for interval in subset
                if interval["reference_chrom"] == REFERENCE_CHROM
            ]
            compare_space = "reference"
        else:
            compare_space = "query"
            candidates = [
                interval
                for interval in candidates
                if str(interval.get("query_chrom", "")) == str(pub["chromosome"])
            ]
        if not candidates:
            validation.append(validation_row(pub, None, compare_space, "no_reconstructed_interval"))
            continue
        best = min(candidates, key=lambda interval: endpoint_distance(pub, interval, compare_space))
        start_diff, end_diff = endpoint_differences(pub, best, compare_space)
        status = "exact_match" if start_diff == 0 and end_diff == 0 else (
            "approximate_match" if abs(start_diff) <= APPROX_MATCH_BP and abs(end_diff) <= APPROX_MATCH_BP else "different"
        )
        validation.append(validation_row(pub, best, compare_space, status))
    return validation


def endpoint_differences(
    published: dict[str, object], interval: dict[str, object], compare_space: str
) -> tuple[int, int]:
    if compare_space == "reference":
        rec_start, rec_end = int(interval["start"]), int(interval["end"])
    else:
        rec_start, rec_end = int(interval["query_start"]), int(interval["query_end"])
    return rec_start - int(published["start"]), rec_end - int(published["end"])


def endpoint_distance(
    published: dict[str, object], interval: dict[str, object], compare_space: str
) -> int:
    start_diff, end_diff = endpoint_differences(published, interval, compare_space)
    return abs(start_diff) + abs(end_diff)


def validation_row(
    published: dict[str, object],
    interval: dict[str, object] | None,
    compare_space: str,
    status: str,
) -> dict[str, object]:
    if interval is None:
        rec_start = rec_end = interval_id = ref_chrom = ""
        start_diff = end_diff = ""
    else:
        interval_id = interval["structural_interval_id"]
        ref_chrom = interval["reference_chrom"]
        if compare_space == "reference":
            rec_start, rec_end = interval["start"], interval["end"]
        else:
            rec_start, rec_end = interval["query_start"], interval["query_end"]
        start_diff, end_diff = endpoint_differences(published, interval, compare_space)
    return {
        "species": published["species"],
        "published_chrom": published["chromosome"],
        "published_start": published["start"],
        "published_end": published["end"],
        "reconstructed_interval_id": interval_id,
        "reconstructed_reference_chrom": ref_chrom,
        "reconstructed_start": rec_start,
        "reconstructed_end": rec_end,
        "comparison_coordinate_space": compare_space,
        "start_difference_bp": start_diff,
        "end_difference_bp": end_diff,
        "status": status,
        "published_description": published["description"],
    }


def validate_outputs(
    blocks: list[dict[str, object]], intervals: list[dict[str, object]], breakpoints: list[dict[str, object]]
) -> dict[str, object]:
    for row in blocks:
        if int(row["reference_start"]) >= int(row["reference_end"]):
            raise ValueError(f"Invalid reference interval for block {row['block_id']}")
        if int(row["query_start"]) >= int(row["query_end"]):
            raise ValueError(f"Invalid query interval for block {row['block_id']}")
        if row["reference_chrom"] != REFERENCE_CHROM:
            raise ValueError(f"Unexpected reference chromosome: {row['reference_chrom']}")
        if "random" in str(row["reference_chrom"]):
            raise ValueError(f"Random scaffold leaked into reference chrom: {row}")
        if row["orientation"] not in {"+", "-"}:
            raise ValueError(f"Unexpected orientation: {row['orientation']}")
        if not row["source_file"]:
            raise ValueError(f"Missing source provenance: {row}")
    ref_min = min(int(row["reference_start"]) for row in blocks)
    ref_max = max(int(row["reference_end"]) for row in blocks)
    for interval in intervals:
        if int(interval["start"]) >= int(interval["end"]):
            raise ValueError(f"Invalid interval: {interval}")
        if not (ref_min <= int(interval["start"]) <= ref_max and ref_min <= int(interval["end"]) <= ref_max):
            raise ValueError(f"Interval outside reference range: {interval}")
    overlap_count = 0
    for species, subset in group_by(intervals, "species").items():
        ordered = sorted(subset, key=lambda row: (row["start"], row["end"]))
        for a, b in zip(ordered, ordered[1:]):
            if int(a["end"]) > int(b["start"]):
                overlap_count += 1
    return {
        "n_structural_blocks": len(blocks),
        "n_structural_intervals": len(intervals),
        "n_breakpoints": len(breakpoints),
        "chr4_coordinate_min": ref_min,
        "chr4_coordinate_max": ref_max,
        "unknown_orientation_blocks": sum(row["orientation"] not in {"+", "-"} for row in blocks),
        "overlapping_structural_intervals": overlap_count,
        "species_with_published_summary_mismatch": "",
    }


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})


def format_value(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_outputs(
    blocks: list[dict[str, object]],
    intervals: list[dict[str, object]],
    breakpoints: list[dict[str, object]],
    states: list[dict[str, object]],
    state_defs: list[dict[str, object]],
    validation: list[dict[str, object]],
    summary: dict[str, object],
) -> None:
    write_tsv(
        BLOCKS_PATH,
        blocks,
        [
            "species",
            "genome_accession",
            "reference_species",
            "reference_assembly",
            "reference_chrom",
            "reference_start",
            "reference_end",
            "reference_midpoint",
            "query_assembly",
            "query_genome",
            "query_chrom",
            "query_start",
            "query_end",
            "query_midpoint",
            "block_id",
            "orientation",
            "reference_orientation",
            "order_index",
            "query_order_index",
            "block_length",
            "reference_block_length",
            "query_original_start",
            "query_original_end",
            "reference_original_start",
            "reference_original_end",
            "query_seq_id",
            "reference_seq_id",
            "query_seq_size",
            "source_file",
            "structural_interval_id",
        ],
    )
    write_tsv(
        INTERVALS_PATH,
        intervals,
        [
            "species",
            "structural_interval_id",
            "reference_chrom",
            "start",
            "end",
            "length_bp",
            "event_type",
            "n_blocks",
            "evidence",
            "source",
            "confidence_or_status",
            "query_chrom",
            "query_start",
            "query_end",
            "query_chrom_count",
            "structural_analysis_reference",
            "main_coordinate_reference",
        ],
    )
    write_tsv(
        BREAKPOINTS_PATH,
        breakpoints,
        [
            "species",
            "structural_interval_id",
            "breakpoint_id",
            "reference_chrom",
            "position",
            "breakpoint_type",
            "source",
        ],
    )
    write_tsv(
        STATES_PATH,
        states,
        [
            "species",
            "structural_interval_id",
            "arrangement_state",
            "block_order_signature",
            "orientation_signature",
        ],
    )
    write_tsv(
        STATE_DEFS_PATH,
        state_defs,
        [
            "arrangement_state",
            "n_blocks",
            "query_chrom_alias_count",
            "order_pattern",
            "orientation_counts",
            "example_intervals",
            "definition",
            "block_order_signature",
            "orientation_signature",
        ],
    )
    write_tsv(
        VALIDATION_PATH,
        validation,
        [
            "species",
            "published_chrom",
            "published_start",
            "published_end",
            "reconstructed_interval_id",
            "reconstructed_reference_chrom",
            "reconstructed_start",
            "reconstructed_end",
            "comparison_coordinate_space",
            "start_difference_bp",
            "end_difference_bp",
            "status",
            "published_description",
        ],
    )
    summary_rows = [{"metric": key, "value": value} for key, value in summary.items()]
    write_tsv(SUMMARY_PATH, summary_rows, ["metric", "value"])
    write_tsv(
        METADATA_PATH,
        [
            {"key": "main_coordinate_reference", "value": REFERENCE_DESCRIPTION},
            {"key": "main_coordinate_reference_assembly", "value": REFERENCE_ASSEMBLY},
            {"key": "structural_analysis_reference", "value": STRUCTURAL_ANALYSIS_REFERENCE},
            {"key": "primary_block_source", "value": str(TOTAB_PATH.relative_to(REPO_ROOT))},
            {"key": "published_structural_block_source", "value": str(PUBLISHED_BLOCKS_PATH.relative_to(REPO_ROOT))},
            {"key": "published_summary_source", "value": str(PUBLISHED_SUMMARY_PATH.relative_to(REPO_ROOT))},
            {"key": "interval_merge_gap_bp", "value": MERGE_GAP_BP},
            {
                "key": "coordinate_caveat",
                "value": "Main coordinates are exact GalGal6 chr4 only; chr4_*_random scaffolds are not used as chr4 coordinates.",
            },
        ],
        ["key", "value"],
    )


def make_plot(blocks: list[dict[str, object]], intervals: list[dict[str, object]], breakpoints: list[dict[str, object]]) -> None:
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chr4avian_mplconfig")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blocks_by_species = group_by(blocks, "species")
    intervals_by_species = group_by(intervals, "species")
    y_positions = {species: i for i, species in enumerate(PLOT_ORDER)}
    chrom_colors = {}
    palette = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#72b7b2", "#e45756", "#9d755d"]

    fig, ax = plt.subplots(figsize=(12, 5.5))
    for species in PLOT_ORDER:
        y = y_positions[species]
        subset = sorted(blocks_by_species.get(species, []), key=lambda r: r["reference_start"])
        for row in subset:
            chrom = str(row["query_chrom"])
            if chrom not in chrom_colors:
                chrom_colors[chrom] = palette[len(chrom_colors) % len(palette)]
            y_offset = 0.12 if row["orientation"] == "+" else -0.12
            ax.plot(
                [int(row["reference_start"]) / 1e6, int(row["reference_end"]) / 1e6],
                [y + y_offset, y + y_offset],
                color=chrom_colors[chrom],
                linewidth=2.0,
                alpha=0.85,
                solid_capstyle="butt",
            )
        for interval in intervals_by_species.get(species, []):
            ax.axvspan(
                int(interval["start"]) / 1e6,
                int(interval["end"]) / 1e6,
                ymin=(y - 0.38 + 0.5) / len(PLOT_ORDER),
                ymax=(y + 0.38 + 0.5) / len(PLOT_ORDER),
                color="#d9d9d9",
                alpha=0.45,
                linewidth=0,
            )
        for bp in [row for row in breakpoints if row["species"] == species]:
            ax.vlines(
                int(bp["position"]) / 1e6,
                y - 0.32,
                y + 0.32,
                color="#111111",
                linewidth=0.7,
                alpha=0.8,
            )
    ax.set_yticks([y_positions[species] for species in PLOT_ORDER], PLOT_ORDER)
    ax.set_xlabel("Chicken GalGal6 chromosome 4 coordinate (Mb)")
    ax.set_ylabel("")
    ax.set_ylim(-0.7, len(PLOT_ORDER) - 0.3)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.6)
    ax.set_title("Structural synteny blocks and candidate rearrangement intervals", loc="left")
    handles = [
        plt.Line2D([0], [0], color=color, lw=3, label=chrom)
        for chrom, color in sorted(chrom_colors.items())
    ]
    if handles:
        ax.legend(handles=handles, title="Query chromosome/scaffold", fontsize=7, title_fontsize=8, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=220)
    plt.close(fig)


def print_report(
    blocks: list[dict[str, object]],
    intervals: list[dict[str, object]],
    breakpoints: list[dict[str, object]],
    states: list[dict[str, object]],
    validation: list[dict[str, object]],
    summary: dict[str, object],
) -> None:
    print("Structural blocks per species:")
    for species in PLOT_ORDER:
        print(f"  {species}: {len([row for row in blocks if row['species'] == species])}")
    print("Inferred structural intervals:")
    for row in intervals:
        print(
            f"  {row['species']} {row['structural_interval_id']}: "
            f"{row['reference_chrom']}:{row['start']}-{row['end']} "
            f"({row['event_type']}, {row['n_blocks']} blocks)"
        )
    print("Breakpoints:")
    for row in breakpoints:
        print(f"  {row['breakpoint_id']}: {row['reference_chrom']}:{row['position']}")
    print("Arrangement states:")
    for row in states:
        print(f"  {row['species']} {row['structural_interval_id']}: {row['arrangement_state']}")
    print("Published summary validation:")
    for row in validation:
        print(
            f"  {row['species']} {row['published_chrom']}:{row['published_start']}-{row['published_end']} "
            f"-> {row['reconstructed_interval_id']} {row['status']}"
        )
    print("Validation summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    print("Files written:")
    for path in (
        BLOCKS_PATH,
        INTERVALS_PATH,
        BREAKPOINTS_PATH,
        STATES_PATH,
        STATE_DEFS_PATH,
        VALIDATION_PATH,
        SUMMARY_PATH,
        METADATA_PATH,
        PDF_PATH,
        PNG_PATH,
    ):
        print(f"  {path.relative_to(REPO_ROOT)}")


class StructuralTests(unittest.TestCase):
    def test_description_parser(self) -> None:
        parsed = parse_description("GCF_000002315.5_GalGal6.chr4")
        self.assertEqual(parsed["assembly"], "GCF_000002315")
        self.assertEqual(parsed["genome"], "5_GalGal6")
        self.assertEqual(parsed["chrom"], "chr4")

    def test_classify_inversion(self) -> None:
        group = [
            {
                "reference_start": 1,
                "reference_end": 10,
                "query_chrom": "chrA",
                "orientation": "-",
                "query_start": 30,
                "query_end": 40,
            },
            {
                "reference_start": 11,
                "reference_end": 20,
                "query_chrom": "chrA",
                "orientation": "-",
                "query_start": 10,
                "query_end": 20,
            },
        ]
        self.assertEqual(classify_event(group), "inversion")

    def test_overlap(self) -> None:
        self.assertTrue(ranges_overlap(10, 20, 20, 30))
        self.assertFalse(ranges_overlap(10, 19, 20, 30))


def run_tests() -> None:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(StructuralTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tests", action="store_true", help="run lightweight structural tests first")
    args = parser.parse_args()

    if args.run_tests:
        run_tests()

    blocks = load_paired_blocks()
    published_blocks = read_published_outlier_blocks()
    intervals = infer_intervals(published_blocks)
    assign_interval_ids_to_blocks(blocks, intervals)
    breakpoints = build_breakpoints(intervals)
    states, state_defs = build_arrangement_states(intervals, blocks)
    validation = validate_against_published(intervals, read_published_summary())
    summary = validate_outputs(blocks, intervals, breakpoints)
    summary["published_exact_matches"] = sum(row["status"] == "exact_match" for row in validation)
    summary["published_approximate_matches"] = sum(row["status"] == "approximate_match" for row in validation)
    summary["published_differences"] = sum(row["status"] == "different" for row in validation)
    summary["published_unmatched"] = sum(row["status"] == "no_reconstructed_interval" for row in validation)
    summary["species_with_published_summary_mismatch"] = ";".join(
        sorted(
            {
                str(row["species"])
                for row in validation
                if row["status"] in {"different", "no_reconstructed_interval"}
            }
            - {REFERENCE_SPECIES}
        )
    )
    write_outputs(blocks, intervals, breakpoints, states, state_defs, validation, summary)
    make_plot(blocks, intervals, breakpoints)
    print_report(blocks, intervals, breakpoints, states, validation, summary)


if __name__ == "__main__":
    main()
