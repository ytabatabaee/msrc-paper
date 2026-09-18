#!/usr/bin/env python3
"""Build Stage 3B-v2 structural-only event-state predictions.

This phase intentionally freezes topology predictions before any chr4 QQS or
genealogy change-point files are read.
"""

from __future__ import annotations

import argparse
import builtins
import csv
import hashlib
import json
import math
import sys
import unittest
import zipfile
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / "results"
BLOCKS = RESULTS / "denovo_chr4_blocks.tsv"
RUNS = RESULTS / "denovo_synteny_runs.tsv"
SETS = RESULTS / "denovo_breakpoint_sets.tsv"
CONSENSUS = RESULTS / "denovo_consensus_structural_breakpoints.tsv"
STAGE27C_MANIFEST = RESULTS / "denovo_breakpoint_manifest.json"
CLADE_DIR = REPO_ROOT / "genetreesupport" / "clade-analysis"
EXAMINED_ZIP = CLADE_DIR / "examined-clades.zip"
NEWANNOTATION = CLADE_DIR / "newannotation.txt"

FORBIDDEN_INPUT_NAMES = {
    "locus_table_chr4.tsv",
    "genealogy_change_points.tsv",
    "genealogy_change_points_raw.tsv",
    "genealogy_change_point_cluster_audit.tsv",
    "shared_genealogy_transitions.tsv",
    "chr4_clade_means.tsv",
    "structural_quartet_predictions.tsv",
    "structural_quartet_predictions_all.tsv",
    "structural_topology_prediction_segment_results.tsv",
    "structural_topology_prediction_summary.tsv",
    "structural_topology_prediction_permutation.tsv",
    "segment_quartet_support.tsv",
    "segment_discrete_arrangement_states.tsv",
    "segment_species_structural_states.tsv",
    "segment_pairwise_structural_similarity.tsv",
    "arrangement_states.tsv",
    "arrangement_state_definitions.tsv",
    "normalized_arrangement_states.tsv",
    "canonical_arrangement_states.tsv",
    "q1",
    "q2",
    "q3",
}

QUERY_SPECIES = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]
REFERENCE_SPECIES = "Chicken"
SPECIES = QUERY_SPECIES + [REFERENCE_SPECIES]
SPECIES_TO_TAXON = {
    "Flamingo": "Phoenicopterus_ruber",
    "Dove": "Columba_livia",
    "Sandgrouse": "Pterocles_gutturalis",
    "Turaco": "Tauraco_erythrolophus",
    "Cuckoo": "Cuculus_canorus",
    "Stork": "Ciconia_maguari",
    "Chicken": "Gallus_gallus",
}
CLADE_TRIPLES = {
    "Columbea": ("Columbimorphae.txt", "Phoenicopteriformes.txt", "Passera.txt"),
    "N61": ("Otidimorphae.txt", "Columbimorphae.txt", "ElementavesTelluraves.txt"),
    "N62": ("Columbiformes.txt", "OtherColumbimorphae.txt", "Otidimorphae.txt"),
}
Q = ["q1", "q2", "q3"]
PRIMARY_WINDOW = 250_000
WINDOWS = [100_000, 250_000, 500_000]
SET_NAMES = ["primary", "stringent", "inclusive"]


def install_read_guard() -> None:
    real_open = builtins.open

    def guarded_open(file, *args, **kwargs):
        path = Path(file) if isinstance(file, (str, Path)) else None
        if path is not None and path.name in FORBIDDEN_INPUT_NAMES:
            raise RuntimeError(f"Stage 3B-v2 Phase A attempted forbidden input read: {path}")
        return real_open(file, *args, **kwargs)

    builtins.open = guarded_open


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


def clade_members(filename: str) -> set[str]:
    with zipfile.ZipFile(EXAMINED_ZIP) as zf:
        return {line.strip() for line in zf.read(filename).decode().splitlines() if line.strip()}


def annotated_taxa() -> set[str]:
    taxa = set()
    for row in read_tsv(NEWANNOTATION):
        if row:
            taxa.add(next(iter(row.values())))
    with NEWANNOTATION.open() as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if parts:
                taxa.add(parts[0])
    return taxa


def taxon_group(taxon: str) -> str:
    with NEWANNOTATION.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] == taxon:
                return parts[1]
    return "unannotated"


def derive_role_mapping() -> tuple[list[dict[str, object]], dict[str, dict[str, list[str]]]]:
    all_taxa = annotated_taxa()
    mapping_rows: list[dict[str, object]] = []
    mapping: dict[str, dict[str, list[str]]] = {}
    for clade, triple in CLADE_TRIPLES.items():
        members = {
            "C1": clade_members(triple[0]),
            "C2": clade_members(triple[1]),
            "S": clade_members(triple[2]),
        }
        members["O"] = all_taxa - members["C1"] - members["C2"] - members["S"]
        mapping[clade] = {}
        for role in ["C1", "C2", "S", "O"]:
            reps = [sp for sp in SPECIES if SPECIES_TO_TAXON[sp] in members[role]]
            mapping[clade][role] = reps
            if role == "O":
                group = "implicit_outgroup"
                source_def = "all newannotation taxa not in C1/C2/S"
            else:
                idx = ["C1", "C2", "S"].index(role)
                group = Path(triple[idx]).stem
                source_def = triple[idx]
            mapping_rows.append(
                {
                    "clade": clade,
                    "role": role,
                    "taxonomic_group": group,
                    "structural_species": ";".join(reps),
                    "structural_taxa": ";".join(SPECIES_TO_TAXON[s] for s in reps),
                    "n_structural_representatives": len(reps),
                    "includes_chicken_reference": REFERENCE_SPECIES in reps,
                    "chicken_taxonomic_group": taxon_group("Gallus_gallus") if REFERENCE_SPECIES in reps else "",
                    "source_file": "genetreesupport/clade-analysis/examined-clades.zip;genetreesupport/clade-analysis/newannotation.txt",
                    "source_definition": source_def,
                    "mapping_status": "resolved" if reps else "unresolved_no_structural_representative",
                    "topology_order": "q1=(C1,C2)|(S,O); q2=(C1,S)|(C2,O); q3=(C2,S)|(C1,O)",
                }
            )
    return mapping_rows, mapping


def load_runs() -> dict[str, list[dict[str, object]]]:
    by_species = {sp: [] for sp in QUERY_SPECIES}
    for row in read_tsv(RUNS):
        sp = row["species"]
        if sp not in by_species:
            continue
        rec = dict(row)
        for key in ["reference_start", "reference_end", "query_start", "query_end", "n_raw_blocks"]:
            rec[key] = int(rec[key])
        by_species[sp].append(rec)
    for sp in by_species:
        by_species[sp].sort(key=lambda r: (int(r["reference_start"]), int(r["reference_end"])))
    return by_species


def load_breakpoints() -> tuple[list[dict[str, object]], dict[str, dict[str, str]]]:
    meta = {r["consensus_breakpoint_id"]: r for r in read_tsv(CONSENSUS)}
    rows = []
    for row in read_tsv(SETS):
        rec = dict(row)
        rec["reference_position"] = int(rec["reference_position"])
        rows.append(rec)
    rows.sort(key=lambda r: int(r["reference_position"]))
    return rows, meta


def endpoint_near_breakpoint(run: dict[str, object], side: str) -> int:
    ori = str(run["relative_orientation"])
    if side == "left":
        return int(run["query_end"]) if ori == "+" else int(run["query_start"])
    return int(run["query_start"]) if ori == "+" else int(run["query_end"])


def gap_class(gap: int | None) -> str:
    if gap is None:
        return "unknown_gap"
    if gap <= 10_000:
        return "adjacent"
    if gap <= 250_000:
        return "near"
    if gap <= 1_000_000:
        return "distant"
    return "very_distant"


def nearest_flanks(runs: list[dict[str, object]], bp: int, window: int) -> tuple[dict[str, object] | None, dict[str, object] | None, int | None, int | None]:
    lefts = [r for r in runs if int(r["reference_start"]) < bp]
    rights = [r for r in runs if int(r["reference_end"]) > bp]
    left = max(lefts, key=lambda r: int(r["reference_end"]), default=None)
    right = min(rights, key=lambda r: int(r["reference_start"]), default=None)
    dl = max(0, bp - int(left["reference_end"])) if left else None
    dr = max(0, int(right["reference_start"]) - bp) if right else None
    if dl is not None and dl > window:
        left = None
    if dr is not None and dr > window:
        right = None
    return left, right, dl if left else None, dr if right else None


def junction_record(species: str, bp: dict[str, object], meta: dict[str, str], window: int, runs_by_species: dict[str, list[dict[str, object]]]) -> dict[str, object]:
    bp_id = str(bp["consensus_breakpoint_id"])
    pos = int(bp["reference_position"])
    if species == REFERENCE_SPECIES:
        return {
            "structural_set": "",
            "flank_window_bp": window,
            "species": species,
            "breakpoint_id": bp_id,
            "breakpoint_position": pos,
            "left_query_chrom": "GalGal6_chr4",
            "right_query_chrom": "GalGal6_chr4",
            "left_query_end_or_boundary": pos,
            "right_query_start_or_boundary": pos + 1,
            "left_relative_orientation": "+",
            "right_relative_orientation": "+",
            "same_query_chromosome": True,
            "query_order_relationship": "forward",
            "query_gap_or_jump": 1,
            "query_gap_class": "adjacent",
            "chromosome_switch": False,
            "orientation_switch": False,
            "order_reversal": False,
            "left_run_length": window,
            "right_run_length": window,
            "distance_left_run_to_breakpoint": 0,
            "distance_right_run_to_breakpoint": 0,
            "junction_signature": "same_chr:+:+:forward:adjacent",
            "confidence": "reference_defined",
        }
    left, right, dl, dr = nearest_flanks(runs_by_species[species], pos, window)
    if left is None or right is None:
        return {
            "structural_set": "",
            "flank_window_bp": window,
            "species": species,
            "breakpoint_id": bp_id,
            "breakpoint_position": pos,
            "left_query_chrom": left["query_chrom"] if left else "",
            "right_query_chrom": right["query_chrom"] if right else "",
            "left_query_end_or_boundary": endpoint_near_breakpoint(left, "left") if left else "",
            "right_query_start_or_boundary": endpoint_near_breakpoint(right, "right") if right else "",
            "left_relative_orientation": left["relative_orientation"] if left else "",
            "right_relative_orientation": right["relative_orientation"] if right else "",
            "same_query_chromosome": "",
            "query_order_relationship": "unresolved",
            "query_gap_or_jump": "",
            "query_gap_class": "unknown_gap",
            "chromosome_switch": "",
            "orientation_switch": "",
            "order_reversal": "",
            "left_run_length": int(left["reference_end"]) - int(left["reference_start"]) + 1 if left else "",
            "right_run_length": int(right["reference_end"]) - int(right["reference_start"]) + 1 if right else "",
            "distance_left_run_to_breakpoint": dl if dl is not None else "",
            "distance_right_run_to_breakpoint": dr if dr is not None else "",
            "junction_signature": "unresolved",
            "confidence": "unresolved_missing_flank",
        }
    lb = endpoint_near_breakpoint(left, "left")
    rb = endpoint_near_breakpoint(right, "right")
    same_chr = left["query_chrom"] == right["query_chrom"]
    jump = abs(rb - lb) if same_chr else None
    if not same_chr:
        order = "chromosome_switch"
    else:
        order = "forward" if rb >= lb else "reverse"
    chrom_switch = not same_chr
    orient_switch = left["relative_orientation"] != right["relative_orientation"]
    order_rev = same_chr and order == "reverse"
    sig = (
        f"{'same_chr' if same_chr else 'different_chr'}:"
        f"{left['relative_orientation']}:{right['relative_orientation']}:"
        f"{order}:{gap_class(jump)}"
    )
    llen = int(left["reference_end"]) - int(left["reference_start"]) + 1
    rlen = int(right["reference_end"]) - int(right["reference_start"]) + 1
    conf = "high"
    if min(llen, rlen) < 50_000:
        conf = "limited_short_flank_run"
    if (dl or 0) > window // 2 or (dr or 0) > window // 2:
        conf = "limited_distant_flank"
    return {
        "structural_set": "",
        "flank_window_bp": window,
        "species": species,
        "breakpoint_id": bp_id,
        "breakpoint_position": pos,
        "left_query_chrom": left["query_chrom"],
        "right_query_chrom": right["query_chrom"],
        "left_query_end_or_boundary": lb,
        "right_query_start_or_boundary": rb,
        "left_relative_orientation": left["relative_orientation"],
        "right_relative_orientation": right["relative_orientation"],
        "same_query_chromosome": same_chr,
        "query_order_relationship": order,
        "query_gap_or_jump": jump if jump is not None else "",
        "query_gap_class": gap_class(jump),
        "chromosome_switch": chrom_switch,
        "orientation_switch": orient_switch,
        "order_reversal": order_rev,
        "left_run_length": llen,
        "right_run_length": rlen,
        "distance_left_run_to_breakpoint": dl,
        "distance_right_run_to_breakpoint": dr,
        "junction_signature": sig,
        "confidence": conf,
    }


def build_junction_rows(bps: list[dict[str, object]], meta: dict[str, str], window: int, runs_by_species: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    rows = []
    for bp in bps:
        for species in SPECIES:
            rows.append(junction_record(species, bp, meta, window, runs_by_species))
    return rows


def group_states(junction_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    by_bp = defaultdict(list)
    for row in junction_rows:
        by_bp[(row["breakpoint_id"], row["flank_window_bp"])].append(row)
    for (_bp_id, _window), rows in sorted(by_bp.items()):
        ref_sig = next(r["junction_signature"] for r in rows if r["species"] == REFERENCE_SPECIES)
        derived_labels: dict[str, str] = {}
        next_idx = 1
        for sig in sorted({r["junction_signature"] for r in rows if r["species"] != REFERENCE_SPECIES and r["junction_signature"] not in {ref_sig, "unresolved"}}):
            derived_labels[sig] = f"A{next_idx}"
            next_idx += 1
        for row in rows:
            sig = row["junction_signature"]
            if sig == "unresolved":
                state = "unresolved"
                refder = "unresolved"
            elif sig == ref_sig:
                state = "A0"
                refder = "reference_like"
            else:
                state = derived_labels[sig]
                refder = state.replace("A", "derived_state_")
            out.append(
                {
                    "breakpoint_id": row["breakpoint_id"],
                    "breakpoint_position": row["breakpoint_position"],
                    "flank_window_bp": row["flank_window_bp"],
                    "species": row["species"],
                    "arrangement_state": state,
                    "junction_signature": sig,
                    "reference_or_derived": refder,
                    "confidence": row["confidence"],
                    "all_pairs_equivalence_rule": "identical normalized junction_signature within breakpoint/window; no single-linkage",
                }
            )
    return out


def event_type_from_pair(left_rows: list[dict[str, object]], right_rows: list[dict[str, object]], support: list[str]) -> str:
    inv = sum(bool(l["orientation_switch"]) and bool(r["orientation_switch"]) and bool(l["same_query_chromosome"]) and bool(r["same_query_chromosome"]) for l, r in paired_species_rows(left_rows, right_rows, support))
    trans = sum(bool(l["chromosome_switch"]) and bool(r["chromosome_switch"]) for l, r in paired_species_rows(left_rows, right_rows, support))
    if inv >= max(1, len(support) // 2):
        return "inversion_like"
    if trans >= max(1, len(support) // 2):
        return "translocation_like"
    return "adjacency_change"


def paired_species_rows(left_rows: list[dict[str, object]], right_rows: list[dict[str, object]], support: list[str]):
    lk = {r["species"]: r for r in left_rows}
    rk = {r["species"]: r for r in right_rows}
    for sp in support:
        yield lk[sp], rk[sp]


def pair_events(bps: list[dict[str, object]], junction_rows: list[dict[str, object]], structural_set: str, window: int) -> list[dict[str, object]]:
    bps = [b for b in bps if str(b.get(f"in_{structural_set}", "0")) in {"1", "true", "True"}]
    by_bp = defaultdict(list)
    for row in junction_rows:
        if int(row["flank_window_bp"]) == window:
            by_bp[row["breakpoint_id"]].append(row)
    used: set[str] = set()
    events = []
    idx = 1
    for a, b in zip(bps, bps[1:]):
        aid, bid = a["consensus_breakpoint_id"], b["consensus_breakpoint_id"]
        if aid in used or bid in used:
            continue
        left_rows = by_bp[aid]
        right_rows = by_bp[bid]
        support = []
        for sp in QUERY_SPECIES:
            lr = next(r for r in left_rows if r["species"] == sp)
            rr = next(r for r in right_rows if r["species"] == sp)
            left_changed = lr["junction_signature"] not in {"unresolved", "same_chr:+:+:forward:adjacent"}
            right_changed = rr["junction_signature"] not in {"unresolved", "same_chr:+:+:forward:adjacent"}
            evidence = (lr["orientation_switch"] and rr["orientation_switch"]) or (lr["chromosome_switch"] and rr["chromosome_switch"]) or (lr["order_reversal"] and rr["order_reversal"])
            if left_changed and right_changed and evidence:
                support.append(sp)
        if len(support) >= 2:
            evtype = event_type_from_pair(left_rows, right_rows, support)
            conf = "high" if len(support) >= 3 else "moderate"
            events.append(
                {
                    "structural_set": structural_set,
                    "flank_window_bp": window,
                    "event_id": f"{structural_set}_W{window}_EV{idx:03d}",
                    "left_breakpoint_id": aid,
                    "right_breakpoint_id": bid,
                    "start": a["reference_position"],
                    "end": b["reference_position"],
                    "event_type": evtype,
                    "species_support": ";".join(support),
                    "pairing_evidence": "adjacent frozen DSBs with same species supporting orientation/chromosome/order change at both boundaries",
                    "pairing_confidence": conf,
                }
            )
            used.update({aid, bid})
            idx += 1
    for bp in bps:
        bid = bp["consensus_breakpoint_id"]
        if bid in used:
            continue
        rows = by_bp[bid]
        support = [
            r["species"]
            for r in rows
            if r["species"] in QUERY_SPECIES and r["junction_signature"] not in {"unresolved", "same_chr:+:+:forward:adjacent"}
        ]
        events.append(
            {
                "structural_set": structural_set,
                "flank_window_bp": window,
                "event_id": f"{structural_set}_W{window}_EV{idx:03d}",
                "left_breakpoint_id": bid,
                "right_breakpoint_id": "",
                "start": bp["reference_position"],
                "end": bp["reference_position"],
                "event_type": "single_boundary_event",
                "species_support": ";".join(support),
                "pairing_evidence": "no second frozen boundary met the paired-event evidence rule",
                "pairing_confidence": "single_boundary",
            }
        )
        idx += 1
    events.sort(key=lambda r: (int(r["start"]), int(r["end"]), r["event_id"]))
    return events


def internal_signature(species: str, event: dict[str, object], runs_by_species: dict[str, list[dict[str, object]]]) -> str:
    if not event["right_breakpoint_id"]:
        return "single_boundary"
    if species == REFERENCE_SPECIES:
        return "GalGal6_chr4:+:reference_interval"
    start, end = int(event["start"]), int(event["end"])
    mid = (start + end) // 2
    hits = [r for r in runs_by_species[species] if int(r["reference_start"]) <= mid <= int(r["reference_end"])]
    if not hits:
        return "unresolved_internal"
    best = max(hits, key=lambda r: int(r["reference_end"]) - int(r["reference_start"]))
    return f"{best['query_chrom']}:{best['relative_orientation']}:coherent_run"


def event_species_states(events: list[dict[str, object]], arrangement_rows: list[dict[str, object]], runs_by_species: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    state_key = {(r["breakpoint_id"], int(r["flank_window_bp"]), r["species"]): r for r in arrangement_rows}
    out = []
    for event in events:
        window = int(event["flank_window_bp"])
        for sp in SPECIES:
            left = state_key[(event["left_breakpoint_id"], window, sp)]
            right = state_key.get((event["right_breakpoint_id"], window, sp)) if event["right_breakpoint_id"] else None
            internal = internal_signature(sp, event, runs_by_species)
            if right is not None:
                combined = f"{left['junction_signature']}|{internal}|{right['junction_signature']}"
                confs = [left["confidence"], right["confidence"]]
                if internal == "unresolved_internal":
                    confs.append("unresolved_internal")
            else:
                combined = str(left["junction_signature"])
                confs = [left["confidence"]]
            if sp == REFERENCE_SPECIES:
                combined_state = "A0"
            elif "unresolved" in combined:
                combined_state = "unresolved"
            else:
                combined_state = combined
            confidence = "high"
            if any(str(c).startswith("unresolved") for c in confs):
                confidence = "unresolved"
            elif any(str(c).startswith("limited") for c in confs):
                confidence = "limited"
            out.append(
                {
                    "structural_set": event["structural_set"],
                    "flank_window_bp": window,
                    "event_id": event["event_id"],
                    "species": sp,
                    "left_state": left["junction_signature"],
                    "internal_state": internal,
                    "right_state": right["junction_signature"] if right else "",
                    "combined_arrangement_state": combined_state,
                    "confidence": confidence,
                }
            )
    relabel_event_states(out)
    return out


def relabel_event_states(rows: list[dict[str, object]]) -> None:
    by_event = defaultdict(list)
    for row in rows:
        by_event[(row["event_id"], row["flank_window_bp"])].append(row)
    for group in by_event.values():
        labels: dict[str, str] = {}
        next_idx = 1
        for sig in sorted({r["combined_arrangement_state"] for r in group if r["species"] != REFERENCE_SPECIES and r["combined_arrangement_state"] not in {"A0", "unresolved"}}):
            labels[sig] = f"A{next_idx}"
            next_idx += 1
        for row in group:
            sig = row["combined_arrangement_state"]
            if sig == "A0":
                row["combined_arrangement_state"] = "A0"
            elif sig == "unresolved":
                row["combined_arrangement_state"] = "unresolved"
            else:
                row["combined_arrangement_state"] = labels[sig]


def role_state(role_species: list[str], state_by_species: dict[str, str]) -> tuple[str, str]:
    if not role_species:
        return "unresolved", "unresolved_no_structural_representative"
    vals = [state_by_species.get(sp, "unresolved") for sp in role_species]
    if any(v == "unresolved" for v in vals):
        return "unresolved", "unresolved_role_species_state"
    uniq = set(vals)
    if len(uniq) == 1:
        return vals[0], "resolved"
    return "role_heterogeneous", "role_heterogeneous"


def strict_prediction(role_states: dict[str, str]) -> tuple[str, str, str]:
    if any(role_states[r] in {"unresolved", "role_heterogeneous"} for r in ["C1", "C2", "S", "O"]):
        return "unresolved_or_heterogeneous", "no_unique_strict_prediction", "ineligible"
    if role_states["C1"] == role_states["C2"] and role_states["S"] == role_states["O"] and role_states["C1"] != role_states["S"]:
        return "C1=C2!=S=O", "q1", "strict_2:2"
    if role_states["C1"] == role_states["S"] and role_states["C2"] == role_states["O"] and role_states["C1"] != role_states["C2"]:
        return "C1=S!=C2=O", "q2", "strict_2:2"
    if role_states["C2"] == role_states["S"] and role_states["C1"] == role_states["O"] and role_states["C2"] != role_states["C1"]:
        return "C2=S!=C1=O", "q3", "strict_2:2"
    counts = sorted([list(role_states.values()).count(v) for v in set(role_states.values())], reverse=True)
    return ":".join(map(str, counts)), "no_unique_strict_prediction", "ineligible"


def prediction_rows(events: list[dict[str, object]], species_state_rows: list[dict[str, object]], mapping: dict[str, dict[str, list[str]]]) -> list[dict[str, object]]:
    by_event_species = {(r["event_id"], r["species"]): r for r in species_state_rows}
    event_meta = {r["event_id"]: r for r in events}
    out = []
    for event in events:
        states = {sp: by_event_species[(event["event_id"], sp)]["combined_arrangement_state"] for sp in SPECIES}
        confidences = {sp: by_event_species[(event["event_id"], sp)]["confidence"] for sp in SPECIES}
        for clade, roles in mapping.items():
            rs, status = {}, {}
            for role in ["C1", "C2", "S", "O"]:
                rs[role], status[role] = role_state(roles[role], states)
            pattern, pred, pred_status = strict_prediction(rs)
            role_confs = [confidences[sp] for role in roles.values() for sp in role if sp in confidences]
            structural_conf = "high"
            if any(c == "unresolved" for c in role_confs):
                structural_conf = "unresolved"
            elif any(c == "limited" for c in role_confs):
                structural_conf = "limited"
            out.append(
                {
                    "structural_set": event["structural_set"],
                    "flank_window_bp": event["flank_window_bp"],
                    "clade": clade,
                    "event_id": event["event_id"],
                    "start": event["start"],
                    "end": event["end"],
                    "event_type": event["event_type"],
                    "pairing_confidence": event["pairing_confidence"],
                    "role_C1_state": rs["C1"],
                    "role_C2_state": rs["C2"],
                    "role_S_state": rs["S"],
                    "role_O_state": rs["O"],
                    "role_status": ";".join(f"{r}:{status[r]}" for r in ["C1", "C2", "S", "O"]),
                    "arrangement_pattern": pattern,
                    "predicted_topology": pred,
                    "prediction_status": pred_status if pred in Q else "no_unique_strict_prediction",
                    "structural_confidence": structural_conf,
                }
            )
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args(argv)
    if args.run_tests:
        return run_tests()
    install_read_guard()
    for path in [BLOCKS, RUNS, SETS, CONSENSUS, STAGE27C_MANIFEST, EXAMINED_ZIP, NEWANNOTATION]:
        if not path.exists():
            raise FileNotFoundError(path)
    mapping_rows, mapping = derive_role_mapping()
    write_tsv(
        RESULTS / "stage3b_v2_quartet_role_mapping.tsv",
        mapping_rows,
        ["clade", "role", "taxonomic_group", "structural_species", "structural_taxa", "n_structural_representatives", "includes_chicken_reference", "chicken_taxonomic_group", "source_file", "source_definition", "mapping_status", "topology_order"],
    )
    bps, meta = load_breakpoints()
    runs = load_runs()
    all_junction_rows = []
    all_arrangement_rows = []
    all_events = []
    all_event_species = []
    all_predictions = []
    for window in WINDOWS:
        junctions = build_junction_rows(bps, meta, window, runs)
        states = group_states(junctions)
        all_junction_rows.extend(junctions)
        all_arrangement_rows.extend(states)
        for structural_set in SET_NAMES:
            set_events = pair_events(bps, junctions, structural_set, window)
            all_events.extend(set_events)
            esp = event_species_states(set_events, states, runs)
            all_event_species.extend(esp)
            all_predictions.extend(prediction_rows(set_events, esp, mapping))
    for rows in [all_junction_rows, all_arrangement_rows]:
        for row in rows:
            in_sets = [s for s in SET_NAMES if any(str(bp["consensus_breakpoint_id"]) == str(row.get("breakpoint_id")) and str(bp.get(f"in_{s}", "0")) in {"1", "true", "True"} for bp in bps)]
            row["structural_set"] = ";".join(in_sets)

    primary_junction = [r for r in all_junction_rows if int(r["flank_window_bp"]) == PRIMARY_WINDOW and "primary" in str(r["structural_set"]).split(";")]
    primary_arrangement = [r for r in all_arrangement_rows if int(r["flank_window_bp"]) == PRIMARY_WINDOW and "primary" in str(r["structural_set"]).split(";")]
    primary_events = [r for r in all_events if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == PRIMARY_WINDOW]
    primary_event_species = [r for r in all_event_species if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == PRIMARY_WINDOW]
    primary_predictions = [r for r in all_predictions if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == PRIMARY_WINDOW]

    write_tsv(RESULTS / "stage3b_v2_breakpoint_junction_states.tsv", primary_junction, ["structural_set", "flank_window_bp", "species", "breakpoint_id", "breakpoint_position", "left_query_chrom", "right_query_chrom", "left_query_end_or_boundary", "right_query_start_or_boundary", "left_relative_orientation", "right_relative_orientation", "same_query_chromosome", "query_order_relationship", "query_gap_or_jump", "query_gap_class", "chromosome_switch", "orientation_switch", "order_reversal", "left_run_length", "right_run_length", "distance_left_run_to_breakpoint", "distance_right_run_to_breakpoint", "junction_signature", "confidence"])
    write_tsv(RESULTS / "stage3b_v2_event_arrangement_states.tsv", primary_arrangement, ["breakpoint_id", "breakpoint_position", "flank_window_bp", "species", "arrangement_state", "junction_signature", "reference_or_derived", "confidence", "all_pairs_equivalence_rule"])
    write_tsv(RESULTS / "stage3b_v2_structural_events.tsv", primary_events, ["structural_set", "flank_window_bp", "event_id", "left_breakpoint_id", "right_breakpoint_id", "start", "end", "event_type", "species_support", "pairing_evidence", "pairing_confidence"])
    write_tsv(RESULTS / "stage3b_v2_event_species_states.tsv", primary_event_species, ["structural_set", "flank_window_bp", "event_id", "species", "left_state", "internal_state", "right_state", "combined_arrangement_state", "confidence"])
    pred_fields = ["structural_set", "flank_window_bp", "clade", "event_id", "start", "end", "event_type", "pairing_confidence", "role_C1_state", "role_C2_state", "role_S_state", "role_O_state", "role_status", "arrangement_pattern", "predicted_topology", "prediction_status", "structural_confidence"]
    write_tsv(RESULTS / "stage3b_v2_structural_quartet_predictions.tsv", primary_predictions, pred_fields)
    write_tsv(RESULTS / "stage3b_v2_structural_quartet_predictions_all.tsv", all_predictions, pred_fields)

    manifest = {
        "phase": "Stage 3B-v2 Phase A structural-only event-state prediction freeze",
        "structural_input_checksums": {str(p.relative_to(REPO_ROOT)): sha256(p) for p in [BLOCKS, RUNS, SETS, CONSENSUS]},
        "stage27c_manifest_checksum": sha256(STAGE27C_MANIFEST),
        "quartet_role_mapping": mapping,
        "quartet_role_mapping_sha256": sha256(RESULTS / "stage3b_v2_quartet_role_mapping.tsv"),
        "chicken_reference_handling": "Chicken/Gallus_gallus/GalGal6 is represented explicitly as a structural taxon with direct reference state A0.",
        "junction_state_definition": "same/different query chromosome, reference-relative left/right orientations, query order relationship, and query discontinuity class; raw query coordinates are audit fields only.",
        "breakpoint_window_sizes_bp": WINDOWS,
        "primary_breakpoint_window_bp": PRIMARY_WINDOW,
        "event_pairing_rule": "Adjacent frozen breakpoints are paired only when the same species support both boundaries and orientation/chromosome/order evidence occurs at both; otherwise single_boundary_event.",
        "all_pairs_equivalence_rule": "Identical normalized state signature within event/breakpoint; no connected-component or single-linkage grouping.",
        "strict_prediction_rule": "C1=C2!=S=O -> q1; C1=S!=C2=O -> q2; C2=S!=C1=O -> q3; all other patterns no_unique_strict_prediction.",
        "frozen_topology_predictions": str((RESULTS / "stage3b_v2_structural_quartet_predictions.tsv").relative_to(REPO_ROOT)),
        "frozen_topology_predictions_sha256": sha256(RESULTS / "stage3b_v2_structural_quartet_predictions.tsv"),
        "all_frozen_prediction_sensitivity_sha256": sha256(RESULTS / "stage3b_v2_structural_quartet_predictions_all.tsv"),
        "script_sha256": sha256(Path(__file__)),
        "anti_circularity": {
            "forbidden_inputs": sorted(FORBIDDEN_INPUT_NAMES),
            "genealogy_or_q_values_loaded": False,
            "old_outlier_or_stage3b_state_files_read": False,
        },
    }
    with (RESULTS / "stage3b_v2_structural_prediction_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


class RegressionTests(unittest.TestCase):
    def test_chicken_reference_signature(self):
        rec = junction_record("Chicken", {"consensus_breakpoint_id": "DSB_X", "reference_position": 1000}, {}, 250, {})
        self.assertEqual(rec["junction_signature"], "same_chr:+:+:forward:adjacent")
        self.assertEqual(rec["confidence"], "reference_defined")

    def test_strict_prediction_only_true_two_two(self):
        self.assertEqual(strict_prediction({"C1": "A", "C2": "A", "S": "B", "O": "B"})[1], "q1")
        self.assertEqual(strict_prediction({"C1": "A", "C2": "A", "S": "A", "O": "B"})[1], "no_unique_strict_prediction")
        self.assertEqual(strict_prediction({"C1": "A", "C2": "B", "S": "C", "O": "D"})[1], "no_unique_strict_prediction")

    def test_gap_class(self):
        self.assertEqual(gap_class(1), "adjacent")
        self.assertEqual(gap_class(100_001), "near")


def run_tests() -> int:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RegressionTests)
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
