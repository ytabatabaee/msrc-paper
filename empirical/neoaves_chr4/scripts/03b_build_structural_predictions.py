#!/usr/bin/env python3
"""Build Stage 3B structural-only quartet topology predictions."""

from __future__ import annotations

import argparse
import builtins
import csv
import hashlib
import json
import math
import statistics
import sys
import zipfile
from collections import defaultdict
from itertools import combinations
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

FORBIDDEN_INPUT_NAMES = {
    "locus_table_chr4.tsv",
    "structural_intervals.tsv",
    "canonical_structural_regions.tsv",
    "canonical_arrangement_states.tsv",
    "normalized_arrangement_states.tsv",
    "outlier-regions-by-maf2synteny.tsv",
    "outlier-regions-by-maf2synteny-summary.tsv",
}
FOCAL_SPECIES = ["Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo", "Stork"]
SPECIES_TO_TAXON = {
    "Flamingo": "Phoenicopterus_ruber",
    "Dove": "Columba_livia",
    "Sandgrouse": "Pterocles_gutturalis",
    "Turaco": "Tauraco_erythrolophus",
    "Cuckoo": "Cuculus_canorus",
    "Stork": "Ciconia_maguari",
}
CLADE_TRIPLES = {
    "Columbea": ("Columbimorphae.txt", "Phoenicopteriformes.txt", "Passera.txt"),
    "N61": ("Otidimorphae.txt", "Columbimorphae.txt", "ElementavesTelluraves.txt"),
    "N62": ("Columbiformes.txt", "OtherColumbimorphae.txt", "Otidimorphae.txt"),
}
TOPOLOGY_BY_PAIRING = {
    ("C1", "C2", "S", "O"): "q1",
    ("C1", "S", "C2", "O"): "q2",
    ("C2", "S", "C1", "O"): "q3",
}


def install_read_guard() -> None:
    real_open = builtins.open

    def guarded_open(file, *args, **kwargs):
        path = Path(file) if isinstance(file, (str, Path)) else None
        if path is not None and path.name in FORBIDDEN_INPUT_NAMES:
            raise RuntimeError(f"Stage 3B Phase A attempted forbidden input read: {path}")
        return real_open(file, *args, **kwargs)

    builtins.open = guarded_open


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value))
    return str(value)


def clade_members(filename: str) -> set[str]:
    with zipfile.ZipFile(EXAMINED_ZIP) as zf:
        return {line.strip() for line in zf.read(filename).decode().splitlines() if line.strip()}


def derive_role_mapping() -> tuple[list[dict[str, object]], dict[str, dict[str, list[str]]]]:
    all_taxa = set()
    with zipfile.ZipFile(EXAMINED_ZIP) as zf:
        for name in zf.namelist():
            if name.endswith(".txt"):
                all_taxa.update(line.strip() for line in zf.read(name).decode().splitlines() if line.strip())
    mapping_rows: list[dict[str, object]] = []
    mapping: dict[str, dict[str, list[str]]] = {}
    for clade, triple in CLADE_TRIPLES.items():
        role_members = {
            "C1": clade_members(triple[0]),
            "C2": clade_members(triple[1]),
            "S": clade_members(triple[2]),
        }
        role_members["O"] = all_taxa - role_members["C1"] - role_members["C2"] - role_members["S"]
        mapping[clade] = {}
        for role in ["C1", "C2", "S", "O"]:
            reps = [sp for sp, taxon in SPECIES_TO_TAXON.items() if taxon in role_members[role]]
            mapping[clade][role] = reps
            mapping_rows.append(
                {
                    "clade": clade,
                    "role": role,
                    "taxonomic_group": Path(triple[["C1", "C2", "S"].index(role)]).stem if role != "O" else "implicit_outgroup",
                    "structural_species": ";".join(reps),
                    "n_structural_representatives": len(reps),
                    "source_file": "genetreesupport/clade-analysis/score-rec-clade-runs.sh;examined-clades.zip",
                    "source_definition": " ".join(Path(x).stem for x in triple) if role != "O" else "all leaves not in C1/C2/S",
                    "mapping_status": "unambiguous" if reps else "unambiguous_no_focal_structural_representative",
                    "notes": "q1=(C1,C2)|(S,O); q2=(C1,S)|(C2,O); q3=(C2,S)|(C1,O)",
                }
            )
    return mapping_rows, mapping


def build_segments(set_name: str, chr_end: int) -> list[dict[str, object]]:
    set_col = f"in_{set_name}"
    bps = [
        row
        for row in read_tsv(SETS)
        if row.get(set_col) in {"1", "true", "True", "TRUE"}
    ]
    bps.sort(key=lambda r: int(r["reference_position"]))
    positions = [0] + [int(r["reference_position"]) for r in bps] + [chr_end]
    ids = ["chromosome_start"] + [r["consensus_breakpoint_id"] for r in bps] + ["chromosome_end"]
    rows = []
    for i in range(len(positions) - 1):
        start = positions[i] + 1 if i > 0 else 1
        end = positions[i + 1]
        if end < start:
            continue
        rows.append(
            {
                "segment_id": f"{set_name}_SEG_{i+1:03d}",
                "start": start,
                "end": end,
                "length_bp": end - start + 1,
                "left_breakpoint_id": ids[i],
                "right_breakpoint_id": ids[i + 1],
                "structural_set": set_name,
            }
        )
    return rows


def load_blocks() -> dict[str, list[dict[str, object]]]:
    rows: dict[str, list[dict[str, object]]] = {sp: [] for sp in FOCAL_SPECIES}
    for row in read_tsv(BLOCKS):
        if row["species"] not in FOCAL_SPECIES or row["reference_chrom"] != "chr4":
            continue
        rows[row["species"]].append(
            {
                "species": row["species"],
                "reference_start": int(row["reference_start"]),
                "reference_end": int(row["reference_end"]),
                "query_chrom": row["query_chrom"],
                "query_start": int(row["query_start"]),
                "query_end": int(row["query_end"]),
                "relative_orientation": row["relative_orientation"],
            }
        )
    for sp in rows:
        rows[sp].sort(key=lambda r: int(r["reference_start"]))
    return rows


def overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def anchor_projection(blocks: dict[str, list[dict[str, object]]], segment: dict[str, object], anchor_size: int) -> dict[str, dict[int, dict[str, object]]]:
    anchors = list(range(int(segment["start"]), int(segment["end"]) + 1, anchor_size))
    out: dict[str, dict[int, dict[str, object]]] = {sp: {} for sp in FOCAL_SPECIES}
    for idx, astart in enumerate(anchors, start=1):
        aend = min(astart + anchor_size - 1, int(segment["end"]))
        alen = aend - astart + 1
        for sp in FOCAL_SPECIES:
            hits = []
            for block in blocks[sp]:
                if int(block["reference_start"]) > aend:
                    break
                if int(block["reference_end"]) < astart:
                    continue
                ov = overlap(astart, aend, int(block["reference_start"]), int(block["reference_end"]))
                if ov:
                    qmid = (int(block["query_start"]) + int(block["query_end"])) / 2
                    hits.append((ov, qmid, block))
            if not hits:
                out[sp][idx] = {"status": "missing", "coverage_fraction": 0.0}
                continue
            hits.sort(key=lambda x: (-x[0], x[1]))
            best_ov, qmid, best = hits[0]
            ambiguous = len(hits) > 1 and hits[1][0] >= 0.8 * best_ov
            out[sp][idx] = {
                "status": "ambiguous" if ambiguous else "informative",
                "coverage_fraction": best_ov / alen,
                "query_chrom": best["query_chrom"],
                "query_mid": qmid,
                "orientation": best["relative_orientation"],
            }
    return out


def species_states(segment: dict[str, object], projections: dict[str, dict[int, dict[str, object]]], min_shared: int) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    rows = []
    states = {}
    total = len(next(iter(projections.values())))
    for sp, anchors in projections.items():
        informative = {aid: rec for aid, rec in anchors.items() if rec.get("status") == "informative" and rec.get("coverage_fraction", 0) >= 0.25}
        chrom_order = {chrom: i for i, chrom in enumerate(sorted({str(r["query_chrom"]) for r in informative.values()}))}
        ordered = sorted(informative.items(), key=lambda item: (chrom_order[str(item[1]["query_chrom"])], float(item[1]["query_mid"])))
        order_sig = ",".join(f"A{aid}" for aid, _ in ordered)
        orient_sig = ",".join(f"A{aid}:{rec['orientation']}" for aid, rec in sorted(informative.items()))
        confidence = "high"
        if total < min_shared:
            confidence = "short_segment_limited_anchors"
        elif len(informative) < min_shared:
            confidence = "low_informative_anchor_count"
        states[sp] = {"informative": informative, "order": [aid for aid, _ in ordered], "confidence": confidence}
        rows.append(
            {
                "segment_id": segment["segment_id"],
                "structural_set": segment["structural_set"],
                "anchor_size": len(anchors) and max(1, int(segment["length_bp"]) // len(anchors)),
                "species": sp,
                "relative_anchor_order": order_sig,
                "orientation_signature": orient_sig,
                "informative_anchor_count": len(informative),
                "missing_anchor_fraction": sum(1 for r in anchors.values() if r.get("status") == "missing") / total if total else 1,
                "ambiguous_anchor_fraction": sum(1 for r in anchors.values() if r.get("status") == "ambiguous") / total if total else 1,
                "confidence": confidence,
            }
        )
    return rows, states


def pair_similarity(segment: dict[str, object], states: dict[str, dict[str, object]], min_shared: int) -> tuple[list[dict[str, object]], dict[tuple[str, str], float]]:
    rows = []
    sim = {}
    for a, b in combinations(FOCAL_SPECIES, 2):
        ia = states[a]["informative"]
        ib = states[b]["informative"]
        shared = sorted(set(ia) & set(ib))
        possible_short = int(segment["length_bp"]) < min_shared * 100000
        if len(shared) < min_shared and not possible_short:
            order_c = orient_c = overall = None
            exact = False
            status = "insufficient_information"
        else:
            pos_a = {aid: i for i, aid in enumerate(states[a]["order"]) if aid in shared}
            pos_b = {aid: i for i, aid in enumerate(states[b]["order"]) if aid in shared}
            denom = concord = 0
            for x, y in combinations(shared, 2):
                denom += 1
                concord += (pos_a[x] - pos_a[y]) * (pos_b[x] - pos_b[y]) > 0
            order_c = concord / denom if denom else (1.0 if len(shared) >= 1 else None)
            orient_c = sum(ia[aid]["orientation"] == ib[aid]["orientation"] for aid in shared) / len(shared) if shared else None
            overall = (order_c + orient_c) / 2 if order_c is not None and orient_c is not None else None
            exact = order_c == 1.0 and orient_c == 1.0
            status = "exact_multi_anchor_order_orientation_match" if exact else "discordant_order_or_orientation"
        sim[(a, b)] = sim[(b, a)] = overall if overall is not None else float("nan")
        rows.append(
            {
                "segment_id": segment["segment_id"],
                "structural_set": segment["structural_set"],
                "anchor_size": segment["anchor_size"],
                "species1": a,
                "species2": b,
                "shared_informative_anchors": len(shared),
                "order_concordance": order_c,
                "orientation_concordance": orient_c,
                "missing_fraction": 1 - (len(shared) / max(1, len(set(ia) | set(ib)))),
                "overall_structural_similarity": overall,
                "exact_equivalent": exact,
                "status": status,
            }
        )
    return rows, sim


def discrete_states(segment: dict[str, object], states: dict[str, dict[str, object]], pair_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    exact = {(r["species1"], r["species2"]): r["exact_equivalent"] for r in pair_rows}
    exact.update({(r["species2"], r["species1"]): r["exact_equivalent"] for r in pair_rows})
    groups: list[list[str]] = []
    for sp in FOCAL_SPECIES:
        placed = False
        for group in groups:
            if all(exact.get((sp, other), False) for other in group):
                group.append(sp)
                placed = True
                break
        if not placed:
            groups.append([sp])
    rows = []
    for i, group in enumerate(groups, start=1):
        for sp in group:
            rows.append(
                {
                    "segment_id": segment["segment_id"],
                    "structural_set": segment["structural_set"],
                    "anchor_size": segment["anchor_size"],
                    "species": sp,
                    "arrangement_state": f"AS{i}",
                    "informative_anchor_count": len(states[sp]["informative"]),
                    "confidence": states[sp]["confidence"],
                }
            )
    return rows


def role_state(role_species: list[str], state_by_species: dict[str, str]) -> tuple[str | None, str]:
    if not role_species:
        return None, "missing_role_representative"
    vals = {state_by_species[sp] for sp in role_species if sp in state_by_species}
    if len(vals) == 1:
        return next(iter(vals)), "single_representative" if len(role_species) == 1 else "resolved"
    return None, "structurally_heterogeneous"


def strict_prediction(role_states: dict[str, str | None], role_status: dict[str, str]) -> tuple[str, str, str]:
    if any(role_states[r] is None for r in ["C1", "C2", "S", "O"]):
        return "unresolved/missing", "", "no_unique_discrete_prediction"
    if any(role_status[r] == "structurally_heterogeneous" for r in role_status):
        return "structurally_heterogeneous_role", "", "no_unique_discrete_prediction"
    vals = {r: role_states[r] for r in ["C1", "C2", "S", "O"]}
    patterns = [
        ("C1", "C2", "S", "O", "q1"),
        ("C1", "S", "C2", "O", "q2"),
        ("C2", "S", "C1", "O", "q3"),
    ]
    for a, b, c, d, q in patterns:
        if vals[a] == vals[b] and vals[c] == vals[d] and vals[a] != vals[c]:
            return f"{a}={b}!={c}={d}", q, "unique_2:2"
    counts = sorted([list(vals.values()).count(v) for v in set(vals.values())], reverse=True)
    return ":".join(map(str, counts)), "", "no_unique_discrete_prediction"


def role_similarity(role_a: list[str], role_b: list[str], sim: dict[tuple[str, str], float]) -> tuple[float, float, float]:
    vals = [sim.get((a, b), float("nan")) for a in role_a for b in role_b if a != b]
    vals = [v for v in vals if not math.isnan(v)]
    if not vals:
        return float("nan"), float("nan"), float("nan")
    return statistics.mean(vals), (statistics.pvariance(vals) if len(vals) > 1 else 0.0), max(vals) - min(vals)


def generalized_scores(mapping: dict[str, list[str]], sim: dict[tuple[str, str], float]) -> tuple[dict[str, float], str, float]:
    pairs = {}
    for a, b in combinations(["C1", "C2", "S", "O"], 2):
        pairs[(a, b)] = pairs[(b, a)] = role_similarity(mapping[a], mapping[b], sim)[0]
    def valid(*keys: tuple[str, str]) -> bool:
        return all(not math.isnan(pairs[k]) for k in keys)
    scores = {"G1": float("nan"), "G2": float("nan"), "G3": float("nan")}
    if valid(("C1", "C2"), ("S", "O"), ("C1", "S"), ("C1", "O"), ("C2", "S"), ("C2", "O")):
        scores["G1"] = (pairs[("C1", "C2")] + pairs[("S", "O")]) / 2 - (pairs[("C1", "S")] + pairs[("C1", "O")] + pairs[("C2", "S")] + pairs[("C2", "O")]) / 4
    if valid(("C1", "S"), ("C2", "O"), ("C1", "C2"), ("C1", "O"), ("S", "C2"), ("S", "O")):
        scores["G2"] = (pairs[("C1", "S")] + pairs[("C2", "O")]) / 2 - (pairs[("C1", "C2")] + pairs[("C1", "O")] + pairs[("S", "C2")] + pairs[("S", "O")]) / 4
    if valid(("C2", "S"), ("C1", "O"), ("C1", "C2"), ("C1", "S"), ("C2", "O"), ("S", "O")):
        scores["G3"] = (pairs[("C2", "S")] + pairs[("C1", "O")]) / 2 - (pairs[("C1", "C2")] + pairs[("C1", "S")] + pairs[("C2", "O")] + pairs[("S", "O")]) / 4
    finite = [(k, v) for k, v in scores.items() if not math.isnan(v)]
    if len(finite) != 3:
        return scores, "", float("nan")
    finite.sort(key=lambda x: x[1], reverse=True)
    pred = {"G1": "q1", "G2": "q2", "G3": "q3"}[finite[0][0]] if finite[0][1] > finite[1][1] else ""
    return scores, pred, finite[0][1] - finite[1][1]


def build_for(set_name: str, anchor_size: int, blocks: dict[str, list[dict[str, object]]], mapping: dict[str, dict[str, list[str]]], segments: list[dict[str, object]], min_shared: int):
    species_rows = []
    pair_rows_all = []
    disc_rows = []
    pred_rows = []
    for seg0 in segments:
        seg = dict(seg0)
        seg["anchor_size"] = anchor_size
        proj = anchor_projection(blocks, seg, anchor_size)
        srows, state_info = species_states(seg, proj, min_shared)
        prows, sim = pair_similarity(seg, state_info, min_shared)
        drows = discrete_states(seg, state_info, prows)
        state_by_species = {r["species"]: r["arrangement_state"] for r in drows}
        species_rows.extend(srows)
        pair_rows_all.extend(prows)
        disc_rows.extend(drows)
        for clade, rolemap in mapping.items():
            rs, rst = {}, {}
            for role in ["C1", "C2", "S", "O"]:
                rs[role], rst[role] = role_state(rolemap[role], state_by_species)
            pattern, pred, status = strict_prediction(rs, rst)
            scores, gpred, margin = generalized_scores(rolemap, sim)
            confidences = [state_info[sp]["confidence"] for reps in rolemap.values() for sp in reps if sp in state_info]
            pred_rows.append(
                {
                    "clade": clade,
                    "segment_id": seg["segment_id"],
                    "structural_set": set_name,
                    "anchor_size": anchor_size,
                    "start": seg["start"],
                    "end": seg["end"],
                    "strict_pattern": pattern,
                    "predicted_discrete_topology": pred,
                    "discrete_prediction_status": status,
                    "G1": scores["G1"],
                    "G2": scores["G2"],
                    "G3": scores["G3"],
                    "predicted_generalized_topology": gpred,
                    "generalized_margin": margin,
                    "prediction_stable_50_100_200kb": "",
                    "structural_confidence": "low" if any(c.startswith("low") for c in confidences) else "limited" if any(c.startswith("short") for c in confidences) else "high",
                    "role_status": ";".join(f"{r}:{rst[r]}" for r in ["C1", "C2", "S", "O"]),
                }
            )
    return species_rows, pair_rows_all, disc_rows, pred_rows


def annotate_stability(rows: list[dict[str, object]]) -> None:
    by_key = defaultdict(list)
    for row in rows:
        by_key[(row["structural_set"], row["clade"], row["segment_id"])].append(row)
    for group in by_key.values():
        disc = {r["predicted_discrete_topology"] for r in group if r["anchor_size"] in {50000, 100000, 200000}}
        gen = {r["predicted_generalized_topology"] for r in group if r["anchor_size"] in {50000, 100000, 200000}}
        stable = len(group) == 3 and len(disc) == 1 and len(gen) == 1
        for r in group:
            r["prediction_stable_50_100_200kb"] = stable


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-size", type=int, default=100000)
    parser.add_argument("--sensitivity-anchor-sizes", type=int, nargs="+", default=[50000, 100000, 200000])
    parser.add_argument("--min-shared-anchors", type=int, default=5)
    args = parser.parse_args()
    install_read_guard()
    for path in [BLOCKS, RUNS, SETS, CONSENSUS, STAGE27C_MANIFEST, EXAMINED_ZIP]:
        if not path.exists():
            raise FileNotFoundError(path)

    mapping_rows, mapping = derive_role_mapping()
    write_tsv(RESULTS / "quartet_role_mapping.tsv", mapping_rows, ["clade", "role", "taxonomic_group", "structural_species", "n_structural_representatives", "source_file", "source_definition", "mapping_status", "notes"])

    blocks = load_blocks()
    chr_end = max(int(r["reference_end"]) for rows in blocks.values() for r in rows)
    all_segment_rows = {}
    for set_name in ["primary", "stringent", "inclusive"]:
        segs = build_segments(set_name, chr_end)
        all_segment_rows[set_name] = segs
        out = RESULTS / f"stage3b_{set_name}_structural_segments.tsv"
        write_tsv(out, segs, ["segment_id", "start", "end", "length_bp", "left_breakpoint_id", "right_breakpoint_id"])

    all_species: list[dict[str, object]] = []
    all_pairs: list[dict[str, object]] = []
    all_disc: list[dict[str, object]] = []
    all_preds: list[dict[str, object]] = []
    for set_name, segs in all_segment_rows.items():
        for anchor_size in sorted(set(args.sensitivity_anchor_sizes + [args.anchor_size])):
            s, p, d, pr = build_for(set_name, anchor_size, blocks, mapping, segs, args.min_shared_anchors)
            all_species.extend(s)
            all_pairs.extend(p)
            all_disc.extend(d)
            all_preds.extend(pr)
    annotate_stability(all_preds)

    main_preds = [r for r in all_preds if r["structural_set"] == "primary" and int(r["anchor_size"]) == args.anchor_size]
    write_tsv(RESULTS / "segment_species_structural_states.tsv", [r for r in all_species if r["structural_set"] == "primary" and int(r["anchor_size"]) == args.anchor_size], ["segment_id", "structural_set", "anchor_size", "species", "relative_anchor_order", "orientation_signature", "informative_anchor_count", "missing_anchor_fraction", "ambiguous_anchor_fraction", "confidence"])
    write_tsv(RESULTS / "segment_pairwise_structural_similarity.tsv", [r for r in all_pairs if r["structural_set"] == "primary" and int(r["anchor_size"]) == args.anchor_size], ["segment_id", "structural_set", "anchor_size", "species1", "species2", "shared_informative_anchors", "order_concordance", "orientation_concordance", "missing_fraction", "overall_structural_similarity", "exact_equivalent", "status"])
    write_tsv(RESULTS / "segment_discrete_arrangement_states.tsv", [r for r in all_disc if r["structural_set"] == "primary" and int(r["anchor_size"]) == args.anchor_size], ["segment_id", "structural_set", "anchor_size", "species", "arrangement_state", "informative_anchor_count", "confidence"])
    pred_fields = ["clade", "segment_id", "structural_set", "anchor_size", "start", "end", "strict_pattern", "predicted_discrete_topology", "discrete_prediction_status", "G1", "G2", "G3", "predicted_generalized_topology", "generalized_margin", "prediction_stable_50_100_200kb", "structural_confidence", "role_status"]
    write_tsv(RESULTS / "structural_quartet_predictions.tsv", main_preds, pred_fields)
    write_tsv(RESULTS / "structural_quartet_predictions_all.tsv", all_preds, pred_fields)

    manifest = {
        "phase": "Stage 3B Phase A structural-only prediction freeze",
        "structural_input_checksums": {str(p.relative_to(REPO_ROOT)): sha256(p) for p in [BLOCKS, RUNS, SETS, CONSENSUS]},
        "stage27c_manifest_checksum": sha256(STAGE27C_MANIFEST),
        "structural_breakpoint_set": "primary",
        "anchor_size": args.anchor_size,
        "sensitivity_anchor_sizes": args.sensitivity_anchor_sizes,
        "role_mapping": mapping,
        "state_equivalence_rule": "complete-linkage all-pairs exact normalized anchor order and corrected reference-relative orientation on informative anchors",
        "generalized_score_formula": "G1/G2/G3 as predeclared pairwise role structural similarity contrasts",
        "frozen_topology_predictions_sha256": sha256(RESULTS / "structural_quartet_predictions.tsv"),
        "all_frozen_prediction_sensitivity_sha256": sha256(RESULTS / "structural_quartet_predictions_all.tsv"),
        "phase_a_script_sha256": sha256(Path(__file__)),
        "anti_circularity": {
            "forbidden_inputs": sorted(FORBIDDEN_INPUT_NAMES),
            "locus_table_read": False,
            "old_outlier_derived_arrangement_state_files_read": False,
        },
    }
    with (RESULTS / "stage3b_structural_prediction_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
