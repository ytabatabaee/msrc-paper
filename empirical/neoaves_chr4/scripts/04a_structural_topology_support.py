#!/usr/bin/env python3
"""Stage 4A: topology support in frozen chr4 structural neighborhoods.

Structural masks are defined only from frozen Stage 3 structural breakpoints:
the union of +/-250 kb neighborhoods around breakpoints in each frozen set.
No topology data are read while masks are constructed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/neoaves_chr4/processed"
RESULTS = ROOT / "empirical/neoaves_chr4/results"
FIGURES = ROOT / "empirical/neoaves_chr4/figures"

LOCI = DATA / "locus_table_chr4.tsv"
BREAKPOINTS = DATA / "denovo_breakpoint_sets.tsv"
EVENTS = DATA / "stage3b_v2_structural_events.tsv"
STAGE3B_SUMMARY = DATA / "stage3b_v2_summary.tsv"
SEGMENTS = {s: DATA / f"stage3b_{s}_structural_segments.tsv" for s in ("stringent", "primary", "inclusive")}

SETS = ("stringent", "primary", "inclusive")
CLADES = ("Columbea", "N61", "N62")
TOPOLOGIES = ("q1", "q2", "q3")
SPECIES_TOPOLOGY = {c: "q1" for c in CLADES}
MASK_RADIUS = 250_000
CHR_START = 1
CHR_END = 91_310_470
SEED = 1729


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.12g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: fmt(row.get(f, "")) for f in fields})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for start, end in sorted(intervals):
        if not out or start > out[-1][1] + 1:
            out.append([start, end])
        else:
            out[-1][1] = max(out[-1][1], end)
    return [(a, b) for a, b in out]


def build_masks(breakpoint_rows: list[dict[str, str]], radius: int = MASK_RADIUS) -> dict[str, list[tuple[int, int]]]:
    """Build masks using structural columns only; callers need not supply loci."""
    masks = {}
    for set_name in SETS:
        positions = [int(r["reference_position"]) for r in breakpoint_rows if truthy(r[f"in_{set_name}"])]
        masks[set_name] = merge_intervals([
            (max(CHR_START, p - radius), min(CHR_END, p + radius)) for p in positions
        ])
    return masks


def in_intervals(position: float, intervals: list[tuple[int, int]]) -> bool:
    return any(start <= position <= end for start, end in intervals)


def interval_span(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def load_tracks() -> dict[str, list[dict[str, object]]]:
    tracks = {c: [] for c in CLADES}
    for row in read_tsv(LOCI):
        if row["Chromosome"] != "chr4" or row["clade"] not in tracks:
            continue
        rec: dict[str, object] = {
            "gene": row["Gene"], "clade": row["clade"], "start": int(row["start"]),
            "end": int(row["end"]), "midpoint": float(row["midpoint"]),
            "dominant_topology": row["dominant_topology"],
        }
        for q in TOPOLOGIES:
            rec[q] = float(row[q])
        tracks[row["clade"]].append(rec)
    for track in tracks.values():
        track.sort(key=lambda r: (float(r["midpoint"]), int(r["start"]), str(r["gene"])))
    return tracks


def collapse_topology_runs(track: list[dict[str, object]]) -> list[dict[str, object]]:
    runs: list[dict[str, object]] = []
    for row in track:
        state = str(row["dominant_topology"])
        if not runs or runs[-1]["topology"] != state:
            runs.append({"topology": state, "start": int(row["start"]), "end": int(row["end"]), "n_windows": 1})
        else:
            runs[-1]["end"] = int(row["end"])
            runs[-1]["n_windows"] = int(runs[-1]["n_windows"]) + 1
    return runs


def topology_summary(rows: list[dict[str, object]], species: str) -> dict[str, object]:
    counts = Counter(str(r["dominant_topology"]) for r in rows)
    n = len(rows)
    freqs = {q: counts[q] / n if n else math.nan for q in TOPOLOGIES}
    support = {q: sum(float(r[q]) for r in rows) for q in TOPOLOGIES}
    mean_support = {q: support[q] / n if n else math.nan for q in TOPOLOGIES}
    alternatives = [q for q in TOPOLOGIES if q != species]
    best_alt = max(alternatives, key=lambda q: (freqs[q], q)) if n else ""
    dominant = max(TOPOLOGIES, key=lambda q: (freqs[q], q)) if n else ""
    return {
        "n": n, "counts": counts, "freqs": freqs, "support": support,
        "mean_support": mean_support, "dominant": dominant, "best_alt": best_alt,
        "margin": freqs[species] - freqs[best_alt] if n else math.nan,
        "support_margin": mean_support[species] - mean_support[best_alt] if n else math.nan,
    }


def strongest_competitor(track: list[dict[str, object]], species: str) -> str:
    counts = Counter(str(r["dominant_topology"]) for r in track)
    return max((q for q in TOPOLOGIES if q != species), key=lambda q: (counts[q], q))


def structural_support_rows(tracks: dict[str, list[dict[str, object]]], masks: dict[str, list[tuple[int, int]]]) -> list[dict[str, object]]:
    out = []
    genome_span = CHR_END - CHR_START + 1
    for set_name in SETS:
        span = interval_span(masks[set_name])
        for clade in CLADES:
            track = tracks[clade]
            species = SPECIES_TOPOLOGY[clade]
            competitor = strongest_competitor(track, species)
            total_alt = sum(str(r["dominant_topology"]) != species for r in track)
            for category, associated in (("structural-associated", True), ("structural-background", False)):
                rows = [r for r in track if in_intervals(float(r["midpoint"]), masks[set_name]) == associated]
                s = topology_summary(rows, species)
                alt_here = sum(str(r["dominant_topology"]) != species for r in rows)
                out.append({
                    "structural_set": set_name, "mask_definition": "union_of_frozen_breakpoint_neighborhoods",
                    "mask_radius_bp": MASK_RADIUS, "clade": clade, "species_topology": species,
                    "strongest_competing_topology": competitor, "structural_category": category,
                    "n_windows": s["n"],
                    **{f"n_{q}": s["counts"][q] for q in TOPOLOGIES},
                    **{f"frequency_{q}": s["freqs"][q] for q in TOPOLOGIES},
                    **{f"mean_support_{q}": s["mean_support"][q] for q in TOPOLOGIES},
                    "dominant_topology": s["dominant"], "best_alternative": s["best_alt"],
                    "species_vs_best_alternative_margin": s["margin"],
                    "mean_support_species_vs_best_alternative_margin": s["support_margin"],
                    "n_alternative_topology_windows": alt_here,
                    "fraction_all_alternative_windows": alt_here / total_alt if total_alt else math.nan,
                    "structural_span_bp": span if associated else genome_span - span,
                    "fraction_genomic_span": span / genome_span if associated else 1 - span / genome_span,
                })
    return out


def window_block_rows(tracks: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    out = []
    for clade in CLADES:
        species = SPECIES_TOPOLOGY[clade]
        competitor = strongest_competitor(tracks[clade], species)
        window_counts = Counter(str(r["dominant_topology"]) for r in tracks[clade])
        runs = collapse_topology_runs(tracks[clade])
        block_counts = Counter(str(r["topology"]) for r in runs)
        for scheme, counts in (("window-weighted", window_counts), ("block-normalized", block_counts)):
            n = sum(counts.values())
            freqs = {q: counts[q] / n for q in TOPOLOGIES}
            best_alt = max((q for q in TOPOLOGIES if q != species), key=lambda q: (freqs[q], q))
            out.append({
                "clade": clade, "species_topology": species, "strongest_competing_topology": competitor,
                "weighting_scheme": scheme, "n_evidence_units": n,
                **{f"n_{q}": counts[q] for q in TOPOLOGIES},
                **{f"q_{q[1:]}": freqs[q] for q in TOPOLOGIES},
                "dominant_topology": max(TOPOLOGIES, key=lambda q: (freqs[q], q)),
                "best_alternative": best_alt,
                "species_vs_best_alternative_margin": freqs[species] - freqs[best_alt],
                "empirical_competing_topology_fraction": freqs[competitor],
            })
    return out


def enrichment_stat(track: list[dict[str, object]], mask: list[tuple[int, int]], competitor: str, states: list[str] | None = None, assoc: list[int] | None = None, n_comp: int | None = None) -> dict[str, float]:
    states = states or [str(r["dominant_topology"]) for r in track]
    assoc = assoc if assoc is not None else [i for i, r in enumerate(track) if in_intervals(float(r["midpoint"]), mask)]
    n_comp = n_comp if n_comp is not None else sum(s == competitor for s in states)
    assoc_comp = sum(states[i] == competitor for i in assoc)
    fa = assoc_comp / len(assoc) if assoc else math.nan
    n_back = len(track) - len(assoc)
    fb = (n_comp - assoc_comp) / n_back if n_back else math.nan
    return {
        "assoc_fraction": fa, "background_fraction": fb, "difference": fa - fb,
        "fold_enrichment": fa / fb if fb else math.inf,
        "fraction_competing_windows_in_associated": assoc_comp / n_comp if n_comp else math.nan,
        "fraction_windows_associated": len(assoc) / len(track),
    }


def stable_offset(*parts: str) -> int:
    return int(hashlib.sha256("|".join(parts).encode()).hexdigest()[:8], 16)


def circular_enrichment(track: list[dict[str, object]], mask: list[tuple[int, int]], competitor: str, permutations: int, seed: int) -> tuple[dict[str, float], list[float]]:
    assoc = [i for i, r in enumerate(track) if in_intervals(float(r["midpoint"]), mask)]
    states = [str(r["dominant_topology"]) for r in track]
    n_comp = sum(s == competitor for s in states)
    observed = enrichment_stat(track, mask, competitor, states, assoc, n_comp)
    rng = random.Random(seed)
    null = []
    for _ in range(permutations):
        shift = rng.randrange(len(states))
        shifted = states[shift:] + states[:shift]
        null.append(enrichment_stat(track, mask, competitor, shifted, assoc, n_comp)["difference"])
    p = (1 + sum(x >= observed["difference"] for x in null)) / (permutations + 1)
    return observed, null + [p]


def enrichment_rows(tracks: dict[str, list[dict[str, object]]], masks: dict[str, list[tuple[int, int]]], permutations: int) -> list[dict[str, object]]:
    out = []
    genome_span = CHR_END - CHR_START + 1
    for set_name in SETS:
        for clade in CLADES:
            competitor = strongest_competitor(tracks[clade], SPECIES_TOPOLOGY[clade])
            observed, packed = circular_enrichment(tracks[clade], masks[set_name], competitor, permutations, SEED + stable_offset(set_name, clade))
            null, p = packed[:-1], packed[-1]
            null_mean = mean(null)
            out.append({
                "structural_set": set_name, "clade": clade, "species_topology": SPECIES_TOPOLOGY[clade],
                "competing_topology": competitor, "mask_radius_bp": MASK_RADIUS,
                "structural_span_bp": interval_span(masks[set_name]),
                "fraction_genomic_span": interval_span(masks[set_name]) / genome_span,
                "fraction_windows_associated": observed["fraction_windows_associated"],
                "competing_frequency_associated": observed["assoc_fraction"],
                "competing_frequency_background": observed["background_fraction"],
                "observed_frequency_difference": observed["difference"],
                "observed_fold_enrichment": observed["fold_enrichment"],
                "fraction_competing_windows_in_associated": observed["fraction_competing_windows_in_associated"],
                "representation_ratio_vs_genomic_span": observed["fraction_competing_windows_in_associated"] / (interval_span(masks[set_name]) / genome_span),
                "null_mean_frequency_difference": null_mean, "null_median_frequency_difference": median(null),
                "null_sd_frequency_difference": math.sqrt(sum((x - null_mean) ** 2 for x in null) / (len(null) - 1)),
                "circular_shift_pvalue_one_sided": p, "permutations": permutations,
                "null_method": "circular shift of complete ordered topology-state track; mask fixed",
                "replication_note": "focal tracks are correlated views, not independent replicates",
            })
    return out


def event_rows(tracks: dict[str, list[dict[str, object]]]) -> list[dict[str, object]]:
    out = []
    events = [r for r in read_tsv(EVENTS) if r["structural_set"] == "primary" and int(r["flank_window_bp"]) == MASK_RADIUS]
    assert len(events) == 16
    for event in events:
        raw_start, raw_end = int(event["start"]), int(event["end"])
        if raw_start == raw_end:
            start, end = max(CHR_START, raw_start - MASK_RADIUS), min(CHR_END, raw_end + MASK_RADIUS)
            neighborhood = "single breakpoint +/-250 kb"
        else:
            start, end = raw_start, raw_end
            neighborhood = "frozen paired-event interval"
        for clade in CLADES:
            local = [r for r in tracks[clade] if start <= float(r["midpoint"]) <= end]
            s = topology_summary(local, SPECIES_TOPOLOGY[clade])
            competitor = strongest_competitor(tracks[clade], SPECIES_TOPOLOGY[clade])
            out.append({
                "event_id": event["event_id"], "clade": clade, "event_start": raw_start, "event_end": raw_end,
                "analysis_start": start, "analysis_end": end, "breakpoint_neighborhood": neighborhood,
                "structural_classification": event["event_type"], "pairing_confidence": event["pairing_confidence"],
                "species_support": event["species_support"], "species_topology": SPECIES_TOPOLOGY[clade],
                "competing_topology": competitor, "n_windows": s["n"],
                "local_species_topology_count": s["counts"][SPECIES_TOPOLOGY[clade]],
                "local_species_topology_frequency": s["freqs"][SPECIES_TOPOLOGY[clade]],
                "local_competing_topology_count": s["counts"][competitor],
                "local_competing_topology_frequency": s["freqs"][competitor],
                **{f"mean_support_{q}": s["mean_support"][q] for q in TOPOLOGIES},
                "dominant_topology": s["dominant"], "topology_margin_species_vs_competing": s["freqs"][SPECIES_TOPOLOGY[clade]] - s["freqs"][competitor] if local else math.nan,
            })
    return out


def make_figures(tracks, masks, support_rows, wb_rows) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    subprocess.run(["Rscript", str(Path(__file__).with_name("04a_make_figures.R"))], cwd=ROOT, check=True)


def report(masks, support_rows, wb_rows, enrich_rows, ev_rows, permutations: int) -> None:
    lines = [
        "# Neoaves chromosome 4 Stage 4A: structural topology support", "",
        "## Scope", "",
        "This analysis asks whether frozen rearrangement-associated regions disproportionately contribute support for competing quartet relationships. It measures spatial association and does not identify a causal rearrangement effect or demonstrate bias in a species-tree inference.", "",
        "## Frozen inputs and topology reference", "",
        f"- Ordered tracks: `{LOCI.relative_to(ROOT)}` (Columbea, N61, N62; q1/q2/q3 and dominant state).",
        "- Frozen structural inputs: `denovo_breakpoint_sets.tsv`, the three `stage3b_*_structural_segments.tsv` partitions, and `stage3b_v2_structural_events.tsv`.",
        "- Stage 3A provenance: `stage3a_final_manifest.json` and its breakpoint enrichment outputs.",
        "- Stage 3B provenance: `stage3b_v2_final_manifest.json`, summary, event states, and structural-only prediction manifest.",
        "- The reference species topology is q1 for every focal clade: q1=(C1,C2)|(S,O). This follows the frozen quartet role mapping. It is distinct from the chromosome-wide empirically dominant topology.", "",
        "## Structural masks", "",
        "The mask is the union of ±250,000 bp neighborhoods around every breakpoint in a frozen set, clipped to chr4:1-91,310,470. The 250 kb radius is the pre-existing Stage 3B primary flank window. No genealogy data enter mask construction. Stringent, primary, and inclusive are sensitivity sets; background is the complement.", "",
    ]
    for s in SETS:
        lines.append(f"- {s}: {len(masks[s])} merged intervals, {interval_span(masks[s]):,} bp ({interval_span(masks[s])/(CHR_END-CHR_START+1):.3%} of chr4).")
    lines += ["", "## Results", "", "### Structural category frequencies", ""]
    for s in SETS:
        lines.append(f"**{s}**")
        lines.append("")
        lines.append("| clade | associated q1/q2/q3 | background q1/q2/q3 | associated dominant | background dominant |")
        lines.append("|---|---:|---:|---|---|")
        for c in CLADES:
            a = next(r for r in support_rows if r["structural_set"] == s and r["clade"] == c and r["structural_category"] == "structural-associated")
            b = next(r for r in support_rows if r["structural_set"] == s and r["clade"] == c and r["structural_category"] == "structural-background")
            af = "/".join(f"{a[f'frequency_{q}']:.3f}" for q in TOPOLOGIES)
            bf = "/".join(f"{b[f'frequency_{q}']:.3f}" for q in TOPOLOGIES)
            lines.append(f"| {c} | {af} | {bf} | {a['dominant_topology']} | {b['dominant_topology']} |")
        lines.append("")
    lines += ["### Spatial enrichment", "", f"One-sided p-values use {permutations:,} deterministic circular shifts of each complete ordered topology-state track while holding the structural mask fixed. The three tracks are correlated views and are not pooled as independent replicates.", "", "| set | clade | competitor | fold | difference | representation/span ratio | p |", "|---|---|---|---:|---:|---:|---:|"]
    for r in enrich_rows:
        lines.append(f"| {r['structural_set']} | {r['clade']} | {r['competing_topology']} | {r['observed_fold_enrichment']:.3f} | {r['observed_frequency_difference']:.3f} | {r['representation_ratio_vs_genomic_span']:.3f} | {r['circular_shift_pvalue_one_sided']:.4f} |")
    lines += ["", "### Window versus block normalization", "", "| clade | scheme | q1/q2/q3 | dominant | species-best-alt margin | empirical competing-topology fraction |", "|---|---|---:|---|---:|---:|"]
    for r in wb_rows:
        q = "/".join(f"{r[f'q_{i}']:.3f}" for i in (1, 2, 3))
        lines.append(f"| {r['clade']} | {r['weighting_scheme']} | {q} | {r['dominant_topology']} | {r['species_vs_best_alternative_margin']:.3f} | {r['empirical_competing_topology_fraction']:.3f} |")
    changes = [c for c in CLADES if len({r["dominant_topology"] for r in wb_rows if r["clade"] == c}) > 1]
    lines += ["", f"Dominant topology changes under block normalization: {', '.join(changes) if changes else 'none'}.", "", "The competing-topology fractions are descriptive empirical analogues of epsilon_window and epsilon_block. They are not literal estimates of theoretical contamination parameters because windows/runs need not satisfy the model's independence or generative assumptions.", "", "### Event heterogeneity", ""]
    for c in CLADES:
        rows = [r for r in ev_rows if r["clade"] == c]
        vals = [float(r["topology_margin_species_vs_competing"]) for r in rows if not math.isnan(float(r["topology_margin_species_vs_competing"]))]
        top = min(rows, key=lambda r: float(r["topology_margin_species_vs_competing"]) if r["n_windows"] else math.inf)
        lines.append(f"- {c}: event margins range {min(vals):.3f} to {max(vals):.3f}; strongest local tilt toward {top['competing_topology']} occurs at {top['event_id']} (margin {top['topology_margin_species_vs_competing']:.3f}).")
    lines += ["", "No event is treated as a strict 2:2 resolved mechanism unit here. The Stage 3B eligibility result remains zero and unchanged.", "", "## Interpretation and limitations", "",
              "1. Competing-topology enrichment is track-specific. N61 q2 is enriched in every frozen mask (fold 1.37-1.51; p=0.0031-0.0292), whereas Columbea and N62 q2 are depleted rather than enriched. The correlated tracks therefore do not support a chromosome-wide claim that structural neighborhoods generally enrich competing topologies.",
              "2. Dense windows amplify N61 q2: its empirical competing-topology fraction falls from 0.454 by windows to 0.332 by blocks. For Columbea and N62, q2 instead rises after block normalization (0.260 to 0.331 and 0.283 to 0.347), so dense sampling dilutes that competitor in those tracks.",
              "3. Dominance changes for N61 (q2 to q1) and N62 (q1 to q2), but not Columbea (q1 in both summaries).",
              "4. N61's event pattern is broadly distributed: 13 of 16 event units have a negative q1-minus-q2 margin. The comparable counts are 3 of 16 for Columbea and 4 of 16 for N62, with substantial event-to-event heterogeneity in every track.",
              "5. Direction and inference are stable across stringent, primary, and inclusive masks: N61 is enriched in all three; Columbea and N62 are depleted in all three.", "",
              "Enrichment is assessed separately for each correlated focal track and across frozen mask sensitivities. Evidence that remains similar after block normalization is less attributable to dense runs of repeated windows; changes in fractions quantify that amplification directly. Event-level ranges show whether the pattern is concentrated or distributed, but overlapping neighborhoods and shared genealogy data prevent treating events or clades as fully independent biological replicates.", "",
              "These results establish association only. They do not show that rearrangements caused a topology shift, estimate a literal MSRC failure parameter, or establish an effect on an actual species-tree analysis relative to independent evidence.", "", "## Figures", "",
              "- Figure A: chromosome position, primary structural neighborhoods, and topology state/support.",
              "- Figure B: primary associated versus background topology frequencies.",
              "- Figure C: window-weighted versus block-normalized support.", ""]
    (RESULTS / "stage4a_report.md").write_text("\n".join(lines))


def validate_stage3b_unchanged() -> None:
    rows = read_tsv(STAGE3B_SUMMARY)
    strict = [r for r in rows if r["analysis"].startswith("strict_event_state")]
    assert strict and all(int(r["n_eligible_events"]) == 0 for r in strict)
    assert all("underpowered" in r["status"] for r in strict)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10_000)
    args = parser.parse_args()
    validate_stage3b_unchanged()
    bp_rows = read_tsv(BREAKPOINTS)
    masks = build_masks(bp_rows)
    tracks = load_tracks()
    support = structural_support_rows(tracks, masks)
    wb = window_block_rows(tracks)
    enrich = enrichment_rows(tracks, masks, args.permutations)
    events = event_rows(tracks)
    write_tsv(RESULTS / "stage4a_structural_topology_support.tsv", support, list(support[0]))
    write_tsv(RESULTS / "stage4a_window_vs_block_support.tsv", wb, list(wb[0]))
    write_tsv(RESULTS / "stage4a_event_level_summary.tsv", events, list(events[0]))
    write_tsv(RESULTS / "stage4a_structural_enrichment.tsv", enrich, list(enrich[0]))
    make_figures(tracks, masks, support, wb)
    report(masks, support, wb, enrich, events, args.permutations)
    manifest = {
        "stage": "4A", "seed": SEED, "permutations": args.permutations,
        "mask_rule": "union of +/-250 kb neighborhoods around frozen breakpoints",
        "species_topology": SPECIES_TOPOLOGY,
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [LOCI, BREAKPOINTS, EVENTS, STAGE3B_SUMMARY, *SEGMENTS.values()]},
        "stage3b_eligibility_relaxed": False,
    }
    (RESULTS / "stage4a_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
