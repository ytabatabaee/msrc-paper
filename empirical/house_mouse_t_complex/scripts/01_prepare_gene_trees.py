#!/usr/bin/env python3
"""Stage 1 freeze of source-backed house-mouse t-complex ML gene trees."""

from __future__ import annotations

import argparse
import io
import json
import sys
import zipfile
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import (  # noqa: E402
    DATA_ROOT,
    METADATA,
    PRIMARY_ARCHIVE,
    PROCESSED,
    RESULTS,
    TreeRecord,
    find_zip_member,
    git_commit,
    load_mapping,
    make_tree_records,
    mapping_gate,
    nested_zip_bytes,
    read_tsv,
    sha256_path,
    stage1_can_run,
    write_json,
    write_tsv,
)


TREE_OUT = PROCESSED / "house_mouse_t_complex_ml_5kb.tre"
META_OUT = PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv"
SUBSPECIES_MAP_OUT = PROCESSED / "house_mouse_t_complex_astral_subspecies_map.tsv"
POPULATION_MAP_OUT = PROCESSED / "house_mouse_t_complex_astral_population_map.tsv"
SUBSPECIES_ASTRAL_MAP_OUT = PROCESSED / "house_mouse_t_complex_astral_subspecies.map"
POPULATION_ASTRAL_MAP_OUT = PROCESSED / "house_mouse_t_complex_astral_population.map"
STANDARD_SUBSPECIES_ASTRAL_MAP_OUT = PROCESSED / "house_mouse_t_complex_astral_standard_only_subspecies.map"
SUBSET_OUT = PROCESSED / "house_mouse_t_complex_tip_subsets.tsv"
TREATMENT_DIR = PROCESSED / "stage1_aster_inputs"
TREATMENT_MANIFEST = PROCESSED / "house_mouse_t_complex_stage1_treatments.tsv"


def astral_taxon_label(value: str) -> str:
    return value.strip().replace(" ", "_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Path to IST-2017-78-v1+1_Data.zip")
    parser.add_argument(
        "--allow-unresolved-stage0",
        action="store_true",
        help="Prepare no-op blocked reports even when the Stage-0 mapping gate has not passed.",
    )
    return parser.parse_args()


def order_records(records: list[TreeRecord]) -> list[TreeRecord]:
    return sorted(records, key=lambda r: ((r.start_bp is None, r.start_bp or 10**30), r.source_filename))


def missing_tips(record: TreeRecord, all_tips: set[str]) -> str:
    missing = sorted(all_tips - set(record.tips))
    return ",".join(missing)


def write_mapping_files(mapping_rows: list[dict[str, str]], tips: set[str]) -> dict[str, object]:
    gate = mapping_gate(mapping_rows, tips)
    exact_rows = [
        row
        for row in mapping_rows
        if row.get("tree_tip") in tips
        and row.get("mapping_confidence") in {"exact", "strong"}
        and row.get("subspecies") not in {"", "NA", None}
    ]
    subspecies_rows = [
        {"individual": row["tree_tip"], "species": row["subspecies"], "mapping_confidence": row["mapping_confidence"]}
        for row in sorted(exact_rows, key=lambda x: x["tree_tip"])
    ]
    write_tsv(SUBSPECIES_MAP_OUT, subspecies_rows, ["individual", "species", "mapping_confidence"])
    with SUBSPECIES_ASTRAL_MAP_OUT.open("w") as handle:
        for row in subspecies_rows:
            handle.write(f"{row['individual']}\t{astral_taxon_label(row['species'])}\n")
    standard_rows = [
        row
        for row in exact_rows
        if row.get("t_status") in {"standard_noncarrier", "outgroup_not_t_haplotype"}
    ]
    with STANDARD_SUBSPECIES_ASTRAL_MAP_OUT.open("w") as handle:
        for row in sorted(standard_rows, key=lambda x: x["tree_tip"]):
            handle.write(f"{row['tree_tip']}\t{astral_taxon_label(row['subspecies'])}\n")
    pop_rows = [
        {
            "individual": row["tree_tip"],
            "population": row.get("population", "NA"),
            "subspecies": row.get("subspecies", "NA"),
            "mapping_confidence": row.get("mapping_confidence", "unresolved"),
        }
        for row in sorted(exact_rows, key=lambda x: x["tree_tip"])
        if row.get("population") not in {"", "NA", None}
    ]
    write_tsv(POPULATION_MAP_OUT, pop_rows, ["individual", "population", "subspecies", "mapping_confidence"])
    with POPULATION_ASTRAL_MAP_OUT.open("w") as handle:
        for row in pop_rows:
            handle.write(f"{row['individual']}\t{astral_taxon_label(row['population'])}\n")
    subset_rows = []
    for treatment in ("all_available", "standard_only", "t_haplotype_only"):
        for row in sorted(exact_rows, key=lambda x: x["tree_tip"]):
            t_status = row.get("t_status", "NA")
            include = treatment == "all_available"
            if treatment == "standard_only":
                include = t_status in {"standard", "noncarrier", "non-t", "standard/noncarrier"}
            if treatment == "t_haplotype_only":
                include = t_status in {"t-haplotype", "t", "pseudo-t", "pseudo_t"}
            if include:
                subset_rows.append(
                    {
                        "subset": treatment,
                        "tree_tip": row["tree_tip"],
                        "individual": row.get("individual", "NA"),
                        "subspecies": row.get("subspecies", "NA"),
                        "population": row.get("population", "NA"),
                        "t_status": t_status,
                        "mapping_confidence": row.get("mapping_confidence", "unresolved"),
                    }
                )
    write_tsv(
        SUBSET_OUT,
        subset_rows,
        ["subset", "tree_tip", "individual", "subspecies", "population", "t_status", "mapping_confidence"],
    )
    return {
        "mapping_gate": gate,
        "subspecies_mapping_rows": len(subspecies_rows),
        "standard_subspecies_mapping_rows": len(standard_rows),
        "population_mapping_rows": len(pop_rows),
        "subspecies_astral_map": str(SUBSPECIES_ASTRAL_MAP_OUT),
        "standard_subspecies_astral_map": str(STANDARD_SUBSPECIES_ASTRAL_MAP_OUT),
        "population_astral_map": str(POPULATION_ASTRAL_MAP_OUT),
        "subset_rows": len(subset_rows),
    }


def prune_newick_to_tips(newick: str, keep: set[str]) -> str:
    from Bio import Phylo

    tree = Phylo.read(io.StringIO(newick), "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep:
            tree.prune(terminal)
    out = io.StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip()


def write_treatment_inputs(records: list[TreeRecord], mapping_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    TREATMENT_DIR.mkdir(parents=True, exist_ok=True)
    by_tip = {row["tree_tip"]: row for row in mapping_rows}
    all_mapped = {
        tip
        for tip, row in by_tip.items()
        if row.get("mapping_confidence") in {"exact", "strong"} and row.get("subspecies") not in {"", "NA", None}
    }
    standard_tips = {
        tip
        for tip in all_mapped
        if by_tip[tip].get("t_status") in {"standard_noncarrier", "outgroup_not_t_haplotype"}
    }
    rows: list[dict[str, object]] = []

    def write_tree_file(name: str, subset: list[TreeRecord], keep: set[str] | None, notes: str) -> None:
        path = TREATMENT_DIR / f"{name}.tre"
        with path.open("w") as handle:
            for record in subset:
                if keep is None:
                    handle.write(record.newick.rstrip() + "\n")
                else:
                    handle.write(prune_newick_to_tips(record.newick, keep) + "\n")
        groups = sorted({by_tip[t].get("subspecies", "") for t in (keep or all_mapped) if by_tip.get(t, {}).get("subspecies")})
        rows.append(
            {
                "treatment": name,
                "tree_file": str(path),
                "n_trees": len(subset),
                "n_tips": len(keep or all_mapped),
                "n_subspecies_groups": len(groups),
                "subspecies_groups": ",".join(groups),
                "valid_for_astral": len(groups) >= 4,
                "notes": notes,
            }
        )

    ordered = order_records([r for r in records if r.tree_valid])
    write_tree_file("T1_ALL_WINDOWS", ordered, None, "All usable published 5-kb ML windows and all mapped tips.")
    write_tree_file(
        "T0_STANDARD",
        ordered,
        standard_tips,
        "Pseudo-t/t-haplotype tips pruned; standard/noncarrier M. musculus and M. spretus outgroup retained.",
    )
    for interval in (50_000, 100_000):
        bins: dict[int, list[TreeRecord]] = defaultdict(list)
        for record in ordered:
            if record.start_bp is not None:
                bins[record.start_bp // interval].append(record)
        thinned = [sorted(bin_records, key=lambda r: (r.start_bp or 0, r.source_filename))[0] for _bin, bin_records in sorted(bins.items())]
        write_tree_file(
            f"T2_SPATIALLY_THINNED_{interval // 1000}kb",
            thinned,
            None,
            f"One deterministic earliest-start tree per {interval // 1000} kb genomic bin.",
        )
    write_tsv(
        TREATMENT_MANIFEST,
        rows,
        ["treatment", "tree_file", "n_trees", "n_tips", "n_subspecies_groups", "subspecies_groups", "valid_for_astral", "notes"],
    )
    return rows


def blocked_report(reason: str, archive: Path) -> None:
    manifest = {
        "analysis": "house_mouse_t_complex_stage1_gene_tree_freeze",
        "status": "blocked",
        "blocker": reason,
        "archive": str(archive),
        "git_commit": git_commit(),
    }
    write_json(RESULTS / "stage1_manifest.json", manifest)
    (RESULTS / "stage1_report.md").write_text(
        "# House mouse t-complex Stage 1 report\n\n"
        f"Stage 1 was not run because {reason}\n\n"
        "This is an intentional scientific gate. Do not manufacture a species/subspecies mapping.\n"
    )


def main() -> int:
    args = parse_args()
    ok, reason, stage0 = stage1_can_run()
    if not ok and not args.allow_unresolved_stage0:
        blocked_report(reason, args.archive)
        print(f"Stage 1 blocked: {reason}")
        return 0
    with zipfile.ZipFile(args.archive) as outer:
        primary_data, info = nested_zip_bytes(outer, PRIMARY_ARCHIVE)
        if primary_data is None or info is None:
            raise SystemExit(f"Primary ML tree archive not found: {PRIMARY_ARCHIVE}")
    records, coord = make_tree_records(primary_data)
    valid = [r for r in order_records(records) if r.tree_valid]
    if any(r.start_bp is None or r.end_bp is None for r in valid):
        raise SystemExit("Cannot freeze Stage 1 because not all valid trees have recovered coordinates.")
    all_tips = sorted({tip for r in valid for tip in r.tips})
    mapping_rows = load_mapping(METADATA / "tip_mapping.tsv")
    mapping_summary = write_mapping_files(mapping_rows, set(all_tips))
    treatment_rows = write_treatment_inputs(valid, mapping_rows)
    with TREE_OUT.open("w") as handle:
        for record in valid:
            handle.write(record.newick.rstrip() + "\n")
    meta_rows = [
        {
            "locus_id": record.window_id,
            "start_bp": record.start_bp,
            "end_bp": record.end_bp,
            "midpoint_bp": record.midpoint_bp,
            "source_filename": record.source_filename,
            "n_tips": record.n_tips,
            "missing_tips": missing_tips(record, set(all_tips)),
            "tree_sha256": record.sha256,
        }
        for record in valid
    ]
    write_tsv(
        META_OUT,
        meta_rows,
        ["locus_id", "start_bp", "end_bp", "midpoint_bp", "source_filename", "n_tips", "missing_tips", "tree_sha256"],
    )
    outputs = {
        "gene_tree_file": str(TREE_OUT.relative_to(DATA_ROOT.parents[1])),
        "metadata_table": str(META_OUT.relative_to(DATA_ROOT.parents[1])),
        "subspecies_mapping": str(SUBSPECIES_MAP_OUT.relative_to(DATA_ROOT.parents[1])),
        "subspecies_astral_map": str(SUBSPECIES_ASTRAL_MAP_OUT.relative_to(DATA_ROOT.parents[1])),
        "standard_subspecies_astral_map": str(STANDARD_SUBSPECIES_ASTRAL_MAP_OUT.relative_to(DATA_ROOT.parents[1])),
        "population_mapping": str(POPULATION_MAP_OUT.relative_to(DATA_ROOT.parents[1])),
        "population_astral_map": str(POPULATION_ASTRAL_MAP_OUT.relative_to(DATA_ROOT.parents[1])),
        "tip_subsets": str(SUBSET_OUT.relative_to(DATA_ROOT.parents[1])),
        "treatment_manifest": str(TREATMENT_MANIFEST.relative_to(DATA_ROOT.parents[1])),
    }
    manifest = {
        "analysis": "house_mouse_t_complex_stage1_gene_tree_freeze",
        "status": "prepared",
        "command": "01_prepare_gene_trees.py --archive " + str(args.archive),
        "git_commit": git_commit(),
        "source_archive": str(args.archive),
        "source_archive_sha256": sha256_path(args.archive),
        "n_valid_trees_written": len(valid),
        "coordinate_summary": coord,
        "mapping_summary": mapping_summary,
        "treatments": treatment_rows,
        "outputs": outputs,
    }
    write_json(RESULTS / "stage1_manifest.json", manifest)
    (RESULTS / "stage1_report.md").write_text(
        "# House mouse t-complex Stage 1 report\n\n"
        f"Status: prepared\n\nUsable ML trees written: {len(valid)}\n\n"
        f"Coordinates recovered for all valid trees: {coord.get('coordinates_recoverable_for_all_final_trees')}\n\n"
        "ASTRAL mapping files use only exact/strong source-backed mappings from `metadata/tip_mapping.tsv`.\n"
        "Tentative and unresolved mappings are excluded automatically.\n\n"
        "Limitation: these windows represent the published chr17 t-complex interval, not a genome-wide background.\n"
    )
    print(f"Wrote {TREE_OUT}")
    print(f"Valid ML trees: {len(valid)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
