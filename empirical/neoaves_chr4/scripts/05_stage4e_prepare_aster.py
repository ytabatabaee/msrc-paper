#!/usr/bin/env python3
"""Prepare Stage 4E ASTER/ASTRAL4 inputs without launching inference."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage4d_core import (  # noqa: E402
    DATA as STAGE4D_DATA,
    RESULTS,
    ROOT,
    is_chr4,
    parse_newick,
    prune_tree,
    emit_newick,
    representatives,
    sha256,
    read_tsv,
    write_tsv,
)

BASE = ROOT / "empirical/neoaves_chr4"
DATA = BASE / "data/stage4e"
CLUSTER = BASE / "cluster/stage4e_aster"
RESULT_TREES = BASE / "results/stage4e/trees"
SOURCE = STAGE4D_DATA / "named_recollapsed.tre"
LOCUS_MANIFEST = RESULTS / "stage4d_locus_manifest.tsv"
J48_TAXA = STAGE4D_DATA / "jarvis48_taxa.tsv"
MASTER_RANDOM_SEED = 4202404

EXPECTED = {
    "T0": 63430,
    "NONCHR4": 57168,
    "PNAS_EXCLUDED": 1431,
    "STRUCT_STRINGENT_EXCLUDED": 495,
    "STRUCT_PRIMARY_EXCLUDED": 835,
    "STRUCT_INCLUSIVE_EXCLUDED": 1308,
}

FULL_TREATMENTS = [
    "T0",
    "T_PNAS",
    "T_STRUCT_PRIMARY",
    "T_STRUCT_STRINGENT",
    "T_STRUCT_INCLUSIVE",
    "T_BLOCK_250K",
    "T_BLOCK_500K",
    "T_BLOCK_1MB",
    "NONCHR4",
]

J48_TREATMENTS = [
    "J48_T0",
    "J48_T_PNAS",
    "J48_T_STRUCT_PRIMARY",
    "J48_T_BLOCK_250K",
    "J48_T_BLOCK_500K",
    "J48_T_BLOCK_1MB",
    "J48_NONCHR4",
]

PRIMARY_CLUSTER_TREATMENTS = ["T0", "T_PNAS", "T_STRUCT_PRIMARY", "T_BLOCK_500K", "NONCHR4"]
SENSITIVITY_CLUSTER_TREATMENTS = [
    "T_STRUCT_STRINGENT",
    "T_STRUCT_INCLUSIVE",
    "T_BLOCK_250K",
    "T_BLOCK_1MB",
]


def boolish(value: object) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def mkdirs() -> None:
    for path in (DATA, DATA / "random_matched", CLUSTER, RESULT_TREES, BASE / "results/stage4e"):
        path.mkdir(parents=True, exist_ok=True)


def treatment_id_sets(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    nonchr4 = {r["locus_id"] for r in rows if not boolish(r["is_chr4"])}
    sets = {
        "T0": {r["locus_id"] for r in rows},
        "NONCHR4": nonchr4,
        "T_PNAS": {r["locus_id"] for r in rows if not boolish(r["published_outlier_region"])},
        "T_STRUCT_PRIMARY": {r["locus_id"] for r in rows if not boolish(r["frozen_structural_primary"])},
        "T_STRUCT_STRINGENT": {r["locus_id"] for r in rows if not boolish(r["frozen_structural_stringent"])},
        "T_STRUCT_INCLUSIVE": {r["locus_id"] for r in rows if not boolish(r["frozen_structural_inclusive"])},
    }
    for label, width in (("T_BLOCK_250K", 250000), ("T_BLOCK_500K", 500000), ("T_BLOCK_1MB", 1000000)):
        reps = set(representatives(rows, width).values())
        sets[label] = nonchr4 | reps
    return sets


def treatment_description(name: str) -> str:
    return {
        "T0": "all 63,430 ID-preserving re-collapsed loci",
        "T_PNAS": "published PNAS chr4 outlier intervals excluded",
        "T_STRUCT_PRIMARY": "primary frozen structural mask excluded",
        "T_STRUCT_STRINGENT": "stringent frozen structural mask excluded",
        "T_STRUCT_INCLUSIVE": "inclusive frozen structural mask excluded",
        "T_BLOCK_250K": "non-chr4 plus deterministic 250-kb chr4 bin representatives",
        "T_BLOCK_500K": "non-chr4 plus deterministic 500-kb chr4 bin representatives",
        "T_BLOCK_1MB": "non-chr4 plus deterministic 1-Mb chr4 bin representatives",
        "NONCHR4": "all chr4 and chr4_* loci excluded",
    }[name]


def validate_manifest(rows: list[dict[str, str]]) -> dict[str, object]:
    if len(rows) != EXPECTED["T0"]:
        raise RuntimeError(f"Expected 63,430 loci, observed {len(rows)}")
    ids = [r["locus_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise RuntimeError("Stage 4D locus IDs are not unique")
    counts = {
        "n_loci": len(rows),
        "n_unique_locus_ids": len(set(ids)),
        "n_chr4_loci": sum(boolish(r["is_chr4"]) for r in rows),
        "n_nonchr4_loci": sum(not boolish(r["is_chr4"]) for r in rows),
        "n_pnas_outlier_loci": sum(boolish(r["published_outlier_region"]) for r in rows),
        "n_struct_stringent_loci": sum(boolish(r["frozen_structural_stringent"]) for r in rows),
        "n_struct_primary_loci": sum(boolish(r["frozen_structural_primary"]) for r in rows),
        "n_struct_inclusive_loci": sum(boolish(r["frozen_structural_inclusive"]) for r in rows),
    }
    assert counts["n_nonchr4_loci"] == EXPECTED["NONCHR4"], counts
    assert counts["n_pnas_outlier_loci"] == EXPECTED["PNAS_EXCLUDED"], counts
    assert counts["n_struct_stringent_loci"] == EXPECTED["STRUCT_STRINGENT_EXCLUDED"], counts
    assert counts["n_struct_primary_loci"] == EXPECTED["STRUCT_PRIMARY_EXCLUDED"], counts
    assert counts["n_struct_inclusive_loci"] == EXPECTED["STRUCT_INCLUSIVE_EXCLUDED"], counts
    if any(boolish(r["is_chr4"]) != is_chr4(r["chromosome"]) for r in rows):
        raise RuntimeError("Frozen chr4 group disagrees with chromosome-group rule")
    return counts


def verify_source_order(rows: list[dict[str, str]]) -> dict[str, object]:
    n_lines = 0
    with SOURCE.open() as handle:
        for n_lines, (row, line) in enumerate(zip(rows, handle, strict=True), start=1):
            if not line.strip().endswith(";"):
                raise RuntimeError(f"Source tree line {n_lines} is not Newick")
            if n_lines != int(row["tree_index"]):
                raise RuntimeError(f"Locus manifest order mismatch at {row['locus_id']}")
    if n_lines != len(rows):
        raise RuntimeError("Source tree count differs from locus manifest")
    known = BASE / "results/stage4d/j48_existing_input_audit.json"
    known_sha = None
    if known.exists():
        payload = json.loads(known.read_text())
        known_sha = payload.get("source_sha256", {}).get("named_recollapsed.tre")
    source_sha = sha256(SOURCE)
    if known_sha and source_sha != known_sha:
        raise RuntimeError("named_recollapsed.tre checksum differs from Stage 4D audit")
    return {"source": str(SOURCE.relative_to(ROOT)), "source_sha256": source_sha, "source_lines": n_lines}


def write_locus_id_file(path: Path, ids: set[str], order: list[str]) -> None:
    with path.open("w") as handle:
        for locus_id in order:
            if locus_id in ids:
                handle.write(locus_id + "\n")


def random_seed(regime: str, replicate: int) -> int:
    seed_material = f"{MASTER_RANDOM_SEED}:{regime}:{replicate}".encode()
    return int.from_bytes(hashlib.sha256(seed_material).digest()[:8], "big")


def random_replicate_ids(rows: list[dict[str, str]], target_chr4: int, seed: int) -> set[str]:
    nonchr4 = {r["locus_id"] for r in rows if not boolish(r["is_chr4"])}
    chr4 = sorted(r["locus_id"] for r in rows if boolish(r["is_chr4"]))
    return nonchr4 | set(random.Random(seed).sample(chr4, target_chr4))


def write_random_manifests(rows: list[dict[str, str]], sets: dict[str, set[str]], order: list[str]) -> list[dict[str, object]]:
    regimes = {
        "T_RANDOM_MATCHED_PNAS": "T_PNAS",
        "T_RANDOM_MATCHED_STRUCT": "T_STRUCT_PRIMARY",
        "T_RANDOM_MATCHED_BLOCK250": "T_BLOCK_250K",
        "T_RANDOM_MATCHED_BLOCK500": "T_BLOCK_500K",
        "T_RANDOM_MATCHED_BLOCK1M": "T_BLOCK_1MB",
    }
    nonchr4 = {r["locus_id"] for r in rows if not boolish(r["is_chr4"])}
    chr4 = sorted(r["locus_id"] for r in rows if boolish(r["is_chr4"]))
    random_rows: list[dict[str, object]] = []
    index_rows: list[dict[str, object]] = []
    outdir = DATA / "random_matched"
    for regime, treatment in regimes.items():
        target_chr4 = len(sets[treatment] - nonchr4)
        for replicate in range(1, 101):
            seed = random_seed(regime, replicate)
            ids = random_replicate_ids(rows, target_chr4, seed)
            path = outdir / f"{regime}_rep{replicate:03d}.locus_ids"
            write_locus_id_file(path, ids, order)
            row = {
                "regime": regime,
                "replicate": replicate,
                "seed": seed,
                "target_treatment": treatment,
                "n_loci": len(ids),
                "n_chr4_loci": target_chr4,
                "relative_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            random_rows.append(row)
            index_rows.append(row)
    write_tsv(RESULTS / "stage4e_random_matched_manifest.tsv", random_rows)
    return index_rows


def write_treatment_inputs(rows: list[dict[str, str]], sets: dict[str, set[str]]) -> list[dict[str, object]]:
    handles = {name: (DATA / f"{name}.tre").open("w") for name in FULL_TREATMENTS}
    counts = {name: 0 for name in FULL_TREATMENTS}
    try:
        with SOURCE.open() as source:
            for row, line in zip(rows, source, strict=True):
                locus_id = row["locus_id"]
                for name, ids in sets.items():
                    if name in handles and locus_id in ids:
                        handles[name].write(line)
                        counts[name] += 1
    finally:
        for handle in handles.values():
            handle.close()
    return [
        {
            "treatment": name,
            "scope": "full363",
            "relative_path": str((DATA / f"{name}.tre").relative_to(ROOT)),
            "n_loci": counts[name],
            "n_chr4_loci": sum(1 for r in rows if boolish(r["is_chr4"]) and r["locus_id"] in sets[name]),
            "sha256": sha256(DATA / f"{name}.tre"),
            "description": treatment_description(name),
        }
        for name in FULL_TREATMENTS
    ]


def write_j48_inputs(rows: list[dict[str, str]], sets: dict[str, set[str]]) -> list[dict[str, object]]:
    taxa = {r["stiller2024_taxon"] for r in read_tsv(J48_TAXA)}
    if len(taxa) != 48:
        raise RuntimeError(f"Expected 48 Jarvis taxa, observed {len(taxa)}")
    mapping = {
        "J48_T0": "T0",
        "J48_T_PNAS": "T_PNAS",
        "J48_T_STRUCT_PRIMARY": "T_STRUCT_PRIMARY",
        "J48_T_BLOCK_250K": "T_BLOCK_250K",
        "J48_T_BLOCK_500K": "T_BLOCK_500K",
        "J48_T_BLOCK_1MB": "T_BLOCK_1MB",
        "J48_NONCHR4": "NONCHR4",
    }
    handles = {name: (DATA / f"{name}.tre").open("w") for name in J48_TREATMENTS}
    counts = {name: 0 for name in J48_TREATMENTS}
    omitted = {name: 0 for name in J48_TREATMENTS}
    observed: set[str] = set()
    try:
        with SOURCE.open() as source:
            for row, line in zip(rows, source, strict=True):
                tree = prune_tree(parse_newick(line), taxa)
                leaves = set()
                if tree is not None:
                    def collect(node):
                        if not node[0]:
                            leaves.add(node[1])
                        for child in node[0]:
                            collect(child)
                    collect(tree)
                    observed.update(leaves)
                for name, full_name in mapping.items():
                    if row["locus_id"] not in sets[full_name]:
                        continue
                    if tree is None or len(leaves) < 4:
                        omitted[name] += 1
                        continue
                    handles[name].write(emit_newick(tree) + "\n")
                    counts[name] += 1
    finally:
        for handle in handles.values():
            handle.close()
    if observed != taxa:
        raise RuntimeError("Jarvis-48 mapping did not recover all target taxa")
    return [
        {
            "treatment": name,
            "scope": "jarvis48",
            "relative_path": str((DATA / f"{name}.tre").relative_to(ROOT)),
            "n_loci": counts[name],
            "n_chr4_loci": "NA",
            "trees_with_fewer_than_four_taxa_omitted": omitted[name],
            "sha256": sha256(DATA / f"{name}.tre"),
            "description": f"Jarvis-48 pruning of {mapping[name]}",
        }
        for name in J48_TREATMENTS
    ]


def write_preflight(input_rows: list[dict[str, object]], source_info: dict[str, object], counts: dict[str, object]) -> None:
    payload = {
        "stage": "Stage 4E - ASTER/ASTRAL4 full-tree treatment analysis",
        "status": "prepared_not_run",
        "large_inference_run_locally": False,
        "source": source_info,
        "frozen_manifest": str(LOCUS_MANIFEST.relative_to(ROOT)),
        "manifest_sha256": sha256(LOCUS_MANIFEST),
        "counts": counts,
        "inputs": input_rows,
        "note": "No ASTER/ASTRAL4 inference has been run by this preparation script.",
    }
    (RESULTS / "stage4e_preflight.json").write_text(json.dumps(payload, indent=2) + "\n")


def write_input_manifest(input_rows: list[dict[str, object]]) -> None:
    fields = [
        "treatment",
        "scope",
        "relative_path",
        "n_loci",
        "n_chr4_loci",
        "trees_with_fewer_than_four_taxa_omitted",
        "sha256",
        "description",
    ]
    normalized = []
    for row in input_rows:
        normalized.append({field: row.get(field, "") for field in fields})
    write_tsv(RESULTS / "stage4e_input_manifest.tsv", normalized, fields)


def write_transfer_manifest(input_rows: list[dict[str, object]], random_rows: list[dict[str, object]]) -> None:
    rows = []
    for row in input_rows + random_rows:
        path = ROOT / row["relative_path"]
        rows.append(
            {
                "relative_path": row["relative_path"],
                "byte_size": path.stat().st_size,
                "sha256": row["sha256"],
                "required_on_cluster": "yes",
                "description": row.get("description", f"{row.get('regime', 'random')} replicate locus IDs"),
            }
        )
    write_tsv(RESULTS / "stage4e_transfer_manifest.tsv", rows)


def write_analysis_plan() -> None:
    text = """# Stage 4E analysis plan - ASTER/ASTRAL4 full-tree treatment analysis

Stage 4E is a new ASTER/ASTRAL4 stage. Historical Stage 4D ASTRAL-MP
reproduction attempts remain preserved as provenance and are not reinterpreted
as successful gates. ASTRAL-MP reproduction was discontinued as a prerequisite
for the new analysis. New full-tree treatments use ASTER/ASTRAL4 consistently.

The prepared full-363 treatments are T0, T_PNAS, T_STRUCT_PRIMARY,
T_STRUCT_STRINGENT, T_STRUCT_INCLUSIVE, T_BLOCK_250K, T_BLOCK_500K,
T_BLOCK_1MB, and NONCHR4. Primary execution on Bridges-2 should run T0,
T_PNAS, T_STRUCT_PRIMARY, T_BLOCK_500K, and NONCHR4 first.

Random-matched controls are predeclared as 100 replicate locus-ID manifests for
PNAS, STRUCT_PRIMARY, BLOCK250, BLOCK500, and BLOCK1M. They retain all non-chr4
loci and sample the same number of chr4 loci as the corresponding deterministic
treatment using deterministic seeds.

Jarvis-48 inputs are generated from the same ID-preserving Stage 4E source.
Trees with fewer than four retained taxa are omitted and counted.

Predeclared interpretation categories, to be assigned only after cluster output
validation and random-matched comparison, are:

- FULL-TOPOLOGY CORRECTION
- FOCAL-BRANCH CORRECTION
- SUPPORT/BRANCH-LENGTH EFFECT
- GENERIC-THINNING EFFECT
- NO FULL-TREE EFFECT

ASTER T0 versus the published 63K tree is a historical sanity comparison, not a
gate. ASTER T0 remains the baseline for Stage 4E treatment comparisons even if
RF(ASTRAL4_T0, published_63K) is greater than zero.

No Stage 4E full-tree inference is complete until Bridges-2 outputs are
returned, collected, and validated.
"""
    (RESULTS / "stage4e_analysis_plan.md").write_text(text)


def slurm_header(job_name: str, array: str | None = None) -> str:
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=32",
        "#SBATCH --time=24:00:00",
        "#SBATCH --mem=64G",
        "# TODO/user: set --account and --partition if required by PSC allocation",
        "# TODO/user: adjust cpus-per-task, time, and mem after Bridges-2 preflight if needed",
    ]
    if array:
        lines.append(f"#SBATCH --array={array}")
    return "\n".join(lines) + "\n\n"


def write_slurm_scripts() -> None:
    run_py = "${MSRC_REPO}/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py"
    setup = """set -euo pipefail

: "${MSRC_REPO:=$PWD}"
: "${NEOAVES_STAGE4E_DATA:=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e}"
: "${NEOAVES_STAGE4E_OUT:=$MSRC_REPO/empirical/neoaves_chr4/results/stage4e/cluster}"
: "${ASTER_BIN:?Set ASTER_BIN to the Bridges-2 astral4 executable}"
: "${PYTHON:=python3}"

mkdir -p "$NEOAVES_STAGE4E_OUT"
"""
    (CLUSTER / "00_preflight.sbatch").write_text(
        slurm_header("neoaves-aster-preflight")
        + setup
        + f"""
"$PYTHON" "{run_py}" \\
  --preflight \\
  --input-manifest "$MSRC_REPO/empirical/neoaves_chr4/results/stage4e_input_manifest.tsv" \\
  --aster-bin "$ASTER_BIN" \\
  --data-dir "$NEOAVES_STAGE4E_DATA" \\
  --output-root "$NEOAVES_STAGE4E_OUT"
"""
    )
    for filename, job, treatments in (
        ("01_full363_primary.sbatch", "neoaves-aster-primary", PRIMARY_CLUSTER_TREATMENTS),
        ("02_full363_sensitivities.sbatch", "neoaves-aster-sens", SENSITIVITY_CLUSTER_TREATMENTS),
    ):
        body = "\n".join(
            [
                f'TREATMENTS=({" ".join(treatments)})',
                'for treatment in "${TREATMENTS[@]}"; do',
                f'  "$PYTHON" "{run_py}" --execute --treatment "$treatment" --aster-bin "$ASTER_BIN" --input "$NEOAVES_STAGE4E_DATA/$treatment.tre" --output-dir "$NEOAVES_STAGE4E_OUT/$treatment"',
                "done",
                "",
            ]
        )
        (CLUSTER / filename).write_text(slurm_header(job) + setup + body)
    j48 = ["J48_T0", "J48_T_PNAS", "J48_T_STRUCT_PRIMARY", "J48_T_BLOCK_500K", "J48_NONCHR4"]
    (CLUSTER / "03_j48_primary.sbatch").write_text(
        slurm_header("neoaves-aster-j48")
        + setup
        + "\n".join(
            [
                f'TREATMENTS=({" ".join(j48)})',
                'for treatment in "${TREATMENTS[@]}"; do',
                f'  "$PYTHON" "{run_py}" --execute --treatment "$treatment" --aster-bin "$ASTER_BIN" --input "$NEOAVES_STAGE4E_DATA/$treatment.tre" --output-dir "$NEOAVES_STAGE4E_OUT/$treatment"',
                "done",
                "",
            ]
        )
    )
    (CLUSTER / "04_random_matched_array.sbatch").write_text(
        slurm_header("neoaves-aster-random", "1-500")
        + setup
        + """
MANIFEST="$MSRC_REPO/empirical/neoaves_chr4/results/stage4e_random_matched_manifest.tsv"
ROW=$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 1 {print; exit}' "$MANIFEST")
if [[ -z "$ROW" ]]; then
  echo "No random-matched manifest row for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi
REGIME=$(printf '%s\\n' "$ROW" | awk -F '\\t' '{print $1}')
REPLICATE=$(printf '%s\\n' "$ROW" | awk -F '\\t' '{print $2}')
RELATIVE_IDS=$(printf '%s\\n' "$ROW" | awk -F '\\t' '{print $7}')
TREATMENT="${REGIME}_rep$(printf '%03d' "$REPLICATE")"
INPUT="$NEOAVES_STAGE4E_OUT/materialized/${TREATMENT}.tre"
mkdir -p "$(dirname "$INPUT")"
"$PYTHON" "$MSRC_REPO/empirical/neoaves_chr4/scripts/05_stage4e_prepare_aster.py" \\
  --materialize-random "$MSRC_REPO/$RELATIVE_IDS" \\
  --source "$NEOAVES_STAGE4E_DATA/T0.tre" \\
  --output "$INPUT"
"$PYTHON" "$MSRC_REPO/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" \\
  --execute \\
  --treatment "$TREATMENT" \\
  --aster-bin "$ASTER_BIN" \\
  --input "$INPUT" \\
  --output-dir "$NEOAVES_STAGE4E_OUT/random_matched/$TREATMENT"
"""
    )
    (CLUSTER / "README.md").write_text(
        """# Stage 4E ASTER/ASTRAL4 Bridges-2 execution

This directory contains submission scripts only. Do not run ASTER/ASTRAL4 on
the large Neoaves inputs locally.

Suggested cluster layout:

```text
$PROJECT/neoaves_stage4e/
  repo/
  inputs/
  outputs/
  software/
```

Transfer the repository and Stage 4E input files manually. Example:

```bash
rsync -avP empirical/neoaves_chr4/data/stage4e/ PSC_USER@bridges2.psc.edu:$PROJECT/neoaves_stage4e/repo/empirical/neoaves_chr4/data/stage4e/
rsync -avP empirical/neoaves_chr4/results/stage4e_* PSC_USER@bridges2.psc.edu:$PROJECT/neoaves_stage4e/repo/empirical/neoaves_chr4/results/
rsync -avP empirical/neoaves_chr4/cluster/stage4e_aster/ PSC_USER@bridges2.psc.edu:$PROJECT/neoaves_stage4e/repo/empirical/neoaves_chr4/cluster/stage4e_aster/
```

On Bridges-2, set:

```bash
export MSRC_REPO=$PROJECT/neoaves_stage4e/repo
export NEOAVES_STAGE4E_DATA=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e
export NEOAVES_STAGE4E_OUT=$PROJECT/neoaves_stage4e/outputs
export ASTER_BIN=$PROJECT/neoaves_stage4e/software/astral4
export PYTHON=python3
```

Run preflight first:

```bash
sbatch empirical/neoaves_chr4/cluster/stage4e_aster/00_preflight.sbatch
```

Then submit primary jobs manually after preflight passes:

```bash
sbatch empirical/neoaves_chr4/cluster/stage4e_aster/01_full363_primary.sbatch
sbatch empirical/neoaves_chr4/cluster/stage4e_aster/02_full363_sensitivities.sbatch
sbatch empirical/neoaves_chr4/cluster/stage4e_aster/03_j48_primary.sbatch
sbatch empirical/neoaves_chr4/cluster/stage4e_aster/04_random_matched_array.sbatch
```
"""
    )


def materialize_random(locus_ids: Path, source: Path, output: Path) -> None:
    ids = {line.strip() for line in locus_ids.read_text().splitlines() if line.strip()}
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(LOCUS_MANIFEST)
    with source.open() as src, output.open("w") as out:
        for row, line in zip(rows, src, strict=True):
            if row["locus_id"] in ids:
                out.write(line)


def run_internal_tests() -> None:
    rows = [
        {"locus_id": "z", "chromosome": "chr4", "midpoint": "250000", "is_chr4": "True"},
        {"locus_id": "a", "chromosome": "chr4", "midpoint": "250000", "is_chr4": "True"},
        {"locus_id": "b", "chromosome": "chr4", "midpoint": "750000", "is_chr4": "True"},
        {"locus_id": "r", "chromosome": "chr4_random", "midpoint": "250000", "is_chr4": "True"},
    ]
    assert representatives(rows) == {("chr4", 0): "a", ("chr4", 1): "b", ("chr4_random", 0): "r"}
    seed_a = int.from_bytes(hashlib.sha256(b"4202404:x:1").digest()[:8], "big")
    seed_b = int.from_bytes(hashlib.sha256(b"4202404:x:1").digest()[:8], "big")
    assert seed_a == seed_b


def prepare_inputs(write_slurm: bool, run_tests: bool) -> None:
    mkdirs()
    rows = read_tsv(LOCUS_MANIFEST)
    counts = validate_manifest(rows)
    source_info = verify_source_order(rows)
    sets = treatment_id_sets(rows)
    order = [r["locus_id"] for r in rows]
    full_rows = write_treatment_inputs(rows, sets)
    j48_rows = write_j48_inputs(rows, sets)
    random_rows = write_random_manifests(rows, sets, order)
    all_inputs = full_rows + j48_rows
    write_input_manifest(all_inputs)
    write_transfer_manifest(all_inputs, random_rows)
    write_preflight(all_inputs, source_info, counts)
    write_analysis_plan()
    if write_slurm:
        write_slurm_scripts()
    if run_tests:
        run_internal_tests()
    print("PREPARATION COMPLETE.")
    print("No ASTER/ASTRAL4 inference was run.")
    print("Transfer the Stage-4E data directory to Bridges-2 and submit the generated SLURM jobs manually.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-inputs", action="store_true")
    parser.add_argument("--write-slurm", action="store_true")
    parser.add_argument("--run-tests", action="store_true")
    parser.add_argument("--materialize-random", type=Path)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.materialize_random:
        if args.output is None:
            parser.error("--materialize-random requires --output")
        materialize_random(args.materialize_random, args.source, args.output)
        return
    if not args.prepare_inputs and not args.write_slurm and not args.run_tests:
        parser.print_help()
        return
    if args.prepare_inputs:
        prepare_inputs(args.write_slurm, args.run_tests)
    else:
        mkdirs()
        if args.write_slurm:
            write_slurm_scripts()
        if args.run_tests:
            run_internal_tests()


if __name__ == "__main__":
    main()
