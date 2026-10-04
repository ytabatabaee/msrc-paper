#!/usr/bin/env bash
#SBATCH --job-name=neoaves-aster-preflight
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
# TODO/user: add your PSC allocation with:
# #SBATCH -A YOUR_ALLOCATION
# TODO/user: adjust cpus-per-task and time after Bridges-2 preflight if needed

set -euo pipefail

: "${MSRC_REPO:=$PWD}"
: "${NEOAVES_STAGE4E_DATA:=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e}"
: "${NEOAVES_STAGE4E_OUT:=$MSRC_REPO/empirical/neoaves_chr4/results/stage4e/cluster}"
: "${ASTER_BIN:?Set ASTER_BIN to the Bridges-2 astral4 executable}"
: "${PYTHON:=python3}"

mkdir -p "$NEOAVES_STAGE4E_OUT"

# NO LARGE INFERENCE: this checks inputs, hashes, astral4 -h, and a tiny four-taxon smoke test only.
"$PYTHON" "${MSRC_REPO}/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" \
  --preflight \
  --input-manifest "$MSRC_REPO/empirical/neoaves_chr4/results/stage4e_input_manifest.tsv" \
  --aster-bin "$ASTER_BIN" \
  --data-dir "$NEOAVES_STAGE4E_DATA" \
  --output-root "$NEOAVES_STAGE4E_OUT"
