#!/usr/bin/env bash
#SBATCH --job-name=neoaves-aster-random
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
# TODO/user: add your PSC allocation with:
# #SBATCH -A YOUR_ALLOCATION
# TODO/user: adjust cpus-per-task and time after Bridges-2 preflight if needed
#SBATCH --array=1-500

set -euo pipefail

: "${MSRC_REPO:=$PWD}"
: "${NEOAVES_STAGE4E_DATA:=$MSRC_REPO/empirical/neoaves_chr4/data/stage4e}"
: "${NEOAVES_STAGE4E_OUT:=$MSRC_REPO/empirical/neoaves_chr4/results/stage4e/cluster}"
: "${ASTER_BIN:?Set ASTER_BIN to the Bridges-2 astral4 executable}"
: "${PYTHON:=python3}"

mkdir -p "$NEOAVES_STAGE4E_OUT"

MANIFEST="$MSRC_REPO/empirical/neoaves_chr4/results/stage4e_random_matched_manifest.tsv"
ROW=$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 1 {print; exit}' "$MANIFEST")
if [[ -z "$ROW" ]]; then
  echo "No random-matched manifest row for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi
REGIME=$(printf '%s\n' "$ROW" | awk -F '\t' '{print $1}')
REPLICATE=$(printf '%s\n' "$ROW" | awk -F '\t' '{print $2}')
RELATIVE_IDS=$(printf '%s\n' "$ROW" | awk -F '\t' '{print $7}')
TREATMENT="${REGIME}_rep$(printf '%03d' "$REPLICATE")"
INPUT="$NEOAVES_STAGE4E_OUT/materialized/${TREATMENT}.tre"
mkdir -p "$(dirname "$INPUT")"
"$PYTHON" "$MSRC_REPO/empirical/neoaves_chr4/scripts/05_stage4e_prepare_aster.py" \
  --materialize-random "$MSRC_REPO/$RELATIVE_IDS" \
  --source "$NEOAVES_STAGE4E_DATA/T0.tre" \
  --output "$INPUT"
"$PYTHON" "$MSRC_REPO/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" \
  --execute \
  --treatment "$TREATMENT" \
  --aster-bin "$ASTER_BIN" \
  --input "$INPUT" \
  --output-dir "$NEOAVES_STAGE4E_OUT/random_matched/$TREATMENT"
