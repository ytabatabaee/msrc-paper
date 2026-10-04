#!/usr/bin/env bash
#SBATCH --job-name=neoaves-aster-j48
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
TREATMENTS=(J48_T0 J48_T_PNAS J48_T_STRUCT_PRIMARY J48_T_BLOCK_500K J48_NONCHR4)
for treatment in "${TREATMENTS[@]}"; do
  "$PYTHON" "${MSRC_REPO}/empirical/neoaves_chr4/scripts/05_stage4e_run_aster.py" --execute --treatment "$treatment" --aster-bin "$ASTER_BIN" --input "$NEOAVES_STAGE4E_DATA/$treatment.tre" --output-dir "$NEOAVES_STAGE4E_OUT/$treatment"
done
