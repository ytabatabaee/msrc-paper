#!/usr/bin/env bash
#SBATCH --job-name=mouse-tcomplex-astral
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
# TODO/user: add your PSC allocation with:
# #SBATCH -A YOUR_ALLOCATION

set -euo pipefail

: "${MSRC_REPO:=$PWD}"
: "${ASTER_BIN:?Set ASTER_BIN to the Bridges2 astral4 executable}"
: "${HOUSE_MOUSE_STAGE1_OUT:=$MSRC_REPO/empirical/house_mouse_t_complex/results/stage1_aster}"

TREE_FILE="$MSRC_REPO/data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre"
MAP_FILE="$MSRC_REPO/data/house_mouse_t_complex/processed/house_mouse_t_complex_astral_subspecies.map"

mkdir -p "$HOUSE_MOUSE_STAGE1_OUT"

test -s "$TREE_FILE"
test -s "$MAP_FILE"

"$ASTER_BIN" \
  -i "$TREE_FILE" \
  -a "$MAP_FILE" \
  -o "$HOUSE_MOUSE_STAGE1_OUT/T1_ALL_WINDOWS_subspecies.nwk" \
  -u 2 \
  -t "${SLURM_CPUS_PER_TASK:-16}" \
  --length CULength \
  --seed 20261005 \
  > "$HOUSE_MOUSE_STAGE1_OUT/T1_ALL_WINDOWS_subspecies.log" 2>&1

