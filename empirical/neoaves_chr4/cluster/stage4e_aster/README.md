# Stage 4E ASTER/ASTRAL4 Bridges-2 execution

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
