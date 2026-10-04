# Stage 4E ASTER/ASTRAL4 Bridges-2 execution

This directory contains submission scripts only. Do not run ASTER/ASTRAL4 on
the large Neoaves inputs locally.

Default Bridges-2 configuration:

- partition: `RM-shared`
- cores: `32`
- memory: determined by the Bridges-2 `RM-shared` core allocation

Provide your PSC allocation at submission time, or add it to the scripts as:

```bash
#SBATCH -A YOUR_ALLOCATION
```

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
sbatch -A YOUR_ALLOCATION empirical/neoaves_chr4/cluster/stage4e_aster/00_preflight.sbatch
```

After preflight completes, inspect:

```text
$NEOAVES_STAGE4E_OUT/stage4e_preflight_runtime.json
```

Only after preflight passes, submit the first primary deterministic job script:

```bash
sbatch -A YOUR_ALLOCATION empirical/neoaves_chr4/cluster/stage4e_aster/01_full363_primary.sbatch
```

Recommended launch order:

1. Run preflight.
2. Inspect `stage4e_preflight_runtime.json`.
3. Run the primary deterministic job.
4. Inspect T0 runtime, memory behavior, and output.
5. Only then submit sensitivities, J48 jobs, or random-matched jobs.

Do not submit the 500-task random-matched array until the deterministic Stage-4E
treatments have completed and the resulting topology effects have been
inspected.
