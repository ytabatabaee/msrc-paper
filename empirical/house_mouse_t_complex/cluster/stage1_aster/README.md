# House mouse t-complex Stage 1 ASTER/ASTRAL4 wrapper

These scripts are placeholders for exploratory Stage 1 ASTRAL4 runs after the Stage 0 mapping gate passes and Stage 1 freezes the ML trees.

Do not submit them before these files exist:

```text
data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre
data/house_mouse_t_complex/processed/house_mouse_t_complex_astral_subspecies.map
```

The mapping file follows the ASTRAL4 convention used elsewhere in this repository: one gene-tree tip followed by one taxon name, separated by a tab, with no header. Tentative and unresolved mappings are excluded by `01_prepare_gene_trees.py`.

This wrapper does not define a genome-wide background. All treatments are within the published chr17 t-complex window range unless future source coordinates establish otherwise.

On Bridges2:

```bash
export MSRC_REPO=$PROJECT/house_mouse_t_complex/repo
export ASTER_BIN=$PROJECT/house_mouse_t_complex/software/astral4
export HOUSE_MOUSE_STAGE1_OUT=$PROJECT/house_mouse_t_complex/outputs/stage1_aster
sbatch -A YOUR_ALLOCATION empirical/house_mouse_t_complex/cluster/stage1_aster/01_exploratory_subspecies.sbatch
```

