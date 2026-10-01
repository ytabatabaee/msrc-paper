# Fire ants chromosome 16 Stage 1 source provenance

Stage 1 retrieves and normalizes the published local-tree/TWISST dataset. It does not classify topologies into the Stage-0 focal quartet and does not calculate focal support statistics.

## Frozen upstream source

- Repository: `https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression`
- Branch recorded in Stage 0: `master`
- Frozen commit: `bdc7823941680fa560e61b6467af9f77950f64b0`

## Raw downloaded files

- `data/fire_ants_chr16/raw/upstream/Topology weighting/results/2021-08-07-twisst/input.tree`: `0b0d64217170610b2bd3e8dd1481138b032ef033589887b18a883099c1780e66`
- `data/fire_ants_chr16/raw/upstream/Topology weighting/results/2021-08-07-twisst/topologies_output.trees`: `36a58a44eb730972bc5e4c2003fc9cdd9d1f8681c9ba2c14e053daa2706bd57c`
- `data/fire_ants_chr16/raw/upstream/Topology weighting/results/2021-08-07-twisst/weights_output.csv`: `617edfa9a1d0ae35ed7b0203c73e8a4eb78a6bb40bf32cdaae80ddb5b4182c7c`

## Published construction order

The upstream `Topology weighting/README.md` constructs `results/2021-08-07-twisst/input.tree` by appending windows from `chr1`, then `chr16A`, then `chr16B`, then the chromosome-16 supergene windows. Within each region, trees are appended in increasing upstream window order from `results/window_coordinates`.

Stage 1 maps coordinate rows to `input.tree` lines by this documented construction order. It does not infer or optimize the correspondence from topology.

## Object distinction

- `input.tree` contains 213 full individual-level RAxML local trees.
- `topologies_output.trees` contains the 945 possible seven-group topologies used by TWISST.
- `weights_output.csv` gives the TWISST weight of those grouped topologies for each local tree/window.

The 945 topology inventory is not a set of 945 loci or 945 gene trees.

## Validation counts

- Coordinate rows: `213`
- Local trees: `213`
- TWISST grouped topologies: `945`
- Weight rows: `213`
- Weight columns: `945`
- Seven-group vocabulary: `geminata, saevissima, pusillignis, invicta/macdonaghi_Sb, invicta/macdonaghi_SB, richteri_Sb, richteri_SB`
- Minimum local-tree tips: `267`
- Maximum local-tree tips: `267`
- Distinct local-tree tip sets: `1`

If tip sets differ among windows, the canonical table records per-window `n_tips` and `tip_set_sha256`; no pruning or imputation is performed in Stage 1.
