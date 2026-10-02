# Fire ants chromosome 16 social supergene

This directory scaffolds a staged empirical analysis for the chromosome-16 social supergene in fire ants of the genus `Solenopsis`.

Stage 0 complete.

Stage 1 complete.

Stage 2 complete.

Stage 3 complete.

Stage 4A complete.

Stage 4B complete — spatial-null inference. It implements the frozen exact circular-shift null and the frozen physical-coordinate-aware same-width interval sensitivity. It does not optimize supergene boundaries, run ASTRAL/ASTER, fit MSRC parameters, or identify a historical mechanism.

Stage 5 complete — frozen for manuscript use.

## Biological system

The focal system is the fire-ant social chromosome described in:

```text
Stolle, E., Pracana, R., Lopez-Osorio, F. et al.
Recurring adaptive introgression of a supergene variant that determines social organization.
Nat Commun 13, 1180 (2022).
DOI: 10.1038/s41467-022-28806-7
```

The public analysis repository is:

```text
https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression
```

The social supergene lies on chromosome 16. The two structural/social haplotypes are `SB` and `Sb`.

## Conceptual role in MSRC paper

The intended use in the MSRC paper is a clean empirical demonstration that a recombination-suppressed structural haplotype defines a coherent alternative genealogy regime across species.

This dataset may also illustrate the difficulty of distinguishing structural-history and introgression/hybridization explanations from quartet patterns alone. It must not be framed as evidence that introgression was absent.

## Public sources

Stage 0 records these sources in `../../data/fire_ants_chr16/metadata/source_manifest.tsv`:

- Stolle et al. 2022, `Nat Commun 13, 1180`, DOI `10.1038/s41467-022-28806-7`.
- Upstream analysis repository `wurmlab/2021-fire-ant-social-supergene-introgression`.
- Wang et al. 2013, `Nature 493, 664-668`, DOI `10.1038/nature11832`.
- Yan et al. 2020, `Nat Ecol Evol 4, 240-249`, DOI `10.1038/s41559-019-1081-1`.

Stage 0 stores only selected small non-topological upstream files under `../../data/fire_ants_chr16/raw/upstream/`. It does not copy the full upstream repository.

## Repository layout

- `scripts/`: executable stage scripts.
- `results/`: stage reports and frozen design summaries.
- `figures/`: rendered Stage-4A descriptive support figures.
- `config/`: reserved for later stage configuration.
- `tests/`: test fixtures or notes if later stages need them.
- `../../data/fire_ants_chr16/`: data tree with `raw/`, `metadata/`, `intermediate/`, and `processed/`.

## Structural haplotypes

Stage 0 freezes the biological haplotype vocabulary independently of local phylogenetic topology:

- `SB`: one social-chromosome haplotype.
- `Sb`: alternative social-chromosome haplotype carrying major rearrangements/inversions.

The published work reports that `Sb` differs structurally from `SB` by multiple large inversions and has strongly suppressed recombination with `SB`. Haplotype state must not be inferred from local topology.

## Author-defined genomic regions

The author-defined region vocabulary comes from the allowed upstream coordinate file `Topology weighting/results/window_coordinates`:

- `chr1`: non-supergene control chromosome.
- `chr16A`: collinear/recombining chromosome-16 region outside the supergene.
- `chr16B`: collinear/recombining chromosome-16 region outside the supergene.
- `chr16_supergene`: author-designated chromosome-16 supergene region.

The first and last BUSCO-window coordinates are treated as observed analysis-window spans, not nucleotide-resolved inversion breakpoints.

## Focal quartet

Stage 0 freezes one canonical focal quartet:

```text
A = invicta/macdonaghi_SB
B = invicta/macdonaghi_Sb
C = richteri_SB
D = richteri_Sb
```

The three possible unrooted resolutions are:

- `species_split`: `AB|CD`
- `haplotype_split`: `AC|BD`
- `third_split`: `AD|BC`

The canonical table is `../../data/fire_ants_chr16/processed/stage0_focal_quartet.tsv`, with checksum stored in `../../data/fire_ants_chr16/processed/stage0_focal_quartet.sha256`.

## Primary future statistic

Later stages will derive per-window support quantities:

```text
q_S(w): support for species_split
q_H(w): support for haplotype_split
q_3(w): support for third_split
```

For windows with complete resolved topology weight, the predeclared identity is:

```text
q_S + q_H + q_3 = 1
```

The primary scalar contrast will be:

```text
D(w) = q_H(w) - q_S(w)
```

Stage 0 does not calculate these quantities.

## Anti-circularity design

The primary future comparison is within chromosome 16:

```text
chr16_supergene
vs.
chr16A + chr16B
```

The predeclared contrast is:

```text
Delta_D = mean(D_supergene) - mean(D_chr16 outside)
```

The chromosome-1 comparison is a secondary descriptive control:

```text
chr16_supergene
vs.
chr1
```

The later primary null will use exact circular shifts of the ordered chromosome-16 `D` track with the author-defined supergene mask fixed. A later physical-coordinate-aware block-placement sensitivity analysis is also predeclared because four-BUSCO-gene windows are irregularly spaced.

## Introgression caveat

The 2022 source paper concludes that the `Sb` supergene originated in the `S. invicta/macdonaghi` lineage and repeatedly introgressed into other species. This analysis must not claim that MSRC explains fire-ant history without gene flow, that the inversion caused introgression, or that topology alone proves MSRC.

## Exploratory-analysis disclosure

Dataset selection and feasibility were exploratory and informed by the published study and an exploratory feasibility calculation of the public TWISST resources.

The formal focal quartet definition is biologically predefined from species identity and `SB`/`Sb` state. The formal pipeline is frozen before reproduction of local support statistics inside this repository.

## Planned stages

Stage 0: Scaffold, provenance, structural definitions, focal quartet, and analysis-design freeze.

Stage 1: Retrieve and normalize the public 213-window coordinate/local-tree/TWISST datasets without testing the focal prediction.

Stage 2: Freeze sample/group-to-species and `SB`/`Sb` mappings independently of local topology support.

Stage 3: Freeze the background species-history topology using the published chromosome 1-15 ASTRAL result, independently of the chromosome-16 local support track.

Stage 4A: First formal local-topology stage: derive `q_species`, `q_haplotype`, `q_third`, and `D` from TWISST weights.

Stage 4B: Complete — spatial-null testing of supergene enrichment on chromosome 16.

Stage 5: Robustness/manuscript freeze.

Stage 6: Complete — post-freeze grouped TWISST / ASTRAL4 summary-tree sensitivity analysis.

## Reproduction

From the repository root:

```bash
python3 empirical/fire_ants_chr16/scripts/00_freeze_analysis_design.py --run-tests
python3 empirical/fire_ants_chr16/scripts/01_retrieve_normalize_twisst.py --skip-download --run-tests
python3 empirical/fire_ants_chr16/scripts/02_freeze_group_states.py --skip-download --run-tests
python3 empirical/fire_ants_chr16/scripts/03_freeze_background_history.py --skip-download --run-tests
python3 empirical/fire_ants_chr16/scripts/04a_local_quartet_support.py --run-tests
python3 empirical/fire_ants_chr16/scripts/04b_spatial_null_tests.py --run-tests
python3 empirical/fire_ants_chr16/scripts/05_finalize_fire_ant_analysis.py --run-tests
python3 empirical/fire_ants_chr16/scripts/06_summary_tree_sensitivity.py --run-tests
python3 empirical/fire_ants_chr16/scripts/06b_compare_cu_branch_lengths.py --run-tests
```

The Stage-0 script validates the manifests, verifies the upstream coordinate table, counts author-defined regions, validates focal groups and quartet resolutions, writes the focal-quartet checksum, and enforces forbidden-input guards.

The Stage-1 script verifies the frozen Stage-0 checksums, validates the immutable raw local-tree/TWISST files, and regenerates:

- `../../data/fire_ants_chr16/processed/stage1_window_trees.tsv`
- `../../data/fire_ants_chr16/processed/stage1_twisst_topologies.tsv`
- `../../data/fire_ants_chr16/processed/stage1_twisst_weights_sparse.tsv`
- `../../data/fire_ants_chr16/processed/stage1_twisst_weight_summary.tsv`
- `../../data/fire_ants_chr16/processed/stage1_window_index.tsv`
- `results/stage1_topology_crosscheck.tsv`
- `results/stage1_validation_summary.tsv`
- `results/stage1_source_provenance.md`
- `results/stage1_manifest.json`

The Stage-2 script verifies the frozen Stage-0 and Stage-1 checksums, validates the upstream group-construction script, and regenerates:

- `../../data/fire_ants_chr16/metadata/twisst_group_definition_rules.tsv`
- `../../data/fire_ants_chr16/processed/stage2_twisst_groups.tsv`
- `../../data/fire_ants_chr16/processed/stage2_sample_label_inventory.tsv`
- `../../data/fire_ants_chr16/processed/stage2_focal_partition.tsv`
- `../../data/fire_ants_chr16/processed/stage2_focal_partition.sha256`
- `results/stage2_report.md`
- `results/stage2_manifest.json`

The Stage-3 script verifies the frozen Stage-0, Stage-1, and Stage-2 checksums, validates the published chromosome 1-15 ASTRAL background tree, and regenerates:

- `../../data/fire_ants_chr16/metadata/background_tree_provenance.tsv`
- `../../data/fire_ants_chr16/processed/stage3_background_tree.nwk`
- `../../data/fire_ants_chr16/processed/stage3_background_tree.sha256`
- `../../data/fire_ants_chr16/processed/stage3_background_partition.tsv`
- `../../data/fire_ants_chr16/processed/stage3_background_partition.sha256`
- `results/stage3_report.md`
- `results/stage3_manifest.json`

Stage 4A is the first formal unblinded local-topology stage.

The Stage-4A script verifies the frozen Stage-0 through Stage-3 checksum chain, classifies the 945 grouped TWISST topologies by the independently frozen focal quartet, aggregates published weights into local quartet support quantities, and regenerates:

- `../../data/fire_ants_chr16/processed/stage4a_topology_classes.tsv`
- `../../data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv`
- `results/stage4a_region_summary.tsv`
- `results/stage4a_primary_contrast.tsv`
- `results/stage4a_report.md`
- `results/stage4a_figure_caption.md`
- `results/stage4a_manifest.json`
- `figures/fire_ants_chr16_quartet_support.pdf`
- `figures/fire_ants_chr16_quartet_support.png`
- `figures/fire_ants_chr1_quartet_support_control.pdf`
- `figures/fire_ants_chr1_quartet_support_control.png`
- `figures/fire_ants_chr16_D_track.pdf`
- `figures/fire_ants_chr16_D_track.png`

The Stage-4B script verifies the frozen Stage-0 through Stage-4A checksum chain, reconstructs physical chromosome-16 order by midpoint, runs the exact circular-shift null and physical-coordinate same-width interval sensitivity, and regenerates the Stage-4B tables, report, manifest, and spatial-null figures.

## Stage 5 manuscript freeze

- `results/FINAL_ANALYSIS.md`
- `results/fire_ants_main_table.tsv`
- `results/fire_ants_methods_text.md`
- `results/fire_ants_results_text.md`
- `figures/fire_ants_main.pdf`
- `figures/fire_ants_supplement.pdf`

Stage 6 is complete as a post-freeze grouped TWISST / ASTRAL4 summary-tree sensitivity analysis and does not alter the Stage-5 primary empirical result.

Stage 6 does not alter the Stage-5 primary empirical analysis.

Stage 6 CU branch-length extension outputs: `results/stage6_cu_branch_length_comparison.tsv`, `results/stage6_cu_branch_length_summary.txt`, `figures/fire_ants_stage6_cu_background_vs_combined.pdf`, and `figures/fire_ants_stage6_cu_delta.pdf`.
