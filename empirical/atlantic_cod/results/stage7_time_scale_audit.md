# Stage 7 time-scale audit

Stage 7 inspected the canonical 250-kb window trees in `data/atlantic_cod/processed/cod_window_trees.tsv`, the frozen Zenodo metadata record, and the published SNAPP workflow scripts under `data/atlantic_cod/raw/github_supergenes/`.

## Findings

1. Of the 426 published 250-kb MCC trees, 425 passed the Stage-7 tree validation checks and were ultrametric within numerical tolerance. The maximum observed tip-depth deviation among valid trees was 7.62e-10.
2. One window failed the nonnegative branch-length validation and was excluded from MRCA-time calculations: LG07_023000001_023250000 (Tree contains a materially negative branch length: minimum branch length -0.003195445245; clamp tolerance 1e-10). Negative branch lengths with absolute magnitude no greater than 1e-10 would be clamped to zero as numerical noise, but the excluded LG07 window is materially below zero and was not silently repaired.
3. Branch lengths are time-scaled SNAPP tree lengths from the published BEAST/SNAPP workflow, not raw pairwise sequence distances. Across valid MCC trees, root heights range from 0.0150647 to 1.52873 in the workflow's time units.
4. The workflow constrains the population-tree crown age with `lognormal(0,3.83,0.093)` and describes that constraint as coming from the largest tree-topology subset of the AIM analysis with a prior distribution on root age. The 250-kb MCC trees are therefore calibrated relative to that shared root-age prior.
5. The script evidence available locally does not label the Newick branch lengths explicitly as years, Ma, or substitutions. Because the paper reports divergence times in Ma and the SNAPP XML workflow applies a root-age prior, Stage 7 treats the branch lengths as the published calibrated SNAPP time units. It does not convert them to calendar years or Ma beyond that source-calibrated scale.
6. The same XML-generation script and root-age constraint were used for all 250-kb windows in the public workflow, so the same calibration rule applies consistently across the 426 windows.
7. All tips are sampled population labels from extant Atlantic cod populations, with no tip-date metadata in the window Newicks; tips are treated as contemporaneous.
8. Per-window posterior samples are not available for Stage 7: 0/426 windows have public/local posterior files in the inspected inventory. The upstream workflow references per-window replicate `.trees` and combined `.trees` files, but the frozen Zenodo record lists only supergene-wide/outside-supergene SNAPP posterior files, not per-250-kb-window posterior files.

## Source evidence

- `data/atlantic_cod/raw/github_supergenes/make_snapp_xmls_windows.slurm` writes the root-age constraint as `lognormal(0,3.83,0.093)` and applies it to each generated window XML.
- `data/atlantic_cod/raw/github_supergenes/combine_snapp_results_windows.sh` combines two replicate chains, removes 10% burn-in from each replicate, resamples each to 2000 trees, combines them, and writes MCC `.tre` files with `treeannotator -b 0 -heights mean`.
- `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM6_ESM_source_data_fig4.txt` contains the published 250-kb MCC Newick strings used here.
- `data/atlantic_cod/raw/zenodo/zenodo_4560275_record.json` lists these population-level SNAPP posterior files: gadus_morhua_outside_supergenes_snapp.trees, gadus_morhua_supergene_lg01_snapp.trees, gadus_morhua_supergene_lg02_snapp.trees, gadus_morhua_supergene_lg07_snapp.trees, gadus_morhua_supergene_lg12_snapp.trees.
- Per-window Zenodo keys matching `LG##_start_end` found in the frozen record: 0.

Stage 7 therefore performs biological interpretation only as within-workflow relative divergence-time sensitivity for the same population pairs inside versus outside inversion intervals.
