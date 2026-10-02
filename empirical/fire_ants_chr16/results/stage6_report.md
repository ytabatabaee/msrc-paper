# Stage 6 report

## Purpose

Stage 6 measures sensitivity of seven-group summary-tree inference to counting linked supergene windows as separate local genealogical inputs.

## Stage-5 frozen status

Stage 5 remains frozen for manuscript use and its primary empirical conclusions are unchanged.

## Mapping feasibility

The exact 267-to-7 TWISST group mapping file was unavailable; individual-level mapped ASTRAL4 was not performed.

## Grouped TWISST representation

The analysis used 945 seven-group TWISST topologies and 213 window-level weight distributions.

## Exact weighted summary objective

See `stage6_weighted_summary.tsv` and `stage6_weighted_candidate_scores.tsv`.

## Modal-tree ASTRAL4 design

One modal topology per non-tied window was supplied to ASTRAL4 with `-R -u 2`.

## Treatment results

See `stage6_astral4_treatment_summary.tsv`.

## Focal quartet behavior

See `stage6_focal_split_summary.tsv`.

## Downweighting

See `stage6_astral4_downweighting.tsv` and `stage6_downweighting_summary.tsv`.

## Agreement between weighted and ASTRAL4 analyses

The focal split summary records agreement treatment by treatment.

## Interpretation

Stage 6 is a post-freeze sensitivity demonstration. It does not show that ASTRAL is wrong and does not identify the true historical mechanism.

## Limitations

Modal-tree ASTRAL4 discards TWISST weight-distribution information, and the supergene windows are linked rather than independent evolutionary replicates.

## Stage boundary

No MSRC model was fit, no raw sequences were reanalyzed, and Stage 5 was not analytically modified.

## Stage 6 CU branch-length extension

Existing Stage 6 ASTRAL4 inferred trees were parsed for `CULength` annotations to test whether adding chr16/supergene windows changes coalescent-unit branch lengths when topology is unchanged. The primary comparison is `T_background.inferred.nwk` versus `T_all.inferred.nwk`; `T_chr16_all` and `T_supergene` are secondary comparisons. `SULength` values were not analyzed. Outputs are `stage6_cu_branch_length_comparison.tsv`, `stage6_cu_branch_length_summary.txt`, and the two CU figures.
