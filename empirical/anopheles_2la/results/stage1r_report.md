# Stage 1R — MalariaGEN regional SNP-access and homozygote pilot

Stage 1R tested whether the authenticated MalariaGEN sample-level data are technically usable for a 2La arrangement-vs-species analysis before any full 2L scan. It used the frozen Stage 1A sample and karyotype data, the fixed AgamP4 interval `2L:20524058-42165532`, and eight coordinate-predefined 50-kb pilot windows. No raw reads were downloaded, no full-chromosome scan was run, and no final manuscript-wide statistic was calculated.

## Environment

- Python: `3.10.0 (v3.10.0:b494f5935c, Oct  4 2021, 14:59:20) [Clang 12.0.5 (clang-1205.0.22.11)]`
- executable: `/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria310/bin/python`
- Python 3.11 migration: not attempted beyond PATH probe: no python3.11 executable available; Stage 1R used validated Python 3.10 environment
- malariagen-data: `14.0.0`
- Ag3 release: `3.10`
- sample set: `fontaine-2015-rebuild`
- bucket: `gs://vo_agam_release_master_us_central1`

## Primary cohort

The primary Stage 1R cohort contains exactly 53 homozygotes. The 19 heterokaryotypes remain in metadata and are excluded from the primary unphased diploid analysis.

## SNP QC

- outside_left_distal_10Mb: 5349 post-filter SNPs, mean missingness 0.0331, status ok
- outside_left_flank_50kb: 5882 post-filter SNPs, mean missingness 0.0681, status ok
- inside_left_boundary_50kb: 4735 post-filter SNPs, mean missingness 0.1123, status ok
- inside_deep_q1: 5784 post-filter SNPs, mean missingness 0.0727, status ok
- inside_deep_center: 6956 post-filter SNPs, mean missingness 0.0532, status ok
- inside_deep_q3: 5711 post-filter SNPs, mean missingness 0.0209, status ok
- inside_right_boundary_50kb: 6161 post-filter SNPs, mean missingness 0.0649, status ok
- outside_right_flank_50kb: 6673 post-filter SNPs, mean missingness 0.0564, status ok

## Crossed gambiae-coluzzii pilot C(w)

`C(w) = mean same-species/opposite-arrangement distance - mean different-species/same-arrangement distance`. Positive values indicate that arrangement similarity can overcome species identity in this pilot contrast.

- outside_left_distal_10Mb: C=-0.00729983
- outside_left_flank_50kb: C=0.00388264
- inside_left_boundary_50kb: C=0.135879
- inside_deep_q1: C=0.113235
- inside_deep_center: C=0.0644514
- inside_deep_q3: C=0.0616843
- inside_right_boundary_50kb: C=0.127423
- outside_right_flank_50kb: C=0.01526

## Crossed quartet pilot D(w)

`D(w)=q_arrangement-q_species` from all gambiae-standard, gambiae-inverted, coluzzii-standard, coluzzii-inverted individual quartets.

- outside_left_distal_10Mb: D=-0.654948, q_arr=0.144, q_species=0.799
- outside_left_flank_50kb: D=0.68099, q_arr=0.839, q_species=0.158
- inside_left_boundary_50kb: D=1, q_arr=1.000, q_species=0.000
- inside_deep_q1: D=1, q_arr=1.000, q_species=0.000
- inside_deep_center: D=1, q_arr=1.000, q_species=0.000
- inside_deep_q3: D=1, q_arr=1.000, q_species=0.000
- inside_right_boundary_50kb: D=1, q_arr=1.000, q_species=0.000
- outside_right_flank_50kb: D=0.979601, q_arr=0.989, q_species=0.010

Mean C inside/outside: `0.100534` / `0.0039476`.
Mean D inside/outside: `1` / `0.335214`.

## Haplotype membership

Regional haplotype availability was checked but not used as the primary Stage 1R analysis. `72` of 72 Fontaine-associated samples were returned by the regional haplotypes query.

## Interpretation

The primary biological inference for Stage 1R is restricted to the crossed gambiae + coluzzii design, because those two species contain both homozygous arrangements. Arabiensis and quadriannulatus are useful anchors in all-53 visualizations, but arrangement and species are fully confounded within each, so they are not allowed to drive the claim that arrangement overrides species.

The pilot does not claim formal significance. It asks whether direct regional SNP access works, whether the fixed homozygote cohort produces usable distances and NJ visualizations, and whether inside windows show a qualitatively stronger arrangement-associated signal than outside controls.

Decision: **PASS**.

Stage 2R has been frozen in `stage2r_preanalysis_plan.md`.
