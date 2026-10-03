# Stage 2R preanalysis plan — full 2L 50-kb scan

Stage 2R uses the same authenticated MalariaGEN Ag3 release `3.10` data resource and the same frozen 53 homozygous Fontaine-associated samples from Stage 1R.

This plan was amended and re-frozen before querying any new full-scan 2L windows in Stage 2R. The amendment adds a continuous crossed-quartet margin because the Stage-1R discrete quartet statistic `D(w)=q_arrangement-q_species` saturated at or near 1 in several windows, including some flanking controls.

## Fixed inputs

- Ag3 release: `3.10`
- sample set: `fontaine-2015-rebuild`
- bucket: `gs://vo_agam_release_master_us_central1`
- primary samples: the 53 frozen homozygotes in `data/anopheles_2la/processed/stage1r_primary_homozygote_samples.tsv`
- crossed primary classes: gambiae standard, gambiae inverted, coluzzii standard, coluzzii inverted
- frozen inversion interval: `2L:20,524,058-42,165,532`, AgamP4-style 1-based inclusive coordinates
- window grid: non-overlapping 50-kb windows anchored to chromosome coordinate 1: `1-50000`, `50001-100000`, `100001-150000`, ...

Do not add heterokaryotypes to the primary analysis. Do not remove samples based on SNP density, PCA, distance, topology, or missingness unless they violate the frozen sample-level missingness rule. Any exclusion must be documented before downstream topology calculations.

## Window classification

- `inside`: window fully contained in 2La
- `outside`: window fully outside 2La
- `boundary`: window intersects either breakpoint and is excluded from primary inside-vs-outside tests

Boundary classification is determined solely from the frozen coordinates.

## SNP filter

Use the same filter as Stage 1R:

- biallelic SNPs only;
- minor allele count >= 2;
- site missingness <= 0.25;
- sample-level missingness flag threshold <= 0.50;
- no primary LD pruning.

Low-signal windows are defined before topology inspection as windows with fewer than 25 post-filter SNPs or sample-level missingness above the threshold. Low-signal windows remain in the genomic grid but have signal statistics marked missing. They are not replaced by nearby windows.

## Diploid distance metric

Use exactly the Stage-1R unphased diploid dosage distance:

```text
d(i,j) = mean_s |g_i(s)-g_j(s)| / 2
```

over jointly callable retained sites, where `g_i(s)` is alternate-allele dosage 0, 1, or 2.

## Primary crossed distance statistic C(w)

For every usable window calculate, using only crossed gambiae/coluzzii homozygotes:

```text
C(w) = mean d_same_species_opposite_arrangement - mean d_different_species_same_arrangement
```

Positive `C(w)` means same arrangement across species is closer than opposite arrangements within species.

## Primary continuous crossed-quartet margin M(w)

For every crossed quartet containing one sample from each class:

- gambiae standard `G_s`
- gambiae inverted `G_i`
- coluzzii standard `C_s`
- coluzzii inverted `C_i`

calculate:

```text
S_species     = d(G_s,G_i) + d(C_s,C_i)
S_arrangement = d(G_s,C_s) + d(G_i,C_i)
S_third       = d(G_s,C_i) + d(G_i,C_s)
M_q           = S_species - S_arrangement
```

Thus `M_q > 0` means the arrangement partition is closer than the species partition, `M_q < 0` means the species partition is closer, and larger absolute magnitude indicates a stronger preference.

For each window record:

- `mean_M = mean_q M_q`
- `median_M`
- `fraction_M_positive`
- `fraction_M_negative`
- `fraction_M_zero`

The genomic window, not the 2304 quartets, is the primary unit.

## Secondary discrete quartet topology statistic D(w)

Retain the Stage-1R topology-frequency statistic as secondary. For each quartet classify the minimum of `S_species`, `S_arrangement`, and `S_third`, then calculate:

```text
D(w) = q_arrangement - q_species
```

Do not use saturation of `D=1` as evidence of stronger magnitude than another `D=1` window. Use `C(w)` and `M(w)` for magnitude.

## Spatial inference

Primary effect sizes:

```text
Delta_C = mean(C_inside) - mean(C_outside)
Delta_M = mean(M_inside) - mean(M_outside)
```

Median contrasts are also reported. `Delta_D` is descriptive.

Primary spatial null: exact circular shift of each ordered genomic signal track relative to the fixed inversion mask. Include the observed alignment in the exact null enumeration and compute one-sided p-values as `#(Delta_null >= Delta_obs) / N_rotations`.

Secondary robustness test: move an interval with the same physical width as 2La along the observed 2L grid, preserving missing/low-signal windows at their genomic positions. Use finite-sample corrected empirical p-values `(k+1)/(N+1)`.

## Secondary analyses

- all-53 distance class summaries are descriptive because arabiensis is fixed inverted and quadriannulatus is fixed standard;
- NJ trees are descriptive and may be generated for usable 50-kb windows if computationally cheap;
- summary topology analysis should emphasize crossed quartet-frequency summaries rather than forcing ASTRAL or heavy ML tree inference;
- limited haplotype validation uses the 53 homozygotes only, two haplotypes each, in a deterministic 3-inside/3-outside subset and does not run a full 2L haplotype scan;
- no raw-read processing;
- no sample removal based on local topology.

Stage 2R execution begins only after this amended plan is written.
