# Anopheles Stage 0R preanalysis plan

## Biological question

Does local genealogy inside the 2La inversion follow inversion arrangement more strongly than species identity, relative to collinear genomic background?

The intended question is: **Does the recombination-suppressed inversion define a spatially coherent genealogy regime associated with arrangement state, above the broader genome-wide discordance caused by introgression?**

This plan does not claim that all arrangement-associated genealogy is caused by the inversion, that introgression is absent, or that MSRC alone explains 2La history.

## Primary sample rule

Use only samples with authoritative homozygous 2La arrangement states:

- `2La/2La`
- `2L+a/2L+a`

Exclude `2La/2L+a` heterokaryotypes from the primary local-tree analysis because an unphased diploid tree mixes two chromosomal arrangements. Retain heterokaryotypes for later phased or secondary analyses.

Do not select samples based on inferred local tree topology.

## Primary data representation

Prefer phased haplotypes if the target Fontaine-associated or Fontaine-compatible samples are present in MalariaGEN haplotype panels across 2L. If phased haplotypes are unavailable for rare species or Fontaine rebuild samples, use diploid genotype distances as a Stage 1R feasibility fallback, but keep this as a weaker design.

## Primary statistic

If explicit local trees are inferred, predefine a window statistic:

`D(w) = Q_arrangement(w) - Q_species(w)`

where `Q_arrangement(w)` is support for grouping same-arrangement lineages across species, and `Q_species(w)` is support for grouping lineages by species identity. Use quartet support/topology weights where the sampled taxa permit unambiguous arrangement and species partitions.

If distance matrices are used first, use:

`D_dist(w) = mean_distance(same_species_opposite_arrangement) - mean_distance(same_arrangement_cross_species)`

Positive `D_dist(w)` means cross-species same-arrangement pairs are closer than same-species opposite-arrangement pairs. Define all eligible pairs before observing distances.

Do not calculate either statistic in Stage 0R.

## Inside versus outside comparison

Inside: 50-kb windows fully contained within the working 2La interval `2L:20524058-42165532`. On a 1-based 2L grid anchored at coordinate 1, this gives 432 fully contained windows; 434 grid windows intersect the inversion.

Outside: all eligible collinear 2L windows outside 2La, excluding boundary-intersecting windows and excluding windows that fail callable-fraction, SNP-density, and missingness thresholds. Prefer same-chromosome controls. If downsampling is required, use a fixed random seed of `20261003` after applying eligibility filters.

Matching criteria should be frozen before topology inference:

- 50-kb length
- callable fraction
- SNP density
- missingness
- recombination environment if an independent map is available

## Feasibility classification rule

- A: full six-species sample-level analysis feasible only if all six focal species have queryable genotypes, clear 2La states, enough homozygous samples in both arrangement classes as needed, and usable sequence variation.
- B: partial multi-species analysis feasible if mainly `gambiae`, `coluzzii`, and `arabiensis` have sufficient data and rare species lack karyotype/phasing support.
- C: Fontaine-only historical comparison required if rare species lack enough sample-level MalariaGEN data or if Fontaine-associated samples cannot be retrieved with reliable 2La states.

Current Stage 0R classification is provisional B/C: MalariaGEN aggregate data identify the Fontaine rebuild and direct regional resources are documented, but sample-level IDs, karyotypes, and haplotype membership were not recovered locally.

## Two-tier strategy

Primary analysis: Fontaine-compatible multi-species design using the `fontaine-2015-rebuild` samples if exact sample IDs, 2La states, and sequence/haplotype access are recovered.

Secondary replication: larger MalariaGEN sample set for species/populations with abundant modern data, likely `gambiae`, `coluzzii`, and `arabiensis`, testing whether arrangement-state signal generalizes with greater sampling.

Neither analysis should run until sample/karyotype and window eligibility are frozen.
