# Neoaves chromosome 4 Stage 4C: independent validation

## Independent reference search

The repository search identified the following candidates before treatment agreement was evaluated:

| candidate | available | evidence | independence |
|---|---|---|---|
| nonchr4_genomewide_QQS_v1 | True | F (genome-wide quartet summaries): `data/neoaves_chr4/raw/genetreesupport/63K_trees.names_header.txt.xz;data/neoaves_chr4/raw/genetreesupport/clade-rec.stat.xz` | yes after prospective exclusion of every chr4-labelled sequence |
| published_q1_branch_definitions | False | E (not usable as a reference tree): `data/neoaves_chr4/raw/genetreesupport/draw-movingaverage.r;data/neoaves_chr4/raw/genetreesupport/clade-analysis/README.md;empirical/neoaves_chr4/results/stage3b_v2_quartet_role_mapping.tsv` | cannot establish independence or extract a reference topology |
| resolved_genetrees_absent | False | A/C (unavailable): `resolved-genetrees.tre.gz mentioned in clade-analysis/README.md but absent` | potentially, but unavailable |

No per-locus Newick gene trees, non-chr4 species tree, frozen ASTRAL tree, or fixed empirical species-tree estimator configuration is present. The README mentions `resolved-genetrees.tre.gz`, but that file is absent.

## Reference selection rule

The rule was fixed before comparison: choose the first available item in A→F order, preferring an analysis that can exclude chr4. No usable A–E item exists: the repository has no gene trees, species tree, or extractable published reference tree. The q1 quadripartition definitions label tested branches but do not provide a standalone independent topology. The primary reference is therefore the available F item, `nonchr4_genomewide_QQS_v1`, after excluding every metadata label whose chromosome group is 4. Selection did not use agreement with T0–T3.

## Reference topology/support

Support is the frequency of the dominant q1/q2/q3 state per non-chr4 locus, matching the Stage 4B estimand. Intervals are 95% chromosome-cluster bootstrap intervals.

| clade | non-chr4 q1/q2/q3 | topology | margin | 95% margin CI | loci | chromosome units |
|---|---:|---|---:|---:|---:|---:|
| Columbea | 0.356/0.310/0.335 | q1 | 0.021 | [0.014, 0.030] | 55984 | 32 |
| N61 | 0.364/0.329/0.307 | q1 | 0.036 | [0.029, 0.043] | 56247 | 32 |
| N62 | 0.317/0.352/0.331 | q2 | 0.022 | [0.015, 0.031] | 54373 | 32 |

## T0-T3 comparison

| clade | treatment | dominant | matches reference | M_reference | Delta vs T0 | 95% M_reference CI |
|---|---|---|---|---:|---:|---:|
| Columbea | T0 | q1 | True | 0.228 | 0.000 | [0.162, 0.281] |
| Columbea | T1 | q1 | True | 0.021 | -0.207 | [-0.009, 0.047] |
| Columbea | T2 | q1 | True | 0.194 | -0.034 | [0.131, 0.244] |
| Columbea | T3 | q1 | True | 0.212 | -0.016 | [0.151, 0.265] |
| N61 | T0 | q2 | False | -0.166 | 0.000 | [-0.231, -0.099] |
| N61 | T1 | q1 | True | 0.014 | 0.180 | [-0.013, 0.037] |
| N61 | T2 | q2 | False | -0.129 | 0.037 | [-0.185, -0.064] |
| N61 | T3 | q2 | False | -0.149 | 0.017 | [-0.208, -0.082] |
| N62 | T0 | q1 | False | -0.171 | 0.000 | [-0.265, -0.057] |
| N62 | T1 | q2 | True | 0.012 | 0.183 | [-0.018, 0.037] |
| N62 | T2 | q1 | False | -0.139 | 0.031 | [-0.241, -0.033] |
| N62 | T3 | q1 | False | -0.156 | 0.015 | [-0.252, -0.049] |

## N61 focal interpretation

The independent non-chr4 reference is q1. T0 is q2 dominant and disagrees. T1 increases M_reference by 0.180 and becomes slightly q1 dominant. T2 and T3 increase M_reference by 0.037 and 0.017, respectively, while remaining q2 dominant. All three changes are in the direction predicted by the earlier structural q2 enrichment result; the structural-only effects are modest.

## Columbea

The non-chr4 reference and T0–T3 are all q1 dominant. T1 substantially reduces the q1 margin, while T2/T3 have smaller negative effects. T0 already agrees with the independent reference.

## N62

The non-chr4 reference is q2, whereas T0 is q1 dominant. T1 moves strongly toward the reference and changes dominance to q2. T2 and T3 move modestly toward q2 but remain q1 dominant. Thus N62 also shows empirical weighting bias relative to non-chr4 evidence, but its improvement is specifically block-aware: Stage 4A found structural-region q2 depletion, and structural exclusion does not reproduce the topology change.

## Genome-wide/full-tree analysis if available

A full species-tree inference was not run. Per-locus multi-taxon Newick trees and a fixed empirical estimator/configuration are absent. This is topology-reference validation using independent non-chr4 focal quartet summaries.

## Control chromosomes if available

- Columbea: 5/29 non-chr4 main chromosomes change dominant topology after run normalization. Chr4 |Delta M_reference|=0.207, empirical percentile=100.0% among non-chr4 controls.
- N61: 7/29 non-chr4 main chromosomes change dominant topology after run normalization. Chr4 |Delta M_reference|=0.180, empirical percentile=100.0% among non-chr4 controls.
- N62: 6/29 non-chr4 main chromosomes change dominant topology after run normalization. Chr4 |Delta M_reference|=0.183, empirical percentile=96.6% among non-chr4 controls.

## Uncertainty

Reference intervals use 10,000 chromosome-cluster bootstrap replicates. Dense loci are never resampled independently. Treatment intervals and paired deltas are the frozen Stage 4B complete-run bootstrap results. Chromosome jackknifing removes one non-chr4 chromosome unit at a time.

- Columbea: jackknife topology counts {'q1': 32}; M_reference range 0.020 to 0.024.
- N61: jackknife topology counts {'q1': 32}; M_reference range 0.034 to 0.038.
- N62: jackknife topology counts {'q2': 32}; M_reference range 0.019 to 0.024.

## Interpretation category

- Columbea: **C — NAIVE CHR4 AGREES WITH REFERENCE**
- N61: **A — EMPIRICAL BIAS / PARTIAL CORRECTION**
- N62: **A — EMPIRICAL BIAS / PARTIAL CORRECTION**

## What claim is now justified

For N61, naive dense chr4 locus weighting produces a focal relationship that disagrees with stable non-chr4 genome-wide quartet evidence. Block normalization reverses the dominant relationship toward that reference; independently frozen structural exclusion/downweighting moves support in the same direction but does not reverse dominance. This supports focal empirical species-tree bias from dense chr4 weighting and strong block-aware, partial structural-aware correction relative to the independent reference.

N62 independently shows naive chr4 disagreement and block-aware correction toward its q2 non-chr4 reference, while its structural-aware effects are modest and do not reverse dominance. Columbea already agrees with its q1 reference, and all weighting treatments reduce rather than improve its reference margin.

## Remaining limitations

The independent evidence consists of published per-locus focal quartet summaries, not recoverable multi-taxon gene trees or a re-estimated full species tree. Locus states on the same chromosome remain linked, which is why uncertainty uses chromosomes. The reference and chr4 tracks derive from the same original genome-wide gene-tree study, although the primary baseline excludes chr4 completely. The result validates focal relationships and weighting behavior; it does not establish a causal rearrangement mechanism or resolve the underpowered Stage 3B 2:2 test.
