# Atlantic cod Stage-6 report

## A. Local quartet-support profiles

Local support fractions are calculated within each published 250-kb population tree from the predeclared informative quartets. They are not ASTRAL4 branch q1/q2/q3 values. Inside-inversion local support is summarized as: LG01 arrangement-dominant, LG02 arrangement-dominant, LG07 mixed, and LG12 arrangement-dominant.

## B. ASTRAL4 all-windows result

The all-window ASTRAL4 population summary tree uses all 426 published 250-kb local population-tree windows. Its unrooted RF distance to the frozen 12-population baseline is 16 and it recovers 1/9 baseline internal splits.

## C. Collinear-only result

The fully outside treatment uses 245 windows outside all frozen inversion intervals and excludes boundary-overlap windows. Its RF distance to the frozen baseline is 12 and it recovers 3/9 baseline internal splits. This directly tests whether filtering inversion regions changes the inferred population-tree topology.

## D. Inversion-only result

The inversion-only treatment uses 176 fully inside windows as a diagnostic linked-supergene tree set. Its RF distance to the frozen baseline is 14 and it recovers 2/9 baseline internal splits.

## E. Leave-one-inversion-out

Leave-one-inversion-out treatments are reported in `stage6_astral4_treatment_summary.tsv`. They identify which fully inside inversion-window set most changes topology or baseline-split recovery when removed.

## F. Downweighting

Progressive inversion-window subsampling keeps all 245 fully outside windows and retains at most m fully inside windows per LG. The results in `stage6_astral4_downweighting.tsv` describe whether reducing effective inversion representation moves ASTRAL4 population-tree inference toward the collinear baseline topology.

## G. Block-balanced m=1

For m=1, each inversion contributes at most one sampled inside window per replicate, plus all outside windows. Across 100 deterministic replicates, exact baseline matches occurred in 0.000 of replicates, median normalized RF was 0.667 with range 0.667-0.778, and mean baseline-split recovery was 0.288.

## H. Interpretation

Treating many linked inversion windows as separate input trees changes the quartet support available to a summary-tree estimator. The downweighting analysis tests whether reducing their effective representation shifts inference toward the collinear population-history signal. This empirical stress test does not show that ASTRAL4 is wrong, does not prove MSRC, and does not fit MSRC parameters.
