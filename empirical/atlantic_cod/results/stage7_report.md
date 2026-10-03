# Atlantic cod Stage-7 report

Stage 7 is a post-freeze divergence-time sensitivity analysis using the published 250-kb SNAPP MCC window trees. It does not modify or invalidate Stage 5 or Stage 6.

## Primary tests

| LG | n inside | n outside | mean A inside | median A inside | mean A outside | median A outside | Delta A | circular p | physical p | BH p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LG01 | 66 | 43 | 0.755333 | 0.768935 | -0.00169664 | -0.00184131 | 0.75703 | 0.00917431 | 0 | 0.010989 |
| LG02 | 21 | 70 | 0.653045 | 0.680616 | -0.0009494 | -0.000608073 | 0.653994 | 0.010989 | 0 | 0.010989 |
| LG07 | 36 | 81 | 0.438901 | 0.427281 | 0.00339981 | 0.00312807 | 0.435501 | 0.00854701 | 0 | 0.010989 |
| LG12 | 53 | 51 | 0.322829 | 0.326013 | -0.00172302 | -0.00280085 | 0.324552 | 0.00961538 | 0.0727273 | 0.010989 |

## Topology-time comparison

| LG | Delta D | topology p | Delta A | time p | interpretation |
|---|---:|---:|---:|---:|---|
| LG01 | 1.23277 | 0.009174 | 0.75703 | 0.00917431 | topology_shift_and_time_shift |
| LG02 | 0.892332 | 0.010989 | 0.653994 | 0.010989 | topology_shift_and_time_shift |
| LG07 | -0.005981 | 0.555556 | 0.435501 | 0.00854701 | discordant_topology_time_pattern |
| LG12 | 1.04317 | 0.038462 | 0.324552 | 0.00961538 | topology_shift_and_time_shift |

## Posterior uncertainty

Per-window posterior tree files available: 0/426.
Stage 7 therefore reports MCC point-estimate sensitivity only and does not treat posterior samples as genomic replicates.

## Tree validation

Ultrametric MCC trees validated and used: 425. Excluded validation failures: 1. Maximum tip-depth deviation among used trees: 7.62e-10.
