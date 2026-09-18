# Neoaves chromosome 4 Stage 4B: inference treatments

## Inputs / frozen definitions

All Stage 1–4A files are read only. Masks are the unchanged unions of ±250 kb neighborhoods around frozen breakpoints. The 16 leave-one-out units are the frozen primary W250000 event table used in Stage 4A. Reference labels q1/q2/q3 and topology runs are unchanged.

## Treatment definitions

- T0: every resolved window has weight 1.
- T1: consecutive windows with the same dominant topology form a run; each run has total weight 1, shared equally among its windows.
- T2: windows whose midpoint is inside the selected frozen structural mask have weight 0; other windows have weight 1.
- T3: windows whose midpoint is inside the selected frozen structural mask have fixed weight 0.5; other windows have weight 1. This rule is topology neutral and was fixed before Stage 4B execution.

## Quartet-level results

Primary mask results:

| clade | treatment | q1/q2/q3 | dominant | M_species | q1-q2 | 95% block-bootstrap CI for M |
|---|---|---:|---|---:|---:|---:|
| Columbea | T0 | 0.488/0.260/0.253 | q1 | 0.228 | 0.228 | [0.162, 0.281] |
| Columbea | T1 | 0.352/0.331/0.317 | q1 | 0.021 | 0.021 | [-0.009, 0.047] |
| Columbea | T2 | 0.464/0.270/0.265 | q1 | 0.194 | 0.194 | [0.131, 0.244] |
| Columbea | T3 | 0.477/0.265/0.259 | q1 | 0.212 | 0.212 | [0.151, 0.265] |
| N61 | T0 | 0.288/0.454/0.258 | q2 | -0.166 | -0.166 | [-0.231, -0.099] |
| N61 | T1 | 0.346/0.332/0.323 | q1 | 0.014 | 0.014 | [-0.013, 0.037] |
| N61 | T2 | 0.299/0.428/0.273 | q2 | -0.129 | -0.129 | [-0.185, -0.064] |
| N61 | T3 | 0.293/0.442/0.265 | q2 | -0.149 | -0.149 | [-0.208, -0.082] |
| N62 | T0 | 0.454/0.283/0.262 | q1 | 0.171 | 0.171 | [0.057, 0.263] |
| N62 | T1 | 0.318/0.347/0.335 | q2 | -0.030 | -0.030 | [-0.058, -0.005] |
| N62 | T2 | 0.432/0.293/0.275 | q1 | 0.139 | 0.139 | [0.033, 0.241] |
| N62 | T3 | 0.444/0.288/0.268 | q1 | 0.156 | 0.156 | [0.049, 0.252] |

## Margin changes

| clade | DeltaM T1 | DeltaM T2 | DeltaM T3 | Delta q2 T1/T2/T3 |
|---|---:|---:|---:|---:|
| Columbea | -0.207 | -0.034 | -0.016 | 0.072/0.011/0.005 |
| N61 | 0.180 | 0.037 | 0.017 | -0.122/-0.026/-0.012 |
| N62 | -0.200 | -0.031 | -0.015 | 0.064/0.010/0.004 |

## Structural vs generic block decomposition

- **Columbea:** T1 changes M by -0.207; T2 changes it by -0.034; T3 changes it by -0.016. Assigned category: **B. GENERIC BLOCK-LENGTH EFFECT**.
- **N61:** T1 changes M by 0.180; T2 changes it by 0.037; T3 changes it by 0.017. Assigned category: **A. STRUCTURE-LINKED WEIGHTING EFFECT**.
- **N62:** T1 changes M by -0.200; T2 changes it by -0.031; T3 changes it by -0.015. Assigned category: **B. GENERIC BLOCK-LENGTH EFFECT**.

## Mask sensitivity

- Columbea: stringent: T2 DeltaM=-0.023, T3 DeltaM=-0.011; primary: T2 DeltaM=-0.034, T3 DeltaM=-0.016; inclusive: T2 DeltaM=-0.039, T3 DeltaM=-0.017.
- N61: stringent: T2 DeltaM=0.024, T3 DeltaM=0.012; primary: T2 DeltaM=0.037, T3 DeltaM=0.017; inclusive: T2 DeltaM=0.048, T3 DeltaM=0.021.
- N62: stringent: T2 DeltaM=-0.024, T3 DeltaM=-0.011; primary: T2 DeltaM=-0.031, T3 DeltaM=-0.015; inclusive: T2 DeltaM=-0.040, T3 DeltaM=-0.018.

## Leave-one-event-out

- Columbea T2: M_species range 0.190 to 0.203; dominant topologies {'q1': 16}.
- Columbea T3: M_species range 0.210 to 0.216; dominant topologies {'q1': 16}.
- N61 T2: M_species range -0.139 to -0.127; dominant topologies {'q2': 16}.
- N61 T3: M_species range -0.153 to -0.148; dominant topologies {'q2': 16}.
- N62 T2: M_species range 0.135 to 0.149; dominant topologies {'q1': 16}.
- N62 T3: M_species range 0.154 to 0.161; dominant topologies {'q1': 16}.

## Primary questions

1. Block normalization materially redistributes support in all three clades and changes the dominant topology for N61 and N62.
2. Structural exclusion shifts support in the same direction as T1 for each clade, but its effects are much smaller and do not change any dominant topology.
3. For N61, both T1 and T2 reduce q2 dominance; only T1 changes dominance from q2 to q1.
4. Dominant topology changes under T1 for N61 (q2 to q1) and N62 (q1 to q2), and under no other primary treatment.
5. N61 T2 and T3 effects have the same positive direction under stringent, primary, and inclusive masks; neither changes dominance.
6. All 16 N61 leave-one-event-out T2 and T3 analyses remain q2 dominant, with narrow margin ranges, so no single event drives the structural weighting effect.
7. N62's T1 topology change does not persist under T2: T2 remains q1 dominant.
8. N62 is therefore most consistent with a generic block-length effect rather than a frozen structural-region-specific effect.
9. Columbea remains q1 dominant in all treatments despite the large T1 support redistribution.
10. Full species-tree inference was not possible from the repository inputs.

## Full species-tree inference

Not performed. The repository contains chr4 locus coordinates and precomputed focal quartet distance summaries, but no per-locus multi-taxon Newick gene trees and no fixed empirical ASTRAL or other species-tree estimator configuration. A genuine full-tree run requires those gene trees, their locus-coordinate mapping, and the project's fixed estimator executable/version/settings. T3 additionally requires native topology-neutral locus weights; otherwise it remains quartet-summary only.

## Interpretation

- Columbea: **B. GENERIC BLOCK-LENGTH EFFECT**.
- N61: **A. STRUCTURE-LINKED WEIGHTING EFFECT**.
- N62: **B. GENERIC BLOCK-LENGTH EFFECT**.

These comparisons show whether naive dense-window weighting changes the inferred focal relationship relative to block-aware and frozen-structure-aware weighting. They do not identify a true topology, establish species-tree bias, or validate a treatment against independent evidence.

## Limitations

Topology support is the weighted frequency of the frozen dominant quartet state per locus, matching Stage 4A. Focal tracks are correlated. Runs and overlapping structural neighborhoods are evidence clusters, not independent biological replicates. The 95% percentile intervals come from 1,000 paired resamples of complete frozen topology runs; they quantify spatial block uncertainty and are not a chromosome replicate analysis.

## Stage 4C requirements

Stage 4C must compare the predeclared treatment results with independent genome-wide evidence, while preserving the frozen labels and treatment rules. Only that external comparison can assess whether any treatment aligns with independent evidence and support language about species-tree bias.
