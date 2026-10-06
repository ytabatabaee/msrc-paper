# House-mouse t-complex Stage 2B report

Stage 2B corrects split-specific ASTRAL annotations, decomposes quartet
support by the standard/pseudo-t arrangement of the three M. musculus-complex
tips, finishes exact balanced resampling, rescues the third filtering dataset,
and completes the spatial phase sweep. These remain linked 5-kb windows from
the chr17 t-complex region, not 4,046 independent loci or genome-wide
background.

The independently specified focal baseline is
`Q_SPECIES = (musculus, castaneus) | (domesticus, spretus)`. The competing
split is `Q_T_ALT = (domesticus, musculus) | (castaneus, spretus)`.

## Corrected ASTRAL annotations

The previous compact summary copied q values from the Germany–France branch of
the seven-population tree into the T0 focal row. The corrected
split-specific values are:

| treatment | topology | CU | localPP | q1 | q2 | q3 |
|---|---|---:|---:|---:|---:|---:|
| T0_STANDARD | Q_SPECIES | 0.056283 | 0.999985 | 0.369910 | 0.343656 | 0.286434 |
| T1_ALL_WINDOWS | Q_T_ALT | 0.049097 | 1.000000 | 0.365364 | 0.323323 | 0.311313 |

Each q triplet is tied to the focal canonical split and sums to one within
floating-point tolerance. Seven-population values are now reported in the
split-specific population tables rather than reduced to the first annotated
branch.

## Arrangement-state decomposition

The pattern order is domesticus / musculus / castaneus, with S for
standard/noncarrier and T for reconstructed/pseudo-t haplotype.

| pattern | q_species | q_t_alt | q_other | unresolved | delta species-alt |
|---|---:|---:|---:|---:|---:|
| SSS | 0.369803 | 0.343630 | 0.286365 | 0.000202 | 0.026174 |
| SST | 0.231582 | 0.306003 | 0.461962 | 0.000453 | -0.074420 |
| STS | 0.199726 | 0.500016 | 0.299799 | 0.000459 | -0.300289 |
| STT | 0.685425 | 0.157016 | 0.157304 | 0.000256 | 0.528409 |
| TSS | 0.398753 | 0.329165 | 0.271871 | 0.000212 | 0.069588 |
| TST | 0.153106 | 0.153167 | 0.693594 | 0.000133 | -0.000062 |
| TTS | 0.130993 | 0.706643 | 0.162107 | 0.000257 | -0.575649 |
| TTT | 0.455309 | 0.267625 | 0.276911 | 0.000154 | 0.187684 |

Homogeneous SSS and TTT both favor Q_SPECIES. The pooled one-T class favors
Q_T_ALT (q_species 0.257808, q_t_alt 0.398295, delta -0.140487), while the
pooled two-T class is slightly species-history leaning (q_species 0.375194,
q_t_alt 0.340423). The strongest alternative-supporting configurations are
TTS, where domesticus and musculus carry T, and STS, where musculus alone
carries T. The strongest species-history configuration is STT; TSS also
favors Q_SPECIES. TST is nearly balanced but has most support on Q_OTHER.

The ALL_TIPS shift is therefore an arrangement-state-associated quartet
distortion concentrated in particular mixed configurations, especially the
one-T class and TTS/STS patterns. It is not evidence that all
pseudo-t-haplotype trees support Q_T_ALT: TTT favors Q_SPECIES, as does the
aggregate T_ONLY analysis.

## Balanced resampling

Across 1,000 deterministic balanced replicates using all 4,046 windows:

| treatment | Q_SPECIES | Q_T_ALT | Q_OTHER |
|---|---:|---:|---:|
| B0_STANDARD_MATCHED | 0.537 | 0.458 | 0.005 |
| B1_T_MATCHED | 1.000 | 0.000 | 0.000 |
| B2_MIXED_BALANCED | 0.268 | 0.658 | 0.074 |

The mixed-state shift survives balanced tip sampling. The T-only matched
treatment remains Q_SPECIES, reinforcing that joint representation of
standard and pseudo-t arrangements is the relevant sampling contrast.

The exact quartet winner matched the ASTRAL topology in all 60 validation
replicates: 20/20 for B0, 20/20 for B1, and 20/20 for B2.

## Filtering robustness

All three published ML filtering datasets are usable after deterministic label
normalization. The normalization removes archive-specific
`OverallCovFiltered` and `_DelRemoved` artifacts and restores the `.fa`
suffix on `_tHaplSubset` labels; the audit found 55/55 exact one-to-one
matches and no fuzzy matches.

| dataset | usable trees | standard q_species / q_t_alt | all-tip q_species / q_t_alt |
|---|---:|---:|---:|
| PASS_SNPs | 4,047 | 0.439601 / 0.353502 | 0.339790 / 0.382868 |
| Coverage + AlleleRatio RAW | 4,046 | 0.369803 / 0.343630 | 0.323187 / 0.365290 |
| Coverage RAW | 4,046 | 0.368636 / 0.344656 | 0.313129 / 0.369594 |

Each dataset returns STANDARD_ONLY = Q_SPECIES and ALL_TIPS = Q_T_ALT.

## Seven-population split analysis

The P0/P1 comparison has two changed canonical splits. P0 loses:

`Afghanistan|CAST|Czech_Republic|Kazakhstan || France|Germany|SPRE`

P1 gains:

`Afghanistan|Czech_Republic|France|Germany|Kazakhstan || CAST|SPRE`

The shared split
`Afghanistan|CAST|Czech_Republic|Kazakhstan|SPRE || France|Germany`
changes CU from 0.886491 to 0.694980. The shared split
`Afghanistan|CAST|France|Germany|SPRE || Czech_Republic|Kazakhstan`
changes CU from 0.529620 to 0.157049. These are sensitivity results for
population history and gene flow; neither tree is labeled correct or wrong.

## Complete spatial phase sweep

All feasible 5-kb phase offsets were evaluated using exact quartet winners;
ASTRAL was retained for the first ten phases at each spacing.

| spacing | phases | Q_SPECIES | Q_T_ALT | Q_OTHER | median localPP | median CU |
|---:|---:|---:|---:|---:|---:|---:|
| 5 kb | 1 | 0.00 | 1.00 | 0.00 | 1.000 | 0.0491 |
| 10 kb | 2 | 0.00 | 1.00 | 0.00 | 0.995 | 0.0476 |
| 20 kb | 4 | 0.00 | 1.00 | 0.00 | 0.952 | 0.0473 |
| 50 kb | 10 | 0.00 | 1.00 | 0.00 | 0.732 | 0.0388 |
| 100 kb | 20 | 0.10 | 0.90 | 0.00 | 0.650 | 0.0425 |
| 250 kb | 50 | 0.10 | 0.86 | 0.04 | 0.533 | 0.0421 |
| 500 kb | 100 | 0.11 | 0.74 | 0.15 | 0.497 | 0.0421 |

Dense linked windows preserve a modal Q_T_ALT topology while strongly
inflating apparent localPP. Spatial thinning reduces support and increases
phase-dependent topology variation; it does not produce a simple monotonic
change in CU.

## Source-defined non-recombined regions

The archive README documents the concatenated non-recombined-region tree but
does not provide interval coordinates. The archive contains concatenated tree
files and their directories for this purpose; no coordinate list, code, or
machine-readable interval table was recovered. Thus
`source_defined_blocks_available = false`, and no blocks were inferred from
the observed quartet signal.

## Scientific conclusion and limitations

The t-complex exhibits arrangement-state-dependent quartet distortion. Neither
standard-only nor t-only sampling alone produces the aggregate quartet
distribution observed when standard and pseudo-t arrangements are sampled
jointly. Dense linked windows strongly inflate confidence in the resulting
topology. These observations are consistent with an MSRC-like mechanism in
which a persistent structural polymorphism alters genealogical sampling across
lineages, but they do not establish causality or prove MSRC. The remaining
limitation is that these are published local trees from one chr17 region, with
no independent genome-wide control and no direct re-estimation of genealogies.

## Final visualization cleanup

The four-group ASTRAL plot now contains exactly two explicit cladogram panels
at
`figures/house_mouse_t_complex_astral_4group.png` and `.pdf`. It shows the
standard-only `Q_SPECIES` tree beside the all-tip `Q_T_ALT` tree. The
ALL-TIPS annotation now correctly reports the biological values
`q_species = 0.323323`, `q_t_alt = 0.365364`, and `q_other = 0.311313`;
these names are not fixed aliases for ASTRAL q1/q2/q3. ASTRAL inference is
unrooted and *M. spretus* is used only to orient the visualization. The
seven-population plot also contains exactly two explicit tree panels at
`figures/house_mouse_t_complex_astral_7population.png` and `.pdf`; it marks the
one P0 split lost and one P1 split gained, with RF(P0,P1) = 2, and does not
collapse the comparison to one branch length.

Direct q tracks are provided in
`figures/house_mouse_t_complex_q_tracks.png` and `.pdf`, with raw 5-kb values
and fixed 250-kb means for STANDARD_ONLY and ALL_TIPS. The aligned
q/recombination figure is at
`figures/house_mouse_t_complex_q_recombination.png` and `.pdf`; the
recombination-only track has an explicit `Fraction of 5-kb windows` y-axis and
is at
`figures/house_mouse_t_complex_recombination_track.png` and `.pdf`. The
The manuscript-oriented `main_v3` composite is at
`figures/house_mouse_t_complex_main_v3.png` and `.pdf`; the regenerated
compatibility `main_v2` output is at
`figures/house_mouse_t_complex_main_v2.png` and `.pdf`; balanced and thinning
controls are separated into `figures/house_mouse_t_complex_controls.png` and
`.pdf`.

The recombination track directly uses the authors' archived
`Tree_results_dom_ML`, `Tree_results_mus_ML`, and `Tree_results_cas_ML` files.
All 4,046 frozen windows for each subspecies match a source classification row,
and all observed source codes have an explicit mapping to the published
classes. The audit is recorded in
`stage2_recombination_source_code_audit.tsv`; this is a source-row coverage
and code-definition audit rather than an independent reclassification test.
The result is a phylogeny-based
recombination-state classification, not a direct cM/Mb rate estimate. Its
stacked y-axis is the fraction of 5-kb windows per 500-kb bin in each class.
Descriptive quartet summaries by each subspecies-specific state are
in `stage2_quartet_by_recombination_state.tsv`. In these summaries, Q_T_ALT
enrichment is strongest for domesticus `RECENT_OR_OLDER` windows and for
very-recent/extensive states in castaneus and musculus, while no-recent states
are species-history leaning. This is a descriptive comparison among linked
windows and does not establish an association or causal effect.

No published coordinate list for the concatenated non-recombined-region tree
was recovered, so no source-defined block boundaries are drawn. The primary
coordinate system remains the 2017 chr17 window system; newer breakpoint
coordinates are not overlaid.

## Final House-mouse emphasis

The manuscript emphasis is now the arrangement-state-conditioned result in
`figures/house_mouse_t_complex_main_v4.png` and `.pdf`. The spatial q track is
retained as a supplementary/descriptive figure because its 250-kb means are
usually close to one third: strong 5-kb genealogical signals frequently
cancel across neighboring windows. The new main figure shows all eight exact
arrangement classes, their weights, their weighted contributions to the
all-tip quartet distribution, and the 1,000-replicate balanced control.

The eight-class contribution table is
`stage2_arrangement_pattern_contributions.tsv`. It reconstructs the frozen
ALL_TIPS means exactly. SSS and TTT favor Q_SPECIES, STS and TTS favor Q_T_ALT,
and TST favors Q_OTHER; STS gives the largest negative weighted contribution
to q_species minus q_t_alt, followed by TTS. This supports the restrained
claim that the t-complex has arrangement-state-dependent quartet distortion
whose mixed-state composition changes the aggregate ASTRAL summary topology.
It does not support a long near-fixed alternative-topology block.

Window heterogeneity and smoothing summaries are in
`stage2_window_topology_heterogeneity.tsv` and
`stage2_q_smoothing_sensitivity.tsv`; the optional summary figure is
`figures/house_mouse_t_complex_q_smoothing.png` and `.pdf`. Missing retained
window intervals are shaded and labeled in the spatial figures. These gaps
represent absent retained local-tree windows, not zero quartet support or no
recombination.
