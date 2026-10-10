# Stage 2C effective MSRC validation

The symmetric 2:2 MSRC quartet result has one identifiable effective parameter, `Delta`, with `q_arr = (1 + 2 Delta)/3` and each non-arrangement topology `(1 - Delta)/3`. The empirical validation estimates this composite parameter only; it does not estimate biological `t` and `m` separately.

The three clean configurations rotate which biological topology is arrangement-concordant: STT maps to Q_SPECIES, TST maps to Q_OTHER, and TTS maps to Q_T_ALT. Mus spretus is treated as `outgroup_not_t_haplotype` according to the frozen tip metadata.

## Aggregate configuration estimates

| configuration | arrangement topology | q_arr | non-arr 1 | non-arr 2 | unresolved | Delta |
|---|---|---:|---:|---:|---:|---:|
| STT | q_species | 0.685600 | 0.157056 | 0.157344 | 0.000256 | 0.528400 |
| TST | q_other | 0.693686 | 0.153126 | 0.153188 | 0.000133 | 0.540529 |
| TTS | q_t_alt | 0.706824 | 0.162149 | 0.131027 | 0.000257 | 0.560236 |

## Held-out validation

The primary leave-one-configuration-out fit gives the two training configurations equal weight. STT-only predictions are reported separately as a simple demonstration. The two non-arrangement frequencies are expected to be equal under the symmetric model; their observed difference is retained as a model-departure diagnostic.

| analysis | held-out | training | Delta | max abs error | L1 | arrangement-support error | non-arrangement symmetry error |
|---|---|---|---:|---:|---:|---:|---:|
| leave_one_configuration_out | STT | TST,TTS | 0.550383 | 0.014655 | 0.029311 | +0.014655 | 0.000288 |
| leave_one_configuration_out | TST | STT,TTS | 0.544318 | 0.002526 | 0.005052 | +0.002526 | 0.000062 |
| leave_one_configuration_out | TTS | STT,TST | 0.534465 | 0.024151 | 0.048303 | -0.017181 | 0.031122 |
| STT_only_prediction | TST | STT | 0.528400 | 0.008086 | 0.016173 | -0.008086 | 0.000062 |
| STT_only_prediction | TTS | STT | 0.528400 | 0.026173 | 0.052346 | -0.021224 | 0.031122 |

## Linkage-aware block bootstrap

Physical blocks were used as resampling units to quantify uncertainty in the aggregate held-out predictions while preserving local linkage and spatial heterogeneity. Blocks are nonoverlapping intervals aligned to the first retained 5-kb window at 5 Mb. The primary analysis uses 1-Mb blocks and 10,000 bootstrap replicates; fixed 500-kb and 2-Mb analyses use 5,000 replicates each. Each replicate samples complete blocks with replacement until the original number of blocks is restored, aggregates quartet counts, and gives the two training configurations equal weight.

| block size | held-out | n_bootstrap | Delta train median [95% CI] | max error median [95% CI] | L1 median [95% CI] |
|---:|---|---:|---|---|---|
| 500000 | STT | 5000 | 0.550755 [0.473430, 0.624056] | 0.026465 [0.004084, 0.079261] | 0.052929 [0.008167, 0.158522] |
| 500000 | TST | 5000 | 0.543635 [0.471527, 0.613985] | 0.020660 [0.003356, 0.062906] | 0.041319 [0.006712, 0.125812] |
| 500000 | TTS | 5000 | 0.533889 [0.469556, 0.600370] | 0.025311 [0.007993, 0.047180] | 0.050622 [0.015986, 0.094360] |
| 1000000 | STT | 10000 | 0.552144 [0.451225, 0.643882] | 0.032117 [0.005209, 0.093991] | 0.064235 [0.010417, 0.187983] |
| 1000000 | TST | 10000 | 0.544728 [0.454984, 0.634164] | 0.025189 [0.004236, 0.076304] | 0.050378 [0.008471, 0.152608] |
| 1000000 | TTS | 10000 | 0.534338 [0.453902, 0.617635] | 0.025864 [0.006849, 0.050329] | 0.051727 [0.013697, 0.100657] |
| 2000000 | STT | 5000 | 0.552239 [0.420662, 0.661259] | 0.038247 [0.006174, 0.116231] | 0.076494 [0.012348, 0.232461] |
| 2000000 | TST | 5000 | 0.546815 [0.432921, 0.648316] | 0.032043 [0.005804, 0.100481] | 0.064087 [0.011607, 0.200963] |
| 2000000 | TTS | 5000 | 0.534712 [0.441681, 0.629061] | 0.026145 [0.005312, 0.054983] | 0.052290 [0.010623, 0.109966] |

## Local spatial heterogeneity diagnostic

A constant Delta is not intended to predict every linked physical block. For each 1-Mb block, `Delta_block = (3 q_arr - 1)/2` is retained without truncation; negative values therefore indicate local model departure rather than an invalid estimate. The bootstrap is the aggregate robustness analysis, while these block estimates describe local heterogeneity.

| configuration | Delta range | median | IQR | fraction < 0 | fraction > 1 |
|---|---|---:|---|---:|---:|
| STT | [-0.161288, 0.991667] | 0.643633 | [0.421454, 0.772743] | 0.094 | 0.000 |
| TST | [-0.314583, 1.000000] | 0.576655 | [0.384529, 0.815783] | 0.062 | 0.000 |
| TTS | [-0.006614, 1.000000] | 0.582563 | [0.398627, 0.756888] | 0.031 | 0.000 |

Spearman correlations of block Delta estimates: STT/TST = 0.2617; STT/TTS = 0.4883.

No conventional multinomial p-values are reported: induced quartets reuse individuals, windows, and linked genomic segments, so they are not independent replicates. The old leave-one-single-block-out table is retained for provenance, but it is not interpreted as the preferred robustness analysis.

The held-out validation is a low-dimensional effective MSRC demonstration. TTS has a modest non-arrangement asymmetry (~0.031), indicating departure from the simplest symmetric model. The result does not establish that the entire mouse genealogy is generated solely by MSRC, and it does not separate the biological switching time and migration parameters.

Frozen tip metadata check: 8 Mus spretus tips are `outgroup_not_t_haplotype`.
