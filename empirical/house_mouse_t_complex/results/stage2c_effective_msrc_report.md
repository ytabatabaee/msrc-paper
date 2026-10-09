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

## Spatially blocked validation

Blocks are nonoverlapping physical intervals aligned to the first retained 5-kb window at 5 Mb. The primary block size is 1 Mb because it is large relative to the 5-kb window spacing while retaining multiple spatial folds; 500-kb and 2-Mb runs are fixed sensitivity analyses, not choices optimized for accuracy. Complete blocks, never individual windows, are assigned to the held-out fold. For each held-out block and configuration, Delta is fit from the other two configurations in all remaining blocks with equal configuration weight.

| block size | configuration | folds | mean max abs error | mean L1 | mean arrangement-support error | mean non-arr symmetry error |
|---:|---|---:|---:|---:|---:|---:|
| 500000 | STT | 58 | 0.171997 | 0.343994 | -0.000030 | 0.077891 |
| 500000 | TST | 58 | 0.173993 | 0.347986 | +0.006362 | 0.088293 |
| 500000 | TTS | 58 | 0.171216 | 0.342432 | -0.006331 | 0.076548 |
| 1000000 | STT | 32 | 0.165921 | 0.331843 | -0.006585 | 0.069582 |
| 1000000 | TST | 32 | 0.171172 | 0.342344 | +0.003340 | 0.085926 |
| 1000000 | TTS | 32 | 0.166644 | 0.333289 | +0.003246 | 0.071804 |
| 2000000 | STT | 17 | 0.157055 | 0.314110 | +0.004963 | 0.060110 |
| 2000000 | TST | 17 | 0.161400 | 0.322800 | +0.005034 | 0.074765 |
| 2000000 | TTS | 17 | 0.167464 | 0.334929 | -0.009997 | 0.043693 |

No conventional multinomial p-values are reported: induced quartets reuse individuals, windows, and linked genomic segments, so they are not independent replicates. The block-CV results are the preferred robustness analysis.

The held-out validation is a low-dimensional effective MSRC demonstration. It does not establish that the mouse affected genealogy is generated solely by MSRC, and it does not separate the biological switching time and migration parameters.

Frozen tip metadata check: 8 Mus spretus tips are `outgroup_not_t_haplotype`.
