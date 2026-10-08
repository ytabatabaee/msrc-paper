# Cross-dataset quartet-mixture CU distortion

Adding a localized genealogy regime changes concordant quartet support while the pooled species-tree topology remains unchanged in both datasets. The linear mixture exactly reconstructs the observed pooled support from the frozen background and affected-region scores, and the MSC quartet-to-CU transform predicts the direction and approximate magnitude of the observed ASTRAL4 CULength change.

The fire-ant fixed-background analysis provides two strong shortening examples and a small lengthening example. These results support the general branch-length-distortion theorem and do not require an exact 2:2 arrangement configuration. They are a theory-guided empirical demonstration, not independent fitting or validation of all MSRC parameters, because the affected-region support `q_R` is observed rather than predicted from independently estimated MSRC parameters. The source fire-ant study interprets its pattern as adaptive introgression; MSRC without gene flow is not claimed here.

Maximum absolute pooled-support prediction error: `8.750e-07`.

All seven branch transforms were defined (`q_bg` and `q_mix` > 1/3); no branch required an undefined MSC length placeholder.

| dataset | branch | q_bg | q_R | q_pred | q_obs | CU_bg | CU_all | predicted % change | observed % change |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Anopheles 2La | 5 taxa / melas + merus | 0.831347 | 0.366710 | 0.626640 | 0.626641 | 1.36548 | 0.57803 | -57.82% | -57.67% |
| Anopheles 2La | 3 taxa / melas + merus + quadriannulatus | 0.957005 | 0.814320 | 0.894142 | 0.894142 | 2.70133 | 1.83158 | -32.87% | -32.20% |
| Anopheles 2La | 4 taxa / coluzzii + gambiae | 0.871768 | 0.402447 | 0.664997 | 0.664997 | 1.63610 | 0.68612 | -58.25% | -58.06% |
| Fire ants chr16 | richteri SB/Sb | 0.967056 | 0.245955 | 0.791013 | 0.791012 | 2.84094 | 1.14248 | -61.43% | -59.79% |
| Fire ants chr16 | invicta/macdonaghi SB/Sb | 0.882244 | 0.331903 | 0.747888 | 0.747888 | 1.68846 | 0.95865 | -43.91% | -43.22% |
| Fire ants chr16 | four focal groups / outgroups | 0.992805 | 1.000000 | 0.994562 | 0.994562 | 3.91279 | 4.19117 | 6.18% | 7.11% |
| Fire ants chr16 | geminata / pusillignis | 0.963298 | 0.942308 | 0.958174 | 0.958173 | 2.74929 | 2.66706 | -4.51% | -2.99% |
