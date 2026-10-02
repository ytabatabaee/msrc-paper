# Stage 1B Synthetic Validation Report

`SYNTHETIC_ONLY = TRUE`

No real Anopheles sequence-derived topology, local-tree, SNP, haplotype, q1/q2/q3, QQS/BQS, or NJ-tree data were read.

Freeze checksum: `64c6b70f9046ea5b3b7f493150dfd68cd6ae454beb513c3147599e4d204a4274`

| Scenario | Classification | Delta | Block-aware p | Interpretation |
| --- | --- | ---: | ---: | --- |
| null_msc_like | not_positive | -0.014 | 0.760 | not systematically positive expected |
| msrc_positive | positive | 0.422 | 0.005 | positive control should be detected |
| broad_alternative | not_positive | 0.014 | 0.240 | broad alternative should be flagged as not 2La-specific |

The prospective test is whether the independently predicted arrangement topology is enriched inside 2La relative to fixed flanks.
The Design-B validation asks whether replacing arrangement background changes topology in the structurally predicted direction while species identity is controlled.
