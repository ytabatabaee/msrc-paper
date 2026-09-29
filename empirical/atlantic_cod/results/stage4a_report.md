# Atlantic cod Stage-4A report

Stage 4A compares published local 250-kb population trees with the independently frozen arrangement and baseline quartet predictions. It is descriptive only; no spatial-null or permutation test was run.

## Input verification

- Stage-2 structural candidate checksum: `303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c`.
- Stage-3 baseline quartet checksum: `5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1`.
- Local-window table checksum: `26ed195ee26cb860ca0885c76d768d427639a4338b74a787d9421498e10dfb0b`.

## Primary result

`D = fraction_arrangement - fraction_baseline`, computed per window using informative quartets only.

| LG | mean D inside | mean D outside | inside - outside | median D inside | median D outside | fraction positive informative quartets |
|---|---:|---:|---:|---:|---:|---:|
| LG01 | 1.000000 | -0.232766 | 1.232766 | 1.000000 | -0.250000 | 1.000000 |
| LG02 | 0.919604 | 0.027273 | 0.892331 | 1.000000 | -0.019481 | 1.000000 |
| LG07 | 0.046942 | 0.052924 | -0.005981 | 0.015504 | 0.007752 | 0.465116 |
| LG12 | 0.930425 | -0.112745 | 1.043170 | 1.000000 | -0.062500 | 1.000000 |

Inside inversion windows show greater alignment with the independently frozen arrangement partition than outside windows for the LGs with positive inside-minus-outside differences. This is descriptive and not a significance claim.

## Topology composition

| LG | region | windows | mean arrangement fraction | mean baseline fraction | mean third fraction |
|---|---|---:|---:|---:|---:|
| LG01 | inside | 66 | 1.000000 | 0.000000 | 0.000000 |
| LG01 | outside | 43 | 0.258929 | 0.491694 | 0.249377 |
| LG02 | inside | 21 | 0.946197 | 0.026592 | 0.027211 |
| LG02 | outside | 70 | 0.365677 | 0.338404 | 0.295918 |
| LG07 | inside | 36 | 0.342808 | 0.295866 | 0.361326 |
| LG07 | outside | 81 | 0.369605 | 0.316681 | 0.313714 |
| LG12 | inside | 53 | 0.964623 | 0.034198 | 0.001179 |
| LG12 | outside | 51 | 0.323529 | 0.436275 | 0.240196 |

## Quartet consistency

| LG | positive informative quartets | fraction positive | median effect | min effect | max effect |
|---|---:|---:|---:|---:|---:|
| LG01 | 112/112 | 1.000000 | 1.058140 | 0.651163 | 1.930233 |
| LG02 | 77/77 | 1.000000 | 0.895238 | 0.509524 | 1.338095 |
| LG07 | 60/129 | 0.465116 | -0.024691 | -1.061728 | 1.074074 |
| LG12 | 16/16 | 1.000000 | 1.041435 | 0.924528 | 1.198298 |

## Spatial pattern

Positive inside-minus-outside D values were observed for 3 of 4 LGs in the descriptive window-level summary. The topology-track PDFs shade the frozen inversion intervals and mark boundary-overlap windows; visual alignment should be assessed from those tracks before Stage 4B.

## Bornholm positive control

Reported LG12 local-switch behavior recovered: yes. Details are in `results/stage4a_bornholm_check.md`.

## Boundary sensitivity

| LG | boundary treatment | delta D | mean D inside | mean D outside |
|---|---|---:|---:|---:|
| LG01 | boundary_excluded | 1.232766 | 1.000000 | -0.232766 |
| LG01 | boundary_as_inside | 1.232766 | 1.000000 | -0.232766 |
| LG01 | boundary_as_outside | 1.204748 | 1.000000 | -0.204748 |
| LG02 | boundary_excluded | 0.892331 | 0.919604 | 0.027273 |
| LG02 | boundary_as_inside | 0.895986 | 0.923259 | 0.027273 |
| LG02 | boundary_as_outside | 0.878631 | 0.919604 | 0.040973 |
| LG07 | boundary_excluded | -0.005981 | 0.046942 | 0.052924 |
| LG07 | boundary_as_inside | -0.006046 | 0.046877 | 0.052924 |
| LG07 | boundary_as_outside | -0.005808 | 0.046942 | 0.052750 |
| LG12 | boundary_excluded | 1.043170 | 0.930425 | -0.112745 |
| LG12 | boundary_as_inside | 1.023625 | 0.910880 | -0.112745 |
| LG12 | boundary_as_outside | 1.043405 | 0.930425 | -0.112981 |

The primary comparison excludes partial boundary-overlap windows from both inside and outside classes. The sensitivity table records alternative treatments without selecting a favorable definition.

## Caveats

- Adjacent 250-kb windows are spatially linked, so Stage 4A does not use iid tests.
- Multiple quartets from the same local tree are dependent; the primary unit is the genomic window.
- LG12 has only 16 informative quartets and is less powered than LG01, LG02, and LG07.
- LG12 ancestral/derived orientation follows the weaker Stage-2 demographic-inference evidence, unlike the stronger direct outgroup evidence for LG01, LG02, and LG07.
- This is a multi-population Atlantic cod analysis, not a direct multi-species persistence-through-speciation test.

## Parsing issues

No unexpected topology/parsing issues were encountered.

## Stage boundary

No circular-shift tests, block permutations, P-values, formal MSRC model fitting, or SNAPP/BEAST reruns were performed.
