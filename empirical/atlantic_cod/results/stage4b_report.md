# Atlantic cod Stage-4B report

Stage 4B tests whether the Stage-4A arrangement-concordance signal is unusually aligned with the independently frozen inversion intervals under spatial null models. It does not fit MSRC parameters or alter frozen definitions.

## Primary spatial test

The primary statistic is `Delta_D = mean(D_inside) - mean(D_outside)`, using fully inside versus fully outside windows and excluding partial boundary-overlap windows. The primary null exactly circularly shifts the ordered D track while keeping the inversion mask fixed.

| LG | observed Delta_D | p_raw | p_BH | rank | null range | minimum possible p |
|---|---:|---:|---:|---:|---|---:|
| LG01 | 1.232766 | 0.009174 | 0.021978 | 1 | -0.803166 to 1.232766 | 0.009174 |
| LG02 | 0.892332 | 0.010989 | 0.021978 | 1 | -0.444651 to 0.892332 | 0.010989 |
| LG07 | -0.005981 | 0.555556 | 0.555556 | 65 | -0.085606 to 0.081730 | 0.008547 |
| LG12 | 1.043170 | 0.038462 | 0.051282 | 4 | -1.101854 to 1.057598 | 0.009615 |

## Effect decomposition

| LG | Delta_A | p_A | Delta_B | p_B |
|---|---:|---:|---:|---:|
| LG01 | 0.741071 | 0.009174 | -0.491694 | 1.000000 |
| LG02 | 0.580520 | 0.010989 | -0.311812 | 1.000000 |
| LG07 | -0.026797 | 0.658120 | -0.020815 | 0.735043 |
| LG12 | 0.641093 | 0.019231 | -0.402076 | 0.971154 |

## Secondary block-placement null

| LG | observed Delta_D | p_block | rank | null range |
|---|---:|---:|---:|---|
| LG01 | 1.232766 | 0.022727 | 1 | -0.393030 to 1.232766 |
| LG02 | 0.892332 | 0.014085 | 1 | -0.444651 to 0.892332 |
| LG07 | -0.005981 | 0.621951 | 51 | -0.061346 to 0.081730 |
| LG12 | 1.043170 | 0.038462 | 2 | -1.082617 to 1.052789 |

## Boundary sensitivity

| LG | boundary treatment | observed Delta_D | p | rank |
|---|---|---:|---:|---:|
| LG01 | exclude | 1.232766 | 0.009174 | 1 |
| LG01 | inside | 1.232766 | 0.009091 | 1 |
| LG01 | outside | 1.204748 | 0.018182 | 1 |
| LG02 | exclude | 0.892332 | 0.010989 | 1 |
| LG02 | inside | 0.895986 | 0.010870 | 1 |
| LG02 | outside | 0.878631 | 0.021739 | 1 |
| LG07 | exclude | -0.005981 | 0.555556 | 65 |
| LG07 | inside | -0.006046 | 0.512605 | 61 |
| LG07 | outside | -0.005807 | 0.537815 | 64 |
| LG12 | exclude | 1.043170 | 0.038462 | 4 |
| LG12 | inside | 1.023625 | 0.028571 | 3 |
| LG12 | outside | 1.043405 | 0.038095 | 4 |

## Leave-one-population-out robustness

| LG | finite Delta_D range | finite difference from full range | finite sign always positive | undefined exclusions |
|---|---:|---:|---|---:|
| LG01 | 0.955065 to 1.406977 | -0.277701 to 0.174211 | yes | 0 |
| LG02 | 0.807029 to 1.067347 | -0.085302 to 0.175015 | yes | 0 |
| LG07 | -0.361797 to 0.358515 | -0.355816 to 0.364496 | no | 0 |
| LG12 | 1.005133 to 1.081206 | -0.038036 to 0.038036 | yes | 2 |

## Spatial topology structure

| LG | unique topologies inside | unique topologies outside | most common topology | frequency | max run | transitions |
|---|---:|---:|---|---:|---:|---:|
| LG01 | 61 | 43 | topology_0196 | 2 | 1 | 109 |
| LG02 | 21 | 70 | topology_0253 | 1 | 1 | 91 |
| LG07 | 36 | 81 | topology_0360 | 1 | 1 | 118 |
| LG12 | 53 | 51 | topology_0410 | 1 | 1 | 104 |

## LG07 contrast

LG07 remains fully reported and shows no corresponding positive spatial association under the primary circular-shift statistic. Stage 4B does not choose among biological explanations for that contrast.

## Global secondary test

The equal-weight four-LG statistic was 0.790571; the seeded Monte Carlo circular-shift summary used seed 20260929 with 100000 draws and gave p = 0.000010. This is secondary to the per-LG tests.

## Limitations

- Only four genomic supergene regions are tested.
- Local windows are linked along chromosomes.
- Quartet observations within a tree are highly dependent.
- This is multi-population, not a direct multi-species persistence test.
- Published MCC trees summarize posterior genealogy uncertainty.
- Spatial association alone does not establish a causal inversion effect.

## Interpretation

Where supported, local genealogies within the inversion are unusually aligned with the arrangement partition relative to other positions on the same linkage group. This is not a claim that MSRC is proven or that the inversions caused the genealogies.
