# Stage7b MSRC coalescent validation

The predeclared theorem-linked direction is `Delta_coal > 0`: opposite-arrangement pairs have extra residual coalescence depth inside the supergene after pair-specific collinear baselines are removed.

| LG | Delta_coal | circular p | BH p | m_eff |
|---|---:|---:|---:|---:|
| LG01 | 0.75703 | 0.00917431 | 0.010989 | 0.660476113 |
| LG02 | 0.653994 | 0.010989 | 0.010989 | 0.764532738 |
| LG07 | 0.435501 | 0.00854701 | 0.010989 | 1.14810298 |
| LG12 | 0.324552 | 0.00961538 | 0.010989 | 1.54058334 |

m_eff is only the effective long-duration approximation 1/(2 Delta_coal), in inverse published SNAPP time units. It is not a recombination rate and does not identify biological m.

The circular null is chromosome-aware; the jackknife uses populations as dependence units and does not treat population pairs as iid. LG07 is the key example of near-zero topology shift with positive coalescence-time shift. This is a theory-consistent empirical signature, not proof that MSRC is the sole historical process.

## Arrangement-class pair summaries

| LG | class | n pairs | mean inside residual | median inside residual | IQR | fraction positive |
|---|---|---:|---:|---:|---:|---:|
| LG01 | opposite_arrangement | 32 | 0.777633 | 0.78357 | 0.0250635 | 1 |
| LG01 | same_arrangement | 34 | 0.0208625 | 0.00778923 | 0.0498156 | 0.764706 |
| LG02 | opposite_arrangement | 32 | 0.650137 | 0.655868 | 0.0192036 | 1 |
| LG02 | same_arrangement | 34 | -0.00373785 | -0.00726274 | 0.00916079 | 0.205882 |
| LG07 | opposite_arrangement | 35 | 0.794454 | 1.10478 | 1.1222 | 0.771429 |
| LG07 | same_arrangement | 31 | 0.358949 | 6.08793e-05 | 1.10972 | 0.516129 |
| LG12 | opposite_arrangement | 27 | 0.317486 | 0.323574 | 0.0244218 | 1 |
| LG12 | same_arrangement | 39 | -0.00740564 | -0.00625993 | 0.0231868 | 0.307692 |

The mean-outside baseline sensitivity preserved the sign and qualitative LG ordering of Delta_coal; because the same pair set contributes to every eligible window, subtracting a different constant baseline per pair changes absolute residual levels but leaves this inside-minus-outside contrast unchanged up to rounding.

| LG | Delta_D | topology p | Delta_coal | coalescent p | BH p | class |
|---|---:|---:|---:|---:|---:|---|
| LG01 | 1.23277 | 0.009174 | 0.75703 | 0.00917431 | 0.010989 | topology_and_time_shift |
| LG02 | 0.892332 | 0.010989 | 0.653994 | 0.010989 | 0.010989 | topology_and_time_shift |
| LG07 | -0.005981 | 0.555556 | 0.435501 | 0.00854701 | 0.010989 | time_shift_without_topology_shift |
| LG12 | 1.04317 | 0.038462 | 0.324552 | 0.00961538 | 0.010989 | topology_and_time_shift |

The data do not separately identify m, lambda, and tau under the finite-duration model. The analysis uses published per-window SNAPP MCC point estimates because posterior samples are unavailable.
