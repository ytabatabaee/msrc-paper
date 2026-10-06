# Conditional quartet experiment results

The fixed configuration `1010` places taxa 1 and 3 in arrangement state 1,
and taxa 2 and 4 in state 0. In the msrc-sim quartet ordering, `q1` is
`12|34` (`Q_SPECIES`), `q2` is `13|24` (`Q_ALT`), and `q3` is `14|23`
(`Q_OTHER`). Thus the arrangement-associated topology is `Q_ALT`.

Across the exact grid, lowering symmetric switching preserves the initial
arrangement contrast for longer. The exact alternative-topology probability
exceeds the species-topology probability at every tested `m` from 0 through
1 for every duration from 0.02 through 2.0. The exact difference
`Q_ALT - Q_SPECIES` decreases monotonically with `m` at every duration in
this grid. The strongest bias occurs at duration 2.0 and `m=0`, where
`(Q_SPECIES, Q_ALT, Q_OTHER) = (0.00610521, 0.98778957, 0.00610521)` and
the exact difference is 0.98168436.

At `m=0`, arrangement backgrounds are isolated during the interval. The
alternative topology is favored even for the shortest duration and becomes
increasingly dominant as the interval lengthens. At `m=1`, switching reduces
the contrast; for example, at duration 2.0 the exact vector is approximately
`(0.222223, 0.555554, 0.222223)`, so the alternative remains the largest
topology on the tested grid.

Monte Carlo estimates from 100,000 loci per cell agree with the exact
finite-state probabilities: the maximum absolute error is 0.0040061 and all
210 topology comparisons are within 3 Monte Carlo SE. The exact calculation
is used for the heatmap, boundary, and smooth curves.

This is an arrangement-conditioned quartet bias in a conditional
structured-interval experiment. It is a mechanistic validation of the
conditional model, not yet a full species-tree inconsistency or full MSRC
anomalous gene-tree result. Those stronger claims require the mechanistic
species-tree simulations planned next.
