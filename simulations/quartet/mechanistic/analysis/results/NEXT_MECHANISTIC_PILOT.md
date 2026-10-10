# Next mechanistic pilot

The first unconditional pilot found no supported positive
`P(Q_A)-P(Q_S)` and observed the terminal `1010` pattern in only six cells,
at frequencies no higher than 0.0833. The next run should therefore increase
independent Wright–Fisher histories and target persistence rather than simply
adding loci to the same histories.

Recommended design:

- retain matched ordinary MSC controls and the same balanced species tree;
- focus first on initial copy fractions 0.20 and 0.50, both ROOT origin depths,
  `τ ∈ {0.5, 2}`, and cross-arrangement fractions `{1, 0.1}`;
- expand to at least 100 independent histories per cell before interpreting
  a positive point estimate;
- record terminal 1010 and all 2:2 pattern frequencies before selecting any
  conditional diagnostic;
- refine only cells whose independent-history interval is close to or above
  zero, with a pre-specified correction for searching across cells.

If the relevant terminal patterns remain rare, report that mechanistic
averaging suppresses the conditional Experiment 1 effect in this parameter
range and investigate persistence mechanisms already supported by `msrc-sim`.
Do not impose `1010` in the primary run or alter the simulator to create a
topology reversal.
