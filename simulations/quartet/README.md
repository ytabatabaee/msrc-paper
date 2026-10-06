# Experiment 1: conditional four-taxon quartet

This experiment fixes arrangement configuration `1010`, with taxa 1 and 3 in
arrangement state 1 and taxa 2 and 4 in state 0. With msrc-sim's taxon order,
the topology ordering is `q1 = 12|34`, `q2 = 13|24`, and `q3 = 14|23`.

The experiment varies structured-interval duration and symmetric switching
`m01 = m10 = m`, using the exact finite-state `H_m(t)` probabilities and
Monte Carlo draws from `ytabatabaee/msrc-sim`. It measures conditional
alternative-topology dominance for the fixed arrangement state. It does not
yet simulate the full species-tree plus forward rearrangement history, so its
results are not interpreted as a full MSRC species-tree inconsistency result.

The `datasets/` tree contains frozen configs, raw simulator outputs, and
normalized tables. The `analysis/` tree contains drivers, validation results,
figures, and tests. The mechanistic WF/Moran experiment is planned but is not
run by this study.
