# Simulator validation

The conditional quartet study used `ytabatabaee/msrc-sim` at commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9` (package version 0.8.8), resolved
from the local checkout recorded in
`datasets/simulator_provenance.json`. The public `msrc-sim --config` CLI was
used for every raw cell; no simulator source was copied into this repository.

The complete upstream test collection could not be collected because the
existing checkout has an import error in
`tests/test_daughter_sorting_eta_analysis.py`: `scripts.run_daughter_sorting_eta`
is not importable. I did not alter that checkout. With that single unrelated
module excluded, the remaining upstream suite passed: 176 tests passed with
8 warnings. The conditional analytic and CLI tests passed in that run.

For the 70-cell production grid and three topologies per cell (210
comparisons), the maximum absolute empirical-versus-exact error was
0.0040061, the median was 0.0006986, 97.619% were within 2 Monte Carlo SE,
and 100% were within 3 SE. The exact `H_m(t)` calculation is therefore used
as the authoritative source for the dominance boundary.

The raw CLI summary contains a legacy `version` field of `0.4.1`, while the
installed package metadata and `msrcsim.__version__` are `0.8.8`. Provenance
uses the package version and exact commit; the discrepancy is retained as
observed output and does not affect the conditional probabilities.
