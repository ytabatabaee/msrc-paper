# Simulation studies

All MSRC simulation is performed with the authoritative
[`ytabatabaee/msrc-sim`](https://github.com/ytabatabaee/msrc-sim) package. This
directory contains experiment definitions, frozen configurations, simulator
outputs, processed tables, analysis scripts, and figures for the manuscript.

The planned progression is:

1. Conditional exact quartet behavior.
2. Full mechanistic quartet simulations under Wright–Fisher.
3. Wright–Fisher versus Moran robustness.
4. Affected-locus fraction and the epsilon failure threshold.
5. Larger species trees with ASTRAL.
6. Linked-window and pseudo-replication analysis.

The first experiment is under `quartet/` and uses the conditional simulator
mode. `msrc-paper` does not vendor or duplicate the simulator implementation.
