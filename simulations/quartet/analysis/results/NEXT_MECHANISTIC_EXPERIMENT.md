# Next mechanistic experiment

The next study should use the existing mechanistic mode in `ytabatabaee/msrc-sim`
to place the species topology `12|34` in a four-taxon species tree and compare
it with the arrangement-associated topology `13|24`. The conditional result
indicates that the first mechanistic grid should concentrate on forward
rearrangement histories whose derived arrangement places taxa 1 and 3
together, then vary persistence across the relevant speciation boundaries and
the effective cross-arrangement recombination fraction.

The design should cross a small set of forward rearrangement origins and
initial copy counts with persistence histories, recombination suppression,
and the simulator's Wright–Fisher population process. A matched Moran version
should follow after the WF pilot has stable diagnostics. Each cell should
report sampled arrangement state, quartet counts, species-topology support,
arrangement-topology support, affected-locus fraction, and history metadata.

Before execution, inspect the current mechanistic config schema and choose
parameters that preserve the conditional experiment's duration and switching
interpretation where the mechanistic API permits. Run smoke configurations,
validate raw summaries against the public simulator outputs, and keep the
mechanistic grid separate from this completed conditional dataset. Do not run
this experiment as part of Experiment 1.
