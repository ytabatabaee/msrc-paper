# Main simulation figure plan

**Figure 3. Consequences of MSRC for quartet and species-tree inference**

The eventual manuscript figure will use a compact four-panel layout:

```text
┌──────────────────────┬──────────────────────┐
│ A. Quartet           │ B. epsilon threshold │
│ scenario + heatmap   │                      │
├──────────────────────┼──────────────────────┤
│ C. Larger tree       │ D. linkage           │
│ ASTRAL               │ pseudo-replication   │
└──────────────────────┴──────────────────────┘
```

Panel A: controlled conditional quartet, using the biological species-tree /
recombination-suppressed-region schematic, four representative exact/Monte
Carlo quartet-curve panels, four quartet-simplex panels, and the complete
exact `q_ALT-q_SPECIES` heatmap in
`quartet/analysis/figures/conditional_quartet_main_panel.{png,pdf}`.

Panel A is now supported by the completed Experiment 1 conditional figure and
the separate Experiment 2 mechanistic pilot figures. Experiment 2 keeps its
mechanistic scenario, unconditional probability, simplex, and bias-map figures
under `quartet/mechanistic/analysis/figures/`; these are diagnostic and are not
yet substituted into the final four-panel manuscript composition.

Panel B: affected-locus fraction / epsilon threshold, showing expected
aggregate quartet support and the topology-flip threshold `epsilon*`.

Panel C: larger-tree species-tree inference, using an approximately 30-taxon
or larger simulation and ASTRAL recovery/error versus MSRC strength or
affected fraction.

Panel D: linked-window / pseudo-replication behavior, comparing dense windows
with block-aware weighting or thinning and reporting topology/support effects.

Panels B–D are planning targets only and are not implemented by Experiment 1.
