# Cross-dataset theory–empirical bridge

The quartet mixture predicts two linked consequences. When the affected fraction is below the first valid crossing threshold ε*, the pooled topology can remain stable while the concordant support and apparent CU branch length change. When ε exceeds ε*, an alternative topology overtakes the background topology.

Anopheles 2La uses 430 inside windows among 976 usable non-boundary windows (ε = 0.440574). Both included topology-sensitive branches are below threshold, consistent with `T_outside6` versus `T_all6` RF = 0, while the inside-only topology differs from the outside topology (RF = 4). Fire ants use 52 supergene windows among 213 windows (ε = 0.244131). Both focal branches are below threshold; free-topology Stage 6C gives RF = 0 for background versus all and RF = 4 for supergene versus background.

House mouse defines ε as the weight of mixed arrangement-state quartets SST, STS, STT, TSS, TST, and TTS within the t-complex, rather than as an inversion-wide fraction. The observed mixed weight is ε = 0.646429, above ε* = 0.311776. The mixture reconstructs the frozen ALL_TIPS vector and agrees with the observed ASTRAL switch from Q_SPECIES (`T0_STANDARD`) to Q_T_ALT (`T1_ALL_WINDOWS`). SSS and TTT are both Q_SPECIES-favoring; the switch is driven by the composition of mixed classes.

This figure is a theory-guided empirical demonstration, not full parameter-level MSRC validation: the affected-component quartet distributions are measured from the data rather than predicted from independently estimated MSRC parameters. The fire-ant source study’s adaptive-introgression interpretation is retained, and Anopheles has known introgression; neither result establishes a rearrangement-only mechanism.

| dataset | branch/configuration | ε_obs | ε*_alt1 | ε*_alt2 | ε* | relation | background | pooled | affected-only |
|---|---|---:|---:|---:|---:|---|---|---|---|
| Anopheles 2La | 5 taxa | melas + merus | 0.440574 | 0.873778 | NA | 0.873778 | below | Q_SPECIES | Q_SPECIES | Q_ALT |
| Anopheles 2La | 4 taxa | coluzzii + gambiae | 0.440574 | NA | 0.949162 | 0.949162 | below | Q_SPECIES | Q_SPECIES | Q_OTHER |
| Fire ants chr16 | richteri SB/Sb | 0.244131 | NA | 0.764866 | 0.764866 | below | Q_SPECIES | Q_SPECIES | Q_OTHER |
| Fire ants chr16 | invicta/macdonaghi SB/Sb | 0.244131 | 0.793689 | NA | 0.793689 | below | Q_SPECIES | Q_SPECIES | Q_ALT |
| House mouse chr17 t-complex | mixed S/T quartets | 0.646429 | 0.311776 | 0.744071 | 0.311776 | above | Q_SPECIES | Q_T_ALT | Q_T_ALT |

Mouse reconstructed ALL_TIPS q = (0.3231866460, 0.3652898935, 0.3112302550); maximum absolute error = 4.047e-11.

Excluded branches:
- Anopheles 2La: arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus — neither affected alternative exceeds r_S
- Fire ants chr16: focal_four_vs_outgroups — not one of the topology-sensitive focal branches requested for Panel B
- Fire ants chr16: focal_related_internal_branch — not one of the topology-sensitive focal branches requested for Panel B
