# House-mouse t-complex: arrangement-state-conditioned quartet distortion

The primary dataset comprised 4,046 published maximum-likelihood gene trees
from non-overlapping 5-kb windows across chr17. The standard-only
four-subspecies ASTRAL analysis recovered `Q_SPECIES`, whereas the analysis
including standard and reconstructed pseudo-t haplotypes recovered `Q_T_ALT`.
The focal ASTRAL quartet annotations were q_species = 0.369910, q_t_alt =
0.343656, and q_other = 0.286434 for standard-only, compared with q_species =
0.323323, q_t_alt = 0.365364, and q_other = 0.311313 for all tips. These named
biological q values are mapped to the inferred focal split and are not fixed
aliases for ASTRAL q1, q2, and q3.

The local spatial signal was heterogeneous rather than a single coherent
alternative-topology block. Across ALL_TIPS windows, the mean frequencies were
q_species = 0.323187, q_t_alt = 0.365290, and q_other = 0.311230, all close to
one third. Nevertheless, 33.9% of individual windows were Q_SPECIES dominant,
37.3% were Q_T_ALT dominant, and 28.4% were Q_OTHER dominant, with 0.4% ties
or unresolved winners. The absolute species-versus-alternative contrast was at
least 0.25 in 49.8% of windows, at least 0.50 in 21.8%, and at least 0.75 in
7.7%. Mean absolute contrast declined from 0.311 at 5-kb resolution to 0.094
at 500-kb aggregation, demonstrating cancellation among neighboring linked
windows rather than absence of strong local genealogies.

The stronger signal emerged after conditioning on arrangement state. The
homogeneous SSS and TTT classes favored Q_SPECIES, whereas STS favored Q_T_ALT
(0.500) and TTS strongly favored Q_T_ALT (0.707). TST favored Q_OTHER (0.694).
The remaining mixed classes also contributed materially: SST favored Q_OTHER,
STT strongly favored Q_SPECIES, and TSS modestly favored Q_SPECIES. Thus the
ALL_TIPS distribution is not an average of only SSS and TTT; six mixed-state
classes are included. The exact pattern weights were SSS 0.325, SST 0.139,
STS 0.200, STT 0.086, TSS 0.108, TST 0.046, TTS 0.067, and TTT 0.029.
Their weighted contributions reconstruct the all-tip quartet frequencies; the
largest negative contribution to q_species − q_t_alt came from STS (-0.060),
followed by TTS (-0.038).

The mixture effect persisted under deterministic balanced resampling. The
matched standard treatment returned Q_SPECIES in 53.7% of replicates, the
matched t treatment returned Q_SPECIES in 100%, and the mixed treatment
returned Q_T_ALT in 65.8%. The mixed treatment was composition-balanced but
not total-sample-size matched to the two single-state treatments. These
results support arrangement-state-dependent quartet distortion and a changed
aggregate summary topology, while the linked windows remain non-independent
and the analysis has no independent genome-wide local-tree control.

