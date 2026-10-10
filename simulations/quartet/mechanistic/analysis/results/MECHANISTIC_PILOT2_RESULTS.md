# Mechanistic quartet Experiment 2B

## Design and validation

Experiment 2B used `ytabatabaee/msrc-sim` commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9` (version 0.8.8) through the public
`msrc-sim-replicates` CLI. The simulator source was not modified. Experiment
2A outputs remain intact in their original `pilot` directories; this run uses
separate `pilot2` configuration and raw-output directories.

The primary design contained all 32 combinations of:

```text
Ne                 20, 100
daughter length τ  0.5, 2.0 coalescent units
ROOT origin        0.25, 0.75 of the ROOT extension
initial frequency  0.20, 0.50
cross fraction     0.1, 1.0
```

Four additional short-branch controls used `τ=0.1`, both population sizes,
both origin depths, selected intermediate/high initial frequencies, and strong
suppression (`cross fraction=0.1`). Every cell used 100 independent
Wright–Fisher histories and 20 loci per shared history. A 36-cell smoke run
used three histories and ten loci per history before the primary run.

The 36 completed pilot2 simulator jobs consumed approximately 178.9 minutes
of recorded subprocess time across four concurrent workers.

The species tree was the balanced tree used in Experiment 2A. Each daughter
ancestral branch has length `2Neτ` generations. Direct inspection of
`SpeciesTree` and the matched MSC controls confirms that the comparable
unrooted quartet length is `t=2τ`: the two daughter ancestral branches each
contribute `τ` coalescent units to the quartet sorting interval. The control
formula is therefore

```text
P(Q_S) = 1 - 2/3 exp(-2τ)
P(Q_A) = P(Q_O) = 1/3 exp(-2τ).
```

The previous audit contained one stale statement using `g/(2Ne)` for this
balanced-tree control; it has been corrected to `g/Ne = 2τ`.

The replicate API records one history per row, terminal arrangement pattern,
ROOT end status, daughter-branch persistence through `persisted_to_target`,
and quartet counts from loci sharing that history. It does not export one
Newick gene tree per locus; the analysis therefore treats the saved topology
counts as the supported public output and does not claim per-locus tree export.
It also does not expose numeric arrangement frequencies at every speciation
boundary in the replicate summary. Experiment 2B therefore reports the
supported branch-end persistence flags and terminal samples; obtaining full
boundary-frequency trajectories would require a future simulator-side export,
which was not added here.

## Persistence

Persistence was measured in two ways:

1. `root_persistence`: the ROOT branch ended segregating.
2. `daughter_AB_persistence`: both daughter ancestral branches `A` and `B`
   ended segregating, using the simulator’s `persistence_target_branches`
   field.

Terminal 2:2 sampling was reported separately and was not treated as a
persistence proxy. Across the 36 cells, daughter-branch persistence ranged
from 0 to 0.95, while terminal 2:2 frequency ranged from 0.08 to 0.31.
The highest daughter-branch persistence was 0.95 in the short-branch control
`Ne=100, τ=0.1, origin=.75, initial frequency=.5, cross=.1`. The highest
terminal 2:2 frequency was 0.31 at `Ne=100, τ=2, origin=.75, initial
frequency=.5, cross=.1`.

The relevant terminal pattern remained uncommon: `1010` reached at most 0.04
and `0101` at most 0.03 in any cell. Across all 3,600 histories, the
non-overlapping decomposition contained 776 lost histories, 157 fixed
histories, 24 `1010` histories, 17 `0101` histories, 604 other 2:2 histories,
and 2,022 segregating non-2:2 histories.

## Bias decomposition

Each history was assigned exactly one category using this precedence:

```text
lost → fixed → terminal_1010 → terminal_0101 → other_2_2
      → segregating_non_2_2 → remaining
```

Thus `1010` and `0101` are reported separately while the `other_2_2` and
combined terminal 2:2 summaries make their unordered relationship explicit.
For each category, the analysis estimates its weight `w_H`, conditional
quartet difference `Δ_H`, contribution `C_H=w_HΔ_H`, and a history-level
bootstrap interval. The contribution sum reconstructs the direct unconditional
`Δ` with maximum absolute discrepancy `2.22×10⁻¹⁶` across all cells.

The strongest conditional effects occurred in the rare 2:2 subsets, with
conditional `Q_A-Q_S` values as high as approximately 0.85 in cells with only
three or four contributing histories. These rare positive conditional effects
were too infrequent to overcome the large negative contributions from lost,
fixed, and non-2:2 histories.

## Unconditional quartet result

No unconditional anomaly was observed. All 36 cells had `Δ<0`, including the
four short-branch controls. The least negative estimate was

```text
Δ = -0.159
Ne=20, τ=0.1, origin=.75, initial frequency=.5, cross=.1
Q_S=0.445, Q_A=0.286, Q_O=0.269
```

The long-branch settings generally had very strong species-topology support
because ordinary MSC sorting dominates once `τ=2`; persistence and terminal
2:2 frequency increased in some of these settings without producing a
topology reversal. This separates “arrangement changes quartet support” from
“arrangement makes the alternative topology more probable.”

Uncertainty was estimated at the history level. The bootstrap resamples
independent histories, so loci sharing one arrangement history are not treated
as independent evolutionary replicates. The processed summary includes both
normal history-level intervals and 2,000-replicate percentile bootstrap
intervals.

## Suppression and matched controls

The primary paired comparison used cross fractions 0.1 and 1.0 with the same
`Ne`, `τ`, origin depth, and initial frequency. There are 16 such paired
settings; the short-branch controls intentionally lack an unsuppressed partner.

Strong suppression changed species-topology support modestly in this pilot.
The largest observed `P_suppressed(Q_S)-P_unsuppressed(Q_S)` was approximately
0.032, and the most negative was approximately -0.027. The largest change in
`Δ_suppressed-Δ_unsuppressed` was approximately 0.041, while the most negative
was approximately -0.068. The corresponding intervals are in
`mechanistic_pilot2_suppression_comparison.tsv`.

Cross fraction 1.0 is not exactly an ordinary MSC in this implementation. It
removes cross-arrangement suppression in the backward genealogy process, but
the arrangement history still changes state frequencies and same-state
coalescence rates. Across the paired settings, the largest deviation of
unsuppressed `Q_S` from the matched ordinary MSC expectation was approximately
0.0088, and the largest `Q_A` deviation was approximately 0.0129. Therefore
the ordinary MSC controls remain a separate neutral reference rather than an
identity imposed on the cross-fraction-1.0 mechanistic condition.

## Adaptive decision

No cell met the primary refinement criterion of an unconditional positive
estimate or a confidence interval close to a supported topology reversal.
Consequently, no 500-history rerun was silently launched. The rare 2:2
subsets did show large conditional effects, so the next targeted run should
increase independent histories in high-frequency, intermediate-persistence
settings and pre-specify a multiple-comparison rule. The reproducible plan is
in `NEXT_MECHANISTIC_PILOT2.md`.

## Files

- `datasets/processed/mechanistic_pilot2_summary.tsv`
- `datasets/processed/mechanistic_pilot2_persistence.tsv`
- `datasets/processed/mechanistic_pilot2_bias_decomposition.tsv`
- `datasets/processed/mechanistic_pilot2_suppression_comparison.tsv`
- `analysis/figures/pilot2_persistence.{png,pdf}`
- `analysis/figures/pilot2_bias_decomposition.{png,pdf}`
- `analysis/figures/pilot2_conditional_unconditional.{png,pdf}`
- `analysis/figures/pilot2_suppression_effect.{png,pdf}`
- `analysis/figures/pilot2_mechanistic_simplex.{png,pdf}`
