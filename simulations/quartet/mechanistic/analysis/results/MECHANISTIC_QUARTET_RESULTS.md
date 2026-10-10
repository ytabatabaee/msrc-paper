# Mechanistic Wright–Fisher quartet pilot

This pilot used `ytabatabaee/msrc-sim` commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9` (package version 0.8.8) through
the public `msrc-sim-replicates` workflow. No simulator source was copied or
modified. The audit is in `SIMULATOR_AUDIT.md` and the provenance record is in
`datasets/simulator_provenance.json`.

The simulator replicate test `tests/test_replicates.py` passed at this commit.
The paper-side smoke and integrity tests also passed.

The species tree was balanced, `((1,2)A,(3,4)B)ROOT`, with tip branches set
to `2Ne` generations and each daughter ancestral branch set to `2Neτ`
generations. Because the simulator allows the 1–2 and 3–4 pairs to coalesce
on their two separate daughter ancestral branches, the comparable unrooted
MSC quartet length is `t = 2τ`. The pilot varied `Ne ∈ {20,100}`,
`τ ∈ {0.1,0.5,2}`, ROOT origin depth `{0.25,0.75}`, initial copy fraction
`{0.025,0.20,0.50}`, and effective cross-arrangement recombination fraction
`{1,0.1}`. Every mechanistic configuration was unconditional:
`conditioning.mode: none`. Histories that lost, fixed, or retained the
rearrangement were retained in the analysis.

The pilot had 72 cells, 12 independent Wright–Fisher histories per cell, and
20 loci per history. The small smoke run completed before the pilot. The
mechanistic replicate API reports history and quartet summaries but does not
emit an affected-locus count, so `affected_locus_fraction` is explicitly
recorded as unavailable in the processed table rather than inferred from
topology counts. Terminal pattern frequencies and persistence probabilities
are reported separately.

## Validation

Six ordinary MSC controls were simulated with 20,000 independent loci per
cell through `msrcsim.structured_coalescent.simulate_msc_genealogy`. Their
observed quartet frequencies agree with

```text
P(Q_S) = 1 - 2/3 exp(-t)
P(Q_A) = P(Q_O) = 1/3 exp(-t),  t = 2τ.
```

The largest absolute control z score was 1.20 and the largest absolute
probability error was 0.00304. This validates the balanced-tree control and
the branch-length conversion used for the pilot.

## Pilot result

No mechanistic pilot cell had a positive point estimate of
`Δ = P(Q_A) − P(Q_S)`: all 72 estimates were negative, and no cell had a
95% independent-replicate interval whose lower bound was positive. The least
negative point estimate was `Δ = -0.0583`; it occurred at
`Ne=100, τ=0.1, origin=.25, initial copy fraction=.025, cross fraction=.1`
and its 95% interval was `[-0.1870, 0.0704]`. This is a small, noisy pilot and
does not establish that an unconditional reversal is impossible.

The terminal `1010` pattern appeared in 6 of 72 cells, with maximum frequency
0.0833. Thus the very strong conditional 1010 behavior from Experiment 1 was
rare in this initial unconditional design. The mechanistic pilot therefore
does not support an unconditional anomalous quartet distribution. The
negative result is consistent with the distinction between a strong
arrangement-conditioned genealogy bias and its contribution after averaging
over stochastic origin, loss, fixation, switching, and terminal sampling.

The pilot also shows why the Experiment 1 switching parameters cannot be
read directly as mechanistic values. In the mechanistic simulator, backward
switching and same-state coalescence rates depend on the realized frequency
history and local effective population size. Experiment 1 instead fixes
`m01`, `m10`, `lambda0`, and `lambda1` over a structured interval. The two
experiments are mechanistically related but are not parameter-equivalent.

The saved smoke, pilot, and control runs took approximately 40.6 minutes in
total on the local machine.

## Outputs and interpretation

The processed unconditional summaries are in
`datasets/processed/mechanistic_quartet_pilot.tsv` and terminal-pattern
contributions are in `datasets/processed/mechanistic_pattern_contributions.tsv`.
The MSC comparison is in `datasets/processed/msc_control_validation.tsv`.
Figures A–D are separate diagnostic figures: `mechanistic_scenario`,
`mechanistic_unconditional_probabilities`, `mechanistic_quartet_simplex`, and
`mechanistic_quartet_bias_map`, each as PNG and PDF.

The scenario figure uses a saved ROOT-origin pilot history as the illustrative
history and labels the terminal arrangement as sampled rather than imposed.
The probability and bias figures use independent-replicate uncertainty. The
simplex summarizes the pilot distributions; it is not the exact conditional
simplex from Experiment 1.

## Follow-up

The next search should increase the number of independent histories before
increasing loci within a history. It should concentrate on biologically
plausible settings with intermediate or high initial copy fraction, longer
ROOT persistence, and terminal 2:2 patterns, while retaining matched MSC
controls. Any positive `Δ` should be refined with substantially more
independent histories and a pre-specified interval-based criterion. A
conditional 1010 analysis may diagnose the mechanism, but it cannot replace
the unconditional analysis.
