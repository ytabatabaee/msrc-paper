# Mechanistic smoke-test report

The smoke suite completed before the pilot for all 72 generated mechanistic
configurations. Each smoke configuration used three independent
Wright–Fisher histories and ten loci per history with
`conditioning.mode: none`. All outputs contained resolved configurations,
replicate summaries, prevalence summaries, and run metadata. The species-tree
taxa were 1, 2, 3, and 4, and terminal arrangements were sampled by the
simulator rather than imposed.

The matched ordinary MSC controls were also run through the public
`msrcsim.structured_coalescent.simulate_msc_genealogy` API. Their full pilot
validation uses 20,000 loci per cell and is reported in
`datasets/processed/msc_control_validation.tsv`.

Command:

```text
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/02_run_mechanistic_pilot.py --smoke --workers 4
```
