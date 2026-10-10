# Experiment 2: mechanistic Wright–Fisher quartet pilot

This experiment uses the public `msrc-sim-replicates` workflow from the
authoritative `ytabatabaee/msrc-sim` package. Each accepted replicate contains
one independently simulated Wright–Fisher rearrangement history, terminal
arrangement sampling, and a locus set simulated conditional on that history.
The pilot is unconditional (`conditioning.mode: none`): histories in which the
rearrangement is lost, fixed, segregating, or produces another terminal
pattern are retained.

The primary species tree is `((1:tip,2:tip)A:internal,(3:tip,4:tip)B:internal)ROOT`
with true quartet `Q_S = 12|34`, arrangement-associated alternative
`Q_A = 13|24`, and remaining topology `Q_O = 14|23`. Internal branch lengths
are specified in coalescent units and converted to generations as
`internal = 2 Ne t`. The pilot varies `Ne`, internal branch length, ancestral
origin depth on `ROOT`, initial copy count, and effective cross-arrangement
recombination fraction.

The mechanistic simulator uses the realized frequency history. Its effective
backward switching and coalescence rates therefore change with arrangement
frequency and are not identified with the fixed `m01`, `m10`, `lambda0`, and
`lambda1` parameters of Experiment 1. See
`analysis/results/SIMULATOR_AUDIT.md`.

The pilot is deliberately modest. It includes ordinary MSC controls generated
through `msrcsim.structured_coalescent.simulate_msc_genealogy`, rather than a
new simulator. A positive point estimate of `Q_A-Q_S` is treated as a
candidate only; independent-replicate uncertainty and matched MSC controls
are reported before any anomaly claim.

Reproducible commands are:

```text
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/01_generate_mechanistic_pilot.py
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/02_run_mechanistic_pilot.py --smoke --workers 4
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/03_run_msc_controls.py
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/02_run_mechanistic_pilot.py --workers 4
/Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/04_process_mechanistic_pilot.py
MPLCONFIGDIR=/private/tmp/msrc-mpl /Users/ytabatabaee/opt/anaconda3/bin/python simulations/quartet/mechanistic/analysis/scripts/05_plot_mechanistic.py
```

The completed pilot contains 72 configurations, 12 independent histories per
configuration, and 20 loci per history. Matched MSC controls contain 20,000
independent loci. Figures A–D are exported separately as
`mechanistic_scenario`, `mechanistic_unconditional_probabilities`,
`mechanistic_quartet_bias_map`, and `mechanistic_conditional_vs_unconditional`;
`mechanistic_quartet_simplex` is retained as a supplementary diagnostic.
