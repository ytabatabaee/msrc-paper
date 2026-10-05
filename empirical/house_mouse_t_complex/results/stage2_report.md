# House-mouse t-complex Stage 2 report

## Scope and gate

This analysis uses the published Kelemen & Vicoso 5-kb maximum-likelihood
gene trees in the chr17 t-complex window range, approximately
chr17:5,000,000–37,054,999. The primary archive is the coverage plus
allele-ratio filtered RAW SNP dataset. The 4,046 windows are densely linked
local genealogies and are not genome-wide background or an outside-inversion
control. `STANDARD_ONLY` is a within-region control that excludes the
reconstructed/pseudo-t haplotype sequences.

The fixed quartet baseline was specified independently as
`Q_SPECIES = (musculus, castaneus) | (domesticus, spretus)`, matching the
commonly supported rooted house-mouse species history with M. spretus as
outgroup. The competing focal split returned by the t-inclusive exploratory
ASTRAL run is `Q_T_ALT = (domesticus, musculus) | (castaneus, spretus)`.

## Fixed quartet scan

The exact scan evaluates every valid individual combination in each window,
using topology only and treating induced polytomies as unresolved. Aggregate
fractions are:

| treatment | q_species | q_t_alt | q_other | unresolved | mean delta species-alt |
|---|---:|---:|---:|---:|---:|
| STANDARD_ONLY | 0.369803 | 0.343630 | 0.286365 | 0.000202 | 0.026174 |
| T_ONLY | 0.455309 | 0.267625 | 0.276911 | 0.000154 | 0.187684 |
| ALL_TIPS | 0.323187 | 0.365290 | 0.311230 | 0.000293 | -0.042103 |

Thus the all-tip distribution shifts away from the independently specified
species-history quartet toward Q_T_ALT relative to STANDARD_ONLY. The T_ONLY
distribution itself favors Q_SPECIES, so the shift is an interaction between
tip composition and multi-individual sampling across the linked windows; it
should not be described as a simple claim that every pseudo-t local tree has
the alternative topology.

## Balanced sampling

The deterministic sampler selects three castaneus, domesticus, musculus, and
spretus tips per replicate. The local exploratory run completed seven
replicates per treatment before runtime became prohibitive. Those partial
frequencies are:

| treatment | n | Q_SPECIES | Q_T_ALT | Q_OTHER |
|---|---:|---:|---:|---:|
| B0_STANDARD_MATCHED | 7 | 0.7143 | 0.2857 | 0.0000 |
| B1_T_MATCHED | 7 | 1.0000 | 0.0000 | 0.0000 |
| B2_MIXED_BALANCED | 7 | 0.2857 | 0.7143 | 0.0000 |

The matched mixed treatment therefore retains the topology sensitivity in the
available pilot, while the matched T-only treatment does not. These are
diagnostic partial results, not a 200-replicate inferential summary. The
checkpointable command for the requested full run is:

```bash
python empirical/house_mouse_t_complex/scripts/02b_balanced_sampling_astral.py \
  --n-replicates 200 --seed 20261005 --threads 2
```

Existing replicate outputs are reused by the script.

## Filtering robustness

The PASS SNP archive contains 4,047 usable trees and reproduces the same
qualitative sensitivity: STANDARD_ONLY returns Q_SPECIES and ALL_TIPS returns
Q_T_ALT. The primary coverage plus allele-ratio filtered RAW SNP archive has
4,046 usable trees and the same pair of focal topologies. The coverage-filtered
RAW SNP archive has 4,046 coordinate-recoverable trees, but its tip labels are
incompatible with the frozen Stage-1 mapping, so no ASTRAL or quartet result
was produced for it. Trees from filtering datasets were never mixed.

## Population analysis

The seven-population STANDARD_ONLY tree and ALL tree differ by RF distance 2
and two changed canonical splits. The 50-kb and 100-kb all-tip thinned trees
retain that RF distance 2 relative to STANDARD_ONLY. This is a sensitivity
result for recent population history and gene flow; it does not identify one
population topology as biologically correct.

## Linkage and pseudo-replication

The all-tip focal topology is Q_T_ALT for the dense 5-kb input and for every
tested phase at 10, 20, and 50 kb. At 100 kb, 9/10 tested phases return
Q_T_ALT and 1/10 returns Q_SPECIES. At 250 kb, 8/10 return Q_T_ALT and 2/10
return Q_OTHER; at 500 kb, 7/10 return Q_T_ALT and 3/10 return Q_OTHER.
LocalPP decreases from 1.0 at 5 kb to phase-dependent values around 0.40–0.96
at 50–100 kb and around 0.43–0.66 at 250 kb, while CU lengths also vary
substantially. The local run used the first ten deterministic 5-kb phase
offsets for spacings with more than ten feasible phases; a full phase sweep is
the next computational extension.

The archive contains paths referring to concatenated “Non-Recombined Regions”
trees, but no machine-readable genomic interval definitions were recovered.
Therefore `source_defined_blocks_available = false`, and no blocks were
defined from the observed quartet signal.

## Interpretation and limitations

The fixed-quartet scan supports topology sensitivity when all mapped tips are
included, and the matched mixed pilot suggests that the signal is not explained
by unequal counts alone. The reduction in ASTRAL support under spatial thinning
is consistent with linked-window support inflation, so thousands of local
windows must not be treated as thousands of independent biological loci. The
filtering result is concordant for PASS and primary datasets but incomplete for
the incompatible third dataset. Full 200-replicate balanced sampling, a full
phase sweep, and a source-defined block analysis remain outstanding; no causal
claim or proof of MSRC is made here.

The newer 2025 inversion breakpoints were not overlaid because assembly,
reference coordinate system, and lift-over compatibility with the 2017 window
coordinates were not established.

## Reproducibility

The source archive is recorded by basename, SHA256
`f925084a87ac3837661ad99f931da608c129ea3e590c3470030b3ba649d47877`, and
byte size 143,697,895. Generated paths are repository-relative. The exact
commands and output checksums are recorded in `stage2_manifest.json`; the
fixed scan, balanced run, population run, filtering run, linkage run, and
figure validation tables are retained beside this report.
