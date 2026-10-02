# Stage 1B / Stage 2 Prospective Analysis Plan

`SYNTHETIC_ONLY = TRUE` for all current validation outputs. This plan was
written before inspecting any real Anopheles sequence-derived genealogy,
topology, local-tree, SNP, haplotype, q1/q2/q3, QQS/BQS, NJ-tree, or
inside/outside topology data.

## Fixed Regions

Coordinates are AgamP4 chromosome arm `2L`, 1-based inclusive:

| Region | Start | End |
| --- | ---: | ---: |
| left_flank | 15,524,058 | 20,524,057 |
| inversion_2La | 20,524,058 | 42,165,532 |
| right_flank | 42,165,533 | 47,165,532 |

The primary flank width is fixed at 5 Mb. Sensitivity flank widths are fixed at
2 Mb and 10 Mb. Coordinates will be clipped only if chromosome boundaries require
it. Flank widths will not be selected from observed topology signal.

## Windowing

The primary local-window scheme is fixed physical 100 kb non-overlapping
windows. Sensitivity windows are 50 kb and 200 kb. Each window carries
chromosome, start, end, midpoint, region class, distance to the left breakpoint,
and distance to the right breakpoint.

## Frozen Prediction Versus Observed Genealogy

Stage 1B eligibility is design-agnostic. A scientifically complete Stage 1A may
enable Design A, Design B, Design C, any combination, or no analysis if no
prospectively eligible designs exist. Zero Design-A quartets is not itself a
reason to alter the analysis plan.

The real-data gate requires
`data/anopheles_2la/processed/stage1a_freeze_complete.json`, successful
Stage-1A freeze/provenance validation, and at least one nonempty prospectively
frozen eligible design table or repository-equivalent frozen row set. Design A
is represented by `frozen_strict_quartets_stage1a.tsv` or equivalent strict
2:2 rows. Design B is represented by frozen arrangement-replacement contrasts.
Design C is represented by frozen geography/population-controlled contrasts.

The Stage 1B runner records `design_A_enabled`, `design_B_enabled`, and
`design_C_enabled`. It executes only analyses supported by nonempty frozen
designs: Design A enables strict 2:2 inside/outside analysis, spatial
localization, and block-aware analysis; Design B enables arrangement-replacement
paired contrasts and their appropriate spatial summaries; Design C enables
geography/population-controlled summaries. No missing design is fabricated.

If Stage 1A completes scientifically but Design A, Design B, and Design C are
all empty, Stage 1B stops before topology inspection with
`NO_ELIGIBLE_PROSPECTIVE_DESIGNS`:

```text
Stage 1A completed successfully, but no prospectively eligible
structural-state analysis designs were available.
```

This is an inconclusive / underpowered empirical result, not evidence against
MSRC.

Frozen Stage 1A design rows must contain quartet or contrast identifiers, sample
labels as applicable, species labels, arrangement states, predicted arrangement
split or contrast direction, design class, independence group, and provenance
fields. Stage 1B consumes these rows without recomputing predicted splits from
local trees.

Observed local-tree topology is a Stage 2 output. The core prospective test is
not "Are gene trees unusual inside 2La?" It is "Does a topology predicted
independently from structural arrangement state become specifically enriched
inside 2La?"

## Missing Tips And Polytomies

Every tree/window/quartet record receives `usable = true/false` and an exclusion
reason. Missing tips, duplicate tips, malformed Newick, and unresolved quartets
are recorded. Unresolved quartets are excluded from topology-frequency
denominators but their counts are reported. Polytomies are not randomly
resolved.

## Primary Statistic And Test

For each frozen structural prediction:

`f_in = P(observed topology = predicted arrangement split | inside 2La)`

`f_out = P(observed topology = predicted arrangement split | flanks)`

The primary effect is `Delta_arr = f_in - f_out`, with enrichment `f_in / f_out`
when `f_out > 0`. The primary test is a block-aware circular-shift permutation
with 10,000 permutations in final analysis. A naive per-window permutation is
implemented only as a diagnostic because it can be anti-conservative when linked
neighboring windows are treated as independent.

## Pseudo-Replication

Two summaries are reported. Window-weighted counts give each usable genomic
window one observation. Block-normalized counts collapse consecutive windows
with the same quartet topology state into one run and give each run total weight
one.

## Breakpoint Localization

Topology transitions are located along 2L and compared with the fixed 2La
breakpoints at 20,524,058 and 42,165,532. The breakpoints are not moved after
observing tree transitions.

## Tree Inference Policy

The primary tree-estimation wrapper is IQ-TREE 2, pinned as
`bioconda::iqtree=2.4.0` in
`empirical/anopheles_2la/config/tree_inference_environment.yml`. The executable
is `iqtree2`, verified before the first real run with `iqtree2 --version`.

The frozen primary command template is:

```bash
iqtree2 \
  -s WINDOW.fasta \
  -seed 1729 \
  -nt AUTO \
  -pre OUTPUT_PREFIX \
  -m MFP
```

Model selection is IQ-TREE ModelFinder Plus (`-m MFP`) for every window. The
seed is `1729`; the thread policy is `-nt AUTO`. No bootstrap, SH-aLRT, UFBoot,
or support threshold is used in the primary quartet-classification analysis,
and no support-based unresolved threshold is introduced. Expected output files
are `OUTPUT_PREFIX.treefile`, `OUTPUT_PREFIX.iqtree`, `OUTPUT_PREFIX.log`,
`OUTPUT_PREFIX.ckp.gz`, and `OUTPUT_PREFIX.model.gz` when emitted by IQ-TREE.

Tree inference records input window, command, software version, model policy,
seed, thread policy, exit status, and resulting tree. Quartet predictions,
Design-A/B/C topology expectations, and q1/q2/q3 expectations are not inputs to
tree inference; frozen structural prediction enters only after local tree files
are finalized.
