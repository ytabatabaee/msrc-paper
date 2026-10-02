# Stage 1B / Stage 2 Final Prospective Freeze

This analysis specification was finalized before authenticated access to real
Anopheles sequence-derived genealogy data.

- Final prospective freeze checksum:
  `64c6b70f9046ea5b3b7f493150dfd68cd6ae454beb513c3147599e4d204a4274`
- Prior prospective freeze checksum:
  `a4c366a9921ac11d3794bf29d7dbb102a278dda543a9650cd4ede1c985aa8fe2`
- Freeze date: 2026-09-16
- Git commit: `844ed61673fc009e329242203a93dfcaec151396`
- Manifest:
  `empirical/anopheles_2la/results/stage1b_synthetic/stage1b_analysis_freeze.json`

The prior checksum is retained as superseded provenance. The checksum changed
only for pre-data hardening:

- design-agnostic execution gate;
- tree-inference environment pinning;
- test-only network isolation.

This was not a response to observed biological results. No real Anopheles
sequence-derived local trees, topology tracks, SNPs, haplotypes, q1/q2/q3
expectations, quartet frequencies, or inside/outside genealogy patterns were
accessed during this hardening work.

Stage 1B eligibility is design-agnostic. A scientifically complete Stage 1A may
enable Design A, Design B, Design C, any combination, or no analysis if no
prospectively eligible designs exist. Zero Design-A quartets is not itself a
reason to alter the analysis plan. If A/B/C are all empty, Stage 1B must stop
before topology inspection with `NO_ELIGIBLE_PROSPECTIVE_DESIGNS`; this is an
inconclusive / underpowered empirical result, not evidence against MSRC.

The primary tree-inference environment is pinned as `bioconda::iqtree=2.4.0`,
executable `iqtree2`, with command template:

```bash
iqtree2 \
  -s WINDOW.fasta \
  -seed 1729 \
  -nt AUTO \
  -pre OUTPUT_PREFIX \
  -m MFP
```

No bootstrap/support threshold is used for primary quartet classification. Tree
estimation remains blind to structural predictions; frozen structural
prediction enters only after local tree files are finalized.

Future deviations from this specification must be labeled explicitly as
post-freeze deviations.
