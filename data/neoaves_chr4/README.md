Neoaves chromosome-4 analysis data.

- `raw/`: immutable external source data copied from `genetreesupport/` and
  `rearrangements/` for the chr4 analyses.
- `metadata/`: source metadata, manifests, checksums, and provenance records.
- `intermediate/`: generated stage handoff products grouped by stage, including
  locus tables, structural states, denovo Stage-2.7c breakpoints, Stage-3A
  enrichment outputs, and exploratory Stage-3B-v1 outputs.
- `processed/`: processed analysis tables and manifests used by downstream or
  reporting stages. `stage03b_v2_event_state_test/` contains the focused final
  Stage-3B-v2 event-state deliverables. The top level also includes a
  compatibility mirror of the TSV/JSON products from
  `empirical/neoaves_chr4/results/`.

Scripts still write active outputs to `empirical/neoaves_chr4/results/` and
figures to `empirical/neoaves_chr4/figures/`. The data tree is the externalized
copy of inputs and generated data products.
