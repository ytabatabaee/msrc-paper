# Atlantic cod data

Data tree for the Atlantic cod inversion-supergene empirical analysis.

- `raw/`: immutable external source data; large files should remain untracked unless explicitly approved.
- `metadata/`: source manifests, region definitions, population-state tables, provenance, and checksums.
- `intermediate/`: generated stage handoff products.
- `processed/`: canonical analysis-ready tables.

Large VCFs, posterior tree archives, BAMs, FASTQs, and similar raw genomic files should remain external or ignored unless explicitly needed for a later approved stage.

Stage 0 records only source and region manifests. It does not download source data, parse local trees, assign population arrangements, or generate processed analysis tables.
