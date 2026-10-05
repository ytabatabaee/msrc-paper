# House mouse t-complex Stage 1 report

Status: prepared

Usable ML trees written: 4046

Coordinates recovered for all valid trees: True

ASTRAL mapping files use only exact/strong source-backed mappings from `metadata/tip_mapping.tsv`.
Tentative and unresolved mappings are excluded automatically.

Limitation: these windows represent the published chr17 t-complex interval, not a genome-wide background.

Exploratory ASTRAL4 diagnostics were run because a local executable was available at
`/Users/ytabatabaee/Desktop/ASTER/bin/astral4`. These runs are Stage-1
diagnostics only and should not be treated as final biological interpretation.

Treatments prepared and run:

- `T0_STANDARD`: 4,046 trees, pseudo-t/t-haplotype tips pruned, 40 standard/outgroup tips.
- `T1_ALL_WINDOWS`: 4,046 trees, all 55 mapped tips.
- `T2_SPATIALLY_THINNED_50kb`: 467 trees, all 55 mapped tips.
- `T2_SPATIALLY_THINNED_100kb`: 248 trees, all 55 mapped tips.

The ASTRAL4 map files use underscore-safe taxon labels as required by the
ASTER/ASTRAL4 whitespace-delimited mapping interface. Results are summarized in
`results/stage1_aster/stage1_astral_summary.tsv`.
