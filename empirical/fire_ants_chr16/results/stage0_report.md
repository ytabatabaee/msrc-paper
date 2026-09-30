# Fire ants chromosome 16 Stage 0 report

Stage 0 scaffolds the fire-ant chromosome-16 social-supergene analysis and freezes the design before local topology support is reproduced in this repository.

## Upstream repository

- Repository: `https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression`
- Default branch: `master`
- HEAD commit: `bdc7823941680fa560e61b6467af9f77950f64b0`
- Retrieval date: `2026-09-30`

## Region inventory

| region | chromosome | windows | observed BUSCO-window span |
| --- | --- | ---: | --- |
| `chr1` | `chr1` | 117 | 289334-29185097 |
| `chr16A` | `chr16` | 42 | 23916-7386513 |
| `chr16B` | `chr16` | 2 | 28498399-28943581 |
| `chr16_supergene` | `chr16` | 52 | 11680438-27917498 |

Coordinates are author-designated analysis regions and observed BUSCO-window spans. They are not treated as nucleotide-resolved inversion breakpoints.

## Frozen focal quartet

- A: `invicta/macdonaghi_SB` (`inv_mac_SB`)
- B: `invicta/macdonaghi_Sb` (`inv_mac_Sb`)
- C: `richteri_SB` (`richteri_SB`)
- D: `richteri_Sb` (`richteri_Sb`)

Frozen resolutions:

- `species_split`: `AB|CD`
- `haplotype_split`: `AC|BD`
- `third_split`: `AD|BC`

Stage-0 focal-quartet checksum: `5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6`

## Raw source-file checksums

- `data/fire_ants_chr16/raw/upstream/README.md`: `d09c83ab0e53f601ca3147e63847cdd115df6693cdd259da37d98fc56a95a122`
- `data/fire_ants_chr16/raw/upstream/Topology weighting/README.md`: `ad5ffe7130c26cdca45f868639c7b89049264e299ad906a72968231dfd24e078`
- `data/fire_ants_chr16/raw/upstream/Topology weighting/supergene/README.md`: `dcd04151ecc09d61b840192b711c45789e0d78434a67a327a04b812df847a4bf`
- `data/fire_ants_chr16/raw/upstream/Topology weighting/results/window_coordinates`: `c46139b1a3a33c8e22f6b1ee85cbe08d2d7dfd23b47f0263d1fab579e157d791`
- `data/fire_ants_chr16/raw/upstream/Coalescence-based phylogenies/README.md`: `7b5cf496f91d617615550193c94111876da8411f916faa0571269afd7c80b90a`

## Guardrails

The Stage-0 validator rejects topology-bearing analysis inputs such as TWISST weights, topology-output trees, input trees, Newick files, quartet-support files, QQS files, and precomputed q1/q2/q3 support tables.

No TWISST weights, local trees, Newick files, ASTRAL/ASTER outputs, or quartet-support values were parsed in Stage 0.
