# Atlantic cod empirical analysis

This directory scaffolds a future multi-population Atlantic cod analysis focused on the four inversion-associated supergenes reported by Matschiner et al. (2022):

- `LG01`
- `LG02`
- `LG07`
- `LG12`

The biological question is whether local population-tree topologies inside these inversion regions preferentially group differentiated Atlantic cod populations according to chromosomal arrangement, relative to local trees outside the inversions.

Stage 0 was only a repository scaffold and planning/provenance freeze.

Stage 1 is complete. It downloaded small published processed source-data files and selected upstream provenance scripts, then normalized the authors' published 250-kb population-window trees into `../../data/atlantic_cod/processed/cod_window_trees.tsv`. Stage 1 did not run SNAPP or BEAST, download raw reads or whole-genome VCFs, assign arrangement states, enumerate quartets, test arrangement concordance, or draw biological conclusions.

Stage 2 is complete. It freezes population identities, ecotype labels, and ancestral/derived arrangement states in `../../data/atlantic_cod/processed/population_arrangements.tsv`, and freezes the structural-only quartet candidate table in `../../data/atlantic_cod/processed/stage2_structural_quartet_candidates.tsv`.

Frozen Stage-2 structural candidate checksum:

```text
303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c
```

Stage 3 is complete. It extracts the independently inferred outside-supergene population tree from Nature Source Data Fig. 2, writes the canonical baseline tree to `../../data/atlantic_cod/processed/baseline_population_tree.nwk`, and freezes baseline population-history splits for every Stage-2 quartet in `../../data/atlantic_cod/processed/stage3_baseline_quartets.tsv`.

Stage 3 now gives each quartet two independent predictions:

```text
structural arrangement prediction
+
baseline population-history prediction
```

Stage 4A is complete. It is the first unblinded local-tree stage: it compares the published 250-kb window trees with the frozen Stage-2 arrangement prediction and frozen Stage-3 baseline prediction, writes descriptive quartet/window summaries, and generates topology-track figures. Stage 4A is descriptive only; no spatial-null test, permutation test, P-value, formal MSRC model fit, or SNAPP/BEAST rerun has been performed.

Stage 4B is complete. It tests whether the Stage-4A topology-alignment signal is unusually aligned with the independently frozen inversion intervals using exact circular shifts of the local-topology track, with boundary, block-placement, leave-one-population-out, topology-run, and global secondary summaries. Stage 4B does not fit MSRC parameters, infer inversion ages, rerun SNAPP/BEAST, or establish causal inversion effects.

Stage 5 is complete. It adds the final physical-coordinate-aware sensitivity check, writes manuscript-ready figures/tables/text, and freezes the Atlantic cod empirical analysis for manuscript use in `results/FINAL_ANALYSIS.md` and `results/stage5_final_manifest.json`.

Stage 6 is complete as a post-freeze extension. It evaluates ASTRAL4 summary-tree topology sensitivity to treating adjacent linked inversion windows as separate input trees. Stage 6 does not replace or modify the Stage-5 freeze.

Stage 7 is complete as a post-freeze extension. It uses the native time-scaled published SNAPP MCC window trees to test divergence-time sensitivity for the same population pairs inside versus outside the frozen inversion intervals. Stage 7 writes pairwise MRCA-time tables, inside/outside time-shift summaries, topology/time comparison tables, and figures in `results/stage7_*` and `figures/atlantic_cod_stage7_*`. One boundary-overlap LG07 window failed the nonnegative branch-length validation and was excluded from Stage-7 MRCA calculations; the primary fully inside/outside counts were unaffected. Per-window posterior SNAPP trees were not available in the local/public inventory, so Stage 7 reports MCC point-estimate sensitivity only and does not rerun SNAPP.

These analyses still do not prove MSRC. The final Atlantic cod freeze provides spatial-null support for the predeclared topology-alignment signal on some linkage groups.

## Conceptual role

The existing empirical systems play different roles:

```text
Neoaves:
structural transitions <-> genealogy-regime transitions across a multi-species radiation

Anopheles:
arrangement-state partition -> predicted local topology in a species-complex/inversion system

Atlantic cod:
arrangement-state partition -> predicted local population-tree topology across differentiated Atlantic cod populations
```

Atlantic cod is not a direct test of persistence through speciation. It tests the more basic prediction that recombination-partitioned chromosomal backgrounds generate persistent local genealogy structure.

The future theorem-level prediction is:

```text
arrangement-state partition -> elevated support for the corresponding local phylogenetic split
```

inside the inversion relative to collinear regions.

## Primary future test

The key future comparison should be:

```text
P(arrangement-concordant topology | inside inversion)

vs.

P(arrangement-concordant topology | outside inversion)
```

This comparison must not be calculated in Stage 0.

## Anti-circularity design

The Atlantic cod workflow must freeze the following inputs before testing local topology enrichment:

1. inversion coordinates;
2. population arrangement states;
3. baseline population relationships.

The arrangement prediction must not be derived from local window trees themselves. Future arrangement-state assignments should come from independent metadata or documented structural-genotype resources, not from topology results. The baseline population tree should be defined from collinear or non-supergene regions before testing enrichment inside the inversions.

## Data sources

Main source paper:

```text
Matschiner et al. 2022
Supergene origin and maintenance in Atlantic cod
Nature Ecology & Evolution
DOI: 10.1038/s41559-022-01661-x
```

Associated public resources:

```text
GitHub:
https://github.com/mmatschiner/supergenes

Zenodo:
https://doi.org/10.5281/zenodo.4560275
```

The Nature article also provides small source-data files associated with figures, including windowed population-tree results. These are recorded in `../../data/atlantic_cod/metadata/source_manifest.tsv` but have not been downloaded.

## Coordinate convention

The Stage-0 region manifest preserves the authors' original linkage-group coordinate convention as author-reported base-pair positions. No BED-like conversion has been applied. Later scripts that require 0-based half-open intervals must perform and document that conversion explicitly without changing the frozen manifest.

Known inversion/supergene regions currently recorded:

- `LG01: 9,114,741-26,192,386`
- `LG02: 18,489,307-24,048,607`
- `LG07: 13,610,591-23,019,113`
- `LG12: 638,100-14,327,837`

These four intervals are frozen from `mmatschiner/supergenes` `cod_phylogenomics/src/split_vcf.sh` for the cod phylogenomics/window-tree analysis. Other source-repository sub-workflows, especially demography, may use slightly different inversion limits. Those alternative coordinate sets must not be mixed silently with the `cod_phylogenomics` coordinate set.

## Planned future tables

The first canonical future table will be:

```text
data/atlantic_cod/processed/cod_window_trees.tsv
```

Intended columns:

```text
lg
window_id
start
end
midpoint
n_sites
tree_newick
n_taxa
inversion_id
inside_inversion
distance_to_left_boundary
distance_to_right_boundary
```

The second canonical future table will be:

```text
data/atlantic_cod/processed/population_arrangements.tsv
```

Intended columns:

```text
population
LG01_state
LG02_state
LG07_state
LG12_state
source
notes
```

Future state labels should be:

```text
ancestral
derived
unknown
```

No population arrangement states are assigned in Stage 0.

The third canonical future table will be:

```text
data/atlantic_cod/processed/local_quartets.tsv
```

Intended columns:

```text
lg
window_id
start
end
quartet_id
taxon1
taxon2
taxon3
taxon4
arrangement_pattern
baseline_topology
arrangement_topology
local_topology
inside_inversion
informative
classification
```

No synthetic rows are created in Stage 0.

## Future stages

Stage 0
Data/source inventory and freeze structural definitions.

Stage 1
Retrieve published 250-kb population-tree outputs and reconstruct the chromosome-wide tree tracks. Completed in `results/stage1_report.md`.

Stage 2
Freeze population arrangement states independently of local topology. Completed in `results/stage2_report.md`.

Stage 3
Define baseline population tree from collinear/non-supergene regions. Completed in `results/stage3_report.md`.

Stage 4A
Extract local quartets from published 250-kb trees and report descriptive inside/outside topology alignment. Completed in `results/stage4a_report.md`.

Stage 4B
Test the descriptive Stage-4A pattern against spatially valid nulls. Completed in `results/stage4b_report.md`.

Stage 5
Test arrangement-concordant topology enrichment inside versus outside inversions.
Completed as final robustness, analysis freeze, and manuscript-ready output generation in `results/FINAL_ANALYSIS.md`.

Stage 6
Test whether topology transitions align with inversion boundaries using spatially valid nulls.

Stage 7
Inspect special local exchange events, including the reported Bornholm LG12 segment.

Stage 8
If useful, extract posterior quartet probabilities from existing SNAPP posterior tree samples.

Raw-read processing and rerunning SNAPP/BEAST are not first-line tasks. The intended order is:

```text
published processed trees
-> topology-track reconstruction
-> quartet analysis
-> only then deeper raw-data work if needed
```

## Layout

- `config/`: human-readable configuration.
- `scripts/`: future stage scripts; currently documentation only.
- `results/`: Stage-0 planning/provenance notes and future outputs.
- `figures/`: reserved for future figures.
- `tests/`: future test documentation.
- `../../data/atlantic_cod/`: data tree with `raw/`, `metadata/`, `intermediate/`, and `processed/`.
