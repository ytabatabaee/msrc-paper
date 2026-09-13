# Neoaves chromosome 4 local QQS tracks

This directory builds an isolated empirical table and first raw spatial plot for
the three focal Columbea-related hypotheses on chromosome 4. It does not modify
the original `genetreesupport/` analysis files.

## Layout

- `scripts/`: executable stage scripts.
- `results/`: tabular analysis outputs and final manifests.
- `figures/`: rendered PDF/PNG figures.
- `config/`: run profiles or parameter files when defaults are externalized.
- `tests/`: shared regression tests; several stage-local tests remain embedded
  behind each script's `--run-tests` flag.
- `../../data/neoaves_chr4/`: data tree following the `raw/`, `metadata/`,
  `intermediate/`, `processed/` convention used by other empirical analyses.
  Raw source inputs are copied under `raw/`; generated stage handoff products
  are grouped under `intermediate/`; manifests/provenance are copied under
  `metadata/`; final and compatibility TSV/JSON products are under
  `processed/`.

## Inputs

- `genetreesupport/63K_trees.names_header.txt.xz`: locus metadata for 63,430
  loci. The published R script uses `ws` and `we`; here these are copied to
  `start` and `end`, and `midpoint = (ws + we) / 2`.
- `genetreesupport/clade-rec.stat.xz`: QQS source statistics with columns
  `C1 C2 topology_id Gene x`.
- `genetreesupport/draw-movingaverage.r`: source for the published clade
  mapping, chromosome parsing, and QQS formulas reproduced here.
- `genetreesupport/README.md`: published description of the quartet statistics.

Chromosome labels can be parsed as in `draw-movingaverage.r`, but the spatial
chromosome-4 table intentionally keeps only metadata rows whose exact
`Chromosome` value is `chr4`. Labels such as `chr4_AADN03010064_random` and
`chr4_JH375171_random` have local scaffold coordinates and are excluded unless
they can be mapped onto main-chromosome coordinates.

## Clades and topology labels

The focal clade mapping is reproduced from the QQS section of
`draw-movingaverage.r`:

| clade | C1 | C2 |
| --- | --- | --- |
| Columbea | Columbimorphae | Phoenicopteriformes |
| N61 | Otidimorphae | Columbimorphae |
| N62 | Columbiformes | OtherColumbimorphae |

For each row, `S` is the sister of `C1` and `C2`, and `O` is outside the parent
of the tested clade. The three local QQS tracks are:

- `q1 = (C1,C2)|(S,O)`
- `q2 = (C1,S)|(C2,O)`
- `q3 = (C2,S)|(C1,O)`

## QQS calculation

`clade-rec.stat.xz` has topology IDs 1, 2, 3 for the three resolved
quadripartitions and 4 for the unresolved quadripartition. After pivoting to
`x1`, `x2`, `x3`, and `x4`, the script computes:

```text
d1 = x1 - x4
d2 = x2 - x4
d3 = x3 - x4
```

The published formulas are then applied exactly:

```text
q1 = (d2 + d3 - d1) / (d1 + d2 + d3)
q2 = (d1 + d3 - d2) / (d1 + d2 + d3)
q3 = (d1 + d2 - d3) / (d1 + d2 + d3)
```

Rows with a zero denominator have undefined QQS and are dropped with a validation
count.

## Regeneration

From the repository root:

```bash
python3 empirical/neoaves_chr4/scripts/01_build_locus_table.py --run-tests
```

The default `python3` on this machine has a broken pandas/NumPy ABI, so the
script uses standard-library streaming for the table build and matplotlib for
the figure. It writes:

- `empirical/neoaves_chr4/results/locus_table_chr4.tsv`
- `empirical/neoaves_chr4/results/validation_summary.tsv`
- `empirical/neoaves_chr4/results/chr4_clade_means.tsv`
- `empirical/neoaves_chr4/figures/chr4_raw_quartet_tracks.pdf`
- `empirical/neoaves_chr4/figures/chr4_raw_quartet_tracks.png`

## Stage 2: structural states

`02_build_structural_states.py` builds a structural-only representation for the
later MSRC test. It does not read `q1`, `q2`, `q3`, dominant topologies, or
Stage-1 change points when defining intervals or breakpoints.

The main block source is `rearrangements/totab_chr4.csv.xz`, described in
`rearrangements/README.md` as postprocessed maf2synteny output containing all
syntenic blocks between GalGal6 and focal species for chromosome 4. The raw
maf2synteny workflow under `maf2synteny/` used stork
`GCA_017639555.1_bCicMag1.pri` sequence `CM030196.1` as the structural-analysis
reference, while this Stage-2 table uses the paired
`GCF_000002315.5_GalGal6.chr4` block rows as the main physical coordinate axis
needed for merging with the chr4 genealogy table.

The canonical block table pairs each maf2synteny `Block` for a target genome
with its exact GalGal6 `chr4` row. Reference coordinates are sorted into
`reference_start` and `reference_end`; the original start/end fields and the
maf2synteny `Strand` values are also retained. The script requires the reference
chromosome to be exactly `chr4`, so `chr4_*_random` scaffolds are not treated as
main-chromosome coordinates.

Structural intervals are inferred from the detailed published maf2synteny
structural block calls in `rearrangements/outlier-regions-by-maf2synteny.tsv`,
not from quartet support. Blocks are first separated into chicken-coordinate
clusters, then split by query chromosome and target-coordinate continuity using
a 1 Mb gap threshold. This preserves multiple adjacent or overlapping pieces
instead of forcing them into a single rearrangement. Event types are conservative:
specific labels such as `inversion` or `reordered_block` are assigned only from
block orientation/order behavior; otherwise intervals remain `unknown` or
`complex`.

Breakpoints are the left and right boundaries of these inferred structural
intervals. They are recorded as published/inferred structural boundaries and are
not nucleotide-resolved breakpoints.

Arrangement states are deterministic `A0`, `A1`, ... labels assigned from
normalized query-chromosome, query-order, and orientation signatures within each
structural interval. The signatures preserve enough block order/orientation
information for later quartet-level arrangement sharing tests, rather than
reducing the data to an inside/outside flag.

The interval audit compares the parsed intervals with the rough published
summary in `rearrangements/outlier-regions-by-maf2synteny-summary.tsv`.
Discrepancies are reported rather than corrected to match the summary.

Regenerate Stage 2 from the repository root with:

```bash
python3 empirical/neoaves_chr4/scripts/02_build_structural_states.py --run-tests
```

It writes:

- `empirical/neoaves_chr4/results/structural_blocks.tsv`
- `empirical/neoaves_chr4/results/structural_intervals.tsv`
- `empirical/neoaves_chr4/results/structural_breakpoints.tsv`
- `empirical/neoaves_chr4/results/arrangement_states.tsv`
- `empirical/neoaves_chr4/results/arrangement_state_definitions.tsv`
- `empirical/neoaves_chr4/results/structural_interval_validation.tsv`
- `empirical/neoaves_chr4/results/structural_validation_summary.tsv`
- `empirical/neoaves_chr4/results/structural_source_metadata.tsv`
- `empirical/neoaves_chr4/figures/chr4_structural_track.pdf`
- `empirical/neoaves_chr4/figures/chr4_structural_track.png`

## Stage 2.5: normalized arrangement states

`025_normalize_arrangement_states.py` converts the exact maf2synteny block
signatures from Stage 2 into local arrangement states that are comparable across
species. Exact raw block signatures are too brittle for this purpose because
different genome assemblies can split the same large-scale arrangement into
different numbers of maf2synteny blocks.

The normalization uses only Stage-2 structural tables:
`structural_blocks.tsv`, `structural_intervals.tsv`, and
`arrangement_states.tsv`. It checks that these inputs do not contain Stage-1
QQS/topology columns before running.

For each independently defined structural interval, the script builds reference
coordinate anchors on the GalGal6 chr4 axis. Candidate anchors are derived from
structural block boundaries after merging boundaries within 50 kb and retaining
anchors of at least 100 kb; if no such anchor exists, the full interval is used
as one anchor. Each species' synteny blocks are projected onto these anchors,
recording the observed target order and orientation, missing anchors, and
ambiguous anchors explicitly.

Species are assigned local normalized labels `A0`, `A1`, ... separately within
each structural interval. Grouping uses exact observed shared-anchor
order/orientation compatibility only. Missing anchors do not force different
states, but two species with no comparable observed anchors are not merged.
Ambiguous anchors are retained in the signatures and ignored for grouping.
No phylogenetic topology, QQS, or gene-tree signal is used.

This normalization intentionally does not infer a quartet partition yet. It only
produces structural partitions showing which focal species share each normalized
state inside each structural interval.

Regenerate Stage 2.5 from the repository root with:

```bash
python3 empirical/neoaves_chr4/scripts/025_normalize_arrangement_states.py --run-tests
```

It writes:

- `empirical/neoaves_chr4/results/normalized_arrangement_states.tsv`
- `empirical/neoaves_chr4/results/normalized_arrangement_state_definitions.tsv`
- `empirical/neoaves_chr4/results/arrangement_state_pairwise_similarity.tsv`
- `empirical/neoaves_chr4/results/structural_partitions.tsv`
- `empirical/neoaves_chr4/results/normalized_arrangement_validation.tsv`
- `empirical/neoaves_chr4/figures/chr4_normalized_arrangement_states.pdf`
- `empirical/neoaves_chr4/figures/chr4_normalized_arrangement_states.png`

## Stage 2.6: canonical structural regions

`026_refine_structural_regions.py` refines the Stage-2 species-specific
intervals into canonical physical regions and recomputes arrangement states on
a multi-anchor reference framework. This stage still uses only structural
inputs: `structural_intervals.tsv`, `structural_blocks.tsv`, and
`structural_breakpoints.tsv`.

Species-specific intervals are clustered into canonical GalGal6 chr4 regions
when both boundaries are within 200 kb, or when spans are near-identical by
reciprocal overlap and length ratio. Nested and partially overlapping calls are
preserved unless their spans are nearly the same, so the broad 44-57 Mb region
and its nested/subregion calls remain separate canonical region families.

Within each canonical region, anchors are fixed-size GalGal6/reference bins.
The default is 100 kb:

```bash
python3 empirical/neoaves_chr4/scripts/026_refine_structural_regions.py --run-tests --anchor-size 100000 --min-shared-anchors 5
```

For each species and anchor, the script summarizes overlapping synteny blocks,
choosing dominant query chromosome and orientation by covered base pairs.
Coverage below 20% is marked missing. Ties or weak dominance are marked
ambiguous rather than forced into a state. Arrangement equivalence is called
only when species share at least five informative anchors and have exact
relative anchor-order and orientation concordance. Normalized state groups are
formed by all-pairs agreement within each group, not by transitive closure, so
different missing-anchor sets cannot bridge otherwise discordant arrangements.

In `canonical_region_summary.tsv`, `n_source_species` is the number of species
that contributed a source structural interval to the canonical region, while
`n_species_evaluated` is the number of focal species evaluated for arrangement
state assignment.

The script also repeats normalization at 50 kb, 100 kb, and 200 kb anchor sizes
and reports partition stability. No QQS, topology, locus table, or genealogy
files are read.

It writes:

- `empirical/neoaves_chr4/results/canonical_structural_regions.tsv`
- `empirical/neoaves_chr4/results/canonical_anchor_states.tsv`
- `empirical/neoaves_chr4/results/canonical_region_pairwise_similarity.tsv`
- `empirical/neoaves_chr4/results/canonical_arrangement_states.tsv`
- `empirical/neoaves_chr4/results/canonical_arrangement_state_definitions.tsv`
- `empirical/neoaves_chr4/results/canonical_region_summary.tsv`
- `empirical/neoaves_chr4/results/anchor_size_sensitivity.tsv`
- `empirical/neoaves_chr4/results/canonical_region_validation.tsv`
- `empirical/neoaves_chr4/figures/chr4_canonical_structural_regions.pdf`
- `empirical/neoaves_chr4/figures/chr4_canonical_structural_regions.png`

## Stage 2.7: de novo structural breakpoints

`027_denovo_structural_breakpoints.py` replaces the inferential use of the
previous Stage-2 structural boundaries. The earlier structural tables ultimately
depended on `rearrangements/outlier-regions-by-maf2synteny.tsv`, which was
generated after restricting maf2synteny blocks to previously selected chr4
outlier windows. That creates circularity for an MSRC test. Stage 2.7 therefore
freezes an independently inferred structural breakpoint set from the complete
chromosome-4 maf2synteny output only:

```text
all chr4 synteny blocks -> blind structural breakpoints
```

The inference stage reads only `rearrangements/totab_chr4.csv.xz`. During
breakpoint inference it must not read prior outlier windows, Stage-2 interval
tables, the chr4 QQS locus table, or Stage-3 genealogy change points. This is
enforced in the code by an allowed-input guard, and by keeping the published
boundary comparison in the separate `027b_compare_denovo_to_published.py`
script.

The complete maf2synteny table is parsed by paired `Block` and `sp` accession
rows. Each focal query block is paired to its exact GalGal6 chr4 row, preserving
block IDs, original coordinates, query scaffold/chromosome labels, reference
strand, query strand, and reference-relative orientation. The relative
orientation is `+` when reference and query strands match, and `-` when they
differ. Chicken/GalGal6 is used only as the reference coordinate axis and is not
treated as a focal query arrangement. The first frozen Stage 2.7 result used
query strand alone as orientation and is archived as
`denovo_breakpoint_manifest_invalid_orientation.json`; it should not be used.

Adjacent GalGal6-ordered blocks are first summarized with structural-only
features: reference gap, query gap/jump, same query chromosome, relative
orientation switch, relative-orientation-conditional query order consistency,
overlap, and reference/query gap discrepancy. Consecutive maf2synteny fragments
are then collapsed into larger collinear synteny runs when they share a query
chromosome/scaffold, relative orientation, expected query order, and do not
exceed configurable reference/query gap tolerances. For relative `+` runs, query
coordinates must increase with GalGal6 coordinate; for relative `-` runs, query
coordinates must decrease. The default run used for the corrected frozen table
was:

```bash
python3 empirical/neoaves_chr4/scripts/027_denovo_structural_breakpoints.py --run-tests --max-reference-gap 100000 --max-query-gap 250000 --min-flank-bp 200000 --min-flank-blocks 2
```

The corrected stage separates raw local disruptions from event-scale structural
breakpoints. `denovo_raw_structural_disruptions.tsv` keeps every adjacent
synteny-run transition and its score. A disruption is promoted only when the
mapping regime has enough flanking support on both sides (`--min-flank-bp`,
`--min-flank-blocks`) and shows either strong categorical evidence or exceeds
the selected score threshold. Flank base pairs are the main persistence measure;
the block count is raw maf2synteny block support, not a requirement for multiple
synteny runs. This allows one long internally collinear synteny run with
multiple blocks to support an event-scale boundary while still preventing a
single short flipped block surrounded by the same mapping regime from becoming
a Stage-3A structural breakpoint.

Candidate event-scale breakpoints are scored with chromosome/scaffold switches,
relative-orientation switches, persistent query-order reversals, robustly
scaled query jumps, and robustly scaled reference/query gap discrepancy.
Lenient, default, and stringent thresholds are run without using expected
coordinates or published regions. Nearby default event calls are merged within
species, then clustered across species using a complete-span constraint:
`max(position)-min(position) <= --cross-species-cluster-bp`, default 200 kb.
Each consensus breakpoint retains at most one representative call per species,
chosen by breakpoint score and robustness, with rejected same-species candidates
kept in audit columns. Consensus coordinates are median GalGal6 positions from
the retained one-call-per-species representatives. Confidence classes are
deterministic species-support classes: strong consensus has at least four
species, moderate consensus has two or three species, and species-specific has
one species. Structural evidence and robustness are retained separately.

Before any genealogy comparison, Stage 2.7c freezes three predeclared breakpoint
sets in `denovo_breakpoint_sets.tsv`: `primary` is `n_species_support >= 2`,
`stringent` is `strong_consensus`, and `inclusive` is every event-scale
consensus call.

The script also runs a structural-only sensitivity analysis over reference gap,
query gap, within-species merge distance, cross-species clustering distance, and
breakpoint score threshold. Parameter choice is based only on empirical block
and gap distributions, generic structural consistency, and stability across
settings; it is not chosen to recover any expected physical coordinates.

It writes:

- `empirical/neoaves_chr4/results/denovo_chr4_blocks.tsv`
- `empirical/neoaves_chr4/results/denovo_adjacent_block_features.tsv`
- `empirical/neoaves_chr4/results/denovo_synteny_runs.tsv`
- `empirical/neoaves_chr4/results/denovo_raw_structural_disruptions.tsv`
- `empirical/neoaves_chr4/results/denovo_species_breakpoints_raw.tsv`
- `empirical/neoaves_chr4/results/denovo_species_breakpoints.tsv`
- `empirical/neoaves_chr4/results/denovo_consensus_structural_breakpoints.tsv`
- `empirical/neoaves_chr4/results/denovo_breakpoint_sets.tsv`
- `empirical/neoaves_chr4/results/denovo_structural_segments.tsv`
- `empirical/neoaves_chr4/results/denovo_breakpoint_sensitivity.tsv`
- `empirical/neoaves_chr4/results/denovo_breakpoint_manifest.json`
- `empirical/neoaves_chr4/results/stage27b_vs_stage27c_consensus.tsv`
- `empirical/neoaves_chr4/figures/denovo_chr4_synteny_map.pdf`
- `empirical/neoaves_chr4/figures/denovo_chr4_synteny_map.png`
- `empirical/neoaves_chr4/figures/denovo_consensus_breakpoints.pdf`
- `empirical/neoaves_chr4/figures/denovo_consensus_breakpoints.png`

Only after `denovo_breakpoint_manifest.json` exists should the descriptive
comparison be run:

```bash
python3 empirical/neoaves_chr4/scripts/027b_compare_denovo_to_published.py
```

That comparison may read `rearrangements/outlier-regions-by-maf2synteny.tsv`
and writes `empirical/neoaves_chr4/results/denovo_vs_published_boundaries.tsv`.
It does not modify any de novo breakpoint calls. Stage 2.7 deliberately does
not compare against genealogy change points; Stage 3A should not be rerun until
the de novo structural set has been reviewed.

## Stage 3A: breakpoint enrichment

`03a_breakpoint_enrichment.py` tests whether genealogy changes in the local
`q1/q2/q3` tracks occur unusually close to independently inferred structural
boundaries. This stage is blind to arrangement-state/topology interpretation:
it uses the Stage-1 main-`chr4` locus table to detect genealogy change points,
finalizes and writes those change points, and only then loads the frozen
Stage-2.7c structural breakpoint sets from `denovo_breakpoint_sets.tsv`.

The Stage 3A chronology is fixed and must be preserved:

1. Genealogy change points were detected from `q1/q2/q3` without structural
   data.
2. Structural breakpoints were inferred from complete chr4 synteny without
   genealogy data.
3. The Stage-2.7c primary/stringent/inclusive structural breakpoint sets were
   frozen.
4. Only after those independent tracks existed were genealogy transitions
   compared with structural breakpoints.

The primary genealogy detector is a rolling multivariate contrast. For each
candidate locus, it compares the mean `(q1,q2,q3)` vector in windows immediately
left and right of the locus and scores the Euclidean contrast. The default
window is 50 loci, and local maxima above a median/MAD threshold are retained.
A deterministic binary segmentation on squared Euclidean loss is also run as a
sensitivity method. Nearby genealogy change points are merged within 200 kb by
retaining the strongest point in each cluster.

Structural boundaries are the predeclared Stage-2.7c de novo consensus
breakpoints. The inferential analysis uses the `primary` set
(`in_primary == 1`, `n_species_support >= 2`). Sensitivity analyses use the
`stringent` set (`in_stringent == 1`, strong consensus only) and the
`inclusive` set (`in_inclusive == 1`, every event-scale consensus call,
including species-specific calls). These sets are frozen before any genealogy
comparison and must not be tuned based on enrichment results.

Stage 3A must not read the old Stage-2.6 canonical or published-outlier
boundary files: `canonical_structural_regions.tsv`,
`canonical_structural_boundaries.tsv`, `structural_intervals.tsv`,
`outlier-regions-by-maf2synteny.tsv`, or
`outlier-regions-by-maf2synteny-summary.tsv`. Structural coordinates for this
stage come only from the frozen Stage-2.7c de novo breakpoint files.

The primary null is a circular shift of finalized genealogy change-point
positions along the usable chr4 coordinate range, preserving the spatial
clustering of the genealogy signal while randomizing its alignment to fixed
structural breakpoints. A secondary null samples the same number of
pseudo-change-points from valid chr4 locus positions.

This stage supports only the claim that structural boundaries and genealogy
transitions are spatially associated. It does not test whether arrangement
states predict a particular topology, whether MSRC is favored over
hybridization, or whether structural polymorphism caused the genealogy change.

Regenerate Stage 3A from the repository root with:

```bash
python3 empirical/neoaves_chr4/scripts/03a_breakpoint_enrichment.py --run-tests --permutations 10000
```

It reads the frozen structural input:

- `empirical/neoaves_chr4/results/denovo_breakpoint_sets.tsv`
- `empirical/neoaves_chr4/results/denovo_consensus_structural_breakpoints.tsv`
- `empirical/neoaves_chr4/results/denovo_breakpoint_manifest.json`

It writes:

- `empirical/neoaves_chr4/results/genealogy_change_points_raw.tsv`
- `empirical/neoaves_chr4/results/genealogy_change_points.tsv`
- `empirical/neoaves_chr4/results/breakpoint_enrichment_summary_final.tsv`
- `empirical/neoaves_chr4/results/breakpoint_enrichment_by_structural_set.tsv`
- `empirical/neoaves_chr4/results/breakpoint_enrichment_by_genealogy_method.tsv`
- `empirical/neoaves_chr4/results/breakpoint_enrichment_sensitivity.tsv`
- `empirical/neoaves_chr4/results/structural_breakpoint_genealogy_support.tsv`
- `empirical/neoaves_chr4/results/shared_genealogy_transitions.tsv`
- `empirical/neoaves_chr4/figures/chr4_change_points_blind.pdf`
- `empirical/neoaves_chr4/figures/chr4_change_points_blind.png`
- `empirical/neoaves_chr4/figures/chr4_change_points_vs_denovo_structure.pdf`
- `empirical/neoaves_chr4/figures/chr4_change_points_vs_denovo_structure.png`
- `empirical/neoaves_chr4/figures/breakpoint_enrichment_null_primary.pdf`
- `empirical/neoaves_chr4/figures/breakpoint_enrichment_structural_set_sensitivity.pdf`
