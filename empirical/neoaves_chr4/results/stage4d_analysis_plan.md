# Neoaves chromosome 4 Stage 4D analysis plan

## Question and claim boundary

The published PNAS analysis already established that chromosome-4 outlier loci interact with reduced taxon sampling and can alter species-tree inference. Stage 4D tests whether spatial over-representation helps explain that established effect, and whether topology-neutral physical-bin normalization or independently frozen structural exclusion reduces it. A treatment result will not be described as causal proof.

## Sequential gates

Novel treatments remain disabled until all three reproduction gates pass:

1. **FULL363_ORIGINAL_REPRO:** ASTRAL-MP 5.15.1, the archived collapsed 63,430-tree input, default unrooted search and local posterior support. Compare against the supplied `63K.tre` by unrooted RF, normalized RF, differing splits, and frozen focal relationships.
2. **Figure 4:** recover all 36 original removal conditions, the signed statistic `-(main-alt)`, locus mean, median, and sample-SE aggregation, and the exact published outlier predicate. Validate the reproduced table and figure against the authors' `all.stat.xz`, `removed-count.tsv`, script, and PDF.
3. **Jarvis-48:** recover a source-controlled exact taxon list and reproduce all-locus and published-outlier-excluded ASTRAL behavior. A list inferred from a preferred result is inadmissible.

If a gate is materially inconsistent or its inputs cannot be established, Stage 4D stops before the freeze and reports an incomplete/inconclusive extension without assigning an outcome category. Category E requires an observed support effect without a full-tree topology effect; a failed checkpoint is not evidence for that category.

## Frozen source definitions

- Chromosome-4 classification is the Stage-4C chromosome-group rule: `chr4` and any `chr4_*` sequence are chromosome 4. All are absent from NONCHR4.
- Published PNAS outliers use the authors' exact chicken GalGal4 `chr4` **window-start** predicate: 25,030,000–32,670,000; 33,510,000–34,470,000; and 44,130,000–56,810,000, inclusive. Random/unplaced scaffolds are not PNAS outliers.
- Structural masks are the untouched Stage-4A unions of ±250 kb neighborhoods around the frozen stringent, primary, and inclusive breakpoint sets. They apply by 10-kb window midpoint on anchored `chr4`.
- Locus coordinates are the existing Stage-1–4C 10-kb window coordinates. The selected 1-kb alignment offset is retained separately.

## Novel treatments to freeze after the gates

- **T0:** all usable loci.
- **T_PNAS:** T0 minus loci satisfying the published outlier predicate.
- **T_STRUCT:** T0 minus anchored-chr4 loci whose midpoint is in the frozen primary structural mask. Stringent and inclusive masks are sensitivities.
- **T_BLOCK_500K:** every non-chr4 locus plus one locus from every nonempty, zero-anchored, half-open 500-kb bin on each chromosome-4 sequence. Select the midpoint closest to the bin midpoint, breaking ties by lexical locus ID. Sequence-local coordinates on random scaffolds are never pooled with anchored chr4 coordinates.
- **T_BLOCK_250K / T_BLOCK_1MB:** predeclared scale sensitivities only.
- **T_RANDOM_MATCHED:** 100 deterministic replicates per regime and bin scale, retaining every non-chr4 locus and uniformly selecting exactly the corresponding number of chr4 loci. A fixed master seed and replicate-specific derived seeds will be frozen.

No focal clade, quartet state, reference-tree agreement, or observed outcome enters binning or selection.

## Taxon regimes

- full Stiller sampling (363 taxa);
- exact Jarvis-48 sampling after source validation;
- every condition in the authors' 36-condition Figure-4 removal grid, in its original order.

Actual ASTRAL inference is required for FULL363 and Jarvis-48 under NONCHR4 and the four primary treatments. Additional Figure-4 transition conditions are predeclared from the published experiment before viewing novel-treatment output and run only when compute permits.

## Comparisons

Within each taxon regime, NONCHR4 is an independent genome-wide comparator, not biological truth. Record unrooted RF and normalized RF, split gains/losses, local posterior support changes, normalized quartet score when ASTRAL reports it, and frozen Columbea/N61/N62 relationships. A treatment moves toward NONCHR4 only when RF decreases or a named branch changes from the T0 resolution to the NONCHR4 resolution.

For the score surface, record S2024–J2014 difference, mean, median, sample SE (to reproduce the original descriptive statistic), and absolute mean-minus-median. Spatial normalization is compared with the full deterministic random-thinning distribution; locus SE is not treated as independence-valid uncertainty.

## Controls

The supplied chromosome archive uses the 80K locus set, whereas Stage 4D uses 63K. Its chromosome-only trees are reported as external controls with that dataset difference stated. Comparable 63K chromosome-only and leave-one-chromosome-out trees are inferred when resources permit; chromosomes are included by a predeclared sampling threshold rather than their outcome.

## Software and reproducibility

Use the authors' historical ASTRAL-MP 5.15.1 source at commit `887669674518f38c13c108d006f936cd54d77467`, seed 692, unrooted gene trees, default search, and annotation level 3. The authors' full run used 28 threads and four P100 GPUs. Hardware/thread differences are recorded; inference settings remain identical across treatments. Large immutable sources are ignored by git; scripts, SHA256 manifests, plans, tables, tree outputs, and reports are retained.
