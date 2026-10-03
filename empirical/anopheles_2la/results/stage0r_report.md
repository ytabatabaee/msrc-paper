# Anopheles Stage 0R report

Access date: 2026-10-03

## Summary

Stage 0R does not justify replacing the existing Fontaine-only scaffold yet. The strongest path is a MalariaGEN-centered redesign, but the decisive sample-level pieces were not recovered locally: original/sample IDs, authoritative 2La karyotypes, and Fontaine rebuild haplotype membership.

The best current classification is **B/C: partial multi-species analysis likely feasible, but Fontaine-only historical comparison remains required until sample-level MalariaGEN access succeeds**. Ag3.10 clearly contains a `fontaine-2015-rebuild` sample set with 72 literature-curated samples, but Stage 0R could only freeze aggregate counts.

## Answers to Stage-0R questions

1. Can the original Fontaine samples be recovered exactly?

Not yet in this local audit. Ag3.10 explicitly includes `fontaine-2015-rebuild`, described as literature-curated Fontaine samples, but exact original sample IDs were not recovered because the local API/package route failed.

2. How many samples are available per species?

From the Ag3.10 public release summary for `fontaine-2015-rebuild`: `arabiensis` 12, `coluzzii` 11, `gambiae` 26, `melas` 4, `merus` 9, and `quadriannulatus` 10. `A. christyi` is not listed.

3. How many have 2La states?

Zero authoritative sample-level 2La states were recovered in Stage 0R. All focal samples are counted as `unknown` in `stage0r_sample_counts.tsv`.

4. How many are homozygous standard/inverted?

Authoritative sample-level counts are zero for both classes in Stage 0R because the karyotype table was not retrieved. Literature-level fixed states are noted in existing files, but they were not substituted for per-sample calls.

5. Are phased haplotypes available?

Ag3 haplotype reference panels are documented for `gamb_colu_arab`, `gamb_colu`, and `arab`, with two homologues represented by phased haplotypes. Fontaine rebuild membership and rare-species support were not verified locally.

6. Is the full six-species analysis feasible?

Not proven. The six focal species are present in aggregate in `fontaine-2015-rebuild`, but full six-species sample-level analysis requires sample IDs, clear 2La states, homozygous counts, and queryable variation. Those are not yet recovered.

7. Is raw-read processing unnecessary?

Probably yes, if MalariaGEN regional SNP/haplotype access works as documented. Raw-read reprocessing should not be needed for Stage 1R.

8. Can we use 50-kb regional queries directly?

The API and download documentation indicate regional SNP and haplotype access is available. This was not locally verified because `malariagen-data` could not be installed/imported in the current environment.

9. What should the primary arrangement-vs-species statistic be?

Use `D(w) = Q_arrangement(w) - Q_species(w)` for inferred local trees/quartet weights. If starting from distance matrices, use the difference between same-species/opposite-arrangement distances and cross-species/same-arrangement distances, with all pairs frozen before computing distances.

10. What should the outside-inversion control be?

Use eligible collinear 2L 50-kb windows outside 2La, preferably all windows passing callable-fraction, SNP-density, and missingness filters. If matching/downsampling is needed, use fixed seed `20261003`.

11. Should the MalariaGEN analysis replace the existing analysis?

Not yet. It should supersede the Fontaine-only scaffold only after Stage 1R recovers sample-level metadata, authoritative 2La states, and regional SNP/haplotype access for a usable homozygous sample set.

12. What should Stage 1R do?

Stage 1R should run in a working MalariaGEN environment, preferably Google Colab or an environment with compatible wheels and authenticated/public GCS access. It should retrieve `fontaine-2015-rebuild` sample metadata, 2La karyotypes, haplotype membership, and a tiny 2L regional SNP/haplotype smoke test. It should not infer local trees.

## Files created

- `empirical/anopheles_2la/results/stage0r_existing_analysis_audit.md`
- `data/anopheles_2la/metadata/stage0r_malariagen_source_manifest.tsv`
- `data/anopheles_2la/processed/stage0r_sample_inventory.tsv`
- `empirical/anopheles_2la/results/stage0r_sample_counts.tsv`
- `empirical/anopheles_2la/results/stage0r_2la_karyotype_audit.md`
- `empirical/anopheles_2la/results/stage0r_2la_karyotype_samples.tsv`
- `empirical/anopheles_2la/results/stage0r_haplotype_availability.tsv`
- `data/anopheles_2la/processed/stage0r_2la_interval.tsv`
- `empirical/anopheles_2la/results/stage0r_data_availability.tsv`
- `data/anopheles_2la/processed/stage0r_fontaine_window_inventory.tsv`
- `empirical/anopheles_2la/results/stage0r_preanalysis_plan.md`
- `empirical/anopheles_2la/results/stage0r_old_vs_new_design.tsv`
- `empirical/anopheles_2la/results/stage0r_report.md`

## Computational scale

The working 2La interval is `2L:20524058-42165532`, length 21,641,475 bp. On a 1-based 50-kb grid anchored at chromosome coordinate 1, 432 windows are fully inside 2La and 434 intersect it. A breakpoint-anchored tiling would have 433 intervals, with the last interval clipped.

Expected Stage 1R scale is modest for metadata and smoke tests. Full Stage 2R local genealogy inference across 432 inside windows plus matched outside windows may involve dozens to hundreds of individuals, or twice as many haplotypes if phased data are used. Regional SNP queries should avoid raw-read reprocessing. Local execution is plausible for metadata and pilot windows; full tree inference across many windows should be benchmarked and may be better on HPC if using many samples or large haplotype panels.

## Recommendation

Proceed to Stage 1R only after fixing the MalariaGEN access environment. The recommended primary sample set is the homozygous subset of `fontaine-2015-rebuild` if sample-level 2La states are recoverable; otherwise use a two-tier plan with a partial `gambiae`/`coluzzii`/`arabiensis` MalariaGEN replication and retain Fontaine-only historical comparison for the rare species.
