# Anopheles Stage 0R existing analysis audit

Access date: 2026-10-03

## Scope

This audit inspected the existing `empirical/anopheles_2la` and `data/anopheles_2la` files before any new MalariaGEN querying. Existing outputs were not modified.

## Current source data

- Dryad Fontaine et al. 2015 dataset DOI `10.5061/dryad.f4114` is recorded as an original source, but the large source archives have not been downloaded into this repository.
- MalariaGEN Ag3.10 `fontaine-2015-rebuild` is recorded as a candidate source with 72 literature-curated samples.
- The local Stage 1A API provenance records a failed anonymous MalariaGEN API attempt, blocked before sample-level metadata retrieval.
- The repository currently contains synthetic Stage 1B/Stage 2 validation fixtures and outputs, not empirical local tree inference from MalariaGEN sample data.

## Current taxa and sample counts

The only non-empty Anopheles counts currently available in the repository are aggregate public Ag3.10 `fontaine-2015-rebuild` counts:

| taxon | count |
|---|---:|
| arabiensis | 12 |
| coluzzii | 11 |
| gambiae | 26 |
| melas | 4 |
| merus | 9 |
| quadriannulatus | 10 |

The canonical local `sample_manifest.tsv` has zero sample-level rows, so no sample IDs, sex calls, per-sample karyotypes, phased-data flags, or per-sample inclusion decisions are currently frozen.

## Current tree/window inputs

- No empirical Fontaine local tree files were found in the Anopheles data tree.
- No empirical MalariaGEN local trees were inferred.
- `frozen_strict_quartets_stage1a.tsv` is a deterministic empty table with only a header.
- Synthetic validation contains generated tree fixtures and 50-kb-style topology tracks. These are method validation inputs only and are not biological evidence.
- The working 2La interval in the existing config is `2L:20524058-42165532` on AgamP4-style chromosome-arm coordinates, stored as 1-based inclusive.

## Current statistics

Current empirical statistics are metadata-gated and empty:

- strict 2:2 quartet candidates: 0
- direct sample-level 2La karyotype calls: 0
- sample-level sequence availability flags: 0
- sample-level phased availability flags: 0

Synthetic outputs include spatial topology summaries, inside/outside summaries, and power-analysis tables, but these are not empirical Anopheles results.

## Current outputs

Important existing files:

- `empirical/anopheles_2la/results/stage0_data_audit.md`
- `empirical/anopheles_2la/results/stage1a_metadata_karyotype_freeze.md`
- `empirical/anopheles_2la/results/stage1a_quartet_design_summary.md`
- `empirical/anopheles_2la/results/stage1b_analysis_plan.md`
- `data/anopheles_2la/metadata/source_manifest.tsv`
- `data/anopheles_2la/metadata/sample_manifest.tsv`
- `data/anopheles_2la/metadata/region_manifest.tsv`
- `data/anopheles_2la/metadata/ag3_10_fontaine_rebuild_availability.tsv`
- `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`

## Current limitations

- The current analysis is not yet a sample-level MalariaGEN analysis.
- Fontaine-associated sample IDs have not been recovered locally.
- Direct per-sample 2La karyotypes have not been retrieved.
- The existing MalariaGEN API attempt was blocked by access to a GCS config object.
- The local Python environment on 2026-10-03 could not install `malariagen-data` because dependencies attempted source builds and failed against local SDK/LLVM tooling.
- No empirical SNP, genotype, haplotype, or tree data have been queried for 2La.
- Rare-species phased haplotype coverage is not established.

## Use of Fontaine window trees

The repository does not currently contain the original 4,063 Fontaine 50-kb window trees as active empirical inputs. The Anopheles scaffold is based on literature provenance and MalariaGEN metadata targets, with downstream behavior validated only on synthetic tree fixtures.

## Assumptions about 2La arrangement states

Existing files distinguish:

- `A0_homozygous`: standard `2L+a/2L+a`
- `A1_homozygous`: inverted `2La/2La`
- `heterokaryotype`
- `unknown`

Existing aggregate notes provisionally treat `melas` and `quadriannulatus` as fixed standard and `arabiensis` and `merus` as fixed inverted, based on literature-level arrangement status. The scaffold explicitly does not use these aggregate states as direct per-sample karyotype calls.

## Frozen manuscript results

No manuscript-ready empirical Anopheles result is frozen. The existing Anopheles status is a blocked metadata/karyotype freeze plus synthetic validation. It should not be described as a completed empirical 2La analysis.
