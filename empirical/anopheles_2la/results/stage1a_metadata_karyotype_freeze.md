# Stage 1A metadata/karyotype freeze

Stage 1A stops after structural/karyotype prediction freezing. It does not infer or inspect local genealogies.

## API retrieval

- MalariaGEN package version: 12.0.1
- Requested release: 3.10
- Sample-set query: fontaine
- Documented sample set identifier: fontaine-2015-rebuild
- Programmatically discovered sample set identifier: not discovered
- Retrieval status: blocked
- Samples retrieved: 0

## Reconciliation

- Expected Stage-0 aggregate total: 72.
- Retrieved sample-manifest rows: 0.
- See `empirical/anopheles_2la/results/stage1a_sample_reconciliation.tsv`.

## Karyotype provenance

- Preferred method: `malariagen_data.Ag3().karyotype("2La", sample_sets=<fontaine sample set>)`.
- Species-fixed states are allowed only with `state_basis=species_fixed` and are not mislabeled as direct sample-level calls.

## State counts

- A0_homozygous: 0
- A1_homozygous: 0
- heterokaryotype: 0
- unknown: 0

## Counts by species

- none: API retrieval blocked before sample-level metadata were returned

## Evidence basis

- direct_sample_karyotype states: 0
- species_fixed states: 0
- unresolved/excluded samples: 0

## Strict quartet predictions

- strict 2:2 candidate quartets: 0
- four-distinct-species strict quartets: 0
- geography-matched strict quartets: 0
- frozen file: `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`
- sha256: `b19c6ac1987b0f141b5f8c89d51a12cb8f47b125015d105b39d0054cb8850cbd`

## Best Design A candidates

- none available from authoritative sample-level manifest

## Design B arrangement-replacement contrasts

- none available until polymorphic gambiae/coluzzii sample-level 2La states are retrieved.

## Design C geography-matched contrasts

- none available until sample-level geography and 2La states are retrieved.

## Proceed to Stage 1B?

No. Authoritative sample-level metadata/karyotype retrieval remains blocked in this environment, so no clean strict 2:2 predictions have been populated beyond the deterministic empty freeze.

## Retrieval blocker

- HttpError: Anonymous caller does not have storage.objects.get access to the Google Cloud Storage object. Permission 'storage.objects.get' denied on resource '//storage.googleapis.com/projects/_/buckets/vo_agam_release/objects/v3-config.json' (or it may not exist)., 401
