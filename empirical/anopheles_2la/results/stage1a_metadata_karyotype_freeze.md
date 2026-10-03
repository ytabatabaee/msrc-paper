# Stage 1A metadata/karyotype freeze

Stage 1A stops after structural/karyotype prediction freezing. It does not infer or inspect local genealogies.

## API retrieval

- MalariaGEN package version: 14.0.0
- Requested release: 3.10
- Sample-set query: fontaine
- Documented sample set identifier: fontaine-2015-rebuild
- Programmatically discovered sample set identifier: fontaine-2015-rebuild
- Retrieval status: SUCCESS
- Samples retrieved: 72

## Reconciliation

- Expected Stage-0 aggregate total: 72.
- Retrieved sample-manifest rows: 72.
- See `empirical/anopheles_2la/results/stage1a_sample_reconciliation.tsv`.

## Karyotype provenance

- Preferred method: `malariagen_data.Ag3().karyotype("2La", sample_sets=<fontaine sample set>)`.
- Species-fixed states are allowed only with `state_basis=species_fixed` and are not mislabeled as direct sample-level calls.

## State counts

- A0_homozygous: 26
- A1_homozygous: 27
- heterokaryotype: 19
- unknown: 0

## Counts by species

- arabiensis: 12
- coluzzii: 11
- gambiae: 26
- melas: 4
- merus: 9
- quadriannulatus: 10

## Evidence basis

- direct_sample_karyotype states: 53
- species_fixed states: 0
- unresolved/excluded samples: 29

## Strict quartet predictions

- strict 2:2 candidate quartets: 42120
- four-distinct-species strict quartets: 0
- geography-matched strict quartets: 630
- frozen file: `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`
- sha256: `d6477c676432ac7a63bca6565e2a3dc20f8a35e80892e1d81a00752816e3b180`

## Best Design A candidates

- none available from authoritative sample-level manifest

## Design B arrangement-replacement contrasts

- none available until polymorphic gambiae/coluzzii sample-level 2La states are retrieved.

## Design C geography-matched contrasts

- none available until sample-level geography and 2La states are retrieved.

## Proceed to Stage 1B?

No. Authoritative sample-level metadata/karyotype retrieval remains blocked in this environment, so no clean strict 2:2 predictions have been populated beyond the deterministic empty freeze.
