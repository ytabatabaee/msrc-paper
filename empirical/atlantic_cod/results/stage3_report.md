# Atlantic cod Stage-3 report

Stage 3 freezes the outside-supergene baseline population-history split for every frozen Stage-2 structural quartet. It does not inspect local 250-kb window trees.

## Baseline tree

- Source: Nature Source Data Fig. 2, panel b, `Maximum-credibility population tree (b)`.
- Source path: `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt`.
- Source SHA256: `7cb2256e1f48b8e883914569c8da31642f045a31010f232a71c44ebf2ca5ddcf`.
- Canonical Newick path: `data/atlantic_cod/processed/baseline_population_tree.nwk`.
- Tree taxa: Gadcha_spc, Gadmac_spc, Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_bor_spc, Gadmor_icc_spc, Gadmor_ico_spc, Gadmor_kie_spc, Gadmor_lfc_spc, Gadmor_lfo_spc, Gadmor_low_spc, Gadmor_twc_spc, Gadmor_two_spc, Gadoga_spc.
- Polytomies: no.
- Branch support: not present in the source Newick.
- Rooting: stored as supplied, but all induced quartet splits are treated as unrooted.
- Label reconciliation: all 12 Stage-2 Atlantic cod labels are present directly in the baseline tree; three non-Stage-2 outgroup labels are retained only in the canonical tree/taxon inventory.

## Quartet design

| LG | total Stage-2 quartets | arrangement = baseline | arrangement != baseline | unresolved | fraction informative |
|---|---:|---:|---:|---:|---:|
| LG01 | 168 | 56 | 112 | 0 | 0.666667 |
| LG02 | 168 | 91 | 77 | 0 | 0.458333 |
| LG07 | 210 | 81 | 129 | 0 | 0.614286 |
| LG12 | 108 | 92 | 16 | 0 | 0.148148 |

## Manual parser/design sanity examples

These examples were selected deterministically from the first resolved rows of the Stage-3 table, not from local window trees.

| LG | quartet | four taxa | arrangement split | baseline split | informative |
|---|---|---|---|---|---|
| LG01 | LG01_SQ0001 | Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_ico_spc | Gadmor_avc_spc,Gadmor_bat_spc|Gadmor_avo_spc,Gadmor_ico_spc | Gadmor_avc_spc,Gadmor_ico_spc|Gadmor_avo_spc,Gadmor_bat_spc | true |
| LG01 | LG01_SQ0002 | Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_lfo_spc | Gadmor_avc_spc,Gadmor_bat_spc|Gadmor_avo_spc,Gadmor_lfo_spc | Gadmor_avc_spc,Gadmor_lfo_spc|Gadmor_avo_spc,Gadmor_bat_spc | true |
| LG01 | LG01_SQ0003 | Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_two_spc | Gadmor_avc_spc,Gadmor_bat_spc|Gadmor_avo_spc,Gadmor_two_spc | Gadmor_avc_spc,Gadmor_avo_spc|Gadmor_bat_spc,Gadmor_two_spc | true |

## Ambiguities

No Stage-2 quartet was impossible to assign from the baseline tree. No induced Stage-2 quartet was unresolved. LG12 arrangement orientation remains the Stage-2 caveat and is not changed in Stage 3.

## Freeze checks

- Stage-2 structural candidate checksum verified: `303d9eb563fc6b5c0dd953b0d9f4c4b363baf83061f32fa37fb2f4e4a0f1a72c`.
- Stage-3 baseline quartet checksum: `5af7eea33f8c998556a140a4e14facb87322e58043632e7b04c15de3eab64bc1`.

## Anti-circularity statement

No 250-kb local window tree, Source Data Fig. 4 topology, Bornholm local topology, or local quartet support was read or used during Stage 3.
