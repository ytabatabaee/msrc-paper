# Atlantic cod Stage-2 report

Stage 2 freezes population identities, ecotype labels, and chromosomal arrangement states from structural metadata only.

## Population-label resolution

| tree label | canonical population | location | ecotype |
|---|---|---|---|
| `Gadmor_avc_spc` | More_stationary | Møre | stationary |
| `Gadmor_avo_spc` | More_migratory | Møre | migratory |
| `Gadmor_bat_spc` | Labrador | Labrador | stationary |
| `Gadmor_bor_spc` | Bornholm_Basin | Bornholm Basin | stationary |
| `Gadmor_icc_spc` | Iceland_stationary | Iceland | stationary |
| `Gadmor_ico_spc` | Iceland_migratory | Iceland | migratory |
| `Gadmor_kie_spc` | Kiel_Bight | Kiel Bight | stationary |
| `Gadmor_lfc_spc` | Lofoten_stationary | Lofoten | stationary |
| `Gadmor_lfo_spc` | Lofoten_migratory | Lofoten | migratory |
| `Gadmor_low_spc` | Suffolk | Suffolk | stationary |
| `Gadmor_twc_spc` | Newfoundland_stationary | Newfoundland | stationary |
| `Gadmor_two_spc` | Newfoundland_migratory | Newfoundland | migratory |

All 12 Stage-1 tree labels are Atlantic cod population/ecotype labels; no outgroup labels occur among them.

## Arrangement orientation

| LG | gadMor2 | gadMor_Stat | evidence | confidence |
|---|---|---|---|---|
| LG01 | derived | ancestral | outgroup_and_contig_alignment | direct_outgroup_alignment |
| LG02 | ancestral | derived | outgroup_colinearity_and_assembly_comparison | direct_outgroup_alignment |
| LG07 | derived | ancestral | outgroup_and_contig_alignment | direct_outgroup_alignment |
| LG12 | ancestral | derived | demographic_inference_after_contig_mapping_was_uninformative | demographic_inference |

## Population arrangement states

| population | LG01 | LG02 | LG07 | LG12 |
|---|---|---|---|---|
| More_stationary | ancestral | derived | ancestral | ancestral |
| More_migratory | derived | ancestral | ancestral | ancestral |
| Labrador | ancestral | ancestral | derived | ancestral |
| Bornholm_Basin | ancestral | ancestral | ancestral | ancestral |
| Iceland_stationary | ancestral | ancestral | derived | ancestral |
| Iceland_migratory | derived | ancestral | derived | ancestral |
| Kiel_Bight | ancestral | derived | ancestral | derived |
| Lofoten_stationary | ancestral | derived | derived | derived |
| Lofoten_migratory | derived | ancestral | derived | ancestral |
| Suffolk | ancestral | derived | ancestral | derived |
| Newfoundland_stationary | ancestral | ancestral | derived | ancestral |
| Newfoundland_migratory | derived | ancestral | derived | ancestral |

## Strict 2:2 candidates

| LG | candidates |
|---|---:|
| LG01 | 168 |
| LG02 | 168 |
| LG07 | 210 |
| LG12 | 108 |

No topology support, local tree grouping, or inside/outside enrichment was calculated.

## Ambiguities

- LG12 orientation is weaker than LG01, LG02, and LG07 because contig mapping was uninformative; the Stage-2 orientation follows the authors' demographic inference and Supplementary Table 5 footnote.
- Source Data Fig. 3 was downloaded during discovery but excluded from Stage-2 structural inputs because it contains topology.
- No arrangement state is inferred from local tree clustering or from the Bornholm LG12 topology switch.

## Anti-circularity declaration

No local tree topology, Source Data Fig. 4 tree grouping, quartet support, or window topology information was used to construct the Stage-2 arrangement predictions.
