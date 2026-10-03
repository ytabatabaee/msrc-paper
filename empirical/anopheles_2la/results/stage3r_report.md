# Stage 3R — species-level summary-tree sensitivity to 2La-linked local genealogies

Stage 3R asked whether hundreds of linked arrangement-associated local genealogies inside 2La propagate into species-level summary-tree inference relative to the collinear 2L background.

- `T_outside`: `(arabiensis,quadriannulatus)|(coluzzii,gambiae)`
- `T_inside`: `(arabiensis,gambiae)|(coluzzii,quadriannulatus)`
- `T_all`: `(arabiensis,quadriannulatus)|(coluzzii,gambiae)`
- final classification: **SUPPORT/BRANCH-LENGTH EFFECT**

`T_inside` differs from `T_outside`, but `T_all` equals `T_outside`. Adding all 2La windows reduces `q(T_outside)` from approximately `0.931` to `0.702` and reduces the requested CU internal branch length from approximately `2.245` to `0.801`. The full combined topology is therefore stable, but the linked inversion windows weaken the support/branch-length profile. Inside-only linked windows nevertheless yield near-unit support for an alternative topology (`localPP=0.99984`), with `q(T_inside | inside)=0.457` versus `q(T_outside | inside)=0.409`.

Interpretation: the extreme local arrangement signal does not overturn the four-species all-window summary topology, but it materially affects summary support and CU branch-length sensitivity. The Stage-2R crossed MalariaGEN result remains the primary Anopheles biological result.

ASTRAL branch lengths, when present, are treated only as summary-coalescent sensitivity metrics. SU lengths are not interpreted as calibrated substitution lengths or divergence times.
