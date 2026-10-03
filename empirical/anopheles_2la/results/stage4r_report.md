# Stage 4R report — six-species species-tree sensitivity using all Fontaine samples

Stage 4R is additive to Stage 2R and Stage 3R. Stage 2R remains the primary clean arrangement-vs-species test using 53 homozygous individuals and the crossed gambiae/coluzzii design. Stage 4R asks whether the 2La local-history regime alters species-level summary-tree inference when all 72 Fontaine-associated individuals and all six biological species are included.

- individuals: 72
- haplotypes: 144
- biological species: 6
- usable outside windows: 546
- usable inside windows: 430
- total usable non-boundary windows: 976
- Stage 4R classification: **SUPPORT/BRANCH-LENGTH EFFECT**
- finalized correction: Stage 4R downweighting all-inside endpoint is `430`, not the Stage-3R restricted-analysis count of 429

`T_outside6`: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`

`T_inside6`: `arabiensis,coluzzii,gambiae,merus|melas,quadriannulatus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,gambiae|coluzzii,melas,merus,quadriannulatus`

`T_all6`: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`

Branch-level sensitivity is reported in `stage4r_fixed_split_scores.tsv` and `stage4r_branch_sensitivity.tsv`. Stage 4R is the primary Anopheles species-tree sensitivity analysis; Stage 3R remains the restricted four-species homozygote sensitivity analysis matched to the Stage-2R cohort. After correcting the Stage-4R downweighting endpoint label from the Stage-3R carryover value 429 to the dynamic all-inside count of 430, no further Anopheles empirical analysis is currently required.
