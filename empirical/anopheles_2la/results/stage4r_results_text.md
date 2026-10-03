# Stage 4R results text

The complete Fontaine/MalariaGEN cohort contained 72 individuals and 144 phased haplotypes across six biological species. Usable local haplotype NJ trees were available for 976 non-boundary windows: 546 outside 2La and 430 inside 2La.

`T_outside6` splits: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`.

`T_inside6` splits: `arabiensis,coluzzii,gambiae,merus|melas,quadriannulatus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,gambiae|coluzzii,melas,merus,quadriannulatus`.

`T_all6` splits: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`.

The RF comparison between `T_outside6` and `T_all6` was 0 (normalized RF 0). Stage 4R classification: **SUPPORT/BRANCH-LENGTH EFFECT**.

The most inversion-sensitive outside-tree branch by outside-to-all support decrease was `arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`. Its baseline quartet support decreased by 0.2068 from outside-only to all windows, and its requested CU length changed by 0.95. Inside-only local histories produced the strongest conflict/weakening for the branches listed in `stage4r_branch_sensitivity.tsv`.

Using all 72 Fontaine-associated individuals represented by 144 phased haplotypes, the inside-2La summary tree differed substantially from the outside-2La tree (RF=4; normalized RF=0.667), changing two of three internal splits. When all 430 usable inside windows were combined with the 546 outside windows, the global six-species topology returned to the outside topology, but all three outside-tree internal branches showed reduced quartet support and shorter requested CU branch lengths.

For the strongest branch, `arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`, quartet support changed from `0.871768` outside to `0.664997` in all windows, and requested CU length changed from `1.63610` to `0.686122`. This is not interpreted as a whole-tree topology failure because `T_all6=T_outside6`.

The downweighting experiment retained all usable outside windows and added deterministic subsets of inside-2La windows. The final endpoint is dynamically labelled as all 430 usable inside windows, giving 546+430=976 total non-boundary windows. It is summarized in `stage4r_downweighting.tsv` and `stage4r_downweighting_summary.tsv`. Dense linked 2La windows are interpreted as local genomic genealogies rather than independent loci.
