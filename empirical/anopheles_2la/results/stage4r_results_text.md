# Stage 4R results text

The complete Fontaine/MalariaGEN cohort contained 72 individuals and 144 phased haplotypes across six biological species. Usable local haplotype NJ trees were available for 976 non-boundary windows: 546 outside 2La and 430 inside 2La.

`T_outside6` splits: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`.

`T_inside6` splits: `arabiensis,coluzzii,gambiae,merus|melas,quadriannulatus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,gambiae|coluzzii,melas,merus,quadriannulatus`.

`T_all6` splits: `arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus; arabiensis,coluzzii,gambiae|melas,merus,quadriannulatus; arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`.

The RF comparison between `T_outside6` and `T_all6` was 0 (normalized RF 0). Stage 4R classification: **SUPPORT/BRANCH-LENGTH EFFECT**.

The most inversion-sensitive outside-tree branch by outside-to-all support decrease was `arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`. Its baseline quartet support decreased by 0.2068 from outside-only to all windows, and its requested CU length changed by 0.95. Inside-only local histories produced the strongest conflict/weakening for the branches listed in `stage4r_branch_sensitivity.tsv`.

The downweighting experiment retained all usable outside windows and added deterministic subsets of inside-2La windows. It is summarized in `stage4r_downweighting.tsv` and `stage4r_downweighting_summary.tsv`. Dense linked 2La windows are interpreted as local genomic genealogies rather than independent loci.
