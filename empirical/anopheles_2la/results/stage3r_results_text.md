# Stage 3R results text

The empirical collinear baseline `T_outside` was `(arabiensis,quadriannulatus)|(coluzzii,gambiae)`. The inside-only topology was `(arabiensis,gambiae)|(coluzzii,quadriannulatus)`, and the all-window topology was `(arabiensis,quadriannulatus)|(coluzzii,gambiae)`. `T_all` matched `T_outside`, whereas `T_inside` differed from `T_outside`.

Quartet support for the outside-baseline topology declined from `0.931104` in the outside-only analysis to `0.701540` when all linked 2La windows were added, a relative reduction of `24.7%`. Combined alternative quartet support increased from `0.0689` outside to `0.2985` in the all-window analysis. The requested CU internal branch length also dropped from `2.24532` to `0.801251`, a reduction of `64.3%`. CU lengths are interpreted only as summary-coalescent sensitivity metrics, not calibrated times.

Using explicit split-aware notation, `q(T_outside | outside)= 0.931104`, `q(T_outside | inside)= 0.408833`, and `q(T_outside | all)= 0.701540`. The inside-only analysis inferred the alternative topology `(arabiensis,gambiae)|(coluzzii,quadriannulatus)` with `q(T_inside | inside)= 0.456944` and near-unit ASTRAL local posterior support (`localPP=0.99984`); the remaining third quartet had support `0.134223`. Thus the inside-winning quartet frequency is only modestly higher than the outside-baseline quartet frequency on the same inside-only data. This illustrates sensitivity to linked-window replication: many physically linked 50-kb windows inside 2La can yield near-unit summary support even when the leading raw quartet-frequency advantage is modest.

The downweighting experiment measured whether linked 2La windows could shift species-level summary-tree inference when all outside windows were retained. The inferred topology matched `T_outside` for every mixed treatment, including all inside windows: 1: 1.00, 2: 1.00, 5: 1.00, 10: 1.00, 20: 1.00, 40: 1.00, 80: 1.00, 160: 1.00, all: 1.00. The support trend is summarized in `stage3r_downweighting_summary.tsv` and plotted in the Stage 3R figure.

Stage 3R classification: **SUPPORT/BRANCH-LENGTH EFFECT**.

This does not replace the Stage-2R primary result. Stage 2R remains the main Anopheles result: within 2La, same-arrangement samples across species are closer than same-species samples carrying opposite arrangements.
