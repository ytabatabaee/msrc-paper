# Stage 2R preanalysis plan — full 2L 50-kb scan, not executed

Stage 2R will use the same authenticated MalariaGEN Ag3 release `3.10` data resource and the same frozen 53 homozygous Fontaine-associated samples from Stage 1R.

Predeclared design:

- construct non-overlapping 50-kb windows across chromosome arm 2L in AgamP4 coordinates;
- classify windows fully inside frozen 2La interval `2L:20524058-42165532` as inside;
- classify windows fully outside the interval as collinear outside controls;
- exclude any boundary-overlap window from inside/outside tests;
- use the same biallelic SNP filter as Stage 1R unless a change is documented before analysis: minor allele count >= 2, site missingness <= 0.25, no LD pruning for the primary distance/tree scan;
- mark low-SNP windows as missing rather than replacing them after inspecting topology;
- primary statistic: crossed gambiae/coluzzii `C(w)` and quartet `D(w)` using only samples from the four crossed classes;
- secondary analysis: all-53 descriptive distances and NJ visualization summaries;
- spatial null/control: compare inside 2La windows against all eligible collinear 2L windows, with sensitivity to matched controls by SNP density/callability only if the matching rule is fixed before topology inspection;
- fixed seed `20261003` for any deterministic quartet subsampling, although exhaustive crossed quartets are preferred;
- no raw-read processing;
- no sample removal based on local topology.

Stage 2R is not executed in Stage 1R.
