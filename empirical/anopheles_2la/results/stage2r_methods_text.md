# Stage 2R methods text

We analyzed MalariaGEN Ag3 release 3.10 `fontaine-2015-rebuild` processed SNP calls for the 53 frozen 2La-homozygous Fontaine-associated samples. We used direct regional `snp_calls` with local biallelic filtering from `gs://vo_agam_release_master_us_central1`; no raw reads were downloaded or reprocessed.

Chromosome arm 2L was divided into non-overlapping 50-kb windows anchored to coordinate 1. The frozen 2La interval was `2L:20524058-42165532`. Windows fully contained in the interval were classified as inside, windows fully outside as outside, and breakpoint-overlap windows as boundary and excluded from primary inside-vs-outside tests.

For each window we retained biallelic SNPs with minor allele count >= 2 and site missingness <= 0.25. Windows with fewer than 25 retained SNPs or maximum sample missingness > 0.5 were retained in the grid but marked low signal.

The primary crossed analysis used only gambiae and coluzzii because both species contain both homozygous arrangements. Pairwise distances used unphased diploid alternate-allele dosage: `d(i,j)=mean_s |g_i(s)-g_j(s)|/2`. Primary statistics were the crossed mean-distance contrast `C(w)` and the continuous quartet margin `M(w)=mean_q[S_species-S_arrangement]`. Discrete topology support `D(w)=q_arrangement-q_species` was retained as secondary.
