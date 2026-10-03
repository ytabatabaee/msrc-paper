# Stage 1R pilot SNP filtering plan

This filter was frozen for the Stage 1R pilot before interpreting distance or tree topology results.

Input data are direct MalariaGEN Ag3 regional calls from `biallelic_snp_calls`, restricted to the 53 frozen homozygous Fontaine-associated samples. Raw-read processing is not used.

Filter:

- use biallelic SNPs only;
- exclude invariant sites by requiring minor allele count > 0;
- require minor allele count >= 2;
- require site missingness <= 0.25;
- record sample missingness and flag windows if any sample has missingness > 0.50;
- do not LD-prune for this initial within-window genealogy/distance pilot;
- do not choose, remove, or reorder windows based on topology.

Pairwise distance formula:

For each retained biallelic SNP, convert an unphased diploid genotype to alternate-allele dosage 0, 1, or 2. For samples `i` and `j`, compare only sites callable in both samples and compute

```text
d(i,j) = mean_s |dosage_i(s) - dosage_j(s)| / 2
```

The denominator 2 scales the per-site distance to [0,1].

Pre-filter and post-filter counts by window are recorded in `stage1r_window_qc.tsv`.
