# Atlantic cod quantitative linkage-validation report

    Stage 9 parsed Matschiner et al. 2022 Source Data Fig. 1 to construct a real physical-coordinate per-SNP linkage track for LG01, LG02, LG07, and LG12. The source linkage statistic is the per-SNP sum of physical distances to strongly linked SNPs (`R^2 > 0.8`) within 250 kb. It is used here as an LD linkage / recombination-suppression proxy, not as a recombination rate.

    ## Source and coordinates

    Source file: `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM3_ESM_source_data_fig1.txt`; SHA256 `5f31d7ebaa30f4f7f145abed365786252aded4c701c4d76ca10d77c057c9b8af`. The linkage block has columns `lg`, `position`, and `linkage`. Coordinates are interpreted as gadMor2 bp positions and are compatible with the frozen cod Stage-1 inversion and 250-kb window coordinates.

    SNP counts: {"LG01": 298, "LG02": 230, "LG07": 249, "LG12": 197}.

    ## SNP-level inside/outside linkage

    | LG | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    | LG01 | 219 | 79 | 6.57408e+07 | 2169 | 6.57387e+07 | 30309.3 |
| LG02 | 199 | 31 | 4.00963e+07 | 508 | 4.00958e+07 | 78929.8 |
| LG07 | 226 | 23 | 5.49202e+07 | 35 | 5.49201e+07 | 1.56915e+06 |
| LG12 | 52 | 145 | 5.32087e+06 | 53371 | 5.2675e+06 | 99.696 |

    ## 250-kb window-level inside/outside linkage

    | LG | n inside windows | n outside windows | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    | LG01 | 0 | 0 | nan | nan | nan | nan |
| LG02 | 1 | 0 | 4.26973e+07 | nan | nan | nan |
| LG07 | 1 | 0 | 5.43937e+07 | nan | nan | nan |
| LG12 | 1 | 0 | 5.32645e+06 | nan | nan | nan |

    ## Spatial linkage-genealogy relationships

    Spearman correlations use 250-kb windows and exclude boundary windows. The Source Data Fig. 1 linkage block is sparse on the chromosome-wide frozen SNAPP grid: most linkage SNPs fall in or near inversion intervals rather than in a dense track across all outside windows. Consequently, linkage-D and linkage-A correlations are reported only when at least three non-boundary 250-kb bins contain linkage SNPs. Circular-shift p-values keep the linkage track fixed, shift D(w) or A(w), and include the identity rotation.

    | LG | comparison | n windows | rho | circular p |
    |---|---|---:|---:|---:|
    | LG01 | D_topology | 0 | nan | nan |
| LG01 | A_time | 0 | nan | nan |
| LG02 | D_topology | 1 | nan | nan |
| LG02 | A_time | 1 | nan | nan |
| LG07 | D_topology | 1 | nan | nan |
| LG07 | A_time | 1 | nan | nan |
| LG12 | D_topology | 1 | nan | nan |
| LG12 | A_time | 1 | nan | nan |

    ## Interpretation

    At the per-SNP level, all four LGs show much higher median linkage scores inside the frozen inversion intervals than outside the intervals represented in Source Data Fig. 1. LG01 and LG02 combine strong quantitative long-range linkage with the previously frozen topology and divergence-time sensitivity signals. LG07 combines strong linkage and divergence-time sensitivity without a corresponding strong topology shift, making it a topology-time discordant inversion rather than a failed control. LG12 shows elevated linkage and supportive topology/time behavior, but Stage-7 physical-coordinate sensitivity was weaker.

    ## Limitations

    The linkage score is an LD-based distance-sum statistic influenced by recombination suppression, population structure, selection, haplotype frequencies, and demography. It is not a direct cM/Mb recombination rate. Individual SNPs are not treated as independent genomic replicates. The 250-kb window summaries are used for physical-coordinate alignment with frozen D(w) and A(w), but the source linkage rows are too sparse on the full SNAPP grid for a meaningful linkage-genealogy correlation test. Spatial correlations are therefore descriptive only and unavailable where fewer than three non-boundary bins have linkage data.
