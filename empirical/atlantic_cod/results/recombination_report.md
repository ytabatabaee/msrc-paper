# Atlantic cod quantitative linkage-validation report

    Stage 9 parsed Matschiner et al. 2022 Source Data Fig. 1 to construct a real physical-coordinate per-SNP linkage track for LG01, LG02, LG07, and LG12. The source linkage statistic is the per-SNP sum of physical distances to strongly linked SNPs (`R^2 > 0.8`) within 250 kb. It is used here as an LD linkage / recombination-suppression proxy, not as a recombination rate.

    ## Source and coordinates

    Source file: `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM3_ESM_source_data_fig1.txt`; SHA256 `5f31d7ebaa30f4f7f145abed365786252aded4c701c4d76ca10d77c057c9b8af`. The linkage block has columns `lg`, `position`, and `linkage`. Coordinates are interpreted as gadMor2 bp positions and are compatible with the frozen cod Stage-1 inversion and 250-kb window coordinates.

    SNP counts: {"LG01": 298, "LG02": 230, "LG07": 249, "LG12": 197}.

    ## Source sampling design

    The `Linkage per SNP (a)` block is boundary-focused. It captures the sharp linkage contrast around the supergene boundaries rather than providing a dense chromosome-wide SNP track across all frozen 250-kb SNAPP windows. The archived 250-kb bin table is therefore retained for provenance only and is not interpreted as a chromosome-wide linkage track.

    ## Whole-source SNP-level inside/outside linkage

    These source-reproduction summaries demonstrate the abrupt linkage contrast represented in Source Data Fig. 1. They should not be interpreted as unbiased chromosome-wide inside/outside linkage effect sizes because the published source sampling is concentrated around the supergene boundaries.

    | LG | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---:|---:|---:|---:|---:|---:|
    | LG01 | 219 | 79 | 6.57408e+07 | 2169 | 6.57387e+07 | 30309.3 |
| LG02 | 199 | 31 | 4.00963e+07 | 508 | 4.00958e+07 | 78929.8 |
| LG07 | 226 | 23 | 5.49202e+07 | 35 | 5.49201e+07 | 1.56915e+06 |
| LG12 | 52 | 145 | 5.32087e+06 | 53371 | 5.2675e+06 | 99.696 |

    ## Boundary-relative analysis

    The primary revised analysis uses a prospective ±250 kb window around each frozen inversion boundary. Relative positions are `SNP position - boundary position`. For left boundaries, negative positions are outside and positive or zero positions are inside. For right boundaries, negative or zero positions are inside and positive positions are outside.

    | LG | boundary | n inside | n outside | median inside | median outside | difference | ratio |
    |---|---|---:|---:|---:|---:|---:|---:|
    | LG01 | left | 204 | 33 | 6.62701e+07 | 7584 | 6.62625e+07 | 8738.14801 |
| LG01 | right | 15 | 46 | 1.87488e+07 | 1205.5 | 1.87476e+07 | 15552.6985 |
| LG02 | left | 137 | 18 | 4.26228e+07 | 604.5 | 4.26222e+07 | 70509.2142 |
| LG02 | right | 62 | 13 | 2.50221e+07 | 423 | 2.50216e+07 | 59153.7943 |
| LG07 | left | 102 | 11 | 6.14563e+07 | 88 | 6.14562e+07 | 698366.784 |
| LG07 | right | 124 | 12 | 5.43085e+07 | 18.5 | 5.43085e+07 | 2935596.14 |
| LG12 | left | 2 | 145 | 5501 | 53371 | -47870 | 0.103070956 |
| LG12 | right | 0 | 0 | nan | nan | nan | NA |

    Boundary coverage is complete and strongly contrasting for LG01, LG02, and LG07. LG12 is less completely sampled at the frozen boundaries: the left boundary has very few inside-side SNPs within ±250 kb, and the published linkage block has no SNPs within ±250 kb of the frozen right boundary. LG12 is therefore interpreted primarily from the whole-source inside/outside linkage contrast, not from a two-sided frozen-boundary transition.

    ## Relationship to frozen topology and divergence-time results

    The revised manuscript-facing summary does not use sparse 250-kb linkage-D or linkage-A correlations. Instead, it presents source-level SNP linkage validation alongside already-frozen linkage-group-level topology and divergence-time summaries. LG01 and LG02 combine strong linkage with robust topology shifts and positive divergence-time shifts. LG07 shows strong published linkage evidence and a divergence-time shift despite little topology enrichment, making it a topology-time discordant inversion. LG12 shows strong whole-source linkage contrast and supportive topology/time effects, while retaining the Stage-7 caveat that physical-coordinate temporal sensitivity is weaker and the boundary-coverage caveat noted above.

    | LG | whole-source SNP linkage ratio | ΔD | topology p | ΔA | time circular p | time physical p |
    |---|---:|---:|---:|---:|---:|---:|
    | LG01 | 30309.3 | 1.23277 | 0.009174 | 0.75703 | 0.00917431 | 0.0217391 |
| LG02 | 78929.8 | 0.892332 | 0.010989 | 0.653994 | 0.010989 | 0.0133333 |
| LG07 | 1.56915e+06 | -0.005981 | 0.555556 | 0.435501 | 0.00854701 | 0.011236 |
| LG12 | 99.696 | 1.04317 | 0.038462 | 0.324552 | 0.00961538 | 0.0892857 |

    ## Retired sparse-grid correlations

    `recombination_spatial_correlations.tsv` remains archived for provenance, but it is uninformative because the source linkage rows populate too few non-boundary windows on the full 250-kb SNAPP grid. NA correlations are not interpreted as negative results.

    ## Interpretation

    Published SNP-level linkage data show abrupt, strong increases in long-range linkage at well-covered Atlantic cod supergene boundaries, independently validating substantial recombination suppression. LG12 has incomplete right-boundary coverage in the source linkage block, so its strongest quantitative support is the whole-source inside/outside linkage contrast rather than a two-sided right-boundary transition. The phylogenetic consequences differ among supergenes: LG01 and LG02 show strong topology and divergence-time effects, LG07 shows a divergence-time effect without a strong topology shift, and LG12 shows supportive topology/time effects. Strong recombination suppression therefore does not imply a single uniform phylogenetic outcome.

    ## Limitations

    The linkage score is an LD-based distance-sum statistic influenced by recombination suppression, population structure, selection, haplotype frequencies, and demography. It is not a direct cM/Mb recombination rate. Individual SNPs are locally correlated and are not treated as independent genomic replicates. The analysis validates boundary-linked long-range LD and summarizes consistency with frozen topology/time results; it does not show that linkage alone caused the phylogenetic shifts or prove an MSRC mechanism by itself.
