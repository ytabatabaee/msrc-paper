# Atlantic cod linkage source audit

- Paper: Matschiner et al. 2022, Supergene origin and maintenance in Atlantic cod, Nature Ecology & Evolution 6:469-481
- DOI: 10.1038/s41559-022-01661-x
- Source-data URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-022-01661-x/MediaObjects/41559_2022_1661_MOESM3_ESM.txt
- Download date: 2026-10-04 UTC
- Source file: `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM3_ESM_source_data_fig1.txt`
- SHA256: `5f31d7ebaa30f4f7f145abed365786252aded4c701c4d76ca10d77c057c9b8af`
- Column names by block: `{"Genetic distance (a)": ["lg", "window_center", "window_start", "window_end", "distance_between_coastal_and_first", "distance_between_haddock_and_first"], "Linkage per SNP (a)": ["lg", "position", "linkage"]}`
- Linkage groups present in linkage block: LG01, LG02, LG07, LG12
- SNP rows per LG: {"LG01": 298, "LG02": 230, "LG07": 249, "LG12": 197}
- Coordinate assembly: gadMor2, as stated for the published Source Data Fig. 1 / Matschiner et al. cod workflow and consistent with the frozen cod Stage-1 coordinate system.
- Linkage-score units: sum of physical distances in bp to SNPs with pairwise R^2 > 0.8 within 250 kb, following the published description. This is an LD linkage score and recombination-suppression proxy, not a meiotic recombination rate.
- Missing values across parsed source blocks: 0.

Parsing decision: the file has two explicit blocks. `Genetic distance (a)` contains window-level genetic-distance columns and is retained for audit context. `Linkage per SNP (a)` has explicit columns `lg`, `position`, and `linkage`; Stage 9 uses this block for the quantitative per-SNP linkage track.
