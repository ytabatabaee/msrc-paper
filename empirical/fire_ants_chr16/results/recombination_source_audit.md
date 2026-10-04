# Recombination/linkage source audit

## Route A: Wang et al. 2013 direct linkage map

- Paper: Wang J., Wurm Y., Nipitwattanaphon M. et al. A Y-like social chromosome causes alternative colony organization in fire ants. Nature 493, 664-668 (2013).
- DOI: 10.1038/nature11832
- Article URL: https://www.nature.com/articles/nature11832
- Supplementary data URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM98_ESM.zip
- Supplementary data file: `data/fire_ants_chr16/raw/recombination_sources/wang_2013_nature/41586_2013_BFnature11832_MOESM98_ESM.zip`
- Supplementary data SHA256: `155f1a9576cf11085c206fcc706ba828683e6a20c18cdfe413ef4a5059ef8613`
- Supplementary information URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM97_ESM.pdf
- Supplementary information file: `data/fire_ants_chr16/raw/recombination_sources/wang_2013_nature/41586_2013_BFnature11832_MOESM97_ESM.pdf`
- Supplementary information SHA256: `70b09fed688ece03cd0d84dc45223f5f6bba64fbb8d1805d57cbe2ff758f2e0b`
- Data type: direct linkage map.
- Marker-level numerical data available: yes.
- Parsed marker rows: 27934 across 7 families (M013, M047, M173, P008, P016, P033, P034).
- Original coordinate system: `Si_gnF.scaffold..._nt...` plus family-specific linkage-group cM.

A documented Wang `Si_gnF` scaffold to frozen Stolle `Si_gnGA`/`gng20170922wFex.fa` chromosome-16 placement table was not recovered. The upstream Stolle README mentions `linkage_map_supergene.txt` and `2018-05-11-linkage-map/results/linkage_map_supergene.txt`, but that file is not committed in the local frozen snapshot or in the public upstream GitHub tree inspected for this analysis. Therefore Route A did not produce target-coordinate cM/Mb values.

## Route B: Yan et al. 2020 physical-coordinate LD

- Paper: Yan Z., Martin S.H., Gotzek D. et al. Evolution of a supergene that regulates a trans-species social polymorphism. Nature Ecology & Evolution 4, 240-249 (2020).
- DOI: 10.1038/s41559-019-1081-1
- Article URL: https://www.nature.com/articles/s41559-019-1081-1
- Supplementary table URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-019-1081-1/MediaObjects/41559_2019_1081_MOESM2_ESM.xls
- Supplementary table file: `data/fire_ants_chr16/raw/recombination_sources/yan_2020/41559_2019_1081_MOESM2_ESM.xls`
- Supplementary table SHA256: `2cc864b2ba12c9be8605b0761aa5ca436ae5954fd8445ba972a85e505a8152e4`
- Supplementary information URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-019-1081-1/MediaObjects/41559_2019_1081_MOESM1_ESM.pdf
- Supplementary information file: `data/fire_ants_chr16/raw/recombination_sources/yan_2020/41559_2019_1081_MOESM1_ESM.pdf`
- Supplementary information SHA256: `18e829cfa04af4597f1b1c0826f449a03508b3aaf7276723f3bc1ef4428bc4de`
- Public raw sequence BioProject: PRJNA421367.

Yan et al. report physical-coordinate LD across chr16 in Extended Data Fig. 5 and provide exact SB-reference inversion breakpoints in the article. The downloadable Supplementary Tables 1-5 contain gene lists, expression tables, and sample metadata, but not the numeric LD matrix or a one-dimensional LD track underlying Extended Data Fig. 5. The BioProject contains large raw sequence data; no small indexed genotype/LD file sufficient to reconstruct chr16 LD was recovered. The Stolle/Wurmlab genome-wide VCF endpoint is approximately 1,503,782,377 bytes and is not the Yan Extended Data Fig. 5 source; it was not downloaded blindly for this figure.

## Evidence hierarchy conclusion

`QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED`. Wang 2013 remains direct experimental linkage support in original coordinates, and Yan 2020 remains published physical-coordinate inversion/LD evidence, but neither yielded a reproducible quantitative chr16 recombination/LD track aligned to the frozen MSRC coordinate axis.
