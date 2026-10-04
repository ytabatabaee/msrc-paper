# Recombination/linkage source audit

## Selected source

- Paper: Wang J., Wurm Y., Nipitwattanaphon M. et al. A Y-like social chromosome causes alternative colony organization in fire ants. Nature 493, 664-668 (2013).
- DOI: 10.1038/nature11832
- Article URL: https://www.nature.com/articles/nature11832
- Supplementary data URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM98_ESM.zip
- Supplementary information URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11832/MediaObjects/41586_2013_BFnature11832_MOESM97_ESM.pdf
- Download/access date: 2026-10-04
- Supplementary data file: `data/fire_ants_chr16/raw/recombination_sources/wang_2013_nature/41586_2013_BFnature11832_MOESM98_ESM.zip`
- Supplementary data SHA256: `155f1a9576cf11085c206fcc706ba828683e6a20c18cdfe413ef4a5059ef8613`
- Supplementary information file: `data/fire_ants_chr16/raw/recombination_sources/wang_2013_nature/41586_2013_BFnature11832_MOESM97_ESM.pdf`
- Supplementary information SHA256: `70b09fed688ece03cd0d84dc45223f5f6bba64fbb8d1805d57cbe2ff758f2e0b`

## Evidence hierarchy result

Direct linkage-map evidence was found. Wang et al. provide RADtag linkage-map tables with marker identifiers, original scaffold positions, linkage groups, and cM positions for seven mapping families. The parsed marker inventory contains 27934 marker rows across 7 families (M013, M047, M173, P008, P016, P033, P034). Marker-level numerical data are therefore available in the source assembly.

## Coordinate status

The marker positions are reported as `Si_gnF.scaffold..._nt...` scaffold coordinates with linkage-map cM positions. The frozen MSRC fire-ant genealogy track uses Stolle/TWISST chromosome coordinates (`chr16`, bp). I did not find a reliable committed scaffold-to-target-chromosome conversion in the frozen analysis inputs, so the marker-level map is preserved in `data/fire_ants_chr16/processed/recombination_map.tsv` but is not projected onto chromosome 16.

## Data type

- Data type: direct linkage map, represented here as a published regional summary for target-coordinate visualization.
- Coordinate system: Wang 2013 original `Si_gnF` scaffold positions plus linkage-group cM positions.
- Marker-level numerical data available: yes.
- Marker-level target-coordinate data available: no.
