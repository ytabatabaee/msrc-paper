# Recombination-suppression validation report

## Source study

The validation uses an independent direct linkage-map source: Wang J., Wurm Y., Nipitwattanaphon M. et al. A Y-like social chromosome causes alternative colony organization in fire ants. Nature 493, 664-668 (2013). DOI:10.1038/nature11832. This source is distinct from the Stolle et al. local-tree/TWISST analysis used for the frozen MSRC fire-ant genealogy results.

## Data type

Direct linkage-map marker data were recovered from the Wang et al. supplementary data archive. The parsed source tables contain 27934 RADtag marker rows across 7 mapping families. Each row includes a marker, an original scaffold position, a linkage group, and a cM coordinate. These marker-level data are written to `data/fire_ants_chr16/processed/recombination_map.tsv` in their original coordinate system.

## Coordinate compatibility

The linkage-map marker positions use original `Si_gnF` scaffold coordinates. The frozen MSRC fire-ant topology track uses Stolle/TWISST chromosome coordinates. Because no reliable scaffold-to-chromosome conversion was found in the committed inputs, no marker-level cM/Mb curve was projected onto chromosome 16. The figure therefore uses the published regional conclusion as a schematic validation track over the frozen author-designated supergene interval `11680438-27917498` bp.

## Recombination/linkage result

Wang et al. report a large social-chromosome region of approximately 13 Mb, about 55% of the chromosome, in which recombination is completely suppressed between the SB and Sb social chromosomes. This is direct linkage-map evidence for suppressed recombination across the social-supergene region.

## Relation to the frozen genealogy signal

The frozen MSRC fire-ant analysis shows that chr16 windows outside the supergene are species-history dominated, whereas windows inside the author-designated supergene interval shift strongly toward the cross-species SB/Sb haplotype quartet. The independent linkage-map evidence supports the biological consistency of this result: the genomic interval with the strong social-haplotype genealogy is also the known recombination-suppressed social chromosome region.

## Relation to Stage 6C branch-length effects

Stage 6C remains unchanged. Its fixed-topology individual-level ASTRAL4/CASTLES-II comparison found strong CULength reductions for the two focal SB/Sb species-pair branches when supergene windows were added, with SULength responses that differed by branch. The recombination validation does not estimate branch lengths and does not reinterpret CU or SU values as recombination rates. It supports the narrative that linked histories in a recombination-suppressed region can influence summary-tree branch estimates even when the global topology remains stable.

## Limitations

The marker-level linkage map could not be placed onto the frozen chromosome-16 coordinate axis without a documented coordinate conversion. The integrated figure is therefore a regional validation figure, not a new recombination-rate map. It should not be read as estimating local cM/Mb values across the Stolle/TWISST windows.

## Introgression caveat

The fire-ant supergene literature invokes recurrent adaptive introgression among socially polymorphic species. This validation supports the role of recombination suppression in maintaining a long linked genealogy, but it does not show that recombination suppression alone caused the observed genealogy or that MSRC without gene flow explains the system.
