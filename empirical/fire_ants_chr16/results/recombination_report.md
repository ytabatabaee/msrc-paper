# Recombination-suppression validation report

## Purpose

The requested update was to replace the earlier schematic recombination bar with a quantitative chromosome-16 recombination or LD track on a real physical coordinate axis. The target layout was physical chr16 position, quantitative recombination/LD evidence, and the existing frozen `q_species` / `q_haplotype` genealogy track.

## Route A: direct Wang linkage map

Wang et al. 2013 provide direct linkage-map marker tables. I parsed 27934 RADtag marker rows across 7 mapping families and preserved them in `data/fire_ants_chr16/processed/recombination_map.tsv`. These rows contain original `Si_gnF` scaffold positions and family-specific genetic positions in cM.

Route A did not produce a chr16 cM/Mb track because no documented `Si_gnF` scaffold-to-Stolle `Si_gnGA`/`gng20170922wFex.fa` chromosome-placement file was recovered. The important candidate file named `linkage_map_supergene.txt` is referenced by the upstream README but was not present in the local frozen inputs or the public upstream repository tree. I did not guess scaffold placements.

## Route B: Yan physical-coordinate LD

Yan et al. 2020 report LD r² across physical chr16 and exact SB-reference inversion breakpoints. The published inversion union is chr16:12612565-24031576. I downloaded the minimal journal supplementary PDF and Supplementary Tables 1-5 and inspected them. They do not contain the numeric Extended Data Fig. 5 LD matrix or a one-dimensional LD track. Reconstructing LD from raw PRJNA421367 reads would require large-scale read/genotype processing, and no small public chr16 genotype/LD file was recovered. I did not digitize the published heatmap or treat pixels as quantitative data.

## Result

`QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED`

No quantitative recombination or LD track was generated. The previous constant schematic bar was copied to archive filenames for provenance, but it should not be used as the final quantitative validation figure requested here.

## Relation to frozen genealogy and Stage 6C

The frozen topology/TWISST/ASTRAL/CASTLES-II results remain unchanged. The local-genealogy result still shows a pronounced switch from species-history support outside the frozen supergene analysis span to SB/Sb haplotype support inside it. Stage 6C branch-length results remain unchanged. This audit only addresses whether an independent quantitative recombination/LD track can be reproducibly aligned to those frozen coordinates.

## Introgression caveat

The fire-ant supergene is known to have experienced recurrent adaptive introgression. Recombination suppression helps preserve a long linked haplotype after such events, but this audit does not imply that MSRC without gene flow fully explains the system.

## Smallest missing objects

A quantitative figure would require one of the following: (1) a documented Wang `Si_gnF` scaffold-to-Stolle `Si_gnGA` chr16 placement table such as `linkage_map_supergene.txt`, or (2) the numeric Yan Extended Data Fig. 5 SNP/genotype/LD source data in physical chr16 coordinates, preferably as a chr16 VCF/genotype matrix or precomputed r² table.
