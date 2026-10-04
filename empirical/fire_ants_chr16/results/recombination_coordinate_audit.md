# Recombination coordinate audit

## Frozen MSRC / Stolle coordinate system

The frozen fire-ant local-genealogy analysis uses Stolle et al. 2022 workflow coordinates from the `Si_gnGA` / `gng20170922wFex.fa` chromosome-level reference. The frozen chromosome-16 supergene analysis interval is `chr16:11680438-27917498`. This is an author-defined BUSCO-window analysis span, not an exact inversion boundary.

## Wang et al. 2013 coordinate system

Wang linkage-map markers are reported on old `Si_gnF` scaffolds, for example `Si_gnF.scaffold00759_nt19793`, with family-specific genetic positions in cM. The Stolle upstream README states that a file named `linkage_map_supergene.txt` / `input/gngs_linkage_map.txt` should contain fields such as `scaffold`, `scaffold_start`, `scaffold_end`, `chr`, `chr_start`, `chr_end`, `orientation`, and possibly `region`, but this placement file was not recovered from the local frozen inputs or public upstream GitHub tree. No Wang marker was accepted as mapped to chr16.

## Yan et al. 2020 coordinate system

Yan et al. report physical chr16 coordinates on the S. invicta SB reference and provide these inversion breakpoints:

- In(16)1: chr16:14549064-24031576
- In(16)2: chr16:13705210-24030990
- In(16)3: chr16:12612565-13683100

The union is chr16:12612565-24031576. These coordinates are useful biological annotations, but without the underlying LD values/genotypes or a documented coordinate equivalence/liftover to the Stolle `Si_gnGA` axis, they are not sufficient to generate a quantitative LD track aligned to the frozen MSRC windows.

## Overlay decision

Assembly compatibility checks did not pass for a quantitative overlay. No direct cM/Mb track and no LD r² track were overlaid on the frozen genealogy track. The previous schematic figure was archived and should not be treated as the requested quantitative validation figure.
