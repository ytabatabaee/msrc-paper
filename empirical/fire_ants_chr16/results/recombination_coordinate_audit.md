# Recombination coordinate audit

## Frozen MSRC fire-ant coordinate system

The frozen local-genealogy analysis uses the author-labeled Stolle/TWISST window-coordinate file recorded in `data/fire_ants_chr16/metadata/region_manifest.tsv`. The chromosome-16 social-supergene analysis interval is `11680438-27917498` bp on `chr16`. These coordinates are an observed BUSCO-window span / author-designated analysis interval, not exact inversion breakpoints.

## Linkage-map coordinate system

Wang et al. 2013 Supplementary Tables 8-14 report RADtag positions in the original assembly as values such as `Si_gnF.scaffold00759_nt19793`, together with linkage groups and cM coordinates.

## Conversion decision

No reliable conversion from the Wang `Si_gnF` scaffold coordinates to the frozen Stolle/TWISST chromosome-16 coordinate system was found in the committed fire-ant analysis inputs. I therefore did not approximate, scale, or otherwise project the marker positions onto chr16. The integrated figure uses a schematic regional recombination-suppression track based on Wang et al.'s published direct linkage-map conclusion rather than a target-coordinate cM/Mb curve.

## Conversion method

None. Marker-level original positions are preserved in `data/fire_ants_chr16/processed/recombination_map.tsv`; target-coordinate visualization is schematic.
