# Fire ants chromosome 16 Stage 4B report

## Purpose

Stage 4B tests whether the Stage-4A shift from frozen species-history quartet support toward SB/Sb haplotype-partition support is unusually aligned with the independently frozen chromosome-16 supergene region.

## Frozen inputs

- Stage-0 focal quartet checksum: `5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6`
- Stage-1 manifest checksum: `fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b`
- Stage-2 focal partition checksum: `52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002`
- Stage-3 background partition checksum: `fa4cc04afcc40deb5e6b3d5ca76dc55b3720dfe070f7d90c3939e253b8006286`
- Stage-4A manifest checksum: `d1fa79815bbe4d1f45a2af7f5a1d5487d39fe3041ccac4f79f961a998ab99bf1`
- Stage-4A window-support checksum: `b8ff3ea95ba31ae5553f56f0aeb45a393bdb0ba2b19081fc1aca822dfdd70758`

All recorded Stage-1 raw and processed checksums, Stage-2 outputs, Stage-3 outputs, and Stage-4A outputs were rechecked before inference.

## Physical chromosome-16 ordering

The Stage-1/TWISST source files concatenate chr16 regions as `chr16A`, `chr16B`, then `chr16_supergene`. Stage 4B does not use that row order for spatial inference. It reconstructs the chromosome-16 track by `chrom == chr16` and increasing physical midpoint.

The verified physical order is `chr16A -> chr16_supergene -> chr16B`, with 96 total chr16 windows: 42 `chr16A`, 52 `chr16_supergene`, and 2 `chr16B`.

## Observed Stage-4A contrast

The predeclared response is `D(w)=q_H(w)-q_S(w)`. The observed contrast is `Delta_D = mean(D_supergene) - mean(D_chr16 outside) = 1.610109863495346`.

Component changes are descriptive: `Delta q_H = 0.800247347783683`, `Delta q_S = -0.809862515711663`, and `Delta q_3 = 0.009615167927980`.

## Primary exact circular-shift null

The primary null rotates the complete ordered chr16 `D` vector through all 96 circular alignments while keeping the frozen 52-window supergene mask fixed. This preserves 52 inside and 44 outside windows, all observed `D` values, and the ordered local genealogy track. It destroys only the alignment between that track and the frozen supergene annotation. Ties to the observed statistic are counted with tolerance `1e-12`.

Observed rank among 96 shifts: `2`. `n_ge_observed = 2`. Exact one-sided `p = 0.0208333333`.

## Physical-coordinate sensitivity

The same-width interval uses the frozen Stage-0 supergene span from `data/fire_ants_chr16/metadata/region_manifest.tsv`: start `11680438`, end `27917498`, width `16237060` bp. The chr16 analysis domain is defined conservatively as the minimum source window start and maximum source window end among frozen chr16 Stage-4A windows: `23916` to `28943581`.

The coordinate null enumerates starts induced by fixed-width interval boundary events `s=x_i` and `s=x_i-L`, plus domain boundaries, the observed interval, and deterministic between-event representatives, then deduplicates exact inside-window membership sets. No arbitrary sliding grid is used.

Unique coordinate placements: `97`. Observed rank: `3`. `n_ge_observed = 3`. Coordinate sensitivity `p = 0.0309278351`.

This sensitivity preserves physical interval width and irregularly spaced window coordinates, but the number of sampled windows inside candidate intervals may vary. It is not the primary P-value.

## Boundary context

Nearest-window context around the independently frozen supergene boundaries is reported without fitting a change point or moving boundaries:

- left left: rank 38, window 155, region chr16A, mid 6855778, D -0.959099779552600, dominant species
- left left: rank 39, window 156, region chr16A, mid 6885696, D -1.000000000000000, dominant species
- left left: rank 40, window 157, region chr16A, mid 6914078, D -1.000000000000000, dominant species
- left left: rank 41, window 158, region chr16A, mid 6945212, D -0.960761605110043, dominant species
- left left: rank 42, window 159, region chr16A, mid 7172472, D -1.000000000000000, dominant species
- left right: rank 43, window 162, region chr16_supergene, mid 11700079, D -1.000000000000000, dominant species
- left right: rank 44, window 163, region chr16_supergene, mid 16702415, D -0.181836353944563, dominant species
- left right: rank 45, window 164, region chr16_supergene, mid 17261275, D 0.995681399298905, dominant haplotype
- left right: rank 46, window 165, region chr16_supergene, mid 17679736, D 0.026483050847458, dominant haplotype
- left right: rank 47, window 166, region chr16_supergene, mid 17859820, D 0.966101694915254, dominant haplotype
- right left: rank 90, window 209, region chr16_supergene, mid 27538229, D 0.999873513786997, dominant haplotype
- right left: rank 91, window 210, region chr16_supergene, mid 27612993, D -0.078389830508475, dominant species
- right left: rank 92, window 211, region chr16_supergene, mid 27820687, D 0.471034657222363, dominant haplotype
- right left: rank 93, window 212, region chr16_supergene, mid 27873246, D 1.000000000000000, dominant haplotype
- right left: rank 94, window 213, region chr16_supergene, mid 27904961, D 0.966101694915254, dominant haplotype
- right right: rank 95, window 160, region chr16B, mid 28669799, D -0.970989483574862, dominant species
- right right: rank 96, window 161, region chr16B, mid 28898377, D -0.958462831122836, dominant species

## Results

The primary circular null gives exact one-sided `p = 0.0208333333` with observed `Delta_D = 1.610109863495346`. The coordinate-aware sensitivity gives `p = 0.0309278351` over `97` unique same-width physical placements.

Dominance-pattern summary is descriptive only. Chr16 outside windows are 44/44 species-dominant. Chr16 supergene windows are species=8, haplotype=43, third=1.

## Interpretation

The shift from species-history quartet support toward SB/Sb haplotype-partition support is unusually aligned with the independently defined chromosome-16 supergene region. The supergene region is spatially associated with a coherent alternative genealogy regime relative to collinear chromosome-16 windows.

## Introgression caveat

The source study's recurrent-introgression interpretation remains an explicit caveat. These spatial nulls do not prove MSRC, do not show that inversions caused the genealogy, and do not rule out introgression.

## Limitations

The primary null is in window space. The coordinate-aware sensitivity preserves physical width but changes the number of sampled windows in candidate intervals. Both analyses use the Stage-4A TWISST-derived local support table and do not re-infer trees.

## Stage boundary

Stage 4B stops after spatial-null inference. No ASTRAL/ASTER run was performed, no MSRC model was fit, and no supergene boundary was optimized.
