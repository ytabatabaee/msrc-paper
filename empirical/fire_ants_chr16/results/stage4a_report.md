# Fire ants chromosome 16 Stage 4A report

## Purpose

Stage 4A is the first formal unblinded local-topology stage. It classifies the 945 published TWISST group topologies according to the frozen focal quartet and aggregates published TWISST weights into per-window `q_S`, `q_H`, `q_3`, and `D=q_H-q_S` values.

Stage 4A is descriptive only. No permutation test, circular shift, block-placement test, bootstrap hypothesis test, P-value, smoothing, fitted boundary, or ASTRAL/ASTER run was performed.

## Frozen inputs

- Stage 0 focal quartet checksum: `5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6`
- Stage 1 manifest checksum: `fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b`
- Stage 2 focal partition checksum: `52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002`
- Stage 3 background partition checksum: `fa4cc04afcc40deb5e6b3d5ca76dc55b3720dfe070f7d90c3939e253b8006286`

## Topology classification

Each seven-group TWISST topology was treated as unrooted. The three non-focal groups were ignored for the focal quartet, and the induced split among A/B/C/D was classified by unweighted graph-distance four-point sums.

Class counts: species=315, haplotype=315, third=315.

## Window-level support calculation

For each of 213 windows, raw TWISST weights were summed within the three topology classes to give `W_species`, `W_haplotype`, and `W_third`. Normalized supports were calculated as class weight divided by total class weight; every window satisfied `q_S + q_H + q_3 = 1` within `1e-12`.

Dominant classes use a tie tolerance of `1e-12`.

## Descriptive region summaries

- `chr1`: n=117, mean q_S=0.983232450588138, mean q_H=0.008234827754813, mean q_3=0.008532721657050, mean D=-0.974997622833325
- `chr16A`: n=42, mean q_S=0.977689760645278, mean q_H=0.011166451742283, mean q_3=0.011143787612439, mean D=-0.966523308902995
- `chr16B`: n=2, mean q_S=0.976552618264609, mean q_H=0.011826460915760, mean q_3=0.011620920819631, mean D=-0.964726157348849
- `chr16_outside`: n=44, mean q_S=0.977638072355247, mean q_H=0.011196452159259, mean q_3=0.011165475485494, mean D=-0.966441620195988
- `chr16_supergene`: n=52, mean q_S=0.167775556643584, mean q_H=0.811443799942942, mean q_3=0.020780643413474, mean D=0.643668243299358

Dominant-class counts:

- `chr1`: species=117, haplotype=0, third=0, ties=0
- `chr16A`: species=42, haplotype=0, third=0, ties=0
- `chr16B`: species=2, haplotype=0, third=0, ties=0
- `chr16_outside`: species=44, haplotype=0, third=0, ties=0
- `chr16_supergene`: species=8, haplotype=43, third=1, ties=0

## Primary descriptive contrast

Primary chr16 contrast `chr16_supergene_vs_chr16_outside`: `delta_D = 1.610109863495346` (`n_inside=52`, `n_outside=44`).

Secondary descriptive control `chr16_supergene_vs_chr1`: `delta_D = 1.618665866132683` (`n_inside=52`, `n_outside=117`).

## Figures

- `empirical/fire_ants_chr16/figures/fire_ants_chr16_quartet_support.pdf` and `empirical/fire_ants_chr16/figures/fire_ants_chr16_quartet_support.png`
- `empirical/fire_ants_chr16/figures/fire_ants_chr1_quartet_support_control.pdf` and `empirical/fire_ants_chr16/figures/fire_ants_chr1_quartet_support_control.png`
- `empirical/fire_ants_chr16/figures/fire_ants_chr16_D_track.pdf` and `empirical/fire_ants_chr16/figures/fire_ants_chr16_D_track.png`

## Interpretation

The formal Stage-4A question is whether local quartet support shifts from the independently frozen background species partition toward the independently frozen SB/Sb haplotype partition within the author-defined supergene region. This descriptive analysis may demonstrate an association between a recombination-suppressed structural haplotype and a coherent alternative genealogy regime, but it does not identify the historical mechanism responsible for sharing that haplotype across species.

The source paper's introgression caveat remains in force: Stage 4A does not claim that introgression is absent, that MSRC without gene flow explains the fire-ant history, or that topology proves MSRC.

## Limitations

`chr16B` contains only two windows and is reported separately for transparency. The primary outside comparator pools the actual `chr16A + chr16B` windows. Formal spatial inference is reserved for Stage 4B.

## Stage boundary

Stop after Stage 4A. Stage 4B will perform the predeclared spatial-null inference.
