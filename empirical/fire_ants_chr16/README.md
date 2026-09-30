# Fire ants chromosome 16 social supergene

This directory scaffolds a staged empirical analysis for the chromosome-16 social supergene in fire ants of the genus `Solenopsis`.

Stage 0 is a design and provenance freeze only. It does not parse TWISST weights, inspect local trees, run ASTRAL/ASTER, calculate quartet support, compare topology frequencies, fit MSRC parameters, or draw new biological conclusions.

## Biological system

The focal system is the fire-ant social chromosome described in:

```text
Stolle, E., Pracana, R., Lopez-Osorio, F. et al.
Recurring adaptive introgression of a supergene variant that determines social organization.
Nat Commun 13, 1180 (2022).
DOI: 10.1038/s41467-022-28806-7
```

The public analysis repository is:

```text
https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression
```

The social supergene lies on chromosome 16. The two structural/social haplotypes are `SB` and `Sb`.

## Conceptual role in MSRC paper

The intended use in the MSRC paper is a clean empirical demonstration that a recombination-suppressed structural haplotype defines a coherent alternative genealogy regime across species.

This dataset may also illustrate the difficulty of distinguishing structural-history and introgression/hybridization explanations from quartet patterns alone. It must not be framed as evidence that introgression was absent.

## Public sources

Stage 0 records these sources in `../../data/fire_ants_chr16/metadata/source_manifest.tsv`:

- Stolle et al. 2022, `Nat Commun 13, 1180`, DOI `10.1038/s41467-022-28806-7`.
- Upstream analysis repository `wurmlab/2021-fire-ant-social-supergene-introgression`.
- Wang et al. 2013, `Nature 493, 664-668`, DOI `10.1038/nature11832`.
- Yan et al. 2020, `Nat Ecol Evol 4, 240-249`, DOI `10.1038/s41559-019-1081-1`.

Stage 0 stores only selected small non-topological upstream files under `../../data/fire_ants_chr16/raw/upstream/`. It does not copy the full upstream repository.

## Repository layout

- `scripts/`: executable stage scripts.
- `results/`: stage reports and frozen design summaries.
- `figures/`: reserved for later rendered figures.
- `config/`: reserved for later stage configuration.
- `tests/`: test fixtures or notes if later stages need them.
- `../../data/fire_ants_chr16/`: data tree with `raw/`, `metadata/`, `intermediate/`, and `processed/`.

## Structural haplotypes

Stage 0 freezes the biological haplotype vocabulary independently of local phylogenetic topology:

- `SB`: one social-chromosome haplotype.
- `Sb`: alternative social-chromosome haplotype carrying major rearrangements/inversions.

The published work reports that `Sb` differs structurally from `SB` by multiple large inversions and has strongly suppressed recombination with `SB`. Haplotype state must not be inferred from local topology.

## Author-defined genomic regions

The author-defined region vocabulary comes from the allowed upstream coordinate file `Topology weighting/results/window_coordinates`:

- `chr1`: non-supergene control chromosome.
- `chr16A`: collinear/recombining chromosome-16 region outside the supergene.
- `chr16B`: collinear/recombining chromosome-16 region outside the supergene.
- `chr16_supergene`: author-designated chromosome-16 supergene region.

The first and last BUSCO-window coordinates are treated as observed analysis-window spans, not nucleotide-resolved inversion breakpoints.

## Focal quartet

Stage 0 freezes one canonical focal quartet:

```text
A = invicta/macdonaghi_SB
B = invicta/macdonaghi_Sb
C = richteri_SB
D = richteri_Sb
```

The three possible unrooted resolutions are:

- `species_split`: `AB|CD`
- `haplotype_split`: `AC|BD`
- `third_split`: `AD|BC`

The canonical table is `../../data/fire_ants_chr16/processed/stage0_focal_quartet.tsv`, with checksum stored in `../../data/fire_ants_chr16/processed/stage0_focal_quartet.sha256`.

## Primary future statistic

Later stages will derive per-window support quantities:

```text
q_S(w): support for species_split
q_H(w): support for haplotype_split
q_3(w): support for third_split
```

For windows with complete resolved topology weight, the predeclared identity is:

```text
q_S + q_H + q_3 = 1
```

The primary scalar contrast will be:

```text
D(w) = q_H(w) - q_S(w)
```

Stage 0 does not calculate these quantities.

## Anti-circularity design

The primary future comparison is within chromosome 16:

```text
chr16_supergene
vs.
chr16A + chr16B
```

The predeclared contrast is:

```text
Delta_D = mean(D_supergene) - mean(D_chr16 outside)
```

The chromosome-1 comparison is a secondary descriptive control:

```text
chr16_supergene
vs.
chr1
```

The later primary null will use exact circular shifts of the ordered chromosome-16 `D` track with the author-defined supergene mask fixed. A later physical-coordinate-aware block-placement sensitivity analysis is also predeclared because four-BUSCO-gene windows are irregularly spaced.

## Introgression caveat

The 2022 source paper concludes that the `Sb` supergene originated in the `S. invicta/macdonaghi` lineage and repeatedly introgressed into other species. This analysis must not claim that MSRC explains fire-ant history without gene flow, that the inversion caused introgression, or that topology alone proves MSRC.

## Exploratory-analysis disclosure

Dataset selection and feasibility were exploratory and informed by the published study and an exploratory feasibility calculation of the public TWISST resources.

The formal focal quartet definition is biologically predefined from species identity and `SB`/`Sb` state. The formal pipeline is frozen before reproduction of local support statistics inside this repository.

## Planned stages

Stage 0: Scaffold, provenance, structural definitions, focal quartet, and analysis-design freeze.

Stage 1: Retrieve and normalize the public 213-window coordinate/local-tree/TWISST datasets without testing the focal prediction.

Stage 2: Freeze sample/group-to-species and `SB`/`Sb` mappings independently of local topology support.

Stage 3: Freeze the background species-history topology using the published chromosome 1-15 ASTRAL result, independently of the chromosome-16 local support track.

Stage 4A: First formal local-topology stage: derive `q_species`, `q_haplotype`, `q_third`, and `D` from TWISST weights.

Stage 4B: Spatial-null testing of supergene enrichment on chromosome 16.

Stage 5: Robustness/final manuscript figure and analysis freeze.

Stage 6: ASTER/ASTRAL4 inference-sensitivity analysis: background versus supergene versus all local windows and progressive supergene downweighting.

## Reproduction

From the repository root:

```bash
python3 empirical/fire_ants_chr16/scripts/00_freeze_analysis_design.py --run-tests
```

The script validates the Stage-0 manifests, verifies the upstream coordinate table, counts author-defined regions, validates focal groups and quartet resolutions, writes the focal-quartet checksum, and enforces forbidden-input guards.
