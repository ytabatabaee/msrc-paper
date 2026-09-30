# Fire ants chromosome 16 Stage 0 analysis design

Stage 0 freezes the empirical design before reproducing any local support statistics inside this repository.

## Biological and structural definitions

The focal system is the chromosome-16 social supergene in fire ants of the genus `Solenopsis`.

Primary source:

```text
Stolle, E., Pracana, R., Lopez-Osorio, F. et al.
Recurring adaptive introgression of a supergene variant that determines social organization.
Nat Commun 13, 1180 (2022).
DOI: 10.1038/s41467-022-28806-7
```

Background structural references:

```text
Wang, J., Wurm, Y., Nipitwattanaphon, M. et al.
A Y-like social chromosome causes alternative colony organization in fire ants.
Nature 493, 664-668 (2013).
DOI: 10.1038/nature11832

Yan, Z., Martin, S.H., Gotzek, D. et al.
Evolution of a supergene that regulates a trans-species social polymorphism.
Nat Ecol Evol 4, 240-249 (2020).
DOI: 10.1038/s41559-019-1081-1
```

The structural/social haplotypes are frozen as:

```text
SB = one social-chromosome haplotype
Sb = alternative social-chromosome haplotype carrying major rearrangements/inversions
```

The published work reports that `Sb` differs structurally from `SB` by multiple large inversions and has strongly suppressed recombination with `SB`. Haplotype state is not inferred from local phylogenetic topology.

## Upstream provenance

The upstream public repository is frozen in `../../data/fire_ants_chr16/metadata/upstream_repository_provenance.json`:

```text
repository: https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression
default branch: master
HEAD commit: bdc7823941680fa560e61b6467af9f77950f64b0
retrieval date: 2026-09-30
```

Only small non-topological files were copied. Stage 0 did not retrieve TWISST weights, topology-output trees, input trees, Newick files, local-tree tracks, quartet-support tables, ASTRAL outputs, or sequencing files.

## Frozen region vocabulary

The author-defined region vocabulary is:

```text
chr1
chr16A
chr16B
chr16_supergene
```

Interpretation:

```text
chr1
    non-supergene control chromosome

chr16A
    collinear/recombining chromosome-16 region outside the supergene

chr16B
    collinear/recombining chromosome-16 region outside the supergene

chr16_supergene
    author-designated chromosome-16 supergene region
```

The coordinate table is `../../data/fire_ants_chr16/raw/upstream/Topology weighting/results/window_coordinates`. Its first and last BUSCO-window coordinates are recorded as an author-designated analysis region / observed BUSCO-window span. They are not treated as nucleotide-resolved inversion breakpoints.

The verified Stage-0 counts are recorded in `stage0_region_inventory.tsv`:

```text
chr1              117 windows
chr16A             42 windows
chr16B              2 windows
chr16_supergene     52 windows
```

## Frozen focal quartet

The biologically predefined focal groups are:

```text
inv_mac_SB     = invicta/macdonaghi_SB
inv_mac_Sb     = invicta/macdonaghi_Sb
richteri_SB    = richteri_SB
richteri_Sb    = richteri_Sb
```

The canonical role assignment is:

```text
A = invicta/macdonaghi_SB
B = invicta/macdonaghi_Sb
C = richteri_SB
D = richteri_Sb
```

The three unrooted resolutions are frozen as:

```text
species_split    = AB|CD
haplotype_split  = AC|BD
third_split      = AD|BC
```

The species-history split is:

```text
(invicta/macdonaghi_SB, invicta/macdonaghi_Sb)
|
(richteri_SB, richteri_Sb)
```

The structural/haplotype split is:

```text
(invicta/macdonaghi_SB, richteri_SB)
|
(invicta/macdonaghi_Sb, richteri_Sb)
```

The third split is:

```text
(invicta/macdonaghi_SB, richteri_Sb)
|
(invicta/macdonaghi_Sb, richteri_SB)
```

## Predeclared support quantities

Later stages will derive three support quantities for every local window:

```text
q_S(w): species-history split support
q_H(w): haplotype split support
q_3(w): third split support
```

For windows with complete resolved topology weight:

```text
q_S + q_H + q_3 = 1
```

The predeclared scalar contrast is:

```text
D(w) = q_H(w) - q_S(w)
```

Stage 0 does not calculate `q_S`, `q_H`, `q_3`, or `D`.

## Primary future comparison

The primary comparison is within chromosome 16:

```text
chr16_supergene
vs.
chr16A + chr16B
```

This controls for chromosome-specific effects.

The predeclared primary contrast is:

```text
Delta_D = mean(D_supergene) - mean(D_chr16 outside)
```

The chromosome-1 comparison is secondary descriptive control evidence:

```text
chr16_supergene
vs.
chr1
```

It is not the primary test.

## Later spatial null

The later primary null is:

```text
exact circular shifts of the ordered chromosome-16 D track
with the author-defined supergene mask fixed
```

This null uses all chromosome-16 windows in physical order and preserves the spatial dependence of neighboring local genealogies.

The predeclared sensitivity analysis is:

```text
physical-coordinate-aware block-placement sensitivity
```

This is included because four-BUSCO-gene windows are irregularly spaced. The null design must not be chosen or tuned after seeing the formal Stage-4 support track.

## Anti-circularity guards

Stage 0 freezes:

- the source/provenance record;
- the region vocabulary and observed BUSCO-window spans;
- the `SB`/`Sb` group labels;
- the focal quartet;
- the species, haplotype, and third unrooted splits;
- the primary chromosome-16 contrast and secondary chromosome-1 control;
- the later circular-shift null and coordinate-aware sensitivity.

The Stage-0 script rejects topology-bearing analysis inputs whose names indicate TWISST weights, topology outputs, input trees, Newick files, quartet-support files, QQS files, or precomputed `q1`/`q2`/`q3` support tables. The guard applies to analysis input files and allows documentation text.

## Introgression caveat

This dataset must not be framed as evidence that introgression was absent. The 2022 source paper concludes that the `Sb` supergene originated in the `S. invicta/macdonaghi` lineage and repeatedly introgressed into other species.

The intended MSRC-paper use is:

```text
a clean empirical demonstration that a recombination-suppressed
structural haplotype defines a coherent alternative genealogy regime
across species
```

and potentially:

```text
an empirical illustration of the difficulty of distinguishing
structural-history and introgression/hybridization explanations
from quartet patterns alone
```

Do not claim that MSRC explains the fire-ant history without gene flow, that the inversion caused introgression, or that topology proves MSRC.

## Exploratory-analysis disclosure

Dataset selection and feasibility were exploratory and informed by published results and an exploratory feasibility calculation of the public TWISST resources.

The formal focal quartet definition is biologically predefined from species identity and `SB`/`Sb` state.

The formal pipeline is frozen before reproduction of the local support statistics inside this repository.

## Stage boundary

Stop after Stage 0. Stages 1-6 are documented in the README but are not implemented here.
