# Fire ants chromosome 16 Stage 3 report

## Purpose

Stage 3 freezes the independent chromosome 1-15 background species-history relationship before the first formal local TWISST support test.

## Prior-stage checksum validation

- Stage 0 focal quartet: `5ea280c1f035d539d82156f201022b2cb4ff1c90bff1a63f93c91a31fcf9d4f6`
- Stage 1 manifest: `fd62bc80d707e89bafc9afe139259b063f6fd48b1435f69e58ffc83c11a14e2b`
- Stage 2 manifest: `2b56ac7ea4a4d0a31588b803b30896522af4e90bc1e91bf23db6d0b9db80e427`
- Stage 2 focal partition: `52c0cea7ccf1622cc7ed3edd85dda5f4c88c253ca782866f7f5a729bab454002`

All recorded prior-stage raw and processed checksums were verified before Stage 3 outputs were written.

## Published background-tree provenance

The frozen upstream `Coalescence-based phylogenies/README.md` states that the authors constructed a tree representing the supergene region and a tree representing the species phylogeny from chromosomes 1 to 15. It records single-copy-gene alignments, RAxML-NG per-gene phylogenies, ASTRAL `5.14.3`, and 100 ASTRAL bootstrap replicates.

Stage 3 uses only the published chromosome 1-15 ASTRAL tree:

```text
data/fire_ants_chr16/raw/upstream/Trees in newick format/Astral.10SNP-genes.chr1-15.nwk
```

Tree tips: `368`.

## Background species relationship

The independent background relationship frozen for the focal comparison is:

```text
invicta/macdonaghi | richteri
```

The chr1-15 tree contains `95` literal `_inv-`/`_mac-` tips and `56` literal `_ric-` tips. These counts are a sanity check only; the Stage-3 partition is frozen from the published chromosome 1-15 species/background interpretation plus Stage-2 group membership.

This is a species-level background statement, not a new representative-individual quartet analysis.

## Mapping to frozen focal quartet

Combining the background species relationship with frozen Stage-2 group membership gives:

```text
(inv_mac_SB, inv_mac_Sb) | (richteri_SB, richteri_Sb)
```

Resulting focal split: `AB|CD`.

Stage-3 background split matches frozen Stage-0 `species_split`: `true`.

The independently published chromosome 1-15 species history supports the predeclared focal species partition AB|CD.

## Anti-circularity

Background species tree: derived from the published chromosome 1-15 ASTRAL analysis.

Supergene/local genealogy: not inspected in Stage 3.

Stage 3 did not read the supergene ASTRAL tree, read TWISST weights, classify grouped TWISST topologies, calculate q-support statistics, select representatives by tree position, compare regions, or run ASTRAL/ASTER.

## Stage boundary

Stage 4A is the first formal unblinded local-topology stage.
