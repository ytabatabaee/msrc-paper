# Atlantic cod Stage-3 baseline tree provenance

The Stage-3 baseline population tree was extracted programmatically from `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt`.

- Source: Nature Source Data Fig. 2 for Matschiner et al. (2022), panel b, `Maximum-credibility population tree (b)`.
- Source SHA256: `7cb2256e1f48b8e883914569c8da31642f045a31010f232a71c44ebf2ca5ddcf`.
- Canonical Newick path: `data/atlantic_cod/processed/baseline_population_tree.nwk`.
- Canonical Newick SHA256: `e54b7102080b9d36c21ec0c7d003dccbdf353bcefb718cc2abfa4bc0bd7b666a`.
- Method reported by source heading: maximum-credibility population tree.
- Stage-1 provenance identifies this as the outside-supergene / collinear baseline population-tree source.
- Tree summary: fully bifurcating for the taxa present.
- Branch lengths are present in the Newick. Branch support annotations are not present in this source Newick.
- The Newick is represented with a rooted nesting, but Stage 3 treats quartet topology as unrooted.

## Taxa present

Mapped Stage-2 Atlantic cod labels: Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_bor_spc, Gadmor_icc_spc, Gadmor_ico_spc, Gadmor_kie_spc, Gadmor_lfc_spc, Gadmor_lfo_spc, Gadmor_low_spc, Gadmor_twc_spc, Gadmor_two_spc.

Non-Stage-2/outgroup labels: Gadcha_spc, Gadmac_spc, Gadoga_spc.

## Label reconciliation

Baseline tree labels match the Stage-2 `source_label` values for all 12 Atlantic cod populations. Outgroup labels are retained in the canonical Newick and listed in `data/atlantic_cod/metadata/baseline_tree_taxa.tsv`, but they are not used for Stage-2 candidate quartets.

## Anti-circularity

No 250-kb local window tree, Source Data Fig. 4 topology, Bornholm local topology, or local quartet support was read or used during Stage 3.
