# Atlantic cod Stage-1 data provenance

Stage 1 uses already-computed published outputs only. No raw reads, whole-genome VCFs, BEAST, or SNAPP reruns were used.

## Upstream workflow observations

- `split_vcf.sh` freezes the cod phylogenomics inversion intervals as `LG01:9114741-26192386`, `LG02:18489307-24048607`, `LG07:13610591-23019113`, and `LG12:638100-14327837`.
- Other source-repository sub-workflows, especially demography, may use different inversion limits; Stage 1 does not harmonize them.
- The paper also notes assembly-placement uncertainty around part of the LG02 region when comparing gadMor2 and gadMor3. Stage 1 does not resolve this; it uses the `cod_phylogenomics` coordinate set because that is the coordinate system used for the published window-tree workflow normalized here.
- `make_snapp_xmls_windows.sh` sets `window_size=250000`, `min_n_sites=500`, `max_n_sites=1000`, and `min_site_dist=50`.
- `make_snapp_xmls_windows.slurm` encodes window IDs as `LG##_000000001_000250000`, using 1-based inclusive coordinates with zero-padded starts and ends.
- Windows are non-overlapping: each next start is the previous end plus one.
- Incomplete terminal behavior follows the authors' loop over starts while `window_start < last_pos`; Source Data Fig. 4 contains only windows that produced published tree rows.
- For `LG12`, the XML-generation script excludes `Gadmor_lfc1` and `Gadmor_lfc2` before population-tree inference.
- `combine_snapp_results_windows.sh` combines replicate `.trees` files and writes maximum-clade-credibility `.tre` files with TreeAnnotator.
- Source Data Fig. 4 supplies the published window tree Newick strings directly; per-window `n_sites` is not included there.

## Downloaded files

| local path | source | size bytes | SHA256 | purpose |
|---|---|---:|---|---|
| `data/atlantic_cod/raw/github_supergenes/split_vcf.sh` | GitHub mmatschiner/supergenes cod_phylogenomics source | 3473 | `9bd8e6d24869dc298b7cccab8c5eb7e03da6b90864af51801e728688da33bb0a` | Cod phylogenomics inversion-coordinate provenance. |
| `data/atlantic_cod/raw/github_supergenes/make_snapp_xmls_windows.sh` | GitHub mmatschiner/supergenes cod_phylogenomics source | 894 | `a4c10861553dd4a0ecfc8a1f30cdffd6dc7a0b3f6c1c847417bdd134b461a734` | Window-design provenance. |
| `data/atlantic_cod/raw/github_supergenes/make_snapp_xmls_windows.slurm` | GitHub mmatschiner/supergenes cod_phylogenomics source | 3654 | `4548e2cc207863d587304c46c89640e22e349e5b3dd5b13abfb19523d0a4a8c4` | Window naming, coordinate, sample, and filtering provenance. |
| `data/atlantic_cod/raw/github_supergenes/combine_snapp_results_windows.sh` | GitHub mmatschiner/supergenes cod_phylogenomics source | 2578 | `9b7c31a02c9e80988e9437b426adf4345d40eb0c427b59e00a9bec680ed0b9ef` | MCC and posterior tree-file provenance. |
| `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM6_ESM_source_data_fig4.txt` | Nature Source Data Fig. 4 | 230593 | `bd6dfeb9d5f7f97091627fcdbca8db9c01440fefd794bb1e591513fd9a29957c` | Canonical 250-kb published population-window trees. |
| `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt` | Nature Source Data Fig. 2 | 85987 | `7cb2256e1f48b8e883914569c8da31642f045a31010f232a71c44ebf2ca5ddcf` | Baseline population-tree source discovery only. |
| `data/atlantic_cod/raw/zenodo/zenodo_4560275_record.json` | Zenodo record metadata | 25136 | `e706bd38ab872387505db91b33d6b42ff036fda43ec746de173ebe3912cd3167` | Inventory of public MCC/posterior files without archive download. |
