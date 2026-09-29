# Atlantic cod Stage-1 report

Stage 1 converted published processed population-tree source data into a canonical 250-kb window table. No topology enrichment, quartet classification, arrangement-state assignment, SNAPP rerun, or raw-data processing was performed.

## Data acquired

Total downloaded disk size recorded in `download_manifest.tsv`: 352315 bytes.

| file | size bytes | purpose |
|---|---:|---|
| `data/atlantic_cod/raw/github_supergenes/split_vcf.sh` | 3473 | Cod phylogenomics inversion-coordinate provenance. |
| `data/atlantic_cod/raw/github_supergenes/make_snapp_xmls_windows.sh` | 894 | Window-design provenance. |
| `data/atlantic_cod/raw/github_supergenes/make_snapp_xmls_windows.slurm` | 3654 | Window naming, coordinate, sample, and filtering provenance. |
| `data/atlantic_cod/raw/github_supergenes/combine_snapp_results_windows.sh` | 2578 | MCC and posterior tree-file provenance. |
| `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM6_ESM_source_data_fig4.txt` | 230593 | Canonical 250-kb published population-window trees. |
| `data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM4_ESM_source_data_fig2.txt` | 85987 | Baseline population-tree source discovery only. |
| `data/atlantic_cod/raw/zenodo/zenodo_4560275_record.json` | 25136 | Inventory of public MCC/posterior files without archive download. |

## Provenance

Primary tree source: Nature Source Data Fig. 4. Baseline-tree discovery source: Nature Source Data Fig. 2. Zenodo was inspected through record metadata only; no archive or large genotype file was downloaded.

## Window inventory

| LG | windows | fully inside | outside | boundary-overlapping | missing/excluded |
|---|---:|---:|---:|---:|---:|
| LG01 | 110 | 66 | 43 | 1 | 3 |
| LG02 | 92 | 21 | 70 | 1 | 4 |
| LG07 | 119 | 36 | 81 | 2 | 6 |
| LG12 | 105 | 53 | 51 | 1 | 4 |

## Inversion membership

Canonical coordinates are 1-based inclusive. `inside_inversion=true` only when the full window lies within the frozen inversion interval. Boundary-overlapping windows are tracked separately with `overlaps_inversion=true` and `inside_inversion=false`.

## Tree parsing

Parsed successfully: 426 trees. Failures: 0.

## Populations

Canonical labels preserve source labels in Stage 1:

Gadmor_avc_spc, Gadmor_avo_spc, Gadmor_bat_spc, Gadmor_bor_spc, Gadmor_icc_spc, Gadmor_ico_spc, Gadmor_kie_spc, Gadmor_lfc_spc, Gadmor_lfo_spc, Gadmor_low_spc, Gadmor_twc_spc, Gadmor_two_spc

## Baseline tree

Located: true. The tree is recorded for Stage-3 source discovery only.

## Posterior trees

Supergene-wide SNAPP posterior `.trees` files are listed in Zenodo: true. Per-window posterior files are not publicly listed in the Zenodo record inspected during Stage 1.

## Bornholm check

LG12 7.50-7.75 Mb window located: true. Window ID: LG12_007500001_007750000.

## Blockers

- Source Data Fig. 4 does not include per-window `n_sites`, so that field is blank in the canonical table.
- Biological population descriptions for abbreviated source labels remain unresolved in Stage 1.
- Per-window posterior SNAPP samples were not identified in the public Zenodo record; only published source-data trees and supergene-wide posterior files were inventoried.
- Coordinate caution: Stage 1 freezes the `cod_phylogenomics/src/split_vcf.sh` intervals for window-tree analysis. Other source-repository workflows and paper notes, especially around LG02 assembly placement, must not be mixed into this coordinate set silently.
