# Stage 6C mapping recovery report

## Supplementary Data 1 source

Source URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-28806-7/MediaObjects/41467_2022_28806_MOESM4_ESM.xlsx
Local file: `data/fire_ants_chr16/raw/upstream/nature_communications_2022/41467_2022_28806_MOESM4_ESM_supplementary_data_1.xlsx`
SHA256: `603578cc5c701d9e5955ee091ddded1f65db9895ccd9033894ea2f2e441f45c6`
Worksheet names: `Sheet1`
Rows after header: `386`
Columns: `Sample name`, `Country`, `Latitude`, `Longitude`, `Species`, `Colony Identity`, `Caste`, `Supergene Variant`, `Source`, `BioProject`, `SRA identifier`, `Number of sequence reads (M)`, `Read length`, `Gigabase`, `Map Reads %`, `Genome coverage`, `RFLP result (colony)`, `Gp-9 variant 1`, `Gp-9 variant 2`, `Gp-9 variant 3`, `Gp-9 variant 4`, `Notes`, `Fraction of sites missing`, `Used for Twisst`, `Used for Quibl (invicta and richteri)`, `Used for Quibl (all species)`, `Used for BBBAA-BBABA (invicta and richteri)`

## Upstream grouping rule

The authoritative grouping rule is the rule in `Topology weighting/supergene/filter_samples_by_missing.R`: remove the `S.` species prefix and spaces, merge `invicta` and `macdonaghi` into `invicta/macdonaghi`, and append `Supergene.Variant` for `invicta/macdonaghi` and `richteri`. The other retained species remain species-only groups.

## Matching result

Tree tips inspected: `267`
Matched tips: `267`
Unmatched tips: `0`
Match-type counts: `{'exact': 267}`
Documented name corrections used: `[]`
Species assignments: `267`
Supergene-variant assignments: `267`
Complete exact individual-to-group mapping available: `true`

Group sizes:
- `geminata`: `1`
- `saevissima`: `11`
- `pusillignis`: `2`
- `invicta/macdonaghi_SB`: `134`
- `invicta/macdonaghi_Sb`: `59`
- `richteri_SB`: `28`
- `richteri_Sb`: `32`

## Validation against tree-label bigB/littleb strings

The tree-label strings were used only as a validation check after metadata-based assignment, not as the mapping source.
bigB/littleb validation mismatches: `0`

## Output files

- Tip inventory: `data/fire_ants_chr16/processed/stage6c_tree_tip_inventory.tsv`
- Match audit: `data/fire_ants_chr16/processed/stage6c_tree_tip_supplement_match_audit.tsv`
- Individual-to-group table: `data/fire_ants_chr16/processed/stage6c_individual_to_group.tsv`
- ASTRAL mapping file: `data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt`
