# Fire ants chromosome 16 Stage 2 report

Stage 2 freezes the biological TWISST group definitions and SB/Sb state provenance. It does not classify topologies, calculate local quartet support, compare regions, or run ASTRAL/ASTER.

## Biological group provenance

The frozen upstream script `Topology weighting/supergene/filter_samples_by_missing.R` reads `tmp/supergene.missing`, reads the first 22 columns of `samples_overview.csv`, renames the missingness sample column to `Sample.name`, and merges the two tables by `Sample.name`.

For the TWISST sample table, the script filters to `F_MISS < 0.05`, excludes `S. interrupta`, `xAdR`, and `S.megergates`, removes `S.` and spaces from `Species`, merges `invicta` and `macdonaghi` into `invicta/macdonaghi`, and constructs `clade` as `paste(Species, Supergene.Variant, sep="_")` only for `invicta/macdonaghi` and `richteri`. Other species use `Species` alone as `clade`.

Thus the four focal group labels are determined by species membership plus the `Supergene.Variant` metadata field. They are not derived from local RAxML clustering, TWISST weights, ASTRAL trees, q-values, or dominant topology.

## Frozen TWISST groups

```text
geminata
saevissima
pusillignis
invicta/macdonaghi_Sb
invicta/macdonaghi_SB
richteri_Sb
richteri_SB
```

## Focal species partition

```text
(invicta/macdonaghi_SB, invicta/macdonaghi_Sb)
|
(richteri_SB, richteri_Sb)
```

## Focal haplotype partition

```text
(invicta/macdonaghi_SB, richteri_SB)
|
(invicta/macdonaghi_Sb, richteri_Sb)
```

## Metadata availability

The exact published group-construction algorithm is available in the frozen upstream script.

The exact 267-row source `samples_overview.csv` / grouped `results/twisst_samples_low_missingness` file is not deposited in the audited repository path at the frozen commit. Stage 2 therefore does not fabricate a sample-to-TWISST-group table.

Stage 4 can still operate on the already-published grouped TWISST weights from Stage 1.

## Sample-label audit

Unique Stage-1 sample labels audited: `267`.

Label parse status counts: `{'literal_only': 155, 'unknown': 112}`.

Literal variant-string counts: `{'bigB': 176, 'littleb': 91}`.

Strings such as `bigB` and `littleb` are retained as literal sample-label metadata only. Stage 2 uses the upstream `Supergene.Variant` field as the authoritative SB/Sb provenance and does not silently equate sample-name strings to SB/Sb.

## Anti-circularity

No TWISST weight row, grouped topology tree, focal quartet support value, regional enrichment result, ASTRAL tree, or local topology classification was used to assign group identity or SB/Sb state.
