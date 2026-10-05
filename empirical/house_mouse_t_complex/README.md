# House mouse t-complex empirical analysis

This directory stages an empirical MSRC analysis of the house-mouse t-haplotype / t-complex on chromosome 17 using the published Kelemen & Vicoso local-gene-tree archive.

The current implementation covers Stage 0 and Stage 1 scaffolding only. Stage 0 audits the source archive, inventories the nested ML tree files, parses the final Newick trees, and creates unresolved tip-mapping tables when the archive does not directly establish biological identities. Stage 1 is gated: it only freezes ASTRAL-ready gene-tree inputs when Stage 0 has source-backed, mostly resolved tip-to-subspecies mappings.

Run Stage 0 with a local copy of the upstream archive:

```bash
python empirical/house_mouse_t_complex/scripts/00_audit_house_mouse_data.py \
  --archive /path/to/IST-2017-78-v1+1_Data.zip
```

The primary dataset is:

```text
Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/
  1-Trees_for_all_5kb_windows/
    ML_trees.zip
```

The PASS-SNP and coverage-filtered RAW-SNP ML archives are inventoried only as possible robustness datasets.

Important limitations are enforced in the scripts:

- The 2017 chr17 window trees are not treated as a genome-wide or outside-inversion background.
- No tip-to-subspecies, population, or t-status mapping is inferred without source-backed evidence.
- Tentative mappings are excluded automatically from ASTRAL mapping files.
- Stage 1 stops unless the Stage 0 mapping gate passes.

Synthetic tests do not require the 144 MB source archive:

```bash
python empirical/house_mouse_t_complex/scripts/00_audit_house_mouse_data.py --run-tests
pytest empirical/house_mouse_t_complex/tests
```

