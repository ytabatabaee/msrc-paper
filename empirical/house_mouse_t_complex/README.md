# House mouse t-complex empirical analysis

This directory stages an empirical MSRC analysis of the house-mouse t-haplotype / t-complex on chromosome 17 using the published Kelemen & Vicoso local-gene-tree archive.

Stages 0 and 1 audit and freeze the published trees. Stage 2 adds exact fixed-quartet scoring, balanced-sampling ASTRAL controls, seven-population sensitivity, filtering robustness, spatial thinning, and a validated main figure. Stage 2B adds split-specific ASTRAL annotation correction, arrangement-state quartet decomposition, 1,000-replicate exact balanced resampling, and the complete spatial phase sweep.

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

Stage 2B outputs are in `results/stage2b_report.md`. The visualization and
recombination extension is generated with:

```bash
python empirical/house_mouse_t_complex/scripts/02l_reconstruct_recombination_track.py \
  --archive /path/to/IST-2017-78-v1+1_Data.zip
python empirical/house_mouse_t_complex/scripts/02m_make_visualizations.py
```

It produces four-group and seven-population ASTRAL drawings, direct quartet
tracks, the source-defined phylogeny-based recombination-state track, the
aligned q/recombination figure, controls, and
`figures/house_mouse_t_complex_main_v2.{png,pdf}`. The recombination track is
an inferred local-tree classification, not a direct cM/Mb rate estimate. The
Stage 2 scripts use repository-relative output paths and retain source archive
provenance by basename, SHA256, and byte size.
