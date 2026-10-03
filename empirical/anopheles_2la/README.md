# Anopheles 2La empirical analysis

This directory scaffolds an empirical MSRC analysis for the Anopheles gambiae
species complex, focused on the 2La inversion on chromosome arm 2L.

The purpose differs from the frozen Neoaves chromosome-4 analysis. Neoaves
currently supports:

```text
structural transitions <-> genealogy-regime transitions
```

The Anopheles 2La analysis is intended to test the theorem-level prediction:

```text
arrangement-state partition -> specific favored local topology
```

Stage 0 does not perform topology tests, fit MSRC, inspect local trees, or draw
biological conclusions.

## Biological target

The focal structural polymorphism is the 2La inversion on chromosome arm `2L`.
The first-pass state vocabulary is:

- `A0`: standard `2L+a` arrangement.
- `A1`: inverted `2La` arrangement.

The primary theorem test should initially use homozygous/homokaryotypic
lineages whenever possible:

- `2L+a / 2L+a` -> `A0_homozygous`
- `2La / 2La` -> `A1_homozygous`

Heterokaryotypes are retained in metadata but excluded from the strict two-state
test unless a later stage explicitly analyzes them.

For four lineages with a clean `A0 A0 A1 A1` pattern, the future prediction is
that the local quartet topology joining the two `A0` lineages and the two `A1`
lineages should receive elevated support inside 2La. Equivalently, if the
arrangement-state roles imply split `AB|CD`, the later statistic will be based
on `q_arrangement = P(AB|CD)` versus the other two quartet topologies.

That signal should eventually be compared among the 2La inversion, the left
flanking region, and the right flanking region. Stage 0 only freezes the
structural/karyotype-derived predictions.

## Layout

- `scripts/`: executable stage scripts.
- `results/`: markdown summaries and final stage reports.
- `figures/`: reserved for later rendered figures.
- `config/`: human-readable configuration.
- `tests/`: shared test documentation; Stage-0 tests are embedded behind
  `--run-tests`.
- `../../data/anopheles_2la/`: data tree with `raw/`, `metadata/`,
  `intermediate/`, and `processed/`.

This mirrors the Neoaves convention while keeping source data outside the code
directory.

## Reference and regions

Stage 0 treated assembly equivalence conservatively because the original
breakpoint structure was reported on PEST/AgamP3 scaffolds. Stage 1A records
explicit Ag1000G/AgamP4-based use of the same 2La interval:

```text
2L:20524058-42165532
```

The coordinates are recorded as 1-based inclusive in TSV manifests. Any later
conversion to 0-based half-open intervals must be explicit and must not silently
shift the stored coordinates.

The left and right flanking regions are represented in
`data/anopheles_2la/metadata/region_manifest.tsv`, but no flank length is frozen
in Stage 0.

## Data sources

Candidate sources identified for Stage 0:

- Dryad Fontaine et al. 2015 dataset, DOI `10.5061/dryad.f4114`.
- MalariaGEN Ag3.10 `fontaine-2015-rebuild`, 72 literature-curated samples.
- MalariaGEN Ag3 sample metadata and `ag3.karyotype("2La")` API for 2La calls.

No large genomic source files are downloaded by Stage 0. Large raw FASTQ, BAM,
CRAM, VCF/BCF, Zarr, and tar archives under `data/anopheles_2la/raw/` are
ignored by `.gitignore`.

## Sample and karyotype definitions

`data/anopheles_2la/metadata/sample_manifest.tsv` is the canonical Stage-0
sample table. It uses one row per biological sample once sample-level metadata
are retrieved.

Controlled `2La_state` values:

- `A0_homozygous`
- `A1_homozygous`
- `heterokaryotype`
- `unknown`

Strict theorem-test candidates require documented `A0_homozygous` or
`A1_homozygous` calls. Do not infer or fabricate karyotypes without a source or
a documented computational procedure.

## Anti-circularity design

`scripts/00_audit_samples.py` reads only sample, species, geography, sequence
availability, and 2La karyotype metadata. It rejects input columns whose names
indicate topology, local-tree, gene-tree, quartet-support, QQS, or Fontaine
topology data. The implied arrangement split is generated solely from
`2La_state`.

This freezes structural predictions before genealogy data are examined,
following the safeguards used in the Neoaves analysis.

## Reproduce Stage 0

From the repository root:

```bash
python3 empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests
```

Outputs:

- `data/anopheles_2la/processed/candidate_strict_quartets.tsv`
- `empirical/anopheles_2la/results/stage0_data_audit.md`

Current Stage 0 reports zero strict quartets because sample-level karyotype
metadata have not yet been retrieved into the sample manifest.

## Stage 1 scope

Stage 1 should retrieve small sample metadata and karyotype tables first,
populate and freeze `sample_manifest.tsv`, rerun Stage 0, and only then inspect
sequence/haplotype/local-tree resources. The first topology test should use the
cleanest frozen strict 2:2 design, not a design chosen after inspecting local
tree support.

## Stage 1A metadata freeze

Stage 1A uses the official MalariaGEN Python API as the authoritative route for
sample identifiers, sample metadata, and 2La karyotypes:

```bash
python3 empirical/anopheles_2la/scripts/01_check_malariagen_access.py
python3 empirical/anopheles_2la/scripts/01_fetch_sample_karyotypes.py --run-tests
python3 empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests --freeze-stage1a
```

In the current environment, `malariagen-data==12.0.1` imports after dependency
repair, but `malariagen_data.Ag3()` cannot read the package config object
`gs://vo_agam_release/v3-config.json` anonymously. The failure is recorded in
`data/anopheles_2la/metadata/stage1a_malariagen_api_provenance.json`.

No sample-level states were fabricated from aggregate counts. The Stage-1A
frozen quartet table is therefore a deterministic empty table until official API
access succeeds.

Stage-1A outputs:

- `data/anopheles_2la/metadata/stage1a_malariagen_api_provenance.json`
- `empirical/anopheles_2la/results/stage1a_sample_reconciliation.tsv`
- `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv`
- `data/anopheles_2la/processed/frozen_strict_quartets_stage1a.sha256`
- `empirical/anopheles_2la/results/stage1a_quartet_design_summary.tsv`
- `empirical/anopheles_2la/results/stage1a_quartet_design_summary.md`
- `empirical/anopheles_2la/results/stage1a_metadata_karyotype_freeze.md`

## Running Stage 1A with MalariaGEN authentication

The blocked local run was an authentication failure, not a biological result and
not a reason to change package versions. Keep `malariagen-data==12.0.1` for the
first authenticated rerun.

From a local machine with Google Cloud SDK installed:

```bash
gcloud auth application-default login
python empirical/anopheles_2la/scripts/01_check_malariagen_access.py
```

Only if the check reaches Ag3 sample-set metadata and reports `SUCCESS`, run:

```bash
python empirical/anopheles_2la/scripts/01_fetch_sample_karyotypes.py --run-tests
python empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests --freeze-stage1a
```

`01_fetch_sample_karyotypes.py` writes timestamped attempt summaries under
`empirical/anopheles_2la/results/provenance/`. Failed attempts do not replace a
previously successful sample manifest or frozen quartet table. Authentication
messages deliberately avoid printing account identifiers, tokens, or credential
paths.

For Google Colab, use
`empirical/anopheles_2la/notebooks/stage1a_metadata_freeze_colab.ipynb`. The
notebook installs the pinned package, uses standard Colab Google
authentication, runs the same check and Stage-1A freeze commands, and packages
only small metadata/result files for transfer back to this repository. It does
not retrieve chromosome-wide SNPs, haplotypes, trees, or topology statistics.

## Stage 1B / Stage 2 synthetic validation

Stage 1A remains waiting for authenticated MalariaGEN metadata access. During
that waiting period, the downstream Stage 1B/Stage 2 machinery has been
implemented and frozen using synthetic fixtures only. No real Anopheles
sequence-derived genealogy/topology data were accessed.

The prospective analysis plan is fixed in
`empirical/anopheles_2la/config/stage1b_analysis_plan.yaml` and summarized in
`empirical/anopheles_2la/results/stage1b_analysis_plan.md`.

Fixed primary regions, AgamP4 2L, 1-based inclusive:

- `left_flank`: `15,524,058-20,524,057`
- `inversion_2La`: `20,524,058-42,165,532`
- `right_flank`: `42,165,533-47,165,532`

The primary flank width is 5 Mb. Sensitivity flank widths are 2 Mb and 10 Mb.
The primary window scheme is fixed 100 kb non-overlapping physical windows, with
50 kb and 200 kb sensitivity analyses.

The core prospective test is not:

```text
Are gene trees unusual inside 2La?
```

It is:

```text
Does a topology predicted independently from structural arrangement state
become specifically enriched inside 2La?
```

The stronger Design-B test is:

```text
When arrangement background changes while species identity is controlled,
does the local genealogy change in the direction predicted from arrangement state?
```

Run the synthetic validation:

```bash
python3 empirical/anopheles_2la/scripts/run_stage1b_synthetic_validation.py --synthetic-validation --permutations 199
python3 empirical/anopheles_2la/scripts/20_power_analysis.py --synthetic-validation --replicates 30 --permutations 49
```

Primary final analysis remains frozen at 10,000 block-aware circular-shift
permutations. The smaller values above are only for repository validation
runtime.

Current synthetic validation summary:

- Null/MSC-like scenario: not positive, `Delta_arr = -0.014`, block-aware
  `p = 0.760`.
- MSRC-positive scenario: positive, `Delta_arr = 0.422`, block-aware `p = 0.005`.
- Broad alternative scenario: not positive for 2La-specific enrichment,
  `Delta_arr = 0.014`, block-aware `p = 0.240`, flagged as broad-alternative-like.
- Design-B synthetic validation: positive replacement contrast changes in the
  predicted direction; null replacement contrast is `null_like`.

Synthetic outputs are under
`empirical/anopheles_2la/results/stage1b_synthetic/` and figures are under
`empirical/anopheles_2la/figures/stage1b_synthetic/`. All reports and power
outputs are labeled `SYNTHETIC_ONLY`.

Freeze manifest:

- `empirical/anopheles_2la/results/stage1b_synthetic/stage1b_analysis_freeze.json`
- Previous SHA256:
  `a4c366a9921ac11d3794bf29d7dbb102a278dda543a9650cd4ede1c985aa8fe2`
- Current SHA256: see
  `64c6b70f9046ea5b3b7f493150dfd68cd6ae454beb513c3147599e4d204a4274`
- Reason for checksum change on 2026-09-16: pre-data hardening:
  design-agnostic execution gate, tree-inference environment pinning, and
  test-only network isolation.

Real Stage 1B is gated. Stage 1B eligibility is design-agnostic. A
scientifically complete Stage 1A may enable Design A, Design B, Design C, any
combination, or no analysis if no prospectively eligible designs exist. Zero
Design-A quartets is not itself a reason to alter the analysis plan.

The local-tree ingestion and tree-inference scripts refuse real-mode execution
until `data/anopheles_2la/processed/stage1a_freeze_complete.json` exists,
Stage-1A freeze/provenance validation succeeds, and at least one prospectively
frozen eligible Design A, Design B, or Design C table or repository-equivalent
row set is nonempty. If all three designs are empty, Stage 1B stops before
topology inspection with `NO_ELIGIBLE_PROSPECTIVE_DESIGNS`; this is
inconclusive / underpowered, not evidence against MSRC.

Future real execution should use the same interfaces:

```bash
python3 empirical/anopheles_2la/scripts/10_ingest_local_trees.py \
  --real-mode \
  --windows <frozen_window_table.tsv> \
  --quartets data/anopheles_2la/processed/frozen_strict_quartets_stage1a.tsv \
  --tree-dir <immutable_local_tree_dir> \
  --out <stage2_spatial_topology_track.tsv>
```

Tree inference, if needed, is kept separate:

```bash
python3 empirical/anopheles_2la/scripts/11_infer_local_trees.py \
  --real-mode \
  --fasta-dir <window_alignment_dir> \
  --windows <frozen_window_table.tsv> \
  --outdir <immutable_tree_dir> \
  --manifest <tree_inference_manifest.tsv>
```

The tree-inference environment is pinned before real sequence access:
`bioconda::iqtree=2.4.0`, executable `iqtree2`, command template
`iqtree2 -s WINDOW.fasta -seed 1729 -nt AUTO -pre OUTPUT_PREFIX -m MFP`. No
bootstrap/support threshold is used for the primary quartet-classification
analysis. Tree estimation reads only alignments/windows and fixed execution
metadata; structural predictions enter only after local tree files are
finalized.

## Stage 1R — MalariaGEN regional SNP-access and homozygote pilot

Stage 1R uses authenticated MalariaGEN Ag3 release 3.10 data from `fontaine-2015-rebuild` via `gs://vo_agam_release_master_us_central1`. The primary pilot cohort is the frozen set of 53 2La homozygotes; 19 heterokaryotypes remain in metadata but are excluded from the unphased diploid pilot.

The fixed pilot queried eight predefined 50-kb windows on 2L, used direct regional `snp_calls` and `biallelic_snp_calls`, built IBS-style diploid distances, generated NJ trees for visualization, and scored crossed gambiae/coluzzii quartet support. No raw reads were downloaded and no full 2L scan was performed.

Stage 1R decision: **PASS**. Mean pilot C inside/outside: 0.100534 / 0.0039476. Mean pilot D inside/outside: 1 / 0.335214.
