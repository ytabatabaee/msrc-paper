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

Stage 0 records MalariaGEN Ag3 as AgamP4-aligned data, but the exact breakpoint
coordinates trace back to PEST/AgamP3/anoGam3-style reporting. The working 2La
interval is recorded provisionally as:

```text
2L:20524058-42165532
```

This is based on White et al. 2007/Fontaine et al. 2015 table S11 lineage of
coordinates and reuse in Ag1000G examples. Exact AgamP4 equivalence should be
confirmed from a machine-readable source or lift-over before Stage 1.

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
