# msrc
Multi-species rearrangement coalescent model

## Running Anopheles 2La Stage 1A with MalariaGEN authentication

The Stage 1A metadata freeze keeps `malariagen-data==12.0.1` for the first
authenticated rerun. From the repository root:

```bash
gcloud auth application-default login
python empirical/anopheles_2la/scripts/01_check_malariagen_access.py
```

Only after the check reaches Ag3 sample-set metadata and reports `SUCCESS`:

```bash
python empirical/anopheles_2la/scripts/01_fetch_sample_karyotypes.py --run-tests
python empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests --freeze-stage1a
```

Google Colab fallback instructions are in
`empirical/anopheles_2la/notebooks/stage1a_metadata_freeze_colab.ipynb`.

## Anopheles 2La Stage 1B synthetic validation

The Stage 1B/Stage 2 downstream analysis has been prospectively implemented
using synthetic inputs only, while Stage 1A waits for authenticated MalariaGEN
metadata access. No real Anopheles sequence-derived genealogy/topology data were
accessed.

Run:

```bash
python3 empirical/anopheles_2la/scripts/run_stage1b_synthetic_validation.py --synthetic-validation --permutations 199
python3 empirical/anopheles_2la/scripts/20_power_analysis.py --synthetic-validation --replicates 30 --permutations 49
```

Previous freeze checksum:
`a4c366a9921ac11d3794bf29d7dbb102a278dda543a9650cd4ede1c985aa8fe2`.
Current freeze checksum:
`64c6b70f9046ea5b3b7f493150dfd68cd6ae454beb513c3147599e4d204a4274`.
The checksum changed on 2026-09-16 only for pre-data hardening:
design-agnostic execution gate, tree-inference environment pinning, and
test-only network isolation.

Real Stage 1B is blocked until
`data/anopheles_2la/processed/stage1a_freeze_complete.json` exists after a
successful Stage 1A metadata/karyotype freeze and the Stage-1A
freeze/provenance validation succeeds. Stage 1B eligibility is design-agnostic.
A scientifically complete Stage 1A may enable Design A, Design B, Design C, any
combination, or no analysis if no prospectively eligible designs exist. Zero
Design-A quartets is not itself a reason to alter the analysis plan. If all
three designs are empty, Stage 1B returns `NO_ELIGIBLE_PROSPECTIVE_DESIGNS`
before topology inspection; this is inconclusive / underpowered, not evidence
against MSRC.

The prospective tree-inference environment is pinned before real sequence
access as `bioconda::iqtree=2.4.0`, executable `iqtree2`, command template
`iqtree2 -s WINDOW.fasta -seed 1729 -nt AUTO -pre OUTPUT_PREFIX -m MFP`.
