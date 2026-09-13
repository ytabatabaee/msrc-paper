Regression tests for the Anopheles 2La empirical analysis.

Stage-0 tests are embedded in `scripts/00_audit_samples.py` and run with:

```bash
python3 empirical/anopheles_2la/scripts/00_audit_samples.py --run-tests
```

The current tests validate allowed 2La states, duplicate sample IDs, missing
species labels, coordinate sanity, impossible 2:2 quartet classifications, and
topology/gene-tree column leakage into the structural-only prediction step.
