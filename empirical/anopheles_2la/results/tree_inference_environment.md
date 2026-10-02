# Tree-Inference Environment Freeze

This policy is frozen before authenticated access to real Anopheles
sequence-derived genealogy or topology data.

- Software: IQ-TREE 2.
- Exact pinned package: `bioconda::iqtree=2.4.0`.
- Executable name: `iqtree2`.
- Environment file: `empirical/anopheles_2la/config/tree_inference_environment.yml`.
- Version verification before first real run: run `iqtree2 --version` inside the
  pinned environment and require the reported version to begin with
  `IQ-TREE multicore version 2.4.0` or the package-equivalent `2.4.0` string.
- Primary command template:

```bash
iqtree2 \
  -s WINDOW.fasta \
  -seed 1729 \
  -nt AUTO \
  -pre OUTPUT_PREFIX \
  -m MFP
```

- Model-selection policy: IQ-TREE ModelFinder Plus (`-m MFP`) for every window,
  fixed prospectively and not chosen from Anopheles topology results.
- Seed: `1729`.
- Thread policy: `-nt AUTO`; wall-clock parallelism is an execution detail, but
  the per-window command records the resolved command and version.
- Bootstrap/support policy: no bootstrap, SH-aLRT, UFBoot, or support threshold
  is used in the primary quartet-classification analysis. Gene trees are not
  converted to unresolved based on branch support in the primary analysis.
- Expected output files per window: `OUTPUT_PREFIX.treefile`,
  `OUTPUT_PREFIX.iqtree`, `OUTPUT_PREFIX.log`, `OUTPUT_PREFIX.ckp.gz`, and
  `OUTPUT_PREFIX.model.gz` when IQ-TREE emits the model checkpoint.
- Prediction blindness: tree estimation consumes only the alignment/window path,
  fixed software configuration, seed, thread policy, and output prefix. It does
  not read predicted arrangement splits, Design-A/B/C topology predictions, or
  q1/q2/q3 expectations.
