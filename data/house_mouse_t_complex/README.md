# House mouse t-complex data

Place the upstream Kelemen & Vicoso data archive locally when running the audit:

```text
IST-2017-78-v1+1_Data.zip
```

The archive is intentionally not committed to git. The Stage 0 script accepts the archive path with `--archive` and inspects nested ZIP contents directly, including:

```text
Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/
  1-Trees_for_all_5kb_windows/
    ML_trees.zip
```

Processed outputs written here are lightweight inventories, mapping tables,
frozen gene-tree files, Stage 2 thinning/sampling inputs, and manifests. The
`raw/` directory is for local user-supplied source material only. Do not add
the upstream archive or nested binary ZIPs to git.
