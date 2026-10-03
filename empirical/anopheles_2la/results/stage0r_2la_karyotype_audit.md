# Stage 0R 2La karyotype audit

Access date: 2026-10-03

## Authoritative route

The intended authoritative source is the MalariaGEN `malariagen_data.Ag3().karyotype("2La", sample_sets=...)` API, paired with `sample_metadata()` for sample alignment. The Ag3 API documentation lists `karyotype` as the route for inferring inversion karyotype from tag SNPs and lists sample metadata access through `sample_metadata()`.

## Stage 0R access result

No direct sample-level 2La karyotype table was retrieved in this environment.

Two installation attempts were made in the configured project interpreter `/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python`:

- `malariagen-data` latest candidate `15.10.0`: failed before installation because `numba` and `llvmlite` attempted source builds and local SDK/LLVM linking failed.
- `malariagen-data==12.0.1`: failed before installation because `llvmlite`, `numba`, and `biotite` attempted source builds and local SDK/LLVM linking failed.

The repository also contains an earlier Stage 1A attempt recording `malariagen-data==12.0.1`; that run reached package import but was blocked by anonymous GCS access to `gs://vo_agam_release/v3-config.json`. No sample-level karyotype calls were produced there either.

## Representation and normalization

The planned normalized arrangement classes are:

- `2La/2La`
- `2La/2L+a`
- `2L+a/2L+a`
- `unknown`

Stage 0R does not apply this normalization to any sample because the raw source representation was not retrieved. Aggregate literature-fixed states for species are not substituted for direct per-sample karyotypes in `stage0r_sample_counts.tsv`.

## Availability and limitations

- Output categories: not verified from a live Stage-0R API return.
- Confidence/quality fields: not verified from a live Stage-0R API return.
- Missingness: all focal samples are treated as `unknown` for authoritative sample-level counts until the API call succeeds.
- Species support: not verified from a live Stage-0R API return; the API is documented for Ag, but the phasing and site-filter documentation is centered on `gambiae`, `coluzzii`, and `arabiensis`.
- Known species limitations: rare species such as `melas`, `merus`, and `quadriannulatus` are represented in `fontaine-2015-rebuild` aggregate counts, but Stage 0R did not verify whether tag-SNP karyotyping is calibrated or emitted for them.

## Per-sample storage

No per-sample 2La karyotype rows were recovered. The companion table `stage0r_2la_karyotype_samples.tsv` therefore stores only aggregate `not_recovered` rows and does not silently substitute inferred states.
