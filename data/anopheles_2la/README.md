Anopheles 2La analysis data.

- `raw/`: immutable external source metadata or downloaded source inputs. Large
  genomic files should remain external unless explicitly approved.
- `metadata/`: region, source, sample, checksum, and provenance manifests.
- `intermediate/`: generated stage handoff products grouped by stage.
- `processed/`: processed Stage-0 tables used by downstream stages.

Stage-0 currently stores source and region manifests plus a schema-complete
sample manifest. Sample-level MalariaGEN/Fontaine metadata and karyotype calls
still need retrieval before strict quartet candidates can be populated.
