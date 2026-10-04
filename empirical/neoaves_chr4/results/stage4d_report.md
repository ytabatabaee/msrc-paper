# Stage 4D checkpoint report — incomplete

No novel Stage-4D treatment has been run. No interpretation category or
full-tree correction claim is assigned. The original paper's chr4/taxon-sampling
effect is prior evidence, not a new MSRC finding.

## Inputs and completed preparation

Existing downloads are reused without downloading raw alignments. Exact URLs,
byte sizes, SHA256 values, and integrity results are recorded in
`stage4d_download_manifest.tsv`. The named and collapsed tree archives contain
63,430 trees. The locus manifest has 63,430 unique IDs and deterministic window
coordinates, including 6,247 anchored chr4 loci and 15 chr4 random-scaffold loci.
NONCHR4 therefore has 57,168 loci. There are 1,431 published-outlier loci and
495/835/1,308 stringent/primary/inclusive structural-mask loci.

The existing exhaustive audit found 31,066 same-row disagreements between the
named and collapsed files. Collapsed IDs must not be assigned by row. An
ID-preserving support<0.95 re-collapse of named trees already exists, but its
full-tree reproduction remains required before novel treatments.

## Historical estimator and current checkpoints

The archived log establishes ASTRAL-MP 5.15.1, unrooted input, 28 CPU threads,
and four Tesla P100 GPUs. It reports total runtime 112,915.644 seconds
(approximately 31.4 hours). The historical source supplies default seed 692
and annotation level 3. The original complete shell command and Java heap
options have not been recovered; source defaults are not evidence of an exact
historical command line.

The final tree printed in that archived log and the supplied `63K.tre` have
unrooted RF=0, normalized RF=0, and zero differing bipartitions. **This validates
archive consistency, not a fresh reproduction.**

Applying the frozen focal role definitions to the supplied published tree gives
Columbea q2, N61 q1, and N62 q1 (each corresponding split has supplied support
label 1). Stage 4C's independent non-chr4 focal summaries instead favor q1,
q1, and q2. Thus the published full-tree topology and the non-chr4 focal
summaries disagree for Columbea and N62. This is an important distinction:
S2024-like and Stage-4C-reference-like are not interchangeable labels.
Neither the frozen definitions nor the Stage-4C results were changed.
Exact group sizes and source checksum are in
`stage4d/published_focal_relationships.json`. This is an extraction from a
published reference, not a new treatment or a new inferred NONCHR4 tree.

The resumed full-363 attempt is recorded under
`stage4d/attempts/full363_resume_01/`. Its provenance includes input and runtime
checksums. It used the historical implementation on CPU, four threads, seed
692, annotation level 3, and an 8 GB heap on a 16 GB host. After 429.57 seconds,
the run was stopped during memory pressure before finishing input loading.
A process sample shows the main thread waiting for G1 heap allocation and
full-GC marking workers. Graceful termination did not respond, so SIGKILL was
used. This was an operator stop, not a reported Java OutOfMemoryError or a
topology mismatch. No species tree was produced. A higher-memory execution
environment is needed to establish full-data feasibility; no exact minimum
memory requirement has been established.

The earlier FULL363 and J48 attempts have empty output trees and do not pass
reproduction. Jarvis-48 all-locus and outlier-excluded input files exist, with
62,945 and 61,524 usable trees respectively; trees with fewer than four sampled
taxa are omitted. Neither Jarvis-48 topology reproduction has passed. The new
`J48_ALL_20260919T121120529946Z` attempt loaded all 62,945 usable trees and
entered search, but was stopped after a process sample showed full-GC
allocation failure at its 4 GB heap limit. The retry
`J48_ALL_20260919T121536445811Z` changes only the heap limit to 8 GB.
Both attempts retain their logs and provenance. Their intermediate search-space
cluster counts differ (1,406,493 versus 1,404,805) despite identical seed and
scientific settings. Thread scheduling is a possible explanation, not a
verified cause; no completed topology is available to assess the consequence.

## Figure 4

The existing reproduction table and Figure A aggregate the authors' per-locus
scores for all 36 conditions. The exact grid is in
`stage4d/figure4_original_grid.tsv`, including expanded removal names and axis
order. Aggregation uses `-(main-alt)`, arithmetic mean, median, and sample
SD/sqrt(n), retaining the nine zero-score loci as the original script did.
These are descriptive locus standard errors, not independence-valid MSRC
uncertainty estimates.

Taxon sets were expanded from the original condition labels and count columns;
this does not independently validate original pruned input trees. De novo score
verification remains pending. The published outlier rule is inclusive on chr4
window starts: 25,030,000–32,670,000; 33,510,000–34,470,000; and
44,130,000–56,810,000.

## Outstanding analyses and claims

The final analysis freeze and its checksum do not yet exist. T_STRUCT,
T_BLOCK500, 250-kb/1-Mb sensitivities, and count-matched random distributions
have not been evaluated. There are no new full-tree RF distances to NONCHR4,
changed-branch results, N61/N62/Columbea outcomes, chromosome controls, or
chromosome jackknife results. Structural-versus-PNAS correction and mean/median
correction cannot be interpreted from mask counts or source summaries alone.

No manuscript-ready positive or negative treatment claim is justified. The
appropriate current statement is: “Full-tree tests of topology-neutral block
normalization and independent structural exclusion remain pending successful
reproduction of the published analyses.”

## Preservation and tests

The authorized Stage-4B typo restoration matches its original frozen checksum.
No scientific Stage 1–4C output was altered. The pre-restoration audit is retained.
The complete Neoaves suite passes 34 tests, including additional checks
that interrupted ASTRAL attempts are not accepted as completed and that
checkpoint inference settings remain fixed, and focal extraction handles
alternative splits, missing groups, polytomies, and dispersed group members.
Passing these software tests does
not imply that scientific reproduction gates have passed.
