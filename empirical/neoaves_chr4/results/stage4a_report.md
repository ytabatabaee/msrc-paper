# Neoaves chromosome 4 Stage 4A: structural topology support

## Scope

This analysis asks whether frozen rearrangement-associated regions disproportionately contribute support for competing quartet relationships. It measures spatial association and does not identify a causal rearrangement effect or demonstrate bias in a species-tree inference.

## Frozen inputs and topology reference

- Ordered tracks: `data/neoaves_chr4/processed/locus_table_chr4.tsv` (Columbea, N61, N62; q1/q2/q3 and dominant state).
- Frozen structural inputs: `denovo_breakpoint_sets.tsv`, the three `stage3b_*_structural_segments.tsv` partitions, and `stage3b_v2_structural_events.tsv`.
- Stage 3A provenance: `stage3a_final_manifest.json` and its breakpoint enrichment outputs.
- Stage 3B provenance: `stage3b_v2_final_manifest.json`, summary, event states, and structural-only prediction manifest.
- The reference species topology is q1 for every focal clade: q1=(C1,C2)|(S,O). This follows the frozen quartet role mapping. It is distinct from the chromosome-wide empirically dominant topology.

## Structural masks

The mask is the union of ±250,000 bp neighborhoods around every breakpoint in a frozen set, clipped to chr4:1-91,310,470. The 250 kb radius is the pre-existing Stage 3B primary flank window. No genealogy data enter mask construction. Stringent, primary, and inclusive are sensitivity sets; background is the complement.

- stringent: 15 merged intervals, 7,846,083 bp (8.593% of chr4).
- primary: 23 merged intervals, 12,521,582 bp (13.713% of chr4).
- inclusive: 30 merged intervals, 19,247,660 bp (21.079% of chr4).

## Results

### Structural category frequencies

**stringent**

| clade | associated q1/q2/q3 | background q1/q2/q3 | associated dominant | background dominant |
|---|---:|---:|---|---|
| Columbea | 0.672/0.175/0.154 | 0.472/0.267/0.261 | q1 | q1 |
| N61 | 0.203/0.658/0.140 | 0.295/0.437/0.268 | q2 | q2 |
| N62 | 0.654/0.193/0.153 | 0.438/0.291/0.271 | q1 | q1 |

**primary**

| clade | associated q1/q2/q3 | background q1/q2/q3 | associated dominant | background dominant |
|---|---:|---:|---|---|
| Columbea | 0.640/0.187/0.173 | 0.464/0.270/0.265 | q1 | q1 |
| N61 | 0.212/0.624/0.164 | 0.299/0.428/0.273 | q2 | q2 |
| N62 | 0.604/0.218/0.178 | 0.432/0.293/0.275 | q1 | q1 |

**inclusive**

| clade | associated q1/q2/q3 | background q1/q2/q3 | associated dominant | background dominant |
|---|---:|---:|---|---|
| Columbea | 0.593/0.216/0.191 | 0.460/0.271/0.269 | q1 | q1 |
| N61 | 0.226/0.576/0.198 | 0.304/0.422/0.274 | q2 | q2 |
| N62 | 0.563/0.235/0.202 | 0.427/0.296/0.278 | q1 | q1 |

### Spatial enrichment

One-sided p-values use 10,000 deterministic circular shifts of each complete ordered topology-state track while holding the structural mask fixed. The three tracks are correlated views and are not pooled as independent replicates.

| set | clade | competitor | fold | difference | representation/span ratio | p |
|---|---|---|---:|---:|---:|---:|
| stringent | Columbea | q2 | 0.655 | -0.092 | 0.618 | 0.9671 |
| stringent | N61 | q2 | 1.506 | 0.221 | 1.317 | 0.0031 |
| stringent | N62 | q2 | 0.664 | -0.098 | 0.600 | 0.9480 |
| primary | Columbea | q2 | 0.691 | -0.083 | 0.691 | 0.9307 |
| primary | N61 | q2 | 1.456 | 0.195 | 1.315 | 0.0166 |
| primary | N62 | q2 | 0.744 | -0.075 | 0.717 | 0.8382 |
| inclusive | Columbea | q2 | 0.799 | -0.054 | 0.819 | 0.8763 |
| inclusive | N61 | q2 | 1.365 | 0.154 | 1.246 | 0.0292 |
| inclusive | N62 | q2 | 0.795 | -0.061 | 0.797 | 0.8470 |

### Window versus block normalization

| clade | scheme | q1/q2/q3 | dominant | species-best-alt margin | empirical competing-topology fraction |
|---|---|---:|---|---:|---:|
| Columbea | window-weighted | 0.488/0.260/0.253 | q1 | 0.228 | 0.260 |
| Columbea | block-normalized | 0.352/0.331/0.317 | q1 | 0.021 | 0.331 |
| N61 | window-weighted | 0.288/0.454/0.258 | q2 | -0.166 | 0.454 |
| N61 | block-normalized | 0.346/0.332/0.323 | q1 | 0.014 | 0.332 |
| N62 | window-weighted | 0.454/0.283/0.262 | q1 | 0.171 | 0.283 |
| N62 | block-normalized | 0.318/0.347/0.335 | q2 | -0.030 | 0.347 |

Dominant topology changes under block normalization: N61, N62.

The competing-topology fractions are descriptive empirical analogues of epsilon_window and epsilon_block. They are not literal estimates of theoretical contamination parameters because windows/runs need not satisfy the model's independence or generative assumptions.

### Event heterogeneity

- Columbea: event margins range -0.154 to 1.000; strongest local tilt toward q2 occurs at primary_W250000_EV001 (margin -0.154).
- N61: event margins range -1.000 to 0.121; strongest local tilt toward q2 occurs at primary_W250000_EV008 (margin -1.000).
- N62: event margins range -0.467 to 1.000; strongest local tilt toward q2 occurs at primary_W250000_EV011 (margin -0.467).

No event is treated as a strict 2:2 resolved mechanism unit here. The Stage 3B eligibility result remains zero and unchanged.

## Interpretation and limitations

1. Competing-topology enrichment is track-specific. N61 q2 is enriched in every frozen mask (fold 1.37-1.51; p=0.0031-0.0292), whereas Columbea and N62 q2 are depleted rather than enriched. The correlated tracks therefore do not support a chromosome-wide claim that structural neighborhoods generally enrich competing topologies.
2. Dense windows amplify N61 q2: its empirical competing-topology fraction falls from 0.454 by windows to 0.332 by blocks. For Columbea and N62, q2 instead rises after block normalization (0.260 to 0.331 and 0.283 to 0.347), so dense sampling dilutes that competitor in those tracks.
3. Dominance changes for N61 (q2 to q1) and N62 (q1 to q2), but not Columbea (q1 in both summaries).
4. N61's event pattern is broadly distributed: 13 of 16 event units have a negative q1-minus-q2 margin. The comparable counts are 3 of 16 for Columbea and 4 of 16 for N62, with substantial event-to-event heterogeneity in every track.
5. Direction and inference are stable across stringent, primary, and inclusive masks: N61 is enriched in all three; Columbea and N62 are depleted in all three.

Enrichment is assessed separately for each correlated focal track and across frozen mask sensitivities. Evidence that remains similar after block normalization is less attributable to dense runs of repeated windows; changes in fractions quantify that amplification directly. Event-level ranges show whether the pattern is concentrated or distributed, but overlapping neighborhoods and shared genealogy data prevent treating events or clades as fully independent biological replicates.

These results establish association only. They do not show that rearrangements caused a topology shift, estimate a literal MSRC failure parameter, or establish an effect on an actual species-tree analysis relative to independent evidence.

## Figures

- Figure A: chromosome position, primary structural neighborhoods, and topology state/support.
- Figure B: primary associated versus background topology frequencies.
- Figure C: window-weighted versus block-normalized support.
