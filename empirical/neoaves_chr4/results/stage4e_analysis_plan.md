# Stage 4E analysis plan - ASTER/ASTRAL4 full-tree treatment analysis

Stage 4E is a new ASTER/ASTRAL4 stage. Historical Stage 4D ASTRAL-MP
reproduction attempts remain preserved as provenance and are not reinterpreted
as successful gates. ASTRAL-MP reproduction was discontinued as a prerequisite
for the new analysis. New full-tree treatments use ASTER/ASTRAL4 consistently.

The prepared full-363 treatments are T0, T_PNAS, T_STRUCT_PRIMARY,
T_STRUCT_STRINGENT, T_STRUCT_INCLUSIVE, T_BLOCK_250K, T_BLOCK_500K,
T_BLOCK_1MB, and NONCHR4. Primary execution on Bridges-2 should run T0,
T_PNAS, T_STRUCT_PRIMARY, T_BLOCK_500K, and NONCHR4 first.

Random-matched controls are predeclared as 100 replicate locus-ID manifests for
PNAS, STRUCT_PRIMARY, BLOCK250, BLOCK500, and BLOCK1M. They retain all non-chr4
loci and sample the same number of chr4 loci as the corresponding deterministic
treatment using deterministic seeds.

Jarvis-48 inputs are generated from the same ID-preserving Stage 4E source.
Trees with fewer than four retained taxa are omitted and counted.

Predeclared interpretation categories, to be assigned only after cluster output
validation and random-matched comparison, are:

- FULL-TOPOLOGY CORRECTION
- FOCAL-BRANCH CORRECTION
- SUPPORT/BRANCH-LENGTH EFFECT
- GENERIC-THINNING EFFECT
- NO FULL-TREE EFFECT

ASTER T0 versus the published 63K tree is a historical sanity comparison, not a
gate. ASTER T0 remains the baseline for Stage 4E treatment comparisons even if
RF(ASTRAL4_T0, published_63K) is greater than zero.

No Stage 4E full-tree inference is complete until Bridges-2 outputs are
returned, collected, and validated.
