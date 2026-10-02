# Songbirds PIPs Stage-0 feasibility and provenance audit

This Stage-0 audit freezes source provenance, PIP/sample metadata, and topology-independent candidate-selection criteria for Pegan and Winger (2025). No local trees were inferred, and no candidate was selected by inspecting Stage-0 tree outcomes.

## Dataset size and released data

- Primary article: Pegan and Winger, Genome Biology and Evolution 17(11):evaf205, DOI 10.1093/gbe/evaf205.
- Public sequence data: NCBI SRA BioProjects PRJNA1043688 and PRJNA1130443.
- Supplementary data downloaded locally: `data/songbirds_pips/raw/oup_supplement/`.
- Species examined: 35.
- Individuals in article Table 1: 1660.
- Individuals with PIP genotype rows in Dataset S4: 1411 across 28 PIP-containing species.
- Total PIPs in Dataset S1: 174; LD-supported/high-confidence PIPs: 139; PCA-only/LD-ambiguous PIPs: 35.
- PCA cluster classes: 133 three-cluster PIPs, 34 two-cluster PIPs, 7 four-to-six-cluster PIPs.

## Feasibility answer

Local genealogies can be reconstructed, but not directly from packaged analysis-ready phylogenetic inputs. The cheapest reliable route is to request or recover the authors' intermediate BAM/ANGSD genotype-likelihood/SAF files. Without those, Stage 1 must download SRA reads for selected samples only, regenerate BAMs/genotype likelihoods for candidate regions and matched controls, and then infer local genealogies or distance trees. PCA output alone is label metadata and is not a phylogenetic input.

## Data availability summary

- Available now: PIP coordinates, reference assemblies, LD-support flags, PCA cluster counts, per-individual PIP genotype assignments, sample coordinates/SRA accessions for PIP species, and representative sample lists for shared-PIP author trees.
- Not found as public packaged analysis-ready files: VCF/BCF, phased haplotypes, SNP matrices, BAM/CRAM, alignments, reusable Newick local trees, or whole-genome trees.
- Raw reads are public but terabyte-scale; Stage 0 did not download them.

## Candidate-selection criteria

Candidates were scored using only structural/data-quality criteria: PCA support, LD support, genotype-class representation, sample count, clear coordinates, COMP/control feasibility, SRA availability, and author-reported shared status. No topology-derived statistic from local trees was used. The one cross-species candidate uses author-reported overlap/shared-PIP metadata; Stage 0 did not inspect tree topology plots.

Anti-circularity audit: Stage 0 read the primary article's prose summary of the authors' neighbor-joining/shared-PIP analysis while inventorying data availability and cross-species claims. No Dataset S3 tree plots, Newick trees, distance matrices, or newly inferred local trees were inspected. Candidate selection used PIP coordinates, PCA/LD/genotype metadata, sample availability, and author-declared shared-region membership only; it did not use whether a candidate's tree looked unusual or favorable for MSRC.

## Selected initial PIPs

- `Empidonax_flaviventris|CM020533.1` (A_strong_within_species): CM020533.1:39775054-99074287, 59.3 Mb, n=56, AA/AB/BB=28/22/6, LD=yes, cross-species=no.
- `Empidonax_minimus|CM020533.1` (B_independent_within_species): CM020533.1:98876212-138925002, 40 Mb, n=43, AA/AB/BB=11/22/10, LD=yes, cross-species=no.
- `Geothlypis_philadelphia|CM019906.1` (independent_warbler_within_species): CM019906.1:5075009-24574663, 19.5 Mb, n=50, AA/AB/BB=6/28/16, LD=yes, cross-species=no.
- `Zonotrichia_albicollis|CM042600.1` (D_small_LD_supported_control): CM042600.1:4625010-5024956, 0.4 Mb, n=66, AA/AB/BB=44/19/3, LD=yes, cross-species=no.
- `shared_CM027530.1` (C_cross_species_warbler_shared_region): 11 warbler species, DatasetS5 representative n=22, author shared-PIP region on CM027530.1; exact shared-overlap coordinates should be reconstructed from constituent PIPs before Stage 1.

## Cross-species/shared PIPs

Dataset S5 contains 26 shared-PIP tree groups and 256 representative rows. The article reports 24 chromosomes with PIPs in multiple species, 26 overlapping PIP cases, and 12 author-interpreted trans-specific cases. The processed cross-species inventory records all DatasetS5 shared groups without inferring introgression or ancient balanced polymorphism beyond the authors' framing.

## Matched control strategy

For within-species candidates, use the authors' COMP intervals where available, then add matched same-chromosome windows outside PIPs with comparable length/SNP density/missingness after Stage-1 genotyping. For the cross-species candidate, reconstruct the authors' shared-PIP and shared-COMP design, then add matched non-PIP regions on the same reference chromosome only after confirming all species use the same reference coordinate frame.

## Coordinate and assembly cautions

Coordinates are species/reference-specific. Dataset S1 gives the reference genome assembly and reference species for each PIP. Do not merge regions across species unless the species were mapped to the same reference genome and overlapping coordinates are explicitly reconstructed. No liftover resources were identified in Stage 0.

## Computational estimate

A selected within-species PIP would require tens of individuals and one focal interval plus matched controls; local execution may be reasonable if starting from BAM/CRAM or small regional read extractions, but full SRA-to-BAM processing is HPC-scale. The full public read archive is terabyte-scale. The cross-species candidate spans many species and should be run on HPC unless author intermediates are obtained.

## Blockers

- No public analysis-ready VCF/BAM/CRAM/phased haplotype/local-tree files were identified.
- Dataset S4 lacks sex and explicit population labels; only coordinates and museums are provided.
- Seven no-PIP species are not included in Dataset S4; their full metadata are in the related Pegan et al. 2025 materials.
- Shared-PIP exact overlap coordinates for cross-species tests need reconstruction from constituent PIPs and/or author metadata before analysis.

## Stage 1 recommendation

Do not launch genome-wide processing. First contact/check for author-provided BAM, ANGSD genotype-likelihood, SAF, or regional distance-matrix intermediates. If unavailable, build a small SRA-reprocessing pilot for the selected PIPs only, starting with one within-species candidate and one matched COMP/control region. Freeze the exact sample sets and control intervals before any local-tree inference.
