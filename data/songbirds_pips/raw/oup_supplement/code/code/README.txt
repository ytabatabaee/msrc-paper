Organization of scripts and data in this repository

We include here csv versions of three of our supplementary datasets: Dataset S1, Dataset S4, and Dataset S5.

The computational and analytical steps used to prepare data for this manuscript can be 
divided into 4 parts.

The first part comprises bioinformatic processing to filter raw data, align it to 
reference genomes, and call SNPs in a genotype likelihood framework. These steps were 
conducted during analyses for Pegan et al. 2025, doi 10.1038/s41559-025-02699-3
Relevant scripts are shared in a public figshare repository associated with that manuscript with 
doi 10.6084/m9.figshare.27284553
Within the above figshare repository, the following scripts contain the code that was used to
pre-process data also used in the present manuscript:
1.Sample_Processing
2.Realign_around_indels
3.Print_SNP_GLs_with_ANGSD
4.Filter_SNPs_with_ngsParalog

The next three sections correspond to code provided here, 
described below:

1) Bioinformatic code for detecting and analyzing PIPs
Code in this section requires analysis of very large sequence datasets which are not 
included here. Sequence data are publically archived on the NCBI SRA with accession 
numbers listed in Dataset S4. Prior to use in this section, sequence data from SRA must 
be processed using code in the four files listed above from the Pegan et al. 2025 
figshare doi listed above. 
The code in this repository is not designed to be ready-to-run, but rather to clearly 
indicate which steps were taken and which arguments were used with each command. 
Each file contains detailed notes about how each step was run,how files could be 
organized, and what a user would need to do to replicate these analyses.
In the folder 1_Bioinformatic_code, scripts labeled A1-A4 in this folder correspond to 
the steps for detecting PIPs described in Fig. S1. Scripts labeled B1-B3 contain code for 
genotyping individuals based on PCA (B1), evaluating heterozygosity (B2), and preparing 
genetic distances for gene trees (B3). Scripts that do not start with one of these 
prefixes are R or python scripts that are called from within the main scripts.

2) R script for generating summary statistics
This script generates summary statistics about PIPs as described in the Result section. 
This script is ready to be run using information provided in the supplementary 
information of the present manuscript.

3) Bioinformatic code for using centroAnno to annotate HORs on reference genomes.
Code in this section requires analysis of very large sequence datasets which are not 
included here. Specifically, the necessary input is published reference genomes listed in 
Dataset S1. The code in this file is not designed to be ready-to-run, but rather to 
clearly indicate which steps were taken and which arguments were used with each command.
