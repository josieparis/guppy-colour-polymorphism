# Guppy Colour Polymorphism

### Code, scripts and additional data pertaining to "A large and diverse autosomal haplotype is associated with sex-linked colour polymorphism in the guppy"

### These data relate to the Nature Communications Publication which is available here:

https://www.nature.com/articles/s41467-022-28895-4 

Please contact me at parisjosephine<@>gmail.com with any queries

## Raw data
DNA sequencing data are available at the European Nucleotide Archive (ENA) under the Study Accession PRJEB36506: Pool-seq Iso-Y data (SAMEA6512722-SAMEA6512725); whole-genome sequencing data for Paria (SAMEA8750557-SAMEA8750565); long-read pacbio data for Iso-Y6 (SAMEA8795870-SAMEA8795872). Whole-genome sequencing data for Upper Marianne individuals are available from the ENA under the Study Accession PRJEB10680 (SAMEA3649957-SAMEA3649973).

Please refer to the file `sample_metadata.txt` for more information on read names and sample accessions

## SNP calling
SNP calling was performed using the pipeline described here: https://github.com/josieparis/gatk-snp-calling

## Colour phenotype analysis
Colormesh can be found here: https://github.com/J0vid/Colormesh

Publication with example data is published in https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7992 

##### NB these analyses and scripts were performed and written by Mitchel Daniel

### R script for performing DAPC and fish image DAPC loading heatmaps:
#### This script also creates the components for Figure 1 and Supplementary Figure 1
01_DAPC_heatmaps.R

### R script for plotting phenotype data PCAs:
#### This script also creates the components for Supplementary Figure 2
02_phenotype_data_PCA_plots.R

## R script for performing phenotype data permutations
#### This script also creates the components for Supplementary Figure 3
03_phenotype_perm.R

## Genomic data analysis
## LG12 coordinate liftover
The LG12 coordinates were lifted over from the genome available in Fraser et al 2020 GBE (https://doi.org/10.1093/gbe/evaa187) to add additional contigs placed by Deborah Charlesworth's genetic maps:

The function for liftover of chr12 coordinates is:
update_chr12_liftover.R

The updates to the files are run in an R script (and bash to sort the vcf files afterwards)

chr12_liftover.R

## Iso-Y Pool-seq data analysis

### Poolfstat to calculate FST and allele frequencies:
Poolfstat is available here: https://cran.r-project.org/web/packages/poolfstat/index.html 

Script to calculate pairwise FST and allele freqs using poolfstat:
poolfstat_FST_AFs.R

### R script for perfoming Z_FST_PCA analysis and figures:
#### This script also creates figure components for Figure 2a & d, and Supplementary Figures 4,6 and 8
Z_FST_PCA_analysis.R

### R script for polarising / folding AFs analysis and figures:
#### This script also creates figure components for Figure 2b and e
polarise_plot_AFs.R

### R script for polarising and plotting density of AFs:
#### This script also creates figure components for Figure 3, and Supplementary Figures 12, 13 & 14
polarised_AF_densities.R

### R script for change point detection analysis:
CPD_analysis.R

### R script for pairwise FST:
#### This script also creates figure components for Supplementary Figures 5 and 7
raw_pairwise_fst.R

### R script for performing Z_pi_PCA analysis and figures:
#### This script also creates figure components for Supplementary Figures 9 & 11
To calculate pi in 10kb windows we use a custom function (written by Jim Whiting)
pool_pi.R

To run the function and create the figures:
Z_pi_PCA_analysis.R

## Natural WGS data analysis
### R script for LD analysis:
##### LDheatmap package can be found here: https://cran.r-project.org/web/packages/LDheatmap/index.html
#### This script creates figure components for Figure 4a & e, and Supplementary Figure 15
LD_heatmap.R

### R script for Lostruct analysis:
##### Lostruct can be found here: https://github.com/petrelharp/local_pca
#### This script creates figure components for Figure 4b & f, and Supplementary Figure 21
lostruct_localPCA.R

### R script for perfoming popuation genetics stats using PopGenome:
##### PoopGenome can be found here: https://cran.r-project.org/web/packages/PopGenome/index.html
#### This script creates figure components for Figure 4c,g,d & h
popgenome_nat_data.R

### Scripts for creating inter-chrom LD can be found here:
https://github.com/josieparis/interchromLD

### Scripts for creating a genotype plot and haplotype plot as per Figure 5 can be found here:
https://github.com/josieparis/genotype_plot

### R Script for performing GO enrichment and KEGG mapping:
GO_KEGG_analysis_LG1.R

## Interchrom LD can be found in the file:
interchromLD

