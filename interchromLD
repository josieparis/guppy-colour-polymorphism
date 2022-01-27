# interchromLD

## Methods and scripts for creating inputs and plotting inter-chromosomal linkage

#### Make the interchr linkage LD pairwise info from plink:

##### Variables
`VCF=two_chr`

##### Add SNP_ID
`bcftools annotate --set-id +'%CHROM\_%POS' $VCF.vcf > ${VCF}_IDs.vcf`

##### Keep only poly sites
`bcftools view -g het ${VCF}_IDs.vcf -o ${VCF}_POLY.vcf`

##### maf filter poly sites
`vcftools --vcf ${VCF}_POLY.vcf --maf 0.1 --recode --out ${VCF}_POLY.maf0.1`

##### Thin
`vcftools --vcf ${VCF}_POLY.maf0.1.recode.vcf --thin 5000 --recode --out ${VCF}_POLY.maf0.1.5kTHIN`

##### Run plink with --inter-chr flag
`plink --vcf ${VCF}_POLY.maf0.1.5kTHIN.recode.vcf --r2 inter-chr --allow-extra-chr --double-id --out ${VCF}_POLY.maf0.1.5kTHIN.inter-chrLD`

## Plot in R

```
######################################################################
############ Script to plot inter chromosomal linkage ################
######################################################################

rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

# Get libs
lib<-as.vector(c("cowplot","gridExtra","viridis","scales", "grid","data.table","ggplot2","scales","dplyr","tidyr", "stringr", "ggpubr"))
lapply(lib,library,character.only=TRUE)

setwd("<path>")

## read input from plink:
dd <- read.table("./inputs/chr1_chr12_interLD.ld", header=TRUE)

## filter each chromosome:
chr1 <- dd %>% filter(CHR_A == "1" & CHR_B == "1")
chr1$SNP_A <- reorder(chr1$SNP_A, chr1$BP_A)
chr1$SNP_B <- reorder(chr1$SNP_B, chr1$BP_B)
## complete the matrix
chr1_comp <- complete(chr1,SNP_A,SNP_B)

chr12 <- dd %>% filter(CHR_A == "12" & CHR_B == "12")
chr12$SNP_A <- reorder(chr12$SNP_A, chr12$BP_A)
chr12$SNP_B <- reorder(chr12$SNP_B, chr12$BP_B)
chr12_comp <- complete(chr12,SNP_A,SNP_B)

# and linkage between both
chr1_chr12 <- dd %>% filter(CHR_A == "1" & CHR_B == "12")
chr1_chr12$SNP_A <- reorder(chr1_chr12$SNP_A, chr1_chr12$BP_A)
chr1_chr12$SNP_B <- reorder(chr1_chr12$SNP_B, chr1_chr12$BP_B)
chr1_chr12_comp <- complete(chr1_chr12,SNP_A,SNP_B)


## first chr
pchr1 <- ggplot(chr1_comp)+
  geom_tile(aes(SNP_A,SNP_B,fill=R2))+
  scale_fill_viridis(option="plasma", begin = 0.1, end = 1,
                     direction=-1, na.value = "#FCFDBF",
                     oob = scales::squish)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# second chr
pchr12 <- ggplot(chr12_comp)+
  geom_tile(aes(SNP_A,SNP_B,fill=R2))+
  scale_fill_viridis(option="plasma", begin = 0.1, end = 1,
                     direction=-1, na.value = "#FCFDBF",
                     oob = scales::squish)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# both chr
pchr112 <- ggplot(chr1_chr12_comp)+
  geom_tile(aes(SNP_A,SNP_B,fill=R2))+
  scale_fill_viridis(option="plasma", 
                     begin = 0.1, end = 1,
                     direction=-1, na.value = "#FCFDBF",
                     oob = scales::squish)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# put together
pchr <- ggarrange(pchr112,pchr12,pchr1,ncol=2,nrow=2)

ggsave("<path>/Figures/LG1_LG12_linkage.png", pchr, width = 40, height = 40, units = "cm")
```

