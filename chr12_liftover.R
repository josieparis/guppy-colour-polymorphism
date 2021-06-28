### chromosome 12 liftover
# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

source("~/Dropbox/Sussex_Guppies/Analyses/update_chr12_liftover.R")

lib<-c("data.table","vcfR","tidyverse","parallel","ggplot2")
lapply(lib,library,character.only=T)

############################################################################################
#### poolseq data ######
############################################################################################

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/vcf_files/")
vcf_file <- "poolseq_variant_filtered.vcf.gz"

vcf<-read.vcfR(vcf_file)
meta<-data.frame(vcf@fix)

# Run the liftover function for the new chr12 STAR positions:
meta[meta$CHROM == "chr12","POS"]<-update_STAR(as.integer(meta[meta$CHROM == "chr12","POS"]))
# Replace metadata
vcf@fix<-as.matrix(meta)

# Write new chr12 VCF
write.vcf(file = "poolseq_variant_filtered_update_STAR.vcf.gz", vcf)

## need to sort the vcf file afterwards using:
# bcftools sort poolseq_variant_filtered_update_STAR.vcf.gz > poolseq_variant_filtered_update_STAR.sorted.vcf.gz 
# Then rename to poolseq_variant_filtered_update_STAR.vcf.gz


############################################################################################
#### natural data ######
############################################################################################

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data")
vcf_file <- "chr12.nat_data.recode.vcf.gz"

vcf<-read.vcfR(vcf_file)
meta<-data.frame(vcf@fix)

# Run the liftover function for the new chr12 STAR positions:
meta[meta$CHROM == "chr12","POS"]<-update_STAR(as.integer(meta[meta$CHROM == "chr12","POS"]))
# Replace metadata
vcf@fix<-as.matrix(meta)

# Write new chr12 VCF
write.vcf(file = "chr12_update_STAR_natural.vcf.gz", vcf)

## Sort vcf file using code above

#########################################################################################################################

## liftover for annotation file
#########################################################################################################################

setwd("/Users/jrp228/Dropbox/Sussex_Guppies/Analyses/pool-seq/chr12_interest_examination")
gff_file <- read_tsv("chr12_genes_OHR95_validated_gene_annotations.gff3", col_names = TRUE)

gff <- as.data.frame(gff_file)

gff$start <- as.character(gff$start)
gff$end <- as.character(gff$end)

# Run the liftover function for the new chr12 STAR positions:
gff[gff$seqID == "chr12","start"]<-update_STAR(as.integer(gff[gff$seqID == "chr12","start"]))
gff[gff$seqID == "chr12","end"]<-update_STAR(as.integer(gff[gff$seqID == "chr12","end"]))


# Find rows that have been inverted
inverted <- which(gff$start > gff$end)
# Save the start pos somewhere
new_start <- gff$start[inverted]
# Now replace the end and start
gff$start[inverted] <- gff$end[inverted]
gff$end[inverted] <- new_start
# And finally fix the strands
inverted_strands <- gff$strand[inverted]
new_strands <- sapply(inverted_strands,function(x){ifelse(x == "+",return("-"),return("+"))})
gff$strand[inverted] <- new_strands

## now sort based on start position:
gff$start <- as.numeric(gff$start)
gff$end <- as.numeric(gff$end)

gff <- gff[order(gff$start),]

## write the new gff:
write.table(gff, "chr12_genes_OHR95_validated_gene_annotations_STAR_liftover.gff3", sep ="\t", quote=F, row.names = F)