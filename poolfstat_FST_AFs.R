# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat")

lib<-c("poolfstat","ggplot2", "cowplot", "readr", "tidyr", "stringr", "data.table")
lapply(lib,library,character.only=T)


pool_data <- vcf2pooldata(vcf.file = "../vcf_files/poolseq_variant_filtered_update_STAR.vcf.gz",
                          poolsizes=rep(96,4),
                          poolnames =c("IF10", "IF6", "IF8", "IF9"))
min.cov.per.pool = 30,
max.cov.per.pool = 500,
min.rc = 10)


## check it's been read properly
is.pooldata(pool_data)

## Calculate pairwise FST
pairwise_FST <- compute.pairwiseFST(pool_data, output.snp.values = TRUE, method = "Anova")

## Write:
write.table(pairwise_FST, "./outputs/poolseq_lifted_pairwise_fst.tsv", sep ="\t", quote=F, row.names = F)

## Allele frequencies 

# take info we want from the pool_data object
dd <- cbind(pool_data@snp.info, pool_data@refallele.readcount, pool_data@readcoverage)

## Add column titles
colnames(dd) <- c("chr","pos","REF","ALT",
                  "refallelereadcount_IF10",
                  "refallelereadcount_IF6",
                  "refallelereadcount_IF8",
                  "refallelereadcount_IF9",
                  "readcoverage_IF10",
                  "readcoverage_IF6",
                  "readcoverage_IF8",
                  "readcoverage_IF9")


write.table(dd, "./outputs/whole_genome_raw_read_counts.tsv", sep ="\t", quote=F, row.names = F)

dd <- data.frame(fread("./outputs/whole_genome_raw_read_counts.tsv", header=TRUE))

## calculate ref allele frequencies
AF_dd <- dd %>% mutate(IF10_AF =  refallelereadcount_IF10 / readcoverage_IF10) %>%
  mutate(IF6_AF =  refallelereadcount_IF6 / readcoverage_IF6) %>%
  mutate(IF8_AF =  refallelereadcount_IF8 / readcoverage_IF8) %>%
  mutate(IF9_AF =  refallelereadcount_IF9 / readcoverage_IF9) %>%
  select(chr, pos, REF, ALT, IF10_AF, IF6_AF, IF8_AF, IF9_AF)

write.table(AF_dd, "./outputs/whole_genome_AFs_final.tsv", sep ="\t", quote=F, row.names = F)
