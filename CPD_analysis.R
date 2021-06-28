# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

#install.packages("changepoint")

lib<-as.vector(c("changepoint","data.table","devtools","rlang","tidyr","ggplot2"))
lapply(lib,library,character.only=TRUE)

## working directory
setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/change-point/")

### Read in the chr12 freqs
### from here:
chr12 <-  data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr12_AFs_final.tsv", header=TRUE))

## Polarise to IF9:
## remove pos
tmp <- chr12 %>% select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-tmp$IF9_AF

for(i in colnames(tmp)){
  tmp[IF9_freqs < 0.5,i]<-(1-tmp[IF9_freqs < 0.5,i])
}

## add IF9 polarised back in back in:
chr12_dd <- cbind(chr12,tmp)

colnames(chr12_dd) <- c("chr", "pos", "REF", "ALT", "IF10_AF", "IF6_AF", "IF8_AF", "IF9_AF",
                        "IF10_polar", "IF6_polar", "IF8_polar", "IF9_polar")


## take out polarised
chr12 <- chr12_dd %>% select(pos,IF10_polar,IF9_polar,IF8_polar,IF6_polar)

## Perform CPD
IF6_chr12 <- chr12 %>% select(IF6_polar)
IF6_chr12 <- as.numeric(IF6_chr12$IF6_polar)

IF8_chr12 <- chr12 %>% select(IF8_polar)
IF8_chr12 <- as.numeric(IF8_chr12$IF8_polar)

IF9_chr12 <- chr12 %>% select(IF9_polar)
IF9_chr12 <- as.numeric(IF9_chr12$IF9_polar)

IF10_chr12 <- chr12 %>% select(IF10_polar)
IF10_chr12 <- as.numeric(IF10_chr12$IF10_polar)

# change in mean allele frequencies
IF6mean=cpt.mean(IF6_chr12, method = "BinSeg", penalty = "SIC", Q=10)

# change in mean allele frequencies
IF8mean=cpt.mean(IF8_chr12, method = "BinSeg", penalty = "SIC", Q=10)

# change in mean allele frequencies
IF9mean=cpt.mean(IF9_chr12, method = "BinSeg", penalty = "SIC", Q=10)

# change in mean allele frequencies
IF10mean=cpt.mean(IF10_chr12, method = "BinSeg", penalty = "SIC", Q=10)

IF6mean@cpts
IF8mean@cpts
IF9mean@cpts
IF10mean@cpts
