rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/")

# Get libs
lib<-as.vector(c("vcfR","pbapply","PopGenome","cowplot","gridExtra","viridis","grid","data.table","ggplot2","parallel","dplyr","tidyr"))
lapply(lib,library,character.only=TRUE)

# Get VCF (you need to put the full path here otherwise it freaks out)
vcf_path <- "/Users/jrp228/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/chr12.nat_data_31893_updateSTAR.DEFAULT.vcf.gz"

vcf_path <- "/Users/jrp228/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/chr1.nat_data_39606.DEFAULT.vcf.gz"

# Get chrs
vcf<-read.vcfR(vcf_path)
chrs<-unique(vcf@fix[,1])
fai<-read.table("~/Dropbox/Sussex_Guppies/genomes/STAR.chromosomes.release.fasta.fai", header=F)
fai<-fai[fai$V1 %in% chrs,]

# Get males and females
inds<-colnames(vcf@gt)[2:ncol(vcf@gt)]
females<-c("UM10","UM11","UM12","UM13","UM14","UM15", "UM16", "UM17", "UM18", "UM8", "UM9", inds[grep("F",inds)])
males<-inds[!(inds %in% females)]

# Run popgenome over them in different window sizes
wind_sizes=c(1000)
wind_sizes=c(10000)

i<-1
# Set window size
wind_size<-wind_sizes[i]

# Run over chrs for comparisons
all_chr<-data.frame(rbindlist(pblapply(c(1:length(chrs)),function(chr){
  #all_chr<-data.frame(rbindlist(lapply(141:length(chrs),function(chr){
  
  # chr name  
  tid = chrs[chr]
  print(tid)
  # chr length
  topos = fai[fai$V1 == tid,]$V2
  
  # Read in the VCF
  snp <- readVCF(file = vcf_path, tid=tid, frompos = 1, topos = topos, numcols=1000000, include.unknown=TRUE)
  
  # Set up males and females
  snp<-set.populations(snp,do.call("list",list(males,females)),diploid=TRUE)
  snp<-set.outgroup(snp,FALSE)
  
  #split into windows
  win_SNP<-sliding.window.transform(snp,width=wind_size,jump=wind_size,type=2)
  
  # Do stats
  win_SNP <-F_ST.stats(win_SNP,mode="nucleotide")
  win_SNP <-neutrality.stats(win_SNP,FAST=FALSE)
  win_SNP <-diversity.stats.between(win_SNP,nucleotide.mode=TRUE)
  
  # Get centre of the window
  genome.pos_1K <- sapply(win_SNP@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val   <- mean(as.numeric(split))
    return(val)
  })
  
  #output results matrix
  PG_out<-data.frame(chr = tid,
                     FST=win_SNP@nucleotide.F_ST,
                     n.seg.M=win_SNP@n.segregating.sites[,1],
                     n.seg.F=win_SNP@n.segregating.sites[,2],
                     dxy=(win_SNP@nuc.diversity.between/wind_size), 
                     pi.M=(win_SNP@nuc.diversity.within/wind_size)[,1],
                     pi.F=(win_SNP@nuc.diversity.within/wind_size)[,2],
                     Taj_all=win_SNP@Tajima.D,
                     Taj.M=win_SNP@Tajima.D[,1],
                     Taj.F=win_SNP@Tajima.D[,2],
                     Watt.M=win_SNP@theta_Watterson[,1],
                     start=genome.pos_1K -(wind_size/2)-1.5,
                     end=genome.pos_1K +(wind_size/2)-1.5,
                     pop="Paria/Marianne")
  colnames(PG_out)[c(2,5)]<-c("FST","DXY")
  #PG_out<-na.omit(PG_out)
  PG_out[PG_out$FST < 0,"FST"]<-0
  rownames(PG_out)<-NULL
  
  return(na.omit(PG_out))
})))

# Calc resid Dxy
all_chr$Da<-all_chr$DXY-all_chr$pi.F


# Write to output
write.table(all_chr,
            "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/popgen_stats/chr12_nat_data_popgen_1kb_STAR_lift.txt",
            row.names = F,quote=F,sep="\t")