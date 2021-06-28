rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/linkage/")
lib<-c("ggplot2","data.table","viridis","parallel", "LDheatmap", "reshape2", "cowplot", "grid")
lapply(lib,library,character.only=T)

##############################
## LG1
##############################

## read in the LD matrix from plink for LG1:
chr1<-read.table("./inputs/chr1.nat_data.POLY.maf0.1.5kTHIN.SQUARE.ld",header = F)

# read in a sites table for the matrix
sites<-read.table("inputs/chr1.nat_data.POLY.maf0.1.5kTHIN_sites.tsv",header = F)

## add the positions to the LD matrix (colnames and rownames)
colnames(chr1)<-sites$V2
rownames(chr1)<-sites$V2

## remove NAs
chr1[chr1 =='NaN']<- NA
## remove Na columns and rows
chr1 <-chr1[rowSums(is.na(chr1)) != ncol(chr1), ]
chr1 <-chr1[,colSums(is.na(chr1)) != nrow(chr1)]

## get a list of the snps
snps<-colnames(chr1)

## list of SNPs every 100 SNPs
snps2<-snps[seq(1,length(snps),100)]

## Keep it as a matrix format:
chr1_ld<-as.matrix(chr1)


chr1_map<-LDheatmap(chr1_ld,add.map = F,color=plasma(100), add.key = TRUE, flip=TRUE,
                    title = "LG1")

png("./figs/LG1_LD.png",res=300, width=600, height=600)
grid.newpage()
grid.draw(chr1_map$LDheatmapGrob)
dev.off()

##############################
## LG12
##############################

## read in the LD matrix from plink
chr12<-read.table("inputs/chr12.nat_data_STAR_update.POLY.maf.5kTHIN.SQUARE.ld",header = F)

# read in a sites table for the matrix
sites<-read.table("./inputs/chr12.nat_data_STAR_update.POLY.maf.5kTHIN.sites.tsv",header = F)

## add sites to the matrix
colnames(chr12)<-sites$V1
rownames(chr12)<-sites$V1

## remove NAs
chr12[chr12 =='NaN']<- NA

## remove Na columns and rows
chr12<-chr12[rowSums(is.na(chr12)) != ncol(chr12), ]
chr12<-chr12[,colSums(is.na(chr12)) != nrow(chr12)]

## get a list of the snps every 100 SNPs
snps<-colnames(chr12)
snps2<-snps[seq(1,length(snps),100)]

chr12_ld<-as.matrix(chr12)

## Plot
chr12_map<-LDheatmap(chr12_ld,add.map = F,color=plasma(100), SNP.name = snps2, add.key=FALSE,title = "Chromosome 12")

png("./figs/LG12_LD.png",res=300, width=600, height=600)
grid.newpage()
grid.draw(chr12_map$LDheatmapGrob)
dev.off()


##############################
## LG1 Males only
##############################

## read in male only linkage
chr1_males <-read.table("./inputs/chr1.nat_data.POLY.maf0.1.5kTHIN.males.SQUARE.ld",header = F)
# read in a sites table for the matrix
males_sites<-read.table("inputs/chr1.nat_data.POLY.maf0.1.5kTHIN.males.sites.tsv",header = F)

## add the positions to the LD matrix (colnames and rownames)
colnames(chr1_males)<-males_sites$V2
rownames(chr1_males)<-males_sites$V2

## remove NAs
chr1_males[chr1_males =='NaN']<- NA
## remove Na columns and rows
chr1_males <-chr1_males[rowSums(is.na(chr1_males)) != ncol(chr1_males), ]
chr1_males <-chr1_males[,colSums(is.na(chr1_males)) != nrow(chr1_males)]

## get a list of the snps
snps<-colnames(chr1_males)

## list of SNPs every 10
snps2<-snps[seq(1,length(snps),10)]

## Keep it as a matrix format:
chr1_males_ld<-as.matrix(chr1_males)

chr1_males_map<-LDheatmap(chr1_males_ld,add.map = T,color=plasma(10), add.key = FALSE, SNP.name = snps2,
                    title = "LG1 Males")

png("./figs/LG1_males_only.png",res=300, width=600, height=600)
grid.newpage()
grid.draw(high_map$LDheatmapGrob)
dev.off()



##############################
## LG1 Females only
##############################

## read in male only linkage
chr1_females <-read.table("./inputs/chr1.nat_data.POLY.maf0.1.5kTHIN.females.SQUARE.ld",header = F)
# read in a sites table for the matrix
females_sites<-read.table("inputs/chr1.nat_data.POLY.maf0.1.5kTHIN.females.sites.tsv",header = F)

## add the positions to the LD matrix (colnames and rownames)
colnames(chr1_females)<-females_sites$V2
rownames(chr1_females)<-females_sites$V2

## remove NAs
chr1_females[chr1_females =='NaN']<- NA
## remove Na columns and rows
chr1_females <-chr1_females[rowSums(is.na(chr1_females)) != ncol(chr1_females), ]
chr1_females <-chr1_females[,colSums(is.na(chr1_females)) != nrow(chr1_females)]

## get a list of the snps
snps<-colnames(chr1_females)

## list of SNPs every 10
snps2<-snps[seq(1,length(snps),10)]

## Keep it as a matrix format:
chr1_females_ld<-as.matrix(chr1_females)

chr1_females_map<-LDheatmap(chr1_females_ld,add.map = T,color=plasma(10), add.key = FALSE, SNP.name = snps2,
                          title = "LG1 Males")

png("./figs/LG1_females_only.png",res=300, width=600, height=600)
grid.newpage()
grid.draw(high_map$LDheatmapGrob)
dev.off()


#############################################################################################################
## For the high linkage region ###
#############################################################################################################

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/linkage/")

high<-read.table("inputs/chr1.10-16.highLD_region.SQUARE.ld",header = F)

# read in a sites table for the matrix
sites<-read.table("inputs/chr1.10-16.highLD_region.sites.tsv",header = F)

## add sites to the matrix
colnames(high)<-sites$V1
rownames(high)<-sites$V1

## remove NAs
high[high =='NaN']<- NA

## remove Na columns and rows
high<-high[rowSums(is.na(high)) != ncol(high), ]
high<-high[,colSums(is.na(high)) != nrow(high)]

## get a list of the snps
snps<-colnames(high)
snps2<-snps[seq(1,length(snps),100)]

high_ld<-as.matrix(high)
high_map<-LDheatmap(high_ld,add.map = F,color=plasma(100), add.key=FALSE, title=NULL)

png("./figs/highLD_LG1.png",res=300, width=600, height=600)
grid.newpage()
grid.draw(high_map$LDheatmapGrob)
dev.off()
