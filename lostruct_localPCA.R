##########################################################################################
# Script runs local PCA across a chromosomal PCA to infer windows of shared ancestry

# Version 2 of this script runs off a whole VCF and subsets in-line
#########################################################################################

#install.packages("ggpubr")
#install.packages("data.table")
#devtools::install_github("petrelharp/local_pca/lostruct")
#library(lostruct)

# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

# Set working directory
setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/lostruct")

## Load packages
lib = c("parallel","data.table","tidyr","dplyr","ggplot2","lostruct","ggpubr","vcfR", "gridExtra", "extrafont")
lapply(lib, library, character.only=TRUE)

# Define various groupings, can either be split into males and females, or all 
female_vec<-c("Paria_F12",  "Paria_F16",  "Paria_F2",  "Paria_F2_RNA",  "Paria_F8",  "UM10",  "UM11",  "UM12",  "UM13",  "UM14",  "UM15",  "UM16",  "UM17",  "UM18",  "UM8",  "UM9")
male_vec<-c("Paria_M12",  "Paria_M2",  "Paria_M5",  "Paria_M9",  "UM1",  "UM2",  "UM3",  "UM5",  "UM6",  "UM7")


# Get parameters from commandline to find outputs

# This is the name of the original vcf
vcf<- "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/chr1.nat_data_39606.DEFAULT.vcf.gz"
## Use both 10 and 100bp windows
window_N<-10
chr<-"chr1"
pops<-c("Paria", "UM")
## change output depending on whether or not you're doing 10bp or 100bp windows
output<-"chr1_nat_data_10bp"

# Make the input
system(paste0("bcftools view -r ",chr," ",vcf," > localPCA_tmp.vcf"),wait=TRUE)

# Read in VCF in popGenome for info on BP locations for unphased vcf
vcf_PG<-read.vcfR("localPCA_tmp.vcf")
BP1<-as.integer(vcf_PG@fix[,2])
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]
inds_vec<-colnames(vcf_PG@gt)[2:length(colnames(vcf_PG@gt))]

# Filter for inds we want
#inds_to_keep<-unlist(lapply(pops,function(x){
#  inds2<-inds_vec[inds_vec %in% male_vec]
#return(inds2[grep(x,inds2)])
#}))


# vcf_PG@gt<-vcf_PG@gt[,c(colnames(vcf_PG@gt)[1],inds_to_keep)]
# Get males and female columns
#females<-match(female_vec,inds_vec)[is.na(match(female_vec,inds_vec))==FALSE]
#males<-setdiff(1:length(inds_vec),females)

# Return logical vector for polymorphic sites
invariant_filter<-is.polymorphic(vcf_PG,na.omit=T)

# Read in VCF for localPCA
vcf_in<-read_vcf("localPCA_tmp.vcf")

# And tidy up...
system("rm -f localPCA_tmp.vcf")

# Filter for invariants and individuals
#vcf_in<-vcf_in[invariant_filter,(inds_vec %in% male_vec)]
BP1<-BP1[invariant_filter]
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]

# Read in eigen_windows, win=number of rows of matrix, k = number of eigenvector/value pairs to return, mc.cores compatible but does seem to crash sometimes with mc.cores
eigenstuff <- eigen_windows(vcf_in, win=window_N, k=3)

# Calculate the distances, npc=N of PCs computed, mc.core compatible
windist <- pc_dist(eigenstuff, npc=3,mc.cores=detectCores()-1)

# Remove any windows for which we are getting NA values
to_keep_logical<-rowSums(is.na(windist)) != ncol(windist)
to_keep<-(1:ncol(windist))[to_keep_logical]
windist2<-windist[to_keep,to_keep]

# Principal Coordinates Analysis, k = dimension of space for data to represented in, eig = return eigenvalues
fit2d <- cmdscale(windist2, eig=TRUE, k=3)

eigenvals <- fit2d$eig/sum(fit2d$eig)
## check out sig eigvalues
eigenvals_10 <- eigenvals[eigenvals > 0.01] ## first 3 MDS

# Plot along a chromosome
plot_dd<-data.frame(mds1=fit2d$points[,1],
                    mds2=fit2d$points[,2],
                    mds3=fit2d$points[,3],
                    window_N=1:length(fit2d$points[,1]),
                    window_pos=window_BP[to_keep])


# Calculate cut-offs for outliers by trimming distributions and taking 3*SD of trimmed distribution
trim_mds1<-plot_dd[plot_dd$mds1 > quantile(plot_dd$mds1,probs=0.05) & plot_dd$mds1 < quantile(plot_dd$mds1,probs=0.95),"mds1"]
mds1_cutoffs<-c(mean(trim_mds1)+3*sd(trim_mds1),
                mean(trim_mds1)-3*sd(trim_mds1))

trim_mds2<-plot_dd[plot_dd$mds2 > quantile(plot_dd$mds2,probs=0.05) & plot_dd$mds2 < quantile(plot_dd$mds2,probs=0.95),"mds2"]
mds2_cutoffs<-c(mean(trim_mds2)+3*sd(trim_mds2),
                mean(trim_mds2)-3*sd(trim_mds2))

trim_mds3<-plot_dd[plot_dd$mds3 > quantile(plot_dd$mds3,probs=0.05) & plot_dd$mds3 < quantile(plot_dd$mds3,probs=0.95),"mds3"]
mds3_cutoffs<-c(mean(trim_mds3)+3*sd(trim_mds3),
                mean(trim_mds3)-3*sd(trim_mds3))

##################
# We also want to return the windows at the extremes of the distribution
##################
# Give window IDs
window_ids<-rep(0)
for(j in 1:length(to_keep)){
  i <- to_keep[j]
  window_ids[j]<-paste0(chr,":",((BP1[(i-1)*window_N])+1),"-",BP1[(i)*window_N])
}
######## This was shortening the vector by filtering again on "to_keep"
plot_dd$window_id <- window_ids
######################################################## 

windows_to_keep1<-unique(c(plot_dd[plot_dd$mds1 < mds1_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds1 > mds1_cutoffs[1],"window_id"]))

windows_to_keep2<-unique(c(plot_dd[plot_dd$mds2 < mds2_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds2 > mds2_cutoffs[1],"window_id"]))

windows_to_keep3<-unique(c(plot_dd[plot_dd$mds3 < mds3_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds3 > mds3_cutoffs[1],"window_id"]))

# Reassign to plot_dd
plot_dd$outlier<-rep("N")
plot_dd[plot_dd$window_id %in% windows_to_keep1,"outlier"]<-"mds1"
plot_dd[plot_dd$window_id %in% windows_to_keep2,"outlier"]<-"mds2"
plot_dd[plot_dd$window_id %in% windows_to_keep3,"outlier"]<-"mds3"

# Throw in the colour vals
cols1 <- c("mds1" = "#FCD225", "mds2" = "#C92D59", "mds3" = "#300060", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")
cols2 <- c("mds2" = "#C92D59", "mds1" = "#FCD225", "mds3" = "#300060", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")
cols3 <- c("mds3" = "#300060", "mds2" = "#C92D59", "mds1" = "#FCD225", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")
cols4 <- c("mds4" = "#a52a2a", "mds3" = "#300060", "mds2" = "#C92D59", "mds1" = "#FCD225", "N" = "black", "mds5" = "green", "mds6" = "blue")
cols5 <- c("mds5" = "green", "mds4" = "#a52a2a", "mds3" = "#300060", "mds2" = "#C92D59", "mds1" = "#FCD225", "N" = "black", "mds6" = "blue")
cols6 <- c("mds6" = "blue", "mds4" = "#a52a2a", "mds3" = "#300060", "mds2" = "#C92D59", "mds1" = "#FCD225", "N" = "black", "mds6" = "blue")

plot_dd<-na.omit(plot_dd)

## p1 is coordinates 1 and 2
p1<-ggplot(plot_dd,aes(x=mds1,y=mds2, colour=outlier))+
  geom_point(alpha=0.7, size=2)+
  xlab("Coordinate 1")+
  ylab("Coordinate 2")+
  scale_colour_manual(values=cols1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, family = "Avenir"),
        legend.position="none",
        strip.text = element_text(size=20))+
  geom_hline(yintercept=mds2_cutoffs,linetype="dashed")+
  geom_vline(xintercept=mds1_cutoffs,linetype="dashed")


## p1.1 is coordinates 1 and 3
p1.1<-ggplot(plot_dd,aes(x=mds1,y=mds3, colour=outlier))+
  geom_point(alpha=0.7, size=2)+
  xlab("Coordinate 1")+
  ylab("Coordinate 3")+
  scale_colour_manual(values=cols1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, family = "Avenir"),
        legend.position="none",
        strip.text = element_text(size=20))+
  geom_hline(yintercept=mds3_cutoffs,linetype="dashed")+
  geom_vline(xintercept=mds1_cutoffs,linetype="dashed")



## p2 is mds1 coordinates along the window position
p2<-ggplot(plot_dd,aes(x=window_pos,y=mds1,colour=outlier))+
  geom_point(alpha=0.7, size=1.5)+
  scale_x_continuous(name="Position (Mb)",
                     breaks=seq(0,46000000, by = 1000000),expand=c(0,0),
                     labels=c("0", "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12", "13", "14",
                              "15", "16", "17", "18", "19", "20", "21",
                              "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
                              "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
                              "42", "43", "44", "45", "46"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylab("MDS 1")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(family="Avenir", size=30, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20),
        axis.text.x = element_blank(),
        #     strip.text = element_text(size=24),
        legend.position="none")+
  scale_colour_manual(values=cols1)


## p2.2 is MDS2 coordinates along window position
p2.2<-ggplot(plot_dd,aes(x=window_pos,y=mds2,colour=outlier))+
  geom_point(alpha=0.7, size=1.5)+
  scale_x_continuous(name="Position (Mb)",
                     breaks=seq(0,46000000, by = 1000000),expand=c(0,0),
                     labels=c("0", "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12", "13", "14",
                              "15", "16", "17", "18", "19", "20", "21",
                              "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
                              "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
                              "42", "43", "44", "45", "46"))+
  ylab("MDS 2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(family="Avenir", size=30, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20),
        axis.text.x = element_blank(),
        strip.text = element_text(size=20),
        legend.position="none")+
  scale_colour_manual(values=cols2)


## p2.3 is MDS3 along window position
p2.3<-ggplot(plot_dd,aes(x=window_pos,y=mds3,colour=outlier))+
  geom_point(alpha=0.7, size=1.5)+
  scale_x_continuous(name="Position (Mb)",
                     breaks=seq(0,46000000, by = 1000000),expand=c(0,0),
                     labels=c("0", "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12", "13", "14",
                              "15", "16", "17", "18", "19", "20", "21",
                              "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
                              "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
                              "42", "43", "44", "45", "46"))+
  ylab("MDS 3")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(family="Avenir", size=30, vjust=0.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20),
        axis.text.x=element_text(family = "Avenir", size=12),
        strip.text = element_text(size=20),
        legend.position="none")+
  scale_colour_manual(values=cols3)


## Try out various plot combinations
plot_1 <- grid.arrange(ggarrange(p1,p1.1,nrow=2), ggarrange(p2,p2.2,p2.3,nrow=3),ncol=2,widths=c(2,3))

ggsave(plot_1, width = 40, height = 30, units = "cm", file = "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/lostruct_LG1_10bp_windows.png")


############################### Chromosome 12 ###################################
## new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/lostruct")

## Load packages
lib = c("parallel","data.table","tidyr","dplyr","ggplot2","lostruct","ggpubr","vcfR", "gridExtra", "extrafont")
lapply(lib, library, character.only=TRUE)

# Define various groupings, can either be split into males and females, or all 
female_vec<-c("Paria_F12",  "Paria_F16",  "Paria_F2",  "Paria_F2_RNA",  "Paria_F8",  "UM10",  "UM11",  "UM12",  "UM13",  "UM14",  "UM15",  "UM16",  "UM17",  "UM18",  "UM8",  "UM9")
male_vec<-c("Paria_M12",  "Paria_M2",  "Paria_M5",  "Paria_M9",  "UM1",  "UM2",  "UM3",  "UM5",  "UM6",  "UM7")

vcf<- "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/natural_data/chr12.nat_data_31893_updateSTAR.DEFAULT.vcf.gz"
window_N<-10
chr<-"chr12"
output<-"chr12_nat_data_10bp"

# Make the input
system(paste0("bcftools view -r ",chr," ",vcf," > localPCA_tmp.vcf"),wait=TRUE)

# Read in VCF in popGenome for info on BP locations for unphased vcf
vcf_PG<-read.vcfR("localPCA_tmp.vcf")
BP1<-as.integer(vcf_PG@fix[,2])
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]
inds_vec<-colnames(vcf_PG@gt)[2:length(colnames(vcf_PG@gt))]

# Filter for inds we want
#inds_to_keep<-unlist(lapply(pops,function(x){
#  inds2<-inds_vec[inds_vec %in% male_vec]
#return(inds2[grep(x,inds2)])
#}))


# vcf_PG@gt<-vcf_PG@gt[,c(colnames(vcf_PG@gt)[1],inds_to_keep)]
# Get males and female columns
#females<-match(female_vec,inds_vec)[is.na(match(female_vec,inds_vec))==FALSE]
#males<-setdiff(1:length(inds_vec),females)

# Return logical vector for polymorphic sites
invariant_filter<-is.polymorphic(vcf_PG,na.omit=T)

# Read in VCF for localPCA
vcf_in<-read_vcf("localPCA_tmp.vcf")

# And tidy up...
system("rm -f localPCA_tmp.vcf")

# Filter for invariants and individuals
#vcf_in<-vcf_in[invariant_filter,(inds_vec %in% male_vec)]
BP1<-BP1[invariant_filter]
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]

# Read in eigen_windows, win=number of rows of matrix, k = number of eigenvector/value pairs to return, mc.cores compatible but does seem to crash sometimes with mc.cores
eigenstuff <- eigen_windows(vcf_in, win=window_N, k=3)

# Calculate the distances, npc=N of PCs computed, mc.core compatible
windist <- pc_dist(eigenstuff, npc=3,mc.cores=detectCores()-1)

# Remove any windows for which we are getting NA values
to_keep_logical<-rowSums(is.na(windist)) != ncol(windist)
to_keep<-(1:ncol(windist))[to_keep_logical]
windist2<-windist[to_keep,to_keep]

# Principal Coordinates Analysis, k = dimension of space for data to represented in, eig = return eigenvalues
fit2d <- cmdscale(windist2, eig=TRUE, k=3)

## Look for drop in eigvals
eigenvals <- fit2d$eig/sum(fit2d$eig)
eigenvals_10 <- eigenvals[eigenvals > 0.01]

# Plot along a chromosome
plot_dd<-data.frame(mds1=fit2d$points[,1],
                    mds2=fit2d$points[,2],
                    mds3=fit2d$points[,3],
                    window_N=1:length(fit2d$points[,1]),
                    window_pos=window_BP[to_keep])


# Calculate cut-offs for outliers by trimming distributions and taking 3*SD of trimmed distribution
trim_mds1<-plot_dd[plot_dd$mds1 > quantile(plot_dd$mds1,probs=0.05) & plot_dd$mds1 < quantile(plot_dd$mds1,probs=0.95),"mds1"]
mds1_cutoffs<-c(mean(trim_mds1)+3*sd(trim_mds1),
                mean(trim_mds1)-3*sd(trim_mds1))

trim_mds2<-plot_dd[plot_dd$mds2 > quantile(plot_dd$mds2,probs=0.05) & plot_dd$mds2 < quantile(plot_dd$mds2,probs=0.95),"mds2"]
mds2_cutoffs<-c(mean(trim_mds2)+3*sd(trim_mds2),
                mean(trim_mds2)-3*sd(trim_mds2))

trim_mds3<-plot_dd[plot_dd$mds3 > quantile(plot_dd$mds3,probs=0.05) & plot_dd$mds3 < quantile(plot_dd$mds3,probs=0.95),"mds3"]
mds3_cutoffs<-c(mean(trim_mds3)+3*sd(trim_mds3),
                mean(trim_mds3)-3*sd(trim_mds3))

##################
# We also want to return the windows at the extremes of the distribution
##################

# Give window IDs
window_ids<-rep(0)
for(j in 1:length(to_keep)){
  i <- to_keep[j]
  window_ids[j]<-paste0(chr,":",((BP1[(i-1)*window_N])+1),"-",BP1[(i)*window_N])
}
######## This was shortening the vector by filtering again on "to_keep"
plot_dd$window_id <- window_ids
######################################################## 

windows_to_keep1<-unique(c(plot_dd[plot_dd$mds1 < mds1_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds1 > mds1_cutoffs[1],"window_id"]))

windows_to_keep2<-unique(c(plot_dd[plot_dd$mds2 < mds2_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds2 > mds2_cutoffs[1],"window_id"]))

windows_to_keep3<-unique(c(plot_dd[plot_dd$mds3 < mds3_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds3 > mds3_cutoffs[1],"window_id"]))

# Reassign to plot_dd
plot_dd$outlier<-rep("N")
plot_dd[plot_dd$window_id %in% windows_to_keep1,"outlier"]<-"mds1"
plot_dd[plot_dd$window_id %in% windows_to_keep2,"outlier"]<-"mds2"
plot_dd[plot_dd$window_id %in% windows_to_keep3,"outlier"]<-"mds3"

# Throw in colour values
cols1 <- c("mds1" = "#FCD225", "mds2" = "#C92D59", "mds3" = "black", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")
cols2 <- c("mds2" = "#C92D59", "mds1" = "#FCD225", "mds3" = "black", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")
cols3 <- c("mds3" = "#300060", "mds2" = "#C92D59", "mds1" = "#FCD225", "N" = "black", "mds4" = "#a52a2a", "mds5" = "green", "mds6" = "blue")

plot_dd<-na.omit(plot_dd)

## Plotting follows LG1 above
p1<-ggplot(plot_dd,aes(x=mds1,y=mds2, colour=outlier))+
  geom_point(alpha=0.7, size=2)+
  xlab("Coordinate 1")+
  ylab("Coordinate 2")+
  scale_colour_manual(values=cols1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, family = "Avenir"),
        legend.position="none",
        strip.text = element_text(size=20))+
  geom_hline(yintercept=mds2_cutoffs,linetype="dashed")+
  geom_vline(xintercept=mds1_cutoffs,linetype="dashed")

p1.1<-ggplot(plot_dd,aes(x=mds1,y=mds3, colour=outlier))+
  geom_point(alpha=0.7, size=2)+
  xlab("Coordinate 1")+
  ylab("Coordinate 3")+
  scale_colour_manual(values=cols1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, family = "Avenir"),
        legend.position="none",
        strip.text = element_text(size=20))+
  geom_hline(yintercept=mds3_cutoffs,linetype="dashed")+
  geom_vline(xintercept=mds1_cutoffs,linetype="dashed")  				

p2<-ggplot(plot_dd,aes(x=window_pos,y=mds1,colour=outlier))+
  geom_point(alpha=0.7, size=1.5)+
  scale_x_continuous(name="Position (Mb)",
                     breaks=seq(0,46000000, by = 1000000),expand=c(0,0),
                     labels=c("0", "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12", "13", "14",
                              "15", "16", "17", "18", "19", "20", "21",
                              "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
                              "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
                              "42", "43", "44", "45", "46"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylab("MDS 1")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(family="Avenir", size=30, vjust =0.5),
        axis.title.x = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20),
        axis.text.x = element_blank(),
        #     strip.text = element_text(size=24),
        legend.position="none")+
  scale_colour_manual(values=cols1)

p2.2<-ggplot(plot_dd,aes(x=window_pos,y=mds2,colour=outlier))+
  geom_point(alpha=0.7, size=1.5)+
  scale_x_continuous(name="Position (Mb)",
                     breaks=seq(0,46000000, by = 1000000),expand=c(0,0),
                     labels=c("0", "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12", "13", "14",
                              "15", "16", "17", "18", "19", "20", "21",
                              "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
                              "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
                              "42", "43", "44", "45", "46"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylab("MDS 2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(family="Avenir", size=30, vjust =0.5),
        axis.title.x = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20),
        axis.text.x = element_text(family = "Avenir", size=14),
        strip.text = element_text(size=20),
        legend.position="none")+
  scale_colour_manual(values=cols2)
