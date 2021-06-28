rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-as.vector(c("vcfR","pbapply","PopGenome","cowplot","gridExtra","viridis","grid","data.table","ggplot2","parallel","dplyr","tidyr", "stringr"))
lapply(lib,library,character.only=TRUE)

source("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/R_scripts/pool_pi.R")

#########################################
#### LG12 #####
##########################################

vcf_path <- "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/vcf_files/chr1_variant_filtered.recode.vcf.gz"
vcf_in <- read.vcfR(vcf_path)

# Usage
chr1_pi <- pool_pi(pool_vcf = vcf_in,wind_size=10000,n_cores=2,min_snp = 10)
chr1_pi <- na.omit(chr1_pi)

# Keep only pi in all pools
wind_counts <- table(chr1_pi$start)
to_keep <- names(wind_counts[wind_counts==4])
# Make a new dd
pca_dd_LG1 <- data.frame(IF10=chr1_pi[chr1_pi$pool=="IF10" & chr1_pi$start %in% to_keep,"pi"],
                         IF9=chr1_pi[chr1_pi$pool=="IF9" & chr1_pi$start %in% to_keep,"pi"],
                         IF8=chr1_pi[chr1_pi$pool=="IF8" & chr1_pi$start %in% to_keep,"pi"],
                         IF6=chr1_pi[chr1_pi$pool=="IF6" & chr1_pi$start %in% to_keep,"pi"])

pi_pca_LG1 <- prcomp(pca_dd_LG1,scale. = T,center = T)
summary(pi_pca_LG1)
pi_pca_LG1

# Scores
LG1_scores <- data.frame(pi_pca_LG1$x)
LG1_scores$start <- as.integer(to_keep)

PC1 <- ggplot(LG1_scores,aes(y=PC1,x=start))+
  geom_line()+
  geom_point(data=LG1_PC1,aes(x=start,y=10), colour="red")+
  ylab(expression(paste("Z-",pi," PC1 (74%)")))+
  xlab("LG1 Position (Mb)")+
  theme_classic()+
  theme(plot.title = element_text(family = "Avenir", size=20),
        axis.text.x=element_text(family = "Avenir", size=18),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_text(size=30, family = "Avenir"),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-4, 10, by = 1), labels=seq(-4, 10,1), limits=c(-4,10)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))

## save
ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/z-pi_PC1_LG1.png", PC1, width = 40, height = 20, units = "cm", dpi = 400)


## calculate alpha PC1
my_threshold_lower <- quantile(LG1_scores$PC1, 0.01) # 5% = -2.419399
my_threshold_upper <- quantile(LG1_scores$PC1, 0.99) #95% = 2.821794

## use the upper 95%, i.e. 3 as your top limit

LG1_PC1 <- LG1_scores %>% mutate(PC1_outlier = ifelse(PC1 > my_threshold_upper, "outlier", "no")) %>%
  dplyr::filter(PC1_outlier == "outlier")


#########################################
#### LG12 #####
##########################################

vcf_path <- "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/vcf_files/chr12_variant_filtered_STAR_updated.vcf.gz"
vcf_in <- read.vcfR(vcf_path)

# Usage eg
chr12_pi <- pool_pi(pool_vcf = vcf_in,wind_size=10000,n_cores=4,min_snp = 1)
chr12_pi <- na.omit(chr12_pi)

# Keep only pi in all pools
wind_counts <- table(chr12_pi$start)
to_keep <- names(wind_counts[wind_counts==4])
# Make a new dd
pca_dd_LG12 <- data.frame(IF10=chr12_pi[chr12_pi$pool=="IF10" & chr12_pi$start %in% to_keep,"pi"],
                          IF9=chr12_pi[chr12_pi$pool=="IF9" & chr12_pi$start %in% to_keep,"pi"],
                          IF8=chr12_pi[chr12_pi$pool=="IF8" & chr12_pi$start %in% to_keep,"pi"],
                          IF6=chr12_pi[chr12_pi$pool=="IF6" & chr12_pi$start %in% to_keep,"pi"])

pi_pca_LG12 <- prcomp(pca_dd_LG12,scale. = T,center = T)
summary(pi_pca_LG12)
pi_pca_LG12

# Scores
LG12_scores <- data.frame(pi_pca_LG12$x)
LG12_scores$start <- as.integer(to_keep)


## top PC2 in positive direction
head(LG12_scores[order(-LG12_scores$PC2),]) # 24840000
## top PC2 in negative direction
head(LG12_scores[order(LG12_scores$PC2),]) # 24280000

## top PC1 in positive direction
head(LG12_scores[order(-LG12_scores$PC1),]) # 24270000, also 6180000
## top PC1 in negative direction
head(LG12_scores[order(LG12_scores$PC1),]) # 880000

## calculate alpha PC1
my_threshold_lower <- quantile(LG12_scores$PC1, 0.01) # 1% = -3.361235
my_threshold_upper <- quantile(LG12_scores$PC1, 0.99) #99% = 4.583852

LG12_PC1 <- LG12_scores %>% mutate(PC1_outlier = ifelse(PC1 > my_threshold_upper, "outlier", "no")) %>%
  dplyr::filter(PC1_outlier == "outlier")

## calculate alpha PC2
my_threshold_lower <- quantile(LG12_scores$PC2, 0.01) # 1% = -1.356175 
my_threshold_upper <- quantile(LG12_scores$PC2, 0.99) #99% = 1.308372

LG12_PC2 <- LG12_scores %>% mutate(PC2_outlier = ifelse(PC2 < my_threshold_lower, "outlier", "no")) %>%
  dplyr::filter(PC2_outlier == "outlier")

# Plot the scores - tells us general diversity (PC1), and IF9 diversity (PC2)


PC1 <- ggplot(LG12_scores,aes(y=PC1,x=start))+
  geom_line()+
  ylab(expression(paste("Z-",pi," PC1 (88%)")))+
  xlab("LG12 Position (Mb)")+
  geom_point(data=LG12_PC1,aes(x=start, y=10.5), colour="red")+
  theme_classic()+
  theme(plot.title = element_text(family = "Avenir", size=20),
        axis.text.x=element_text(family = "Avenir", size=18),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-4, 11, by = 1), labels=seq(-4, 11,1), limits=c(-4,11)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))


PC2 <- ggplot(LG12_scores,aes(y=PC2,x=start))+
  geom_line()+
  ylab(expression(paste("Z-",pi," PC2 (7%)")))+
  xlab("LG12 Position (Mb)")+
  geom_point(data=LG12_PC2,aes(x=start, y=-3.5), colour="red")+
  theme_classic()+
  theme(plot.title = element_text(family = "Avenir", size=20),
        axis.text.x=element_text(family = "Avenir", size=18),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_text(size=30, family = "Avenir"),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-3, 3, by = 1), labels=seq(-3, 3,1), limits=c(-4,3)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))


PC1_PC2 <- cowplot::plot_grid(PC1,PC2, ncol=1)

## save
ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/z-pi_PC1_PC2LG12.png", PC1_PC2, width = 40, height = 20, units = "cm", dpi = 400)

cowplot::plot_grid(
  ggplot(scores,aes(y=PC1,x=start))+geom_line()+ggtitle("PC1 (87.79%) - General Pi"),
  ggplot(scores,aes(y=PC2,x=start))+geom_line()+ggtitle("PC2 (7.05%) - IF9 Pi"),
  ncol=1)



########################
# If you want to look at the Z-transformation analysis as well..
# Get PC1
pi_PC1 <- pi_pca$x[,1]
IF9_piZ <- scale(pca_dd$IF9)
plot_dd <- data.frame(start=as.integer(to_keep),
                      IF9_diff=scale(pca_dd$IF9)-pi_PC1,
                      IF10_diff=scale(pca_dd$IF10)-pi_PC1,
                      IF8_diff=scale(pca_dd$IF8)-pi_PC1,
                      IF6_diff=scale(pca_dd$IF6)-pi_PC1)
ggplot(plot_dd,aes(start,IF9_diff))+geom_line()
