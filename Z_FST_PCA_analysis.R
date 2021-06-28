# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat")

lib<-c("poolfstat","ggplot2", "cowplot", "readr", "tidyr", "stringr", "data.table", "ggridges")
lapply(lib,library,character.only=T)


## Read in so we don't have to repeat the above:
dd <- data.frame(fread("./outputs/poolseq_lifted_pairwise_fst.tsv", header = TRUE))

## remove NA
dd <- drop_na(dd)

# PCA with scale and center (same as Zscores)
Z_PCA<-prcomp(dd[,5:10],scale=T,center=T)

# Add PC scores to dd
dd$Z_PC1<-Z_PCA$x[,1]
dd$Z_PC2<-Z_PCA$x[,2]
dd$Z_PC3<-Z_PCA$x[,3]
dd$chr<-factor(dd$chr,levels=unique(dd$chr))

# Summary
summary(Z_PCA)

# Summarise by chr (PC1)
FST_sum <- dd %>% group_by(chr) %>% dplyr::summarise(Z_FST=mean(Z_PC1))

write.table(file = "./outputs/Z_FST_1_SUM_whole_genome.tsv", FST_sum, sep = "\t", quote = FALSE, col.names = TRUE)


# Set up bits for plotting (don't need all of this)
dd <- dd %>%
  mutate(chr_num = str_extract(chr, "[^r]+$")) %>%
  mutate(chr_name = ifelse(str_detect(chr, "chr"), chr, "unplaced_chr0")) %>%
  mutate(color = ifelse(chr_num %in% seq(0, 23, 2), "even", "odd")) %>%
  mutate(number = 1) %>% # this is a column to be used latter for a running number by cumsum
  filter(str_detect(chr,"chr")) %>%
  group_by(chr) %>% mutate(rownum = cumsum(number))

dd$chr_num <- factor(dd$chr_num,
                      levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                 "9", "10", "11", "12", "13", "14", "15",
                                 "16", "17", "18", "19", "20", "21", "22", "23"))

## plot PCA scores
g1 <-   ggplot(dd)+
  geom_density_ridges_gradient(aes(x=`Z_PC1`,y=`chr_num`, group=chr_num, fill = stat(x)),
                               rel_min_height = 0.01, scale=3, quantile_lines = TRUE, quantiles = 0.5)+
  geom_vline(xintercept = 3, linetype="dashed", alpha=0.8)+
  scale_fill_viridis_c(name = expression(paste("Z-F" [ST], " PC1")), option = "plasma", direction=-1) +
  scale_x_continuous(breaks=seq(-2, 10, by = 1), labels=seq(-2, 10,1), limits=c(-2,10)) +
  xlab(expression(paste("Z-F" [ST], " PC1")))+
  ylab("Chromosome")+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title.y = element_text(family="Avenir",size=40,vjust=0.5),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(family = "Avenir", size = 18),
        legend.key.size = unit(0.8, "cm"),
        axis.ticks.x = element_line(),
        axis.text.y = element_text(family="Avenir", size = 24),
        axis.text.x = element_text(family="Avenir", size = 24),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=28, family="Avenir"))

ggsave(g1, width = 40, height = 20, units = "cm", dpi = 400,
       file = "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/Figure_2A_density_ZPC.png")


## pull out LG1
chr1 <- dd %>% filter(chr == "chr1")

chr1 <- drop_na(chr1)

# PCA with scale and center (same as Zscores) ## try with no scale
Z_PCA<-prcomp(chr1[,5:10],scale=T,center=T)

chr1$Z_PC1<-Z_PCA$x[,1]
chr1$Z_PC2<-Z_PCA$x[,2]
chr1$Z_PC3<-Z_PCA$x[,3]

summary(Z_PCA)

Z_PCA$rotation

## set up spline
x=chr1$pos
y=chr1$Z_PC1
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr1$spline_PC1 <- spline

x=chr1$pos
y=chr1$Z_PC2
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr1$spline_PC2 <- spline

x=chr1$pos
y=chr1$Z_PC3
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr1$spline_PC3 <- spline

## Add in shading for regions 1, 2 and 3 where regions have been calculated using changepoint detection
chr1_regions <- read.table("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/chr1_regions/chr1_region_coords.txt", header=T)
chr1_regions <- as_tibble(chr1_regions)


LG1_PC1 <- ggplot(chr1)+
  geom_point(aes(x=pos,y=Z_PC1),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(aes(x=pos,y=spline_PC1),colour="#f4b41a", size = 2.2, alpha=0.8)+
  geom_rect(data=chr1_regions,alpha = 0.3, fill=c("grey", "grey", "grey"), colour = "white",
            linetype=1, size = 0,
            mapping=aes(xmin=start,xmax=end,
                        ymin=c(-2,-2,-2), ymax=c(10, 10, 10)))+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC1")))+
  xlab("LG1 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-2, 10, by = 1), labels=seq(-2,10,1), limits=c(-2,10)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))


LG1_PC2 <- ggplot(chr1)+
  geom_point(aes(x=pos,y=Z_PC2),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(aes(x=pos,y=spline_PC2),colour="#f4b41a", size = 2.2, alpha=0.8)+
  geom_rect(data=chr1_regions,alpha = 0.3, fill=c("grey", "grey", "grey"), colour = "white",
            linetype=1, size = 0,
            mapping=aes(xmin=start,xmax=end,
                        ymin=c(-6,-6,-6), ymax=c(3, 3, 3)))+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC2")))+
  xlab("LG1 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-6, 3, by = 1), labels=seq(-6,3,1), limits=c(-6,3)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))


LG1_PC3 <- ggplot(chr1)+
  geom_point(aes(x=pos,y=Z_PC3),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(aes(x=pos,y=spline_PC3),colour="#f4b41a", size = 2.2, alpha=0.8)+
  geom_rect(data=chr1_regions,alpha = 0.3, fill=c("grey", "grey", "grey"), colour = "white",
            linetype=1, size = 0,
            mapping=aes(xmin=start,xmax=end,
                        ymin=c(-7,-7,-7), ymax=c(5, 5, 5)))+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC3")))+
  xlab("LG1 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_text(family = "Avenir", size=18),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-7, 5, by = 1), labels=seq(-7,5,1), limits=c(-7,5)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))


LG1_PC2_PC3 <- cowplot::plot_grid(LG1_PC2,LG1_PC3, ncol=1)

ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/PC2_PC3_LG1_loadings.png", LG1_PC2_PC3 , width = 40, height = 20, units = "cm", dpi = 400)


################################################
# pull out LG12
################################################
chr12 <- dd %>% filter(chr == "chr12")
chr12 <- drop_na(chr12)

# PCA with scale and center (same as Zscores) ## try with no scale
Z_PCA<-prcomp(chr12[,5:10],scale=T,center=T)

chr12$Z_PC1<-Z_PCA$x[,1]
chr12$Z_PC2<-Z_PCA$x[,2]
chr12$Z_PC3<-Z_PCA$x[,3]

summary(Z_PCA)

Z_PCA$rotation

## set up spline
x=chr12$pos
y=chr12$Z_PC1
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr12$spline_PC1 <- spline

## set up spline
x=chr12$pos
y=chr12$Z_PC2
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr12$spline_PC2 <- spline

## set up spline
x=chr12$pos
y=chr12$Z_PC3
spline=predict(smooth.spline(x,y,spar=0.5),x)$y
chr12$spline_PC3 <- spline


LG12_PC1 <- ggplot(chr12,aes(x=pos,y=spline_PC1))+
  geom_point(aes(y=Z_PC1),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(colour="#f4b41a", size = 2.2, alpha=0.8)+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC2")))+
  xlab("LG12 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-2, 10, by = 1), labels=seq(-2,10,1), limits=c(-2,10)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))

LG12_PC2 <- ggplot(chr12,aes(x=pos,y=spline_PC2))+
  geom_point(aes(y=Z_PC2),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(colour="#f4b41a", size = 2.2, alpha=0.8)+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC2")))+
  xlab("LG12 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-3, 5, by = 1), labels=seq(-3,5,1), limits=c(-3,5)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))

LG12_PC3 <- ggplot(chr12,aes(x=pos,y=spline_PC3))+
  geom_point(aes(y=Z_PC3),alpha=0.4, size = 0.1,colour = "#1e272c")+
  geom_line(colour="#f4b41a", size = 2.2, alpha=0.8)+
  theme_classic()+
  ylab(expression(paste("Z-F" [ST], " PC3")))+
  xlab("LG12 (Mb)")+
  theme(plot.title = element_blank(),
        axis.text.x=element_text(size=14, family="Avenir"),
        axis.title.y=element_text(size=30, family="Avenir"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, family = "Avenir"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(-6, 2, by = 1), labels=seq(-6,2,1), limits=c(-6,2)) +
  scale_x_continuous(breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))

LG12_PC2_PC3 <- cowplot::plot_grid(LG12_PC2,LG12_PC3, ncol=1)


ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/PC2_PC3_LG12_loadings.png", LG12_PC2_PC3 , width = 40, height = 20, units = "cm", dpi = 400)

############################################################
####### Calculate alphas for sig SNPs #####################
############################################################

## calculate alpha PC1
my_threshold_lower <- quantile(dd$Z_PC1, 0.05)
my_threshold_upper <- quantile(dd$Z_PC1, 0.95)

## use the upper 95%, 2.94, rounded to 3 as your top limit

PC1 <- dd %>% mutate(PC1_outlier = ifelse(Z_PC1 > 3, "outlier", "no")) %>%
  dplyr::filter(PC1_outlier == "outlier")

## count per chromosomes:
PC1_sum <- PC1 %>% group_by(chr) %>% dplyr::summarise(Z_FST=dplyr::count(PC1_outlier))

## count per chromosomes:
PC1_sum <- PC1 %>% group_by(chr) %>% count(Z_FST=(PC1_outlier))

write.table(file = "./outputs/PC1_loadings_per_chr_scaff_new.tsv", PC1_sum, sep = "\t", quote = FALSE, col.names = TRUE)

## Get sum of total outliers:
tot <- sum(PC1_sum$n) # 190259

options(scipen = 999)
PC1_sum <- PC1_sum %>% mutate(perct_r=round((n/tot)*100))

write.table(file = "./outputs/PC1_loadings_per_chr_scaff_new.tsv", PC1_sum, sep = "\t", quote = FALSE, col.names = TRUE)

### Read in table with each percentage and expected

dd <- read.table(file="./outputs/PC_zscores_percentage_chrs.txt", header = TRUE, sep="\t")

## lock in chr order:
dd$chr <- factor(dd$chr, levels = dd$chr)

dd2 <- tidyr::gather(dd,key = type, value = PC,PC1:PC3)

ggplot(dd2)+
  geom_col(aes(x=chr,y=PC, fill = type), position = "dodge")


### Make a plot with all the scaffolds and chrs:

dd <- read.table("LGs_scaffolds_outilers.txt", header=TRUE, sep="\t")

chr <- dd %>% filter(grepl("LG",Chromosome))

# lock in factor level order
chr$Chromosome <- factor(chr$Chromosome, levels = chr$Chromosome)


## chrs:
chrs <- ggplot(chr)+
  geom_col(aes(x=Chromosome,y=Proportion.of.high.Z.FST.PC1.above.3))+
  ylab(expression(paste("% SNPs > critical ", alpha)))+
  xlab("Chromosome")+
  theme_bw()+
  theme(axis.text.x=element_text(family="Avenir", size=12),
        axis.text.y=element_text(family="Avenir", size=16),
        axis.title = element_text(family = "Avenir", size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_y_continuous(breaks = seq(0,30, by=1), labels = seq(0,30, by=1))

ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/chromosomes_proportion_alpha.png", chrs, width = 40, height = 20, units = "cm", dpi = 400)


scaff <- dd %>% filter(!grepl("LG",Chromosome))

# lock in factor level order
chr$Chromosome <- factor(chr$Chromosome, levels = chr$Chromosome)

## scafs:
scaf <- ggplot(scaff)+
  geom_col(aes(x=Chromosome,y=Proportion.of.high.Z.FST.PC1.above.3))+
  ylab(expression(paste("% SNPs > critical ", alpha)))+
  xlab("Unplaced genomic scaffolds")+
  theme_bw()+
  theme(axis.text.x=element_text(family="Avenir", size=4, angle = 90),
        axis.text.y=element_text(family="Avenir", size=16),
        axis.title = element_text(family = "Avenir", size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_y_continuous(breaks = seq(0,1, by=0.1), labels = seq(0,1, by=0.1))

ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/scaffolds_proportion_alpha.png", scaf, width = 40, height = 20, units = "cm", dpi = 400)
