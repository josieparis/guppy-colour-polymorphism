rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-as.vector(c("ggpmisc","readr","tidyverse","datatable","gridExtra","dplyr","tidyr","stringr"))
lapply(lib,library,character.only=TRUE)

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/polarising_freqs/")

AF_dd <- data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr1_AFs_final.tsv"))

#############################################################################################################
## Region 2:  9,627,619 - 17,074,870 bp
#############################################################################################################

dd2 <- AF_dd %>% dplyr::filter(pos >= 9627619 & pos <= 17074870)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-dd3$IF9_AF

for(i in colnames(dd3)){
  dd3[IF9_freqs < 0.5,i]<-(1-dd3[IF9_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)

dd <- dd3

## IF6
densityCurve <- ggplot(dd, aes(x=IF6_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

local <- localMins[2:3]

IF6 <- ggplot(dd, aes(x=IF6_AF))+
  geom_density(fill="#B22222", alpha=0.5, lwd=2) + 
  geom_vline(xintercept = local, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("AF")+
  ggtitle("Iso-Y6")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=24),
        axis.title = element_text(family="Avenir", size=34))

ggsave("./figs/IF6_AF_densities_Region2.png", IF6, width = 15, height = 20, units = "cm", dpi = 400)


## IF8
densityCurve <- ggplot(dd, aes(x=IF8_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

IF8 <- ggplot(dd, aes(x=IF8_AF))+
  geom_density(fill="#2d708e", alpha=0.5, lwd=2) + 
  geom_vline(xintercept = 0.6973848, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("AF")+
  ggtitle("Iso-Y8")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=24),
        axis.title = element_text(family="Avenir", size=34))

ggsave("./figs/IF8_AF_densities_Region2.png", IF8, width = 15, height = 20, units = "cm", dpi = 400)


## IF10
densityCurve <- ggplot(dd, aes(x=IF10_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

local <- localMins[5]

IF10 <- ggplot(dd, aes(x=IF10_AF))+
  geom_density(fill="#3cbb75", alpha=0.5, lwd=2) + 
  geom_vline(xintercept = local, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("AF")+
  ggtitle("Iso-Y10")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=24),
        axis.title = element_text(family="Avenir", size=34))

ggsave("./figs/IF10_AF_densities_Region2.png", IF10, width = 15, height = 20, units = "cm", dpi = 400)

## IF9
densityCurve <- ggplot(dd, aes(x=IF9_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

IF9 <- ggplot(dd, aes(x=IF9_AF))+
  geom_density(fill="#481567", alpha=0.5, lwd=1) + 
  #geom_vline(xintercept = localMins, color="grey3", linetype = "dashed")+
  ylab("Density")+
  xlab("AF")+
  ggtitle("Iso-Y9")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=24),
        axis.title = element_text(family="Avenir", size=34))

ggsave("./figs/IF9_AF_densities_Region2.png", IF9, width = 15, height = 20, units = "cm", dpi = 400)

## plot together
plots <-list(IF6,IF8,IF9,IF10)

arranged <- grid.arrange(grobs = plots, nrow = 1)

ggsave("./figs/All_density_plots_Region2.png", arranged, width = 60, height = 20, units = "cm", dpi = 400)

#############################################################################################################
## Region 3:21,944,840 - 24,959,750 bp
#############################################################################################################
AF_dd <- data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr1_AFs_final.tsv"))

dd2 <- AF_dd %>% filter(pos >= 21944840 & pos <= 24959750)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-dd3$IF9_AF

for(i in colnames(dd3)){
  dd3[IF9_freqs < 0.5,i]<-(1-dd3[IF9_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)

dd <- dd3

## IF6
densityCurve <- ggplot(dd, aes(x=IF6_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

local <- localMins[2:3]

IF6 <- ggplot(dd, aes(x=IF6_AF))+
  geom_density(fill="#B22222", alpha=0.5, lwd=2) + 
  geom_vline(xintercept = 0.59, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("Allele Frequency")+
  ggtitle("Iso-Y6")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=20),
        axis.title = element_text(family="Avenir", size=30))


ggsave("./figs/IF6_AF_densities_region3.png", IF6, width = 20, height = 20, units = "cm", dpi = 400)

## IF8
densityCurve <- ggplot(dd, aes(x=IF8_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

IF8 <- ggplot(dd, aes(x=IF8_AF))+
  geom_density(fill="#2d708e", alpha=0.5, lwd=2) + 
  #  geom_vline(xintercept = localMins, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("Allele Frequency")+
  ggtitle("Iso-Y8")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=20),
        axis.title = element_text(family="Avenir", size=30))

ggsave("./figs/IF8_AF_densities_region3.png", IF8, width = 20, height = 20, units = "cm", dpi = 400)


## IF10
densityCurve <- ggplot(dd, aes(x=IF10_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

local <- localMins[5]

IF10 <- ggplot(dd, aes(x=IF10_AF))+
  geom_density(fill="#3cbb75", alpha=0.5, lwd=2) + 
  # geom_vline(xintercept = local, color="grey3", linetype = "dashed", lwd=3)+
  ylab("Density")+
  xlab("Allele Frequency")+
  ggtitle("Iso-Y10")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=20),
        axis.title = element_text(family="Avenir", size=30))

ggsave("./figs/IF10_AF_densities_Region3.png", IF10, width = 20, height = 20, units = "cm", dpi = 400)


## IF9
densityCurve <- ggplot(dd, aes(x=IF9_AF)) + geom_density()

# extract the data from the graph
densityCurveData <- ggplot_build(densityCurve)
# get the indices of the local minima
localMins <- which(ggpmisc:::find_peaks(-densityCurveData$data[[1]]$density) == TRUE)

# get the value of the local minima
localMins <- densityCurveData$data[[1]]$x[localMins]
localMins <- c(-Inf, localMins, +Inf)

IF9 <- ggplot(dd, aes(x=IF9_AF))+
  geom_density(fill="#481567", alpha=0.5, lwd=1) + 
  #geom_vline(xintercept = localMins, color="grey3", linetype = "dashed")+
  ylab("Density")+
  xlab("Allele Frequency")+
  ggtitle("Iso-Y9")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, family="Avenir", size=50),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_line(),
        axis.text = element_text(family="Avenir", size=20),
        axis.title = element_text(family="Avenir", size=30))

ggsave("./figs/IF9_AF_densities_region3.png", IF9, width = 20, height = 20, units = "cm", dpi = 400)

## plot together
plots <-list(IF6,IF8,IF9,IF10)

arranged <- grid.arrange(grobs = plots, nrow = 1)

ggsave("./figs/All_density_plots_region3.png", arranged, width = 60, height = 20, units = "cm", dpi = 400)

