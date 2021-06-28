# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat")

lib<-c("ggplot2", "cowplot", "readr", "tidyr", "stringr", "data.table")
lapply(lib,library,character.only=T)

############
## Chr1 ##
############

## read in data:
chr1 <- data.frame(fread("./outputs/final/chr1_raw_pairwise_fst.tsv", header=TRUE))

chr1[is.na(chr1)] <- 0

## now subset for each line:
IF10_vs_IF6 <- chr1 %>% dplyr::select(chr,pos,IF10_vs_IF6) %>% set_colnames(c("chr","pos","fst"))
IF10_vs_IF8 <- chr1 %>% dplyr::select(chr,pos,IF10_vs_IF8) %>% set_colnames(c("chr","pos","fst"))
IF10_vs_IF9 <- chr1 %>% dplyr::select(chr,pos,IF10_vs_IF9) %>% set_colnames(c("chr","pos","fst"))
IF6_vs_IF8 <- chr1 %>% dplyr::select(chr,pos,IF6_vs_IF8) %>% set_colnames(c("chr","pos","fst"))
IF6_vs_IF9 <- chr1 %>% dplyr::select(chr,pos,IF6_vs_IF9) %>% set_colnames(c("chr","pos","fst"))
IF8_vs_IF9 <- chr1 %>% dplyr::select(chr,pos,IF8_vs_IF9) %>% set_colnames(c("chr","pos","fst"))

## add spline chr1
x=IF10_vs_IF6$pos
y=IF10_vs_IF6$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF6$spline <- spline_fst

x=IF10_vs_IF8$pos
y=IF10_vs_IF8$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF8$spline <- spline_fst

x=IF10_vs_IF9$pos
y=IF10_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF9$spline <- spline_fst

x=IF6_vs_IF8$pos
y=IF6_vs_IF8$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF6_vs_IF8$spline <- spline_fst

x=IF6_vs_IF9$pos
y=IF6_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF6_vs_IF9$spline <- spline_fst

x=IF8_vs_IF9$pos
y=IF8_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF8_vs_IF9$spline <- spline_fst

## Load chr1 regions of differentation
chr1_regions <- read.table("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/chr1_regions/chr1_region_coords.txt", header=T)
chr1_regions <- as_tibble(chr1_regions)

# Function for plotting all pairwise
plot_fst_chr1 <- function(dataset1,pop){
  ggplot(dataset1) +
    geom_point(mapping=aes(x =pos, y = fst), colour = "#1e272c",size = 0.1, alpha = 0.2) +
    geom_line(mapping=aes(x=pos, y=spline),colour="#fab201", size = 1.8)+
    # geom_smooth(mapping=aes(x=pos, y=fst),colour="gold2", size = 1, span=0.000001)+
    geom_rect(alpha = 0.3, fill=c("grey", "grey", "grey"), colour = "white",
              linetype=1, size = 0,
              data=chr1_regions,
              mapping=aes(xmin=start,
                          xmax=end,ymin=ymin,
                          ymax=c(0.99, 0.99, 0.99)))+
    ylab(expression(F[ST]))+
    xlab("Position") +
    scale_y_continuous(breaks=seq(-0.2, 1, by = 0.2), limits = c(0,1)) +
    scale_color_manual(values=c("#f4b41a")) + guides(colour=FALSE)+
    scale_x_continuous(name="Position (Mb)",
                       breaks=seq(0,46000000, by = 5000000),expand=c(0,0),
                       labels=seq(0,46, by = 5))+
    theme_bw() +
    theme(plot.title = element_text(family ="Avenir", size=30, hjust =0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.4),
          axis.title.y = element_blank(),
          axis.text.y = element_text(family = "Avenir", size=30),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(family = "Avenir", size=30),
          axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"))
}


IF10_IF6_plot1 <- plot_fst_chr1(IF10_vs_IF6) + ggtitle ("Iso-Y6 vs Iso-Y10") + theme(axis.text.x = element_blank())
IF10_IF8_plot1 <- plot_fst_chr1(IF10_vs_IF8) + ggtitle ("Iso-Y10 vs Iso-Y8")
IF10_IF9_plot1 <- plot_fst_chr1(IF10_vs_IF9) + ggtitle ("Iso-Y10 vs Iso-Y9") + theme(axis.text.x = element_blank())
IF6_IF8_plot1 <- plot_fst_chr1(IF6_vs_IF8) + ggtitle ("Iso-Y6 vs Iso-Y8")
IF6_IF9_plot1 <- plot_fst_chr1(IF6_vs_IF9) + ggtitle ("Iso-Y6 vs Iso-Y9") + theme(axis.text.x = element_blank())
IF8_IF9_plot1 <- plot_fst_chr1(IF8_vs_IF9) + ggtitle ("Iso-Y8 vs Iso-Y9")


plots_chr1 <-list(IF10_IF6_plot1,IF6_IF9_plot1,IF6_IF8_plot1,IF10_IF9_plot1, IF10_IF8_plot1, IF8_IF9_plot1)

## Create empty matrix for chr1 to plot in pairwise triangles
m1 <- matrix(NA, 3, 3)
m1[lower.tri(m1, diag = T)] <- 1:6
m1
tri_fst_plot_chr1 <- grid.arrange(grobs = plots_chr1, layout_matrix = m1)

ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/chr1_pairwise_plot_nonwindowed.png", tri_fst_plot_chr1, width = 40, height = 40, units = "cm", dpi = 400)

############
## Chr12 ##
############

## read in again:
chr12 <- data.frame(fread("./outputs/final/chr12_raw_pairwise_fst.tsv", header=TRUE))

chr12[is.na(chr12)] <- 0

## now subset for each line:
IF10_vs_IF6 <- chr12 %>% dplyr::select(chr,pos,IF10_vs_IF6) %>% set_colnames(c("chr","pos","fst"))
IF10_vs_IF8 <- chr12 %>% dplyr::select(chr,pos,IF10_vs_IF8) %>% set_colnames(c("chr","pos","fst"))
IF10_vs_IF9 <- chr12 %>% dplyr::select(chr,pos,IF10_vs_IF9) %>% set_colnames(c("chr","pos","fst"))
IF6_vs_IF8 <- chr12 %>% dplyr::select(chr,pos,IF6_vs_IF8) %>% set_colnames(c("chr","pos","fst"))
IF6_vs_IF9 <- chr12 %>% dplyr::select(chr,pos,IF6_vs_IF9) %>% set_colnames(c("chr","pos","fst"))
IF8_vs_IF9 <- chr12 %>% dplyr::select(chr,pos,IF8_vs_IF9) %>% set_colnames(c("chr","pos","fst"))

## add spline chr12
x=IF10_vs_IF6$pos
y=IF10_vs_IF6$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF6$spline <- spline_fst

x=IF10_vs_IF8$pos
y=IF10_vs_IF8$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF8$spline <- spline_fst

x=IF10_vs_IF9$pos
y=IF10_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF10_vs_IF9$spline <- spline_fst

x=IF6_vs_IF8$pos
y=IF6_vs_IF8$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF6_vs_IF8$spline <- spline_fst

x=IF6_vs_IF9$pos
y=IF6_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF6_vs_IF9$spline <- spline_fst

x=IF8_vs_IF9$pos
y=IF8_vs_IF9$fst
spline_fst=predict(smooth.spline(x,y,spar=0.5),x)$y
IF8_vs_IF9$spline <- spline_fst

# Function for plotting all pairwise
plot_fst_chr12 <- function(dataset1,pop){
  ggplot(dataset1) +
    geom_point(mapping=aes(x =pos, y = fst), colour = "#1e272c",size = 0.1, alpha = 0.2) +
    geom_line(mapping=aes(x=pos, y=spline),colour="#fab201", size = 1.8)+
    ylab(expression(F[ST]))+
    xlab("Position") +
    scale_y_continuous(breaks=seq(-0.2, 1, by = 0.2), limits = c(0,1)) +
    scale_color_manual(values=c("#f4b41a")) + guides(colour=FALSE)+
    scale_x_continuous(name="Position (Mb)",
                       breaks=seq(0,46000000, by = 5000000),expand=c(0,0),
                       labels=seq(0,46, by = 5))+
    theme_bw() +
    theme(plot.title = element_text(family ="Avenir", size=30, hjust =0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.4),
          axis.title.y = element_blank(),
          axis.text.y = element_text(family = "Avenir", size=30),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(family = "Avenir", size=30),
          axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"))
}



IF10_IF6_plot1 <- plot_fst_chr12(IF10_vs_IF6) + ggtitle ("Iso-Y6 vs Iso-Y10") + theme(axis.text.x = element_blank())
IF10_IF8_plot1 <- plot_fst_chr12(IF10_vs_IF8) + ggtitle ("Iso-Y10 vs Iso-Y8")
IF10_IF9_plot1 <- plot_fst_chr12(IF10_vs_IF9) + ggtitle ("Iso-Y10 vs Iso-Y9") + theme(axis.text.x = element_blank())
IF6_IF8_plot1 <- plot_fst_chr12(IF6_vs_IF8) + ggtitle ("Iso-Y6 vs Iso-Y8")
IF6_IF9_plot1 <- plot_fst_chr12(IF6_vs_IF9) + ggtitle ("Iso-Y6 vs Iso-Y9") + theme(axis.text.x = element_blank())
IF8_IF9_plot1 <- plot_fst_chr12(IF8_vs_IF9) + ggtitle ("Iso-Y8 vs Iso-Y9")


plots_chr12 <-list(IF10_IF6_plot1,IF6_IF9_plot1,IF6_IF8_plot1,IF10_IF9_plot1, IF10_IF8_plot1, IF8_IF9_plot1)

## Create empty matrix for chr12 in pairwise triangle
m1 <- matrix(NA, 3, 3)
m1[lower.tri(m1, diag = T)] <- 1:6
m1
tri_fst_plot_chr12 <- grid.arrange(grobs = plots_chr12, layout_matrix = m1)

ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/chr12_pairwise_plot_nonwindowed.png", tri_fst_plot_chr12, width = 40, height = 40, units = "cm", dpi = 400)

