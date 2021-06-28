rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-as.vector(c("cowplot","data.table","ggplot2","parallel","dplyr","tidyr", "stringr"))
lapply(lib,library,character.only=TRUE)

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/polarising_freqs/")

AF_dd <- data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr1_AFs_final.tsv"))

## read in CPD regions file
chr1_regions <- read.table("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/chr1_regions/chr1_region_coords.txt", header=T)
chr1_regions <- as_tibble(chr1_regions)

## extract region 1 (IF6 fixed) from frequencies:
dd2 <- AF_dd %>% filter(pos >= 4079988 & pos <= 5984584)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF6 freqs and polarise to this
IF6_freqs<-dd3$IF6_AF

for(i in colnames(dd3)){
  dd3[IF6_freqs < 0.5,i]<-(1-dd3[IF6_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)

## make long
dd3 <- dd3 %>% gather(pool, AF,IF10_AF:IF9_AF)

## read in gtf for gene locations
genes <- read.table("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/gene_density/chr1_OHR95_genes.tsv")
chr1_genes <- genes[genes$V1 == "chr1",]
chr1_genes <- as_tibble(chr1_genes)
region1_genes <- filter(chr1_genes, V4 >= 4079988 & V5 <= 5984584)
region1_genes <- add_column(region1_genes, y = 1.02)


## plot region 1:
af_plot_region1 <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=AF, colour=pool), size=1.5, alpha=0.2)+
  geom_rug(data = region1_genes, aes(V4), sides = "top", length = unit(0.06, "npc"), size=0.5)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir", size=30),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_colour_manual(values = c("#3CBB75", "#B22222", "#2D708E", "#481567"))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.1),
                     limits =c(0,1.02))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(4000000,6000000, by = 500000),
                     labels=c("4.0", "4.5", "5.0", "5.5", "6.0"))


ggsave("./figs/whole_region1_polarised.png", af_plot_region1, width = 20, height = 20, units = "cm", dpi = 400)

##############################################################################################
### pull out individual AFs for each line for region 2:
head(AF_dd)

## and now the same for region 2:
## extract region 1 (IF9 fixed) from frequencies:
dd2 <- AF_dd %>% filter(pos >= 9627619 & pos <= 17074870)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs and polarise to this:
IF9_freqs<-dd3$IF9_AF

for(i in colnames(dd3)){
  dd3[IF9_freqs < 0.5,i]<-(1-dd3[IF9_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)

## make long
dd3 <- dd3 %>% gather(pool, AF,IF10_AF:IF9_AF)

## genes
region2_genes <- filter(chr1_genes, V4 >= 9627619 & V5 <= 17074870)
region2_genes <- add_column(region2_genes, y = 1.02)


## plot region 2:
af_plot_region2 <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=AF, colour=pool), size=1.2, alpha=0.2)+
  geom_rug(data = region2_genes, aes(V4), sides = "top", length = unit(0.06, "npc"), size=0.5)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir", size=30),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_colour_manual(values = c("#3CBB75", "#B22222", "#2D708E", "#481567"))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.1),
                     limits =c(0,1.02))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(9500000,17500000, by = 1000000),
                     labels=c("9.5", "10.5", "11.5", "12.5", "13.5", "14.5", "15.5",
                              "16.5", "17.5"))

ggsave("./figs/whole_region2_polarised.png", af_plot_region2, width = 20, height = 20, units = "cm", dpi = 400)


###########################################################################################################
## and now the same for region 3:
## extract region 3 (IF9 fixed) from frequencies:
dd2 <- AF_dd %>% filter(pos >= 21944840 & pos <= 24959750)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs and polarise to this:
IF9_freqs<-dd3$IF9_AF

for(i in colnames(dd3)){
  dd3[IF9_freqs < 0.5,i]<-(1-dd3[IF9_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)

## make long
dd3 <- dd3 %>% gather(pool, AF,IF10_AF:IF9_AF)

## genes
region3_genes <- filter(chr1_genes, V4 >= 21944840 & V5 <= 24959750)
region3_genes <- add_column(region3_genes, y = 1.02)

## plot region 3:
af_plot_region3 <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=AF, colour=pool), size=1.5, alpha=0.2)+
  geom_rug(data = region3_genes, aes(V4), sides = "top", length = unit(0.06, "npc"), size=0.5)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir", size=30),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_colour_manual(values = c("#3CBB75", "#B22222", "#2D708E", "#481567"))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.1),
                     limits =c(0,1.02))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(22000000,25000000, by = 1000000),
                     labels=c("22.0", "23.0", "24.0", "25.0"))

ggsave("./figs/whole_region3_polarised.png", af_plot_region3, width = 20, height = 20, units = "cm", dpi = 400)

plot_grid(af_plot_region1,af_plot_region2,af_plot_region3)

##################################################################################################################
### pull out individual AFs for each line for region 2 or region 3:
head(AF_dd)

## and now the same for region 2:
## extract region 1 (IF9 fixed) from frequencies:
dd2 <- AF_dd %>% filter(pos >= 9627619 & pos <= 17074870)
dd2 <- AF_dd %>% filter(pos >= 21944840 & pos <= 24959750)

## remove pos
dd3 <- dd2 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF6 freqs:
IF9_freqs<-dd3$IF9_AF

for(i in colnames(dd3)){
  dd3[IF9_freqs < 0.5,i]<-(1-dd3[IF9_freqs < 0.5,i])
}

## add pos back in:
dd3 <- add_column(dd3,pos=dd2$pos)


## genes
region2_genes <- filter(chr1_genes, V4 >= 9627619 & V5 <= 17074870)
region2_genes <- add_column(region2_genes, y = 1.02)


## plot IF9:
IF9_plot <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=IF9_AF), colour="#481567", size=0.5, alpha=0.6)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir",size=20),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.05),
                     limits =c(0,1))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(22000000,25000000, by = 1000000),
                     labels=c("22", "23", "24", "25"))

## plot IF6:
IF6_plot <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=IF6_AF), colour="#B22222", size=0.5, alpha=0.6)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir",size=20),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.05),
                     limits =c(0,1))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(22000000,25000000, by = 1000000),
                     labels=c("22", "23", "24", "25"))

## plot region:
IF8_plot <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=IF8_AF), colour="#2d708eff", size=0.5, alpha=0.6)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir",size=20),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.05),
                     limits =c(0,1))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(22000000,25000000, by = 1000000),
                     labels=c("22", "23", "24", "25"))
## plot region:
IF10_plot <- ggplot(dd3)+
  geom_jitter(aes(x=pos, y=IF10_AF), colour="#3cbb75ff", size=0.5, alpha=0.6)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(family="Avenir",size=20),
        axis.title.x = element_text(family="Avenir",size=20),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(family = "Avenir", size = 20),
        axis.title.y = element_text(family = "Avenir", size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.4))+
  scale_y_continuous(name = "Allele Frequency",
                     breaks=seq(0,1, by=0.05),
                     limits =c(0,1))+
  scale_x_continuous(name="LG1 Position (Mb)",
                     breaks=seq(22000000,25000000, by = 1000000),
                     labels=c("22", "23", "24", "25"))


## plot together
plots <-list(IF6_plot,IF8_plot,IF9_plot,IF10_plot)

arranged <- grid.arrange(grobs = plots, nrow = 1)

ggsave("./figs/region3_AFs_seperataely_plots.png", arranged, width = 60, height = 20, units = "cm", dpi = 400)


ggsave("./figs/IF9_ONLY_polarised_region3.png", IF9_plot, width = 20, height = 20, units = "cm", dpi = 400)
ggsave("./figs/IF6_ONLY_polarised_region3.png", IF6_plot, width = 20, height = 20, units = "cm", dpi = 400)
ggsave("./figs/IF8_ONLY_polarised._region3.png", IF8_plot, width = 20, height = 20, units = "cm", dpi = 400)
ggsave("./figs/IF10_ONLY_polarised_region3.png", IF10_plot, width = 20, height = 20, units = "cm", dpi = 400)

