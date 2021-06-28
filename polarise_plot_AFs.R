# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-as.vector(c("ggplot2","dplyr","tidyr","stringr","data.table","plyr","viridis"))
lapply(lib,library,character.only=TRUE)

## working directory
setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/geno_plots/")

#########################################
## LG1 plot ##
#########################################

## read in the frequencies 
chr1 <-  data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr1_AFs_final.tsv", header=TRUE))

## Polarise to IF9:
## remove pos
tmp <- chr1 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-tmp$IF9_AF

# Polarise to IF9

for(i in colnames(tmp)){
  tmp[IF9_freqs < 0.5,i]<-(1-tmp[IF9_freqs < 0.5,i])
}

## add IF9 polarised back in back in:
chr1_dd <- cbind(chr1,tmp)

colnames(chr1_dd) <- c("chr", "pos", "REF", "ALT", "IF10_AF", "IF6_AF", "IF8_AF", "IF9_AF",
                       "IF10_polar", "IF6_polar", "IF8_polar", "IF9_polar")


## take out normal
# chr1 <- chr1_dd %>% select(pos,IF10_AF,IF9_AF,IF8_AF,IF6_AF)

## take out polarised
chr1 <- chr1_dd %>% dplyr::select(pos,IF10_polar,IF9_polar,IF8_polar,IF6_polar)

#################################
# Within each pool, create the segments
max_pos<-max(chr1$pos)+10
long_AFs<-data.frame(rbindlist(lapply(2:5,function(x){
  
  # Get the pool
  tmp<-data.frame(chr1[,c(1,x)])
  
  # Bin the 0s
  tmp<-tmp[tmp[,2] >= 0,]
  
  # Set up xend
  tmp$xend<-c(tmp$pos[2:nrow(tmp)],max_pos)
  
  # Set up y/yend
  tmp$y<-x-1.5
  tmp$yend<-x-0.5
  
  # Name pool
  tmp$pool<-colnames(chr1)[x]
  
  # Tidy
  colnames(tmp)[2]<-"AFs"
  return(tmp)
})))

## convert AFs to factor for plotting 
long_AFs$chr1<-as.factor(long_AFs$AFs)

# rename for ease
chr1_long_AFs<- long_AFs

# Plot
LG1_geno <- ggplot(chr1_long_AFs,aes(x=pos,y=y))+
  geom_segment(aes(xend=xend,yend=yend,colour=AFs))+
  ## CP detection
  #   geom_segment(data = CP,aes(x =x, y = y, xend = xend, yend = yend), colour="black", size = 1)+
  viridis::scale_color_viridis(name = 'AF', option="magma")+
  # ggtitle("IF9 polarised AFs")+
  theme_classic()+
  theme(axis.text.x=element_text(family = "Avenir", size=18),
        axis.ticks.x=element_line(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20, angle=90, hjust=0.5,
                                 colour=c("#3CBB75", "#481567", "#2D708E", "#B22222")),
        axis.title.x=element_text(family="Avenir", size=30, vjust =-2),
        axis.title.y=element_text(size=30, family="Avenir"),
        legend.title = element_blank(),
        #    legend.key.size = unit(2,"line"),
        #    legend.text = element_text("Avenir", size=20),
        legend.position = "none")+
  scale_x_continuous(name="LG1 (Mb)",
                     breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))+
  scale_y_continuous(name = "AFs", 
                     breaks = seq(1,4, by = 1),
                     labels=c("Iso\nY10", "Iso\nY9", "Iso\nY8", "Iso\nY6"))

ggsave(width = 30, height = 12, units = "cm", dpi=800, "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/LG1_AF_plot_Isoy9-polarised.png", LG1_geno)


#########################################
## LG12 plot ##
#########################################

## read in LG12 data:
chr12 <-  data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr12_AFs_final.tsv", header=TRUE))

## Polarise to IF9:
## remove pos
tmp <- chr12 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-tmp$IF9_AF

## polarise to Iso-Y9
for(i in colnames(tmp)){
  tmp[IF9_freqs < 0.5,i]<-(1-tmp[IF9_freqs < 0.5,i])
}

## add IF9 polarised back in back in:
chr12_dd <- cbind(chr12,tmp)

colnames(chr12_dd) <- c("chr", "pos", "REF", "ALT", "IF10_AF", "IF6_AF", "IF8_AF", "IF9_AF",
                        "IF10_polar", "IF6_polar", "IF8_polar", "IF9_polar")

## take out polarised
chr12 <- chr12_dd %>% dplyr::select(pos,IF10_polar,IF9_polar,IF8_polar,IF6_polar)

#################################
# create segments
max_pos<-max(chr12$pos)+10
long_AFs<-data.frame(rbindlist(lapply(2:5,function(x){
  
  # Get the pool
  tmp<-data.frame(chr12[,c(1,x)])
  
  # Bin the 0s
  tmp<-tmp[tmp[,2] >= 0,]
  
  # Set up xend
  tmp$xend<-c(tmp$pos[2:nrow(tmp)],max_pos)
  
  # Set up y/yend
  tmp$y<-x-1.5
  tmp$yend<-x-0.5
  
  # Name pool
  tmp$pool<-colnames(chr12)[x]
  
  # Tidy
  colnames(tmp)[2]<-"AFs"
  return(tmp)
})))

## convert AFs to factor for plotting 
long_AFs$chr12<-as.factor(long_AFs$AFs)

# rename for plotting
chr12_long_AFs <- long_AFs


## plot
LG12_geno <- ggplot(chr12_long_AFs,aes(x=pos,y=y))+
  geom_segment(aes(xend=xend,yend=yend,colour=AFs))+
  viridis::scale_color_viridis(name = 'AF', option="magma")+
  theme_classic()+
  theme(axis.text.x=element_text(family = "Avenir", size=18),
        axis.ticks.x=element_line(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20, angle=90, hjust=0.5,
                                 colour=c("#3CBB75", "#481567", "#2D708E", "#B22222")),
        axis.title.x=element_text(family="Avenir", size=30, vjust =-2),
        axis.title.y=element_text(size=30, family="Avenir"),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_x_continuous(name="LG12 (Mb)",
                     breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))+
  scale_y_continuous(name = "AFs", 
                     breaks = seq(1,4, by = 1),
                     labels=c("Iso\nY10", "Iso\nY9", "Iso\nY8", "Iso\nY6"))

ggsave(width = 30, height = 12, units = "cm", dpi=800, "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/LG12_AF_plot_Isoy9-polarised.png", LG1_geno)


## Can also look at folded AFs, e.g.:

chr12 <-  data.frame(fread("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/poolfstat/outputs/final/chr12_AFs_final.tsv", header=TRUE))

## fold AFs
chr12$IF10_folded = ifelse(chr12$IF10_AF < 0.5, abs(chr12$IF10_AF - 1), chr12$IF10_AF )
chr12$IF6_folded = ifelse(chr12$IF6_AF < 0.5, abs(chr12$IF6_AF - 1), chr12$IF6_AF )
chr12$IF8_folded = ifelse(chr12$IF8_AF < 0.5, abs(chr12$IF8_AF - 1), chr12$IF8_AF )
chr12$IF9_folded = ifelse(chr12$IF9_AF < 0.5, abs(chr12$IF9_AF - 1), chr12$IF9_AF )


## Polarise to IF9:
## remove pos
tmp <- chr12 %>% dplyr::select(IF10_AF, IF6_AF, IF8_AF, IF9_AF)

## save the IF9 freqs:
IF9_freqs<-tmp$IF9_AF

for(i in colnames(tmp)){
  tmp[IF9_freqs < 0.5,i]<-(1-tmp[IF9_freqs < 0.5,i])
}

## add IF9 polarised back in back in:
chr12_dd <- cbind(chr12,tmp)

colnames(chr12_dd) <- c("chr", "pos", "REF", "ALT", "IF10_AF", "IF6_AF", "IF8_AF", "IF9_AF",
                        "IF10_folded", "IF6_folded", "IF8_folded", "IF9_folded",
                        "IF10_polar", "IF6_polar", "IF8_polar", "IF9_polar")


## take out folded
chr12 <- chr12_dd %>% dplyr::select(pos,IF10_folded,IF9_folded,IF8_folded,IF6_folded)


#################################
# create segments
max_pos<-max(chr12$pos)+10
long_AFs<-data.frame(rbindlist(lapply(2:5,function(x){
  
  # Get the pool
  tmp<-data.frame(chr12[,c(1,x)])
  
  # Bin the 0s
  tmp<-tmp[tmp[,2] >= 0,]
  
  # Set up xend
  tmp$xend<-c(tmp$pos[2:nrow(tmp)],max_pos)
  
  # Set up y/yend
  tmp$y<-x-1.5
  tmp$yend<-x-0.5
  
  # Name pool
  tmp$pool<-colnames(chr12)[x]
  
  # Tidy
  colnames(tmp)[2]<-"AFs"
  return(tmp)
})))


long_AFs$chr12<-as.factor(long_AFs$AFs)
chr12_long_AFs <- long_AFs


## plot
LG12_AFs_folded <- ggplot(chr12_long_AFs,aes(x=pos,y=y))+
  geom_segment(aes(xend=xend,yend=yend,colour=AFs))+
  viridis::scale_color_viridis(name = 'AF', option="magma")+
  theme_classic()+
  theme(axis.text.x=element_text(family = "Avenir", size=18),
        axis.ticks.x=element_line(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_text(family = "Avenir", size=20, angle=90, hjust=0.5,
                                 colour=c("#3CBB75", "#481567", "#2D708E", "#B22222")),
        axis.title.x=element_text(family="Avenir", size=30, vjust =-2),
        axis.title.y=element_text(size=30, family="Avenir"),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_x_continuous(name="LG12 (Mb)",
                     breaks=seq(0,34000000, by = 1000000),
                     labels=seq(0,34, by = 1))+
  scale_y_continuous(name = "AFs", 
                     breaks = seq(1,4, by = 1),
                     labels=c("Iso\nY10", "Iso\nY9", "Iso\nY8", "Iso\nY6"))


