rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-c("adagenet","ggplot2", "cowplot", "factoextra", "readr", "tidyr", "stringr", "data.table")
lapply(lib,library,character.only=T)

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/phenotype_data_current/")

output_folder = "figs/"
extracted_color = "extracted_colour/"

#bind the UV csv files with the Vis ones and then revise code belows to analyze 4 channels instead of 3

IF6_Vis_4T_1R = read.csv(paste(extracted_color, "IF6_Vis_4T_1R.csv", sep=""))

IF6_UV_4T_1R = read.csv(paste(extracted_color, "IF6_UV_4T_1R.csv", sep=""))

IF6_Vis_4T_1R = cbind(IF6_Vis_4T_1R, IF6_UV_4T_1R[,4:2479])

# write.csv(IF6_Vis_4T_1R, paste(extracted_color, "IF6_Vis_UV_4T_1R.csv"), row.names = FALSE)


IF8_Vis_4T_1R = read.csv(paste(extracted_color, "IF8_Vis_4T_1R.csv", sep=""))

IF8_UV_4T_1R = read.csv(paste(extracted_color, "IF8_UV_4T_1R.csv", sep=""))

IF8_Vis_4T_1R = cbind(IF8_Vis_4T_1R, IF8_UV_4T_1R[,4:2479]) #specifying columns is to not include the column names in the UV file

# write.csv(IF8_Vis_4T_1R, paste(extracted_color, "IF8_Vis_UV_4T_1R.csv"), row.names = FALSE)



IF9_Vis_4T_1R = read.csv(paste(extracted_color, "IF9_Vis_4T_1R.csv", sep="")) #waiting on landmarking for IF9 UV and Vis to do this one

IF9_UV_4T_1R = read.csv(paste(extracted_color, "IF9_UV_4T_1R.csv", sep=""))

IF9_Vis_4T_1R = cbind(IF9_Vis_4T_1R, IF9_UV_4T_1R[,4:2479])

# write.csv(write.csv(IF9_Vis_4T_1R, paste(extracted_color, "IF9_Vis_UV_4T_1R.csv"), row.names = FALSE))


IF10_Vis_4T_1R = read.csv(paste(extracted_color, "IF10_Vis_4T_1R.csv", sep="")) #waiting on landmarking for IF10 UV and Vis to do this one

IF10_UV_4T_1R = read.csv(paste(extracted_color, "IF10_UV_4T_1R.csv", sep=""))

IF10_Vis_4T_1R = cbind(IF10_Vis_4T_1R, IF10_UV_4T_1R[,4:2479])

# write.csv(write.csv(IF10_Vis_4T_1R, paste(extracted_color, "IF10_Vis_UV_4T_1R.csv"), row.names = FALSE))



# write.csv(IF10_Vis_4T_1R, paste(extracted_color, "IF10_Vis_UV_4T_1R.csv"), row.names = FALSE)

# # Just read in the clean files
# IF10_Vis_4T_1R = read.csv("data/Iso-Y_phenotype_data-main/ IF10_Vis_UV_4T_1R.csv", sep="")
# IF6_Vis_4T_1R = read.csv("data/Iso-Y_phenotype_data-main/ IF6_Vis_UV_4T_1R.csv", sep="")
# IF8_Vis_4T_1R = read.csv("data/Iso-Y_phenotype_data-main/ IF8_Vis_UV_4T_1R.csv", sep="")
# IF9_Vis_4T_1R = read.csv("data/Iso-Y_phenotype_data-main/IF9_Vis_UV_4T_1R.csv", sep="")


#unwarp all populations to the same, global consensus so that there is only 1 set of coords for below - need to update the Sampling_locations.csv file to be global consensus once I know how to make that file


## X Y coords of sampling points, saved from image processing/sample collection custom code
## The X Y coords from the image processing may need to have the y axis flipped. If the fish shown in the plots at the end of this code are upside-down, just multiply your y coords by -1
X_Y_coords = read.csv("extracted_colour//Sampling_locations.csv")

IF6_name = rep("Iso-Y6", times = 41)
IF8_name = rep("Iso-Y8", times = 48)
IF9_name = rep("Iso-Y9", times = 42)
IF10_name = rep("Iso-Y10", times = 42)


## Combine the above vectors into one vector. This vector is in the appropriate order so that the population label aligns with the dataframe rows
#dapc.all_0p_4DT.pops = c(IF6_name, IF8_name, IF10_name)
dapc.all_0p_4DT.pops = c(IF6_name, IF8_name, IF9_name, IF10_name)

dapc.all_0p_4DT = rbind(IF6_Vis_4T_1R, IF8_Vis_4T_1R, IF9_Vis_4T_1R, IF10_Vis_4T_1R)


dapc_4DT.2_noPRLP17 = dapc(dapc.all_0p_4DT[,c(4:9907)], dapc.all_0p_4DT.pops, var.contrib = TRUE, scale = TRUE, n.pca = 10, n.da = 10, var.loadings = TRUE) 
#n.da is the number of discriminant functions to use, which defaults to 1 less than number of groups


IF_line = c(rep(6,41),rep(8,48),rep(10,42),rep(10,42))

dapc_4DT.2_noPRLP17$ind.coord_lin_names = cbind(IF_line, dapc_4DT.2_noPRLP17$ind.coord)




fish_dapc_loading_heatmap = function(dapc_object, axis, filename, ylim, X_Y_coords, color.channel = "sum",plot.legend=T){
  
  
  variable.loadings <- dapc_object$var.load[,axis]
  
  
  variable.names <- names(dapc_object$var.load[,axis])
  
  
  
  variable.loadings.data.frame = as.data.frame(variable.loadings)
  variable.names.data.frame = as.data.frame(variable.names)
  
  
  loadings.and.names = cbind(variable.loadings.data.frame, variable.names.data.frame)
  
  
  colnames(loadings.and.names) = c("corrs", "point")
  
  
  red.loadings.and.names = loadings.and.names[1:2476, ] #the number of rows in each needs to match number of rows in X_Y_coords - if I want to use thresholds >0 above, need to find a way to workaround this issue
  green.loadings.and.names = loadings.and.names[2477:4952, ]
  blue.loadings.and.names = loadings.and.names[4953:7428, ]
  UV.loadings.and.names = loadings.and.names[7429:9904, ]
  
  all.color.loadings = as.data.frame(cbind(red.loadings.and.names$corrs, green.loadings.and.names$corrs, blue.loadings.and.names$corrs, UV.loadings.and.names$corrs))
  colnames(all.color.loadings) = c("red", "gre", "blu", "UV")
  all.color.loadings[1:5, ]
  
  
  coordinates.and.loadings = cbind(X_Y_coords, all.color.loadings)
  
  coordinates.and.loadings$mean = rowMeans(coordinates.and.loadings[,c("red","gre","blu", "UV")])
  
  if(color.channel == "mean"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings[7])
    low.color = "white"
    high.color = "black"
  } else if(color.channel == "red"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$red)
    low.color = "black"
    high.color = "red"
  }else if(color.channel == "green"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$gre)
    low.color = "black"
    high.color = "green"
  }else if(color.channel == "blue"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$blu)
    low.color = "black"
    high.color = "blue"
  }else {   #default is UV
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$UV)
    low.color = "black"
    high.color = "white"
  }
  names(all.color.loadings) = c("V1","V2", "response")
  
  if(plot.legend){
    ggplot(all.color.loadings, aes(V1,-V2, color = response)) +
      ggtitle("") +
      geom_point(size=1.5) +
      coord_fixed() +
      scale_color_gradient2(low = low.color, mid = "grey", high = high.color,breaks = c(-0.006,0,0.006),labels = c(-0.006,0,0.006), limits = ylim, name = toupper(color.channel)) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            plot.title = element_text(size = 24, family = "Avenir", hjust = 0.5, vjust=5),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),
            legend.position = "right",
            legend.text = element_text(size=10, family = "Avenir"),
            legend.title = element_text(size=12, family = "Avenir"))
    
  } else {
    ggplot(all.color.loadings, aes(V1,-V2, color = response)) +
      ggtitle("") +
      geom_point(size=1.5) +
      coord_fixed() +
      scale_color_gradient2(low = low.color, mid = "grey", high = high.color, limits = ylim, name = toupper(color.channel)) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            plot.title = element_text(size = 24, family = "Avenir", hjust = 0.5, vjust=5),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),
            legend.position = "none") 
  }
  
}

# List colours 
cols <- c("red","green","blue","UV")
dapc1_fish <- lapply(cols,function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_4DT.2_noPRLP17, axis = 1, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = X_Y_coords, color.channel = colour,plot.legend = F)
})
dapc1_fish[[1]] <- dapc1_fish[[1]]+ggtitle("DF1")
dapc1_fish_fig <- cowplot::plot_grid(plotlist = dapc1_fish,ncol=1)

dapc2_fish <- lapply(cols,function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_4DT.2_noPRLP17, axis = 2, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = X_Y_coords, color.channel = colour,plot.legend = T)
})
dapc2_fish[[1]] <- dapc2_fish[[1]]+ggtitle("DF2")
dapc2_fish_fig <- cowplot::plot_grid(plotlist = dapc2_fish,ncol=1)

dapc22_fish <- lapply(cols,function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_4DT.2_noPRLP17, axis = 2, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = X_Y_coords, color.channel = colour,plot.legend = F)
})
dapc22_fish[[1]] <- dapc22_fish[[1]]+ggtitle("DF2")
dapc22_fish_fig <- cowplot::plot_grid(plotlist = dapc22_fish,ncol=1)

dapc3_fish <- lapply(cols,function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_4DT.2_noPRLP17, axis = 3, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = X_Y_coords, color.channel = colour,plot.legend = T)
})
dapc3_fish[[1]] <- dapc3_fish[[1]]+ggtitle("DF3")
dapc3_fish_fig <- cowplot::plot_grid(plotlist = dapc3_fish,ncol=1)

# Read in DAPC figure and plot alongside...
# dapc_fig <- readRDS("data/complete_dapc_fig.rds")

###### MAKE DAPC FIGURE HERE ########
library(adegenet)
library(ggplot2)
DAPC_data <- readRDS("DAPC_results.rds")

# Plot with ggplot
plot_dd <- data.frame(DAPC_data$ind.coord)
plot_dd$IF_line <- DAPC_data$grp
plot_dd$IF_line <- gsub('\\s+', '', plot_dd$IF_line)

# Get all the label info...
cols <- data.frame(line=paste0("Iso-Y",c(10,6,8,9)),
                   col=c("#3CBB75","#B22222","#2D708E","#481567"))

plot_labs <- data.frame(DAPC_data$grp.coord)
plot_labs$IF_line <- rownames(plot_labs)
plot_labs$IF_line <- c("Iso-Y10", "Iso-Y6", "Iso-Y8", "Iso-Y9")

# And get the eigenvals
eigenvals <- round((DAPC_data$eig/sum(DAPC_data$eig)) * 100,2)


dapc_fig <- ggplot(plot_dd,aes(LD1,LD2,colour=IF_line))+
  theme_classic()+
  geom_vline(xintercept = 0,alpha=0.3)+
  geom_hline(yintercept = 0,alpha=0.3)+
  geom_point(show.legend = F,size=3, alpha=0.7)+
  scale_colour_manual(breaks = cols$line,
                      values = cols$col)+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y6",],aes(x=LD1,xend=plot_labs[plot_labs$IF_line == "Iso-Y6","LD1"],
                                                              y=LD2,yend=plot_labs[plot_labs$IF_line == "Iso-Y6","LD2"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y8",],aes(x=LD1,xend=plot_labs[plot_labs$IF_line == "Iso-Y8","LD1"],
                                                              y=LD2,yend=plot_labs[plot_labs$IF_line == "Iso-Y8","LD2"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y9",],aes(x=LD1,xend=plot_labs[plot_labs$IF_line == "Iso-Y9","LD1"],
                                                              y=LD2,yend=plot_labs[plot_labs$IF_line == "Iso-Y9","LD2"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y10",],aes(x=LD1,xend=plot_labs[plot_labs$IF_line == "Iso-Y10","LD1"],
                                                               y=LD2,yend=plot_labs[plot_labs$IF_line == "Iso-Y10","LD2"]))+
  
  geom_label(data=plot_labs,aes(x=LD1,y=LD2,label=IF_line),show.legend = F,size=8, family = "Avenir", alpha=0.8)+
  theme(axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, famil = "Avenir"),
        legend.position = "none")+
  scale_x_continuous(breaks=c(-6,0,6))+
  scale_y_continuous(breaks=c(-6,0,6))+
  labs(x=paste0("Discriminant Function 1 (",eigenvals[1],"%)"),
       y=paste0("Discriminant Function 2 (",eigenvals[2],"%)"))+
  stat_ellipse()

png("figs/Figure1_full_dapc_only.png",width = 40, height = 20, units = "cm", res = 400)
# plot_grid(dapc_fig,dapc1_fish_fig,dapc2_fish_fig,ncol = 3,rel_widths = c(2.5,0.75,1.05),labels=c("a","b"),label_size = 20)
plot_grid(dapc_fig,ncol = 1,rel_widths = c(2.5))
dev.off()

## Make figure for DAPC Axes 2 & 3
dapc23_fig <- ggplot(plot_dd,aes(LD2,LD3,colour=IF_line))+
  theme_classic()+
  geom_vline(xintercept = 0,alpha=0.5)+
  geom_hline(yintercept = 0,alpha=0.5)+
  geom_point(show.legend = F,size=2)+
  scale_colour_manual(breaks = cols$line,
                      values = cols$col)+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y6",],aes(x=LD2,xend=plot_labs[plot_labs$IF_line == "Iso-Y6","LD2"],
                                                              y=LD3,yend=plot_labs[plot_labs$IF_line == "Iso-Y6","LD3"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y8",],aes(x=LD2,xend=plot_labs[plot_labs$IF_line == "Iso-Y8","LD2"],
                                                              y=LD3,yend=plot_labs[plot_labs$IF_line == "Iso-Y8","LD3"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y9",],aes(x=LD2,xend=plot_labs[plot_labs$IF_line == "Iso-Y9","LD2"],
                                                              y=LD3,yend=plot_labs[plot_labs$IF_line == "Iso-Y9","LD3"]))+
  geom_segment(data=plot_dd[plot_dd$IF_line == "Iso-Y10",],aes(x=LD2,xend=plot_labs[plot_labs$IF_line == "Iso-Y10","LD2"],
                                                               y=LD3,yend=plot_labs[plot_labs$IF_line == "Iso-Y10","LD3"]))+
  
  geom_label(data=plot_labs,aes(x=LD2,y=LD3,label=IF_line),show.legend = F,size=8,alpha=0.8, family = "Avenir")+
  theme(axis.title = element_text(size=20, family = "Avenir"),
        axis.text = element_text(size=18, family = "Avenir"),
        legend.position = "none")+
  labs(x=paste0("Discriminant Function 2 (",eigenvals[2],"%)"),
       y=paste0("Discriminant Function 3 (",eigenvals[3],"%)"))+
  stat_ellipse()




#####################################
empty <- ggplot()+theme_void()

# Compose
library(cowplot)
# pdf("figs/full_dapc_heatmaps.pdf",width = 14,height=10)
# plot_grid(plot_grid(dapc2_fish_fig,dapc_fig,ncol=2,rel_widths = c(1,3)),
#           plot_grid(empty,dapc1_fish_fig,ncol=2,rel_widths = c(1,3)),
#           nrow = 2,ncol = 1,rel_heights = c(3,1))
# dev.off()

png("figs/Figure1_full_dapc_heatmaps.png",width = 40, height = 20, units = "cm", res = 400)
# plot_grid(dapc_fig,dapc1_fish_fig,dapc2_fish_fig,ncol = 3,rel_widths = c(2.5,0.75,1.05),labels=c("a","b"),label_size = 20)
plot_grid(dapc_fig,dapc1_fish_fig,dapc2_fish_fig,ncol = 3,rel_widths = c(2.5,0.75,1.05))
dev.off()

pdf("figs/FigureSX_full_dapc_DF2_DF3_heatmaps.pdf",width = 14,height=8)
plot_grid(dapc23_fig,dapc22_fish_fig,dapc3_fish_fig,ncol = 3,rel_widths = c(2.5,0.75,1.15),labels=c("A","B"),label_size = 20)
dev.off()




