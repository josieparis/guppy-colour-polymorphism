rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-c("dplyr","ggplot2", "devtools", "magrittr", "gridExtra", "readr", "tidyr", "stringr", "data.table")
lapply(lib,library,character.only=T)

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/phenotype_data_current/")

D_all_lines1 <- read.csv(file = "./extracted_colour/IF_lines_trace.csv")
extracted_color = "extracted_colour/"

#bind the UV csv files with the Vis ones and then revise code belows to analyze 4 channels instead of 3

IF6_Vis_4T_1R = read.csv(paste(extracted_color, "IF6_Vis_4T_1R.csv", sep=""))

IF6_UV_4T_1R = read.csv(paste(extracted_color, "IF6_UV_4T_1R.csv", sep=""))

IF6_Vis_4T_1R = cbind(IF6_Vis_4T_1R, IF6_UV_4T_1R[,4:ncol(IF6_UV_4T_1R)])

IF8_Vis_4T_1R = read.csv(paste(extracted_color, "IF8_Vis_4T_1R.csv", sep=""))

IF8_UV_4T_1R = read.csv(paste(extracted_color, "IF8_UV_4T_1R.csv", sep=""))

IF8_Vis_4T_1R = cbind(IF8_Vis_4T_1R, IF8_UV_4T_1R[,4:ncol(IF8_UV_4T_1R)]) #specifying columns is to not include the column names in the UV file


IF9_Vis_4T_1R = read.csv(paste(extracted_color, "IF9_Vis_4T_1R.csv", sep="")) #waiting on landmarking for IF9 UV and Vis to do this one

IF9_UV_4T_1R = read.csv(paste(extracted_color, "IF9_UV_4T_1R.csv", sep=""))

IF9_Vis_4T_1R = cbind(IF9_Vis_4T_1R, IF9_UV_4T_1R[,4:ncol(IF9_UV_4T_1R)])

IF10_Vis_4T_1R = read.csv(paste(extracted_color, "IF10_Vis_4T_1R.csv", sep=""))

IF10_UV_4T_1R = read.csv(paste(extracted_color, "IF10_UV_4T_1R.csv", sep=""))

IF10_Vis_4T_1R = cbind(IF10_Vis_4T_1R, IF10_UV_4T_1R[,4:ncol(IF10_UV_4T_1R)])

IF6_vars = apply(IF6_Vis_4T_1R[c(4:9904)], 2, var) 

#write.csv(data.frame(names(IF6_vars), IF6_vars), paste(extracted_color, "IF6_vars.csv"), row.names = FALSE, col.names =c("color channel and sampling location", "variance "))

write.table(data.frame(names(IF6_vars), IF6_vars), paste(extracted_color, "IF6_vars.csv"), sep = ",", row.names = FALSE, col.names =c("color channel and sampling location", "variance "))


IF8_vars = apply(IF8_Vis_4T_1R[c(4:9904)], 2, var) 

#write.csv(IF8_vars, paste(extracted_color, "IF8_vars.csv"), col.names = c("color channel and sampling location", "variance "))

write.table(data.frame(names(IF8_vars), IF8_vars), paste(extracted_color, "IF8_vars.csv"), sep = ",", row.names = FALSE, col.names =c("color channel and sampling location", "variance "))


IF9_vars = apply(IF9_Vis_4T_1R[c(4:9904)], 2, var) 

IF10_vars = apply(IF10_Vis_4T_1R[c(4:9904)], 2, var) 

#write.csv(IF10_vars, paste(extracted_color, "IF10_vars.csv"), col.names = c("color channel and sampling location", "variance "))

write.table(data.frame(names(IF10_vars), IF10_vars), paste(extracted_color, "IF10_vars.csv"), sep = ",", row.names = FALSE, col.names =c("color channel and sampling location", "variance "))

#red variances for each line
IF6_r_vars = sum(IF6_vars[1:2476])

IF8_r_vars = sum(IF8_vars[1:2476])

IF9_r_vars = sum(IF9_vars[1:2476])

IF10_r_vars = sum(IF10_vars[1:2476])


#blue variances for each line
IF6_blue_vars = sum(IF6_vars[2477:4952])

IF8_blue_vars = sum(IF8_vars[2477:4952])

IF9_blue_vars = sum(IF9_vars[2477:4952])

IF10_blue_vars = sum(IF10_vars[2477:4952])

#green variances for each line
IF6_green_vars = sum(IF6_vars[4953:7428])

IF8_green_vars = sum(IF8_vars[4953:7428])

IF9_green_vars = sum(IF9_vars[4953:7428])

IF10_green_vars = sum(IF10_vars[4953:7428])

#UV variances for each line
IF6_UV_vars = sum(IF6_vars[7429:9904])

IF8_UV_vars = sum(IF8_vars[7429:9904])

IF9_UV_vars = sum(IF9_vars[7429:9904])

IF10_UV_vars = sum(IF10_vars[7429:9904])

#total variance for each line
IF6_total_vars = sum(IF6_vars)

IF8_total_vars = sum(IF8_vars)

IF9_total_vars = sum(IF9_vars)

IF10_total_vars = sum(IF10_vars)


IF9_IF10_diff = IF9_total_vars - IF10_total_vars

IF8_IF10_diff = IF8_total_vars - IF10_total_vars

IF6_IF10_diff = IF6_total_vars - IF10_total_vars

IF8_IF9_diff = IF8_total_vars - IF9_total_vars

IF6_IF9_diff = IF6_total_vars - IF9_total_vars

IF6_IF8_diff = IF6_total_vars - IF8_total_vars


## Load comparisons
IF9_IF10 = D_all_lines1[,7]
IF8_IF10 = D_all_lines1[,6]
IF8_IF9 = D_all_lines1[,5]
IF6_IF10 = D_all_lines1[,4]
IF6_IF9 = D_all_lines1[,3]
IF6_IF8 = D_all_lines1[,2]

## Create all the permutation plots

AR_PERM_plot6 = ggplot(D_all_lines1, aes(x = IF6_IF8)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y6 - Iso-Y8", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF6_IF8_diff, linetype="dashed", size=1.5) +
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5) 
  ggtitle("a")
AR_PERM_plot6


IF6_IF8_prop = sum(abs(IF6_IF8) >= abs(IF6_IF8_diff))  #count the number of traces that were larger than observed
IF6_IF8_P = IF6_IF8_prop/1000


AR_PERM_plot5 = ggplot(D_all_lines1, aes(x = IF6_IF9)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y6 - Iso-Y9", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF6_IF9_diff, linetype="dashed", size=1.5) +
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5)
  ggtitle("b")
AR_PERM_plot5



IF6_IF9_prop = sum(abs(IF6_IF9) >= abs(IF6_IF9_diff))  #count the number of traces that were larger than observed
IF6_IF9_P = IF6_IF9_prop/1000


AR_PERM_plot4 = ggplot(D_all_lines1, aes(x = IF8_IF9)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y8 - Iso-Y9", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF8_IF9_diff, linetype="dashed", size=1.5) + #need to make line show up
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5)
  ggtitle("c")
AR_PERM_plot4



IF8_IF9_prop = sum(abs(IF8_IF9) >= abs(IF8_IF9_diff))  #count the number of traces that were larger than observed
IF8_IF9_P = IF8_IF9_prop/1000


AR_PERM_plot2 = ggplot(D_all_lines1, aes(x = IF6_IF10)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y6 - Iso-Y10", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF6_IF10_diff, linetype="dashed", size=1.5) +
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5)
  ggtitle("d")
AR_PERM_plot2



IF6_IF10_prop = sum(abs(IF6_IF10) >= abs(IF6_IF10_diff))  #count the number of traces that were larger than observed
IF6_IF10_P = IF6_IF10_prop/1000


AR_PERM_plot = ggplot(D_all_lines1, aes(x = IF8_IF10)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y8 - Iso-Y10", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF8_IF10_diff, linetype="dashed", size=1.5) +
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5)
  ggtitle("e")
AR_PERM_plot



IF8_IF10_prop = sum(abs(IF8_IF10) >= abs(IF8_IF10_diff))  #count the number of traces that were larger than observed
IF8_IF10_P = IF8_IF10_prop/1000 #turns this into a p-value


AR_PERM_plot3 = ggplot(D_all_lines1, aes(x = IF9_IF10)) + 
  geom_density(fill="gray30", alpha = 0.8)+
  scale_x_continuous(name = "Iso-Y9 - Iso-Y10", limits=c(-16.5,16.5)) +
  scale_y_continuous(name = "")+
  #ggtitle("Proportion of permutations less than the observed")+
  theme_classic()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0,size = 12, family = "Avenir"),
        text=element_text(family="Avenir"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12)) +
  geom_vline(xintercept=IF9_IF10_diff, linetype="dashed", size=1.5) +
  #geom_rug(data=observ_diff_trace_ar, color = "red", lwd = 2.5)
  ggtitle("f")
AR_PERM_plot3



IF9_IF10_prop = sum(abs(IF9_IF10) >= abs(IF9_IF10_diff))  #count the number of traces that were larger than observed
IF9_IF10_P = IF9_IF10_prop/1000
#P = 0, indicating that Iso-Y 9 is significantly more variable than Iso-Y 10

## Put them all together
figure_grid <- grid.arrange(AR_PERM_plot6, AR_PERM_plot5, AR_PERM_plot4, AR_PERM_plot2, AR_PERM_plot, AR_PERM_plot3, ncol = 3, nrow = 3,
                            layout_matrix = rbind(c(1, NA, NA),
                                                  c(2, 3, NA),
                                                  c(4, 5, 6)
                                                  #C(AR_PERM_plot2, AR_PERM_plot, AR_PERM_plot3, NA))
                            ),
                            bottom = "Difference in phenotypic variance between lines",
                            left = "Density",
                            heights = unit(c(2,2,2),c("in","in","in")),
                            widths = unit(c(2,2,2),c("in","in","in")))

ggsave(figure_grid, file = "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/perm_tests_phenotypes.tiff",  dpi = 1200, height=160, width=200, units="mm") 
