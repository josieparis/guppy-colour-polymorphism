# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-c("lemon", "grid", "gridExtra", "ggplot2", "cowplot", "factoextra", "readr", "tidyr", "stringr", "data.table")
lapply(lib,library,character.only=T)

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/phenotype_data_current/extracted_colour/")

PC.for.ADONIS <- read.csv(file = "PC_for_ADONIS.csv", header = TRUE, sep = ",")

#adding line labels
## create a vector with the line ID; the number is how many males in that group
IF6_name = rep("Iso-Y 6", times = 41)
IF8_name = rep("Iso-Y 8", times = 48)
IF9_name = rep("Iso-Y 9", times = 42)
IF10_name = rep("Iso-Y 10", times = 42)

Line = factor(c(IF6_name, IF8_name, IF9_name, IF10_name), levels = c("Iso-Y 6", "Iso-Y 8", "Iso-Y 9", "Iso-Y 10") )

PC.for.ADONIS <- cbind(Line, PC.for.ADONIS)


element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

scaleFUN <- function(x) sprintf("%.3f", x)

myCol <- c("#B22222", "#2D708E", "#481567", "#3CBB75") # 6,8,9,10

c("#3CBB75", "#B22222", "#2D708E", "#481567")


pc_plotter = function(dataset, column, title, ylabel = " ", legend.position = "none") {
  
  
  
  a <- ggplot(dataset, aes(x = dataset[,column]))
  a = a + geom_density(aes(fill = Line), alpha = 0.4) +
    labs(title = title) +
    ylab(ylabel) + 
    theme(plot.title = element_textbox(
      family="Avenir", size = 24, hjust = 0.5, margin = margin(t = 5, b = 5)
    ), axis.title.x=element_blank(), axis.text.x = element_text(family="Avenir", size = 20), axis.text.y = element_text(family="Avenir", size = 24), axis.title.y = element_text(family="Avenir", size = 24), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = legend.position) +
    scale_y_continuous(labels=scaleFUN) +
    scale_fill_manual( values = myCol)
  a
}



a = pc_plotter(PC.for.ADONIS, 3, "PC1", legend.position = "top")

b = pc_plotter(PC.for.ADONIS, 4, "PC2")

c = pc_plotter(PC.for.ADONIS, 5, "PC3")

d = pc_plotter(PC.for.ADONIS, 6, "PC4")

e = pc_plotter(PC.for.ADONIS, 7, "PC5")

f = pc_plotter(PC.for.ADONIS, 8, "PC6")

g = pc_plotter(PC.for.ADONIS, 9, "PC7")

h = pc_plotter(PC.for.ADONIS, 10, "PC8")

i = pc_plotter(PC.for.ADONIS, 11, "PC9")

j = pc_plotter(PC.for.ADONIS, 12, "PC10")

k = pc_plotter(PC.for.ADONIS, 13, "PC11")

l = pc_plotter(PC.for.ADONIS, 14, "PC12")

m = pc_plotter(PC.for.ADONIS, 15, "PC13")

n = pc_plotter(PC.for.ADONIS, 16, "PC14")

o = pc_plotter(PC.for.ADONIS, 17, "PC15")

p = pc_plotter(PC.for.ADONIS, 18, "PC16")

q = pc_plotter(PC.for.ADONIS, 19, "PC17")




pp <- grid_arrange_shared_legend(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q, nrow = 6, ncol = 3, position = c("top")) #ncol determines number of columns
#nrow determines number of rows


ggsave("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/Figures/phenotype_PCs_SIFigure1.png", pp, width = 40, height = 40, units = "cm", dpi = 400)

