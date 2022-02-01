# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/recomb_fst_correlations/")

# Analysis of poolseq AFs
lib <- c("data.table","ggplot2","vcfR","tidyverse","ggridges","viridis","pbmcapply")
sapply(lib,library,character.only=T)

# # Read in poolvcf 
# pool_vcf <- read.vcfR("~/Exeter/VCFs/poolseq_variant_filtered.vcf.gz")

# # Go through each chromosome and build up allele frequencies...
# ADtoAF <- function(AD_vector,haploid_N=FALSE){
#   
#   af_dd <- data.frame(AD_vector) %>% separate("AD_vector",sep=",",into=c("ref","alt"))
#   af_dd$ref <- as.integer(af_dd$ref)
#   af_dd$alt <- as.integer(af_dd$alt)
#   af_dd$AF <- af_dd$ref / rowSums(af_dd)
#   
#   if(haploid_N){
#     return(rowSums(af_dd[,c("ref","alt")]))
#   } else {
#     return(af_dd$AF)
#   }
# }
# 
# if(!(file.exists("~/Exeter/josie_poolseq/outputs/all_AF_change.rds"))){
#   chr_AFs <- lapply(paste0("chr",1:23),function(x){
#     print(x)
#     
#     # Subset chrom and read it in
#     # vcf_tempfile <- tempfile(pattern = "gt_plot", fileext = '.vcf')
#     # on.exit({ unlink(vcf_tempfile) })  
#     # system(paste0("bcftools view -r ",x," ~/Exeter/VCFs/poolseq_variant_filtered.vcf.gz > ",vcf_tempfile))
#     
#     # Fetch the genotypes
#     pool_vcf <- read.vcfR(paste0("~/Exeter/VCFs/poolseq_variant_filtered_",x,".vcf.gz"))
#     chr_geno <- extract.gt(pool_vcf[pool_vcf@fix[,1] == x,],element = "AD")
#     
#     # Convert to AF
#     chr_AFs <- apply(chr_geno,2,ADtoAF)
#     # chr_haploids <-  apply(chr_geno,2,ADtoAF,haploid_N=TRUE)
#     
#     # Calculate pairwise AF change among all pools...
#     chr_AF_change <- matrix(ncol=6,nrow=nrow(chr_AFs))
#     comps_to_make <- combn(colnames(chr_AFs),2)
#     colnames(chr_AF_change) <- sapply(1:ncol(comps_to_make),function(i) paste0(comps_to_make[1,i],"_",comps_to_make[2,i]))
#     for(i in 1:ncol(comps_to_make)){
#       chr_AF_change[,i] <- chr_AFs[,comps_to_make[1,i]] - chr_AFs[,comps_to_make[2,i]]
#     }
#     chr_AF_change <- data.frame(chr_AF_change)
#     chr_AF_change$chr <- x
#     
#     # Return
#     return(list(chr_AF_change,chr_haploids))
#   })
#   saveRDS(chr_AFs,"~/Exeter/josie_poolseq/outputs/all_AF_change.rds")
# } else {
#   chr_AFs <- readRDS("~/Exeter/josie_poolseq/outputs/all_AF_change.rds")  
# }
# 
# #chr_haploidN <- lapply(chr_AFs,'[[',2)
# chr_AF_change <- lapply(chr_AFs,'[[',1)
# 
# # Condense to a df and plot
# #chr_haploidN <- data.frame(rbindlist(lapply(chr_haploidN,function(x) return(data.frame(x)))))
# AF_dd <- data.frame(rbindlist(chr_AF_change))
# AF_dd_melt <- melt(AF_dd)
# colnames(AF_dd_melt) <- c("chr","comparison","AF_change")
# 
# # Plot ridges...
# AF_dd_melt$chr_F <- factor(AF_dd_melt$chr,levels = paste0("chr",23:1))
# ggplot(AF_dd_melt,aes(x=abs(AF_change),y=chr_F))+
#   stat_density_ridges()+
#   facet_wrap(~comparison,nrow=1)
# 
# # Now report mean, median and sd for all chroms...
# AF_dd_melt$AFD <- abs(AF_dd_melt$AF_change)
# avg_AF <- data.frame(AF_dd_melt %>% group_by(chr,comparison) %>% summarise(mean_AF = mean(AFD,na.rm = T),
#                                                                            median_AF = median(AFD,na.rm = T),
#                                                                            sd_AF = sd(AFD,na.rm = T)))
# 
# avg_AF$chr_F <- factor(avg_AF$chr,levels=paste0("chr",1:23))
# ggplot(avg_AF,aes(fill=comparison,x=chr_F,y=median_AF))+
#   geom_bar(stat = "identity",position="dodge")

# Read in the Fst and show distributions for each chromosome --------------
all_fst <- data.frame(fread("../poolfstat/outputs/final_pairwise_fst.tsv"))

# Transform to long-form and plot as ridges...
all_fst_melt <- melt(all_fst[,c("chr","pos",grep("IF",colnames(all_fst),value = T))],
                     id.vars=c("chr","pos"))

all_fst_melt <- all_fst_melt[grep("chr",all_fst_melt$chr),]
all_fst_melt$chr_F <- factor(all_fst_melt$chr,levels=paste0("chr",1:23))
ggplot(all_fst_melt,aes(x=value,y=chr_F))+
  stat_density_ridges(quantile_lines = 2)+
  facet_wrap(~variable,nrow=1)

chr_avgs <- data.frame(all_fst_melt %>% group_by(chr,variable) %>% dplyr::summarise(avg_Fst=median(value,na.rm = T)))
chr_avgs$chr_F <- factor(chr_avgs$chr,paste0("chr",1:23))
ggplot(chr_avgs,aes(fill=variable,x=chr_F,y=avg_Fst))+
  geom_bar(stat = "identity",position="dodge")


# Get our recombination map -----------------------------------------------
rec_maps <- lapply(paste0("../recomb_maps/chr",1:23,"_smoothed.map"),read.table,header=T)
chrom_lengths <- read.table("~/Dropbox/Sussex_Guppies/genomes/STAR.chromosomes.release.fasta.fai")[1:23,1:2]

# What is average recombination per chromosomes
chr_rec_avg <- data.frame(chr=paste0("chr",1:23),
                          avg_rec=sapply(1:23,function(chr){
  rec_tmp <- rec_maps[[chr]]
  colnames(rec_tmp)[1] <- "start"
  rec_tmp$end <- c(rec_tmp$start[2:nrow(rec_tmp)],chrom_lengths[chr,2])
  rec_tmp$coverage <- rec_tmp$end - rec_tmp$start
  weighted.mean(rec_tmp$COMBINED_rate,rec_tmp$coverage)
}))

chr_rec_avg_fst <- merge(chr_avgs,chr_rec_avg,by = "chr")
ggplot(chr_rec_avg_fst,aes(x=avg_rec,y=avg_Fst,colour=variable))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Mean Rec Rate (per chrom)",y="Mean Fst")


# chr-level fst-rec correlations
chr_rec_fst <- data.frame(rbindlist(pbmclapply(1:23,function(chr){
  
  rec_tmp <- rec_maps[[chr]]
  colnames(rec_tmp)[1] <- "start"
  rec_tmp$end <- c(rec_tmp$start[2:nrow(rec_tmp)],chrom_lengths[chr,2])
  
  # Windowise rec_rate
  rec_winds <- seq(0,max(rec_tmp$end),100000)
  rec_avgs <- sapply(1:length(rec_winds),function(x){
    rec_sub <- rec_tmp[rec_tmp$start < rec_winds[x]+100000 &
                         rec_tmp$end > rec_winds[x],]
    if(nrow(rec_sub) > 0){
      rec_sub[rec_sub$end > rec_winds[x]+100000,"end"] <- rec_winds[x]+100000
      rec_sub[rec_sub$start < rec_winds[x],"start"] <- rec_winds[x]
      
      rec_sub$window_coverage <- rec_sub$end-rec_sub$start
      return(weighted.mean(rec_sub$COMBINED_rate,weights=rec_sub$window_coverage))
    } else {
      return(NA)
    }
  })
  rec_avg_winds <- data.frame(mid=rec_winds+50000,
                              rec=rec_avgs)
  rec_avg_winds$window <- cut_interval(rec_avg_winds$mid,length=100000)
  
  # Fetch the fst
  fst_tmp <- all_fst_melt[all_fst_melt$chr == paste0("chr",chr),]

  # Windowise fst
  fst_tmp$window <- cut_interval(fst_tmp$pos,length=100000)
  
  # Get avg fst
  mean_fst_windows <- data.frame(fst_tmp %>% group_by(window,variable) %>% dplyr::summarise(mean_fst=mean(value)))
  
  # Merge
  mean_fst_windows <- merge(mean_fst_windows,rec_avg_winds[,c("window","rec")])
  
  comp_corrs <- sapply(unique(mean_fst_windows$variable),function(comp){
    cor_res <- cor.test(mean_fst_windows[mean_fst_windows$variable == comp,"mean_fst"],
             mean_fst_windows[mean_fst_windows$variable == comp,"rec"],method = "kendall")
    cor_res$estimate
  })
  
  comp_corrs2 <- sapply(unique(mean_fst_windows$variable),function(comp){
    cor_res <- cor.test(mean_fst_windows[mean_fst_windows$variable == comp,"mean_fst"],
                        mean_fst_windows[mean_fst_windows$variable == comp,"rec"],method = "kendall")
    cor_res$estimate
  })
  
  data.frame(comps=unique(mean_fst_windows$variable),
             rec_corr=comp_corrs,
             index_corr=comp_corrs2,
             chr=paste0("chr",chr))
},mc.cores=4)))

# Plot fst-rec heats
chr_rec_fst$comps <- gsub("IF", "Iso-Y",chr_rec_fst$comps)
chr_rec_fst$comps <- gsub("_vs_", " vs ",chr_rec_fst$comps)

chr_rec_fst$chr <- sub("chr","LG",chr_rec_fst$chr)
chr_rec_fst$chr_F <- factor(chr_rec_fst$chr,levels=paste0("LG",23:1))

rec_corr_heats <- ggplot(chr_rec_fst,aes(y=chr_F,x=comps,fill=rec_corr))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle=45,hjust=1,family="Avenir"),
        axis.text.y = element_text(size=14, family="Avenir"),
        axis.title = element_blank())
  #labs(fill=expression('F'[ST]*'-Recombination Correlation'))
  #labs(fill="Fst-Recombination\nCorrelation")

# Get chromosome averages...
# chrom_rec_avgs <- data.frame(chr_rec_fst %>% group_by(chr) %>% summarise(mean_corr=mean(rec_corr)))
chrom_rec_avgs <- data.frame(chr_rec_fst %>% group_by(chr) %>% dplyr::summarise(median_corr=median(rec_corr)))
chrom_rec_avgs <- chrom_rec_avgs[order(chrom_rec_avgs$median_corr),]
chr_rec_fst$chr_F2 <- factor(chr_rec_fst$chr,levels=chrom_rec_avgs$chr)


chr_rec_corr_sum <- ggplot(chr_rec_fst,aes(y=chr_F2,x=rec_corr))+
  geom_boxplot(fill="skyblue",alpha=0.4)+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=12,family="Avenir"),
        axis.title.x = element_text(size=14,family="Avenir"))+
  labs(x=expression('F'[ST]*'-Recombination Correlation'))

# Plot histogram...
rec_corr_hist <- ggplot(chr_rec_fst,aes(x=rec_corr))+
  geom_histogram()+
  theme_bw()+
  labs(y="Counts",x=expression('F'[ST]*'-Recombination Correlation'))+
 # scale_y_continuous(labels=seq(0,15,5),breaks=seq(0,15,5))+
  theme(axis.text = element_text(size=14,family="Avenir"),
        axis.title = element_text(size=16,family="Avenir"))

library(cowplot)
p1 <- plot_grid(plot_grid(rec_corr_heats,chr_rec_corr_sum,
          ncol=2,rel_widths=c(2,2)),
          rec_corr_hist,
          ncol=1,rel_heights = c(2,1))


ggsave("figs/FST_Recombination_corrs.png",p1,width=20,height=20,units="cm",dpi=400)


# Fetch the lostruct results and compile  ---------------------------------
lostruct_res <- readRDS("data_localPCA_results_windsize100.rds")
start_end <- lostruct_res[[1]]
mds_scores <- lostruct_res[[2]]

mds_res <- data.frame(rbindlist(lapply(1:23,function(x){
  data.frame(start_end[[x]][1:nrow(mds_scores[[x]]$points),],
             mds1=mds_scores[[x]]$points[,1])
})))

# Calculate FST-PC1
fst_pca <- prcomp(na.omit(all_fst[,grep("IF",colnames(all_fst),value=T)]),scale=T,center = T)
clean_all_fst <- na.omit(all_fst)
fst_pc1_snp_scores <- data.frame(chr=clean_all_fst$chr,
                                 pos=clean_all_fst$pos,
                                 PC1=fst_pca$x[,1])

# Merge these together
# Set up lostruct as genomic regions
mds_res_bed <- data.table(mds_res[,c("chr","start","end")])
setkey(mds_res_bed, chr, start, end)

# Set these fst as genomic regions
fst_pc1_snp_scores$end <- fst_pc1_snp_scores$pos+1
fst_pc1_snp_scores_bed <- data.table(fst_pc1_snp_scores)
setkey(fst_pc1_snp_scores_bed, chr, pos, end)

# Find overlaps
fst_lostruct_overlap <- data.frame(foverlaps(fst_pc1_snp_scores_bed, mds_res_bed, minoverlap=1L,nomatch=NA))

# Use overlaps to group genes
fst_lostruct_overlap$lostruct_window <- paste0(fst_lostruct_overlap$chr,":",fst_lostruct_overlap$start,"-",fst_lostruct_overlap$end)
fst_lostruct_overlap <- na.omit(fst_lostruct_overlap)

# Merge back with lostruct
mds_res$lostruct_window <- paste0(mds_res$chr,":",mds_res$start,"-",mds_res$end)
mds_res_merge <- merge(mds_res,fst_lostruct_overlap[,c("PC1","lostruct_window")],by="lostruct_window")

# Fetch average PC1 per window...
fst_mds1_windows <- data.frame(mds_res_merge %>% group_by(lostruct_window) %>% summarise(mean_fst_pc1=mean(PC1),
                                                                                         mds1=mean(mds1),
                                                                                         chr=unique(chr),
                                                                                         start=unique(start),
                                                                                         end=unique(end)))

# Highlight the LG1 co-ords
fst_mds1_windows$LG1_outlier <- "No"
fst_mds1_windows[fst_mds1_windows$chr == "chr1" &
                   fst_mds1_windows$start > 11500000 &
                   fst_mds1_windows$end < 15800000,"LG1_outlier"] <- "LG1_outlier"

fst_mds1_windows2 <- fst_mds1_windows %>% mutate(LG1_outlier = ifelse(chr == "chr3" && fst_mds1_windows$start > )) 

# And plot
ggplot(fst_mds1_windows,aes(x=abs(mds1),y=mean_fst_pc1,colour=LG1_outlier))+
  geom_point(alpha=0.5)+
  scale_colour_manual(breaks=c("Yes","No"),values = c("red2","black"))+
  theme_bw()+
  theme(axis.text = element_text(size=14, family="Avenir"),
        axis.title = element_text(size=16, family="Avenir"),
        legend.position="top")+
  labs(colour="LG1 Region 2",x="abs(MDS1)",y=expression('Mean Z-F'[ST]*' PC1'))+
  geom_smooth(method="lm",colour="blue2")
 

expression('F'[ST]*'-Recombination Correlation')