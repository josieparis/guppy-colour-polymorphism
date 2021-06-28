# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("~/Dropbox/Sussex_Guppies/Analyses/pool-seq/colour_gene_exploration/")

lib<-c("clusterProfiler", "enrichplot", "parallel","ggplot2", "cowplot", "ggrepel", "gridExtra", "biomaRt", "AnnotationHub")
lapply(lib,library,character.only=T)

## use ensembl annotations
ensembl=useMart("ensembl")

### list species available:
# ensembl_dd <- as.data.frame(listDatasets(ensembl))

## get guppy
guppy <- useDataset("preticulata_gene_ensembl",mart=ensembl)

## check it out if you fancy
# listAttributes(guppy)

## Grab the guppy universe
guppy_universe<-getBM(attributes = c("ensembl_gene_id","entrezgene_id","kegg_enzyme"),
                      mart=guppy)

guppy_entrez_universe<-as.character(unique(guppy_universe$entrezgene_id))


# Fetch the Uniprot IDs...
ensembl_digger<-function(outliers=NULL,genome1=NULL,genome2=NULL,aln_block_length=1000,aln_qual=40,
                         biomart_attributes=c("ensembl_gene_id","entrezgene_id","uniprot_gn_id", "description", "external_gene_name", "kegg_enzyme","chromosome_name","start_position","end_position", "go_id", "name_1006", "definition_1006"),n.cores=1){
  # Firstly, take all the outlier regions and make a multi-fasta
  lapply(outliers,function(x){
    system(paste0("samtools faidx ",genome1," ",x," >> outputs/tmp_multi.fa"))
  })
  # Align to old genome
  system2('/usr/local/bin/minimap2',
          args=c(paste0("-t ",n.cores),genome2,"outputs/tmp_multi.fa"),
          stdout ="outputs/tmp_multi.paf",wait=T)
  # Fetch regions
  aln<-read.table("outputs/tmp_multi.paf",fill=T)
  # Keep tidy!
  system(paste0("rm -f outputs/tmp_multi*"))
  # Only keep "high support alignment regions"
  aln<-aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]
  regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))
  # Pull uniprot genes from biomaRt for each region
  tmp_biomart<-getBM(attributes = biomart_attributes,
                     filters= "chromosomal_region",
                     values=regions,
                     mart=guppy)
  # Return
  return(tmp_biomart)
}

## Coordinates for analysis
## Region 1 genes: 4,079,988 - 5,984,584 bp
### Region 2 genes 9,627,619 - 17,074,870 bp
### Region 3 genes:  21,944,840 - 24,959,750 bp
## Region of interest: chr1:12210000-12220000
## LG12 MDS2 coordinates: 5,608,703 - 7,063,435 bp
## LG12 MDS1 coordinates: 21,913,633 - 25,664,533 bp

outliers<- "chr1"

male_genome <- "~/Dropbox/Sussex_Guppies/genomes/STAR.chromosomes.release.fasta"

female_genome <- "~/Dropbox/Sussex_Guppies/genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa"

LG1 <- ensembl_digger(outliers=outliers, genome1 = male_genome, genome2 = female_genome)

## write results:
write.table(file = "./outputs/LG12_MDS1_genes.tsv", LG12_MDS1, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## Read results in 
LG1_region1 <- read_tsv("outputs/LG1_region1_genes.tsv")
LG1_region2 <- read_tsv("outputs/LG1_region2_genes.tsv")
LG1_region3 <- read_tsv("outputs/LG1_region3_genes.tsv")


### First try KEGG enrichment:
LG1_region3_enrich<-enrichKEGG(gene=na.omit(unique(LG1_region3$entrezgene_id)),
                          keyType = 'kegg',
                          organism = 'pret',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "fdr",
                          universe = guppy_entrez_universe)

## Add to a table, and write results:
LG1_region1_enrich_results <- LG1_region1_enrich@result
write.table(file="./outputs/LG1_region1_KEGG_results.tsv", LG1_region1_enrich_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

LG1_region2_enrich_results <- LG1_region2_enrich@result
write.table(file="./outputs/LG1_region2_KEGG_results.tsv", LG1_region2_enrich_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

LG1_region3_enrich_results <- LG1_region3_enrich@result
write.table(file="./outputs/LG1_region3_KEGG_results.tsv", LG1_region3_enrich_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

# Try with all genes in the identified regions on LG1:

region1 <- LG1_region1$entrezgene_id
region2 <- LG1_region2$entrezgene_id
region3 <- LG1_region3$entrezgene_id
  
## append
LG1_regions <- c(region1,region2,region3)

LG1_regions_enrich<-enrichKEGG(gene=na.omit(unique(LG1_regions)),
                               keyType = 'kegg',
                               organism = 'pret',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "fdr",
                               universe = guppy_entrez_universe)

# Write results
LG1_regions_enrich_results <- LG1_regions_enrich@result
write.table(file="./outputs/LG1_regions_KEGG_results.tsv", LG1_regions_enrich_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

## Ok and now GO analysis ... 

library(AnnotationHub)

## Get the hub
hub <- AnnotationHub()

## Take a look for latest OrgDB
d <- display(hub) 

## on 14/06/2021, this was the latest hub for PRET
PRET <- hub[["AH86018"]]

## You can look at what info the annotation has (if you so desire):
# biomaRt::keytypes(PRET)

## Now run for each region
LG1_region3_enrich_GO <-enrichGO(gene=na.omit(unique(LG1_region3$entrezgene_id)),
                             keyType = 'ENTREZID',
                             OrgDb = PRET,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "fdr",
                             universe = guppy_entrez_universe)

## Add to a table, and write results:
LG1_region1_enrich_GO_results <- LG1_region1_enrich_GO@result
write.table(file="./outputs/LG1_region1_GO_results.tsv", LG1_region1_enrich_GO_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

LG1_region2_enrich_GO_results <- LG1_region2_enrich_GO@result
write.table(file="./outputs/LG1_region2_GO_results.tsv", LG1_region2_enrich_GO_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)

LG1_region3_enrich_GO_results <- LG1_region3_enrich_GO@result
write.table(file="./outputs/LG1_region3_GO_results.tsv", LG1_region3_enrich_GO_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)


## And finally run for all the regions 
LG1_regions_enrich_GO <-enrichGO(gene=na.omit(unique(LG1_regions)),
                                 keyType = 'ENTREZID',
                                 OrgDb = PRET,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "fdr",
                                 universe = guppy_entrez_universe)

LG1_regions_enrich_GO_results <- LG1_regions_enrich_GO@result
write.table(file="./outputs/LG1_regions_GO_results.tsv", LG1_regions_enrich_GO_results, sep="\t", quote=FALSE,col.names=TRUE,row.names=FALSE)



