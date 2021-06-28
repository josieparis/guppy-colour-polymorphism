pool_pi <- function(pool_vcf,wind_size=10000,n_cores=1,min_snp=10){
  # Add a chr error
  if(length(unique(pool_vcf@fix[,1])) > 1){
    stop("Too many chromosomes in VCF")
  }
  # Get all the window starts and end
  bp <- as.integer(pool_vcf@fix[,2])
  winds <- seq(0,max(bp),wind_size)
  winds2 <- winds+wind_size
  winds2[length(winds2)] <- max(bp)
  # Get the pool info
  pools <- colnames(pool_vcf@gt)[2:ncol(pool_vcf@gt)]
  # Work through windows
  winds_out <- data.frame(rbindlist(mclapply(1:length(winds),function(x){
    # Subset VCF
    vcf_sub <- pool_vcf[bp < winds2[x] & bp >= winds[x],]
    # Get the AD
    allele_counts <- extract.gt(vcf_sub,element="AD")
    # Get counts for each pool
    pool_pi <- sapply(pools,function(pool){
      if(nrow(vcf_sub@fix) < min_snp){
        return(NA)
      } else {
        # Take column
        tmp_counts <- allele_counts[,pool]
        # Separate 
        split_counts <- str_split(tmp_counts,",")
        # Merge to data.frame
        count_dd <- data.frame(count1 = as.integer(sapply(split_counts,function(y){return(y[1])})),
                               count2 = as.integer(sapply(split_counts,function(y){return(y[2])})))
        # Sum
        count_dd$sum <- rowSums(count_dd)
        # Get Pi
        count_dd$pi <- (2*(count_dd$count1 * count_dd$count2))/(count_dd$sum*(count_dd$sum-1))
        # Avg Pi
        pi_avg <- sum(count_dd$pi)/(winds2[x]-winds[x])
        return(pi_avg)
      }
    })
    # Out format
    out <- data.frame(pool = pools,
                      pi = pool_pi,
                      chr = unique(pool_vcf@fix[,1]),
                      start = as.integer(winds[x]),
                      end = as.integer(winds2[x]))
    rownames(out) <- NULL
    return(out)
  },mc.cores=n_cores)))
  return(winds_out)
}



# Usage eg:

#vcf_path <- "~/Dropbox/Sussex_Guppies/Analyses/pool-seq/vcf_files/chr12_variant_filtered_STAR_updated.vcf.gz"

#vcf<-read.vcfR(vcf_path)

#chr12_pi <- pool_pi(pool_vcf = vcf,wind_size=10000,n_cores=2,min_snp = 10)

# Libraries required!
#lib<-as.vector(c("vcfR","pbapply","PopGenome","cowplot","gridExtra","viridis","grid","data.table","ggplot2","parallel","dplyr","tidyr", "stringr"))
#lapply(lib,library,character.only=TRUE)