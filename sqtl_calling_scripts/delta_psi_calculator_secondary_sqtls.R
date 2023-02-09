################################
# This script gets delta PSI values for all exons in a significant exon group
# (This is a fork of the script 'delta_psi_calculator_all_exons.R')
# 
# Jonah Einson
# 5/3/2021
# 
# USAGE: Rscript --vanilla delta_psi_calculator_secondary_sqtls.R <tissue_prefix> <chrom>
# 
# This script will read in the output of qtltools, select sGene - snp
# pairs, flip cases so 1 represents the higher included sQTL genotype, 
# then write delta-psi values. 
################################

source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/")
source("sQTL/QTL_plotter.R")

cn <- read_lines("sQTL/files/qtltools_nominal_header.txt") 

# Parse the command line for a tissue to use
if(interactive()){
  prefix = "Brain_Cortex"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  prefix = args[1]
  print(prefix)
}

completed_tissues <- read_lines("sQTL_updated/completed_tissues.txt")
individual_delta_psi_by_tiss <- list()

for(prefix in completed_tissues[8:9]){
  psi <- read.csv(sprintf("sQTL_updated/bed/%s_psi_v8_protocol.bed.gz", prefix),
                  sep = "\t", 
                  row.names = 4, 
                  check.names = F)
  psi <- psi[,-(1:3)]
  
  print(chrom)
  
  # Read in sQTLs and psi info
  
  # This is the set of secondary sQTLs that includes all exons with matching
  # signs to the top exon
  sQTLs <- read_tsv("sQTL_updated/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs.tsv")
  sQTLs <- filter(sQTLs, tiss == prefix)
  
  
  # Use the homemade get genotype function to read the VCF
  VCF = "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz" 
  
  gt <- c("0/0", "0/1", "1/1")
  n <- nrow(sQTLs)
  
  out <- matrix(nrow = n, ncol = 3)
  out_n <- matrix(nrow = n, ncol = 3)
  
  colnames(out) <- gt
  colnames(out_n) <- gt
  
  # Try doing this with a geometric mean, which uses the log-sum-exp trick
  # gm_mean = function(x, na.rm=TRUE){
  #   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  # }
  
  snp <- "starting_placeholder"
  for(i in 1:n){
    sExon <- sQTLs$pid[i]
    
    if(snp != sQTLs$vid[i]){
      snp = sQTLs$vid[i]
      gt_i <- enframe(get_gt_info(snp, VCF), name = "name", value = "gt")
    }
    
    psi_i <- enframe(unlist(psi[sExon,]), name = "name", value = "psi")
    
    df <- inner_join(gt_i, psi_i) 
    df <- df %>% filter(gt == "1/1" | gt == "0/1" | gt == "0/0")
    
    # ggplot(df, aes(gt, psi)) + 
    #   geom_boxplot()
    
    x <-
      df %>% 
      split(.$gt) %>%
      map("psi") %>%
      #map_dbl(mean, na.rm = T)
      map_dbl(median, na.rm = T)
    out[i,] <- x[gt]
    
    y <- 
      df %>% 
      split(.$gt) %>%
      map("psi") %>%
      map(is.na) %>%
      map(not) %>%
      map_dbl(sum)
    
    out_n[i,] <- y[gt]
  }
  
  
  psi_flipper <- function(beta, h1, h2, h3){
    if(beta > 0){
      return(c(h1, h2, h3))
    } else {
      return(c(h3, h2, h1))
    }
  }
  
  out_beta_aware <- matrix(nrow = n, ncol = 3)
  colnames(out_beta_aware) <- colnames(out)
  
  out_n_beta_aware <- matrix(nrow = n, ncol = 3)
  colnames(out_n_beta_aware) <- colnames(out)
  
  for(i in 1:n){
    out_beta_aware[i,] <- psi_flipper(sQTLs$slope[i], out[i, 1], out[i, 2], out[i, 3])
    out_n_beta_aware[i,] <- psi_flipper(sQTLs$slope[i], out_n[i,1], out_n[i,2], out_n[i,3])
  }
  
  colnames(out_beta_aware) <- c("mean_00_psi", "mean_01_psi", "mean_11_psi")
  colnames(out_n_beta_aware) <- c("n00", "n01", "n11")
  sQTL_af <- (2*out_n_beta_aware[,3] + out_n_beta_aware[,2]) / (rowSums(out_n_beta_aware) * 2)
  
  sqtl_psi <- data.frame(out_beta_aware, out_n_beta_aware, sQTL_af, check.names = F)
  
  
  # Calculate delta psi
  sqtl_psi$delta_psi <- apply(sqtl_psi, 1, function(x) mean(c(x[2]-x[1], x[3]-x[2]), na.rm = T))
  
  # Combine info for output
  sqtl_psi <- cbind(sQTLs, sqtl_psi)
  
  individual_delta_psi_by_tiss[[prefix]] <- sqtl_psi
}

# To run this whole thing in parallel
#tisss=(Adipose_Subcutaneous_gencode_v26 Cells_Cultured_Fibroblasts_gencode_v26 Muscle_skeletal_gencode_v26 Whole_Blood_gencode_v26 Pituitary Brain_Cortex Skin_Sun_Exposed)
#for t in ${tisss[@]}; do bash run_delta_psi_calculator_batch.sh -t $t; done
