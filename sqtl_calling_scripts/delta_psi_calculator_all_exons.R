################################
# This script gets delta PSI values for each sQTL, even the non-significant ones
# 
# Jonah Einson
# 3/3/2020
#
# Updated 5/27/21 to use the directory with GTEx v8 annotation files
# 
# USAGE: Rscript --vanilla delta_psi_calculator.R <tissue_prefix> <chrom>
# 
# This script will read in the output of qtltools, select sGene - snp
# pairs, flip cases so 1 represents the higher included sQTL genotype, 
# then write delta-psi values. 
################################

source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/")
source("sQTL/QTL_plotter.R")

cn <- read_lines("sQTL/files/qtltools_group_header.txt") 

# Parse the command line for a tissue to use
args <- commandArgs(trailingOnly = TRUE)
prefix = as.numeric(args[1])

completed_tissues <- read_lines("sQTL_v8_anno/completed_tissues.txt")
prefix <- completed_tissues[prefix]

print(prefix)
psi <- read.csv(sprintf("sQTL_v8_anno/bed/%s_psi_v8_protocol.bed.gz", prefix),
                sep = "\t", 
                row.names = 4, 
                check.names = F)
psi <- psi[,-(1:3)]

for(chrom in 1:22){
  print(chrom)
  
  # Read in sQTLs and psi info
  #fn <- list.files(sprintf("sQTL_updated/qtltools_results/%s_MAF10", prefix), full.names = T)
  fn <- list.files(sprintf("sQTL_v8_anno/qtltools_results/%s", prefix), full.names = T)
  fn <- fn[order(nchar(fn), fn)]
  fn <- fn[chrom]
  
  sQTLs <- suppressMessages(read_delim(fn, col_names = cn, delim = " ", progress = F))
  
  # Limit to significant sQTLs
  sQTLs <- filter(sQTLs, bpval < .05)
  sQTLs <- filter(sQTLs, !duplicated(vid))
  
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
  
  # Make a dictionary to convert GTEx Exon IDs to to their coordinates, to 
  # identify them. 
  exon_id_map <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/gtex_v8_exon_id_map.tsv")
  exon_id_map <- deframe(select(exon_id_map, GTEx_ID, ID))
  
  for(i in 1:n){
    snp <- sQTLs$vid[i]
    sExon <- sQTLs$top_pid[i]
    sExon <- exon_id_map[sExon]
    
    # Use our handy dicionary to get this in the right format
    
    gt_i <- enframe(get_gt_info(snp, VCF), name = "name", value = "gt")
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
  
  #### scratch zone ####
  #Plot the significant PSI values
  # sig <- which(sQTLs$bpval < .05)
  # plot(NULL, ylim = c(0, 1), xlim = c(1,3)); k = 1
  # for(i in sig){
  #   lines(out_beta_aware[i,], col = k, type = "p", pch = 8);
  #   lines(out_beta_aware[i,], col = k, type = "l")
  #   k <- k + 1
  # }
  # 
  # x <- which((out_beta_aware[,3] - out_beta_aware[,2]) < .01 &
  #              (out_beta_aware[,2] - out_beta_aware[,1]) < .01)
  # weird <- intersect(x, sig)
  # 
  # 
  # plot(NULL, ylim = c(0, 1), xlim = c(1,3)); k = 1
  # for(i in weird){
  #   lines(out_beta_aware[i,], col = k); k <- k + 1
  # }
  # 
  # out_n[weird,]
  # 
  # snp <- sQTLs$vid[50]
  # sExon <- sQTLs$top_pid[50]
  # 
  # gt_i <- enframe(get_gt_info(snp, VCF), name = "name", value = "gt")
  # psi_i <- enframe(unlist(psi[sExon,]), name = "name", value = "psi")
  # 
  # df <- inner_join(gt_i, psi_i) 
  # df <- df %>% filter(gt == "1/1" | gt == "0/1" | gt == "0/0")
  # 
  # library(ggplot2)
  # ggplot(df, aes(gt, psi)) +
  #   geom_boxplot()+
  #   stat_summary(
  #     fun.data = stat_box_data, 
  #     geom = "text"
  #   )
  #############################
  
  colnames(out_beta_aware) <- c("mean_00_psi", "mean_01_psi", "mean_11_psi")
  colnames(out_n_beta_aware) <- c("n00", "n01", "n11")
  sQTL_af <- (2*out_n_beta_aware[,3] + out_n_beta_aware[,2]) / (rowSums(out_n_beta_aware) * 2)
  
  sqtl_psi <- data.frame(out_beta_aware, out_n_beta_aware, sQTL_af, check.names = F)
  rownames(sqtl_psi) <- sQTLs$vid
  
  # Calculate delta psi
  sqtl_psi$delta_psi <- apply(sqtl_psi, 1, function(x) mean(c(x[2]-x[1], x[3]-x[2]), na.rm = T))
  sqtl_psi <- rownames_to_column(sqtl_psi, "ID")
  
  # Combine info for output
  sqtl_psi <- cbind(sQTLs, sqtl_psi)
  
  # Save the output
  #out_path <- sprintf("sQTL_updated/combined_sQTLs_w_delta_PSI/%s.chr%i_MAF10_all_exon_delta_psi.tsv", prefix, chrom)
  out_path <- sprintf("sQTL_v8_anno/combined_sQTLs_w_delta_PSI/%s.chr%i_MAF05_all_exon_delta_psi.tsv", prefix, chrom)
  write_tsv(sqtl_psi, out_path)
}

# To run this whole thing in parallel
#tisss=(Adipose_Subcutaneous_gencode_v26 Cells_Cultured_Fibroblasts_gencode_v26 Muscle_skeletal_gencode_v26 Whole_Blood_gencode_v26 Pituitary Brain_Cortex Skin_Sun_Exposed)
#for t in ${tisss[@]}; do bash run_delta_psi_calculator_batch.sh -t $t; done
