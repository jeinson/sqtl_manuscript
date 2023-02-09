# This script will take the top sQTLs that include secondary exons, and 
# calculate the delta PSI score, including all secondary exons. 

# It is similar to delta_psi_calculator_secondary_sQTLs, but does everything
# at once, as opposed to be run one tissue at a time. 


source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/")

top_sQTLs <- read_tsv("sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs_exons_labeled_by_coord.tsv")

# Read in PSI scores from all tissues. (This gonna take a lot of RAM!)
tissues <- read_lines("sQTL_v8_anno/completed_tissues.txt")
cross_tiss_psi <- list()
for(tiss in tissues){
  message(tiss)
  psi_fp <- sprintf("sQTL_v8_anno/bed/%s_psi_v8_protocol.bed.gz", tiss)
  cross_tiss_psi[[tiss]] <- read.csv(file = psi_fp, 
                                     sep = "\t",
                                     row.names = 4, 
                                     check.names = F)
}

# VCF with GTEx genotype info
VCF = "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz" 
psi_flipper <- function(beta, h1, h2, h3){
  if(beta > 0){
    return(c(h1, h2, h3))
  } else {
    return(c(h3, h2, h1))
  }
}

# Split up the top sQTLs by gene, and process each gene group at a time
top_sQTLs_bygene <- split(top_sQTLs, top_sQTLs$group)

top_sQTLs_bygene_w_dpsi <- list()
for(gene_sqtl in top_sQTLs_bygene){
  # Grab PSI info
  tiss <- gene_sqtl$tiss[1]
  exons <- gene_sqtl$pid
  psi <- cross_tiss_psi[[tiss]][exons,-(1:3)]
  
  # Grab genotype info
  snp = gene_sqtl$vid[1]
  gt_i <- enframe(get_gt_info(snp, VCF), name = "name", value = "gt")
  
  # Loop through all exons
  gt <- c("0/0", "0/1", "1/1")
  n = nrow(gene_sqtl)
  out <- matrix(nrow = n, ncol = 3)
  out_n <- matrix(nrow = n, ncol = 3)
  for(i in 1:n){
    psi_i <- enframe(unlist(psi[i,]), name = "name", value = "psi")
    df <- inner_join(gt_i, psi_i) 
    df <- df %>% filter(gt == "1/1" | gt == "0/1" | gt == "0/0")
    
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
    
    out_beta_aware <- matrix(nrow = n, ncol = 3)
    colnames(out_beta_aware) <- colnames(out)
    
    out_n_beta_aware <- matrix(nrow = n, ncol = 3)
    colnames(out_n_beta_aware) <- colnames(out)
    
    for(i in 1:n){
      out_beta_aware[i,] <- psi_flipper(gene_sqtl$slope[i], out[i, 1], out[i, 2], out[i, 3])
      out_n_beta_aware[i,] <- psi_flipper(gene_sqtl$slope[i], out_n[i,1], out_n[i,2], out_n[i,3])
    }
    
    colnames(out_beta_aware) <- c("mean_00_psi", "mean_01_psi", "mean_11_psi")
    colnames(out_n_beta_aware) <- c("n00", "n01", "n11")
    
    sQTL_af <- (2*out_n_beta_aware[,3] + out_n_beta_aware[,2]) / (rowSums(out_n_beta_aware) * 2)
    
    # Compile everytyhing
    sqtl_psi <- data.frame(out_beta_aware, out_n_beta_aware, sQTL_af, check.names = F)
    
    # Calculate delta psi
    sqtl_psi$delta_psi <- apply(sqtl_psi, 1, function(x) mean(c(x[2]-x[1], x[3]-x[2]), na.rm = T))
    
    # Put it back together
    top_sQTLs_bygene_w_dpsi[[gene_sqtl$group[1]]] <- cbind(gene_sqtl, sqtl_psi)
  }
}

top_sQTLs_w_dpsi <- bind_rows(top_sQTLs_bygene_w_dpsi)
write_tsv(top_sQTLs_w_dpsi, "sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_w_dpsi.tsv")
