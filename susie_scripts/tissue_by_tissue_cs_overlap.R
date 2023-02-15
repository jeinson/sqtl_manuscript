# Explore tissue-by-tissue overlap between psi-sQTL and eQTL credible sets. 
# We would like to look for cases where splicing in a non-symmetric exon 
# triggers an NMD event. Most interestingly, the strongest eQTLs could be 
# triggered by this type of mechasnism. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

# Look at all tissues
tissues = read_lines("../sQTL_v8_anno/completed_tissues.txt")
tissue_key <- read_tsv("gtex_qtl_catalogue_tissue_key.txt", col_names = F) %>% 
  deframe()
tissues = tissues[-which(tissues == "Colon_Transverse")]

for(tiss in tissues){
  tiss_lc = tissue_key[tiss]
  
  ########### Load and process Credible Sets ##########
  # Read in CS data
  sqtl_cs <- read_tsv(
    sprintf(
      "sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt.gz", 
      tiss))
  
  eqtl_cs <- read_tsv(
    sprintf(
      "GTEx_ge_cred_sets/GTEx_ge_%s.purity_filtered.txt.gz",
      tiss_lc
    )
  )
  
  # Fix column names in sqtl_cs to match what's in eqtl_cs
  sqtl_cs <- rename(sqtl_cs, 
                    "molecular_trait_id" = "gene_id", 
                    "variant" = "variant_id", 
                    "chromosome" = "chr", 
                    "position"= "pos") %>%
    mutate(finemapped_region = str_remove(finemapped_region, "chr")) %>%
    mutate(molecular_trait_id = remove_trailing_digit(molecular_trait_id))
  
  # Apply filtering used in Kerimov et. al
  eqtl_cs <- 
    eqtl_cs %>%
    group_by(cs_id) %>%
    filter(n() <= 30 & max(abs(z)) > 3) %>%
    ungroup
  
  ######### Analysis ##########
  
  # 1) How much does this credible set overlap in general?
  # Side note. This looks like shit. How does this overlap compare to other 
  # splicing phenotypes from the qtl catalogue?
  # txqtl_cs <- read_tsv("GTEx_txrevise_cred_sets/GTEx_txrev_nerve_tibial.purity_filtered.txt.gz")
  # txqtl_cs$gene <- txqtl_cs$molecular_trait_id %>% str_split(pattern = "\\.") %>% map_chr(1)
  # 
  genes_per_assay_list <-
    list(sGenes = unique(sqtl_cs$molecular_trait_id),
         eGenes = unique(eqtl_cs$molecular_trait_id)
    )
  
  # Get overlaps over genes with an eQTL and an sQTL
  esGenes <- intersect(eqtl_cs$molecular_trait_id, sqtl_cs$molecular_trait_id)
  
  library(UpSetR)
  upset_fp <- paste0("figures/", tiss, "_overlapping_esGenes.png")
  png(filename = upset_fp)
  upset(fromList(genes_per_assay_list))
  dev.off()
  
  # 2) For genes that have both an eQTL and an sQTL, how many actually share a 
  # causal variant?
  
  # For every esGene, check if there is an overlap of variants in the credible
  # set. Label the genes where there is an overlap
  overlaps <- 
    tibble(eCS_id = character(), 
           sCS_id = character(), 
           variants = character()
    )
  
  for(gene in esGenes){
    eCS <- filter(eqtl_cs, molecular_trait_id == gene) %>% split(.$cs_id)
    sCS <- filter(sqtl_cs, molecular_trait_id == gene) %>% split(.$collapsed_cs_id)
    
    # Do a double for loop to check for overlaps 
    for(this_eCS in eCS){
      for(this_sCS in sCS){
        bicausal_vars <- intersect(this_eCS$variant, this_sCS$variant)
        if(length(bicausal_vars) > 0){
          overlaps <- add_row(overlaps, 
                              eCS_id = this_eCS$cs_id[1],
                              sCS_id = this_sCS$collapsed_cs_id[1],
                              variants = paste(bicausal_vars, collapse = ","))
        }
      }
    }
  }
  
  # Add the exons involved
  overlaps <- 
    left_join(x = overlaps, 
              y = distinct(select(sqtl_cs, collapsed_cs_id, exons_in_cs)), 
              by = c("sCS_id" = "collapsed_cs_id"))
  
  # Plot some summary results
  n_eGenes <- length(unique(eqtl_cs$molecular_trait_id))
  n_eGene_CSs <- length(unique(eqtl_cs$cs_id))
  
  n_sGenes <- length(unique(sqtl_cs$molecular_trait_id))
  n_sGene_CSs <- length(unique(sqtl_cs$collapsed_cs_id))
  
  n_overlapping_CSs <- nrow(overlaps)
  
  # Save the result
  out_name <- paste0("figures/", tiss, "_credible_set_overlaps.png")
  
  png(out_name, width = 850, height = 600, res = 100)
  barplot(c(`eGenes` = n_eGenes, `eQTL\nCredible sets` = n_eGene_CSs, 
            `sGenes` = n_sGenes, `sQTL Multi-Exon\nCredible Sets` = n_sGene_CSs, 
            `esGenes` = length(esGenes), 
            `Overlapping\nCredible Sets` = n_overlapping_CSs), 
          col = c("dodgerblue4", "dodgerblue1", "springgreen4", 
                  "springgreen1", "tomato4", "tomato1"), 
          space = c(.2, 0, .2, 0, .2, 0),
          main = paste("Number of credible sets in", tiss))
  grid()
  dev.off()
  
  # Save the output file
  overlaps_filename <- paste0("credible_set_overlaps/", tiss, 
                              "_cs_overlaps.tsv")
  write_tsv(overlaps, overlaps_filename)
}

