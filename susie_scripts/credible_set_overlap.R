# This script will look for patterns of overlap between psi-sQTL and eQTL 
# credible sets. This is still a pilot analysis. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

# start with Skeletal muscle tissue
tissue = "Lung"

# Read in sQTL cred sets
in_path <- sprintf("sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt.gz", tissue)
sqtl_cs <- read_tsv(in_path)

# Read in eQTL cred sets
in_path_eqtl <- sprintf("GTEx_ge_cred_sets/GTEx_ge_%s.purity_filtered.txt.gz", 
                        tolower(tissue))
eqtl_cs <- read_tsv(in_path_eqtl)

# Process data so it can be overlapped. 
sqtl_cs$gene_id <- remove_trailing_digit(sqtl_cs$gene_id)
shared_genes <- intersect(sqtl_cs$gene_id, eqtl_cs$molecular_trait_id)

# For each gene, see how many variant intersect per credible set, per gene
overlapping_set_list <- list(); i = 1
for(gene in shared_genes){
  s_gene <- filter(sqtl_cs, gene_id == gene)
  e_gene <- filter(eqtl_cs, molecular_trait_id == gene)
  
  # Split the lists by credible set, and see how many overlaps there are
  s_gene_cs <- split(s_gene, s_gene$collapsed_cs_id)
  e_gene_cs <- split(e_gene, e_gene$cs_id)
  
  for(s_cs in s_gene_cs){
    for(e_cs in e_gene_cs){
      overlaps <- intersect(s_cs$variant_id, e_cs$variant)
      if(length(overlaps) > 0) {
        overlapping_set_list[[i]] <- 
          tibble(gene_id = e_cs$molecular_trait_id[1], 
                 eqtl_cred_set = e_cs$cs_id[1], 
                 eqtl_cred_set_size = e_cs$cs_size[1],
                 sqtl_cred_set = s_cs$collapsed_cs_id[1], 
                 sqtl_cred_set_size = s_cs$cs_size[1],
                 shared_vars = overlaps) 
        i = i+1
      } else {
        overlapping_set_list[[i]] <- 
          tibble(gene_id = e_cs$molecular_trait_id[1], 
                 eqtl_cred_set = e_cs$cs_id[1], 
                 eqtl_cred_set_size = e_cs$cs_size[1],
                 sqtl_cred_set = s_cs$collapsed_cs_id[1], 
                 sqtl_cred_set_size = s_cs$cs_size[1],
                 shared_vars = NULL)
        i = i+1
      }
    }
  }
}

# This table contains all variants that between two eQTL and sQTL credible
# sets. 
overlapping_set <- bind_rows(overlapping_set_list)
write_tsv(overlapping_set, "eqtl_sqtl_credible_set_overlaps_v1.tsv")



