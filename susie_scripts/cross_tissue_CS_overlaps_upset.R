# Make and UPset plot of shared causal variants in e/s-Genes across all tissues

rm(list = ls())
source("~/myPackages.R")

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/")


# # This file was made in analysis/collapse_eQTLs.R
# eqtls <- read_tsv("GTEx_ge_cross_tissue_collapsed_cs.tsv.gz")
# 
# # This file comes from credible_set_Kerimov.R
# sqtls <- read_tsv("sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_cross_tissue_collapsed_sqtl_credible_set.tsv")
# 
# ####### Analysis ######
# nrow(eqtls)
# nrow(sqtls)
# 
# # Reorganize the data by nesting for simpler comparisons
# eqtls_nest <- eqtls %>%
#   select(-molecular_trait_id, -chromosome, -position, -ref, -alt) %>%
#   group_by(gene, cs_id, finemapped_region, cross_tissue_cs_id, tissue_in_cs, cs_index) %>%
#   summarise(variant = list(variant))
# 
# sqtls_nest <- 
#   sqtls %>%
#   select(-chr, -pos, -ref, -alt) %>%
#   group_by(gene_id, finemapped_region, collapsed_cs_id, cross_tissue_cs_id) %>%
#   summarise(variant_id = list(variant_id), exons_in_cs = list(exons_in_cs))
# 
# # Add some information about the number of variants per set
# sqtls_nest$n_var_per_cs <- sapply(sqtls_nest$variant_id, length)
# eqtls_nest$n_var_per_cs <- sapply(eqtls_nest$variant, length)
# 
# ggplot() + 
#   geom_histogram(aes(sqtls_nest$n_var_per_cs), fill = "blue", alpha = .5) + 
#   geom_histogram(aes(eqtls_nest$n_var_per_cs), fill = "red", alpha = .5)
# 
# library(UpSetR)
# 
# split(sqtls, 'cross_tissue_cs_id')


# You know what, I don't really now what I'm doing... So let's redo this 
# from scratch on both datasets simultaneously, and hopefully get them 
# as similar looking as possible. 

tissues <- read_tsv("gtex_qtl_catalogue_tissue_key.txt", 
                    col_names = c("gtex_id", "qtl_catalogue_id"))

# Read in collapsed sQTl credible sets
sqtl_cs <- 
  map_dfr(tissues$gtex_id, function(x){
  fn <- sprintf("sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt.gz", tiss)
  read_tsv(fn)
}
)
