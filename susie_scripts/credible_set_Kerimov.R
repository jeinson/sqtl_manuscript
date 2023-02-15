# This script will apply the same method used to select credible sets
# from sQTLs with multiple correlated phenotypes (i.e. multiple exons per gene)
# using the same selection method as in Kerimov et. al. (the QTL catalogue paper)

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

# Load all credible sets combined into one data frame
all_cs <- readRDS("all_credible_sets_12_2.rds")
all_cs <- bind_rows(all_cs, .id = "tissue")
all_cs$cs_id_tiss <- paste(all_cs$cs_id, all_cs$tissue, sep = "_")

### Quick summary stats

# The number of credible set variants per tissue (not a very useful statistic)
table(all_cs$tissue)

# The number of credible sets per tissue
n_cs_per_tiss <- 
  tapply(all_cs$cs_id, all_cs$tissue, function(x) length(unique(x))) %>%
  enframe(name = "Tissue", value = "n_cs") %>%
  arrange(n_cs) %>%
  mutate(Tissue = factor(Tissue, levels = Tissue))

library(ggplot2)
ggplot(n_cs_per_tiss, aes(y = Tissue, x = n_cs)) + 
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  xlab("Number of credible sets per tissue")

# 1) Select CSs with < 30 variants and max z-score > 3.
all_cs_f1 <- 
  all_cs %>%
  group_by(cs_id_tiss) %>%
  filter(n() < 30 & max(abs(z)) > 3) %>%
  ungroup

# 2) For each gene, select the smallest credible set in each dataset. (I think
# this means in each tissue). If there are ties, use the one with the highest 
# maximal posterior inclusion probability
# all_cs_f1 %>%
#   group_by(gene_id, cs_id_tiss) %>%
#   mutate(n_in_set = n()) %>%
#   ungroup(cs_id_tiss) %>%
#   filter(n_in_set == min(n_in_set)) %>%
#   slice(which.max(posterior_mean))

##### Collapse across genes, within tissues ----
# Alternatively, just get the highest level gene sets across
# all genes, in genes expressed per tissue. 
all_cs_bytiss <- all_cs_f1 %>% split(.$tissue)
tissues <- names(all_cs_bytiss)

# Loop across genes per tissue

# Save the output here, cause this is an lapply statement
gene_level_cs <- 
  lapply(all_cs_bytiss, function(x){
    message(x$tissue[1])
    y <- x %>% split(.$gene_id)
    for(i in seq_along(y)){
      
      # Find overlaps between all credible sets to get a finalized list 
      # "subgraphs" in this network of credible sets 
      cred_sets <- split(y[[i]]$variant_id, y[[i]]$cs_id)
      collapsed_cred_sets <- collapse(cred_sets)
      names(collapsed_cred_sets) <- paste0("FCS", seq_along(collapsed_cred_sets))
      
      # Now assign the variants to the final credible sets. 
      cs_key <- stack(collapsed_cred_sets)
      cs_key$ind <- paste0(y[[i]]$gene_id[1], "_", cs_key$ind)
      cs_key <- deframe(cs_key)
      y[[i]]$collapsed_cs_id <- cs_key[y[[i]]$variant_id]
    }
    output <- bind_rows(y)
    
    output %>%
      group_by(collapsed_cs_id) %>%
      mutate(exons_in_cs = paste(unique(phenotype_id), collapse = ",")) %>%
      ungroup %>%
      
      # Remove the exon IDs, since we're at gene level right now
      #select(-phenotype_id, -pip, -z, -cs_min_r2, -cs_avg_r2, -posterior_mean, 
      #       -posterior_sd) %>%
      select(-phenotype_id, -cs_id_tiss, -cs_id, -cs_index, -(pip:cs_log10bf)) %>%
      distinct
  })

# Save everything to their own files
for(tiss in tissues){
  out_path <- 
    sprintf(
      "sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt", 
      tiss
    )
  
  write_tsv(gene_level_cs[[tiss]], out_path)
  system(paste("gzip -f", out_path))
}

### Collapse across tissues ----
# Finally, collapse all credible sets across all tissues, to get one huge
# set of credible sets genome wide. 
cs_cross_tiss <- bind_rows(gene_level_cs)
all_cs_bygene <- split(cs_cross_tiss, cs_cross_tiss$gene_id)

cross_tissue_gene_cs <- lapply(all_cs_bygene, function(x){
  x$collapsed_cs_id_tiss <- with(x, paste0(collapsed_cs_id, "_", tissue))
  cred_sets <- split(x$variant_id, x$collapsed_cs_id_tiss)
  
  # Get overlaps
  cross_tissue_credible_set <- collapse(cred_sets)
  names(cross_tissue_credible_set) <- paste0("CT_CS", seq_along(cross_tissue_credible_set))
  
  # Now assign the variants to the final credible sets. 
  cs_key <- stack(cross_tissue_credible_set)
  
  cs_key$ind <- paste0(x$gene_id[1], "_", cs_key$ind)
  cs_key <- deframe(cs_key)
  x$cross_tissue_cs_id <- cs_key[x$variant_id]
  
  x
})

# combine rows and remove tissue information, to get a succinct list of 
# variants in each credible set. 
cross_tissue_gene_cs_full <- 
  bind_rows(cross_tissue_gene_cs) %>%
  select(-collapsed_cs_id, -collapsed_cs_id_tiss, -chr, -pos, -ref, -alt) %>%
  distinct %>%
  group_by(gene_id, finemapped_region, cross_tissue_cs_id, exons_in_cs) %>%
  summarise(tissue = paste0(unique(tissue), collapse = ","), 
            variant_id = paste0(variant_id, collapse = ","))

write_tsv(cross_tissue_gene_cs_full,
          "sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_cross_tissue_collapsed_sqtl_credible_set.tsv")
