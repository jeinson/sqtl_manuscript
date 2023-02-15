# In this script, I apply the same collapsing procedure on GTEx eQTLs, so I can 
# compare them more easily. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

fps <- list.files("GTEx_ge_cred_sets", full.names = T, pattern = "GTEx")

tissues <- fps %>% 
  str_split("/") %>% 
  map_chr(2) %>% 
  str_remove("GTEx_ge_") %>%
  str_remove(".purity_filtered.txt.gz")
names(fps) <- tissues

# Read in everything at once
all_eqtl_cs <- lapply(fps, read_tsv)
all_eqtl_cs <- bind_rows(all_eqtl_cs, .id = 'tissue')

all_eqtl_cs$cs_id_tiss <- 
  paste(all_eqtl_cs$cs_id, all_eqtl_cs$tissue, sep = "_")

# Run the same filter they wrote about in Kerimov et. al, but apparently didn't
# apply themselves here... (Max 30 vars in size, max Z score > 3)
all_eqtl_cs <- 
  all_eqtl_cs %>%
  group_by(cs_id_tiss) %>%
  filter(n() <= 30 & max(abs(z)) > 3) %>%
  ungroup

### Collapse across tissues ----
# Run roughly the same process as in the cross tissue sQTL set collapsing 
all_cs_bygene <- split(all_eqtl_cs, all_eqtl_cs$molecular_trait_id)

cross_tissue_gene_cs <- lapply(all_cs_bygene, function(x){
  cred_sets <- split(x$variant, x$cs_id_tiss)
  
  # Get overlaps
  cross_tissue_credible_set <- collapse(cred_sets)
  names(cross_tissue_credible_set) <- paste0("CT_CS", seq_along(cross_tissue_credible_set))
  
  # Now assign the variants to the final credible sets. 
  cs_key <- stack(cross_tissue_credible_set)
  
  cs_key$ind <- paste0(x$molecular_trait_id[1], "_", cs_key$ind)
  cs_key <- deframe(cs_key)
  x$cross_tissue_cs_id <- cs_key[x$variant]
  
  # Add a label for tissues in the cs
  x$tissue_in_cs <- paste(unique(x$tissue), collapse = ",")
  
  # Remove duplicates variants that originated from different credible sets
  x <- x %>% select(-tissue, -(pip:cs_id_tiss)) %>% distinct
  
  # Return final list of variants per credible set
  x
})

# Collapse everythig down into one big data frame
all_eqtl_cs_collapsed <- bind_rows(cross_tissue_gene_cs, .id = "gene")

# Save to disk
write_tsv(all_eqtl_cs_collapsed, "GTEx_ge_cross_tissue_collapsed_cs.tsv")
system("gzip -f GTEx_ge_cross_tissue_collapsed_cs.tsv")
