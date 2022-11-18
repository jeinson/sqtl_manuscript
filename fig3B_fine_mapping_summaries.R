# Scripts for Figure 2B, which looks at fine-mapping of sQTLs. 
# 

library(here)
source(here("myPackages.R"))

cs_path <- "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/sqtl_v8_anno_collapsed_credible_sets/"
tissues <- read_lines(here('data/completed_tissues.txt'))
tissues <- tissues[-which(tissues == "Colon_Transverse")]
collapsed_cs_list <- list()
for(tiss in tissues){
  message(tiss)
  cs <- read_tsv(paste0(cs_path, "GTEx_psi_", tiss, ".collapsed.txt.gz"))
  cs$exons_in_cs <- str_split(cs$exons_in_cs, ",")
  
  cs %<>% 
    select(-c(tissue, pos, ref, alt, chr)) %>% 
    group_by(gene_id, finemapped_region, collapsed_cs_id, exons_in_cs) %>% 
    summarise(variant_id = list(variant_id)) 
  
  collapsed_cs_list[[tiss]] <- cs
}

# Summarize across datasets
n_variants_in_cs <- 
  sapply(collapsed_cs_list, function(x){
    map_dbl(x$exons_in_cs, length)
  })

mean(unlist(n_variants_in_cs))

cs_per_gene_per_tiss <- 
  sapply(collapsed_cs_list, function(x){
    table(x$gene_id) 
  }
  )

mean(unlist(cs_per_gene_per_tiss))
