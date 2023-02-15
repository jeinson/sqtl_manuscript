# Make some summary plots for credible collapsed credible sets from PSI sQTLs

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

tissues = read_lines("../sQTL_v8_anno/completed_tissues.txt")
tissues = tissues[-which(tissues == "Colon_Transverse")]
collapsed_cs <- list()

for(tiss in tissues){
  ########### Load and process Credible Sets ##########
  # Read in CS data
  collapsed_cs[[tiss]] <- read_tsv(
    sprintf(
      "sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt.gz", 
      tiss))
}

# How many credible sets do we get per tissue group?
n_CS_per_tiss <- sapply(collapsed_cs, function(x) length(unique(x$collapsed_cs_id)))
n_sGenes_per_Tiss <- sapply(collapsed_cs, function(x) length(unique(x$gene_id)))
                            
n_cs_per_df <- 
  cbind(enframe(n_CS_per_tiss, name = "Tissue", value = "n_CS_per_tiss"), n_sGenes_per_Tiss) %>%
  arrange(desc(n_CS_per_tiss)) %>%
  mutate(Tissue = fct_inorder(Tissue)) 
mean(n_cs_per_df$n_CS_per_tiss / n_cs_per_df$n_sGenes_per_Tiss)

n_cs_per_df %>%
  pivot_longer(cols = starts_with("n"), names_to = "count") %>%
  ggplot(aes(Tissue, value, fill = count)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ggtitle("Number of sGenes and sQTL credible sets across tissues")
