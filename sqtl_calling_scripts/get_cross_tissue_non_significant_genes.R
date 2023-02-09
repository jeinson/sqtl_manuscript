#!/nfs/sw/R/R-3.4.1/bin/R

##########################
# This script takes the sQTLs and selects genes that are not significant in 
# any tissue. It will save the exon with the most significant sQTL that is still
# far away from being significant. 
# 
# Jonah Einson jeinson@nygenome.org
# 10/19/21
##########################

rm(list = ls())
source("~/myPackages.R")
library("ggplot2")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/")

fp <- list.files("combined_qtltools_results", full.names = T)

# 7/27/20 Update: Remove Testis, since it causes weird issues
fp <- fp[!grepl("Testis", fp)]

combined_sQTLs <- data.frame()
for(p in fp){
  tiss = str_split(p, "/")[[1]][2]
  tiss = str_remove(tiss, "_combined_sQTLs.tsv")
  message(tiss, "\n")
  input <- stfu(read_tsv(p))
  input$tiss = tiss
  
  combined_sQTLs <- rbind(combined_sQTLs, input)
}

# Check that the allele frequency is still correct
#combined_sQTLs %<>% filter(sQTL_af > .1 & sQTL_af < .9)
combined_sQTLs %<>% filter(sQTL_af > .05 & sQTL_af < .95)

# Get the most significant exon/tissue per gene
most_sig_sQTL_per_gene <-
  combined_sQTLs %>%
  group_by(group) %>%
  filter(bpval == min(bpval)) 

# interesting to notice that almost all variants are nominally significant, 
# but the beta p-value gives something a bit more reasonable. 

plot(ecdf(most_sig_sQTL_per_gene$bpval), xlab = "Beta P-Value", ylab = "pct total")
abline(v = 0.05, lty = 3, col = "red")
abline(h = .5, lty = 3, col = "blue")

non_significant_genes <- filter(most_sig_sQTL_per_gene, bpval > .2)

# Save the results
write_tsv(non_significant_genes, "cross_tissue_top_sQTLs/cross_tissue_nonsignificant_genes.tsv")
