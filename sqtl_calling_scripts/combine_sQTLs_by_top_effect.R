#!/nfs/sw/R/R-3.4.1/bin/R

##########################
# This script takes the sQTLs with delta PSI calculated, and selects the 'top sQTL' across all tissues where the 
# sQTL is significant. It draws some plots to help understand the shape of the data as well. 
# 
# Jonah Einson (This header was written 9/8/20. Script itself was written way earlier(
##########################

rm(list = ls())
source("~/myPackages.R")
library("ggplot2")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/")

fp <- list.files("combined_sQTLs_w_delta_PSI", full.names = T, pattern = "MAF05")

# 7/27/20 Update: Remove Testis, since it causes weird issues
fp <- fp[!grepl("Testis", fp)]

combined_sQTLs <- data.frame()
for(p in fp){
  tiss = str_split(p, "/")[[1]][2]
  tiss = gsub(".chr.+", "", tiss)
  
  input <- stfu(read_tsv(p))
  input$tiss = tiss
  
  combined_sQTLs <- rbind(combined_sQTLs, input)
}

# Check that the allele frequency is still correct
#combined_sQTLs %<>% filter(sQTL_af > .1 & sQTL_af < .9)
combined_sQTLs %<>% filter(sQTL_af > .05 & sQTL_af < .95)

# Make an upset plot showing the overlap of different sQTLs
sig_gene_list <- 
  combined_sQTLs %>% 
  split(.$tiss) %>%
  map("group")

library(UpSetR)
svg("figures/tissue_sQTL_upset.svg", width = 10, height = 6)
upset(fromList(sig_gene_list), nsets = length(sig_gene_list), order.by = "freq", 
      mainbar.y.label = "Significant sQTLs")
dev.off()

# Work some dplyr magic to get this looking good
combined_sQTLs %>%
  # group_by(group) %>%
  # mutate(n_tiss = n()) %>%
  # mutate(n_unique_top_exons = length(unique(top_pid))) %>%
  # arrange(desc(delta_psi)) %>%
  # slice(1) %>%
  split(.$group) %>%
  map(~ mutate(.x, n_tiss = nrow(.x))) %>%
  map(~ filter(.x, .x$delta_psi == max(.x$delta_psi))) %>%
  bind_rows() %>%
  filter(delta_psi > 0) ->
  top_sQTLs

# Randomly choose the winner in the few cases where there's a tie
top_sQTLs <- 
  top_sQTLs %>%
  group_by(group) %>%
  filter(row_number() == 1) %>%
  ungroup

# How many genes now have a significant sQTL?
nrow(top_sQTLs)

top_sQTLs

# What percentage of tissues accounts for the most sQTLS?
# (Make some summary stats plots)
ggplot(combined_sQTLs, aes(fct_infreq(tiss))) + 
  geom_bar(fill = "cornflowerblue") + 
  ggtitle("Number of significant sQTLs per tissue") + 
  xlab("GTEx Tissue") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(top_sQTLs, aes(tiss)) + 
  geom_bar(fill = "darkgoldenrod") +
  ggtitle("Number of genes for which each tissue has the top sQTL")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(top_sQTLs, aes(n_tiss)) +
  geom_bar(fill = "darkgreen") +
  ggtitle("Number of tissues covering each gene") 

ggplot(top_sQTLs, aes(delta_psi)) + 
  geom_histogram(fill = "chocolate1", col = "black") + 
  ggtitle("Distribution of ΔPSI scores across top sQTLs") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(top_sQTLs, aes(factor(n_tiss), delta_psi)) + 
  geom_boxplot() + 
  ggtitle("ΔPSI by the number of tissues in which the gene is covered")

# Visualize a few significant sQTLs
source("../pipeline/scripts/QTL_plotting_function.R")


# This sQTL's top effect is in muscle
# muscle_psi <- read_tsv("bed/Muscle_skeletal_gencode_v26_psi_v8_protocol.bed.gz")
# adipose_psi <- read_tsv("bed/Adipose_Subcutaneous_gencode_v26_psi_v8_protocol.bed.gz")
# cortex_psi <- read_tsv("bed/Brain_Cortex_psi_v8_protocol.bed.gz")
# pituitary_psi <- read_tsv("bed/Pituitary_psi_v8_protocol.bed.gz")
# blood_psi <- read_tsv("bed/Whole_Blood_gencode_v26_psi_v8_protocol.bed.gz")
# testis_psi <- read_tsv("bed/Testis_psi_v8_protocol.bed.gz")
# vcf_path = "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
# 
# top_sQTLs %>% filter(n_tiss == 7) %>% .[1,] %>% t
# sid = "chr6_29956612_A_G_b38"
# gene = "chr6_30061917_30062017_+"
# 
# this_sqtl <- top_sQTLs %>% filter(n_tiss == 7) %>% .[1,] %>% unlist
# 
# gene = this_sqtl['top_pid']
# sid = this_sqtl['ID']
# 
# 
# muscle_plt = plot_QTL(sid = sid, 
#                       gene = gene, 
#                       vcf = vcf_path,
#                       phenotypes = muscle_psi
# ) + ggtitle("Muscle")
# 
# adipose_plt = plot_QTL(sid = sid, 
#                        gene = gene,  
#                        vcf = vcf_path,
#                        phenotypes = adipose_psi
# ) + ggtitle("Adipose")
# 
# cortex_plt = plot_QTL(sid = sid, 
#                       gene = gene, 
#                       vcf = vcf_path,
#                       phenotypes = cortex_psi
# ) + ggtitle("Cortex")
# 
# pituitary_plt = plot_QTL(sid = sid, 
#                          gene = gene, 
#                          vcf = vcf_path,
#                          phenotypes = pituitary_psi
# ) + ggtitle("Pituitary")
# 
# library(patchwork)
# (muscle_plt | adipose_plt) / (cortex_plt | pituitary_plt)

# Save the results, and make sure to include the date!
top_sQTLs <- top_sQTLs[with(top_sQTLs, order(nchar(chr), chr)),]
write_tsv(top_sQTLs, "cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")
