# Figure 2A: Colocalization Analysis. 
# 
# This script will deal with everything colocalization. 
# 
# Jonah Einson
# 6/8/22

rm(list = ls())
library("here")
dr_here()

source(here("myPackages.R"))
source(here("sqtl_manuscript_functions.R"))
library(ggplot2)


# Read in colocalization results for downstream processing
cluster_path = "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/"
results <- read_tsv(paste0(cluster_path, "combined_coloc_results_full.tsv"))

# Add gene information to each co-localizing exon
exon_gene_map <- map_dfr(
  list.files(here("../sQTL_v8_anno/exon_gene_maps"), full.names = T), 
  read_delim, " ") 
exon_gene_map <- distinct(exon_gene_map)
results <- left_join(results, exon_gene_map, by = c("phenotype_id" = "ExonID"))

# Add columns for pp_coloc and pp_power, and filter based on this
results$pp_power <- with(results, PP.H3.abf + PP.H4.abf)
results$pp_coloc <- with(results, PP.H4.abf / (PP.H3.abf + PP.H4.abf))

# Only include genes that are in the top sQTL set
top_sQTLs <- read_tsv(here("data/top_sQTLs_MAF05.tsv"))

results_top_coloc <- 
  results %>%
  filter(phenotype_id %in% top_sQTLs$top_pid) %>%
  group_by(GeneID) %>%
  filter(pp_power + pp_coloc == max(pp_power + pp_coloc)) %>%
  dplyr::slice(1) %>%
  #slice(1) %>%
  ungroup

# Add the Euclidean distance to (1,1)
results_top_coloc %<>%
  mutate(dist = sqrt((1-pp_power)^2+(1-pp_coloc)^2))

# Plot PP-coloc vs. PP-power
ggplot(results_top_coloc, aes(pp_coloc, pp_power, color = dist < .25)) + 
  geom_point() + 
  scale_color_manual(values = c("red", "blue"))

# Test what Euclidean distance optimizes the number of things we call colocalized
pct <- function(x) sum(x) / length(x)
x = seq(0, 1.2, by = .05)
y = sapply(x , function(x) pct(results_top_coloc$dist < x))
plot(x, y, type = 'o')

# Define colocalization as having a Euclidean distance < .25 from 1,1
results_top_coloc$has_coloc <- results_top_coloc$dist < .25
barplot(table(results_top_coloc$has_coloc), col = "cornflowerblue")
text(x = c(1,2), y = table(results_top_coloc$has_coloc)-100, labels = table(results_top_coloc$has_coloc))

# Save this file for downstream analyses. I think this will make our lives easier. 
write_tsv(results_top_coloc, here("data/top_sQTLs_with_top_coloc_event.tsv"))

