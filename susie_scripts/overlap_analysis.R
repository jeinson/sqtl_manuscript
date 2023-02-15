# Analyze the overlaps found between credible sets on a tissue by 
# tissue basis. This script will start with the output from 
# tissue_by_tissue_cs_overlap.R


rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")


# Load in some files with overlaps

completed_tissues <- read_lines("../sQTL_v8_anno/completed_tissues.txt")

fps <- list.files("credible_set_overlaps", full.names = T)
names(fps) <- fps %>% str_remove("credible_set_overlaps/") %>% str_remove("_cs_overlaps.tsv")

overlaps <- lapply(fps, read_tsv)
names(overlaps) <- names(fps)

sapply(overlaps, nrow) %>%
  sort() %>%
  enframe("Tissue", "Overlapping Credible Sets") %>%
  mutate(Tissue = factor(Tissue, levels = Tissue)) %>%
  ggplot(aes(Tissue, `Overlapping Credible Sets`)) + 
  geom_bar(stat = "identity", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  coord_flip()

mean(sapply(overlaps, nrow))
# How do I even go about analyzing this? Can I look at the effect sizes of 
# eQTLs that overlap? 
