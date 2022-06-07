# This script just re-plots the analyses from Mariia, for consistency with 
# the rest of the manuscript. 

rm(list = ls())

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sqtl_manuscript")
source("sqtl_manuscript_functions.R")
source("~/myPackages.R")

# Read in data
fps <- list.files("stuff_from_mariia", full.names = T)
# [1] "stuff_from_mariia/to_draw_inclusion_groups_overall.csv"       
# [2] "stuff_from_mariia/to_draw_inclusion_groups_random_shuffle.csv"
# [3] "stuff_from_mariia/to_draw_inclusion_groups.csv"               
# [4] "stuff_from_mariia/to_draw_psi_groups.csv"                     
# [5] "stuff_from_mariia/to_draw_sQTL.csv"  

inclusion_group_overall <- read_csv(fps[1])
inclusion_groups_random <- read_csv(fps[2])
inclusion_groups <- read_csv(fps[3])
psi_groups <- read_csv(fps[4])
sQTL_res <- read_csv(fps[5])

results <- map_dfr(fps[-c(2,4)], read_csv)

# Figure 3
library(ggplot2)
ggplot(results %>% filter(test_type == "ks"), 
       aes(statistics, trait, xmin = conf_low, xmax = conf_top)) + 
  geom_pointrange(aes(color = pair), position = position_dodge(.5)) + 
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Average difference between normalized means") + 
  ylab("Exon Domain Feature") + 
  sqtl_manuscript_theme()

# Supplements to Figure 3
# What do distributions of individual features look like?
sqtl_details <- read_csv("../splicing_project/Data/combined_sQTL_data.csv")
sqtl_details