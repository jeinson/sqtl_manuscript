################################################################################
#
# Script for combining sQTL results from qtltools. It saves a dataframe in the 
# format for haplotype calling using the haplotype_config python script
# 
# Jonah Einson
# 1/3/2020
#
################################################################################

args = commandArgs(trailingOnly=TRUE)
#tissue = "Adipose_Subcutaneous" # This is the name of the directory containing info
tissue = args[1]

source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno")

# Read in sQTL results

fp <- list.files(sprintf("qtltools_results/%s_MAF10", tissue), full.names = T)
fp <- fp[order(nchar(fp), fp)]

group_header <- read_lines("files/qtltools_group_header.txt")

out <- data.frame()
for(p in fp){
  x <- read_delim(p, delim = " ", col_names = F)
  out <- rbind(out, x)
}

names(out) <- group_header

# Save the output
write_tsv(
  out, 
  #sprintf("combined_qtltools_results/%s_combined_sQTLs.tsv", tissue)
  sprintf("combined_qtltools_results/%s_MAF10_combed_sQTLs.tsv", tissue)
)
