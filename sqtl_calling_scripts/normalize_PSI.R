#!/nfs/sw/R/R-3.4.1/bin/R
##########################################################################################
#
# Preprocessing and normalization of PSI data for sQTL calling
# Based on the GTEx pipeline for eQTL analysis
# It will save the "violently" normalized PSI file for downstream analysis. 
# Jonah Einson: jeinson@nygenome.org
#
##########################################################################################

# This is now on a new branch used to annotate exons using the OFFICIAL GTEx
# v8 release. 

args = commandArgs(trailingOnly=TRUE)
#tissue = "Adipose_Subcutaneous" # This is the name of the directory containing info
tissue = args[1]

# In this script, I take the output from IPSA-nf PSI calling and apply the standard normal transformation
# and run PEER to prepare for sQTL analysis

source("~/myPackages.R")
#library(tidyverse)
#library(magrittr)
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/")

# Get the paths for psi scores
psi_fp = sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/v8_psi_updated/tissues_v8_anno/%s_v8/data/all.A.psi.filtered.sorted.tsv", tissue)
psi_detail_fp = sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/v8_psi_updated/tissues_v8_anno/%s/data/psi_detail.tsv", tissue)

psi <- read_tsv(psi_fp, col_names = F, skip  = 1)
psi_colnames <- strsplit(system(sprintf("head -1 %s", psi_fp), intern = T), "\t")[[1]]
psi_colnames <- c("Exon", psi_colnames)
# Remove the tissue identifiers...
psi_colnames <- psi_colnames %>% str_split("-") %>% map(~.x[1:2]) %>% map_chr(~paste(.x, collapse="-"))
colnames(psi) <- psi_colnames
psi <- column_to_rownames(psi, "Exon-NA")

# Make sure to only keep annotated exons. This is relevant when PSI was called
# using a different annotation file which contains more exons than the gold 
# standard GTEx gencode v26 annotation. 
# exons_tabix <- 
#   rownames(psi) %>%
#   str_split("_") %>%
#   map_chr(~ paste0(.x[1], ":", .x[2], "-", .x[3]))
# keep <- (exons_tabix %in% psi_detail$X1)
# psi <- psi[keep,]
# 
# # Make sure they line up!!!
# exons_tabix <- 
#   rownames(psi) %>%
#   str_split("_") %>%
#   map_chr(~ paste0(.x[1], ":", .x[2], "-", .x[3]))
# psi_detail <- psi_detail[match(exons_tabix, psi_detail$X1),]

## Only use significantly included genes for sQTL analysis
# This is the same protocols as in the v8 data prep, as to avoid numerical 
# issues with the beta distribution. 
sample_frac_threshold = .5 # Must be present in over 50% of samples

# Exon must be present in 252 or more individuals
min_exons_required = sample_frac_threshold*ncol(psi) 

# Each exon must have more than 10 unique values or have more than
# 10% of their values be unique, whichever is greater. 
has_enough_unique_values <- function(x){
  cutoff <- max(10, 0.1 * sum(!is.na(x)))
  length(unique(x)) > cutoff
}

# Apply these two filters
exons_to_keep <- apply(psi, 1, function(x){
  enough_avail <- sum(!is.na(x)) > min_exons_required
  enough_unique_values <- has_enough_unique_values(x)
  
  enough_avail & enough_unique_values
}
)

sprintf("Number of exons total: %i", length(exons_to_keep)) # 54,366
sprintf("Number of good exons: %i", sum(exons_to_keep)) #78,230

# Only include rows which pass filtering, 
psi_filtered <- psi[exons_to_keep,]

# get rid of psi and psi_detail to save some ram
rm(psi)

# Make sure to only include individuals who have genotype data.
# it will break fastQTL if individuals not in the VCF are included. 
vcf_sample_names <- 
  read_lines("files/vcf_sample_names.txt")
tissue_indvs <- 
  names(psi_filtered) %>% 
  str_split("-") %>%
  map(~.x[1:2]) %>%
  map_chr(~ paste(.x, collapse = "-"))

psi_filtered <- 
  psi_filtered[,tissue_indvs %in% vcf_sample_names]

# Save the individuals used for downstream analysis:
write_lines(colnames(psi_filtered), sprintf("indvs_per_tiss/%s_indvs.txt", tissue))

########################## Process PSI detail more #############################
# message("Process PSI Detail")
# vec_to_named_vec <- function(x) {
#   x <- str_split(x, " ")
#   nm <- unlist(map(x, 1))
#   val <- unlist(map(x, 2))
#   names(val) <- nm
#   ln <- length(val)
#   val[ln] <- str_remove(val[ln], ";")
#   val
# }
# 
# psi_detail_filtered$X2 <- 
#   psi_detail_filtered$X2 %>%
#   str_remove_all("\"") %>%
#   str_split("; ") %>%
#   map(vec_to_named_vec) 

###################### Get Exon Annotations from GTEx ##########################
exon_id_map <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/gtex_v8_exon_id_map.tsv",
                        col_names = c("chr", "start", "end", "exon_id", "gtex_exon_id"), 
                        skip = 1)
exon_id_map$gene_id <- remove_trailing_digit(exon_id_map$gtex_exon_id)

psi_detail_filtered <- 
  exon_id_map[match(rownames(psi_filtered), exon_id_map$exon_id),]

################# Make and save a phenotype group column #######################
# We can do grouped permutations with fastQTL, where exons are assigned to a 
# gene group. 
GeneID <- psi_detail_filtered$gene_id
ExonID <- psi_detail_filtered$gtex_exon_id

# as always, remove that pesky trailing digit!
GeneID <- str_extract(GeneID, "[^\\\\.]+")
Ensembl_ExonID <- str_split(ExonID, ",") %>% 
  map(~ str_extract(.x, "[^\\\\.]+")) %>%
  map(paste, collapse = ",") %>%
  unlist

psi_filtered_w_GeneID <- 
  psi_filtered %>%
  rownames_to_column("Exon_coord") %>%
  cbind(GeneID, .) %>%
  cbind(ExonID, .) 

groups <- psi_filtered_w_GeneID[,c("ExonID", "GeneID", "Exon_coord")]
write.table(groups, 
            sep = " ", row.names = F, quote = F, 
            col.names = c("ExonID", "GeneID", "Exon_coord"), 
            file = sprintf("exon_gene_maps/%s_exon_gene_map.txt", tissue))

psi_filtered <- psi_filtered_w_GeneID %>%
  select(-"GeneID", -"ExonID")  %>%
  column_to_rownames("Exon_coord")

###################### Standard normal transforamtion ##########################
# Add a bit of noise, so there is an actual ranking and the data will be 
# normally distributed. This is of course an extremely violent process, but 
# shouldn't matter too much in the long run
message("Standard Normal Transform")
n_indvs = ncol(psi_filtered) 
stochastic_norm = apply(psi_filtered, 1, function(x){
  x = x + runif(n_indvs, min = 0, max = 0.000000001)
  qnorm((rank(x, na.last='keep') - 0.5)/sum(!is.na(x)))
})

########### Save the Results in BED format for use in fastQTL ##################
# Raw PSI values
psi_filtered_out <- 
  str_split(rownames(psi_filtered), "_") %>%
  do.call("rbind", .) %>% .[,1:3] %>%
  set_colnames(c("#Chr", "start", "end")) %>%
  cbind(rownames_to_column(psi_filtered, var = "ID"))

out_path = sprintf("bed/%s_psi_v8_protocol.bed", tissue)
write_tsv(psi_filtered_out, out_path, na = "")
system(paste("/nfs/sw/htslib/htslib-1.1/bgzip -f", out_path))
system(paste0("/nfs/sw/htslib/htslib-1.1/tabix -p bed ", out_path, ".gz"))

# clean up
rm(psi_filtered_out)

# Violently normalized PSI values (used for calculating PEER factors)
stochastic_norm <- apply(stochastic_norm, 2, round, 5)
violently_normalized_psi_out <- 
  str_split(rownames(psi_filtered), "_") %>%
  do.call("rbind", .) %>% .[,1:3] %>%
  set_colnames(c("#Chr", "start", "end")) %>%
  cbind(rownames_to_column(as.data.frame(t(stochastic_norm)), var = "ID"))

out_path = sprintf("bed/violently_normalized_%s_psi_v8_protocol.bed", tissue)
write_tsv(violently_normalized_psi_out, out_path, na = "")
system(paste("/nfs/sw/htslib/htslib-1.1/bgzip -f", out_path))
system(paste0("/nfs/sw/htslib/htslib-1.1/tabix -p bed ", out_path, ".gz"))

# clean up
rm(violently_normalized_psi_out, psi_filtered, psi_filtered_w_GeneID)

######### Save the results in a format for QTLtools, group mode ################
# Dump exons which appear in a secondary gene, or are from non-protein-coding
# or lncRNA genes :-/
TSSs <-
  read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/hg19_pc_gene_tsss.tsv")
is_pc_gene <- groups$GeneID %in% TSSs$GeneID

#exons_to_keep <- !grepl("[0-9]$", groups$ExonID)
exons_to_keep <- is_pc_gene
  #exons_to_keep & is_pc_gene

# Take care of that
qtltools_input_info_columns <- groups[exons_to_keep,]
qtltools_input <- t(stochastic_norm)[exons_to_keep,]

# Add each gene's TSS
qtltools_input_info_columns <- 
  qtltools_input_info_columns %>%
  left_join(TSSs) %>% 
  mutate(end = TSS+1) %>% 
  mutate(strand = str_split(Exon_coord, "_")) %>%
  mutate(strand = map_chr(strand, 4)) %>% 
  .[,c(4, 5, 6, 1, 2, 7)] %>% 
  set_colnames(c("#chr", "start", "end", "exon", "geneID", "strand")) 
qtltools_input <- 
  cbind(qtltools_input_info_columns, qtltools_input)

#Sanity check
sum(!(rownames(qtltools_input) == qtltools_input$exon))

# Sort the bed file for tabix indexing
qtltools_input <- 
  qtltools_input %>% 
  arrange(`#chr`, start)

qtltools_out_path <- 
  sprintf("bed/%s_normalized_psi_grouped_exons_qtltools.bed", tissue)
write_tsv(qtltools_input, qtltools_out_path, na = "")
system(paste("/nfs/sw/htslib/htslib-1.1/bgzip -f", qtltools_out_path))
system(paste0("/nfs/sw/htslib/htslib-1.1/tabix -p bed -f ", qtltools_out_path, ".gz"))

##### Save another file for a standard non-grouped pass using QTLtools ########
# I want to see if by considering high quality exons individually, without the 
# grouped permuatation thing, how it will effect the number of interpretable
# sQTLs from PSI scores. 

# qtltools_input <- t(stochastic_norm)[exons_to_keep,]
# qtltools_input <- apply(qtltools_input, 2, round, 5)
# qtltools_input <- as.data.frame(qtltools_input)
# 
# 
# qtltools_input_info_columns <- 
#   rownames(qtltools_input) %>%
#   sapply(str_split, "_") %>%
#   do.call("rbind", .) %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   mutate(length = as.numeric(V3)-as.numeric(V2)) %>%
#   .[,c(2,3,4,1,6,5)] %>%
#   set_colnames(c("#chr", "start", "end", "gene", "length", "strand")) 
# 
# # Filter out short or super long exons, by eyeballing the cutoffs
# exons_to_keep_by_length <- (abs(qtltools_input_info_columns$length) > 10) &
#   (abs(qtltools_input_info_columns$length) < 1000)
#   
# qtltools_input <- qtltools_input[exons_to_keep_by_length,]
# qtltools_input_info_columns <- qtltools_input_info_columns[exons_to_keep_by_length,]
# 
# #Sanity check
# sum(!(rownames(qtltools_input) == qtltools_input$gene))
# 
# qtltools_input <- cbind(qtltools_input_info_columns, qtltools_input)
# 
# # Save the output
# qtltools_out_path <- 
#   sprintf("bed/%s_normalized_psi_single_exons_qtltools.bed", tissue)
# write_tsv(qtltools_input, qtltools_out_path, na = "")
# system(paste("~/.linuxbrew/bin/bgzip -f", qtltools_out_path))
# system(paste0("~/.linuxbrew/bin/tabix -p bed -f ", qtltools_out_path, ".gz"))
# 
# # Just to double check
# hist(qtltools_input$length, breaks = 50)
