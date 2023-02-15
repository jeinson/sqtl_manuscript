# Look at the results in "tissue_by_tissue_cs_overlap.R" to investigate 
# the exons that seem to be involved with an eQTL result as well. 
# I'll tear apart these results to see if there any interesting patterns that
# pop out. 

rm(list = ls())
source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")

tissues <- 
  list.files("credible_set_overlaps/") %>% str_remove("_cs_overlaps.tsv")

overlaps <- lapply(tissues, function(x) read_tsv(paste0("credible_set_overlaps/", x, "_cs_overlaps.tsv")))
names(overlaps) <- tissues
overlaps <- bind_rows(overlaps, .id = "tissue")

overlaps$gene_id <- str_split(overlaps$eCS_id, "_") %>% map_chr(1)

length(unique(overlaps$gene_id))
nrow(overlaps)

# So there are about 600 distinct genes that have a shared eQTL / sQTL causal
# variant. That's actually not terrible, but still not great either. 

# Convert variants and affected exons into a list
overlaps$exons_in_cs <- str_split(overlaps$exons_in_cs, ",")
overlaps$variants <- str_split(overlaps$variants, ",")
overlaps$n_exons_in_cs <- sapply(overlaps$exons_in_cs, length)

gene = "ENSG00000273136"
tissue = "Adipose_Subcutaneous"
variant = "chr1_121011431_C_T"
library(seqminer)
aFC_file = "/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/afc/ct_afcs.union_independent.gtex.v8.txt.gz"
afc_col_names <- unlist(str_split(read_lines(aFC_file, n_max = 1), "\t"))

get_eQTL_aFC <- function(variant, gene, tissue){
  # parse variant id
  y <- unlist(str_split(variant, "_"))
  variant_tbx <- paste0(y[1], ":", y[2], "-", y[2])
  
  afc <- tabix.read.table(tabixFile = aFC_file, 
                          tabixRange = variant_tbx)
  if(nrow(afc) > 0){
    colnames(afc) <- afc_col_names
    afc$gene_id <- remove_trailing_digit(afc$gene_id)
    
    if(gene %in% afc$gene_id) {
      # Return the aFC score from the tissue and gene of interest
      afc[afc$gene_id == gene,tissue]
    } else {
      return("Not_a_top_eSNP")
    }
  } else {
    return("Not_a_top_eSNP")
  }
}


overlaps_un <- unnest(overlaps, cols = c(variants))

# Get everything! ## Warning this take a while
top_aFC <- with(overlaps_un, mapply(get_eQTL_aFC, variants, gene_id, tissue))

# Now process the results
overlaps_un$eqtl_afc <- unlist(top_aFC)
overlaps_un$eqtl_afc <- as.numeric(overlaps_un$eqtl_afc)

# Now add information about exon symmetry
exon_id_map <- readRDS("../../../data/gtex_stuff/gtex_v8_exon_id_map.rds")

exon_is_non_symmetric <- function(x){ 
  # Takes a list of exons as input. 
  # Returns a list of Booleans
  x <- exon_id_map[x]
  y <- unlist(str_split(x, "_"))
  
  # Test if the length of the exon is divisible by 3
  # If it isn't symmetric, return TRUE (1)
  (as.numeric(y[3]) - as.numeric(y[2]) + 1) %% 3 != 0
}

# As an alternative approach, calculate the total length of spliced exons. 
# When there are multiple spliced together, if the total length induces a 
# frameshift, that could also cause a stop gain variant. 
get_exon_length <- function(x){
  # Takes a list of exons as input. 
  # Returns a list of Booleans
  x <- exon_id_map[x]
  y <- unlist(str_split(x, "_"))
  
  # Test if the length of the exon is divisible by 3
  # If it isn't symmetric, return TRUE (1)
  (as.numeric(y[3]) - as.numeric(y[2]) + 1)
}

# Calculate exon length and symmetry on all exons
exon_symmetry <- lapply(overlaps_un$exons_in_cs, function(x) sapply(x, exon_is_non_symmetric))
exon_lengths <- lapply(overlaps_un$exons_in_cs, function(x) sapply(x, get_exon_length))
total_exon_length <- sapply(exon_lengths, sum)

# Add these metrics back to the full dataframe
overlaps_un$exon_symmetry_list <- exon_symmetry
overlaps_un$CS_has_NS_exon  <- as.logical(sapply(overlaps_un$exon_symmetry_list, sum))
overlaps_un$exon_lengths <- exon_lengths
overlaps_un$total_exon_length <- total_exon_length
overlaps_un$total_exon_symmmetry <- (overlaps_un$total_exon_length %% 3) == 0

#### Analysis ####
overlaps_afc <- 
  overlaps_un %>% 
  filter(!is.na(eqtl_afc))

# Select the gene with the higher aFC in cases where the gene is significant
# in multiple tissues. 
overlaps_afc %>%
  group_by(gene_id) %>%
  filter(abs(eqtl_afc) == max(abs(eqtl_afc))) %>%
  ungroup
