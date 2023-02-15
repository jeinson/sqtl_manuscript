# Get aFC of eQTLs that share a causal variant with a PSI-sQTL. 
# 
# This script does the same thing as the first part of "analyze_overlap_signals.R"
# but uses the aFC stored on the exchange, which Elise recommended. 
#
# The processing doesn't take too long and output format would be a nested list, 
# so let's just do the analysis in this script. 

#### COPY PASTE FROM OLD SCRIPT ####
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

overlaps_un <- unnest(overlaps, cols = c(variants), )
## END COPY PASTE ##

#### Function for retrieving aFC ####
aFC_base = "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_aFC_top_eQTL/"
aFC_bytiss <- list()

# Read in all aFC data first
for(tiss in tissues){
  aFC_bytiss[[tiss]] <- 
    read_tsv(paste0(aFC_base, tiss, ".aFC.txt.gz")) %>%
    mutate(pid = remove_trailing_digit(pid)) %>%
    mutate(sid = str_remove(sid, "_b38"))
}

# Define the function
# variants = list(), gene = 'character', tissue = 'character'
get_eQTL_aFC <- function(variants, gene, tissue) {
  x = filter(aFC_bytiss[[tissue]], pid == gene & sid %in% variants)
  if(nrow(x) > 0){
    x$log2_aFC
  } else {
    NA
  }
}

# Run the function for all overlapping CSs
overlap_aFC <- with(overlaps, 
                    mapply(get_eQTL_aFC, variants, gene_id, tissue, SIMPLIFY = F))
overlap_aFC <- unlist(overlap_aFC)

# Add the result back to the list
overlaps$eQTL_log2_aFC <- overlap_aFC

#### Analysis ####
# Filter down the overlap list to include each gene once. Select the tissue with
# the top aFC to represent the gene
nrow(overlaps)

overlaps_og <- 
  overlaps %>%
  group_by(gene_id) %>%
  arrange(desc(abs(eQTL_log2_aFC))) %>%
  #slice(1) %>% # incase of ties or all NA
  filter(row_number() == 1) %>%
  ungroup

barplot(table(!is.na(overlaps_og$eQTL_log2_aFC)))

 # What does the distribution of aFC look like?
boxplot(overlaps_og$eQTL_log2_aFC, horizontal = T, 
        main = "log2(aFC) of eQTLs that overlap an sQTL")
abline(v = 0, lty = 2, col = "red")


# Add symmetry information to the overlapping CSs
overlaps_afc <- filter(overlaps_og, !is.na(eQTL_log2_aFC))
exon_id_map <- read_rds("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/gtex_v8_exon_id_map.rds")

get_exon_lengths <- function(x){
  # Takes a list of exons as input. 
  # Returns a list of Booleans
  x <- exon_id_map[x]
  y <- str_split(x, "_")
  
  end <- as.integer(unlist(map(y, 3)))
  start <- as.integer(unlist(map(y, 2)))
  
  end - start + 1
}

exon_lengths <- map(overlaps_afc$exons_in_cs, get_exon_lengths)

overlaps_afc$exon_lengths <- exon_lengths
overlaps_afc$total_exon_length <- map_dbl(overlaps_afc$exon_lengths, sum)
overlaps_afc$exon_mod3 <- map(exon_lengths, ~ magrittr::mod(.x, 3))

# Look at aFC outliers
overlaps_afc$aFC_outlier <- flag_outliers(overlaps_afc$eQTL_log2_aFC)
