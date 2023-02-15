# Perform some analyses to help me understand what susie is actually doing. 
# I'm still confused about what genes get a credible set. 
# 
# I think Susie is run on all genes that have a top sQTL. But are there are 
# many genes that have a top-sQTL that don't have any credible set in the 
# final file. Are these filtered out for some reason? Or am I missing something?


setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")
source("~/myPackages.R")

# Work with one tissue at a time here
tiss = "Nerve_Tibial"

# Read in CS
all_cs <- readRDS("all_credible_sets_12_2.rds")
cred_set <- all_cs[[tiss]]

# and sQTL calls
sqtls <- read_tsv(sprintf(
  "../sQTL_v8_anno/combined_qtltools_results/%s_combined_sQTLs.tsv", tiss)
)

sqtls <- filter(sqtls, bpval < .05)

# Are all the genes in here at least?
remove_trailing_digit(unique(cred_set$gene_id)) %in% sqtls$group %>% not %>% sum

# Ok well that's good. How about the other way around?
sqtls$group %in% remove_trailing_digit(unique(cred_set$gene_id)) %>% table

# about 2/3 are there. But happened to the other 1/3?
missing_sqtls <- 
  sqtls[!(sqtls$group %in% remove_trailing_digit(unique(cred_set$gene_id))), ]

# Ah, ok. it looks like the susie pipeline is applying FDR correction to 
# the beta QTL p-values, and only using the ones that are significant at this 
# threshold. I think I should probably be doing this too.... 
# 
# I don't like just throwing FDR at something because "that's what you're 
# supposed to do." But I don't have any real reason not to do it either. 