# Perform summary analyses on SuSiE psi-sQTL fine mapped credible sets. 
# 
# Also do some comparison analyses with eQTL fine mapped variants. How much
# overlap do we see? Does this correspond with an sQTL including non-symmetric
# exon? 

rm(list = ls())
source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")

# # Read in data
# fps <- list.files("sqtl_v8_anno_susie_output", pattern = "22.txt.gz", full.names = T)
# tissues <- read_lines("../sQTL_v8_anno/completed_tissues.txt")
# # Read and process everything in one go!
# tissue_level_cs <- lapply(tissues[-8], function(x) {
#   fps <- map_dfr(grep(x, fps, value = T), read_tsv)
#   out <- fps[!duplicated(fps),] %>%
#     mutate(gene_id = str_remove(phenotype_id, "_.+")) %>%
#     relocate(gene_id, .before = phenotype_id) %>%
#     mutate(variant_id = variant_id %>% str_remove("chr"))
#   out$tissue <- x
#   out %>% relocate(tissue, .before = gene_id)
# })
# 
# # # 
# names(tissue_level_cs) <- tissues[-8]
# write_rds(tissue_level_cs, "all_credible_sets_12_2.rds")

if(FALSE) tissue_level_cs$Artery_Tibial %>% View

# Just read in the RDS file
tissue_level_cs <- read_rds("all_credible_sets_12_2.rds")

# Define a function to apply credible set filtering, on a tissue by tissue
# basis
# x: Tissue
# y: Gene
# z: 
# 
# Thanks Paul for figuring out how to do this. 
collapse <- function(x) {
  # Perform unions if there's an intersection
  y <- lapply(
    X = seq_along(along.with = x),
    FUN = function(i) {
      return(Reduce(
        f = function(a, b) {
          if (length(x = intersect(x = a, y = b))) {
            return(union(x = a, y = b))
          }
          return(a)
        }, 
        x = x[seq.int(from = i, to = length(x = x))]
      ))
    }
  )
  # Remove values that had been added previously
  ni <- vector()
  for (i in seq_along(along.with = y)) {
    if (i == 1L) {
      next
    }
    for (j in seq.int(from = 1L, to = i - 1L)) {
      if (all(x[[i]] %in% x[[j]])) {
        ni <- c(ni, i)
        break
      }
    }
  }
  if (length(x = ni)) {
    y <- y[-ni]
  }
  # Run collapse again
  if (!identical(x = y, y = x)) {
    return(collapse(x = y))
  }
  return(y)
}

# Collapse everything down
tissue_level_collapsed_cs <- 
  lapply(tissue_level_cs, function(x){
    y <- x %>% split(.$gene_id)
    for(i in seq_along(y)){
      
      # Find overlaps between all credible sets to get a finalized list 
      # "subgraphs" in this network of credible sets 
      cred_sets <- split(y[[i]]$variant_id, y[[i]]$cs_id)
      collapsed_cred_sets <- collapse(cred_sets)
      names(collapsed_cred_sets) <- paste0("FCS", seq_along(collapsed_cred_sets))
      
      # 5/3/22
      # I think 'FCS' was supposed to stand for Final Credible Set, though
      # I'm not sure how final this really is... 
      
      # Now assign the variants to the final credible sets. 
      cs_key <- stack(collapsed_cred_sets)
      cs_key$ind <- paste0(y[[i]]$gene_id[1], "_", cs_key$ind)
      cs_key <- deframe(cs_key)
      y[[i]]$collapsed_cs_id <- cs_key[y[[i]]$variant_id]
    }
    bind_rows(y)
  })

collapsed_credible_sets <- bind_rows(tissue_level_collapsed_cs, .id = "tissue")

# Now collapse credible sets again across tissues
more_collapse_credible_sets <- 
  collapsed_credible_sets %>% 
  split(.$gene_id) %>%
  lapply(function(x){
    collapsed_cred_sets_2 <- collapse(split(x$variant_id, x$collapsed_cs_id))
    names(collapsed_cred_sets_2) <- paste0("CS", seq_along(collapsed_cred_sets_2))
    
    # Now assign the variants to the final credible sets. 
    cs_key <- stack(collapsed_cred_sets_2)
    cs_key$ind <- paste0(x$gene_id[1], "_", cs_key$ind)
    cs_key <- deframe(cs_key)
    x$collapsed_cs_id_cross_tiss <- cs_key[x$variant_id]
    
    x
  })

more_collapse_credible_sets <- bind_rows(more_collapse_credible_sets)

# Just for my own understanding, save this as a VCF for visualization in IGV
var_list <- select(more_collapse_credible_sets, variant_id, gene_id, collapsed_cs_id_cross_tiss) %>%
  distinct
var_df <- data.frame(
  do.call("rbind", str_split(var_list$variant_id, "_")))
colnames(var_df) <- c("#CHROM", "POS", "REF", "ALT")
# Fill in missing info with placeholders
var_df$ID <- "."

var_vcf <- var_df[,c("#CHROM", "POS", "ID", "REF", "ALT")]

var_vcf$QUAL <- "."
var_vcf$FILTER <- "."
var_vcf$INFO <- 
  paste(
    paste0("TG=", var_list$gene_id),
    #paste0("TISS=", var_list$tissue),
    paste0("CS_ID=", var_list$collapsed_cs_id_cross_tiss), 
    sep = ";"
  )

# Sort for indexing
var_vcf <- var_vcf[with(var_vcf, order(nchar(`#CHROM`), `#CHROM`, nchar(POS), POS)),]

# write this thing
out_name = "gtex_v8_sqtl_collapsed_credible_sets.vcf"
system(paste("touch", out_name))
write_lines("##fileformat=VCFv4.2", out_name, append = T)
write_lines('##INFO=<ID=TG,Number=A,Type=String,Description="Target gene of the sQTL variant">', out_name, append = T)
write_lines('##INFO=<ID=TISS,Number=A,Type=String,Description="Tissue where the sQTL is active">', out_name, append = T)
write_lines('##INFO=<ID=CS_ID,Number=A,Type=String,Description="The collapsed credible set ID (assigned by Jonah)">', out_name, append = T)
write_tsv(var_vcf, out_name, append = T, col_names = T)
system(paste("/nfs/sw/htslib/htslib-1.9/bin/bgzip -f", out_name))
#system("/nfs/sw/htslib/htslib-1.9/bin/tabix -p vcf topmed/topmed_csnps_v8_anno.vcf.gz")

# Try plotting some GTEx GE Css
ge_cs <- read_tsv("GTEx_ge_cred_sets/GTEx_ge_adipose_subcutaneous.purity_filtered.txt.gz")
x <- filter(ge_cs, molecular_trait_id == unique(molecular_trait_id)[1])
