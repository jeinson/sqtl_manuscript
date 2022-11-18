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

# How many splicing events colocalize with something?
# Add the Euclidean distance to (1,1)
results %<>%
  mutate(dist = sqrt((1-pp_power)^2+(1-pp_coloc)^2))
results$has_coloc <- results$dist < .25

exons_that_colocalize <- 
  tapply(results$has_coloc, results$phenotype_id, function(x) as.logical(sum(x)))
sum(exons_that_colocalize) / length(exons_that_colocalize)

# Coloc heat map ----
# Make a heatmap showing colocalization events across tissues and traits
n_hits_per_group <-
  results %>%
  filter(has_coloc) %>%
  count(qtl, gwas) %>%
  rename("n" = "n_hits")

n_total_hits <-
  results %>%
  count(qtl, gwas) %>%
  rename("n" = "n_total")

percentage_of_hits_per_tests <-
  inner_join(n_hits_per_group, n_total_hits) %>%
  mutate(percent_of_hits = n_hits / n_total)

# Assign a ceiling for overlap percentage so
outlier_cutoff <- with(percentage_of_hits_per_tests,
                       quantile(percent_of_hits)[4] + 1.5 * IQR(percent_of_hits) )

percentage_of_hits_per_tests$percent_of_hits[percentage_of_hits_per_tests$percent_of_hits > outlier_cutoff] <- outlier_cutoff

tiss_order <- tapply(percentage_of_hits_per_tests$percent_of_hits, percentage_of_hits_per_tests$qtl, sum) %>% sort() %>% names
trait_order <- tapply(percentage_of_hits_per_tests$percent_of_hits, percentage_of_hits_per_tests$gwas, sum) %>% sort() %>% names

percentage_of_hits_per_tests$qtl <- factor(percentage_of_hits_per_tests$qtl,
                                           levels = tiss_order)
percentage_of_hits_per_tests$gwas <- factor(percentage_of_hits_per_tests$gwas,
                                            levels = trait_order)

library(pheatmap)
n_hits_per_group_plt <-
  spread(n_hits_per_group, key = "gwas", value = "n_hits") %>%
  as.data.frame %>%
  column_to_rownames("qtl")
n_hits_per_group_plt[is.na(n_hits_per_group_plt)] <- 0

trait_annotations <- read_tsv("/gpfs/commons/groups/lappalainen_lab/dglinos/projects/GTEx_eVariants/gwas_metadata.txt")
trait_key <- deframe(trait_annotations %>% select(Tag, Category))
trait_key <- sort(trait_key)
trait_key <- trait_key[names(trait_key) %in% colnames(n_hits_per_group_plt)]
n_hits_per_group_plt <- n_hits_per_group_plt[,names(trait_key)]

# Save this for plotting in a separate script
write.table(n_hits_per_group_plt, "data/n_gwas_hits_per_group_plt.tsv", sep = "\t")
trait_key <- data.frame(GWAS_trait = trait_key)
pheatmap(log2(n_hits_per_group_plt + 1),
         show_colnames = F, cluster_cols = F,
         annotation_col = trait_key
)

# Filter down to top genes per trait ----
# Only include genes that are in the top sQTL set
# Include the terminal exon filter. 
top_sQTLs <- read_tsv(here("data/top_sQTLs_MAF05.tsv"))
n_exons_per_gene <- read_tsv(here("data/gtex_v8_n_exons_per_gene.tsv"))
top_sQTLs <- 
  top_sQTLs %>%
  mutate(exon_number = str_split(top_pid, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("group" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons))
top_sQTLs %<>% filter(!terminal_exon_flag)

results_top_coloc <- 
  results %>%
  filter(phenotype_id %in% top_sQTLs$top_pid) %>%
  group_by(GeneID) %>%
  filter(pp_power + pp_coloc == max(pp_power + pp_coloc)) %>%
  dplyr::slice(1) %>%
  #slice(1) %>%
  ungroup

# Test what Euclidean distance optimizes the number of things we call colocalized
pct <- function(x) sum(x) / length(x)
x = seq(0, 1.2, by = .05)
y = sapply(x , function(x) pct(results_top_coloc$dist < x))
plot(x, y, type = 'o')

# Save this file for downstream analyses. I think this will make our lives easier. 
write_tsv(results_top_coloc, here("data/top_sQTLs_with_top_coloc_event.tsv"))


# Retrieve the derived and ancestral alleles for the top colocalizing alleles
# Let's just do this super inefficiently but whatever. 
x_subset <- results_top_coloc %>%
  select(qtl_lead_snv, qtl, phenotype_id, gwas)


# For each row, grab the qtl beta value and the derived allele. 
x <- str_split(x_subset$qtl_lead_snv, "_")
x_subset$chr <- map_chr(x, 1)
x_subset$pos <- as.numeric(map_chr(x, 2))
x_subset$ref <- map_chr(x, 3)
x_subset$alt <- map_chr(x, 4)

x_subset$anc_allele <- character(length = nrow(x_subset))
x_subset$beta <- numeric(length = nrow(x_subset))

get_qtl_path <- function(tiss, chr){
  sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/qtltools_results/%s_nominal/%s.qtltools.nominal.MAF05.%s.06.02.21.txt.gz", 
          tiss, 
          tiss,
          chr)
}

library(seqminer)
tabix.read.table.custom <- 
  function (tabixFile, tabixRange, col.names = TRUE, stringsAsFactors = FALSE) {
    stopifnot(seqminer:::local.file.exists(tabixFile))
    stopifnot(all(isTabixRange(tabixRange)))
    tabixFile <- path.expand(tabixFile)
    storage.mode(tabixFile) <- "character"
    storage.mode(tabixRange) <- "character"
    header <- .Call("readTabixHeader", tabixFile, PACKAGE = "seqminer")
    body <- .Call("readTabixByRange", tabixFile, tabixRange, 
                  PACKAGE = "seqminer")
    body <- do.call(rbind, strsplit(body, "\t"))
    body <- as.data.frame(body, stringsAsFactors = FALSE)
    # Don't do the type conversion. that screw everything
    if (ncol(body) > 0) {
      #   for (i in 1:ncol(body)) {
      #     body[, i] <- utils::type.convert(body[, i], as.is = !stringsAsFactors)
      #   }
      num.col <- ncol(body)
      header <- header[nchar(header) > 0]
      if (length(header) == 0 || !col.names) {
        colNames <- paste0("V", 1L:num.col)
      }
      else {
        hdrLine <- header[length(header)]
        hdrLine <- sub("^#", "", hdrLine)
        colNames <- make.names(strsplit(hdrLine, "\t")[[1]])
        if (length(colNames) > ncol(body)) {
          colNames <- colNames[1:ncol(body)]
        }
        else if (length(colNames) < ncol(body)) {
          tmpNames <- paste0("V", 1L:num.col)
          tmpNames[1:length(colNames)] <- colNames
          colNames <- tmpNames
        }
      }
      colnames(body) <- colNames
    }
    body
  }

anc_allele_file = "/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/analysis/release_91_homo_sapiens.txt.gz"
get_ancestral_allele <- function(chr, pos, ref, alt){
  chr_abs <- str_remove(chr, "chr")
  trange <- paste0(chr_abs, ":", pos, "-", pos)
  tabix.result <- tabix.read.table.custom(anc_allele_file, trange)
  
  if(nrow(tabix.result) > 0){
    anc_allele <- tabix.result$V7[tabix.result$V4 == ref & tabix.result$V5 == alt]
  } else {
    anc_allele <- NA
  }

  anc_allele
}

qtl_path = ""

# Sort to avoid having to read in the data frame many many times
x_subset %<>% arrange(qtl, chr)

# Do the whole shebang
for(i in 1:nrow(x_subset)){
  attach(x_subset[i,], warn.conflicts = F)
  message(paste0(i, ", "), appendLF = F)
  qtl_path_new <- get_qtl_path(qtl, chr)
  
  # Prevent having to read in the file every time I do this
  if(qtl_path != qtl_path_new){
    qtl_path <- qtl_path_new
    qtl_res <- read_delim(qtl_path, delim = " ", col_names = F)
  }
  
  anc_allele <- get_ancestral_allele(chr, pos, ref, alt)
  qtl_beta <- filter(qtl_res, X1 == phenotype_id & X10 == pos)$X13
  
  x_subset$anc_allele[i] <- anc_allele
  x_subset$beta[i] <- qtl_beta
  
}

write_tsv(x_subset, here("data/top_sQTLs_with_top_coloc_event_aa_and_beta.tsv"))
