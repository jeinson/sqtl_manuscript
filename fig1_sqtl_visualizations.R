# This script just makes a few summary plots for visualizing sQTLs mapped 
# across 18 tissues.
# 
# Wrapped into Figure 1 of the Manuscript
rm(list = ls())
source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno")
source("../sqtl_manuscript/sqtl_manuscript_functions.R")
library(ggplot2)
sqtls <- read_tsv("cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")

# Make a barplot of the number of times a gene is significantg in a tissue. 
# svg("sqtl_analyses/tissue_coverage_plot.svg", width = 3, height = 2)
# ggplot(sqtls, aes(n_tiss)) + 
#   geom_bar(fill = "cornflowerblue") + 
#   xlab("N tissues where a gene is significant") + 
#   ylab("N genes") +
#   gtex_v8_figure_theme()
# dev.off()

# Find the number of genes tested vs the number significant. 
tissues <- read_lines("completed_tissues.txt")
qtl_lists <- tissues %>% map(~ read_tsv(paste0("combined_qtltools_results/", 
                                               .x, "_combined_sQTLs.tsv")))
names(qtl_lists) <- tissues

sig_genes <- sapply(qtl_lists, function(x) sum(x$bpval < .05))
tested_genes <- sapply(qtl_lists, nrow)
percent_sig <- sig_genes / tested_genes
n_per_tiss <- sapply(tissues, function(x) length(read_lines(paste0("indvs_per_tiss/", x, "_indvs.txt"))))

summary_table <- data.frame(
  "Tissue" = tissues, 
  "Tested Genes" = tested_genes, 
  "Significant Genes" = sig_genes, 
  "Percent Significant" = percent_sig, 
  "Total Individuals" = n_per_tiss, 
  check.names = F
)

summary_table$Tissue <- tissues

library(ggrepel)
save_plot("fig1_sqtl_totals_summary.svg", width = 5, height = 5)
ggplot(summary_table, aes(`Total Individuals`, `Significant Genes`)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_label_repel(aes(label = Tissue), 
                   box.padding = .4,
                   point.padding = .7,
                   size = 3) + 
  sqtl_manuscript_theme()
dev.off()
 summary(lm(`Significant Genes`~`Total Individuals`, data = summary_table))

# Do a similar thing for leafcutter sQTLs
lc_root = "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/sqtl/GTEx_Analysis_v8_sQTL/"
leafcutter_lists <- tissues %>% map(~ read_tsv(paste0(lc_root, .x, ".v8.sgenes.txt.gz")))
names(leafcutter_lists) <- tissues

qtl_lists <- leafcutter_lists
sig_genes <- sapply(qtl_lists, function(x) sum(x$pval_beta < .05))
tested_genes <- sapply(qtl_lists, nrow)
percent_sig <- sig_genes / tested_genes

summary_table_lc <- data.frame(
  "Tissue" = tissues, 
  "Tested Genes" = tested_genes, 
  "Significant Genes" = sig_genes, 
  "Percent Significant" = percent_sig, 
  "Total Individuals" = n_per_tiss, 
  check.names = F
)

summary_table_lc$Tissue <- tissues
# Plot them side by side

save_plot("fig1S_leafcutter_psi_pct_sig_comparison.svg", width = 8, height = 4)
list(PSI = select(summary_table, Tissue, `Percent Significant`), 
     Leafcutter = select(summary_table_lc, Tissue, `Percent Significant`)
) %>% bind_rows(.id = "Splicing quantification method") %>%
  ggplot(aes(Tissue, `Percent Significant`, fill = `Splicing quantification method`)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  sqtl_manuscript_theme() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
dev.off()

# What does the distribution of higher inclusion sQTL allele frequencies 
# look like?

save_plot("fig1_high_inclusion_allele_frequency_distribution.svg",
    width = 4, height = 3)
ggplot(sqtls, aes(sQTL_af)) + 
  geom_histogram(fill = "cornflowerblue") + 
  xlab("Higher inclusion sQTL allele frequency") + 
  ggtitle("Distribution of high exon inclusion sQTL allele frequencies") +
  gtex_v8_figure_theme()
dev.off()

