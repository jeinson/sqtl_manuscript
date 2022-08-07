# Plot the colocalization heatmap figure
# The data is saved from fig2A_colocalization_analysis_prelim.R script


library(here)
source(here("myPackages.R"))
source(here("sqtl_manuscript_functions.R"))


n_hits_per_group_plt <- read.csv("data/n_gwas_hits_per_group_plt.tsv", 
                                 sep = "\t")

trait_annotations <- read_tsv(here("data/gwas_metadata.txt"))
trait_key <- deframe(trait_annotations %>% select(Tag, Category))
trait_key <- sort(trait_key)
trait_key <- trait_key[names(trait_key) %in% colnames(n_hits_per_group_plt)]
n_hits_per_group_plt <- n_hits_per_group_plt[,names(trait_key)]
trait_key <- data.frame(GWAS_trait = trait_key)

fig_S3 <- pheatmap(log2(n_hits_per_group_plt + 1),
         show_colnames = F, cluster_cols = F,
         annotation_col = trait_key
)

save_plot("figS3_GWAS_coloc_heatmap.svg", width = 12, height = 2.5)
fig_S3
dev.off()
