# Make Figure 1 diagrams for sQTL Manuscript. 
# 
# This combines pieces of code from other scripts, only taking what will
# actually be used downstream. 
# 
# Jonah Einson
# 4/12/22
library(here)
source(here("myPackages.R"))
source(here("sqtl_manuscript_functions.R"))
library(ggplot2)


#### 1B. N Significant vs. N individuals in GTEx ####
# Find the number of genes tested vs the number significant. 
tissues <- read_lines(here("data/completed_tissues.txt"))
qtl_lists <- tissues %>% map(~ read_tsv(paste0(here("../sQTL_v8_anno/combined_qtltools_results/"), 
                                               .x, "_combined_sQTLs.tsv")))
names(qtl_lists) <- tissues

sig_genes <- sapply(qtl_lists, function(x) sum(x$bpval < .05))
tested_genes <- sapply(qtl_lists, nrow)
percent_sig <- sig_genes / tested_genes
n_per_tiss <- sapply(tissues, function(x) length(read_lines(paste0(here("data/indvs_per_tiss/"), x, "_indvs.txt"))))

summary_table <- data.frame(
  "Tissue" = tissues, 
  "Tested Genes" = tested_genes, 
  "Significant Genes" = sig_genes, 
  "Percent Significant" = percent_sig, 
  "Total Individuals" = n_per_tiss, 
  check.names = F
)

summary_table$Tissue <- tissues

# Supplementary figure, indvs to genes
library(ggrepel)
figS1A <- 
  ggplot(summary_table, aes(`Total Individuals`, `Significant Genes`)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_label_repel(aes(label = Tissue), 
                   box.padding = .4,
                   point.padding = .7,
                   size = 3) + 
  sqtl_manuscript_theme()
summary(
  lm(`Significant Genes` ~ `Total Individuals`, data = summary_table)
)

# Make a barplot of the number of significant sQTLs mapped
# per tissue. Keep it simple. 
gtex_colors <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/gtex_colors_alt.txt")
summary_table_w_color <- 
  summary_table %>%
  left_join(gtex_colors %>% select(tissue_id, tissue_site_detail, color_hex), 
            by = c("Tissue" = "tissue_id")) %>%
  mutate(color_hex = paste0("#", color_hex)) %>%
 # arrange(`Significant Genes`) %>%
  mutate(tissue_site_detail = as_factor(tissue_site_detail))
fig1B <- 
ggplot(summary_table_w_color, aes(tissue_site_detail, `Significant Genes`, fill = Tissue)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = summary_table_w_color$color_hex) + 
  sqtl_manuscript_theme() + 
  theme(legend.position = "none") + 
  coord_flip()

### Do the same thing with Leafcutter sQTLs for a comparison plot. 
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

# plot the plot
figS1F <- list(PSI = select(summary_table, Tissue, `Percent Significant`),
     Leafcutter = select(summary_table_lc, Tissue, `Percent Significant`)
) %>% bind_rows(.id = "Splicing quantification method") %>%
  ggplot(aes(Tissue, `Percent Significant`, fill = `Splicing quantification method`)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))


# All sQTLs
sqtls <- read_tsv(here("data/top_sQTLs_MAF05.tsv"))

# Remove terminal exons
n_exons_per_gene <- read_tsv(here("data/gtex_v8_n_exons_per_gene.tsv"))
sqtls %<>%
  mutate(exon_number = str_split(top_pid, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("group" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons))
sqtls %<>% filter(!terminal_exon_flag)

## Fig 1S1: sQTL length comparison ----
# Calculate the lengths of exons
exon_id_map <- read_rds(here("data/gtex_v8_exon_id_map.rds"))
sqtls$pid <- exon_id_map[sqtls$top_pid]

calculate_exon_length <- function(x){
  k = unlist(str_split(x, "_"))
  start = as.integer(k[2])
  end = as.integer(k[3])
  end - start + 1
}


sqtls$exon_length <- sapply(
  sqtls$pid, calculate_exon_length
)

hist(log10(sqtls$exon_length), col = "cornflowerblue", 
     main = "Distribution of sExon lengths", 
     xlab = "log10(Exon Length) (bp)")
lim <- quantile(sqtls$exon_length)[4] + (1.5 * IQR(sqtls$exon_length))
abline(v = log10(lim), col = "red", lty = 2)

# Also grab non-sQTLS
non_sqtls <- read_tsv(here("data/cross_tissue_nonsignificant_genes.tsv")) %>%
  mutate(pid = exon_id_map[top_pid]) %>%
  mutate(exon_number = str_split(top_pid, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("group" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons)) %>%
  filter(!terminal_exon_flag) 

non_sqtls$exon_length <- sapply(non_sqtls$pid, calculate_exon_length)

exon_lengths <- 
  rbind(
    non_sqtls %>% 
      select(exon_length) %>%
      mutate(tag = "Variable exons\nwith no sQTL") %>%
      mutate(outlier = flag_outliers(exon_length)),
    sqtls %>% 
      select(exon_length) %>% 
      mutate(tag = "sExons") %>%
      mutate(outlier = flag_outliers(exon_length))
  )
exon_lengths$tag <- as_factor(exon_lengths$tag)

fig_S1B <- ggplot(exon_lengths %>% filter(!outlier), aes(tag, exon_length)) + 
  geom_violin(aes(fill = tag)) + 
  scale_fill_brewer(palette = "Set3") + 
  geom_boxplot(width = .2)  + 
  sqtl_manuscript_theme() + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("Exon Length (bp)")

wilcox.test(exon_length ~ tag, data = exon_lengths %>% filter(!outlier))

#### 1C. sExon symmetry vs. all exons ####
sExon_cds_location <- sqtls$exon_number / sqtls$n_exons
non_sExon_cds_location <- non_sqtls$exon_number / non_sqtls$n_exons

# fig_S1D <- 
#   ggplot() + 
#   geom_histogram(aes(x = sExon_cds_location), bins = 10, 
#                  color = "black", fill = "darkolivegreen3") + 
#   xlab("Normalized location in CDS") + 
#   ylab("sExons") + 
#   sqtl_manuscript_theme()
# library(spgs)
chisq.unif.test(sExon_cds_location)

exon_cds_locations <- 
  rbind(
    enframe(sExon_cds_location) %>% mutate(tag = 'sExons'),
    enframe(non_sExon_cds_location) %>% mutate(tag = 'non sExons')
  )

fig_S1C<- 
  ggplot(exon_cds_locations, aes(value, fill = tag)) + 
  geom_density(alpha = .5) + 
  scale_fill_brewer(palette = "Set3") + 
  xlab("Exon's relative position in transcript") +
  sqtl_manuscript_theme()

ks.test(sExon_cds_location, non_sExon_cds_location)
wilcox.test(sExon_cds_location, non_sExon_cds_location)

### Save Figure S2
# library(cowplot)
# save_plot("figS2_sexon_length_and_positions.svg", width = 7, height = 3)
# plot_grid(fig_S2A, fig_S2B, align = "hv", labels = "AUTO", 
#           rel_widths = c(.4, .6))
# dev.off()

#### All exon symmetry
# Get information about ALL exons
library(rtracklayer)
gtf <- rtracklayer::import("../data/gencode.v26.GRCh38.GTExV8.genes.gtf")
gtf_df <- data.frame(gtf)
all_exons <- filter(gtf_df, transcript_type == "protein_coding" & type == "exon") %>%
  select(gene_id, transcript_id, exon_id, seqnames, start, end, strand, width)

all_exons$exon_number <- 
  str_split(all_exons$exon_id, " ") %>%
  map_chr(3) %>%
  as.integer()

# Add a BOOL column for being a terminal exon
all_exons <- 
  all_exons %>%
  group_by(gene_id) %>%
  mutate(terminal = (exon_number == 1 | exon_number == n())) %>%
  ungroup

# Remove terminal exons
all_exons <- filter(all_exons, !terminal)

# Check symmetry of both groups, run a fisher test, and make table for plotting
sqtls$is_symmetric <- (sqtls$exon_length %% 3 == 0)
all_exons$is_symmetric <- (all_exons$width %% 3 == 0)

is.symmetric.c.table <- rbind(
  table(all_exons$is_symmetric),
  table(sqtls$is_symmetric)
)

rownames(is.symmetric.c.table) <- c("All Exons", "sExons")
fisher.test(is.symmetric.c.table)

c.table <- rbind(table(all_exons$width %% 3), table(sqtls$exon_length %% 3))
rownames(c.table) <- c("All Exons", "sExons")

p.table <- t(apply(c.table, 1, function(x) x / sum(x)))

fig1C <- 
  as.data.frame(p.table) %>% 
  rownames_to_column("group") %>% 
  pivot_longer(cols = c('0', '1', '2'), names_to = "Exon_length", values_to = "Percent") %>%
  ggplot(aes(Exon_length, Percent, fill = group)) + 
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) + 
  geom_bar(stat = "identity", position = "dodge") + 
  sqtl_manuscript_theme() + 
  theme(legend.position = c(.8, .9)) + 
  xlab("Exon Length mod 3") + ylab("Percent of All Exons")

#### 1D High Inclusion Allele Frequency Histogram ####
# This one should be quite easy. 
fig1D <- 
  ggplot(sqtls, aes(sQTL_af)) + 
  geom_histogram(fill = "cornflowerblue") + 
  xlab("Higher inclusion sQTL allele frequency") + 
  sqtl_manuscript_theme()

#### 1E Derived allele effect direction bars ####
sqtls_da <- read_tsv(here("data/top_sQTLs_MAF05_w_anc_allele.tsv"))
sqtls_da <- filter(sqtls_da, group %in% sqtls$group)

# Add if derived is lower or higher allele
sqtls_da$der_allele_hl <- 
  with(sqtls_da, 
       # This is based on the ancestral allele, so flip it since we're interested
       # in the ancestral allele
       ifelse(anc_allele == hi_allele,
              "Decreases\nPSI", 
              ifelse(anc_allele == li_allele, 
                     "Increases\nPSI", NA
              )
       )
  )

# Save the results of this for downstream analysis. 
sqtls_da %>%
  select(group:top_pid, vid, anc_allele, slope,
         anc_allele_freq, hi_allele, li_allele,
         der_allele_hl) %>%
  write_tsv(here("data/top_sQTLs_MAF05_derived_allele_effect_directions.tsv"))


fig1E <- 
  ggplot(sqtls_da, aes(der_allele_hl, fill = der_allele_hl)) + 
  geom_bar() + 
  xlab("Derived allele PSI effect") + 
  ylab("N Genes") + 
  theme(legend.position = "none") + 
  sqtl_manuscript_theme()

#### 1F Derived allele histogram ####
sqtls_da$der_allele_freq <- 1 - sqtls_da$anc_allele_freq
sqtl_subset <- subset(sqtls_da, !is.na(der_allele_hl) & delta_psi > .035) 

fig1F <- 
  ggplot(sqtl_subset, aes(der_allele_freq, fill = der_allele_hl)) +
  geom_density(alpha = .5) + 
  labs(fill = "Derived allele\neffect\n(ΔPSI > 0.035)") + 
  xlab("Derived Allele Frequency") + 
  ylab("Density") + 
  sqtl_manuscript_theme() + 
  theme(legend.position = c(.8, .8))

KS_test <- sqtl_subset %>% 
  split(.$der_allele_hl) %>%
  ks.test(x = .[['Increases PSI']]$der_allele_freq, y = .[['Decreases PSI']]$der_allele_freq) 

table(sqtl_subset$der_allele_hl)

# 1G Difference in PSI distributions across dPSI cutoffs ----
# Do this across several cutoffs to plot. 
delta_psi <- seq(0, .45, by = .05)
KS_test_results <- 
  sapply(X = delta_psi, FUN = function(x) {
  subset(sqtls_da, !is.na(der_allele_hl) & delta_psi > x) %>% 
      split(.$der_allele_hl) %>%
      ks.test(x = .[['Increases\nPSI']]$der_allele_freq, y = .[['Decreases\nPSI']]$der_allele_freq)
  }) %>%
  t %>% as.data.frame
KS_test_results$delta_psi <- delta_psi       

figS1D <- 
  ggplot(KS_test_results, aes(delta_psi, unlist(statistic))) + 
  geom_point() + 
  geom_line() + 
  xlab("Delta PSI Cutoff") + 
  ylab("D statistic") + 
  ylim(0, .5) +
  sqtl_manuscript_theme() + 
  theme(panel.grid.major = element_line(colour = "grey", size = .5, linetype = 1), 
        panel.grid.minor = element_line(colour = "grey", size = .25, linetype = 2))

# 1H Annotions of high vs low inclusion sQTLs ----
# This was a point brought up. What do the annotations look like among different
# classes of top sQTLs?

# Get the variants in VEP format for a VEP run
# sqtls_vep <- sqtls_da %>% select(chr, start, end)
# sqtls_vep$chr <- str_remove(sqtls_vep$chr, "chr")
# sqtls_vep$allele <- paste0(sqtls_da$ref_allele, "/", sqtls_da$alt_allele)
sqtls_vep <- 
  sqtls_da$vid %>% 
  str_split("_") %>% 
  map(as.data.frame.character) %>% 
  bind_cols() %>% t %>% 
  as_tibble()

sqtls_vep %<>% 
  mutate(V6 = paste0(V3, "/", V4)) %>%
  select(V1, V2, V6) %>%
  mutate(V3 = V2) %>%
  relocate(V3, .before = V6) %>%
  mutate(V7 = sqtls_da$strand) %>%
  arrange(nchar(V1), V1, as.numeric(V2)) %>%
  mutate(V1 = str_remove(V1, "chr"))

# Save this for running VEP
write_tsv(sqtls_vep, "cross_tissue_top_sQTLs/top_sQTLs_MAF05.var", col_names = F)

### After running:
library(vcfR)
vcf_fp = "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/cross_tissue_top_sQTLs/top_sQTLs_MAF05.annotated.vcf"
sqtls_annotated <- read.vcfR(vcf_fp)
sqtls_annotated <- vcfR::extract_info_tidy(sqtls_annotated) %>%
  mutate(genes = str_split(genes, ",")) %>%
  mutate(annotations = str_split(annotations, ",")) %>%
  mutate(transcripts = str_split(transcripts, ",")) %>%
  mutate(distances = str_split(distances, ",")) %>% 
  unnest(cols = c(genes, transcripts, distances, annotations)) %>%
  mutate(transcripts = remove_trailing_digit(transcripts))

sqtls_annotated_subset <- sqtls_da %>% 
  left_join(sqtls_annotated, by = c("group" = "transcripts")) %>%
  select(group, ID, annotations, der_allele_hl, delta_psi) %>%
  distinct()
sqtls_annotated_subset$annotations[is.na(sqtls_annotated_subset$annotations)] <- 'intergenic'

# Calculate percentage of each
sqtl_anno_table <- 
  with(sqtls_annotated_subset 
       #%>% filter(delta_psi > .05)
       , table(annotations, der_allele_hl))

normalize <- function(x) x / sum(x)
sqtl_anno_table %<>% apply(2, normalize)

sqtl_anno_table <- sqtl_anno_table %>%
  data.frame %>%
  rownames_to_column("annotation") %>%
  pivot_longer(names_to = "exonic_DA_effect", 
               cols = c("Decreases.PSI", "Increases.PSI"), 
               names_repair = "minimal")

figS1E <- 
  ggplot(sqtl_anno_table, aes(x = exonic_DA_effect, y = value, fill = annotation)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_brewer(palette = "Set2") + 
  sqtl_manuscript_theme() + 
  labs(fill = "sVariant\nAnnotation") + 
  xlab("sVariant derived allele effect") + 
  ylab("Proportion in group")

#### Plot the Plots! ####
# fig_w = 5
# fig_h = 5
# save_plot("fig1B_sqtl_totals_summary.svg", width = fig_w, height = fig_h)
# fig1B
# dev.off()
# 
# save_plot("fig1B_sqtl_counts_per_tiss.svg", width = fig_w, height = fig_h) 
# fig1B_alt
# dev.off()
# 
# save_plot("fig1C_sqtl_symmetry_barplot.svg", width = 2, height = fig_h)
# fig1C
# dev.off()
# 
# save_plot("fig1D_high_inclusion_allele_frequency_distribution.svg",
#           width = fig_w, height = fig_h)
# fig1D
# dev.off()
# 
# save_plot("fig1E_n_high_vs_low_effect_sqtls.svg", width = 2, height = fig_h)
# fig1E
# dev.off()
# 
# save_plot("fig1F_derived_allele_effect_density.svg", width = fig_w, height = fig_h)
# fig1F
# dev.off()
# 
# save_plot("fig1G_derived_allele_effect_density_across_cutoffs.svg", width = fig_w, height = fig_h)
# fig1G
# dev.off()
# 
# save_plot("fig1S_sexon_cds_location.svg", width = 3, height = 2)
# fig_1S
# dev.off()

# Cowplot everything
library(cowplot)
f1_top <- plot_grid(NULL, fig1B, fig1C, rel_widths = c(1.5, 2, 1), nrow = 1, 
                    labels = "AUTO")
f1_bot <- plot_grid(fig1D, fig1E, fig1F, rel_widths = c(1.5, 1, 1.5), nrow = 1, 
                    labels = c("D", "E", "F"))
save_plot("fig1_full.svg", width = 12, height = 8)
plot_grid(f1_top, f1_bot, nrow = 2)
dev.off()

# Cowplot the Supplemental figures Part 1
S1_topR <- plot_grid(fig_S1B, fig_S1C, figS1D, figS1E, nrow = 2,
                     rel_widths = c(1, 1.7, 1, 1.7), labels = list("B", "C", "D", "E"))
S1_top <- plot_grid(figS1A, S1_topR, rel_widths = c(1, 1.3))
fig_S1 <- plot_grid(S1_top, figS1F, nrow = 2, rel_heights = c(1.2, 1), 
                    labels = list("A", "F"))

save_plot("fig_S1_sQTL_supplements.svg", width = 12, height = 9)
fig_S1
dev.off()

