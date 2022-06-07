# Supplemental figure for the ΨQTL manuscript, looking at some properties
# of symmetric / non-symmetric exons
# In this script, I dig a bit deeper into symmetric and non-symmetric sQTLs.  
# What is going on with these? 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/sqtl_analyses")
suppressPackageStartupMessages(source("~/myPackages.R"))
source("../../sqtl_manuscript/sqtl_manuscript_functions.R")

sqtls <- read_tsv("../multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs_exons_labeled_by_coord.tsv")

# Relabel exons by their ID (instead of their coordinates)
id_exon_map <- read_rds("../../../../data/gtex_stuff/gtex_v8_id_exon_map.rds")
sqtls$exon_id <- id_exon_map[sqtls$pid]

# Limit to sQTLs where only one exon is significant. This will probably simplify
# things moving forward. 
genes_with_one_sig_exon <- 
  sqtls %>%
  select(group, n_sig) %>%
  filter(n_sig == 1) %>%
  distinct %>% 
  .$group

sqtls <- filter(sqtls, group %in% genes_with_one_sig_exon)

# Calculate the lengths of exons
calculate_exon_length <- function(x){
  k = unlist(str_split(x, "_"))
  start = as.integer(k[2])
  end = as.integer(k[3])
  end - start + 1
}

sqtls$exon_length <- sapply(
  sqtls$pid, calculate_exon_length
)

# Get rid of sQTLs where the sExon is a terminal exon.
n_exons_per_gene <- read_tsv("../../../../data/gtex_stuff/gtex_v8_n_exons_per_gene.tsv")
sqtls %<>%
  mutate(exon_number = str_split(exon_id, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("group" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons.x))

   # Real quick, how is the significant exon a terminal exon
   barplot(table(sqtls$terminal_exon_flag), main = "Significant exon is the terminal exon")

sqtls %<>% filter(!terminal_exon_flag)

#### Plot 1 - Lenghts of all sExons ####
hist(log10(sqtls$exon_length), col = "cornflowerblue", 
     main = "Distribution of sExon lengths", 
     xlab = "log10(Exon Length) (bp)")
lim <- quantile(sqtls$exon_length)[4] + (1.5 * IQR(sqtls$exon_length))
abline(v = log10(lim), col = "red", lty = 2)

# Dump abnormally long exons
#     Actually nevermind, don't do
# sqtls <- sqtls[!flag_outliers(sqtls$exon_length),]
# barplot(table(sqtls$exon_length %% 3), 
#         main = "sExon Lengths mod 3")
# 
# hist(sqtls$exon_length)
# 


# 1) Is being symmetric affected by exon length? ----
sqtls$is_symmetric <- (sqtls$exon_length %% 3) == 0
barplot(table(sqtls$is_symmetric))

library(ggplot2)
ggplot(sqtls, aes(exon_length, fill = is_symmetric)) + 
  geom_density(alpha = .5)

# 2) How about overall PSI? ----
# Are non symmetric less spliced in in general?
med_psi <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/median_psi/median_psi_all_tissues.tsv")
med_psi <- dplyr::rename(med_psi, "tiss" = "tissue")
sqtls <- left_join(sqtls, med_psi)

ggplot(sqtls, aes(mean_psi, fill = is_symmetric)) + 
  geom_density(alpha = .5)

ggplot(sqtls, aes(sd_psi, fill = is_symmetric)) + 
  geom_density(alpha = .5)

# 3) (non)-symmetric status all exons? (Not just sExons?) ----
# Let's just read in all protein coding genes from the GTEx gtf, and 
# perform the same calculations. 
library(rtracklayer)
gtf <- rtracklayer::import("../../data/gencode.v26.GRCh38.GTExV8.genes.gtf")
#gtf <- rtracklayer::import("../../data/gencode.v26.annotation.gtf.gz")
gtf_df <- data.frame(gtf)
all_exons <- filter(gtf_df, transcript_type == "protein_coding" & type == "exon") %>%
  select(gene_id, transcript_id, exon_id, seqnames, start, end, strand, width)

all_exons$exon_number <- 
  str_split(all_exons$exon_id, " ") %>%
  map_chr(3) %>%
  as.integer()

all_exons <- 
  all_exons %>% 
  split(.$gene_id) %>%
  lapply(function(x) {
    # Actually, don't flip it around. Francois already did this!! :-) 
    x$normalized_cds_loc <- (cumsum(x$width) + (x$width/2)) / sum(x$width)
    x$exon_pct <- x$exon_number / nrow(x)
    
    x
  }) %>%
  bind_rows

# Add a BOOL column for being a terminal exon
all_exons <- 
  all_exons %>%
  group_by(gene_id) %>%
  mutate(terminal = (exon_number == 1 | exon_number == n())) %>%
  ungroup

all_exons$is_symmetric <- (all_exons$width %% 3 == 0)

# Plot these side by side. Are non-terminal exons 
# more likely to be symmetric?
all_exons_sym_ct <- table(all_exons$width %% 3)
nt_exons_sym_ct <- table(all_exons[!all_exons$terminal, ]$width %% 3)

barplot(rbind(
  prop.table(all_exons_sym_ct),
  prop.table(nt_exons_sym_ct)
), beside = T, col = c("salmon", "lightgreen"), 
xlab = "Exon length (bp) mod 3"
)
legend("topright", 
       legend = c("All exons", "All non-terminal-exons"), 
       col = c("salmon", "lightgreen"), 
       pch = 15
       )

# Now, can I replicate the result from the Magen and Ast Paper?

# First, get a dictionary of cross tissue median PSI, to find constitutive vs
# alternative exons
cross_tiss_med_psi <- med_psi %>% split(.$pid) %>% map_dbl(~ median(.$median_psi))
cross_tiss_mean_psi <- med_psi %>% split(.$pid) %>% map_dbl(~ mean(.$mean_psi))

all_exons$exon_coord <- with(all_exons, paste(seqnames, start, end, strand, sep = "_"))
all_exons$cross_tissue_med_psi <- cross_tiss_med_psi[all_exons$exon_coord]

all_exons$constitutive <- with(all_exons, cross_tissue_med_psi == 1)
all_exons <- all_exons %>% filter(cross_tiss_med_psi > 0 & !is.na(cross_tissue_med_psi))

# Now try to make the plots that look really awful in the original paper
all_exons$normalized_cds_loc_bin <- 
  cut(all_exons$normalized_cds_loc, breaks = seq(0, 1, by = .1))
all_exons$exon_pct_bin <- 
  cut(all_exons$exon_pct, breaks = seq(0, 1, by = .1))

par(mfrow = c(2,2))
all_exons %>%
  filter(is_symmetric & !constitutive) %$%
  table(exon_pct_bin) %>%
  prop.table %>%
  barplot(main = "Symmetrical alternative exons", 
          xlab = "Normalized location in CDS", 
          ylab = "Percentage of exons", 
          col = "grey3")

all_exons %>%
  filter(is_symmetric & constitutive) %$%
  table(exon_pct_bin) %>%
  prop.table %>%
  barplot(main = "Symmetrical constituve exons", 
          xlab = "Normalized location in CDS", 
          ylab = "Percentage of exons", 
          col = "grey3")

all_exons %>%
  filter(!is_symmetric & !constitutive) %$%
  table(exon_pct_bin) %>%
  prop.table %>%
  barplot(main = "Non-Symmetrical alternative exons", 
          xlab = "Normalized location in CDS", 
          ylab = "Percentage of exons")

all_exons %>%
  filter(!is_symmetric & constitutive) %$%
  table(exon_pct_bin) %>%
  prop.table %>%
  barplot(main = "Non-Symmetrical constitutive exons", 
          xlab = "Normalized location in CDS", 
          ylab = "Percentage of exons")

# Are non-symmetric exons more likely to pop up in the later parts of the 
# transcript? I think I could just plot this. 
ggplot(all_exons, aes(normalized_cds_loc, color = is_symmetric)) + 
  geom_density()

all_exons$cross_tiss_median_psi <- all_exons$exon_id

# 4) Statistical comparison of sExons to all exons ----
c.table <- rbind(table(all_exons$width %% 3), table(sqtls$exon_length %% 3))
rownames(c.table) <- c("all_exons", "s_exons")

# Make sure to only use non-terminal exons!
all_exons <- all_exons %>% filter(!terminal)

is.symmetric.c.table <- rbind(
  table(all_exons$is_symmetric),
  table(sqtls$is_symmetric)
)
rownames(is.symmetric.c.table) <- c("all_exons", "s_exons")
fisher.test(is.symmetric.c.table)

par(mfrow = c(1,1))
p.table <- t(apply(c.table, 1, function(x) x / sum(x)))

# save_plot("fig1_sqtl_symmetry_barplot.svg", width = 5, height = 5)
# barplot(p.table, beside = T, 
#         col = c("aquamarine", "coral"), 
#         main = expression('Exon length ' *italic(mod)*' 3'), 
#         ylab = "Percent of all exons", space = c(.2, 0, .2, 0, .2, 0))
# legend("topright", c("All exons", "sExons"), col = c("aquamarine", "coral"), 
#        pch = 15)
# dev.off()

# Remake this in ggplot
save_plot("fig1_sqtl_symmetry_barplot.svg", width = 5, height = 5)
as.data.frame(p.table) %>% 
  rownames_to_column("group") %>% 
  pivot_longer(cols = c('0', '1', '2'), names_to = "Exon_length", values_to = "Percent") %>%
  ggplot(aes(Exon_length, Percent, fill = group)) + 
  scale_fill_manual(values = c("aquamarine", "coral")) + 
  geom_bar(stat = "identity", position = "dodge") + 
  sqtl_manuscript_theme() + 
  theme(legend.position = c(.8, .9)) + 
  xlab("Exon Length mod 3") + ylab("Percent of All Exons")
dev.off()

# 5) Comparison to full gencode v26 set ----
gencode_v26 <- rtracklayer::import("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/v8_psi/psi_annotation/gencode.v26.annotation.gtf")
all_exon_widths <- width(gencode_v26[gencode_v26$type == "exon" & 
                                       gencode_v26$gene_type == "protein_coding"])

v26_exons_sym_ct <- table(all_exon_widths %% 3)
barplot(rbind(
  prop.table(all_exons_sym_ct),
  prop.table(nt_exons_sym_ct),
  prop.table(v26_exons_sym_ct)
), beside = T, col = c("salmon", "lightgreen", "violet")
)
