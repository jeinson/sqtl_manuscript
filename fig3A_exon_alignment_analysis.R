# Read in the results from sExon alignment to MANE, and plot some summary
# stats and basic results. 

library(here)
source(here("myPackages.R"))
source(here("sqtl_manuscript_functions.R"))

# Read in data
sqtls <- read_tsv(here("data/top_sQTLs_MAF05.tsv"))
alignments <- read_csv(here("data/top_sQTLs_top_coloc_mapped.csv"))
alignments <- distinct(alignments)

# Remove terminal exons
n_exons_per_gene <- read_tsv(here("data/gtex_v8_n_exons_per_gene.tsv"))
sqtls %<>%
  mutate(exon_number = str_split(top_pid, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("group" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons))
sqtls %<>% filter(!terminal_exon_flag)

# Align the top sQTLs with the alignments
alignments <- 
  alignments %>%
  mutate(group = remove_trailing_digit(GENE.ID)) %>%
  inner_join(
    sqtls
  )

# Is there an alignments?
alignments$exon_in_MANE <- 
  !(alignments$MAPPED.TO.OTHER | is.na(alignments$MAPPED.TO.OTHER)) | 
  !is.na(alignments$MAPPED.TO.OTHER) # yikes!

# Add median PSI
median_psi <- read_tsv(here("data/top_sQTLs_median_psi.tsv"))
alignments <- left_join(alignments, median_psi, 
                        by = c('tiss', 'top_pid' = 'exon_id'))

# How does this match up?
library(ggplot2)
ggplot(alignments, aes(median_psi, fill = exon_in_MANE)) + 
  geom_density(alpha = .5) + 
  xlab("Median PSI") +
  sqtl_manuscript_theme()

