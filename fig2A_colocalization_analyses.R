# Plot the results from fig2A_colocalization-analyses_prelim.R to get this
# in the same place. 
rm(list = ls())
library(here)
source(here("myPackages.R"))
source(here("sqtl_manuscript_functions.R"))
library(ggplot2)

sqtls_coloc <- read_tsv(here("data/top_sQTLs_with_top_coloc_event.tsv"))

# Remove terminal exons
n_exons_per_gene <- read_tsv(here("data/gtex_v8_n_exons_per_gene.tsv"))
sqtls_coloc %<>%
  mutate(exon_number = str_split(phenotype_id, "_") %>% map(2) %>% as.integer) %>%
  left_join(n_exons_per_gene, by = c("GeneID" = "gene")) %>% 
  mutate(terminal_exon_flag = (exon_number == 1) | (exon_number == n_exons))
sqtls_coloc %<>% filter(!terminal_exon_flag)

# Plot PP.power vs. PP.coloc for all genes
fig2S_1 <- 
  ggplot(sqtls_coloc, aes(pp_power, pp_coloc, color = has_coloc)) + 
  geom_point() + 
  ylab("PP coloc (PP.H4 / (PP.H3 + PP.H4))") +
  xlab("PP power (PP.H3 + PP.H4)") +
  sqtl_manuscript_theme()

fig2S_1

## 1) How many sQTLs colocalize with something? ----
fig2A <- 
  ggplot(sqtls_coloc %>%
           mutate(has_coloc_long = 
                    ifelse(has_coloc, "At least one\ncolocalized trait", "No colocalization")), 
         aes(has_coloc_long)) + 
  geom_bar(fill = "cornflowerblue") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  ylab("ψQTL genes") + 
  ylim(0, 4500) + 
  xlab("") +
  sqtl_manuscript_theme()

fig2A

## 2) Contribution of derived alleles and exon symmetry ----
# Calculate lengths of exons

# Symmetry
calculate_exon_length <- function(x){
  k = unlist(str_split(x, "_"))
  start = as.integer(k[2])
  end = as.integer(k[3])
  end - start + 1
}

sqtls_coloc$exon_length <- sapply(qtls_coloc$Exon_coord, calculate_exon_length)
sqtls_coloc %<>% mutate(symmetric = (exon_length %% 3) == 0)

symmetry_coloc_table <- with(sqtls_coloc, table(has_coloc, symmetric))
symmetry_fisher_test <- fisher.test(symmetry_coloc_table)

symmetry_fisher_test <- 
  fisher.test.to.data.frame(symmetry_fisher_test)


# Derived allele effect direction
ancestral_alleles <- read_tsv(here("data/top_sQTLs_MAF05_w_anc_allele.tsv"))
top_sQTLs_top_coloc_anc_alleles <-
  ancestral_alleles %>%
  select(top_pid, li_allele, hi_allele, anc_allele) %>%
  right_join(top_sQTLs_top_coloc)

top_sQTLs_top_coloc_anc_alleles$der_allele_hl <-
  with(top_sQTLs_top_coloc_anc_alleles,
       # This is based on the ancestral allele, so flip it since we're interested
       # in the ancestral allele
       ifelse(anc_allele == hi_allele,
              "lower",
              ifelse(anc_allele == li_allele,
                     "higher", NA
              )
       )
  )
