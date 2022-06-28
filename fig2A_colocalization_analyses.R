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

# Add has_coloc flag
sqtls_coloc$has_coloc <- sqtls_coloc$dist < .25

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
# Symmetry ----
calculate_exon_length <- function(x){
  k = unlist(str_split(x, "_"))
  start = as.integer(k[2])
  end = as.integer(k[3])
  end - start + 1
}

sqtls_coloc$exon_length <- sapply(sqtls_coloc$Exon_coord, calculate_exon_length)
sqtls_coloc %<>% mutate(symmetric = (exon_length %% 3) == 0)

symmetry_coloc_table <- with(sqtls_coloc, table(symmetric, has_coloc))
symmetry_coloc_table
symmetry_fisher_test <- fisher.test(symmetry_coloc_table)
symmetry_fisher_test

symmetry_fisher_test <- 
  fisher.test.to.data.frame(symmetry_fisher_test)


# Derived allele effect direction ----
# Read in the file compiled from the previous script. 
ancestral_alleles <- read_tsv(here("data/top_sQTLs_with_top_coloc_event_aa_and_beta.tsv"))

sqtls_coloc_aa_beta <- 
  sqtls_coloc %>%
  left_join(
    ancestral_alleles
  )

allele_caller <- function(ref, alt, beta, hi = T){
  if(!is.na(beta)){
    if(beta > 0){
      li_allele = ref; hi_allele = alt
    } else if(beta < 0){
      hi_allele = ref; li_allele = alt
    }
    ifelse(hi, hi_allele, li_allele)
  }
  else {
    NA
  }
}

sqtls_coloc_aa_beta$hi_allele <- 
  with(sqtls_coloc_aa_beta, mapply(allele_caller, ref, alt, beta, hi = T))
sqtls_coloc_aa_beta$li_allele <- 
  with(sqtls_coloc_aa_beta, mapply(allele_caller, ref, alt, beta, hi = F))

sqtls_coloc_aa_beta$der_allele_hl <-
  with(sqtls_coloc_aa_beta,
       # This is based on the ancestral allele, so flip it since we're interested
       # in the ancestral allele
       ifelse(anc_allele == hi_allele,
              "lower",
              ifelse(anc_allele == li_allele,
                     "higher", NA
              )
       )
  )

da_coloc_table <- with(sqtls_coloc_aa_beta, table(der_allele_hl, has_coloc))
da_coloc_table
da_fisher_test <- fisher.test(da_coloc_table)
da_fisher_test

da_fisher_test <- 
  fisher.test.to.data.frame(da_fisher_test)

## 3) Distance from Exon to sSNP ----
top_sqtls <- read_tsv(here("data/top_sQTLs_MAF05.tsv"))
sqtls_coloc %<>%
  left_join(
    top_sqtls %>% select(top_pid, dist) %>%
      rename("snp_dist" = "dist", 
             "phenotype_id" = "top_pid")
  )

fig2B_coloc_dist <- 
  ggplot(sqtls_coloc, aes(abs(snp_dist), fill = has_coloc)) + 
  geom_density(alpha = .5) + 
  scale_x_log10() + 
  sqtl_manuscript_theme()
fig2B_coloc_dist
## 4) Delta PSI of top exon association ----
top_sqtls_anc_alleles <- read_tsv(here("data/top_sQTLs_MAF05_w_anc_allele.tsv"))
sqtls_coloc %<>%
  left_join(
    top_sqtls_anc_alleles %>%
      select(top_pid, delta_psi), 
    by = c("phenotype_id" = "top_pid")
  )

ggplot(sqtls_coloc, aes(delta_psi, fill = has_coloc)) + 
  geom_density(alpha = .5) + 
  scale_x_log10()
