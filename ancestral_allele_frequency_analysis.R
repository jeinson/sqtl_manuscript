# This script will look at the ancestral allele in top sQTLs across multiple 
# tissues, to see if the higher higher or lower included allele is more likely
# to be ancestral. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/anc_allele_analysis")
source("~/myPackages.R")
library("ggplot2")

top_sQTLs <- read_tsv("top_sQTLs_MAF05_w_anc_allele.tsv", comment = "[")

hist(top_sQTLs$hi_allele_freq, breaks = 20)
top_sQTLs$hi_allele_freq %>% range(na.rm = T)

with(top_sQTLs, sum(hi_allele == anc_allele, na.rm = T))
with(top_sQTLs, sum(li_allele == anc_allele, na.rm = T))
with(top_sQTLs, sum(is.na(anc_allele)))

# Add a column for the derived allele frequency, vs the anecestral allele freq
top_sQTLs$der_allele_freq <- 1 - top_sQTLs$anc_allele_freq

# Make a histogram of the derived allele frequency
ggplot(top_sQTLs, aes(der_allele_freq)) + 
  geom_histogram(fill = "firebrick") + 
  xlab("Derived allele frequency") + 
  ggtitle("Derived allele frequencies (gnomad) of sQTL variants")

# Add if derived is lower or higher allele
top_sQTLs$der_allele_hl <- 
  with(top_sQTLs, 
       # This is based on the ancestral allele, so flip it since we're interested
       # in the ancestral allele
       ifelse(anc_allele == hi_allele,
              "lower", 
              ifelse(anc_allele == li_allele, 
                     "higher", NA
              )
       )
  )

# Save the plot, paying attention to it looking nice
svg("../sqtl_analyses/ancestral_derived_allele_counts.svg", width = 2, height = 3)
ggplot(top_sQTLs, aes(der_allele_hl, fill = der_allele_hl)) + 
  geom_bar() + 
  xlab("Derived allele PSI effect") + 
  ylab("N Genes") + 
  ggtitle("Counts of derived higher vs.\nlower exon inclusion alleles") + 
  theme(legend.position = "none") #+
  ggpubr::add_summary() +
  gtex_v8_figure_theme()
dev.off()

der_allele_counts <- table(top_sQTLs$der_allele_hl)
binom.test(der_allele_counts[1], sum(der_allele_counts))

#### Plot higher/lower overlaps and run KS Tests ####
library(ggplot2)
svg("../sqtl_analyses/ancestral_derived_allele_frequency.svg", width = 4, height = 3)
ggplot(subset(top_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = der_allele_hl)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies for high and low inclusion sQTL haplotypes") + 
  labs(fill = "Derived allele\neffect direction") + 
  gtex_v8_figure_theme()
dev.off()

base_p <- top_sQTLs %>% 
  split(.$der_allele_hl) %>%
  ks.test(x = .[['higher']]$der_allele_freq, y = .[['lower']]$der_allele_freq) %>% 
  .[['p.value']]

# what happens if we crank up the delta-psi values
le_sQTLs <- top_sQTLs %>% filter(delta_psi > .05)

## delta-PSI > .05
svg("../sqtl_analyses/ancestral_derived_allele_frequency_dPSI_.05.svg", width = 4, height = 3)
ggplot(subset(le_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = der_allele_hl)) +
geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies for high and low inclusion sQTL haplotypes\nwith delta PSI > .05") + 
  gtex_v8_figure_theme()
dev.off()

le_sQTLs %>% 
  split(.$der_allele_hl) %>%
  ks.test(x = .[['higher']]$der_allele_freq, y = .[['lower']]$der_allele_freq) %>% 
  .[['p.value']]

## delta-PSI > 0.1
me_sQTLs <- top_sQTLs %>% filter(delta_psi > .1)
svg("../sqtl_analyses/ancestral_derived_allele_frequency_dPSI_.1.svg", width = 4, height = 3)
ggplot(subset(me_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = der_allele_hl)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies for high and low inclusion sQTL haplotypes\nwith delta PSI > .1") + 
  gtex_v8_figure_theme()
dev.off()

me_sQTLs %>%
  split(.$der_allele_hl) %>%
  ks.test(x = .[['higher']]$der_allele_freq, y = .[['lower']]$der_allele_freq) %>% 
  .[['p.value']]

## delta-PSI > 0.15
he_sQTLs <- top_sQTLs %>% filter(delta_psi > .15)
he_sQTLs <- subset(he_sQTLs, !is.na(der_allele_hl)) 

svg("../sqtl_analyses/ancestral_derived_allele_frequency_dPSI_.15.svg", width = 4, height = 3)
ggplot(he_sQTLs, aes(der_allele_freq, fill = der_allele_hl)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies for high and low inclusion sQTL haplotypes\nwith delta PSI > .15") + 
  gtex_v8_figure_theme()
dev.off()

ks_test_output <- 
  he_sQTLs %>%
  split(.$der_allele_hl) %>%
  ks.test(x = .[['higher']]$der_allele_freq, y = .[['lower']]$der_allele_freq) %>% 
  .[['p.value']]

## Plot all sQTLS and their delta psi scores
svg("../figures/sQTL_dPSI_distribution.svg", width = 5, height = 2)
ggplot(top_sQTLs, aes(delta_psi)) + 
  geom_histogram() + 
  geom_vline(xintercept = .05, col = "green") + 
  geom_vline(xintercept = .1, col = "orange") + 
  geom_vline(xintercept = .15, col = "red") + 
  xlab("ΔPSI") +
  ylab("sQTL count") +
  ggtitle("Distribition of ΔPSI across top sQTLs") + 
  gtex_v8_figure_theme()
dev.off()

# What happens if we permute everything?
permuted_qtls <- le_sQTLs %>% mutate(
  der_allele_hl = sample(der_allele_hl,nrow(.), replace = F)
  ) %>%
  filter(!is.na(der_allele_hl))

ggplot(permuted_qtls, aes(der_allele_freq, fill = der_allele_hl)) +
    geom_density(alpha = .5) + 
    ggtitle("Derived allele frequencies with permuted H/L designation")

permuted_qtls %>%
  split(.$der_allele_hl) %>%
  ks.test(x = .[['higher']]$der_allele_freq, y = .[['lower']]$der_allele_freq)

#### Follow Ups ####
# This raises a number of questions about the differences between the high inclusion
# derived allele and the low inclusion derived allele. Do the genes differ
# in length, for example? What are the pLI scores of the differences?

# Look at the length of exons between the two. 
#### Exon Length ####
exon_id_map <- read_rds("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gtex_stuff/gtex_v8_exon_id_map.rds")

top_sQTLs$exon_length <- 
  exon_id_map[top_sQTLs$top_pid] %>% 
  str_split("_") %>%
  sapply(function(x) abs(as.numeric(x[2]) - as.numeric(x[3])) + 1)

ggplot(top_sQTLs, aes(der_allele_hl, exon_length)) + 
  geom_boxplot() + 
  scale_y_log10()

# Just out of curiosity?
barplot(table(top_sQTLs$exon_length %% 3), main = "sQTL Exon lengths mod 3 (bp)")
# There doesn't seem to be a lot of bias for symmetric exons in sQTLs, though 
# it's higher than the baseline of all exons in gencode. 

# Break this up 
top_sQTLs <- 
  top_sQTLs %>%
  mutate(is_mod3 = (exon_length %% 3) == 0)

# Are lower included sQTL haplotypes more likely to not be mod3?
cgy_tbl <- table(top_sQTLs %>% select(der_allele_hl, is_mod3))
fisher.test(cgy_tbl)

#### Overall high inclusion vs. low inclusion ####
hist(top_sQTLs$mean_01_psi, main = "Median PSI of sQTL heterozygotes")
abline(v = .3, col = "red")
abline(v = .7, col = "red")
abline(v = .5, lty = 2, col = "blue")

# Add a column for average exon inclusion level
# top_sQTLs$median_psi_level <- ifelse(top_sQTLs$mean_01_psi < .3, "< 30%", 
#                               ifelse(top_sQTLs$mean_01_psi > .7, "> 70%", "mid"))
# top_sQTLs$median_psi_level <- 
#   top_sQTLs$median_psi_level %>%
#   fct_relevel("< 30%", "mid", "> 70%")

top_sQTLs$median_psi_level <- ifelse(top_sQTLs$mean_01_psi < .5, "< 50%", "> 50%")

barplot(table(top_sQTLs$median_psi_level), 
        xlab = "Median sExon inclusion", 
        ylab = "Number of sQTLs", 
        col = "firebrick")

ggplot(top_sQTLs, aes(median_psi_level, fill = der_allele_hl)) + 
  geom_bar(position = "dodge") + 
  ggtitle("Counts of exon inclusion/exclusion derived sQTL alleles\n in highly and lowly included exons")

hl_table <- table(top_sQTLs %>% select(der_allele_hl, median_psi_level)) 
fisher.test(hl_table)

# Look at the allele frequency distributon of overall higher or overall low
# inclusion exons. 
ggplot(top_sQTLs, aes(der_allele_freq, fill = median_psi_level)) + 
  geom_density(alpha = .5) 

# Look at higher effect size 
ggplot(top_sQTLs %>% filter(delta_psi > .05), aes(der_allele_freq, fill = median_psi_level)) + 
  geom_density(alpha = .5) +
  ggtitle("Delta PSI > .05")
top_sQTLs %>% filter(delta_psi > .05) %>% .$median_psi_level %>% table

# How many high inclusion exons are there as you increase the effect size cutoff?
n_sQTLs <- 
  sapply(seq(0, .49, by = .05), function(x){
    top_sQTLs %>% filter(delta_psi > x) %>% .$median_psi_level %>% table
  }
)
colnames(n_sQTLs) <- seq(0, .49, by = .05)
n_sQTLs <- 
  n_sQTLs %>% t %>% as.data.frame %>%
  rownames_to_column("cutoff") %>%
  pivot_longer(-cutoff, names_to = "median_psi", values_to = "n_sQTLs")

ggplot(n_sQTLs, aes(cutoff, n_sQTLs, fill = median_psi)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Delta PSI >")

#### derived allele distribution by median psi ####
le_sQTLs <- top_sQTLs %>% filter(delta_psi > .05)
ggplot(subset(le_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = median_psi_level)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies overall high and low\nsQTL haplotypes with delta PSI > .05")
le_sQTLs %>% 
  split(.$median_psi_level) %>%
  ks.test(x = .[['< 50%']]$der_allele_freq, y = .[['> 50%']]$der_allele_freq)

me_sQTLs <- top_sQTLs %>% filter(delta_psi > .1)
ggplot(subset(me_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = median_psi_level)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies overall high and low\nsQTL haplotypes with delta PSI > .10")
me_sQTLs %>% 
  split(.$median_psi_level) %>%
  ks.test(x = .[['< 50%']]$der_allele_freq, y = .[['> 50%']]$der_allele_freq)

he_sQTLs <- top_sQTLs %>% filter(delta_psi > .15)
ggplot(subset(me_sQTLs, !is.na(der_allele_hl)), aes(der_allele_freq, fill = median_psi_level)) +
  geom_density(alpha = .5) + 
  ggtitle("Derived allele frequencies overall high and low\nsQTL haplotypes with delta PSI > .15")
he_sQTLs %>% 
  split(.$median_psi_level) %>%
  ks.test(x = .[['< 50%']]$der_allele_freq, y = .[['> 50%']]$der_allele_freq)

#### Combine overall inclusion and effect direction ####
top_sQTLs$psi_level_and_effect_direction <- paste(
  top_sQTLs$median_psi_level,
  top_sQTLs$der_allele_hl, 
  sep = " & "
)

ggplot(top_sQTLs %>% filter(!is.na(der_allele_hl)), 
       aes(psi_level_and_effect_direction, fill = psi_level_and_effect_direction)) + 
  geom_bar()

ggplot(top_sQTLs %>% filter(!is.na(der_allele_hl)), 
       aes(der_allele_freq, fill = psi_level_and_effect_direction)) +
  geom_density() + 
  facet_wrap("psi_level_and_effect_direction") +
  theme(legend.position = "none")

# Crank up the effect size
svg("../figures/sQTL_counts_split_by_effect.svg", width = 5, height = 5)
ggplot(top_sQTLs %>% filter(!is.na(der_allele_hl) & delta_psi > .05), 
       aes(psi_level_and_effect_direction, fill = psi_level_and_effect_direction)) + 
  geom_bar() + 
  theme_classic() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  ggtitle("Counts of top sQTLs with ΔPSI > 0.05") +
  xlab("Starting exon PSI & QTL effect direction") + 
  ylab("sQTL count") +
  ylim(0, 600) +
  theme(legend.position = 'none') 
dev.off()

svg("../figures/sQTL_derived_allele_frequency_density.svg", width = 6, height = 4)
ggplot(top_sQTLs %>% filter(!is.na(der_allele_hl) & delta_psi > .05), 
       aes(der_allele_freq, fill = psi_level_and_effect_direction)) +
  geom_density() + 
  facet_wrap("psi_level_and_effect_direction") +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("Derived allele frequency") + 
  ylab("Density")
dev.off()

