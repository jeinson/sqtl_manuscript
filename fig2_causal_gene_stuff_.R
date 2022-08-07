# Figure 2 of the ψQTL paper, which looks into colocalization and shared
# causal variants with eQTLs. 
# 
# Jonah Einson
# 5/2/22

# Panel 2A looks at shared causal variants between eQTLs and sQTLs
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie")

rm(list = ls())
source("~/myPackages.R")
source("../sqtl_manuscript/sqtl_manuscript_functions.R")
library(ggplot2)


#### Figure 2B: Causal gene overlap between sQTLs and eQTLs
tissues <- read_lines("../sQTL_v8_anno/completed_tissues.txt")
tissues <- tissues[-which(tissues == "Colon_Transverse")]

ct_table_list <- list()
mean_exon_length_mann_test_list <- list()
for(tiss in tissues){
  message(tiss)
  overlaps <- read_tsv(sprintf("credible_set_overlaps/%s_cs_overlaps.tsv", tiss))
  
  all_css_raw <- read_tsv(
    sprintf("sqtl_v8_anno_collapsed_credible_sets/GTEx_psi_%s.collapsed.txt.gz", 
            tiss))
  
  # Do a bit of reformatting so the exons and variants in the credible set 
  # are in their own nested lists
  all_css <- 
    all_css_raw %>%
    select(-ref, -alt, -chr, -pos) %>%
    mutate(exons_in_cs = strsplit(exons_in_cs, split = ",")) %>%
    unnest(c(exons_in_cs)) %>%
    chop(variant_id) %>%
    chop(exons_in_cs)
  
  # Calculate lengths of all exons in each CS
  exon_id_map <- readRDS("../../../data/gtex_stuff/gtex_v8_exon_id_map.rds")
  
  # Takes a vector of exons and returns their lengths
  get_exon_length <- function(x){
    x <- exon_id_map[x]
    y <- str_split(x, "_")
    as.numeric(map_chr(y, 3)) - as.numeric(map_chr(y, 2)) + 1
  }
  
  all_css <- 
    all_css %>%
    mutate(exon_lengths = map(exons_in_cs, get_exon_length))
  
  # Now add info about exon symmetry
  all_css <- 
    all_css %>%
    mutate(indv_exon_symmetry = map(exon_lengths, ~ .x %% 3 == 0))
  all_css$cs_contains_nonsymmetric_exon <- 
    as.logical(map_int(all_css$indv_exon_symmetry, ~ sum(!.x)))
  
  all_css <- 
    all_css %>%
    mutate(total_exon_symmetry = map_lgl(exon_lengths, ~ sum(.x) %% 3 == 0))
  
  # Merge the table with eQTLs. Make a new table with eqtl overlaps
  all_css_ol <- 
    all_css %>%
    left_join(overlaps %>% select(-exons_in_cs), 
              by = c("collapsed_cs_id" = "sCS_id"))
  
  all_css_ol$cs_has_ol <- !is.na(all_css_ol$eCS_id)
  
  # Fisher test if credible sets where the total length of all exons is 
  # divisible by 3. 
  total_symmetry_ctable <- with(all_css_ol, table(cs_has_ol, total_exon_symmetry))
  ct_table_list[[tiss]] <- total_symmetry_ctable
  fisher.test(total_symmetry_ctable)
  
  # Fisher test if credible sets that contain at least one nonsymmetric exon 
  # are more likely to share a causal variant with an eQTL
  at_least_one_nonsymmetric_ctable <- 
    with(all_css_ol, table(cs_has_ol, cs_contains_nonsymmetric_exon))
  at_least_one_nonsymmetric_ctable
  fisher.test(at_least_one_nonsymmetric_ctable)
  
  # Now check if overlapping sets are generally longer. i.e. there's a large
  # splicing effect that's appearing as an eQTL. 
  all_css_ol$total_exon_length <- map_dbl(all_css_ol$exon_lengths, sum)
  all_css_ol$total_exon_length_norm <- inv_norm(all_css_ol$total_exon_length)
  ggplot(all_css_ol, aes(total_exon_length, color = cs_has_ol)) + 
    geom_density()
  # Do a t-test for difference in means
  mean_exon_length_mann_test_list[[tiss]] <- 
    wilcox.test(formula = total_exon_length ~ cs_has_ol, data = all_css_ol, 
                conf.int = T, conf.level = .95)
  
  # Try doing this as a logistic regression
  logistic_regression <- 
    glm(as.numeric(cs_has_ol) ~ total_exon_length, 
        family = binomial(link = "logit"), 
        data = all_css_ol)
}

### Now plot the overall trend across tissues ####
ft_test_list <- lapply(ct_table_list, fisher.test)
n_es_genes_per_tiss <- sapply(ct_table_list, function(x) sum(x[2,]))

extract_data_from_htest <- function(htest){
  tibble(
    p.value = htest$p.value,
    estimate = htest$estimate, 
    conf.int.low = htest$conf.int[1], 
    conf.int.hi = htest$conf.int[2]
  )
}

fisher.test.results <- 
  bind_rows(lapply(ft_test_list, extract_data_from_htest), .id = "tissue")

fisher.test.results$n_es_genes_per_tiss <- n_es_genes_per_tiss

fisher.test.results <- arrange(fisher.test.results, estimate)
fisher.test.results$tissue <- fct_inorder(fisher.test.results$tissue)

fisher.test.results <- 
  fisher.test.results %>%
  mutate(log.estimate = log10(estimate), 
         log.conf.int.low = log10(conf.int.low),
         log.conf.int.hi = log10(conf.int.hi))

library(ggplot2)
fisher.test.result.figure <- 
  ggplot(fisher.test.results, 
         aes(tissue, log.estimate, 
             ymin = log.conf.int.low, ymax = log.conf.int.hi)) + 
  geom_point() + 
  geom_pointrange(color = ifelse(fisher.test.results$p.value > .05, "deepskyblue4", "deepskyblue")) +
  geom_hline(yintercept = 0, lty = 2, color = "red") + 
  ylab("Log(Odds Ratio) of fully\nnon-symmetric exon sets among eQTLs") +
  coord_flip() + 
  sqtl_manuscript_theme()

n_per_tiss_barplot <- 
  ggplot(fisher.test.results, aes(tissue, n_es_genes_per_tiss)) + 
  geom_bar(stat = "identity", fill = "deepskyblue4") + 
  ylab("N [e/s]Genes\nper tissue") + 
  gtex_v8_figure_theme() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  sqtl_manuscript_theme()

save_plot("fig2A_shared_causal_gene_symmetry_enrichment.svg", width = 8, height = 4)
plot_grid(fisher.test.result.figure, n_per_tiss_barplot, 
          ncol = 2, rel_widths = c(4, 1), align = "h")
dev.off()

### Plot the results of a T-test across tissues ####
extract_data_from_mann_test <- function(mann_test){
  tibble(
    p.value = mann_test$p.value,
    estimate = mann_test$estimate, 
    conf.int.low = mann_test$conf.int[1], 
    conf.int.hi = mann_test$conf.int[2]
  )
}

mann.test.results <- 
  bind_rows(lapply(mean_exon_length_mann_test_list, extract_data_from_mann_test), .id = "tissue")

mann.test.results$n_es_genes_per_tiss <- n_es_genes_per_tiss

# Put in the same order as the fisher.test for symmetry
# mann.test.results$tissue <- fct_inorder(fisher.test.results$tissue)
mann.test.results <- 
  mann.test.results[match(fisher.test.results$tissue, mann.test.results$tissue),]


library(ggplot2)
mann.test.results.figure <- 
  ggplot(mann.test.results, 
         aes(tissue, estimate, 
             ymin = conf.int.low, ymax = conf.int.hi)) + 
  geom_point() + 
  geom_pointrange(color = ifelse(fisher.test.results$p.value > .05, "deepskyblue4", "deepskyblue")) +
  geom_hline(yintercept = 0, lty = 2, color = "red") + 
  ylab("Difference in total spliced exon length between sQTLs that share a causal variant with an eQTL\nMann-Whitey U-Test Statistic +/- 95% Confidence Interval") +
  coord_flip() + 
  sqtl_manuscript_theme()

mann.test.results.figure
