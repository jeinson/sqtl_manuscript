# This script compiles a sample sheet for running coloc on sQTLs by tissue
# and GWAS hits. 

# This will have 3 nested for loops. Loop through every tissue, every GWAS
# trait, and every chromosome. 
library(tidyverse)
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno")


tissues=read_lines("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/completed_tissues.txt")
gwas_fps=list.files("/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/gwas", 
                    full.names = T, pattern = ".txt.gz")
n_samples=read_delim("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/files/samples_with_n.txt",
                     delim=" ", col_names = F)

# Grab the number of samples in each GWAS
gwas_n_table <- tibble()
for(gwas_fp in gwas_fps){
  line = system(paste("zcat", gwas_fp, "| head -2 | tail -1"), intern = T)
  line = strsplit(line, "\t", )[[1]]
  gwas_n = as.numeric(line[4])
  gwas_n_table <- rbind(gwas_n_table,
                        c(gwas_fp, gwas_n))
}

row_n <- length(tissues) * nrow(gwas_n_table) * 22
out_list <- list()
for(tiss in tissues){
  out <- tibble()
  message(paste("Processing", tiss))
  # Grab covariates
  cov_fp=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/covariates/%s_covariates.tsv",
                 tiss)
  
  qtl_n=with(n_samples, X2[grepl(tiss, X1)])
  
  # Grab the sdY file
  sdY_fp=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sdY/%s_normalized_PSI_sdY.tsv",
                 tiss)
  
  for(row in 1:nrow(gwas_n_table)){
    # Process gwas info 
    
    gwas_name=sapply(strsplit(gwas_n_table[row, 1], "/"), '[[', 10) %>% str_remove(".txt.gz")
    gwas_fp = gwas_n_table[row, 1]
    gwas_n = gwas_n_table[row, 2]
    
    for(chr in seq(1, 22)){
      qtl_fp=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/qtltools_results/%s_nominal/%s.qtltools.nominal.MAF05.chr%i.06.02.21.txt.gz",
                     tiss, tiss, chr)
      
      out <- rbind(out, 
                   c(tiss, chr, qtl_fp, qtl_n, cov_fp, sdY_fp, gwas_fp, gwas_name, gwas_n))
      
    }
  }
  colnames(out) <- str_c("X", 1:9)
  out_list[[tiss]] <- out
}

out <- bind_rows(out_list)

names(out) <- str_c("X", 1:9)
out$X8 <- str_remove(out$X8, "imputed_")
write_tsv(out, "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sample_sheets/sQTL_gwas_coloc_sample_sheet.tsv", 
          col_names = F)

out_height = out[grep("height", out[,7]),]
write_tsv(out_height, "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sample_sheets/sQTL_gwas_coloc_sample_sheet_height.tsv")


# Grab the samples which are actually used in GTEx v8
library(readxl)
s11 <- read_excel("aaz1776_Suppl.Excel_seq1_July20.xlsx", sheet = 2, skip = 1)

out_gtex <- out[out$X8 %in% s11$Tag,]

# Upate the n-sample info directly from the sample sheet
n_key <- deframe(select(s11, Tag, Sample_Size))
out_gtex$X9 <- n_key[out_gtex$X8]

# Add a column for number of cases in cc studies
case_key <- deframe(select(s11, Tag, Cases))
out_gtex$X10 <- case_key[out_gtex$X8]

# Add a column for quant/cc 
type_key <- deframe(select(s11, Tag, Binary))
out_gtex$X11 <- ifelse(type_key[out_gtex$X8] == 0, "quant", "cc")

write_tsv(out_gtex, "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sample_sheets/sQTL_gwas_coloc_gtex_final_list.tsv")

