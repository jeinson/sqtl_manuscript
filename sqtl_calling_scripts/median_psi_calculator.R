############################
# Calculate Median PSI
# 
# Jonah Einson
# 9/24/21
############################

source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno")

all_tissues <- read_lines("completed_tissues.txt")

exon_id_map <- readRDS("../../../data/gtex_stuff/gtex_v8_exon_id_map.rds")
id_exon_map <- names(exon_id_map)
names(id_exon_map) <- exon_id_map

for(tissue in all_tissues){
  message(tissue)
  # if(interactive()) {
  #   tissue = "Whole_Blood"
  # } else {
  #   args <- commandArgs(trailingOnly = TRUE)
  #   tissue = args[1]
  # }
  # 
  
  fp <- paste0("bed/", tissue, "_psi_v8_protocol.bed.gz")
  dat <- read.csv(fp, sep = "\t", row.names = 4)
  dat <- dat[,-(1:3)]
  
  rowMedians <- apply(dat, 1, median, na.rm = T)
  meanPSI <- apply(dat, 1, mean, na.rm = T)
  rowSDs <- apply(dat, 1, sd, na.rm = T)
  
  out <- enframe(rowMedians, "pid", "median_psi")
  out$mean_psi <- meanPSI
  out$sd_psi <- rowSDs
  
  # Add the exon ID from the official annotation
  out$exon_id <- id_exon_map[out$pid]
  
  out_fp <- paste0("median_psi/", tissue, "_median_psi.tsv")
  write_tsv(out, out_fp)
}

all_med_psi <- lapply(all_tissues, function(x){
  fp <- paste0("median_psi/", x, "_median_psi.tsv")
  read_tsv(fp)
    })
names(all_med_psi) <- all_tissues
all_med_psi <- bind_rows(all_med_psi, .id = "tissue")
write_tsv(all_med_psi, "median_psi/median_psi_all_tissues.tsv")
