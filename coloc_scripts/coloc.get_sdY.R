#!/nfs/sw/R/R-4.0.0/bin/Rscript
#----------------------------------
# Script to get sdY for coloc
# Author: Silva Kasela
# edited by Jonah Einson
#----------------------------------

library(tidyverse)
library(doParallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Need to specify five arguments")

pheno_file <- args[1]
cov_file <- args[2]
interact_file <- ifelse(args[3] %in% "NA", NA, args[3])
cores <- as.numeric(args[4])
out_file <- args[5]

cat("Read in phenotype file: ", pheno_file, fill = T)
pheno <- read.csv(pheno_file, sep = "\t", row.names = 4, check.names = F)
pheno <- pheno[,-c(1:5)]

cat("Covariates file: ", cov_file, fill = T)
cov <- read.csv(cov_file, sep = "\t", check.names = F, row.names = 1)
cov <- t(cov)
stopifnot(rownames(cov) == colnames(pheno))

if (!is.na(interact_file)) {
  cat("Interaction file: ", interact_file, fill = T)
  interact <- fread(here(interact_file), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
  # Add interaction variable to the covariates file
  cov <- cbind(cov, interact[match(rownames(cov), interact$ID), 2])
}

cat("Regressing the phenotypes on the covariates", fill = T)

run_regress <- function(pheno_id) {
  y <- as.numeric(pheno[pheno_id,])
  m <- lm(y ~ cov)
  sdY <- sd(m$residuals)
  names(sdY) <- pheno_id
  return(sdY)
}

if (cores > 2) {
  # Parallel computing
  num_cores <- cores - 1
  cl <- makeCluster(num_cores, outfile = "")
  registerDoParallel(cl)
  result <- foreach(i = rownames(pheno), .combine = c, .export = "cov") %dopar% run_regress(i)
  stopCluster(cl)
} else{
  result <- foreach(i = rownames(pheno), .combine = c, .export = "cov") %do% run_regress(i)
}

write.csv(result, file = out_file, col.names = F, row.names = T, sep = "\t", quote = F)

cat("Done!", fill = T)
