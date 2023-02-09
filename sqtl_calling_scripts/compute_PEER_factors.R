########################
# This script computes PEER factors for sQTL analysis, given "violently 
# normalized" PSI values for a specified tissue
########################

args = commandArgs(trailingOnly=TRUE)
#tissue = "Adipose_Subcutaneous" # This is the name of the directory containing info
tissue = args[1]
gtex_tissue = args[2]

if(is.na(gtex_tissue)) {
  gtex_tissue <- tissue
}

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno")
#source("~/myPackages.R")
library("readr")
library("purrr")
library("stringr")

library(peer, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4/")
num_infer = 15 # go with this from the GTEx sqtl analysis

eqtl_covariates <-
  read_tsv(sprintf("/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/%s.v8.covariates.txt", gtex_tissue))

genotype_PCs <- eqtl_covariates[1:5,]
technical_factors <- eqtl_covariates[(nrow(eqtl_covariates)-2):nrow(eqtl_covariates),]

# Read in the "violently normalized" expression data
expr <- read_tsv(sprintf("bed/violently_normalized_%s_psi_v8_protocol.bed.gz", tissue))
expr <- expr[,-(1:4)]
expr <- t(expr) # This puts it in the format for PEER, with the samples in rows
expr <- as.matrix(expr)

# Fill in NA rows with the mean of each gene's expression
# This seems like it's a problem unique to splicing data, but PEER will 
# crash if you don't do this. 
expr <- 
  apply(expr, 2, function(x){
    x[is.na(x)] <- mean(x, na.rm = T)
    x
  })

################# Run the pipeline with the full dataset #######################
covar <- as.matrix(genotype_PCs[,-1])
samp_names <- rownames(expr) %>% str_split("-") %>% map(~.x[1:2]) %>% map_chr(~paste(.x, collapse="-"))
covar <- covar[,samp_names]

tech_covar <- as.matrix(technical_factors[,-1])
tech_covar <- tech_covar[,samp_names]

covar <- rbind(covar, tech_covar)
# From the Documentation: 
# The C observed covariates are assumed to be in a NxC matrix:
covar <- t(covar)

# This bit is copy-pasted from every script that's used to run PEER....

names.use <- rownames(expr)

message("Running PEER")
model <- peer::PEER()
x = peer::PEER_setNk(self = model, Nk = num_infer)

# Set the mean
x = peer::PEER_setPhenoMean(self = model, matrix = expr)

# Add genetic ancestry covariates
x = peer::PEER_setCovariates(self = model, matrix = covar)

# Just for testing that the dimensions look right
# x = peer::PEER_setNmax_iterations(self = model, Nmax_iterations = 1)

x = peer::PEER_update(self = model)

############## Get and assemble the results from PEER ##########################
message("Assembling results")

peer.covar <- as.data.frame(x = t(x = peer::PEER_getX(self = model)))
peer.alpha <- as.data.frame(x = peer::PEER_getAlpha(self = model))
peer.resid <- as.data.frame(x = peer::PEER_getResiduals(self = model))
covar.names <- c( # Get the names correct
  paste0("PC", 1:5),
  technical_factors$ID,
  paste0('InferredCov', 1:num_infer)
)

# Assemble back together
peer.covar <- cbind(covar.names, peer.covar)
colnames(peer.covar) <- c('ID', names.use)

# Rearrange so it's in the same order
peer.covar <- peer.covar[c(1:5, 9:23, 6:8),]

# Save to disk
write_tsv(peer.covar, sprintf("covariates/%s_covariates.tsv", tissue))

# To submit this as a job to the cluster
# echo "module load R/3.4.1; Rscript compute_PEER_factors.R" | qsub -V -cwd -N PEER_psi -m ae -M jeinson@nygenome.org -o logs -e logs