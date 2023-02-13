#!/nfs/sw/R/R-4.0.0/bin/Rscript

# This script is designed to run colocalization analysis of GTEx psi-sQTLs 
# against GWAS hits, using p-values (as opposed to beta values). It can take a 
# GWAS and sQTL file and input. 

library(argparse)

# Master Parameters
# This is really helpful when test running this :-) thanks me from a year ago!
TEST_QTL_TRAIT="Adipose_Subcutaneous"
TEST_CHR=22
TEST_GWAS_TRAIT="imputed_UKB_50_Standing_height"

# QTL Info
TEST_QTL_FILE=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/qtltools_results/%s_nominal/%s.qtltools.nominal.MAF05.chr%i.06.02.21.txt.gz", TEST_QTL_TRAIT, TEST_QTL_TRAIT, TEST_CHR)
TEST_SDY_FILE=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sdY/%s_normalized_PSI_sdY.tsv", TEST_QTL_TRAIT)
TEST_COVARIATES_FILE=sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/covariates/%s_covariates.tsv", TEST_QTL_TRAIT)

samples <- read.csv("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/files/samples_with_n.txt", header=F, sep = " ")
TEST_QTL_N=samples$V2[grepl(TEST_QTL_TRAIT, samples$V1)]

# GWAS Info
TEST_GWAS_FILE=sprintf("/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/gwas/%s.txt.gz", TEST_GWAS_TRAIT)
TEST_GWAS_N=13550
CORES=4

GENOFILE="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"

parser <- ArgumentParser()

# Use Argparser to help with command line input
parser$add_argument("--qtl_trait", default=TEST_QTL_TRAIT, 
					help="Name of the QTL trait being tested. (Probably a tissue name)")
parser$add_argument("--gwas_trait", default=TEST_GWAS_TRAIT, 
					help="Name of the GWAS trait being tested")
parser$add_argument("--qtl_fp", default=TEST_QTL_FILE,
                    help="Path to file containing the results from qtltools")
parser$add_argument("--sdY", default=TEST_SDY_FILE, 
					help="Path to file containing sdY estimates for QTL data")
parser$add_argument("--covariates", default=TEST_COVARIATES_FILE, 
					help="Path to file containing covariates used for qtl analysis")
parser$add_argument("--gwas_fp", default=TEST_GWAS_FILE,
                    help="Path to file containing GWAS hits")
parser$add_argument("--gwas_type", default="quant", 
					help="Specifies if the GWAS trait is quantitative (quant) or case-control (cc). Defaults to 'quant'")
parser$add_argument("--chr", default=TEST_CHR, type="integer",
                    help="Chromosome number (no 'chr')")
parser$add_argument("--gwas_n", default=TEST_GWAS_N, type = "integer",
                    help="The number of samples in the GWAS")
parser$add_argument("--qtl_n", default=TEST_QTL_N, type="integer",
                    help="The number of samples in the QTL analysis")
parser$add_argument("--min_gwas_p", default = 1e-5, type="double", 
					help="The minimum GWAS p-value cutoff for performing colocalization. Defaults to 1e-5")
parser$add_argument("--genofile", default=GENOFILE, 
					help="Path to the sample genotype phased VCF. Defaults to the GTEx v8 path")
parser$add_argument("--priors", default="1e-4,1e-4,5e-6",
					help="A comma separated list of priors used for p1, p2, and p12 respectively. Defaults to '1e-4,1e-4,5e-6'")
parser$add_argument("--cores", default=CORES, type="integer", 
                    help="The number of cores to use while multi-threading")
args <- parser$parse_args()

# Load the rest of the packages
library(coloc)
library(tidyverse)
library(foreach)
library(doParallel)

# just to test

##### Read in sQTL Data ####
message("Reading QTL File:")
message(args$qtl_fp)

qtl_header <- read_lines("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL/files/qtltools_nominal_header_for_coloc.txt")
qtls <- read_delim(args$qtl_fp, delim = " ", col_names = qtl_header)

# Add the variant MAF data by loading in files containing variant info
MAF_file_path <- sprintf(
	"/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/data/common_variants/common_variants_chr%i.txt",
	args$chr
)
maf <- read_tsv(MAF_file_path, col_names=c("id", "qtl_maf"))
# Fix the MAf so it's actually a MAF
maf$qtl_maf <- sapply(maf$qtl_maf, function(x){min(x, 1-x)})

qtls <- left_join(qtls, maf)

# Calculate SE(beta), from p-value and beta-hat
# This column isn't provided from fastQTL, but can easily be calculated
calculate_se <- function(pval, bhat){
  bhat = abs(bhat)
  bhat / qnorm(pval/2, lower.tail = FALSE)
}

qtls$qtl_se <- calculate_se(qtls$qtl_pval, qtls$qtl_beta)

# Read in the sdY file, which gives the standard deviation of the phenotypes, when 
# accounting or covariates
sdY <- read_csv(args$sdY, 	col_names = c("phenotype_id", "qtl_sdY"), skip = 1)
qtls <- left_join(qtls, sdY) # add them to the table of sQTLs

#### Read in GWAS Data #### 
message("Reading in GWAS File:")
message(args$gwas_fp)

gwas <- read_tsv(args$gwas_fp)
gwas <- filter(gwas, chromosome == paste0("chr", args$chr))

# rename to match the other dataset
gwas <- 
  rename(gwas, 
         "gwas_rsid" = "variant_id",
         "id" = "panel_variant_id",
         "gwas_zscore" = "zscore", 
         "gwas_pval" = "pvalue",
         "gwas_beta" = "effect_size",
         "gwas_se" = "standard_error"
  )

gwas$gwas_maf <- sapply(gwas$frequency, function(x){min(x, 1 - x)})
gwas$gwas_type <- args$gwas_type

#### Merge GWAS and QTL into one table ####
dat <- inner_join(qtls, gwas)

#Add study type and sample size to merged dataset
dat$qtl_n <- args$qtl_n
dat$gwas_n <- args$gwas_n

# remove loci which are missing pvals or betas (I don't have SE...)
rem_na1 <- which(is.na(dat$qtl_beta) | is.na(dat$gwas_beta))
rem_na2 <- which(is.na(dat$qtl_pval) | is.na(dat$gwas_pval))
rem_na3 <- which(is.na(dat$qtl_se) | is.na(dat$gwas_se))

rem <- unique(c(rem_na1, rem_na2, rem_na3))
message("No. of SNPs removed: ", length(rem))
if (length(rem) != 0) {
  dat <- dat[-rem,]
}

# Order by phenotype_id
dat <- dplyr::arrange(dat, phenotype_id, position)

# Run coloc only for region with GWAS lead p-value < min_gwas_p
min_gwas_p <- base::tapply(dat$gwas_pval, dat$phenotype_id, min)

# Do this the dplyr way
dat %>%
  group_by(phenotype_id) %>%
  filter(gwas_pval == min(gwas_pval)) %>%
  mutate(n = n()) -> x

message(paste0("Number of loci with lead GWAS p >= ", args$min_gwas_p,": ", sum(min_gwas_p >= args$min_gwas_p)))
min_gwas_p <- min_gwas_p[min_gwas_p < args$min_gwas_p]
phenotypes <- names(min_gwas_p)

message(sprintf("Number of phenotypes to test for coloc after excluding regions with lead GWAS p > %f: %i",
 args$min_gwas_p, length(phenotypes)))

if (length(phenotypes) == 0) {
  message("No molQTLs and GWAS hits to colocalize! Done!")
  quit("no")
}

browser()

## Specify the priors from the argument parser
priors <- str_split(args$priors, ",")[[1]]
p1used <- as.numeric(priors[1])
p2used <- as.numeric(priors[2])
p12used <- as.numeric(priors[3])

message("-------- Running coloc! --------")
source("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/pipeline/coloc/functions_run_coloc_JE.R")

# Do this in parallel
fig_out="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/figures"
message("Multi-processing using ", args$cores-1, " cores")
cl <- makeCluster(args$cores-1, outfile = "")
registerDoParallel(cl)
result <- foreach(i = phenotypes, .combine = "rbind", .packages = c("coloc", "ggplot2")) %dopar%
run_coloc(phenotype = i, data = dat,
		  qtl_coloc_input = "beta", qtl_coloc_mode = "single",
		  gwas_coloc_input = "beta", gwas_coloc_mode = "single",
		  genofile = args$genofile, covfile = args$covariates,
		  locuscompare = TRUE,
		  locuscompare_fig = paste0(fig_out, "/", args$gwas_trait, "_locuscompare"),
		  locuszoom = TRUE,
		  locuszoom_fig = paste0(fig_out, "/", args$gwas_trait, "_locuszoom"),
		  plot_main = paste0(args$gwas_trait, " GWAS and ", args$qtl_trait),
		  ld_show = FALSE,
		  ld_variant = "qtl",
		  p1 = p1used, p2 = p2used, p12 = p12used,
		  sensitivity = FALSE,
		  sensitivity_fig = NULL,
		  qtl_trait = args$qtl_trait,
		  gwas_trait = args$gwas_trait, 
		  pthr = args$min_gwas_p)
stopCluster(cl)

# Save results
base_dir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/coloc_results"
out_name=sprintf("%s/%s_chr%i_%s_coloc_results.tsv", base_dir, args$qtl_trait, args$chr, args$gwas_trait)

write_tsv(result, out_name)
system(paste("gzip", out_name))
