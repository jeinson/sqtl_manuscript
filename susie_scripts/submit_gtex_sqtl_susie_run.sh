#!/bin/bash
#SBATCH --job-name=run_susie
#SBATCH --mem='24G'
#SBATCH --time=03:00:00 
#SBATCH --output=logs/%A-%a.out
#SBATCH --error=logs/%A-%a.err
#SBATCH --array=181,182,183,187,191,193,194,195,196,197,217,225,227,231,232,233

# load modules
module load R/4.0.0

# This just runs a test run of susie on GTEx PSI-sqtls. It specify the inputs needed
# to run the run_susie.R script from the gtex catalog toolset. 
TISS_FILE="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/completed_tissues_and_chrom.txt"

# Get the info from the job array
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $TISS_FILE)
read -r -a line <<< "$line"

TISS=${line[1]}
BATCH=${line[0]}

BASE_DIR="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/susie_formatted_psi_data"

# Metadata about the phenotypes (exons in this case)
# phenotype_id  group_id    gene_id chromosome  phenotype_pos   strand
PHENOTYPE_META="${BASE_DIR}/${TISS}_phenotype_meta.tsv"

# Information about the samples 
# sample_id genotype_id qtl_group
SAMPLE_META="${BASE_DIR}/${TISS}_sample_meta.tsv"

# The actual normalized expression level of each gene across all samples
# phenotype_id ...all_samples... (RN normalized)
EXPRESSION_MATRIX="${BASE_DIR}/${TISS}_normalized_psi_matrix.tsv"

# The output from QTLtools. (I'm not sure why this is called Phenotype list)
# It seems to have the top variant per gene
# molecular_trait_object_id molecular_trait_id  n_traits    n_variants  variant chromosome  position    pvalue  beta    p_perm  p_beta
QTL_BASE="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno"
PHENOTYPE_LIST="${BASE_DIR}/${TISS}_qtl_results.tsv"

# The genotypes in dosage format (0/1/2). This is the output from the first step in the pipeline.
# CHROM POS REF ALT ...all_samples...
GENOTYPE_MATRIX="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/GTEX_v8_genotypes.dose.tsv.gz"

# Covariates, in a pretty standard format. 
# SampleID  ..all_esamples...
COVARIATES="${QTL_BASE}/covariates/${TISS}_covariates.tsv"

# The name of the directory where the output is going to be saved
OUT_BASE="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/sqtl_v8_anno_susie_output"
OUT_PREFIX="${OUT_BASE}/${TISS}_susie"

# Do this in batches
N_BATCHES=22
OUT_PREFIX="${OUT_PREFIX}_${BATCH}_${N_BATCHES}"

# Other things are 
# PERMUTED (shoudl be true by default, which is what I want)

# eqtl_utils
# This should be left where it is. 


# Now actually run susie!
# Updated 'cisdistance' with the same distance used in eqtl catalogue fine mapping
susie_base="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/susie/susie-workflow"
Rscript ${susie_base}/bin/run_susie.R --expression_matrix ${EXPRESSION_MATRIX}\
    --phenotype_meta ${PHENOTYPE_META}\
    --sample_meta ${SAMPLE_META}\
    --phenotype_list ${PHENOTYPE_LIST}\
    --covariates ${COVARIATES}\
    --genotype_matrix ${GENOTYPE_MATRIX}\
    --chunk "${BATCH} ${N_BATCHES}" \
    --cisdistance 1000000 \
    --out_prefix ${OUT_PREFIX}\
    --eqtlutils "/gpfs/commons/home/jeinson/Downloads/eQTLUtils" \
    --permuted "true"

# Gzip the output
gzip -f ${OUT_PREFIX}.cred.txt
gzip -f ${OUT_PREFIX}.snp.txt
gzip -f ${OUT_PREFIX}.txt