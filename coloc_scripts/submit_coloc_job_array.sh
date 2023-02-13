#!/bin/bash

#SBATCH --array=35001-40000
#SBATCH --job-name=coloc_run
#SBATCH --output=log_files/%A_%a.out
#SBATCH --error=log_files/%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

SAMPLE_SHEET=/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/coloc_v8_anno/sample_sheets/sQTL_gwas_coloc_gtex_final_list.tsv
N_RUNS=$(wc -l $SAMPLE_SHEET)

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET)
#LINE=$(sed -n 9p $SAMPLE_SHEET)
read -r -a array <<< $LINE

TISS=${array[0]}
CHR=${array[1]}
QTL_FP=${array[2]}
QTL_N=${array[3]}
COV_FP=${array[4]}
SDY_FP=${array[5]}
GWAS_FP=${array[6]}
GWAS_NAME=${array[7]}
GWAS_N=${array[8]}



# Print the task id.

# Run the thing!
srun /gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/pipeline/coloc/coloc_wrapper_JE.R \
--qtl_trait $TISS \
--gwas_trait $GWAS_NAME \
--qtl_fp $QTL_FP \
--sdY $SDY_FP \
--covariates $COV_FP \
--gwas_fp $GWAS_FP \
--chr $CHR \
--gwas_n $GWAS_N \
--qtl_n $QTL_N \
--min_gwas_p 1e-6 \
--cores 4

