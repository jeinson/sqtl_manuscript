#!/usr/bin/bash

set -eo pipefail

source /gpfs/commons/home/jeinson/python3_env/bin/activate

function Usage() {
    echo -e "\
Purpose: This script submits coloc.get_sdY.R for a specified tissue. 
Usage: $(basename $0) -t Tissue name
Where: -t|--tiss is the name of the tissue 
" >&2
    exit 1
}

# Janky argument parsing
[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -t|--tiss)
            TISS="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

# Set arguments for coloc.get_sdY.R
basedir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance"
PHENO_FILE="${basedir}/sQTL_v8_anno/bed/${TISS}_normalized_psi_grouped_exons_qtltools.bed.gz"
COV_FILE="${basedir}/sQTL_v8_anno/covariates/${TISS}_covariates.tsv"
INTERACT_FILE="NA" # No interactions in this analysis
CORES=5 # To run things in parallel
OUT_FILE="${basedir}/coloc_v8_anno/sdY/${TISS}_normalized_PSI_sdY.tsv"

# check to make sure these files exist
if [[ ! -f $PHENO_FILE || ! -f $COV_FILE ]]; then
  echo "Make sure $PHENO_FILE and $COV_FILE exist"
  exit
fi


module load R/4.0.0

# Submit slurm job
sbatch --mem=16g --time=0:30:00 --cpus-per-task=${CORES} \
-o ${basedir}/slurm_logs/${TISS}_sdY.out -e ${basedir}/slurm_logs/${TISS}_sdY.err \
${basedir}/pipeline/coloc/coloc.get_sdY.R \
$PHENO_FILE $COV_FILE $INTERACT_FILE $CORES $OUT_FILE