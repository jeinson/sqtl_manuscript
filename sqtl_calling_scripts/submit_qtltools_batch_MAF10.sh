#!/bin/bash
set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the results from sQTL_results_combiner.R, and runs sQTL calling on violently normalized PSI scores
Where: -t|--tiss is the name of the tissue 
" >&2
    exit 1
}

# Janky argument parsing
[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -t|--tiss)
            tiss="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

# Set te paths to important stuff
gtex_vcf="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz" 
date=$(printf '%(%m.%d.%y)T\n' -1)

# Check to make sure everything is where it's supposed to be
if [ ! -f bed/${tiss}_normalized_psi_grouped_exons_qtltools.bed.gz ]
then
    echo "Cannot find the input file bed/${tiss}_normalized_psi_grouped_exons_qtltools.bed.gz. Please check the file path"
    exit
fi

if [ ! -f covariates/${tiss}_covariates.tsv ]
then 
    echo "Cannot find the covariate file covariates/${tiss}_covariates.tsv. Please make sure it exists"
    exit
fi

# Set the output directory, and delete it if it already exists
out_dir="${tiss}_MAF10"

function pause(){
	read -p "$*"
}

if [ ! -d qtltools_results/$out_dir ]
then 
# 	pause 'Output directory already exists. Press [enter] to delete it and continue'
# 	rm -r qtltools_results/$out_dir
#     mkdir qtltools_results/$out_dir
# else
	echo "Creating new output directory qtltools_results/${out_dir}"
	mkdir qtltools_results/$out_dir
fi

# Actually run QTLtools
module load qtltools
echo "Submitting QTLtools jobs"

# Do this whole thing as a slurm grid array

sbatch --job-name sQTLM10_${TISS}_%A -t 1:00:00 --mem=10G \
-e slurm_logs/${tiss}.qtltools.%a.$date.stderr.txt \
-o slurm_logs/${tiss}.qtltools.%a.$date.stdout.txt \
--array=1-22 \
--wrap=" 
qtltools cis \
--vcf $gtex_vcf \
--bed bed/${tiss}_normalized_psi_grouped_exons_qtltools.bed.gz \
--include-site ../data/common_variants/MAF_gt_.1_variants_only_chr\${SLURM_ARRAY_TASK_ID}.txt \
--cov covariates/${tiss}_covariates.tsv  \
--permute 100 1000 \
--region "chr"\${SLURM_ARRAY_TASK_ID} \
--grp-best \
--out qtltools_results/$out_dir/${tiss}.qtltools.grouped.permutation.MAF10.chr\${SLURM_ARRAY_TASK_ID}.$date.txt"

