#!/bin/bash -l
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=10240
#SBATCH -J compute_PEER

#$ -cwd
#$ -V
#$ -l h_rt=4:00:00
#$ -l h_vmem=10G
#$ -N calculate_PEER
#$ -m ea
#$ -M jeinson@nygenome.org

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the results from sQTL_results_combiner.R, and runs the PEER calculation workflow, "compute_PEER_factors.R"
Usage: $(basename $0) -t Tissue name -g GTEx tissue name
Where: -t|--tiss is the name of the tissue 
       -g|--gtex_tiss is the name of the tissue in GTEx, for matching existing technical covariates
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
        -g|--gtex_tiss)
            GTEX_TISS="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

module load R/3.4.1

Rscript --vanilla /gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/pipeline/sQTL/compute_PEER_factors.R $TISS $GTEX_TISS
