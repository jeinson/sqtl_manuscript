# -*- sh -*-
#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=4:00:00
#$ -l h_vmem=1G
#$ -N psi_filter
#$ -m ea
#$ -M jeinson@nygenome.org

source ~/python3_env/bin/activate

set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the results from combine_gffs_batch.sh
Usage: $(basename $0) -d directory
Where: -d|--dir is the directory containing the gff files
" >&2
    exit 1
}

# Janky argument parsing
[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -d|--dir)
            DIR="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done


python psi_filterer.py \
	--input $DIR/all.A.psi.tsv \
	--output $DIR/all.A.psi.filtered.tsv
