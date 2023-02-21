#/usr/bin/bash
# -*- sh -*-

set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the results from combine_gffs_batch.sh, and combines all the tmp annotation files together
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

contigs=(
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY) 

# Set where the output goes
out_dir=$DIR

touch $out_dir/psi_detail.tsv
for c in ${contigs[@]}; do
echo $c;
cat tmp/${c}_* >> $out_dir/psi_detail.tsv;
done

# Always remove your temp files kids!
rm -r tmp