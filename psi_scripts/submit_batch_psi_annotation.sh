#/usr/bin/bash
# -*- sh -*-

set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the results from combine_gffs_batch.sh, after filtering and sorting
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
chr19
chr20
chr21
chr22
chrX
chrY)

mkdir tmp

for contig in ${contigs[@]};
do 
    echo $contig;

    echo "source ~/python3_env/bin/activate; \
    python psi_annotation/psi_annotation_tool.py \
    --psi_file $DIR/all.A.psi.filtered.sorted.tsv \
    --annotation_gfx psi_annotation/Homo_sapiens_GRCh38.97.processed.anno.gtf.gz  \
    --contig $contig > \
    tmp/${contig}_psi_detail.tsv" | 
    qsub -cwd -V -l h_rt=0:5:00 -l h_vmem=1G -N psi_annotation_$contig -o logs -e logs

done
