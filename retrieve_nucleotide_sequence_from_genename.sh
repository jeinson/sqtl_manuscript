#!/bin/bash

set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes the name of a gene and a GTF file, and uses grep to get the features 
from the GTF that contain that gene, and retrieve the nucleotide sequence from hg38. 
Usage: $(basename $0) -g Gene_name
Where: -g|--gene is the name of the gene (without quotes)
" >&2
    exit 1
}

# Janky argument parsing
[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -g|--gene)
            GENE="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

GTEX_GTF=~/splicing_project/gencode.v26.GRCh38.GTExV8.genes.gtf
GENOME=~/splicing_project/Homo_sapiens_assembly38.fasta

# First, get a gtf representation of the file
export GTF=$(cat <(head -6 $GTEX_GTF) <(grep \"$GENE\" $GTEX_GTF | grep exon))

if [[ $(wc -l <<< "$GTF") -eq 6 ]]; then
	echo -e "Gene not found. Please check spelling" >&2; 
fi
# Then use bedtools to pull out the amino acid sequence
#module load bedtools/2.29.0
bedtools getfasta -fi $GENOME -bed <(echo "$GTF") -split | fold -w 60 > output.txt
