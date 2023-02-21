# -*- sh -*-
#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=4:00:00
#$ -l h_vmem=50G
#$ -N combine_gffs
#$ -e logs
#$ -o logs
#$ -m ea
#$ -M jeinson@nygenome.org

set -eo pipefail

function Usage() {
    echo -e "\
Purpose: This script takes a directory containing the gff files from the A07 step of IPSA-nf, then combines them locally
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

# Tests to make sure the inputs are actually there and exist
[[ -z "${DIR}" ]] && Usage
[[ -d "${DIR}" ]] || (echo "Can't find a directory at ${DIR}" >&2; exit 1)

# Start running the thing
module load perl
export PERL5LIB=$(pwd)/Perl

declare -a filenames
for file in $DIR/A07/*; do
	filenames=("${filenames[@]}" "$file")
done

# Process each element of the array to parse out the gtex donor names
declare -a samples_names
for ((i=0; i <= ${#filenames[@]}; i++)); do
	sample_names[$i]=$(echo ${filenames[$i]} | cut -f3 -d "/" | cut -f1 -d ".")
done

# Combine the two lists together 
declare input
for ((j=0; j <= ${#filenames[@]}; j++)) do
	input="$input -i ${filenames[$j]} ${sample_names[$j]}"
done

#echo $input
# Do the thing!
perl ipsa-full/Perl/merge_gff.pl $input -o psi $DIR/all.A.psi.tsv -o cosi $DIR/all.A.cosi.tsv
