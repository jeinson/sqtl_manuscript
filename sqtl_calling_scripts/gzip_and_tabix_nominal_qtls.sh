#!/bin/bash -l
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=5210
#SBATCH -J gzip_and_tabix
#SBATCH --array=396

module load htslib/1.7

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p completed_tissues_and_chrom.txt)
read -r -a line_array <<< "$LINE"

file=qtltools_results/${line_array[1]}_nominal/${line_array[1]}.qtltools.nominal.MAF05.chr${line_array[0]}.06.02.21.txt

gzip $file