#!/usr/bin/bash

gtex_vcf=$1
tiss=$2
i=$3
chr=$4
date=$5

qtltools cis \
	--vcf $gtex_vcf \
	--bed bed/${tiss}_normalized_psi_grouped_exons_qtltools.bed.gz \
	--include-site ../data/common_variants/common_variants_only_chr$i.txt \
	--cov covariates/${tiss}_covariates.tsv  \
	--permute 100 1000 \
	--region $chr \
	--grp-best \
	--out qtltools_results/$out_dir/${tiss}.qtltools.grouped.permutation.MAF05.$chr.$date.txt 
