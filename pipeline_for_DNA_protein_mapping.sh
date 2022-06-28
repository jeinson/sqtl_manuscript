#!/bin/bash

echo Loading modules
module load R_packages/4.1.1
module load bioinfo-tools
module load samtools
module load blast

echo Building index
samtools faidx Homo_sapiens_assembly38.fasta

echo building BLAST DB
makeblastdb -in MANE.GRCh38.v1.0.ensembl_protein.faa -out MANE_prot -parse_seqids -taxid 9606 -dbtype prot

FILE=data/top_sQTLs_with_top_coloc_event.tsv
N=$(wc -l < $FILE)
N=$(($N - 1))

for (( i=1; i<=$N; i++))
# for (( i=75; i<=76; i++))
do
	echo This is round $i
	# Exon extraction
	echo Exon extraction
	Rscript get_needed_exons.R $i $FILE
	Rscript get_needed_exon.R $i $FILE

	# blastx
	echo blastx
	blastx -query exon_seq.fa -db MANE_prot -out try_blastx.txt -outfmt "6 qseqid sacc sstart send sseq length evalue score"

	# Combine final table
	echo Combine final table
	Rscript get_final_table.R

	rm try_blastx.txt
	rm exon_seq.fa
done

