This step of the pipeline is run in the following order:

1. `normalize_psi.R`: Take in the output from PSI calling. Annotate exons, and
   perform inverse normal transformation. Saves a version of the file called
   `{tissue}_normalized_psi_grouped_exons_qtltools.bed` which is used as input
   for QTL calling in grouped permutation mode, where genes define exon groups. 
   	1a. `median_psi_calculator.R`: Gets the median exon PSI across all
	individuals, using the output from (1). This is used in some downstream
	analyses. 
2. `compute_PEER_factors.R`: Use the output from (1) to calculate PEER factors,
   which are used to control for population stratification when performing QTL
   analysis `submit_PEER_slurm.sh` is a wrapper used to running this on the
   SLURM load scheduling system. 
3. `submit_qtltools_batch.sh`: Takes output from (1) and (2), and submits a
   QTLtools run. 
