This step of the pipeline is run in the following order:

1. `normalize_psi.R`: Take in the output from PSI calling. Annotate exons, and
   perform inverse normal transformation. Saves a version of the file called
   `{tissue}_normalized_psi_grouped_exons_qtltools.bed` which is used as input
   for QTL calling in grouped permutation mode, where genes define exon groups. 
   	- `median_psi_calculator.R`: Gets the median exon PSI across all
	individuals, using the output from (1). This is used in some downstream
	analyses. 
2. `compute_PEER_factors.R`: Use the output from (1) to calculate PEER factors,
   which are used to control for population stratification when performing QTL
   analysis `submit_PEER_slurm.sh` is a wrapper used to running this on the
   SLURM load scheduling system. 
3. `submit_qtltools_batch.sh`: Takes output from (1) and (2), and submits a
   QTLtools run. `submit_qtltools_batch_no_cutoff.sh` takes the same input and
   submits a nominal QTL calling pass. In both scripts, QTL calling is run one
   tissue at a time. To combine, use `sQTL_results_combiner.R`
4. `delta_psi_calculator_all_exons.R`: Uses the significant sQTLs output from
   (3), and then reads in the PSI files from (1) to calculate delta PSI. i.e.
   the change in PSI induced by a sQTL. It saves another file called
   `{tissue}.chr{X}.MAF05_all_exon_delta_psi.tsv`, to reference the effect size
   of the QTL. 
5. `combine_sQTLs_by_top_effect.R`: Takes the output files from the previous
   step, run across all tissues, and puts together a single file that represents
   the top sQTLs based on their delta PSIs. 


Other scripts: 
- `QTL_plotter.R`: A script which can plot a simple QTL plot, given a VCF and a
  phenotype file. Useful for making a shiny app. 
