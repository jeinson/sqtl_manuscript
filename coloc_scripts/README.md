Scripts used to run colocalization analysis. Credit for Silva Kasela
(@silvakasela) for many of these scripts. 

This step of the pipeline is run in the following order:
1. `coloc.get_sdY.R`: Takes in 1) a phenotype file (normalized PSI in this case)
   2) the covariate file used for QTL mapping. 3) An interaction file, can be
   NA. 4) The number of cores, if running using `doParallel` 5) an output file
   name. It returns the standard deviation of the residuals a linear model of
   the phenotype ~ covariates + interactions. This is used in coloc. 

2. `compile_samplesheet.R`: Make a sample sheet using for running lots of coloc
   tasks in parallel. It produces a sample sheet with a line for every
   combination of trait, tissue, and chromosome. 

3. `coloc_wrapper.R`: This script is designed to run colocalization analysis of
   GTEx psi-sQTLs against GWAS hits, using p-values (as opposed to beta values).
   It can take a GWAS and QTL result file as input. 
   	- `functions_run_coloc_JE.R`: This is used heavily in the above script.
	  Contains wrappers for the actual coloc function. 

4. `submit_coloc_job_array.sh`: Runs the whole pipeline as an array of jobs on
   SLURM, using the sample sheet and wrapper functions. 
