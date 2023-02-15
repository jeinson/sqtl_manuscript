In this analysis, we use the scripts for running susie from https://github.com/eQTL-Catalogue/susie-workflow. The wrapper script submit_gtex_sqtl_susie_run.sh specifies input files, and submits them as a batch job. 

After running susie on all sQTLs across all GTEx tissues, the steps for collapsing and analyzing credible sets are as follows:

+ credible_set_analysis.R - Reads in raw data and collapses variants across tissues
+ tissue_by_tissue_cs_overlap.R - Find overlaps between eQTLs and sQTLs, and save them
to the disk in the 'credible_set_overlaps' directory. ! Filter to genes with fewer
than 30 variants in the credible set, and with the max(abs(Zscore)) > 3, as suggested
in Kerimov et. al. 
+ analyze_overlap_signals_simple.R - Does some analysis on the overlaps between
eQTLs and sQTLs. This may be the only figure included in the analysis.

+ collapse_eQTLs.R - Use the same collapsing procedure to make a set of cross tissue
credible sets from GTEx eQTLs. 
+ cross_tissue_CS_overlaps_upset.R - Take the results
