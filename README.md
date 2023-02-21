This repo contains code used in [preprint](https://www.biorxiv.org/content/10.1101/2022.12.05.518915v1.article-metrics) *The impact of genetically controlled splicing on exon inclusion and protein structure*. 

The scripts are organized by the figure they generate. All data to run these scripts is available at https://zenodo.org/record/7275062. 

Other directories included are for: 
- `psi_scripts`: Scripts used for handling raw PSI, which starts with GTEx v8
  BAM files. We use a modified version of
  [IPSA-nf](https://github.com/guigolab/ipsa-nf) for this task. 
- `sqtl_calling_scripts`: Scripts used to take PSI and run QTL calling across 18
  GTEx tissues, then merge these QTLs by top effect size. 
- `coloc_scripts`: Scripts used to perform colocalization analysis of sQTLs
  against GWAS summary statistics from [Barbeira et
  al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02252-4)
- `susie_scripts`: Scripts used to run fine mapping on sQTLs with SuSiE.
- `data`: Small data snippets and final output files used to make the plots in
  the root directory. For larger files, refer to the Zenodo repo. 

Disclaimer: this repo is a living document! If anything is missing or unclear,
please feel free to reach out to me, Jonah Einson, at jeinson@nygenome.org. 
