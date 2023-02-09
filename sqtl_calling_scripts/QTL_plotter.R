### This script has a function for plotting QTLs given a VCF, phenotype, 
### a list of samples, and a list of corresponding variants

#' Genotype info function
get_gt_info <- function(sid, vcf, samp_names = NULL){
  # Get header info
  if(is.null(samp_names)){
    samp_names <- 
      system(paste("zcat", vcf, "| head -4000 | grep 'CHROM'"), 
             intern = T)
    samp_names <- strsplit(samp_names, "\t")[[1]][-(1:9)]
  }
  
  # Get genotype info for specified region
  chr = strsplit(sid, "_")[[1]][1]
  start = as.numeric(strsplit(sid, "_")[[1]][2])
  end = start + 1
  region <- paste0(chr, ":", start, "-", end)
  
  out <- system(paste(" /nfs/sw/htslib/htslib-1.7/devel/tabix ", vcf, region), intern = T)
  
  if(length(out) == 0) stop(message("This SNP is not in the VCF file"))
  
  out <- strsplit(out, "\t")[[1]]
  out <- map_chr(strsplit(out[-(1:9)], ":"), 1)
  names(out) <- samp_names
  
  return(out)
}

#' QTL plotting function
plot_QTL <- function(vcf_path, qtlVariant, phenotype_matrix) {
  samp_names <- 
    system(paste("zcat", vcf_path, "| head -4000 | grep 'CHROM'"), 
           intern = T)
  samp_names <- strsplit(samp_names, "\t")[[1]][-(1:9)]
  
  sVariants <- 
    PSI_sqtls %>% 
    filter(group %in% gene_categories[[m]]) %>%
    .$vid
  # p.values <- 
  #   PSI_sqtls %>% 
  #   filter(group %in% gene_categories[[m]]) %>%
  #   .$qval 
  # sExons <- 
  #   PSI_sqtls %>% 
  #   filter(group %in% gene_categories[[m]]) %>%
  #   .$top_pid
  # sGenes <- 
  #   PSI_sqtls %>% 
  #   filter(group %in% gene_categories[[m]]) %>%
  #   .$group
  # leafcutter.p.values <- # Get these in the right order for comparison
  #   J.tools::two_col_2_named_vector(leafcutter_sqtls$qval, 
  #                                   J.tools::remove_trailing_digit(leafcutter_sqtls$group_id)
  #   )
  leafcutter.p.values <- leafcutter.p.values[sGenes]
  
  # Make a nice function for plotting the number of exons available
  stat_box_data <- function(y, upper_limit = 1.07){
    return( 
      data.frame(
        y = 0.95 * upper_limit,
        label = paste(length(y))
      )
    )
  }
  
  # pdf(filenames[m], width=3, height=3)
  # for(k in 1:length(p.values)){
  #   message(k)
  gt <- 
    get_gt_info(qtlVariant, vcf_path, samp_names) %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    set_colnames(c("SampleName", "Genotype")) %>%
    mutate(Genotype = as.factor(Genotype))
  
  psi <- 
    t(PSI[which(PSI$ID == sExons[k]),][-(1:4)]) %>%
    as.data.frame %>%
    rownames_to_column() %>%
    set_colnames(c("SampleName", "PSI"))
  
  plot_tbl <- inner_join(gt, psi)
  plot_tbl <- plot_tbl[complete.cases(plot_tbl),]
  
  # Save the figure
  print(
    ggplot(plot_tbl, aes(Genotype, PSI)) +
      geom_boxplot(fill = "grey") +
      stat_summary(
        fun.data = stat_box_data, 
        geom = "text"
      ) +
      ggtitle(paste0("Exon: ", sExons[k], "\n", 
                     "Gene: ", sGenes[k], "\n",
                     "Variant: ", sVariants[k], 
                     "\nq = ", round(p.values[k], 8),
                     "\nleafcutter q = ", round(leafcutter.p.values[k], 8))) +
      theme_classic() +
      theme(plot.title = element_text(size = 8))
  )
  # }
  # dev.off()
}