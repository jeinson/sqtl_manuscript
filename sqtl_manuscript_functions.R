### This is an R script that has some functionf to help the sqtl manuscript run smoothly


save_plot <- function(fn, save_fn = svg, ...){
  full_path <- paste0(
    "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sqtl_manuscript/figures/",
    fn)
  save_fn(full_path, ...)
  
}
  

sqtl_manuscript_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),
               axis.text=element_text(size=11), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=11), legend.title = element_text(size=12)))
}
