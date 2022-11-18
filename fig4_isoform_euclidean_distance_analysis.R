# Make a density plot of Euclidean distances between spliced in and spliced 
# out isoforms among genes with a sQTL-GWAS colocalization event. 

rm(list = ls())
source("myPackages.R")
source("sqtl_manuscript_functions.R")

dat <- read_csv("data/top_sQTLs_with_top_coloc_with_AF_new.csv")
dat <- filter(dat, !is.na(EUCL.DISTANCE))

fig4A <- ggplot(dat, aes(EUCL.DISTANCE)) + 
  geom_histogram(fill = "cornflowerblue", color = "black") + 
  #geom_density() +
  sqtl_manuscript_theme() + 
  xlab("Euclidean Distance between Isoforms")

save_plot("fig4A_euclidean_distance_histogram.svg", width = 8, height = 5)
fig4A
dev.off()
