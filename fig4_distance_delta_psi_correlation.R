library(ggplot2)
library(ggExtra)
library(stringr)
library(data.table)

source("sqtl_manuscript_functions.R")

gnomad = fread("./data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
head(gnomad)

data2 = read.csv("./data/top_sQTLs_with_top_coloc_with_AF_new.csv")
data2 = data2[data2$has_coloc ==T,]
gnomad_needed = gnomad[gnomad$gene_id %in% data2$GeneID ]
nrow(data2)
nrow(gnomad_needed)

data2 = merge(data2, gnomad_needed, by.x="GeneID", by.y="gene_id", all.x=T)


data = read.csv("./data/top_sQTLs_MAF05_with_AF.csv")
data[data$ALPHAFOLD.NAME %in% data2$ALPHAFOLD.NAME,] -> new_data
new_data = new_data[!is.na(new_data$ALPHAFOLD.NAME),]
delta_psi = c()
data2 = data2[!is.na(data2$ALPHAFOLD.NAME) & data2$ALPHAFOLD.NAME != "",]

for (name in data2$ALPHAFOLD.NAME){
        if (name %in% new_data$ALPHAFOLD.NAME){
                delta_psi = c(delta_psi, new_data[new_data$ALPHAFOLD.NAME == name,]$delta_psi[1])
        } else {
                delta_psi = c(delta_psi, NA)
        }
}

data2$delta_psi = delta_psi
data2$EUCL.DISTANCE = as.numeric(data2$EUCL.DISTANCE)
data_sub = data2[!is.na(data2$EUCL.DISTANCE) & !is.na(data2$delta_psi),]
fit = lm(EUCL.DISTANCE~delta_psi, data_sub)


plot = ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        annotate("text",label=paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                           "Intercept =",signif(fit$coef[[1]],5 ),
                           " Slope =",signif(fit$coef[[2]], 5), "\n",
                           " P =",signif(summary(fit)$coef[2,4], 5),
                           "cor =", signif(cor(data_sub$delta_psi, 
                                               data_sub$EUCL.DISTANCE),5)), x=0.25, y=13.0) +
        sqtl_manuscript_theme()

plot = ggExtra::ggMarginal(plot, type = "histogram")
ggsave(plot = plot, filename = paste0("delta_psi_distance_correlation.png"), path = "~/splicing_project/Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


plot = ggplot(data = data2, aes(x=as.numeric(EUCL.DISTANCE))) + geom_density() +
        sqtl_manuscript_theme() +
        labs(title=paste("min = ",signif(min(data_sub$EUCL.DISTANCE), 5),
                   " median =",signif(median(data_sub$EUCL.DISTANCE),5 ),
                   " mean =",signif(mean(data_sub$EUCL.DISTANCE), 5),
                   " max =",signif(max(data_sub$EUCL.DISTANCE), 5)))


ggsave(plot = plot, filename = paste0("structure_distance_distribution.png"), path = "~/splicing_project/Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)



fit = lm(EUCL.DISTANCE~oe_lof_upper, data_sub)
plot = ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        annotate("text",label=paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                    "Intercept =",signif(fit$coef[[1]],5 ),
                                    " Slope =",signif(fit$coef[[2]], 5), "\n",
                                    " P =",signif(summary(fit)$coef[2,4], 5),
                                    "cor =", signif(cor(data_sub$oe_lof_upper, 
                                                        data_sub$EUCL.DISTANCE),5)), x=1.0, y=13.0) +
        sqtl_manuscript_theme()
plot = ggExtra::ggMarginal(plot, type = "histogram")
plot
ggsave(plot = plot, filename = paste0("loeuf_distance_correlation.png"), path = "~/splicing_project/Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


data_sub$length = sapply(data_sub$COORDS, function(x){
        temp = unlist(str_split(x, pattern="-"))
        print(temp)
        (as.numeric(temp[2]) - as.numeric(temp[1]))/3
})

fit = lm(EUCL.DISTANCE~ALIGN.LENGTH, data_sub)

plot = ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        annotate("text",label=paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                    "Intercept =",signif(fit$coef[[1]],5 ),
                                    " Slope =",signif(fit$coef[[2]], 5), "\n",
                                    " P =",signif(summary(fit)$coef[2,4], 5),
                                    "cor =", signif(cor(data_sub$oe_lof_upper, 
                                                        data_sub$EUCL.DISTANCE),5)), x=50.0, y=13.0) +
        gtex_v8_figure_theme()
plot = ggExtra::ggMarginal(plot, type = "histogram")
plot
ggsave(plot = plot, filename = paste0("align_length_distance_correlation.png"), path = "~/splicing_project/Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


