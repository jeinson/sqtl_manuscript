library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(data.table)

library(here)

source(here("auxillary_scripts/MedianBootstrap.R"))
source(here("auxillary_scripts/FlagOutlier.R"))
source(here("auxillary_scripts/NiceFigures.R"))
source(here("sqtl_manuscript_functions.R"))
## Random Shuffle
data_full = read.csv(here("data/top_sQTLs_MAF05_with_AF.csv"))
data_full_1 = read.csv(here("data/cross_tissue_nonsignificant_genes_with_AF.csv"))
data_full_1 %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full_1
data_full = data_full[!flag_outliers(data_full$ALIGN.LENGTH) & data_full$ALIGN.LENGTH >= 15,]
data_full_1 = data_full_1[!flag_outliers(data_full_1$ALIGN.LENGTH) & data_full_1$ALIGN.LENGTH >= 15,]
data_full %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full
data2 = data_full[data_full$mean_01_psi < 0.5, ]

data = data_full[data_full$mean_01_psi > 0.5, ]

data = data[sample(1:nrow(data), nrow(data2)),]
data_full = rbind(data, data2)
data_full$mean_01_psi = data_full[sample(1:nrow(data_full)),]$mean_01_psi
data2 = data_full[data_full$mean_01_psi < 0.5, ]
data = data_full[data_full$mean_01_psi > 0.5, ]




data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = paste0("highly", "\n", nrow(data))
data2$GROUP = paste0("lowly", "\n", nrow(data2))

to_draw = rbind(data, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p11 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE,  draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40)) +
        ylab("Q1")
p21 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE,  draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25)) + 
        ylab("Q3_pLDDT")




p31 = ggviolin(to_draw, x="GROUP", y="ASN..",  rug = TRUE, draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + ylab("ASN.." ) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p41 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE, draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + ylab("CYS.." ) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p51 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none") +
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..))) + xlab("Random Shuffle") + ylab("ALIGN.LENGTH") + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))



## Highly vs. lowly
data_full = read.csv(here("data/top_sQTLs_MAF05_with_AF.csv"))
data_full_1 = read.csv(here("data/cross_tissue_nonsignificant_genes_with_AF.csv"))


inter = intersect(colnames(data_full), colnames(data_full_1))
data_full  = data_full[, inter]
data_full_1 = data_full_1[,inter]

data_full = rbind(data_full, data_full_1)

data2 = data_full[data_full$mean_01_psi < 0.5, ]

data = data_full[data_full$mean_01_psi > 0.5, ]


data = unique.data.frame(data)
data2 = unique.data.frame(data2)

data = data[!flag_outliers(data$ALIGN.LENGTH) & data$ALIGN.LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$ALIGN.LENGTH) & data2$ALIGN.LENGTH >= 15,]

data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = paste0("highly", "\n", nrow(data))
data2$GROUP = paste0("lowly", "\n", nrow(data2))

to_draw = rbind(data, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p12 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40))
p22 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))
p32 = ggviolin(to_draw, x="GROUP", y="ASN.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p42 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE, color = "lightgray", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p52 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.y = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..)))  + xlab("Highly vs. lowly") + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))


## sQTL vs Non sQTL



data = read.csv(here("data/top_sQTLs_with_top_coloc_with_AF_new.csv"))
data2 = read.csv(here("data/cross_tissue_nonsignificant_genes_with_AF.csv"))

data = unique.data.frame(data)
data2 = unique.data.frame(data2)


data = data[!flag_outliers(data$ALIGN.LENGTH) & data$ALIGN.LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$ALIGN.LENGTH) & data2$ALIGN.LENGTH >= 15,]
nrow(data)
nrow(data2)


data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = paste0("sQTL", "\n", nrow(data))
data2$GROUP = paste0("non_sQTL", "\n", nrow(data2))
to_draw = rbind(data, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p13 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE,  draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40))
p23 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE,  draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))
p33 = ggviolin(to_draw, x="GROUP", y="ASN.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p43 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE,  draw_quantiles = 0.5, color = "GROUP") + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p53 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, draw_quantiles = 0.5, color = "GROUP") + 
        theme(legend.position="none", axis.title.y = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..))) + xlab("sQTLs vs. non_sQTLs") + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))


## Highly sQTLs vs. lowly sQTLs

data_full = read.csv(here("data/top_sQTLs_with_top_coloc_with_AF_new.csv"))

data = data_full[data_full$has_coloc,]
data2 = data_full[!data_full$has_coloc,]


data = data[!flag_outliers(data$ALIGN.LENGTH) & data$ALIGN.LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$ALIGN.LENGTH) & data2$ALIGN.LENGTH >= 15,]


data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = paste0("coloc", "\n", nrow(data))
data2$GROUP = paste0("non_coloc", "\n", nrow(data2))

to_draw = rbind(data, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p14 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40))
p24 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))
p34 = ggviolin(to_draw, x="GROUP", y="ASN.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p44 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.x = element_blank(), 
              axis.title.y = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p54 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, color = "lightgray", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.y = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..))) + xlab("Coloc vs. non_coloc") + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))

## Increasing vs. decreasing

data  = fread(here("data/top_sQTLs_MAF05_w_anc_allele.tsv"))
data_full = read.csv(here("data/top_sQTLs_MAF05_with_AF.csv"))

data_full = merge(data_full, data, by="top_pid", all = F)

data1 = data_full[data_full$anc_allele == data_full$li_allele,]
data2 = data_full[data_full$anc_allele == data_full$hi_allele,]
nrow(data2)
data1 = data1[!flag_outliers(data1$ALIGN.LENGTH) & data1$ALIGN.LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$ALIGN.LENGTH) & data2$ALIGN.LENGTH >= 15,]

data1$GROUP = paste0("Increasing", "\n", nrow(data1))
data2$GROUP =  paste0("Decreasing", "\n", nrow(data2))

data1 = data1[colnames(data1) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data1)]

to_draw = rbind(data1, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p15 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE, draw_quantiles = 0.5, color = "GROUP") + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40))
p25 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE, draw_quantiles = 0.5, color = "GROUP") + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) + 
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))
p35 = ggviolin(to_draw, x="GROUP", y="ASN.." ,  rug = TRUE, draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p45 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE,  draw_quantiles = 0.5, color = "lightgray") + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank())  +
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p55 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        theme(legend.position="none", axis.title.y = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +  
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..))) + xlab("Incr vs. decr") + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))

# Variable vs. Constitutive

data_full = read.csv(here("data/top_sQTLs_MAF05_with_AF.csv"))
data_full_1 = read.csv(here("data/cross_tissue_nonsignificant_genes_with_AF.csv"))
data_full_1 %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full_1
data_full = data_full[!flag_outliers(data_full$ALIGN.LENGTH) & data_full$ALIGN.LENGTH >= 15,]
data_full_1 = data_full_1[!flag_outliers(data_full_1$ALIGN.LENGTH) & data_full_1$ALIGN.LENGTH >= 15,]
data_full %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full
data_full_1 = read.csv("data/cross_tissue_constitutive_exons_with_AF.csv")

data2 = data_full_1

data = data_full


data = unique.data.frame(data)
data2 = unique.data.frame(data2)

data = data[!flag_outliers(data$ALIGN.LENGTH) & data$ALIGN.LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$ALIGN.LENGTH) & data2$ALIGN.LENGTH >= 15,]

data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = paste0("variable", "\n", sum(!is.na(data$ALPHAFOLD.NAME)))
data2$GROUP = paste0("constitutive", "\n", sum(!is.na(data2$ALPHAFOLD.NAME)))

to_draw = rbind(data, data2)
to_draw$ASN.. = to_draw$ASN.. / to_draw$ALIGN.LENGTH
to_draw$CYS.. = to_draw$CYS.. / to_draw$ALIGN.LENGTH

p16 = ggviolin(to_draw, x="GROUP", y="Q1",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        stat_compare_means(label.y = 120, aes(label = paste0("p =", ..p.format..))) + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +   
        scale_y_continuous(limits = c(0, 130), breaks = seq(0, 120, by = 40))
p26 = ggviolin(to_draw, x="GROUP", y="Q3_pLLDT",  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        stat_compare_means(label.y = 0, aes(label = paste0("p =", ..p.format..))) + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) + 
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))
p36 = ggviolin(to_draw, x="GROUP", y="ASN.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) + 
        scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p46 = ggviolin(to_draw, x="GROUP", y="CYS.." ,  rug = TRUE, color = "GROUP", draw_quantiles = 0.5) + 
        stat_compare_means(label.y = 0.27, aes(label = paste0("p =", ..p.format..))) + 
        theme(legend.position="none", axis.title.y = element_blank(),  
              axis.title.x = element_blank(), axis.text.x = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) + 
        scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0.0, 0.3, by = 0.1))
p56 = ggviolin(to_draw, x="GROUP", y="ALIGN.LENGTH",  rug = TRUE, color = "lightgray", draw_quantiles = 0.5) + 
        stat_compare_means(label.y = 110, aes(label = paste0("p =", ..p.format..)))  + xlab("Variab vs. Const") + 
        theme(legend.position="none", axis.title.y = element_blank(), 
              axis.line.y = element_blank(), axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank()) +  
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 25))


save_plot("fig3_exon_feature_comparisons.svg", width = 6, height = 6)
grid.arrange(p11, p12, p13, p14, p15, p16, p21, p22, p23, p24, p25, p26, p31, p32, 
             p33, p34, p35, p36, p41, p42, p43, p44, p45, p46, p51, p52, p53, p54, p55, p56, ncol=6, nrow=5,
             widths = c(1.5,1,1,1,1,1), 
             heights = c(1,1,1,1,1.5))
dev.off()

