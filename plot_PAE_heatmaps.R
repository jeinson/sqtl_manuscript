library(rjson)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(DescTools)
library(stringr)


plot_pae = function(data, data_2, pos1, pos2){
        mat1 = as.data.frame(cbind(data[[1]]$residue1, data[[1]]$residue2, data[[1]]$distance))
        mat1 = matrix(mat1$V3, byrow = T, nrow = max(mat1$V1))
        print(max(data_2[[1]]$residue1))
        print(dim(mat1))
        new_mat1 = matrix(NA, nrow = max(data_2[[1]]$residue1), max(data_2[[1]]$residue1))
        print(dim(new_mat1))
        if (dim(mat1)[1] < dim(new_mat1)[1] & pos2 < dim(new_mat1)[1]){
                new_mat1[1:pos1, 1:pos1] = mat1[1:pos1, 1:pos1]
                new_mat1[1:pos1, pos2:ncol(new_mat1)] = mat1[1:pos1, (pos1+1):ncol(mat1)]
                new_mat1[pos2:nrow(new_mat1), 1:pos1] = mat1[(pos1+1):nrow(mat1), 1:pos1]
                new_mat1[pos2:nrow(new_mat1), pos2:ncol(new_mat1)] = mat1[(pos1+1):nrow(mat1), (pos1+1):ncol(mat1)]
                to_plot = melt(new_mat1, na.rm = FALSE)
                p_heatmap_1 <- ggplot(to_plot, aes(x=Var1, y=Var2)) +
                        geom_tile(aes(fill=value)) +
                        scale_fill_gradient2("pae", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                        theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                        theme(legend.position = "bottom", legend.direction = "horizontal") +
                        guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) + ggtitle("Spliced out")
                #p_heatmap_1
                
                to_plot = as.data.frame(cbind(data_2[[1]]$residue1, data_2[[1]]$residue2, data_2[[1]]$distance))
                p_heatmap_2 <- ggplot(to_plot, aes(x=V1, y=V2)) +
                        geom_tile(aes(fill=V3)) +
                        scale_fill_gradient2("pae", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                        theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                        theme(legend.position = "bottom", legend.direction = "horizontal") +
                        guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) + ggtitle("Spliced in")
                #p_heatmap_2
                
                
                
                lay <- rbind(c(1,2),
                             c(1,2))
                return(grid.arrange(grobs = list(p_heatmap_1, p_heatmap_2), layout_matrix = lay))
        }
}


plot_pae_corr = function(data, data_2, pos1, pos2){
        mat1 = as.data.frame(cbind(data[[1]]$residue1, data[[1]]$residue2, data[[1]]$distance))
        mat2 = as.data.frame(cbind(data_2[[1]]$residue1, data_2[[1]]$residue2, data_2[[1]]$distance))
        mat1 = matrix(mat1$V3, byrow = T, nrow = max(mat1$V1))
        mat2 = matrix(mat2$V3, byrow = T, nrow = max(mat2$V1))
        new_mat1 = matrix(NA, nrow = nrow(mat2), ncol=ncol(mat2))
        if (dim(mat1)[1] < dim(new_mat1)[1] & pos2 < dim(new_mat1)[1]){
                new_mat1[1:pos1, 1:pos1] = mat1[1:pos1, 1:pos1]
                new_mat1[1:pos1, pos2:ncol(new_mat1)] = mat1[1:pos1, (pos1+1):ncol(mat1)]
                new_mat1[pos2:nrow(new_mat1), 1:pos1] = mat1[(pos1+1):nrow(mat1), 1:pos1]
                new_mat1[pos2:nrow(new_mat1), pos2:ncol(new_mat1)] = mat1[(pos1+1):nrow(mat1), (pos1+1):ncol(mat1)]
                cor_mat = new_mat1 - mat2
                dist = sqrt(sum(cor_mat^2, na.rm = T))/ncol(mat1)
                print(dist)
                melted_cor_mat = melt(cor_mat)
                p_heatmap <- ggplot(melted_cor_mat, aes(x=Var1, y=Var2)) +
                        geom_tile(aes(fill=value)) +
                        scale_fill_gradient2("Entropy", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                        theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                        theme(legend.position = "bottom", legend.direction = "horizontal") +
                        guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5))
                p_heatmap
                return(list(p_heatmap, dist))
        }
}


setwd("~/splicing_project/Coloc_predictions")
dirs = list.dirs(path=".")
dirs = dirs[2:length(dirs)]
dirs = unique(sapply(dirs, function(x){
        unlist(str_split(unlist(str_split(unlist(str_split(x, pattern="./"))[2], pattern="[.]"))[1], pattern="_"))[1]
}))
dirs
metadata = read.csv("./data/top_sQTLs_with_top_coloc_with_AF.csv")
metadata$EUCL.DISTANCE = NA
for (dir in dirs){
        start = metadata[metadata$ALPHAFOLD.NAME == dir,]$START
        end = metadata[metadata$ALPHAFOLD.NAME == dir,]$END
        files = c(list.files(path=paste0("./", dir, ".result/"), pattern="predicted_aligned_error_v1.json", full.names = TRUE)[1],
                  list.files(path=paste0("./", dir, "_mutant.result/"), pattern="predicted_aligned_error_v1.json", full.names = TRUE)[1])
        print(files)
        print(start)
        print(end)
        if (!is.na(files[1]) & !is.na(files[2]) & metadata[metadata$ALPHAFOLD.NAME == dir,]$MAPPED.TO.OTHER == "False"){
                data_2 <- fromJSON(file = files[1])
                data = fromJSON(file=files[2])
                plot = plot_pae(data, data_2, start[1], end[1])
                #plot saving
                # ggsave(plot = plot, filename = paste0(dir, "_no_", start, "_", end, "_PAE_new.png"), path = "./AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
                plot = plot_pae_corr(data, data_2, start[1], end[1])
                #plot saving
                # ggsave(plot = plot[[1]], filename = paste0(dir, "_no_", start, "_", end, "_PAE_corr_new.png"), path = "./AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
                metadata[metadata$ALPHAFOLD.NAME == dir,]$EUCL.DISTANCE = plot[[2]]
        }
        
}

write.csv(metadata, "./data/top_sQTLs_with_top_coloc_with_AF_new.csv")

metadata[which(metadata$EUCL.DISTANCE == max(metadata$EUCL.DISTANCE, na.rm=T)),]

