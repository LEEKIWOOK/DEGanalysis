#heatmap
#https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html
####################################################################################
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
# install.packages("circlize")
####################################################################################
#0. Data load
library(ComplexHeatmap)
library(circlize)
library(gtools)

setwd("/home/kwlee/Projects_gflas/Team_CrisprCas/RNAexpression")
file_list = list.files(path = './data', pattern = "ppt*", full.names = TRUE)

data <- lapply(file_list,function(i){read.csv(i, header=TRUE, sep='\t')})
for(i in 1:length(data)){
  
  #!!!
  data[[i]] <- data[[i]][which(data[[i]]$X.log2FC. > 3.0),]
  #!!!
  rownames(data[[i]]) <- data[[i]]$GeneSymbol
  lenc <- length(colnames(data[[i]]))
  
  if(i <= 2){
    sample_name = c(sprintf("N-RO_%d",1:3), sprintf("M-RBO_%d",1:3))
    sample_col = list(type = c(
      "N-RO_1" = "#b2e2e2", "N-RO_2" = "#66c2a4", "N-RO_3" = "#238b45",
      "M-RBO_1" = "#cbc9e2", "M-RBO_2" = "#9e9ac8", "M-RBO_3" = "#6a51a3"))
    title = "N-RO / M-BRO"
  }else{
    sample_name = c(sprintf("RSC_",1:4), sprintf("Y79_%d",1:5))
    title = "RSC / Y79"
    sample_col = list(type = c( 
      "RSC_1" = "#b2e2e2", "RSC_2" = "#66c2a4", "RSC_3" = "#2ca25f", "RSC_4" = "#006d2c",
      "Y79_1" = "#dadaeb", "Y79_2" = "#bcbddc", "Y79_3" = "#9e9ac8", "Y79_4" = "#756bb1", "Y79_5" = "#54278f"))
  }
  
  mat <- as.matrix(data[[i]][8:lenc])
  mat_scaled <- t(apply(mat, 1, scale))
  colnames(mat_scaled)<-sample_name
  info<-data[[i]][c(2,3,6,7)]
  signif <- stars.pval(info$P.Value)
  
  #print(summary(mat_scaled))
  #htop = HeatmapAnnotation(type = sample_name, col = sample_col, show_annotation_name = FALSE, show_legend = FALSE)
   
  hright = rowAnnotation(
    signif = anno_text(signif, location = 0.5, show_name = TRUE),
    log2FC = anno_barplot(info$log2FC, bar_width = 1, gp = gpar(col = "white", fill = "#FFE200"), border = FALSE, width = unit(2, "cm")),
    show_annotation_name = TRUE, gap = unit(8, "mm"), annotation_name_rot = 0
  )
  
  col_main_func = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  hmain <- Heatmap(mat_scaled, col = col_main_func,
                   cluster_columns = FALSE, show_row_dend = FALSE, 
                   show_column_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                   column_title = title, 
                   column_names_side = "top",column_names_rot = 360, column_names_centered = TRUE, column_names_gp = gpar(fontsize = 8),
                   #top_annotation = htop, 
                   right_annotation = hright,
                   row_title = "... cell markers",
                   heatmap_legend_param = list(at = c(-2, 0, 2), labels = c(-2, 0, 2), 
                                               color_bar = "continous", direction = "horizontal", 
                                               title_position = "lefttop", title = "log2 expression", 
                                               legend_width = unit(40, "mm"))
  ) # + hright
  
  draw(hmain, heatmap_legend_side = "bottom")
  
  #column_title_gp = gpar(fontsize = 6), column_title_rot = 90,
  
  
}
####################################################################################