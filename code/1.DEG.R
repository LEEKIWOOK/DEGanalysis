# cran.packages <- c('knitr', 'mixOmics', 'xtable', 'ggplot2', 'vegan', 'cluster', 'gridExtra', 'pheatmap', 'ruv', 'lmerTest', 'bapred')
# install.packages(cran.packages)
# bioconductor.packages <- c('sva', 'limma', 'AgiMicroRna', 'variancePartition', 'pvca')
# if (!requireNamespace('BiocManager', quietly = TRUE))
#     install.packages('BiocManager')
# BiocManager::install(bioconductor.packages)

library(knitr)
library(xtable) # table
library(mixOmics)
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(gridExtra) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(AgiMicroRna) # RLE plot
library(cluster) # silhouette coefficient
library(variancePartition) # variance calculation
library(pvca) # PVCA
library(pheatmap) # heatmap
library(ruv) # RUVIII
library(lmerTest) # lmer
library(bapred) # FAbatch

#0. Data load
setwd("/home/kwlee/Projects_gflas/Team_CrisprCas/RNAexpression")
fpkm<-read.delim("./data/1.deg/fpkm.tsv", sep='\t')

#1.Log2 transform
log2col<-paste("log2",colnames(fpkm[-1]), sep="_")
log_fpkm = log2col = log(fpkm[-1]+1e-6, 2)
rownames(log_fpkm)<-fpkm$GeneSymbol
log_fpkm = log_fpkm[!is.infinite(rowSums(log_fpkm[-1])),]

#2.Batch correction


#2-1.define conditions, library methods, and replicates
sample_list = colnames(fpkm)[-1]
condition = c("N-RO","N-RO","N-RO","M-RBO","M-RBO","M-RBO","RSC","RSC","RSC","RSC","Y79","Y79","Y79","Y79","Y79")
library_methods = c(rep("inhouse",10),rep("publish_data",5))
replicates = c(1:3,1:3,1:4,1:5)

#2-2.calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(log_fpkm[,sample_list])

