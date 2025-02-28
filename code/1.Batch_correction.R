# cran.packages <- c('knitr', 'xtable', 'ggplot2', 'vegan', 'cluster', 'gridExtra', 'pheatmap', 'ruv', 'lmerTest', 'bapred')
# install.packages(cran.packages)
# 
# bioconductor.packages <- c('sva', 'limma', 'AgiMicroRna', 'variancePartition', 'pvca', 'mixOmics')
# if (!requireNamespace('BiocManager', quietly = TRUE))
#     install.packages('BiocManager')
# BiocManager::install(bioconductor.packages)
# 
# install.packages("remotes")
# remotes::install_github("ajabadi/mixOmics2")

####################################################################################

library(knitr)
library(xtable) # table
library(mixOmics)
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(gridExtra) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)
library(pheatmap) # heatmap

#0. Data load
setwd("/home/kwlee/Projects_gflas/Team_CrisprCas/RNAexpression")
source(file = './code/Functions.R')

fpkm<-read.delim("./data/fpkm.tsv", sep='\t')
rownames(fpkm)<-fpkm$GeneSymbol
fpkm<-as.matrix(fpkm[-1])

####################################################################################
#1.Log transform
log_fpkm = logratio.transfo(fpkm, logratio = "CLR", offset=1e-6) #CLR (centered log ratio transformation)
class(log_fpkm) = 'matrix'

####################################################################################
#2.Batch effect detection
replicates = c(rep("N-RO",3),rep("M-RBO",3),rep("RSC",4),rep("Y79",5))
batch<-factor(replicates)
names(batch) = colnames(fpkm)

experiment = rep("inhouse",15)
library<-factor(experiment)
names(library) = colnames(fpkm)

fpkm.pca.before <- pca(t(log_fpkm), ncomp = ncol(log_fpkm), scale=TRUE, multilevel = replicates, center = TRUE)


Scatter_Density(mat = fpkm.pca.before$variates$X, B = batch, T = library, E = fpkm.pca.before$prop_expl_var$X,
                xlim = c(-150,150), ylim = c(-150,150),
                batch.legend.title = 'Cell type (batch)',
                trt.legend.title = 'Experiment (trt)',
                title = 'Before batch effect correction')

####################################################################################
#3.Batch effect correction

d = data.frame(batch = as.integer(batch))
rownames(d) = colnames(fpkm)

#fpkm.mod <- model.matrix( ~ d) # full model
fpkm.mod <- model.matrix(~batch, data = d) # full model
fpkm.mod0 <- model.matrix( ~ 1, data = d) # null model

##################################################################################################################
#Adjusting for surrogate variables using the limma package

#fpkm.sva.n <- num.sv(dat = log_fpkm, mod = fpkm.mod, method = "leek") #variables in rows and samples in columns
fpkm.sva.n <- num.sv(dat = log_fpkm, mod = fpkm.mod, method = "be")

#To estimate the surrogate variables with both full and null models
fpkm.sva <- sva(dat = log_fpkm, mod = fpkm.mod, mod0 = fpkm.mod0, n.sv = fpkm.sva.n)
#Number of significant surrogate variables is:  13

#Estimated surrogate variables in both the null and full models
fpkm.mod.bat <- cbind(fpkm.mod, fpkm.sva$sv)
fpkm.mod0.bat <- cbind(fpkm.mod0, fpkm.sva$sv)

#Calculate parametric F-test P-values and Q-values (adjusted P-values) 
fpkm.sva.trt_p <- f.pvalue(log_fpkm, fpkm.mod.bat, fpkm.mod0.bat)
fpkm.sva.trt_adjp <- p.adjust(fpkm.sva.trt_p, method='fdr')
##################################################################################################################
#removeBatchEffect is a function implemented in the LIMMA package that fits a linear model 
# + for each variable given a series of conditions as explanatory variables, 
# + including the batch effect and treatment effect.
fpkm.limma <- removeBatchEffect(log_fpkm, batch = batch, design = fpkm.mod)

##################################################################################################################
#Applying the ComBat function to adjust for known batches
combat_res = ComBat(dat = log_fpkm, batch = batch, mod = fpkm.mod0, par.prior = TRUE, prior.plots = TRUE)

####################################################################################
#4. Result of batch effect correction
fpkm.pca.combat = pca(t(combat_res), ncomp = ncol(log_fpkm), scale=TRUE, multilevel = replicates, center = TRUE)
fpkm.pca.limma = pca(t(fpkm.limma), ncomp = ncol(log_fpkm), scale=TRUE, multilevel = replicates, center = TRUE)

fpkm.pca.plot.before <- Scatter_Density(mat = fpkm.pca.before$variates$X, B = batch, T = library, E = fpkm.pca.before$prop_expl_var$X, xlim = c(-150,150), ylim = c(-150,150), batch.legend.title = 'Cell type (batch)', trt.legend.title = 'Experiment (trt)', title = 'Before batch effect correction')
fpkm.pca.plot.limma <- Scatter_Density(mat = fpkm.pca.limma$variates$X, B = batch, T = library, E = fpkm.pca.limma$prop_expl_var$X, xlim = c(-100,150), ylim = c(-100,150), batch.legend.title = 'Cell type (batch)', trt.legend.title = 'Experiment (trt)', title = 'Batch correction with rBE')
fpkm.pca.plot.combat <- Scatter_Density(mat = fpkm.pca.combat$variates$X, B = batch, T = library, E = fpkm.pca.combat$prop_expl_var$X, xlim = c(-100,150), ylim = c(-100,100), batch.legend.title = 'Cell type (batch)', trt.legend.title = 'Experiment (trt)', title = 'Batch correction with Combat')


grid.arrange(fpkm.pca.plot.before, fpkm.pca.plot.limma, fpkm.pca.plot.combat, ncol=2)

####################################################################################
#5. boxplot
gene_name <- "ASCC3"
gene_idx = which(rownames(fpkm) %in% gene_name)

fpkm.before.df = data.frame(value = log_fpkm[gene_idx, ], batch = batch)
fpkm.boxplot.before = box_plot_fun(data = fpkm.before.df, x = fpkm.before.df$batch, y = fpkm.before.df$value, title = paste(gene_name, "- before", sep=' '), batch.legend.title = 'Cell type (batch)')

fpkm.limma.df = data.frame(value = fpkm.limma[gene_idx, ], batch = batch)
fpkm.boxplot.limma = box_plot_fun(data = fpkm.limma.df, x = fpkm.limma.df$batch, y = fpkm.limma.df$value, title = paste(gene_name, "- limma", sep=' '), batch.legend.title = 'Cell type (batch)')

grid.arrange(fpkm.boxplot.before, fpkm.boxplot.limma, ncol=2)