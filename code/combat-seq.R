setwd("/home/kwlee/Projects_gflas/Team_CrisprCas/RNAexpression/data/temp")

library("sva") #Note this exercise requires sva (>= v3.36.0) which is only available for R (>= v4.x)
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")

uncorrected_data = read.table("GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv", header=TRUE, sep="\t", as.is=c(1,2))
#simplify the names of the data columns
# (A = Universal Human Reference RNA and B = Human Brain Reference RNA)
# RNA = polyA enrichment and RIBO = ribosomal RNA depletion
# 1, 2, 3, 4 are replicates
names(uncorrected_data) = c("Gene", "Chr", "UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4", "HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4", 
                            "UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4", "HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")
sample_names = names(uncorrected_data)[3:length(names(uncorrected_data))]

#review data structure
head(uncorrected_data)
dim(uncorrected_data)

#define conditions, library methods, and replicates
conditions = c("UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR", "UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR")
library_methods = c("Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly")
replicates = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)

#calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(uncorrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

#assign labels to the data frame
pca_uncorrected[,"condition"] = conditions
pca_uncorrected[,"library_method"] = library_methods
pca_uncorrected[,"replicate"] = replicates

#plot the PCA
#create a classic 2-dimension PCA plot (first two principal components) with conditions and library methods indicated
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (uncorrected data)", color="Condition", shape="Library Method")
p1 = p1 + scale_colour_manual(values = cols)



#perform the batch correction

#first we need to transform the format of our groups and batches from names (e.g. "UHR", "HBR", etc.) to numbers (e.g. 1, 2, etc.)
#in the command below "sapply" is used to apply the "switch" command to each element and convert names to numbers as we define
groups = sapply(as.character(conditions), switch, "UHR" = 1, "HBR" = 2, USE.NAMES = F)
batches = sapply(as.character(library_methods), switch, "Ribo" = 1, "Poly" = 2, USE.NAMES = F)

#now run ComBat_seq
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = batches, group = groups)

#join the gene and chromosome names onto the now corrected counts from ComBat_seq
corrected_data = cbind(uncorrected_data[,c("Gene","Chr")], corrected_data)

#compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
dim(corrected_data)

#visually compare values of corrected and uncorrected data sets
head(uncorrected_data)
head(corrected_data)


#calculate principal components for the uncorrected data
pca_corrected_obj = prcomp(corrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

#assign labels to the data frame
pca_corrected[,"condition"] = conditions
pca_corrected[,"library_method"] = library_methods
pca_corrected[,"replicate"] = replicates

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
p2 = p2 + scale_colour_manual(values = cols)

pdf(file="Uncorrected-vs-BatchCorrected-PCA.pdf")
grid.arrange(p1, p2, nrow = 2)
dev.off()


#perform differential expression analysis on the uncorrected data and batch corrected data sets

#first define the sets of samples to be compared to each other
uhr_ribo_samples = c("UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4")
uhr_poly_samples = c("UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4")
hbr_ribo_samples = c("HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4")
hbr_poly_samples = c("HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")
uhr_samples = c(uhr_ribo_samples, uhr_poly_samples)
hbr_samples = c(hbr_ribo_samples, hbr_poly_samples)

#create a function that will run edgeR (DE analysis) for a particular pair of sample sets
run_edgeR = function(data, group_a_name, group_a_samples, group_b_samples, group_b_name){
  #create a list of all samples for this current comparison
  samples_for_comparison = c(group_a_samples, group_b_samples)
  
  #define the class factor for this pair of sample sets
  class = factor(c(rep(group_a_name,length(group_a_samples)), rep(group_b_name,length(group_b_samples))))
  
  #create a simplified data matrix for only these samples
  rawdata = data[,samples_for_comparison]
  
  #store gene names for later
  genes = rownames(data)
  gene_names = data[,"Gene"]
  
  #make DGElist object
  y = DGEList(counts=rawdata, genes=genes, group=class)
  
  #perform TMM normalization
  y <- calcNormFactors(y)
  
  #estimate dispersion
  y <- estimateCommonDisp(y, verbose=TRUE)
  y <- estimateTagwiseDisp(y)
  
  #differential expression test
  et <- exactTest(y)
  
  #print number of up/down significant genes at FDR = 0.05 significance level and store the DE status in a new variable (de)
  summary(de <- decideTestsDGE(et, p=.05))
  summary(de <- decideTestsDGE(et, adjust.method="fdr", p=.05))
  
  #create a matrix of significant DE genes
  mat <- cbind(
    genes, gene_names,
    sprintf('%0.3f', log10(et$table$PValue)),
    sprintf('%0.3f', et$table$logFC)
  )[as.logical(de),]
  colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")
  
  #order by log fold change
  o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
  mat <- mat[o,]
  
  #fix the issue where corrected p-values that are 0 become -Inf upon log10 conversion
  x = as.numeric(mat[,"Log10_Pvalue"])
  lowest_pvalue = min(x[which(!x == -Inf)])
  i = which(mat[,"Log10_Pvalue"] == -Inf)
  mat[i,"Log10_Pvalue"] = lowest_pvalue
  
  return(mat)
}

#run the five comparisons through edgeR using the *uncorrected data*
uhr_ribo_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

#run the same five comparisons through edgeR using the *batch corrected data*
uhr_ribo_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

#how much of a difference does batch correction make when doing the comparison of all UHR vs all HBR samples?
dim(uhr_vs_hbr_uncorrected)
dim(uhr_vs_hbr_corrected)

#create upset plots to summarize the overlap between the comparisons performed above

#first create upset plot from the uncorrected data
pdf(file="Uncorrected-UpSet.pdf")
listInput = list("4 UHR Ribo vs 4 HBR Ribo" = uhr_ribo_vs_hbr_ribo_uncorrected[,"Gene"], 
                 "4 UHR Poly vs 4HBR Poly" = uhr_poly_vs_hbr_poly_uncorrected[,"Gene"],
                 "4 UHR Ribo vs 4 HBR Poly" = uhr_ribo_vs_hbr_poly_uncorrected[,"Gene"],
                 "4 UHR Poly vs 4 HBR Ribo" = uhr_poly_vs_hbr_ribo_uncorrected[,"Gene"],
                 "8 UHR vs 8 HBR" = uhr_vs_hbr_uncorrected[,"Gene"])
upset(fromList(listInput), order.by = "freq", number.angles=45, point.size=3)
dev.off()

#now create an upset plot from the batch corrected data
pdf(file="BatchCorrected-UpSet.pdf")
listInput = list("4 UHR Ribo vs 4 HBR Ribo" = uhr_ribo_vs_hbr_ribo_corrected[,"Gene"], 
                 "4 UHR Poly vs 4 HBR Poly" = uhr_poly_vs_hbr_poly_corrected[,"Gene"],
                 "4 UHR Ribo vs 4 HBR Poly" = uhr_ribo_vs_hbr_poly_corrected[,"Gene"],
                 "4 UHR Poly vs 4 HBR Ribo" = uhr_poly_vs_hbr_ribo_corrected[,"Gene"],
                 "8 UHR vs 8 HBR" = uhr_vs_hbr_corrected[,"Gene"])
upset(fromList(listInput), order.by = "freq", number.angles=45, point.size=3)
dev.off()

#write out the final set of DE genes where all UHR and HBR samples were compared using the corrected data
write.table(uhr_vs_hbr_corrected, file="DE_genes_uhr_vs_hbr_corrected.tsv", quote=FALSE, row.names=FALSE, sep="\t")

#To exit R type the following
quit(save="no")


