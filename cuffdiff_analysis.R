# setwd(dir = "C:\\Users/sujay/Desktop/USC Assignments and Material/TRGN515/RNAseq Project/")
diff_exp_genes = read.table(file = "../cuffdiff results/Cuffdiff_wStats/gene_exp.diff", header = TRUE, sep = "\t")
sig_gene_data = subset(diff_exp_genes, p_value < 0.01)

genes_fpkm = read.table(file = "../cuffdiff results/Cuffdiff_NoStats/genes.fpkm_tracking", header = TRUE, sep = "\t")
genes_fpkm = genes_fpkm[genes_fpkm$gene_id %in% sig_gene_data$gene_id,]

library(dplyr)
genes_fpkm_matrix = select(genes_fpkm, 
                           gene_short_name, 
                           Ctrl1_FPKM, Ctrl2_FPKM, Ctrl3_FPKM, 
                           MidAD1_FPKM, MidAD2_FPKM, MidAD3_FPKM, 
                           LateAD1_FPKM, LateAD2_FPKM, LateAD3_FPKM)

genes_fpkm_matrix[,2:10] = log10(genes_fpkm_matrix[,2:10] + 1)

genes_fpkm_matrix_clean = genes_fpkm_matrix[,-1]
rownames(genes_fpkm_matrix_clean) = genes_fpkm_matrix[,1]

library(gplots)
library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(8, "Paired"))(25)
heatmap.2(x = as.matrix(genes_fpkm_matrix_clean), col = palette, margins = c(12, 9), trace = "none")

genes_fpkm_matrix_clean = as.data.frame(lapply(genes_fpkm_matrix_clean, as.double))
genes_fpkm_matrix_clean_t = as.data.frame(t(genes_fpkm_matrix_clean))

library(ggfortify)
pca_components = prcomp(genes_fpkm_matrix_clean_t)
genes_fpkm_matrix_clean_t$Samples = rownames(genes_fpkm_matrix_clean_t)
genes_fpkm_matrix_clean_t$Group = ifelse(grepl("Ctrl", genes_fpkm_matrix_clean_t$Samples), "CTRL",
                                         ifelse(grepl("Mid", genes_fpkm_matrix_clean_t$Samples), "Mid AD", "Late AD"))

autoplot(pca_components, data = genes_fpkm_matrix_clean_t, label = TRUE, label.size = 2, colour = "Group")

library(factoextra)
genes_fpkm_matrix_clean_t[is.na(genes_fpkm_matrix_clean_t)] = 2
dist.eucl = dist(genes_fpkm_matrix_clean_t, method = "euclidean")
hclust.ward.eucl = hclust(d = dist.eucl, method = "ward.D2")
plot(hclust.ward.eucl, main = "Dendrogram (Euclidean - Ward's)")

upregulated_genes = subset(sig_gene_data, log2.fold_change. > 0)
# write.table(upregulated_genes, "upregulated_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_genes = subset(sig_gene_data, log2.fold_change. < 0)
# write.table(downregulated_genes, "downregulated_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)