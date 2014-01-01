library(gplots)
library(RColorBrewer)
data=read.table ("HMMstatecount_allcluster.bed", header=T)
row.names(data) <- data$gene
data <- data[,2:14]
data <- data.matrix(data)
data = data +1
data <- log2(data)
row.dist <- as.dist(1 - cor(t(data)))
col.dist <- as.dist(1 - cor(data))
hmcol <- colorRampPalette(brewer.pal (10, "RdBu"))(256)
pdf (file="HMM_fig7.pdf")
heatmap.2(hmm_exp,col=hmcol,Colv=as.dendrogram(hclust(col.dist, method = "complete")),RowV = as.dendrogram(hclust(row.dist, method = "complete")), key = TRUE,keysize = 1.2,  xlab = "HMM States", ylab = "Genes", trace="none", cexRow=0.6, cexCol=1.0)
dev.off()
