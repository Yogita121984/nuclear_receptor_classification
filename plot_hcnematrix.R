library(gplots)
library(RColorBrewer)
rawdata<-read.table("HCNE_4GRID.allSP_allPC_2mb.bed", header=T)
data_matrix<-rawdata[,2:6]
row.names(data_play)<-rawdata$Gene
data_matrix=data_matrix+0.0001
data_matrix_log<-log2(data_matrix)
data_dist_m<-dist(data_matrix_log, method = "euclidean")
data_dist.matrix_m<-as.matrix(data_dist_m)
hmcol <- colorRampPalette(brewer.pal (10, "RdBu"))(256)
colnames(data_dist.matrix_m)<-rawdata$Gene
row.names(data_dist.matrix_m)<-rawdata$Gene
pdf (file="Fig2_hcne.pdf")
heatmap.2(data_dist.matrix_m, trace="none", dendrogram="row", col=hmcol, key=TRUE, keysize=1.0, xlab="Genes", ylab="Genes", cexRow=0.6, cexCol=0.6, density.info="none")
dev.off()
