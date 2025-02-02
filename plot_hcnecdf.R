hcne_cl1 <- read.table("HCNE_cl1_1kbwindow_allSP_allPC.2mb", header=T, sep="\t")
hcne_cl2 <- read.table("HCNE_cl2_1kbwindow_allSP_allPC.2mb", header=T, sep="\t")
cl1.mm9 <- ecdf(hcne_cl1$mouse)
cl1.galgal <- ecdf(hcne_cl1$chicken)
cl1.xentro <- ecdf(hcne_cl1$frog)
cl1.fr <- ecdf(hcne_cl1$fugu)
cl2.danrer <- ecdf(hcne_cl1$zebrafish)
cl2.mm9 <- ecdf(hcne_cl2$mouse)
cl2.galgal <- ecdf(hcne_cl2$chicken)
cl2.xentro <- ecdf(hcne_cl2$frog)
cl2.fr <- ecdf(hcne_cl2$fugu)
cl2.danrer <- ecdf(hcne_cl2$zebrafish)
names <- c("Cluster1", "Cluster2")
plot_colors <- c("red","blue")
pdf (file = "Figure3_cdf.pdf")
par(mfrow=c(2,3))
plot(cl1.mm9, ylab="Fraction", xlab = "HCNE density (1 kb window)", main = "Human : Mouse", verticals = FALSE, col.01line = "gray70", pch = 19, col = plot_colors[1])
lines (cl2.mm9, col = plot_colors[2])
box()
legend(2, 0.2, names, cex=0.9, col=plot_colors, pch=21:21, lty=1:1)
plot(cl1.galgal, ylab="Fraction", xlab = "HCNE density (1 kb window)", main = "Human : Chicken", verticals = FALSE, col.01line = "gray70", pch = 19, col = plot_colors[1])
lines (cl2.galgal, col = plot_colors[2])
box()
legend(2, 0.2, names, cex=0.9, col=plot_colors, pch=21:21, lty=1:1)
plot(cl1.xentro, ylab="Fraction", xlab = "HCNE density (1 kb window)", main = "Human : Frog", verticals = FALSE, col.01line = "gray70", pch = 19, col = plot_colors[1])
lines (cl2.xentro, col = plot_colors[2])
box()
legend(2, 0.2, names, cex=0.9, col=plot_colors, pch=21:21, lty=1:1)
plot(cl1.fr, ylab="Fraction", xlab = "HCNE density (1 kb window)", main = "Human : Fugu", verticals = FALSE, col.01line = "gray70", pch = 19, col = plot_colors[1]) 
lines (cl2.fr, col = plot_colors[2])
box()
legend(1, 0.2, names, cex=0.9, col=plot_colors, pch=21:21, lty=1:1)
plot(cl1.danrer, ylab="Fraction", xlab = "HCNE density (1 kb window)", main = "Human : Zebrafish", verticals = FALSE, col.01line = "gray70", pch = 19, col = plot_colors[1])
lines (cl2.danrer, col = plot_colors[2])
box()
legend(1, 0.2, names, cex=0.9, col=plot_colors, pch=21:21, lty=1:1)
dev.off()
