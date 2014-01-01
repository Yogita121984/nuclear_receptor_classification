library (BSgenome)
library (ShortRead)
library (IRanges)
library (BayesPeak)
library (chipseq)
library (rtracklayer)
library (GenomicRanges)
library (GenomicFeatures)
library (BSgenome.Hsapiens.UCSC.hg18)
librarysize = read.table("TotalLibrarysize_allcellline.bed", header=F)
  
Gm_H3K4me2 = read.table ("CCAT_Significant_region_reads_h3k4me2_Gm", header=F)
H3K4me2.GRanges <- GRanges(seqnames=as.vector(Gm_H3K4me2$V1), IRanges(start=Gm_H3K4me2$V2, end=Gm_H3K4me2$V3), strand=Gm_H3K4me2$V6)
H3K4me2.cov   <- coverage(H3K4me2.GRanges)
H3K4me2.cov <-H3K4me2.cov/librarysize[1]

cl1.nexp.gm<-read.table("genesetcl1.notexp.10kb.gm.bed")
cl1.nexp.gm.sub1 <-cl1.nexp.gm.df[cl1.nexp.gm$V5 == "+", ]
cl1.nexp.gm.sub2 <-cl1.nexp.gm.df[cl1.nexp.gm$V5 == "-", ]
cl1.exp.gm<-read.table("geneset.exp.10kb.gm.bed")
cl1.exp.gm.sub1 <-cl1.exp.gm.df[cl1.exp.gm$V5 == "+", ]
cl1.exp.gm.sub2 <-cl1.exp.gm.df[cl1.exp.gm$V5 == "-", ]
cl2.nexp.gm<-read.table("genesetcl2.notexp.10kb.gm.bed")
cl2.nexp.gm.sub1 <-cl2.nexp.gm.df[cl2.nexp.gm$V5 == "+", ]
cl2.nexp.gm.sub2 <-cl2.nexp.gm.df[cl2.nexp.gm$V5 == "-", ]
cl2.exp.gm<-read.table("genesetcl2.exp.10kb.gm.bed")
cl2.exp.gm.sub1 <-cl2.exp.gm.df[cl2.exp.gm$V5 == "+", ]
cl2.exp.gm.sub2 <-cl2.exp.gm.df[cl2.exp.gm$V5 == "-", ]

cluster1_Nexp.sub1<-data.frame()
for (i in c(1:nrow(cl1_nexp.sub1))) {
chromosome<-as.character(cl1_nexp.sub1$V1[i])
start <-cl1_nexp.sub1$V2[i]
end <-cl1_nexp.sub1$V3[i]
midpoint=(start+end)/2
view <-window (H3K4me2.cov[[chromosome]], start=midpoint-10000, end=midpoint+10000)
a1 = as.vector(view)
b1 <-t(a1)
cluster1_Nexp.sub1 <-rbind(cluster1_Nexp.sub1, b1)
}
cluster1_Nexp.sub2<-data.frame()
for (i in c(1:nrow(cl1_nexp.sub2))) {
chromosome<-as.character(cl1_nexp.sub2$V1[i])
start <-cl1_nexp.sub2$V2[i]
end <-cl1_nexp.sub2$V3[i]
midpoint=(start+end)/2
view <-window (H3K4me2.cov[[chromosome]], start=midpoint-10000, end=midpoint+10000)
a1 = as.vector(view)
b1 <-t(a1)
cluster1_Nexp.sub2 <-rbind(cluster1_Nexp.sub2, b1)
}

cluster1_Nexp.sub2<-rev(cluster1_Nexp.sub2)
cluster1_Nexp.total<-rbind(cluster1_Nexp.sub2, cluster1_Nexp.sub1)

cluster1_exp.sub1<-data.frame()
for (i in c(1:nrow(cl1_exp.sub1))) {
chromosome<-as.character(cl1_exp.sub1$V1[i])
start <-cl1_exp.sub1$V2[i]
end <-cl1_exp.sub1$V3[i]
midpoint=(start+end)/2
view <-window (H3K4me2.cov[[chromosome]], start=midpoint-10000, end=midpoint+10000)
a1 = as.vector(view)
b1 <-t(a1)
cluster1_exp.sub1 <-rbind(cluster1_exp.sub1, b1)
}
cluster1_exp.sub2<-data.frame()
for (i in c(1:nrow(cl1_exp.sub2))) {
chromosome<-as.character(cl1_exp.sub2$V1[i])
start <-cl1_exp.sub2$V2[i]
end <-cl1_exp.sub2$V3[i]
midpoint=(start+end)/2
view <-window (H3K4me2.cov[[chromosome]], start=midpoint-10000, end=midpoint+10000)
a1 = as.vector(view)
b1 <-t(a1)
cluster1_exp.sub2 <-rbind(cluster1_exp.sub2, b1)
}

cluster1_exp.sub2<-rev(cluster1_exp.sub2)
cluster1_exp.total<-rbind(cluster1_exp.sub2, cluster1_exp.sub2)
#Plotting for average coverage for all cell lines.
plot_colors<-rainbow(5)
names <-c("Gm12878", "Hepg2", "Huvec", "H1hesc", "k562")
pdf (file = "H3k4me2_coverage.pdf")
> attach(mtcars)
> par(mfrow=c(1,2))
> plot( colSums(cluster1_exp.total),type="l", axes=FALSE,col=plot_colors[1])
> lines(colSums(Hepg2_H3k4me2.sub1),col=plot_colors[2])
> lines(colSums(Huvec_H3k4me2.sub1),col=plot_colors[3])
> lines(colSums(H1hesc_H3k4me2.sub1),col=plot_colors[4])
> lines(colSums(k562_H3k4me2.sub1),col=plot_colors[5])
> axis (1,at=seq(from=0, to=10000, by=1000), lab=seq(from=-5000,to=5000,by=1000))
> legend("topright", names, cex=0.8, col=plot_colors,
+ pch=21:23, lty=1:3);
> title(main="H3k4me2_Cluster1_exp", col.main="blue", font.main=4)
> plot( colSum(cluster1_Nexp.total),type="l", axes=FALSE,col=plot_colors[1], xlab="Position", ylab="Coverage")
> lines(colSums(Hepg2_H3k4me2.sub2),col=plot_colors[2])
> lines(colSums(Huvec_H3k4me2.sub2),col=plot_colors[3])
> lines(colSums(H1hesc_H3k4me2.sub2),col=plot_colors[4])
> lines(colSums(k562_H3k4me2.sub2),col=plot_colors[5])
> axis (1,at=seq(from=0, to=10000, by=1000), lab=seq(from=-5000,to=5000,by=1000))
> legend("topright", names, cex=0.8, col=plot_colors,pch=21:23, lty=1:3);
> title(main="H3k4me2_Cluster1_Nexp", col.main="blue", font.main=4)
> dev.off()






