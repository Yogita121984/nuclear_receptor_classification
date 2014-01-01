library (BSgenome)
library(biomaRt)
library (IRanges)
library (rtracklayer)
library (GenomicRanges)
library (GenomicFeatures)
library (BSgenome.Hsapiens.UCSC.hg18)

slideNo <- "H3k4me1"
chipData <- read.table( file="H3k4me1_h1hesc.bed", header=T, sep="\t")
chip_Range <- GRanges(seqnames = Rle( chipData$chromosome), c(1:length(chipData$start)), ranges = IRanges(start=chipData$start, end = chipData$end), strand = Rle(strand( NA )),score=chipData$chipSignal  )
chrs = paste("chr", c("1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X"), sep = "")
seqlengths(chip_Range) <-seqlengths(Hsapiens)[chrs]

targetGenes <- read.table( file="target_grb.gcord", header=T)
target_ranges <- GRanges(seqnames = Rle( targetGenes$chromosome), c(1:length(targetGenes$start)), ranges = IRanges(start=targetGenes$start, end =targetGenes$end), strand = Rle(strand( targetGenes$strand )),genes=targetGenes$gene  )
chrs = paste("chr", c("1","14","15","2","3","5","6","9"), sep = "")
seqlengths(target_ranges) <- seqlengths(Hsapiens)[chrs]

nontargetGenes <- read.table( file="nonTar_grb.gcord", header=T)
nontarget_ranges <- GRanges(seqnames = Rle( nontargetGenes$chromosome), c(1:length(nontargetGenes$start)), ranges = IRanges(start=nontargetGenes$start, end =nontargetGenes$end), strand = Rle(strand( nontargetGenes$strand )),genes=nontargetGenes$gene  )
chrs = paste("chr", c("1","11","12","14","15","17","19","20","22","3","4","5","6","8","9","X"), sep = "")
seqlengths(nontarget_ranges) <- seqlengths(Hsapiens)[chrs]


cPgData <- read.table( file="human_chrX_CpG.txt", header=T)
cpg_ranges <- GRanges(seqnames = Rle( cPgData$chromosome), c(1:length(cPgData$start)), ranges = IRanges(start=cPgData$start, end =cPgData$end), strand = Rle(strand( NA )),length=cPgData$length  )
#seqlengths(cpg_ranges) <- chromosomeLength
seqlengths(cpg_ranges) <- seqlengths(Hsapiens)['chrX']
ensembl <-useMart("ensembl_mart_52", dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=TRUE)
gene.ids <- unique(unlist(lapply(as.list(c("X")),
function(this.chr)
getBM(attributes="ensembl_gene_id", filters="chromosome_name",
values=this.chr, mart=ensembl)[,1]), use.names=FALSE))
sel.attributes=c("ensembl_gene_id", "chromosome_name",
"strand", "start_position","end_position", "description")

hggenes <- getBM(attributes=sel.attributes, filters="ensembl_gene_id",value=gene.ids, mart=ensembl)

hggenes$name <- hggenes$"ensembl_gene_id"
hggenes$gene <- hggenes$"ensembl_gene_id"
hggenes$chr <- hggenes$chromosome_name
hggenes$symbol <- hggenes$"mgi_symbol"
hggenes$start <- hggenes$"start_position"
hggenes$end <- hggenes$"end_position"
hggenes$feature <- rep("gene",nrow(hggenes))

ensembl_ranges <-GRanges(seqnames = Rle( hggenes$chromosome_name), c(1:length( hggenes$ensembl_gene_id )), ranges = IRanges(start=hggenes$start_position, end =hggenes$end_position),strand = Rle(strand( hggenes$strand )), hggenes$ensembl_gene_id   )
seqlengths(ensembl_ranges) <- seqlengths(Hsapiens)['chrX']

cpg_ensembl_genes <- ensembl_ranges[as.matrix(findOverlaps( flank(ensembl_ranges, 1000, both=T ), cpg_ranges))[,1], ]
non_cpg_ensembl_genes <- ensembl_ranges[-as.matrix(findOverlaps( flank(ensembl_ranges, 1000, both=T ), cpg_ranges))[,1], ]

#function_to  get upstreem and downstreem coordinates of TSS
getUpstreamAndDownstreamCoordinatesOfTSS <- function( geneCoordinates, flankingDistance )
{
   #get tss sites
   #in GRanges object on + starnd TSS - start coordinate
   #in GRanges object on - starnd TSS - END coordinate
   geneCoordinates$tss <- ifelse( geneCoordinates$strand == "+", (geneCoordinates$start ), geneCoordinates$end )
   geneCoordinates$start <- geneCoordinates$tss - flankingDistance
   geneCoordinates$end <- geneCoordinates$tss + flankingDistance
   #re arrange data frame,seqnames, start, end, genes, strand, tss
   geneCoordinates <- geneCoordinates[ , c(1,2,3,7,5, 8)]
   names(geneCoordinates) <- c( "chromosome", "start", "end", "gene", "strand", "tss")
   return(geneCoordinates)

}
#Function to conververt data frame to GRange object
#dataframe with 6 columns
dataframe2GrangesObject <- function( dataFrameObject)
{
   grangesObject <- GRanges(seqnames = Rle( dataFrameObject$chromosome), c(1:length(dataFrameObject$start)),
   ranges =IRanges(start=dataFrameObject$start, end =dataFrameObject$end),
   strand = Rle(strand( dataFrameObject$strand )),genes=dataFrameObject$gene, tss=dataFrameObject$tss  )
   return(grangesObject)

}


findOverLapsAndExtractElements <- function( genes, chipSignal, flankingDistance )
{
  overLaps <- findOverlaps(genes, chipSignal)
  tss_and_chipsignals <- as.data.frame(genes[ as.matrix(overLaps)[ , 1],])[,c(1,8,7)]
  tss_and_chipsignals$signal <- as.data.frame(chipSignal[ as.matrix(overLaps)[ , 2],])[ , 7]
 tss_and_chipsignals$chipRelativePosition <- (start(genes[as.matrix(overLaps)[ , 1], ]) + flankingDistance )- start(chipSignal[as.matrix(overLaps)[ ,2 ], ])
  tss_and_chipsignals <- tss_and_chipsignals[ , c(1,2,5,4,3)]
  return(tss_and_chipsignals)
}

#Function to sample genes from set of genes and give total number of signals

sampleGenesAndGiveQunatityOfSignals <- function( geneSet, totalNoOfGenes, noOfGenesForSampling, distances,chip_Range )
{
   #geneSet, chip_Range - GRanges objects
   #totalNoOfGenes, noOfGenesForSampling, distances - numaric
   up_and_down_tss <- getUpstreamAndDownstreamCoordinatesOfTSS(as.data.frame( geneSet
   [ sample(1: totalNoOfGenes, noOfGenesForSampling), ]), distances)
   up_and_down_tss <- dataframe2GrangesObject(up_and_down_tss)
   chipEnrichmentData <- findOverLapsAndExtractElements( up_and_down_tss ,chip_Range, distances)
   noOfSignals <- length(chipEnrichmentData$signal)
   return(noOfSignals)
}


#put the flanking distances in the loop say 1KB, 1MB, 2 MB
distances <- as.integer( c( 10000, 1000000, 2000000))

for( i in 1:length( distances ))
{
   #initiate empty vectors to hold values
   ranTargetSigNumbers <- numeric()
   ranCpGSigNumbers <- numeric()
   ranNonCpGSigNumbers <- numeric()
   rannonTargetSigNumbers <- numeric()
   for ( x in 1:1000 )
   {
      up_and_down_tss <- getUpstreamAndDownstreamCoordinatesOfTSS(as.data.frame( target_ranges[ sample(1:25, 24), ]),
      distances[i])
      up_and_down_tss <- dataframe2GrangesObject(up_and_down_tss)
      chipEnrichmentData <- findOverLapsAndExtractElements( up_and_down_tss ,chip_Range, distances[i])
      #how many unique target genes with chip signls?
      noOfGenesForSampling <- length(unique(as.data.frame(chipEnrichmentData)$genes))
      noOfSignals_target <- length(chipEnrichmentData$signal)
      ranTargetSigNumbers[x] <- noOfSignals_target
      #Get signals for random CpG genes
      noOfSignals_RanCpG <- sampleGenesAndGiveQunatityOfSignals( cpg_ensembl_genes, 402, noOfGenesForSampling,
      distances[i], chip_Range )
      ranCpGSigNumbers[x] <- noOfSignals_RanCpG
   #Get signals for random random NonCpG genes
      noOfSignals_RanNonCpG <- sampleGenesAndGiveQunatityOfSignals( non_cpg_ensembl_genes, 1652, noOfGenesForSampling,
      distances[i], chip_Range )
      ranNonCpGSigNumbers[x] <- noOfSignals_RanNonCpG

     #Get signals for random random NonTarget genes

     noOfSignals_nontarget <- sampleGenesAndGiveQunatityOfSignals( nontarget_ranges, 23, noOfGenesForSampling,
     distances[i], chip_Range )
     rannonTargetSigNumbers[x] <-noOfSignals_nontarget


}
   all_signals <- data.frame( ranTargetSigNumbers, ranCpGSigNumbers, ranNonCpGSigNumbers, rannonTargetSigNumbers)
   outfileName <- paste(slideNo,"all_chipSidnals", distances[i], sep="_")
   write.table( all_signals, file=outfileName, col.names=T, row.names=F, quote=F, sep="\t")

}



