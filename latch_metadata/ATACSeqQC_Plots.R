#!/bin/Rscript
library("ATACseqQC")
library("Rsamtools")
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BSgenome.Hsapiens.UCSC.hg38")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("ChIPpeakAnno")
library("GenomeInfoDb")
library("MotifDb")
library("Rsamtools")
library('Cairo')

args <- commandArgs(trailingOnly=TRUE)

merged_bam_file <- args[1]
outPath <- args[2]


if (dir.exists(outPath))
{
  unlink(outPath, recursive = TRUE, force = TRUE)
}
dir.create(outPath)

fragSize <- fragSizeDist(merged_bam_file, "BamFile")
x <- fragSize$BamFile
filt <- x[1:1010]
filt[is.na(filt)] <- 0
y <- filt / sum(filt)
y <- as.numeric(y)
write.table(y, paste(outPath, "Frag_Sizes.txt", sep=""))

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))
bamTop100 <- Rsamtools::scanBam(Rsamtools::BamFile(merged_bam_file,
                                                   yieldSize = 100),
                     param = Rsamtools::ScanBamParam(tag = possibleTag))[[1]]$tag
seqlev <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
            "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
            "chr18","chr19","chr20","chr21","chr22","chrX","chrY")

which <- as(GenomeInfoDb::seqinfo(Hsapiens)[seqlev], "GRanges")
tags <- names(bamTop100)[lengths(bamTop100)>0]

gal <- readBamFile(merged_bam_file, tag=tags, which=which,
                   asMates=TRUE, bigFile=TRUE)
shiftedBamFile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamFile)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs <- txs[seqnames(txs) %in% seqlev]
genome <- Hsapiens
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)

bamFiles <- file.path(outPath,
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
(librarySize <- estLibSize(bamFiles))

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")],
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

out <- featureAlignedDistribution(sigs, reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt = .5, n.tile = NTILE, type = "l",
                                  ylab = "Averaged coverage")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
write.table(out, paste(outPath,"featurealignment_coverage.txt", sep = ""), sep = "\t")

s <- estimateLibComplexity(readsDupFreq(merged_bam_file))
while (!is.null(dev.list()))  dev.off()
write.table(s, paste(outPath, "Saturation_Plots.txt", sep = ""), sep = "\t")

#CTCF <- query(MotifDb, c("CTCF"))
#CTCF <- as.list(CTCF)
#while (!is.null(dev.list()))  dev.off()
#CTCF_Forward <- list()
#CTCF_Reverse <- list()
#for (j in 1: length(CTCF)){
#    sigs <- factorFootprints(shiftedBamFile, pfm=CTCF[[j]],
#                             genome=genome, min.score="90%",
#                             seqlev=c('chrX','chrY'), upstream=100,
#                             downstream=100)
#    for_hits = nrow(sigs$signal$`+`)
#    rev_hits = nrow(sigs$signal$`-`)
#    signal_avg_forward = colSums(sigs$signal$`+`)
#    signal_avg_reverse = colSums(sigs$signal$`-`)
#    while (!is.null(dev.list())) dev.off()
#    for (i in 1: (length(seqlev)-2))
#    {
#        print(paste0(j,"Chromosome--->",seqlev[i]))
#        sigs <- factorFootprints(shiftedBamFile, pfm=CTCF[[j]],
#                                 genome=genome, min.score="90%",
#                                 seqlev=c(seqlev[i]), upstream=100,
#                                 downstream=100)
#        signal_avg_forward = signal_avg_forward + colSums(sigs$signal$`+`)
#        signal_avg_reverse = signal_avg_reverse + colSums(sigs$signal$`-`)
#        for_hits = for_hits + nrow(sigs$signal$`+`)
#        rev_hits = rev_hits + nrow(sigs$signal$`-`)
#        while (!is.null(dev.list())) dev.off()
#    }
#    CTCF_Forward[[j]] <- (signal_avg_forward/for_hits)
#    CTCF_Reverse[[j]] <- (signal_avg_reverse/rev_hits)
#}
