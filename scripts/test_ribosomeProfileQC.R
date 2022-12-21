# ribosomeProfilQC vigenette code
#
# Chap 1
#
library(BiocManager)
BiocManager::install(c("ribosomeProfilingQC", 
                       "AnnotationDbi", "Rsamtools",
                      # "BSgenome.Drerio.UCSC.danRer10",
                       "TxDb.Drerio.UCSC.danRer10.refGene",
                       "motifStack"), force=TRUE)

## load libraries             
library(ribosomeProfilingQC)
library(AnnotationDbi)
library(Rsamtools)
library(BSgenome.Drerio.UCSC.danRer10)
library(GenomicAlignments)
library(GenomicFeatures)
## set genome, Drerio is a shortname for BSgenome.Drerio.UCSC.danRer10
genome <- Drerio

## which is corresponding to BSgenome.Drerio.UCSC.danRer10
library(TxDb.Drerio.UCSC.danRer10.refGene)
txdb <- TxDb.Drerio.UCSC.danRer10.refGene ## give it a short name
CDS <- prepareCDS(txdb)

# get bam
library(Rsamtools)
## input the bamFile from the ribosomeProfilingQC package 
bamfilename <- system.file("extdata", "RPF.WT.1.bam",
                           package="ribosomeProfilingQC")
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)       

estimatePsite(bamfile, CDS, genome)                           
#?ribosomeProfilingQC:::fixSeqlevelsStyle
#what if I separate stand - and stand + and see the coverage of P-sites?

# inside of readsEndPlot.R

readsEndPlot(bamfile, CDS, toStartCodon=TRUE, fiveEnd=TRUE)
toStartCodon = TRUE
fiveEnd = TRUE
shift = 0
window = c(-29, 30)
readLen = 25:30
anchor = "5end"

# code of readEndPlot
    CDS <- CDS[CDS$isFirstExonInCDS]
    which <- promoters(CDS, upstream = abs(window)[1], downstream = abs(window)[2])
    h <- scanBamHeader(bamfile)
    seqs <- h$targets
    which <- ribosomeProfilingQC:::fixSeqlevelsStyle(which, names(seqs))
    which <- which[as.character(seqnames(which)) %in% names(seqs)]
    seqlevels(which) <- seqlevels(which)[seqlevels(which) %in%
        names(seqs)]
    param <- ScanBamParam(what = c("qwidth"), tag = character(0),
        flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE,
            isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE),
        which = which)
    open(bamfile)
    reads <- readGAlignments(bamfile, param = param)
    close(bamfile)

    chr1_CDS=CDS[seqnames(CDS)=="chr1"]
    reads <- reads[njunc(reads) == 0]
    reads <- narrow(reads)
    reads <- reads[qwidth(reads) %in% readLen]
    reads <- as(reads, "GRanges")
    reads <- ribosomeProfilingQC:::fixSeqlevelsStyle(reads, CDS)
    which <- ribosomeProfilingQC:::fixSeqlevelsStyle(which, CDS)
    x <- promoters(reads, upstream = 0, downstream = 1) # first nucleotide
    cvg <- coverage(x)
    w <- split(which, seqnames(which))
    cvg.sub <- unlist(lapply(cvg, sum))
    cvg <- cvg[cvg.sub > 0]
    seq <- intersect(names(cvg), names(w))
    vws <- Views(cvg[seq], w[seq])
    vws <- lapply(vws, function(.ele) {
        viewApply(.ele[width(.ele) == sum(abs(window)[c(1, 2)])],
            as.numeric)
    })
    vws <- do.call(cbind, vws)
    at <- seq(-abs(window[1]), abs(window[2]))
    at <- at[at != 0]
    if (length(dim(vws)) != 2) {
        stop("Not enough data available.")
    }
    height <- rowSums(vws)
    names(height) <- at

#
# separate + and - strand of CDS and observe the coverage
#