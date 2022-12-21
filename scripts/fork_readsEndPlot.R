# fork_readsEndPlot.R
# this function replace ribosomeProfilingQC's readsEndPlot.R

fork_readsEndPlot <- function (bamfile, CDS, toStartCodon = TRUE, fiveEnd = TRUE,
    shift = 0, window = c(-29, 30), readLen = 25:30)
{
    stopifnot(is(bamfile, "BamFile"))
    stopifnot(is(fiveEnd, "logical"))
    stopifnot(is(shift, "numeric"))
    stopifnot(is(window, "numeric"))
    stopifnot(is(readLen, "numeric"))
    stopifnot(is(CDS, "GRanges"))
    if (length(CDS$internalPos) != length(CDS) || length(CDS$isFirstExonInCDS) !=
        length(CDS) || length(CDS$isLastExonInCDS) != length(CDS) ||
        length(CDS$tx_name) != length(CDS) || length(CDS$gene_id) !=
        length(CDS)) {
        stop("CDS must be output of prepareCDS")
    }
    if (toStartCodon) {
        CDS <- CDS[CDS$isFirstExonInCDS]
        which <- promoters(CDS, upstream = abs(window)[1], downstream = abs(window)[2])
    }
    else {
        CDS <- CDS[CDS$isLastExonInCDS]
        CDS <- ribosomeProfilingQC:::switch.strand(CDS)
        which <- promoters(CDS, upstream = abs(window)[1], downstream = abs(window)[2])
        which <- ribosomeProfilingQC:::switch.strand(which)
    }
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
    reads <- GenomicAlignments::readGAlignments(bamfile, param = param)
    close(bamfile)
    reads <- reads[njunc(reads) == 0]
    reads <- narrow(reads)
    reads <- reads[qwidth(reads) %in% readLen]
    reads <- as(reads, "GRanges")
    reads <- ribosomeProfilingQC:::fixSeqlevelsStyle(reads, CDS)
    which <- ribosomeProfilingQC:::fixSeqlevelsStyle(which, CDS)
    if (fiveEnd[1]) {
        x <- promoters(reads, upstream = 0, downstream = 1)
        if (shift[1] != 0) {
            x <- shift(x, shift = shift[1] - 1)
        }
    }
    else {
        x <- ribosomeProfilingQC:::switch.strand(reads)
        x <- promoters(x, upstream = 0, downstream = 1)
        x <- ribosomeProfilingQC:::switch.strand(x)
        if (shift[1] != 0) {
            x <- shift(x, shift = shift[1] + 1)
        }
    }
    cvg <- coverage(x)
    w <- split(which, seqnames(which))
    cvg.sub <- unlist(lapply(cvg, sum))
    cvg <- cvg[cvg.sub > 0]
    seq <- intersect(names(cvg), names(w))
    vws <- Views(cvg[seq], w[seq])
    #vws <- lapply(vws, function(.ele) {
    #    viewApply(.ele[width(.ele) == sum(abs(window)[c(1, 2)])],
    #        as.numeric)
    #})
    vws <- lapply(names(vws), function(seq_name) {
        .ele <- vws[[seq_name]]
        .which <- w[[seq_name]]
        view_matrix <- viewApply(.ele[width(.ele) == sum(abs(window)[c(1, 2)])],
            as.numeric)
        rev_ind <- which(strand(.which)=="-")    
        view_matrix[, rev_ind] <- view_matrix[nrow(view_matrix):1, rev_ind]
        view_matrix
    })

    vws <- do.call(cbind, vws)
    at <- seq(-abs(window[1]), abs(window[2]))
    at <- at[at != 0]
    if (length(dim(vws)) != 2) {
        stop("Not enough data available.")
    }
    height <- rowSums(vws)
    names(height) <- at
    barplot(height, las = 3, space = 0.5, ylab = "counts", xlab = paste("distance from",
        ifelse(fiveEnd, "5'", "3'"), "of reads to", ifelse(toStartCodon,
            "start", "stop"), "codon"))
    at1 <- which(at == -1)
    abline(v = at1 * 1.5 + 0.25, lty = 4)
    at <- seq((at1%%3), to = length(at), by = 3)
    at <- at[at != at1]
    at <- at * 1.5 + 0.25
    ymax <- max(height)
    segments(x0 = at, y0 = 0, x1 = at, y1 = ymax * 0.9, lty = 3,
        col = "gray80")
    return(invisible(height))
}
