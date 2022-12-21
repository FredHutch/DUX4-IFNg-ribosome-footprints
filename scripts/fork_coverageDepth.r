
fork_coverageDepth <- function (RPFs, RNAs, gtf, level = c("tx", "gene"), bestpsite = 13,
    readsLen = c(28, 29), anchor = "5end", region = "cds", ext = 5000,
    ...)
{
    stopifnot(is.character(gtf) || is(gtf, "TxDb"))
    level <- match.arg(level)
    anchor <- match.arg(anchor, choices = c("5end", "3end"))
    region <- region[1]
    if (!region %in% c("cds", "utr5", "utr3", "exon", "transcripts",
        "feature with extension")) {
        stop("region must be cds, utr5, utr3, exon,\n         transcripts, feature with extension")
    }
    cvgs <- list()
    txdb <- NULL
    if (is(gtf, "TxDb")) {
        txdb <- gtf
    }
    else {
        if (is.character(gtf)) {
            gtf <- gtf[1]
            suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtf,
                ...)))
        }
    }
    if (!is(txdb, "TxDb")) {
        stop("Can not determine annotations from gtf parameter.")
    }
    if (!missing(RPFs)) {
        stopifnot(is.character(RPFs))
        stopifnot(is.numeric(readsLen))
        stopifnot(is.numeric(bestpsite))
        cd <- getCvgs(files = RPFs, txdb = txdb, level = level,
            bestpsite = bestpsite, readsLen = readsLen, anchor = anchor,
            region = region, ext = ext)
        cvgs[["RPFs"]] <- cd
    }
    if (!missing(RNAs)) {
        stopifnot(is.character(RNAs))
        cd <- getCvgs(files = RNAs, txdb = txdb, level = level,
            region = region, ext = ext)
        cvgs[["mRNA"]] <- cd
    }
    cvgs
}
<bytecode: 0x3c563860>
<environment: namespace:ribosomeProfilingQC>