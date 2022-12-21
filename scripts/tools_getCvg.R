.getCvgs_2 <- function(reads, anchor, ext, region) {
    # pc per frame per region
    ignore.strand <- FALSE
    features <- switch(region, cds = cds(txdb, columns = c("gene_id",
        "tx_name")), utr5 = {
        f <- fiveUTRsByTranscript(txdb, use.name = TRUE)
        fs <- unlist(f)
        mcols(fs) <- DataFrame(tx_name = rep(names(f), lengths(f)))
        suppressMessages(id_map <- AnnotationDbi::select(txdb, keys = unique(fs$tx_name),
            columns = c("TXNAME", "GENEID"), keytype = "TXNAME"))
        fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME),
            "GENEID"]
        fs
    }, utr3 = {
        f <- threeUTRsByTranscript(txdb, use.name = TRUE)
        fs <- unlist(f)
        mcols(fs) <- DataFrame(tx_name = rep(names(f), lengths(f)))
        suppressMessages(id_map <- AnnotationDbi::select(txdb, keys = unique(fs$tx_name),
            columns = c("TXNAME", "GENEID"), keytype = "TXNAME"))
        fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME),
            "GENEID"]
        fs
    }, exon = exons(txdb, columns = c("gene_id", "tx_name")),
        transcripts = transcripts(txdb, columns = c("gene_id",
            "tx_name")), `feature with extension` = {
            f <- transcripts(txdb, columns = c("gene_id", "tx_name"))
            e <- exons(txdb, columns = c("gene_id", "tx_name"))
            tx_name <- e$tx_name
            e <- rep(e, lengths(tx_name))
            e$tx_name <- unlist(tx_name)
            CDS <- cds(txdb, columns = c("gene_id", "tx_name"))
            tx_name <- CDS$tx_name
            CDS <- rep(CDS, lengths(tx_name))
            CDS$tx_name <- unlist(tx_name)
            getExt <- function(ext, f, e, start = TRUE) {
                exts <- flank(f, ext, start = start)
                diffs <- setdiff(exts, f)
                ols <- findOverlaps(diffs, exts, type = "within")
                es <- diffs[queryHits(ols)]
                mcols(es) <- mcols(exts[subjectHits(ols)])
                dists <- distance(es, f[match(es$tx_name, f$tx_name)])
                es <- es[dists == 0]
                ol <- findOverlaps(es, e, maxgap = 1)
                ol <- ol[es[queryHits(ol)]$tx_name == e[subjectHits(ol)]$tx_name]
                es <- es[queryHits(ol)]
                e1 <- e[subjectHits(ol)]
                dist <- distance(es, e1)
                es <- es[dist == 0]
                e1 <- e1[dist == 0]
                ol <- ol[dist == 0]
                start(es) <- ifelse(start(es) < start(e1), start(es),
                  start(e1))
                end(es) <- ifelse(end(es) < end(e1), end(e1),
                  end(es))
                e[subjectHits(ol)] <- es
                e
            }
            exStart <- getExt(ext, f, e, start = TRUE)
            exEnd <- getExt(ext, f, e, start = FALSE)
            f_ex <- c(f, CDS, exStart, exEnd)
            f_ex <- split(f_ex, f_ex$tx_name)
            f_ex <- disjoin(f_ex)
            f_ex <- unlist(f_ex)
            mcols(f_ex) <- mcols(f)[match(names(f_ex), f$tx_name),
                ]
            f_ex
        })
    features <- switch(level, gene = {
        f <- rep(features, lengths(features$gene_id))
        mcols(f) <- DataFrame(feature_id = unlist(features$gene_id))
        f[!is.na(f$feature_id)]
    }, tx = {
        f <- rep(features, lengths(features$tx_name))
        mcols(f) <- DataFrame(feature_id = unlist(features$tx_name))
        f[!is.na(f$feature_id)]
    })
    # 
    cvgs <- coverage(reads)



}

.getCvgs <- function (reads, txdb, level, anchor, region, ext)
{
    ignore.strand <- FALSE
    features <- switch(region, cds = cds(txdb, columns = c("gene_id",
        "tx_name")), utr5 = {
        f <- fiveUTRsByTranscript(txdb, use.name = TRUE)
        fs <- unlist(f)
        mcols(fs) <- DataFrame(tx_name = rep(names(f), lengths(f)))
        suppressMessages(id_map <- AnnotationDbi::select(txdb, keys = unique(fs$tx_name),
            columns = c("TXNAME", "GENEID"), keytype = "TXNAME"))
        fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME),
            "GENEID"]
        fs
    }, utr3 = {
        f <- threeUTRsByTranscript(txdb, use.name = TRUE)
        fs <- unlist(f)
        mcols(fs) <- DataFrame(tx_name = rep(names(f), lengths(f)))
        suppressMessages(id_map <- AnnotationDbi::select(txdb, keys = unique(fs$tx_name),
            columns = c("TXNAME", "GENEID"), keytype = "TXNAME"))
        fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME),
            "GENEID"]
        fs
    }, exon = exons(txdb, columns = c("gene_id", "tx_name")),
        transcripts = transcripts(txdb, columns = c("gene_id",
            "tx_name")), `feature with extension` = {
            f <- transcripts(txdb, columns = c("gene_id", "tx_name"))
            e <- exons(txdb, columns = c("gene_id", "tx_name"))
            tx_name <- e$tx_name
            e <- rep(e, lengths(tx_name))
            e$tx_name <- unlist(tx_name)
            CDS <- cds(txdb, columns = c("gene_id", "tx_name"))
            tx_name <- CDS$tx_name
            CDS <- rep(CDS, lengths(tx_name))
            CDS$tx_name <- unlist(tx_name)
            getExt <- function(ext, f, e, start = TRUE) {
                exts <- flank(f, ext, start = start)
                diffs <- setdiff(exts, f)
                ols <- findOverlaps(diffs, exts, type = "within")
                es <- diffs[queryHits(ols)]
                mcols(es) <- mcols(exts[subjectHits(ols)])
                dists <- distance(es, f[match(es$tx_name, f$tx_name)])
                es <- es[dists == 0]
                ol <- findOverlaps(es, e, maxgap = 1)
                ol <- ol[es[queryHits(ol)]$tx_name == e[subjectHits(ol)]$tx_name]
                es <- es[queryHits(ol)]
                e1 <- e[subjectHits(ol)]
                dist <- distance(es, e1)
                es <- es[dist == 0]
                e1 <- e1[dist == 0]
                ol <- ol[dist == 0]
                start(es) <- ifelse(start(es) < start(e1), start(es),
                  start(e1))
                end(es) <- ifelse(end(es) < end(e1), end(e1),
                  end(es))
                e[subjectHits(ol)] <- es
                e
            }
            exStart <- getExt(ext, f, e, start = TRUE)
            exEnd <- getExt(ext, f, e, start = FALSE)
            f_ex <- c(f, CDS, exStart, exEnd)
            f_ex <- split(f_ex, f_ex$tx_name)
            f_ex <- disjoin(f_ex)
            f_ex <- unlist(f_ex)
            mcols(f_ex) <- mcols(f)[match(names(f_ex), f$tx_name),
                ]
            f_ex
        })
    features <- switch(level, gene = {
        f <- rep(features, lengths(features$gene_id))
        mcols(f) <- DataFrame(feature_id = unlist(features$gene_id))
        f[!is.na(f$feature_id)]
    }, tx = {
        f <- rep(features, lengths(features$tx_name))
        mcols(f) <- DataFrame(feature_id = unlist(features$tx_name))
        f[!is.na(f$feature_id)]
    })
    cvgs <- lapply(reads, GenomicRanges::coverage)
    seqs <- lapply(cvgs, function(.ele) names(.ele)[vapply(.ele,
        sum, FUN.VALUE = 0) > 0])
    seqs <- unique(unlist(seqs))
    seqs <- intersect(seqlevels(features), seqs)
    features <- features[as.character(seqnames(features)) %in%
        seqs]
    seqlevels(features) <- seqlevels(features)[seqlevels(features) %in%
        seqs]
    features1 <- disjoin(features, with.revmap = TRUE)
    f <- rep(features1, lengths(features1$revmap))
    f$feature_id <- features$feature_id[unlist(features1$revmap)]
    f$revmap <- NULL
    f <- f[!duplicated(f) | !duplicated(f$feature_id)]
    f <- f[order(f$feature_id, as.character(seqnames(f)), ifelse(as.character(strand(f)) ==
        "-", -1, 1) * start(f))]
    f <- split(f, f$feature_id)
    coverages <- lapply(reads, coverageByTranscript, transcripts = f,
        ignore.strand = ignore.strand)
    cd <- cvgd(coverage = coverages, granges = f)
    return(cd)
}
