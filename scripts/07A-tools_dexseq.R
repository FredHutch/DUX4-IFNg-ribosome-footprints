.cds_parts <- function(txdb, linked.to.single.gene.only = TRUE) {
  if (!isTRUEorFALSE(linked.to.single.gene.only))
        stop("'linked.to.single.gene.only' must be TRUE or FALSE")
    ex <- .tidy_cds(txdb, drop.geneless = linked.to.single.gene.only)
    GenomicFeatures:::.break_in_parts(ex, linked.to.single.gene.only, extra_mcol = "exonic_part")
}

.tidy_cds <- function (txdb, drop.geneless = FALSE)
{
    tx <- GenomicFeatures:::tidyTranscripts(txdb, drop.geneless = drop.geneless)
    cds_by_tx <- cdsBy(txdb, by="tx")
    names(tx) <- mcols(tx)$tx_id
    tx <- tx[names(cds_by_tx)]
    ans <- unlist(cds_by_tx, use.names = FALSE)
    idx <- rep(seq_along(tx), lengths(cds_by_tx))
    mcols(ans) <- cbind(mcols(tx)[idx, , drop = FALSE], mcols(ans))
    names(mcols(ans))[c(4, 5)] <- c("exon_id", "exon_name")
    ans
}

.do_dexseq <- function(dxd, treatments) {
    # subset specific treatments; filter should be done before calling this function
    sub <- dxd[, dxd$condition %in% treatments]
    sub$condition <- factor(sub$condition, levels=treatments)
    # sub <- sub[rowSums(counts(sub)) >= 12, ]
    # test for differential exon usage and exon fold changes
    sub <- DEXSeq::estimateSizeFactors(sub)
    sub <- DEXSeq::estimateDispersions(sub)
    sub <-  DEXSeq::testForDEU(sub)
    sub <- estimateExonFoldChanges(sub, fitExpToVar="condition")
    return(sub)
}

.get_exons_by_quantile <- function(.x, n_quantile = c(1, 2, 3, 4), 
                                   type = c("any", "start", "end", "within", "equal")) {
  require(plyranges)                                       
  # n_quantile: 1 - first, 2 - second, 3 - third, 4 - fourth                                       
  n_quantile <- match.arg(as.character(n_quantile), choices=c("1", "2", "3", "4"))
  n_quantile <- as.integer(n_quantile)        
  type <- match.arg(type)                     

  strand <- as.character(.x$genomicData.strand)[1]
  l <- sum(.x$genomicData.width)
  w <- cumsum(.x$genomicData.width)
  quarter_length <- sum(.x$genomicData.width) /4
  ir_q <- data.frame(start=c(1, w[-length(w)]+1), 
                     width=.x$genomicData.width) %>% 
            plyranges::as_iranges()

  # if "-" reverse n_quantile
  if (strand == "-") {
      q <- c(1, 2, 3, 4)
      n_quantile <- rev(q)[n_quantile]
  }

  ir_s <- data.frame(start=max(1, (round(quarter_length * (n_quantile - 1)) - 1 )), 
                     end=(ceiling(quarter_length * n_quantile) + 1)) %>% 
            plyranges::as_iranges()

  ol <- findOverlaps(ir_q, ir_s, type=type)
  selected_idx <- queryHits(ol)

  # if type="within", and n_quantile = 1 and queryHitls(ol) == ingeter(0) => and ol is empty, then change type to "any"
  # need unit test
  if (identical(selected_idx, integer(0)) & type == "within") {
     selected_idx <- ifelse(strand=="+", 1, length(ir_q))
  }

  .x[selected_idx, ]
}