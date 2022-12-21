.assign_utr_frame <- function(pc, CDS, txdb) {
  # input: 
  #     pc is the p_sites made by ribosomeProfilingQC::getPC..
  #     CDS must be the product of ribosomeProfilingQC::prepareCDS
  # output: 
  #     three vector of p-sites in CDS, 5' UTR and 3' UTR in which readingFrame 
  #     are assigned to 5'UTR and 3'UTR p-sites

  # get start codon
  first_exon <- CDS[CDS$isFirstExonInCDS]
  start_codon <- promoters(first_exon, upstream=0, downstream=1)

  ov_cds <- findOverlaps(pc, CDS, ignore.strand=FALSE)
  utr_candidates <- pc[-unique(queryHits(ov_cds))]

  utr_regions <- c("utr5", "utr3")
  utr_pc <- lapply(utr_regions, function(utr_region) {
    # define features
    if (utr_region == "utr5")
      features <- fiveUTRsByTranscript(txdb, use.name = TRUE)
    if (utr_region == "utr3")   
      features <- threeUTRsByTranscript(txdb, use.name = TRUE)  
      
    features <- unlist(features)
    features$tx_name <- names(features)
   
    # get overlaps
    ov_utr <- findOverlaps(utr_candidates, features, ignore.strand=FALSE) %>% 
      as.data.frame() %>%
      dplyr::distinct(queryHits, .keep_all=TRUE) %>%
      dplyr::mutate(tx_name = features[subjectHits]$tx_name) %>% 
      left_join(as.data.frame(start_codon), by="tx_name") # add start codon of the subjectHits
    
    start_codon_gr <-  ov_utr %>%
      dplyr::select(seqnames, start, end, strand, width, tx_name, queryHits, subjectHits) %>% 
      plyranges::as_granges()
    width(start_codon_gr) <- 0  
    utr_pc <- utr_candidates[start_codon_gr$queryHits] # keep width 0
    width(utr_pc) <- 0
    dist <- distance(utr_pc, start_codon_gr, ignore.strand = FALSE) 
    utr_pc$readingFrame <- dist %% 3 
    width(utr_pc) <- 1

    return(utr_pc)   
  })
  names(utr_pc) <- utr_regions
  return(list(cds=pc[unique(queryHits(ov_cds))], utr5=utr_pc$utr5, utr3=utr_pc$utr3))
}