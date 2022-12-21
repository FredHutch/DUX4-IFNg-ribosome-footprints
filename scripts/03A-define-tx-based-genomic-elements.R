# define transcript-based genomic elements (features) such as 
# 5'UTR, 1st coding exons, and 3' UTR, round translation start site: TSS [-13, 13], 
# basically up/down 13 nucleotides from the start of the 1st coding exons.

# load library and workers
library(DESeq2)
library(GenomicAlignments)
library(plyranges)
library(tidyverse)
library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35

library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"

#
# define features
#
feature_5p  <- fiveUTRsByTranscript(txdb, use.name=TRUE)
feature_3p  <- threeUTRsByTranscript(txdb, use.name=TRUE)
feature_cds <- cdsBy(txdb, by="tx", use.name=TRUE)

.unique_UTRs <- function(utrs) {
    # exclude UTRs that are not unique
    exons_names_by_tx <- bplapply(utrs, function(gr) paste(gr$exon_name, collapse=","))
    keep_tx <- as.data.frame(unlist(exons_names_by_tx)) %>% rownames_to_column(var="tx_name") %>%
                 dplyr::rename(exons_names=`unlist(exons_names_by_tx)`) %>%
                 dplyr::distinct(exons_names, .keep_all=TRUE)
    utrs[keep_tx$tx_name]                 
}

unique_feature_5p <- .unique_UTRs(utrs=feature_5p)
unique_feature_3p <- .unique_UTRs(utrs=feature_3p)

# first coding exon and TSS [-13, 13]
first_exon_per_tx <- bplapply(feature_cds, function(gr) {
  st <- as.character(strand(gr[1]))
  if (st %in% c("+", "*")) x <- gr[1]
  if (st == "-") x <- gr[length(gr)]
  return(x)
})
first_exon_per_tx <- GRangesList(first_exon_per_tx)

.distinct_ranges <- function(gr) { # keeping distinct 1st coding exons only
    mcols(gr)$rng <- as.character(gr)
    gr_mcols <- as.data.frame(mcols(gr)) %>% dplyr::distinct(rng)
    gr <- gr[names(gr) %in% rownames(gr_mcols)]
}


first_exon <- .distinct_ranges(unlist(first_exon_per_tx))
#around_TSS <- promoters(first_exon, upstream=2, downstream=3)
around_TSS <- promoters(first_exon, upstream=13, downstream=14) # 27 nucleotides, 13 up/downstream of TSS

tx_based_features <- list(feature_5p=unique_feature_5p, feature_3p = unique_feature_3p,
                          first_exon = first_exon, around_TSS = around_TSS)
save(tx_based_features, file=file.path(pkg_dir, "data", "tx_based_features.rda"))                          
