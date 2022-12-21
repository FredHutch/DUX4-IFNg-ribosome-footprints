# random_codes.R

getPsiteCoordinates <- function (bamfile, bestpsite, anchor = "5end")
{
    stopifnot(is(bamfile, "BamFile"))
    stopifnot(is.numeric(bestpsite))
    anchor <- match.arg(anchor, choices = c("5end", "3end"))
    offset <- round(bestpsite[1]) - 1
    param <- ScanBamParam(what = c("qwidth"), tag = character(0),
        flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE,
            isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE))
    reads <- GAlignments()
    open(bamfile)
    while (length(chunk <- readGAlignments(bamfile, param = param))) {
        reads <- c(reads, shiftReads(chunk, shift = offset, anchor = anchor))
    }
    close(bamfile)
    reads <- as(reads, "GRanges")
    reads <- promoters(reads, upstream = 0, downstream = 1)
    reads$Psite <- bestpsite
    if (anchor == "3end")
        reads$Psite <- reads$qwidth - reads$Psite
    metadata(reads) <- list(totalReads = length(reads))
    reads
}



# 
load(file.path("~/tapscott/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2/data/p_sites.rda"))
freq = lapply(p_sites[1:2], function(pc) {
    length_freq = table(pc$qwidth)
})

library(ribosomeProfilingQC)
library(tidyverse)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
fig_dir <- file.path(pkg_dir, "figures", "QC")
source(file.path(pkg_dir, "scripts", "tools.R"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

#
# sample information and CDS
#
bam_dir <- file.path(scratch_dir, "bam", "merged_bam_runs")
bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)
sample_info <- data.frame(
  bam_files = bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)) %>%
    dplyr::mutate(sample_name = str_replace(basename(bam_files), ".bam", ""),
                  treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
    dplyr::mutate(treatment = factor(treatment, levels=c("untreated", "DOX-pulse", "IFNg", "DOX-pulse_IFNg")))

CDS <- prepareCDS(txdb)
read_length_freq <- bplapply(sample_info$bam_files[1:3], function(x) {
  bam_file <- BamFile(x)
  param <- ScanBamParam(what = c("qwidth"), tag = character(0),
  flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE,
           isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE))
  #open(bam_file)
  reads <- readGAlignments(bam_file, param = param)
  #close(bamfile)
  reads <- reads[njunc(reads) == 0]
  
  #p_site <- estimatePsite(bam_file, CDS, genome)
  #pc <- getPsiteCoordinates(bam_file, bestpsite = p_site)
  #read_length_freq <- summaryReadsLength(pc, widthRange = c(20:39), plot=FALSE)
})
names(read_length_freq) <- sample_info$sample_name
