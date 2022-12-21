#
# load library
#
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(DEXSeq)
library(ribosomeProfilingQC)
library(hg38.HomoSapiens.Gencode.v35)
library(tidyverse)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35

library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# Define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))
source(file.path(pkg_dir, "scripts", "07A-tools_dexseq.R"))

sample_info <- as(colData(rse_cds_by_gene)[, c("sample_name", "treatment", "bam_files")], "data.frame") %>%
  dplyr::rename(condition = "treatment") %>%
  mutate(sample_name = str_replace_all(sample_name, "-", "_"))

#
# get exonic/cds parts
#
exonic_parts = exonicParts(txdb, linked.to.single.gene.only = TRUE) # GenomicFeatures
cds_parts <- .cds_parts(txdb, linked.to.single.gene.only = TRUE) # tools_dexseq.R

############################################################
#
# (1) Get P sites coordinates; limited to reads of length 30:33
#  
############################################################
ignore.strand <- FALSE
yieldSize <- 1e+07
best_p_site <- 13
reads_len <- c(26:29)
anchor <- "5end"

# (a) get p site coordinates
p_sites <- bplapply(sample_info$bam_files, function(f) {
  bam_file <- Rsamtools::BamFile(file = f)
  pc <- ribosomeProfilingQC::getPsiteCoordinates(bam_file, bestpsite = best_p_site,
                                                 anchor = anchor)
  pc.sub <- pc[pc$qwidth %in% reads_len]
})
names(p_sites) <- sample_info$sample_name
#save(p_sites, file=file.path(pkg_dir, "data", "p_sites.rda"))

###########################################################
#
# (2) Ribo-seq: count exon bins and make dxd from rse_exon
#
###########################################################
exonic_parts = exonicParts(txdb, linked.to.single.gene.only = TRUE)
cds_parts <- .cds_parts(txdb, linked.to.single.gene.only = TRUE)
# clean up 5'UTR and 3'UTR; keep 
rse <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=cds_parts, 
                    reads=pc, 
                    singelEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=FALSE)
})
rse_exon <- do.call(cbind, rse)
colData(rse_exon) <- as(sample_info, "DataFrame")
dxd_cds_parts <- DEXSeqDataSetFromSE(rse_exon, design= ~ sample + exon + condition:exon)
save(dxd_cds_parts, file=file.path(pkg_dir, "data", "dxd_cds_parts.rda"))

###########################################################
#
# (3) RNA-seq: count exon bins and make dxd from rse_exon
#
###########################################################
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda"))
sample_info <- as(colData(rse_cds_mRNA)[, c("sample_name", "treatment", "bam_files")], "data.frame") %>%
  dplyr::rename(condition = "treatment") %>%
  mutate(sample_name = str_replace_all(sample_name, "-", "_"))
rse <- summarizeOverlaps(features=cds_parts, reads=sample_info$bam_files, 
                         singleEnd=TRUE, inter.feature=FALSE, ignore.strand=FALSE,
                         mode="IntersectionStrict", BPPARAM=bp_param)  
colData(rse) <- as(sample_info, "DataFrame")                         
dxd_cds_parts_mRNA <- DEXSeqDataSetFromSE(rse, design= ~ sample + exon + condition:exon)
save(dxd_cds_parts_mRNA, file=file.path(pkg_dir, "data", "dxd_cds_parts_mRNA.rda"))

