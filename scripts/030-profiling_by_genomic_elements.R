# profiling_by_genomic_elements.R (updated version of ribo_profiling.R)
# profiling by genomic elements (5'UTR by tx, CDS by tx and gene, and 3'UTR by tx)
#
# OUTPUT:
#  (a) data: `rse_cds_by_gene` and `rse_5p_3p_cds_by_tx`     
#  (b) PCA for `dds_cds_by_gene`, `dds_5p and dds_3p` -> ../figures/PCA
#
# NOTE: 
#  (a) CDS; 5p, TSS, 1st exon, and 3p profiling by gene defines size factors for all the profilings. 
#  (b) make genomic element specific profiling by transcripts (tx)
#  (c) CDS profiling by gene will be used to estimate size factors and 
#      translation efficiency 
#  (d) `rse_cds_by_gene.rda` replaces the depricated `dds_r1.rda` or `rse.rda`
#  (e) `rse_5p_3p_cds_by_tx.` replaces `dds_5p_3p_cds.rda` and `ribo_profile.rda`
#

# load library and workers
library(DESeq2)
library(GenomicAlignments)
library(plyranges)
library(tidyverse)
library(ggrepel)
library(writexl)
library(Rsamtools)
library(ribosomeProfilingQC)
library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35

library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

#
# define parameter and sample info
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
fig_dir <- file.path(pkg_dir, "figures", "PCA")
source(file.path(pkg_dir, "scripts", "tools.R"))
bam_dir <- file.path(scratch_dir, "bam", "merged_bam_runs")
bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)
sample_info <- data.frame(
  bam_files = bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)) %>%
    dplyr::mutate(sample_name = str_replace(basename(bam_files), ".bam", ""),
                  treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
    dplyr::mutate(treatment = str_replace(treatment, "-", "_")) %>%                
    dplyr::mutate(treatment = factor(treatment, levels=c("untreated", "DOX_pulse", "IFNg", "DOX_pulse_IFNg")))

#
# tools
#
.get_col_row_data <- function(rse, txdb, dds_cds_by_gene) {
  # colnames and append sample_info
  colnames(rse) <- colnames(dds_cds_by_gene)
  colData(rse) <- colData(dds_cds_by_gene)
  # rowData: tx_name -> gene_id, gene_name, gene_type
  tx_name <- rownames(rse)
  df <- AnnotationDbi::select(txdb, keys=tx_name, columns="GENEID", keytype="TXNAME",
                              multiVals="first") %>% as.data.frame() %>%
    dplyr::distinct(TXNAME, .keep_all=TRUE) %>%
    dplyr::rename(tx_name=TXNAME, gene_id=GENEID) %>%
    dplyr::left_join(as.data.frame(gene.anno), by="gene_id") %>%
    dplyr::select(tx_name, gene_id, gene_type, gene_name, hgnc_id)                       
  rownames(df) <- df$tx_name
  rowData(rse) <- df[rownames(rse), ]
  
  return(rse)
}

.get_row_data_txname <- function(rse, txdb) {
  tx_name <- rownames(rse)
  df <- AnnotationDbi::select(txdb, keys=tx_name, columns="GENEID", keytype="TXNAME",
                              multiVals="first") %>% as.data.frame() %>%
    dplyr::distinct(TXNAME, .keep_all=TRUE) %>%
    dplyr::rename(tx_name=TXNAME, gene_id=GENEID) %>%
    dplyr::left_join(as.data.frame(gene.anno), by="gene_id") %>%
    dplyr::select(tx_name, gene_id, gene_type, gene_name, hgnc_id)                       
  rownames(df) <- df$tx_name
  rowData(rse) <- df[rownames(rse), ]
  
  return(rse)
}


############################################################
#
# (1) Get P sites coordinates; limited to reads of length 26:29
#  save to "p_sites.rda"
############################################################
ignore.strand <- FALSE
yieldSize <- 1e+07
best_p_site <- 13
reads_len <- c(26:29)
anchor <- "5end"

# (a) get p site coordinates (we load from )
reads <- bplapply(sample_info$bam_files, function(f) {
  bam_file <- BamFile(file = f)
  pc <- getPsiteCoordinates(bam_file, bestpsite = best_p_site,
                            anchor = anchor)
  pc.sub <- pc[pc$qwidth %in% reads_len]
})
names(reads) <- sample_info$sample_name
p_sites <- reads
save(p_sites, file=file.path(pkg_dir, "data", "p_sites.rda"))

##############################
#
# (2) CDS by gene profiling; convert to dds
#
##############################
cds_by_gene <- cdsBy(txdb, by="gene")
ignore.strand <- FALSE

# profiling
rse_cds_by_gene <- bplapply(reads, function(pc) {
  summarizeOverlaps(features=cds_by_gene, 
                    reads=pc, 
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_cds_by_gene <- do.call(cbind, rse_cds_by_gene)

# tidy colData and rowData
colnames(rse_cds_by_gene) <- sample_info$sample_name
colData(rse_cds_by_gene) <- append(colData(rse_cds_by_gene), 
                                   as(sample_info, "DataFrame"))
rowData(rse_cds_by_gene) <- gene.anno[rownames(rse_cds_by_gene), 
                                      c("gene_id", "gene_type", "gene_name", "hgnc_id")]
save(rse_cds_by_gene, file=file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))                                      

# dds
dds_cds_by_gene <- DESeqDataSet(rse_cds_by_gene, design = ~treatment)              
dds_cds_by_gene <- dds_cds_by_gene[rowSums(counts(dds_cds_by_gene)) > 12]
dds_cds_by_gene <- estimateSizeFactors(dds_cds_by_gene)
save(dds_cds_by_gene, file=file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

gg_pca <- .make_pca_by_treatment(dds_cds_by_gene) +
  theme(legend.position="bottom") +
  labs(title="PCA on CDS (by gene) profiling")
pdf(file.path(pkg_dir, "figures", "PCA", "PCA_cds_by_gene.pdf"), width=5, height=4)
plot(gg_pca)
dev.off()


############################################################
#
# (3) Profile genomic features defined based on transcripts (tx): 5'UTR/TSS/1st coding exon/3'UTR
#     The normalization of p-site counts on these features is based on size factors for CDS (by genes) with DESeq2
#
############################################################

load(file.path(pkg_dir, "data", "p_sites.rda")) # same as "reads" made previously
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda")) # need the sizeFactor and column data
load(file.path(pkg_dir, "data", "tx_based_features.rda")) # 5'UTR/TSS/1stExon/3'UTR
around_TSS <- tx_based_features$around_TSS
feature_5p <- tx_based_features$feature_5p
feature_3p <- tx_based_features$feature_3p
first_exon <- tx_based_features$first_exon
feature_cds <- cdsBy(txdb, by="tx", use.name=TRUE)
ignore.strand <- FALSE
inter.feature <- FALSE

#
# (3a) TSS [-13, 13] tx_based_features$around TSS
#
rse_TSS_by_tx <- bplapply(p_sites, function(pc) {
    summarizeOverlaps(features=around_TSS, 
                    reads=pc, 
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_TSS_by_tx <- do.call(cbind, rse_TSS_by_tx)
rse_TSS_by_tx <- .get_col_row_data(rse_TSS_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)

save(rse_TSS_by_tx, file=file.path(pkg_dir, "data", "rse_TSS_by_tx.rda"))                            

#
# (3b) 5'UTR profiling by tx; use size factor from dds_cds_by_gene
#
rse_5UTR_by_tx <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=feature_5p, 
                    reads=pc, 
                    singleEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_5UTR_by_tx <- do.call(cbind, rse_5UTR_by_tx)
rse_5UTR_by_tx <- .get_col_row_data(rse_5UTR_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)
save(rse_5UTR_by_tx, file=file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))

# PCA
dds_5p <- DESeqDataSet(rse_5UTR_by_tx, design = ~treatment)              
dds_5p <- dds_5p[rowSums(counts(dds_5p)) > 12]
gg_pca <- .make_pca_by_treatment(dds_5p, sample_labels=TRUE) +
  theme(legend.position="none") +
  labs(title="PCA on 5'UTR (by tx) profiling")
pdf(file.path(pkg_dir, "figures", "PCA", "PCA_5UTR_by_tx.pdf"), width=5, height=4)
plot(gg_pca)
dev.off()

#
# (3c) 3'UTR profiling by tx; use size factor from dds_cds_by_gene
#
rse_3UTR_by_tx <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=feature_3p, 
                    reads=pc, 
                    singleEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_3UTR_by_tx <- do.call(cbind, rse_3UTR_by_tx)
rse_3UTR_by_tx <- .get_col_row_data(rse_3UTR_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)
save(rse_3UTR_by_tx, file=file.path(pkg_dir, "data", "rse_3UTR_by_tx.rda"))

#PCA
dds_3p <- DESeqDataSet(rse_3UTR_by_tx, design = ~treatment)              
dds_3p <- dds_3p[rowSums(counts(dds_3p)) > 12]
gg_pca <- .make_pca_by_treatment(dds_3p, sample_labels=TRUE) +
  theme(legend.position="none") +
  labs(title="PCA on 3'UTR (by tx) profiling")
pdf(file.path(pkg_dir, "figures", "PCA", "PCA_3UTR_by_tx.pdf"), width=5, height=4)
plot(gg_pca)
dev.off()

#
# (3d) 1st coding exon profiling by tx; use size factor from dds_cds_by_gene1st exons
#
rse_1st_exon_by_tx <- bplapply(p_sites, function(pc) {
    summarizeOverlaps(features=first_exon, 
                    reads=pc, 
                    inter.feature=FALSE,
                    ignore.strand=FALSE)
})
rse_1st_exon_by_tx <- do.call(cbind, rse_1st_exon_by_tx)
rse_1st_exon_by_tx <- .get_col_row_data(rse=rse_1st_exon_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)      
save(rse_1st_exon_by_tx, file=file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda")) 

