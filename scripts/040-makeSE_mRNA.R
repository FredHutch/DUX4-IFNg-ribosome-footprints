# makeSE_mRNA.R
# Count reads on different genomic elements: 5'UTR (by tx), around TSS (by tx), first exons (by tx), 
# exons by genes, CDS by genes and transcripts, and 3' UTR (by tx)

library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35
library(GenomicAlignments)
library(tidyverse)

library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

#
# Define parameters and sample info
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
rna_bam_dir <- "/fh/fast/tapscott_s/CompBio/dchamm/DCH_Paired_RiboSeq_RNA-Seq_BAM_files"
sample_info <- data.frame(
  bam_files = bam_files <- list.files(rna_bam_dir, pattern=".bam$", full.names=TRUE)) %>%
    dplyr::mutate(sample_name = str_replace(basename(bam_files), ".bam", ""),
                  treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
    dplyr::mutate(treatment = str_replace(treatment, "-", "_")) %>%
    dplyr::mutate(treatment = recode(treatment, Untreated="untreated")) %>%
    dplyr::mutate(treatment = factor(treatment, levels=c("untreated", "DOX_pulse", "IFNg", "DOX_pulse_IFNg")))
rownames(sample_info) <- sample_info$sample_name                                  

#
# tools for tx_name
#
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

# 
# Define features by transcripts and genes
#

# (a) tx-based features
load(file.path(pkg_dir, "data", "tx_based_features.rda")) # 5'UTR/TSS/1stExon/3'UTR
around_TSS <- tx_based_features$around_TSS
feature_5p <- tx_based_features$feature_5p
feature_3p <- tx_based_features$feature_3p
first_exon <- tx_based_features$first_exon

# (b) gene-based features
cds <- cdsBy(txdb, by="gene")
exons <- exonsBy(txdb, by="gene")

#
# CDS expression by genes
#
rse_cds_mRNA <- summarizeOverlaps(features = cds, 
                                  reads=BamFileList(sample_info$bam_files),
                                  mode = "IntersectionStrict",
                                  inter.feature = TRUE, ignore.strand=TRUE)
colnames(rse_cds_mRNA) <- sample_info$sample_name                                  
colData(rse_cds_mRNA) <- as(sample_info, "DataFrame")
rowData(rse_cds_mRNA) <- gene.anno[rownames(rse_cds_mRNA),
                                   c("gene_id", "gene_type", "gene_name", "hgnc_id")]
save(rse_cds_mRNA, file=file.path(pkg_dir, "data", "rse_cds_mRNA.rda"))                                   

#
# exons expression by genes
#
rse_exons_mRNA <- summarizeOverlaps(features = exons, 
                                  reads=BamFileList(sample_info$bam_files),
                                  mode = "IntersectionStrict",
                                  inter.feature = TRUE, ignore.strand=TRUE)
colnames(rse_exons_mRNA) <- sample_info$sample_name                                  
colData(rse_exons_mRNA) <- as(sample_info, "DataFrame")
rowData(rse_exons_mRNA) <- gene.anno[rownames(rse_exons_mRNA),
                                   c("gene_id", "gene_type", "gene_name", "hgnc_id")]
save(rse_exons_mRNA, file=file.path(pkg_dir, "data", "rse_exons_mRNA.rda"))   

#
# first exon by transcripts (tx) - inherit sizeFactor from rse_cds_mRNA; must turn off inter.feature
#

rse_1st_exon_by_tx_mRNA <- 
  summarizeOverlaps(features = first_exon, 
                    reads=BamFileList(sample_info$bam_files),
                    mode = "IntersectionStrict",
                    inter.feature = FALSE, ignore.strand=TRUE, BPPARAM=bp_param)
colnames(rse_1st_exon_by_tx_mRNA) <- sample_info$sample_name                                  
colData(rse_1st_exon_by_tx_mRNA) <- as(sample_info, "DataFrame")
rse_1st_exon_by_tx_mRNA <- .get_row_data_txname(rse_1st_exon_by_tx_mRNA, txdb)
save(rse_1st_exon_by_tx_mRNA, file=file.path(pkg_dir, "data", "rse_1st_exon_by_tx_mRNA.rda")) 


# sanity check
psmb9_inx <- which(rowData(rse_1st_exon_by_tx_mRNA)$gene_name == "PSMB9")
psmb9_tx <- rownames(rse_1st_exon_by_tx_mRNA[psmb9_inx])
assays(rse_1st_exon_by_tx_mRNA[psmb9_tx])[[1]]

#
# around TSS
#
rse_TSS_by_tx_mRNA <- 
    summarizeOverlaps(features = around_TSS, 
                      reads=BamFileList(sample_info$bam_files),
                      mode = "Union",
                      inter.feature = FALSE, ignore.strand=TRUE, BPPARAM=bp_param)
colnames(rse_TSS_by_tx_mRNA) <- sample_info$sample_name                                  
colData(rse_TSS_by_tx_mRNA) <- as(sample_info, "DataFrame")
rse_TSS_by_tx_mRNA <- .get_row_data_txname(rse_TSS_by_tx_mRNA, txdb)
save(rse_TSS_by_tx_mRNA, file=file.path(pkg_dir, "data", "rse_TSS_by_tx_mRNA.rda"))                    

assays(rse_TSS_by_tx_mRNA[psmb9_tx])[[1]]

#
# 5' UTR
#
rse_5UTR_by_tx_mRNA <- 
    summarizeOverlaps(features = feature_5p, 
                      reads=BamFileList(sample_info$bam_files),
                      mode = "Union",
                      inter.feature = FALSE, ignore.strand=TRUE, BPPARAM=bp_param)
colnames(rse_5UTR_by_tx_mRNA) <- sample_info$sample_name                                  
colData(rse_5UTR_by_tx_mRNA) <- as(sample_info, "DataFrame")
rse_5UTR_by_tx_mRNA <- .get_row_data_txname(rse_5UTR_by_tx_mRNA, txdb)
save(rse_5UTR_by_tx_mRNA, file=file.path(pkg_dir, "data", "rse_5UTR_by_tx_mRNA.rda"))     

assays(rse_5UTR_by_tx_mRNA[psmb9_tx])[[1]]

#
# 3' UTR
#
rse_3UTR_by_tx_mRNA <- 
    summarizeOverlaps(features = feature_3p, 
                      reads=BamFileList(sample_info$bam_files),
                      mode = "Union",
                      inter.feature = FALSE, ignore.strand=TRUE, BPPARAM=bp_param)
colnames(rse_3UTR_by_tx_mRNA) <- sample_info$sample_name                                  
colData(rse_3UTR_by_tx_mRNA) <- as(sample_info, "DataFrame")
rse_3UTR_by_tx_mRNA <- .get_row_data_txname(rse_3UTR_by_tx_mRNA, txdb)
save(rse_3UTR_by_tx_mRNA, file=file.path(pkg_dir, "data", "rse_3UTR_by_tx_mRNA.rda"))                    
