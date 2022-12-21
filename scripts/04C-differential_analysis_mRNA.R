# differential_analysis.R
#   subject: ribo-seq ribosome protected fragments (RPF)
#   size factor: use dds_cds_by_gene
#   analysis includes: 
#     (1) DESeq2 by CDS by gene (dds_cds_by_gene)
#     (2) 5'UTR by tx (rse_5p)
#   this script is adapted from my previous script, ribo_profiling.R
#  
#   output: * deseq2_cds -   
# 
#   figures:  

library(ribosomeProfilingQC)
library(tidyverse)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(writexl)
library(wesanderson)
library(viridis)  
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# load datasets and set parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
source(file.path(pkg_dir, "scripts", "tools.R"))
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda"))
load(file.path(pkg_dir, "data", "DUX4_altered.rda"))
load(file.path(pkg_dir, "data", "IFNg_altered.rda"))

# list of comparisons
list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))


.subset_by_treatment <- function(treatments, dds, threshold=12) {
    # treatments must include to treatments for comparison
    sub <- dds[, dds$treatment %in% treatments]
    sub <- sub[rowSums(counts(sub)) >= threshold, ]
    sub$treatment <- factor(sub$treatment, levels=treatments)
    sub <- DESeq(sub)
}

#
# (1c) DESeq2: no previously estimated size factors
#
tmp_cds <- DESeqDataSet(rse_cds_mRNA, design = ~treatment)  
deseq2_cds <- lapply(list_comp, .subset_by_treatment, dds=tmp_cds)    

summary <- map_dfr(deseq2_cds, function(x) {
  res <- results(x, lfcThreshold=0, alpha=0.05)
  as.data.frame(res) %>% dplyr::filter(padj < 0.05) %>%
    summarise(`LFC > 1 (up)` = sum(log2FoldChange > 1, na.rm=TRUE),
              `LFC < -1 (down)` = sum(log2FoldChange < -1, na.rm=TRUE))
}, .id="name") %>%
  add_column(comparison=map_chr(list_comp, function(x) paste0(x[2], "-vs-", x[1])), .after="name")

summary

tidy_results <- map(deseq2_cds, function(dds) {
  norm_cnt <- as.data.frame(counts(dds, normalized=TRUE)) %>%
    rownames_to_column(var="ensembl")
  results(dds, lfcThreshold=0, alpha=0.05, tidy=TRUE) %>%
    rename(ensembl="row") %>% #dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(gene_name = rowData(dds[ensembl])$gene_name, .before="baseMean") %>%
    dplyr::mutate(gene_type = rowData(dds[ensembl])$gene_type, .after="gene_name") %>%
    dplyr::mutate(significant = abs(log2FoldChange) > 1  & padj < 0.05) %>%
    dplyr::mutate(DUX4_altered = ensembl %in% DUX4_altered$ensembl) %>%
    dplyr::mutate(IFNg_altered = ensembl %in% IFNg_altered$ensembl) %>%
    left_join(norm_cnt, by="ensembl")
})

# output
file_name <- file.path(pkg_dir, "stats", "differential_analysis", "deseq2_cds_rnaseq.xlsx")
writexl::write_xlsx(
  list(DOX_pulse_vs_untreated      = tidy_results[["S1"]] ,
       IFNg_vs_untreated           = tidy_results[["S2"]],
       DOX_pulse_IFNg_vs_untreated = tidy_results[["S3"]],
       DOX_pulse_vs_IFNg           = tidy_results[["S4"]],
       DOX_pulse_IFNg_vs_DOX_pulse = tidy_results[["S5"]],
       DOX_pulse_IFNg_vs_IFNg      = tidy_results[["S6"]]),
  path=file_name)
