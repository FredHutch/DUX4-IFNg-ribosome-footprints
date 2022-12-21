# identify_DUX4_induced.R:
# This script identifies DUX4-induced by using RNA-seq reads mapping to 
# CDS. The criteria to call DUX4 targets is |lfc(DOX/control)| > 1 and FDR < 0.05

library(DESeq2)
library(readxl)
library(writexl)
library(tidyverse)
library(corrr)
library(plyranges)
library(wesanderson)

library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda")) # cds by genes

#
# (1) DUX4-induced: DOX_pulse vs. untreated
#
rna_S1 <- rse_cds_mRNA[, rse_cds_mRNA$treatment %in% c("untreated", "DOX_pulse")]
rna_S1 <- rna_S1[rowSums(assays(rna_S1)[["counts"]]) >= 12] 
rna_S1$treatment <- factor(rna_S1$treatment, levels=c("untreated", "DOX_pulse"))
dds_rna_S1 <- DESeqDataSet(rna_S1, design = ~ treatment)
dds_rna_S1 <- estimateSizeFactors(dds_rna_S1)
dds_rna_S1 <- DESeq(dds_rna_S1)
cnt <- as.data.frame(counts(dds_rna_S1, normalized=TRUE)) %>%
  rownames_to_column(var="ensembl")

# DUX4_induced with lfc >2 and padj < 0.05
DUX4_induced <- results(dds_rna_S1, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 2 ) %>%
  dplyr::arrange(padj) %>%
  dplyr::rename(ensembl="row") %>%
  dplyr::mutate(gene_name=rowData(dds_rna_S1[ensembl])$gene_name,
                gene_type=rowData(dds_rna_S1[ensembl])$gene_type, .before=2) %>%
  dplyr::mutate(status=if_else(log2FoldChange > 1, "up", "down")) %>%                
  dplyr::left_join(cnt, by="ensembl")
write_xlsx(DUX4_induced, path=file.path(pkg_dir, "stats", "differential_analysis", "deseq2_DUX4_induced_CDS_RNAseq.xlsx"))
save(DUX4_induced, file=file.path(pkg_dir, "data", "DUX4_induced.rda"))

# DUX4_induced with lfc > 1 and padj < 0.05
DUX4_induced_v2 <- results(dds_rna_S1, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 1 ) %>%
  dplyr::arrange(padj) %>%
  dplyr::rename(ensembl="row") %>%
  dplyr::mutate(gene_name=rowData(dds_rna_S1[ensembl])$gene_name,
                gene_type=rowData(dds_rna_S1[ensembl])$gene_type, .before=2) %>%
  dplyr::mutate(status=if_else(log2FoldChange > 1, "up", "down")) %>%                
  dplyr::left_join(cnt, by="ensembl")
#write_xlsx(DUX4_induced, path=file.path(pkg_dir, "stats", "differential_analysis", "deseq2_DUX4_induced_CDS_RNAseq.xlsx"))
save(DUX4_induced_v2, file=file.path(pkg_dir, "data", "DUX4_induced_v2.rda"))


#
# (2) MAplot
#
results(dds_rna_S1, alpha=0.05, tidy=TRUE) %>% 
  dplyr::mutate(log10_baseMean=log10(baseMean)) %>%
  dplyr::mutate(significant = padj < 0.05 & abs(log2FoldChange) > 2) %>%
  ggplot(aes(x=log10_baseMean, y=log2FoldChange, color=significant)) +
    geom_point(alpha=0.5, show.legend=FALSE) +
    theme(legend.position="none") +
    theme_bw() +
    scale_color_manual(values=c("grey", "#F21A00")) +
    labs(title="RNA-seq: DOX_pulse / untreated", 
         x=bquote(~log[10]~"(base mean)"), y="RNA-seq \nlog2FC (DOX_pulse / untreated)") +
         annotate("text", x=4, y=6.5, label="up: 476 \ndown: 171", hjust = 0, vjust=1) 
ggsave(file.path(pkg_dir, "figures", "RNA-S1-MAPlot.pdf"), width=4, height=3)  

#
# (3) S2 - IFNg induced: IFNg vs. untreated
#
rna_S2 <- rse_cds_mRNA[, rse_cds_mRNA$treatment %in% c("untreated", "IFNg")]
rna_S2 <- rna_S2[rowSums(assays(rna_S2)[["counts"]]) >= 12] 
rna_S2$treatement <- factor(rna_S2$treatment, levels=c("untreated", "IFNg"))
dds_rna_S2 <- DESeqDataSet(rna_S2, design = ~ treatment)
dds_rna_S2 <- estimateSizeFactors(dds_rna_S2)
dds_rna_S2 <- DESeq(dds_rna_S2)
IFNg_induced <- results(dds_rna_S2, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 2 ) %>%
  dplyr::arrange(padj) %>%
  dplyr::rename(ensembl="row") %>%
  dplyr::mutate(gene_name=rowData(dds_rna_S6[ensembl])$gene_name,
                gene_type=rowData(dds_rna_S6[ensembl])$gene_type, .before=2) %>%
  dplyr::mutate(counts=counts(dds_rna_S6[ensembl], normalized=TRUE))  

# output                
write_xlsx(IFNg_induced, path=file.path(pkg_dir, "stats", "differential_analysis", "deseq2_IFNg_induced_CDS_RNAseq.xlsx"))
save(IFNg_induced, file=file.path(pkg_dir, "data", "IFNg_induced.rda"))

IFNg_induced_v2 <- results(dds_rna_S2, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 1 ) %>%
  dplyr::arrange(padj) %>%
  dplyr::rename(ensembl="row") %>%
  dplyr::mutate(gene_name=rowData(dds_rna_S6[ensembl])$gene_name,
                gene_type=rowData(dds_rna_S6[ensembl])$gene_type, .before=2) %>%
  dplyr::mutate(counts=counts(dds_rna_S6[ensembl], normalized=TRUE))  
save(IFNg_induced_v2, file=file.path(pkg_dir, "data", "IFNg_induced_v2.rda"))