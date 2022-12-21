# 06C-translation-changes-5UTR-plus-1st-exon.R
library(DESeq2)
library(readxl)
library(writexl)
library(tidyverse)
library(corrr)
library(plyranges)
library(wesanderson)
library(goseq)
library(Biostrings)
library(GenomicAlignments)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)

library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"

load(file.path(pkg_dir, "data", "p_sites.rda"))
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda")) # Ribo
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda")) # RNA
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda")) # Ribo



load(file.path(pkg_dir, "data", "tx_based_features.rda")) # 5'UTR/TSS/1stExon/3'UTR
first_exon <- tx_based_features$first_exon
first_exon <- S4Vectors::split(first_exon,  names(first_exon))
feature_5p <- tx_based_features$feature_5p


#
# (A) combine first coding exons ad 5p UTRs
#
tx_intersect <- intersect(names(feature_5p), names(first_exon))
tx_1stexon_only <- setdiff(names(first_exon), names(feature_5p))
tx_5p_only <- setdiff(names(feature_5p), names(first_exon))

comb <- bplapply(tx_intersect, function(x) {
    GenomicRanges::union(feature_5p[[x]], first_exon[[x]])
})
comb <- GRangesList(comb)
names(comb) <- tx_intersect
comb <- c(comb, first_exon[tx_1stexon_only])
comb <- c(comb, feature_5p[tx_5p_only])

#
# (B) count p-sites
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

rse_5UTR_1stExon_by_tx <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=comb, 
                    reads=pc, 
                    singleEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=FALSE)
}, BPPARAM=bp_param)
rse_5UTR_1stExon_by_tx <- do.call(cbind, rse_5UTR_1stExon_by_tx)
rse_5UTR_1stExon_by_tx <- .get_col_row_data(rse_5UTR_1stExon_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)
save(rse_5UTR_1stExon_by_tx, file=file.path(pkg_dir, "data", "rse_5UTR_1stExon_by_tx.rda"))

#
# (C) RNA expression
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

feature <- rowRanges(rse_5UTR_1stExon_by_tx)
si <- colData(rse_cds_mRNA) # RNA
rse_5UTR_1stExon_mRNA <- summarizeOverlaps(features = feature, 
                                  reads=BamFileList(si$bam_files),
                                  mode = "IntersectionStrict",
                                  inter.feature = TRUE, ignore.strand=TRUE, 
                                  BPPARAM = bp_param)
colnames(rse_5UTR_1stExon_mRNA) <- si$sample_name                                  
colData(rse_5UTR_1stExon_mRNA) <- si
rse_5UTR_1stExon_mRNA <- .get_row_data_txname(rse_5UTR_1stExon_mRNA, txdb)
save(rse_5UTR_1stExon_mRNA, file=file.path(pkg_dir, "data", "rse_5UTR_1stExon_mRNA.rda"))   

#
# (C) translation changes
#
source(file.path(pkg_dir, "scripts", "06A-tools_TE.R"))
source(file.path(pkg_dir, "scripts", "tools.R"))

load(file.path(pkg_dir, "data", "DUX4_induced_v2.rda")) # LFC > 1
load(file.path(pkg_dir, "data", "IFNg_induced_v2.rda")) # LFC > 1
# identify H1 and H2 variants and PSM subunits
histone_variants <- .histone_variants(gene.anno)
exclude_gene <- c(histone_variants$gene_id, DUX4_induced_v2$ensembl)
psm_name <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "Gene Lists_Proteasome_MHC-I.xlsx"),
                              sheet=1, range="A4:B49") 
PSM <- as.data.frame(gene.anno) %>% dplyr::select(gene_name, gene_id) %>%
  dplyr::filter(gene_name %in% psm_name$Symbol)

.do_seq_tidy_results_S6 <- function(dds_S6, exclude_gene, DUX4_induced_v2, IFNg_induced_v2, PSM) {
  dds_S6 <- dds_S6[!rowData(dds_S6)$gene_id %in% exclude_gene]        
  dds_S6 <- DESeq(dds_S6)
  S6_ribo_over_rna <- results(dds_S6, name="protocolRibo.treatmentDOX_pulse_IFNg", alpha=0.05) # translation efficiency changes
  print(summary(S6_ribo_over_rna))
  S6_rna <- results(dds_S6, name="treatment_DOX_pulse_IFNg_vs_IFNg", alpha=0.05)  
  S6_ribo <- results(dds_S6, contrast=list(c("treatment_DOX_pulse_IFNg_vs_IFNg", "protocolRibo.treatmentDOX_pulse_IFNg")), alpha=0.05)
  tidy_S6 <- .tidy_res_ribo_over_rna(dds=dds_S6, inter_res=S6_ribo_over_rna, 
                                   ribo_res=S6_ribo, rna_res=S6_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             rename(tx_name="ensembl") %>%
             dplyr::mutate(gene_id = rowData(dds_S6[tx_name])$gene_id, .after="tx_name") %>%                                   
             dplyr::mutate(IFNg_induced_v2 = gene_name %in% IFNg_induced_v2$gene_name) %>%
             dplyr::mutate(PSM = gene_id %in% PSM$gene_id)
}

# starts
dds_S6 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_5UTR_1stExon_mRNA, rse_ribo_by_tx = rse_5UTR_1stExon_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = c("IFNg", "DOX_pulse_IFNg"),
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S6 <- .do_seq_tidy_results_S6(dds_S6=dds_S6, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 


fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S6")

.scatter_plot_TE(res=tidy_S6, label_x_pos=0.8, label_y_pos=-3.8) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
      title="5'UTR +1st coding exon: Ribo-seq vs. RNA-seq (by tx)", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "S6_scatter_plot_5UTR_1stExon_tx.pdf"), width=6, height=5)

tidy_S6_down <- tidy_S6 %>% dplyr::filter(status=="down")
tidy_S6_up <- tidy_S6 %>% dplyr::filter(status=="up") # 651 transcripts / 325 genes

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR-plus-1st-exon")
write_xlsx(list(`all` = tidy_S6,
                `translation_upregulation` = tidy_S6_up, #
                `translation_downregulation` = tidy_S6_down),
           path=file.path(stats_dir, "S6-5UTR-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"))  