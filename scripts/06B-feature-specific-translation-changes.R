#
# Translation changes analysis part II: transcript-specific (5'UTR, TSS [-2, 2], 1st exon, 3'UTR) translation chanages
# Previously in part I, we analyzed translation changes in CDS. In part II to 
# identify translation changes in other features including TSS, 1st exon, 5'UTR for 
# each of the transcripts, for S6 comparison.
#

#
# load library
#
library(DESeq2)
library(readxl)
library(writexl)
library(tidyverse)
library(corrr)
library(plyranges)
library(wesanderson)
library(goseq)
library(Biostrings)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)

library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S6")
source(file.path(pkg_dir, "scripts", "06A-tools_TE.R"))
source(file.path(pkg_dir, "scripts", "tools.R"))

# rse and dds CDS by genes
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda")) # Ribo
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda")) # RNA
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

# colors
my_color <- wesanderson::wes_palette("Darjeeling1", n=5)[2:5]
names(my_color) <- levels(rse_cds_by_gene$treatment)

list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))

my_color <- wesanderson::wes_palette("Darjeeling1", n=5)[2:5]
names(my_color) <- levels(rse_cds_by_gene$treatment)
#
# load DUX4_altered and INFg_altered: what's the difference between first and secon version?
#
load(file.path(pkg_dir, "data", "DUX4_induced.rda")) # LFC > 2
load(file.path(pkg_dir, "data", "IFNg_induced.rda")) # LFC > 2
load(file.path(pkg_dir, "data", "DUX4_induced_v2.rda")) # LFC > 1
load(file.path(pkg_dir, "data", "IFNg_induced_v2.rda")) # LFC > 1

#
# identify H1 and H2 variants and PSM subunits
#
histone_variants <- .histone_variants(gene.anno)
exclude_gene <- c(histone_variants$gene_id, DUX4_induced_v2$ensembl)
psm_name <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "Gene Lists_Proteasome_MHC-I.xlsx"),
                              sheet=1, range="A4:B49") 
PSM <- as.data.frame(gene.anno) %>% dplyr::select(gene_name, gene_id) %>%
  dplyr::filter(gene_name %in% psm_name$Symbol)

#
# tools
#
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

.do_seq_tidy_results_S2 <- function(dds_S2, exclude_gene, DUX4_induced_v2, IFNg_induced_v2, PSM) {
  dds_S2 <- dds_S2[!rowData(dds_S2)$gene_id %in% exclude_gene]        
  dds_S2 <- DESeq(dds_S2)
  S2_ribo_over_rna <- results(dds_S2, name="protocolRibo.treatmentIFNg", alpha=0.05) # translation efficiency changes
  print(summary(S2_ribo_over_rna))
  S2_rna <- results(dds_S2, name="treatment_IFNg_vs_untreated", alpha=0.05)  
  S2_ribo <- results(dds_S2, contrast=list(c("treatment_IFNg_vs_untreated", "protocolRibo.treatmentIFNg")), alpha=0.05)
  tidy_S2 <- .tidy_res_ribo_over_rna(dds=dds_S2, inter_res=S2_ribo_over_rna, 
                                   ribo_res=S2_ribo, rna_res=S2_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             rename(tx_name="ensembl") %>%
             dplyr::mutate(gene_id = rowData(dds_S2[tx_name])$gene_id, .after="tx_name") %>%                                   
             dplyr::mutate(IFNg_induced_v2 = gene_name %in% IFNg_induced_v2$gene_name) %>%
             dplyr::mutate(PSM = gene_id %in% PSM$gene_id)
}

.do_seq_tidy_results_S1 <- function(dds_S1, exclude_gene, DUX4_induced_v2, IFNg_induced_v2, PSM) {
  dds_S1 <- dds_S1[!rowData(dds_S1)$gene_id %in% exclude_gene]        
  dds_S1 <- DESeq(dds_S1)
  S1_ribo_over_rna <- results(dds_S1, name="protocolRibo.treatmentDOX_pulse", alpha=0.05) # translation efficiency changes
  print(summary(S1_ribo_over_rna))
  S1_rna <- results(dds_S1, name="treatment_DOX_pulse_vs_untreated", alpha=0.05)  
  S1_ribo <- results(dds_S1, contrast=list(c("treatment_DOX_pulse_vs_untreated", "protocolRibo.treatmentDOX_pulse")), alpha=0.05)
  tidy_S1 <- .tidy_res_ribo_over_rna(dds=dds_S1, inter_res=S1_ribo_over_rna, 
                                   ribo_res=S1_ribo, rna_res=S1_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             rename(tx_name="ensembl") %>%
             dplyr::mutate(gene_id = rowData(dds_S1[tx_name])$gene_id, .after="tx_name") %>%                                   
             dplyr::mutate(IFNg_induced_v2 = gene_name %in% IFNg_induced_v2$gene_name) %>%
             dplyr::mutate(PSM = gene_id %in% PSM$gene_id)
}

.scatter_plot_TE <- function(res, label_x_pos=1, label_y_pos=-3.5, gap=0.5) {
  require(scales)
  require(ggrepel)
  msg_up   <- sprintf("log2(Ribo / RNA) > 1: %4.0f (%s)",  sum(res$status == "up"), percent(sum(res$status=="up")/nrow(res), accuracy=0.1))
  msg_down <- sprintf("log2(Ribo / RNA) < -1: %4.0f (%s)", sum(res$status == "down"), percent(sum(res$status=="down")/nrow(res), accuracy=0.1))  
  pearson <- cor(res$rna_lfc, res$ribo_lfc)             

  ggplot(res, aes(x=rna_lfc, y=ribo_lfc)) +
    geom_point(size=1.5, alpha=0.5, aes(color=status, shape=IFNg_induced_v2)) +
    theme_bw() +
    annotate("text", x=label_x_pos, y=label_y_pos, label=paste0("Pearson = ", format(pearson, digit=2)),
             hjust = 0, vjust=1) +
    annotate("text", x=label_x_pos, y=label_y_pos - gap, label=msg_up, color="red", hjust = 0, vjust=1) +
    annotate("text", x=label_x_pos, y=label_y_pos - 2*gap, label=msg_down, color="blue", hjust=0, vjust=1) +
    scale_color_manual(values=c("gray75", "red", "blue"), guide=FALSE) + #guide=FALSE
    scale_shape_manual(values=c(19, 1)) +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", legend.box="vertical") 
}

.summarize_feature <- function(tidy_feature) {
  up <- tidy_feature %>% dplyr::filter(status == "up")
  down <- tidy_feature %>% dplyr::filter(status == "down")
  c(n=nrow(tidy_feature), up=nrow(up), freq_up=nrow(up)/nrow(tidy_feature),
    down=nrow(down), freq_down = nrow(down)/nrow(tidy_feature))
}

#
# (A) S6: 5UTR + 1set coding exon
#
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR-plus-1st-exon")
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))




#
# (1) S6 - 1st exon by tx
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S6")
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx_mRNA.rda"))

dds_S6 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_1st_exon_by_tx_mRNA, rse_ribo_by_tx = rse_1st_exon_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S6"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S6_1stExon <- .do_seq_tidy_results_S6(dds_S6=dds_S6, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

# scatter plot                                    
.scatter_plot_TE(res=tidy_S6_1stExon, label_x_pos=0.8, label_y_pos=-3.8) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
      title="1st coding exon: Ribo-seq vs. RNA-seq (by tx)", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "S6_scatter_plot_1st_coding_exon_tx.pdf"), width=6, height=5)


# GO analysis
tidy_S6_1stExon_down <- tidy_S6_1stExon %>% dplyr::filter(status=="down")
tidy_S6_1stExon_up <- tidy_S6_1stExon %>% dplyr::filter(status=="up") # 651 transcripts / 325 genes
tidy_S6_1stExon_down %>% distinct(gene_name) %>% nrow()
tidy_S6_1stExon_up %>% distinct(gene_name) %>% nrow()

library(goseq)
source(file.path(pkg_dir, "scripts", "tools.R"))
universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")

# down
selected <- tidy_S6_1stExon_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"
# up                         
selected <- tidy_S6_1stExon_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"

write_xlsx(list(`all_1st_exon` = tidy_S6_1stExon,
                `translation_upregulation` = tidy_S6_1stExon_up, #
                `Enriched_GO_for_up` = enriched_go_up,
                `translation_downregulation` = tidy_S6_1stExon_down,
                `Enriched_GO_for_down` = enriched_go_down), 
           path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"))  

#
# (2) S6: Correct TSS 
#
load(file.path(pkg_dir, "data", "rse_TSS_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_TSS_by_tx.rda")) # by tx
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")

dds_S6 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_TSS_by_tx_mRNA, rse_ribo_by_tx = rse_TSS_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S6"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S6_TSS <- .do_seq_tidy_results_S6(dds_S6=dds_S6, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

# scatter plot                                    
.scatter_plot_TE(res=tidy_S6_TSS, label_x_pos=0.8, label_y_pos=-3.8) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
      title="TSS [-13, 13]: Ribo-seq vs. RNA-seq (by tx)", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "S6_scatter_plot_TSS_tx.pdf"), width=6, height=5)

#` GO term
tidy_S6_TSS_down <- tidy_S6_TSS %>% dplyr::filter(status=="down")
tidy_S6_TSS_up <- tidy_S6_TSS %>% dplyr::filter(status=="up")
tidy_S6_TSS_down %>% distinct(gene_name) %>% nrow()
tidy_S6_TSS_up %>% distinct(gene_name) %>% nrow()

universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
# down
selected <- tidy_S6_TSS_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"
# up                         
selected <- tidy_S6_TSS_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"

write_xlsx(list(`all_TSS` = tidy_S6_TSS,
                `translation_upregulation` = tidy_S6_TSS_up, 
                `Enriched_GO_for_up` = enriched_go_up,
                `translation_downregulation` = tidy_S6_TSS_down,
                `Enriched_GO_for_down` = enriched_go_down), 
           path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx")) 

#
# 5'UTR
#           
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")

dds_S6 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_5UTR_by_tx_mRNA, rse_ribo_by_tx = rse_5UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S6"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S6_5UTR <- .do_seq_tidy_results_S6(dds_S6=dds_S6, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

# scatte plot
.scatter_plot_TE(res=tidy_S6_5UTR, label_x_pos=1, label_y_pos=-3.5) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
      title="5' UTR: Ribo-seq vs. RNA-seq (by tx)", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "S6_scatter_plot_5UTR_tx.pdf"), width=6, height=5)

# GO
tidy_S6_5UTR_down <- tidy_S6_5UTR %>% dplyr::filter(status=="down")
tidy_S6_5UTR_down %>% distinct(gene_name) %>% nrow() # 364
tidy_S6_5UTR_up <- tidy_S6_5UTR %>% dplyr::filter(status=="up")
tidy_S6_5UTR_up %>% distinct(gene_name) %>% nrow() 

universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
selected <- tidy_S6_5UTR_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"
selected <- tidy_S6_5UTR_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"

write_xlsx(list(`all_5UTR` = tidy_S6_5UTR,
                `translation_upregulation` = tidy_S6_5UTR_up,
                `Enriched_Go_for_up` = enriched_go_up,
                `translation_downregulation` = tidy_S6_5UTR_down,
                `Enriched_GO_for_down` = enriched_go_down), 
           path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"))                          



#
# 3' UTR
#
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx.rda"))
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "3UTR_v2")

dds_S6 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_3UTR_by_tx_mRNA, rse_ribo_by_tx = rse_3UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S6"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S6_3UTR <- .do_seq_tidy_results_S6(dds_S6=dds_S6, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM)                                            
# scatte plot
.scatter_plot_TE(res=tidy_S6_3UTR, label_x_pos=1, label_y_pos=-3.5) +
 labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
      title="3' UTR: Ribo-seq vs. RNA-seq (by tx)", shape="IFNg-induced")
ggsave(file.path(fig_dir, "S6_scatter_plot_3UTR_tx.pdf"), width=6, height=5)


# GO
tidy_S6_3UTR_down <- tidy_S6_3UTR %>% dplyr::filter(status=="down")
tidy_S6_3UTR_up <- tidy_S6_3UTR %>% dplyr::filter(status=="up")
universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
# down
selected <- tidy_S6_3UTR_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"
# up                         
selected <- tidy_S6_3UTR_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) # dds must be "by gene"
write_xlsx(list(`all_3UTR` = tidy_S6_3UTR,
                `translation_upregulation` = tidy_S6_3UTR_up,
                `Enriched_Go_for_up` = enriched_go_up,
                `translation_downregulation` = tidy_S6_3UTR_down,
                `Enriched_GO_for_down` = enriched_go_down), 
           path=file.path(stats_dir, "S6-3UTR-by_tx-exclude-DUX4-indueced.xlsx"))       
         
#
# save the statistics
#
tidy_TE_S6_by_tx <- list(`5UTR`=tidy_S6_5UTR, `TSS`=tidy_S6_TSS, `1stCodingExon`=tidy_S6_1stExon, `3UTR`=tidy_S6_3UTR)
save(tidy_TE_S6_by_tx, file=file.path(pkg_dir, "data", "tidy_TE_S6_by_tx.rda"))

#
# boxplot of TE logFC for all genomic features
#
load(file.path(pkg_dir, "data", "tidy_TE_CDS.rda"))
fig_dir <- fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S6")
tidy_S6_CDS <- tidy_TE_CDS$S6 %>%
  rename(gene_id = "ensembl")
my_color <- wesanderson::wes_palette("Darjeeling1", n=5)
tmp <- 
  bind_rows(dplyr::select(tidy_S6_5UTR, gene_id, log2FoldChange),
            dplyr::select(tidy_S6_TSS, gene_id, log2FoldChange),
            dplyr::select(tidy_S6_1stExon, gene_id, log2FoldChange),
            dplyr::select(tidy_S6_CDS, gene_id, log2FoldChange), 
            dplyr::select(tidy_S6_3UTR, gene_id, log2FoldChange), .id="feature") %>%
  dplyr::mutate(feature=case_when(
        feature == 1 ~ "5' UTR",
        feature == 2 ~ "TSS [-13, 13]",
        feature == 3 ~ "1st coding exon",
        feature == 4 ~ "CDS", 
        feature == 5 ~ "3' UTR")) %>%
  dplyr::mutate(feature=factor(feature, levels=c("5' UTR", "TSS [-13, 13]", "1st coding exon", "CDS", "3' UTR")))
tmp_S6 <- tmp

ggplot(tmp, aes(x=feature, y=log2FoldChange, fill=feature)) +
  geom_boxplot(outlier.shape=NA, ) +
  theme_bw() +
  scale_fill_manual(values=my_color) +
  labs(y="Translation changes (log2FC)", title="S6") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  coord_cartesian(ylim=c(-4, 3))        
ggsave(file.path(fig_dir, "boxplot_all_genomic_features.pdf"), width=2, height=4)  

# overall summary (numbers)
.summarize_feature(tidy_S6_5UTR)
.summarize_feature(tidy_S6_TSS)
.summarize_feature(tidy_S6_1stExon)
.summarize_feature(tidy_TE_CDS$S6)
.summarize_feature(tidy_S6_3UTR)

#
# S2: IFNg vs untreated
#

# 1st coding exons
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S2")
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx_mRNA.rda"))
dds_S2 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_1st_exon_by_tx_mRNA, rse_ribo_by_tx = rse_1st_exon_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S2"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
exclude_gene_S2 = c(exclude_gene, IFNg_induced_v2$ensembl)                                           
tidy_S2_1stExon <- .do_seq_tidy_results_S2(dds_S2=dds_S2, exclude_gene=exclude_gene_S2, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

pdf(file.path(fig_dir, "S2_scatter_plot_1st_coding_exon_tx.pdf"), width=6, height=5)
gg <-  .scatter_plot_TE(res=tidy_S2_1stExon, label_x_pos=-1, label_y_pos=-3.5, gap=0.7)                                                                         
gg + labs(x=bquote(~"RNA-seq:" ~log[2]~"(IFNg / untreated)"), 
         y=bquote(~"Ribo-seq:" ~log[2]~"(FNg / untreated)"), 
         title="1st coding exon: Ribo-seq vs. RNA-seq") 
dev.off()

# (2) TSS
load(file.path(pkg_dir, "data", "rse_TSS_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_TSS_by_tx.rda")) # by tx
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")

dds_S2 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_TSS_by_tx_mRNA, rse_ribo_by_tx = rse_TSS_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S2"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S2_TSS <- .do_seq_tidy_results_S2(dds_S2=dds_S2, exclude_gene=exclude_gene_S2, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

pdf(file.path(fig_dir, "S2_scatter_plot_TSS_tx.pdf"), width=6, height=5)
gg <-  .scatter_plot_TE(res=tidy_S2_TSS, label_x_pos=-4, label_y_pos=-3.5, gap=0.7)                                                                         
gg + labs(x=bquote(~"RNA-seq:" ~log[2]~"(IFNg / untreated)"), 
         y=bquote(~"Ribo-seq:" ~log[2]~"(FNg / untreated)"), 
         title="TSS: Ribo-seq vs. RNA-seq") 
dev.off()                                   

# (3) 5' UTR
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))
dds_S2 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_5UTR_by_tx_mRNA, rse_ribo_by_tx = rse_5UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S2"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S2_5UTR <- .do_seq_tidy_results_S2(dds_S2=dds_S2, exclude_gene=exclude_gene_S2, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

pdf(file.path(fig_dir, "S2_scatter_plot_5UTR_tx.pdf"), width=6, height=5)
gg <-  .scatter_plot_TE(res=tidy_S2_5UTR, label_x_pos=2, label_y_pos=-4, gap=0.7)                                                                         
gg + labs(x=bquote(~"RNA-seq:" ~log[2]~"(IFNg / untreated)"), 
         y=bquote(~"Ribo-seq:" ~log[2]~"(FNg / untreated)"), 
         title="5' UTR: Ribo-seq vs. RNA-seq") 
dev.off()                 

# (4) 3' UTR
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx.rda"))
dds_S2 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_3UTR_by_tx_mRNA, rse_ribo_by_tx = rse_3UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S2"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S2_3UTR <- .do_seq_tidy_results_S2(dds_S2=dds_S2, exclude_gene=exclude_gene_S2, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM)  

pdf(file.path(fig_dir, "S2_scatter_plot_3UTR_tx.pdf"), width=6, height=5)
gg <-  .scatter_plot_TE(res=tidy_S2_3UTR, label_x_pos=0.5,  label_y_pos=-4.5, gap=0.7)                                                                         
gg + labs(x=bquote(~"RNA-seq:" ~log[2]~"(IFNg / untreated)"), 
         y=bquote(~"Ribo-seq:" ~log[2]~"(FNg / untreated)"), 
         title="3' UTR: Ribo-seq vs. RNA-seq") 
dev.off()   

# GO
tidy_S2_3UTR_down <- tidy_S2_3UTR %>% dplyr::filter(status=="down")
tidy_S2_3UTR_up <- tidy_S2_3UTR %>% dplyr::filter(status=="up")
universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
# down
selected <- tidy_S2_3UTR_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=FALSE, dds=dds_cds_by_gene) # dds must be "by gene"
# up                         
selected <- tidy_S2_3UTR_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=FALSE, dds=dds_cds_by_gene) # dds must be "by gene"

# (5) visualization of the translation chagnes
load(file.path(pkg_dir, "data", "tidy_TE_CDS.rda"))
tidy_S2_CDS <- tidy_TE_CDS$S2 %>%
  rename(gene_id = "ensembl")
#tidy_S2_CDS <- tidy_S2  %>% rename(gene_id = "ensembl")
my_color <- wesanderson::wes_palette("Darjeeling1", n=5)
tmp <- 
  bind_rows(dplyr::select(tidy_S2_5UTR, gene_id, log2FoldChange),
            dplyr::select(tidy_S2_TSS, gene_id, log2FoldChange),
            dplyr::select(tidy_S2_1stExon, gene_id, log2FoldChange),
            dplyr::select(tidy_S2_CDS, gene_id, log2FoldChange), 
            dplyr::select(tidy_S2_3UTR, gene_id, log2FoldChange), .id="feature") %>%
  dplyr::mutate(feature=case_when(
        feature == 1 ~ "5' UTR",
        feature == 2 ~ "TSS [-13, 13]",
        feature == 3 ~ "1st coding exon",
        feature == 4 ~ "CDS", 
        feature == 5 ~ "3' UTR")) %>%
  dplyr::mutate(feature=factor(feature, levels=c("5' UTR", "TSS [-13, 13]", "1st coding exon", "CDS", "3' UTR")))

tmp_S2 <- tmp

my_color <- wesanderson::wes_palette("Darjeeling1", n=5)[1:5]
ggplot(tmp_S2, aes(x=feature, y=log2FoldChange, fill=feature)) +
  #geom_hline(yintercept=hline, size=1, color="#5B1A18", alpha=0.3) +
  geom_boxplot(outlier.shape=NA) +
  theme_bw() +
  scale_fill_manual(values=my_color) +
  labs(y="Translation changes (log2FC)", title="S2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        #legend.position = "none", panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  coord_cartesian(ylim=c(-4, 3))        
ggsave(file.path(fig_dir, "boxplot_all_genomic_features_2.pdf"), width=2, height=4)  

# summary
.summarize_feature(tidy_S2_5UTR)
.summarize_feature(tidy_S2_TSS)
.summarize_feature(tidy_S2_1stExon)
.summarize_feature(tidy_TE_CDS$S2)
.summarize_feature(tidy_S2_3UTR)

# output statistics: 5'UTR, TSS, 
.output_tidy_res_go <- function(tidy_res, universe_dds, output_path, return.DEInCat=FALSE) {
  tidy_down <- tidy_res %>% dplyr::filter(status=="down")
  tidy_up <- tidy_res %>% dplyr::filter(status=="up")
  universal <- rownames(universe_dds) %>% str_replace("\\..*", "")
  # down
  selected <- tidy_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
  enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                            return.DEInCat=return.DEInCat, dds=universe_dds) # dds must be "by gene"
  # up                         
  selected <- tidy_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
  enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                           return.DEInCat=return.DEInCat, dds=universe_dds) # dds must be "by gene"
  writexl::write_xlsx(list(`all` = tidy_res,
                  `translation_upregulation` = tidy_up,
                  `Enriched_Go_for_up` = enriched_go_up,
                  `translation_downregulation` = tidy_down,
                  `Enriched_GO_for_down` = enriched_go_down), 
             path=output_path)  
}

.output_tidy_res_go(tidy_res=tidy_S2_1stExon, universe_dds=dds_cds_by_gene, 
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "1st_exon_v2", "S2-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"))
.output_tidy_res_go(tidy_res=tidy_S2_TSS, universe_dds=dds_cds_by_gene, 
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "around_TSS_v2", "S2-TSS-by_tx-exclude-DUX4-indueced.xlsx"))                                          
.output_tidy_res_go(tidy_res=tidy_S2_5UTR, universe_dds=dds_cds_by_gene, 
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "5UTR_v2", "S2-5UTR-by_tx-exclude-DUX4-indueced.xlsx")) 
.output_tidy_res_go(tidy_res=tidy_S2_3UTR, universe_dds=dds_cds_by_gene, 
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "3UTR_v2", "S2-3UTR-by_tx-exclude-DUX4-indueced.xlsx"))                                          

#
# S1: DOX vs untreated
#


# (1) S1 - 1st coding exon
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S1")
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx_mRNA.rda"))
dds_S1 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_1st_exon_by_tx_mRNA, rse_ribo_by_tx = rse_1st_exon_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S1"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S1_1stExon <- .do_seq_tidy_results_S1(dds_S1=dds_S1, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 



# (2) TSS
load(file.path(pkg_dir, "data", "rse_TSS_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_TSS_by_tx.rda")) # by tx
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")

dds_S1 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_TSS_by_tx_mRNA, rse_ribo_by_tx = rse_TSS_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S1"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S1_TSS <- .do_seq_tidy_results_S1(dds_S1=dds_S1, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

# (3) 5' UTR
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))

dds_S1 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_5UTR_by_tx_mRNA, rse_ribo_by_tx = rse_5UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S1"]],
                                           rna_mean_filter = 10, ribo_mean_filter = 5)
tidy_S1_5UTR <- .do_seq_tidy_results_S1(dds_S1=dds_S1, exclude_gene=exclude_gene, 
                                   DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM) 

# (4) 3' UTR
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx_mRNA.rda"))
load(file.path(pkg_dir, "data", "rse_3UTR_by_tx.rda"))
dds_S1 <- .comb_RNA_Ribo_sizeFactor_by_cds(rse_rna_by_tx = rse_3UTR_by_tx_mRNA, rse_ribo_by_tx = rse_3UTR_by_tx, 
                                           rse_rna_cds_by_gene = rse_cds_mRNA, rse_ribo_cds_by_gene = rse_cds_by_gene,
                                           treatments = list_comp[["S1"]],
                                           rna_mean_filter = 5, ribo_mean_filter = 3)
tidy_S1_3UTR <- .do_seq_tidy_results_S1(dds_S1=dds_S1, exclude_gene=exclude_gene, 
                                        DUX4_induced_v2=DUX4_induced_v2, IFNg_induced_v2=IFNg_induced_v2, PSM=PSM)  

# GO
tidy_S1_3UTR_down <- tidy_S1_3UTR %>% dplyr::filter(status=="down")
tidy_S1_3UTR_up <- tidy_S1_3UTR %>% dplyr::filter(status=="up")
universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
# down
selected <- tidy_S1_3UTR_down %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=FALSE, dds=dds_cds_by_gene) # dds must be "by gene"
# up                         
selected <- tidy_S1_3UTR_up %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.05,
                         return.DEInCat=FALSE, dds=dds_cds_by_gene) # dds must be "by gene"

# (5) visualization of the translation chagnes
load(file.path(pkg_dir, "data", "tidy_TE_CDS.rda"))
tidy_S1_CDS <- tidy_TE_CDS$S1 %>%
  rename(gene_id = "ensembl")
my_color <- wesanderson::wes_palette("Darjeeling1", n=5)
tmp <- 
  bind_rows(dplyr::select(tidy_S1_5UTR, gene_id, log2FoldChange),
            dplyr::select(tidy_S1_TSS, gene_id, log2FoldChange),
            dplyr::select(tidy_S1_1stExon, gene_id, log2FoldChange),
            dplyr::select(tidy_S1_CDS, gene_id, log2FoldChange), 
            dplyr::select(tidy_S1_3UTR, gene_id, log2FoldChange), .id="feature") %>%
  dplyr::mutate(feature=case_when(
        feature == 1 ~ "5' UTR",
        feature == 2 ~ "TSS [-13, 13]",
        feature == 3 ~ "1st coding exon",
        feature == 4 ~ "CDS", 
        feature == 5 ~ "3' UTR")) %>%
  dplyr::mutate(feature=factor(feature, levels=c("5' UTR", "TSS [-13, 13]", "1st coding exon", "CDS", "3' UTR")))

hline <- tmp %>% dplyr::filter(feature == "CDS") %>% summarize(median=median(log2FoldChange)) %>% pull(median)
tmp_S1 <- tmp

ggplot(tmp_S1, aes(x=feature, y=log2FoldChange)) +
  geom_hline(yintercept=hline, size=1, color="#5B1A18", alpha=0.3) +
  geom_boxplot(outlier.shape=NA, fill="gray75", alpha=0.5) +
  theme_minimal() +
  #scale_fill_manual(values=my_color) +
  labs(y="Translation changes (log2FC)", title="S1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "none", panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  coord_cartesian(ylim=c(-4, 3))        
ggsave(file.path(fig_dir, "test1.pdf"), width=2, height=4)
#ggsave(file.path(fig_dir, "boxplot_all_genomic_features.pdf"), width=2, height=4)  

.summarize_feature(tidy_S1_5UTR)
.summarize_feature(tidy_S1_TSS)
.summarize_feature(tidy_S1_1stExon)
.summarize_feature(tidy_TE_CDS$S1)
.summarize_feature(tidy_S1_3UTR)

.output_tidy_res_go(tidy_res=tidy_S1_1stExon, universe_dds=dds_cds_by_gene, return.DEInCat=TRUE,
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "1st_exon_v2", "S1-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"))
.output_tidy_res_go(tidy_res=tidy_S1_TSS, universe_dds=dds_cds_by_gene, return.DEInCat=TRUE,
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "around_TSS_v2", "S1-TSS-by_tx-exclude-DUX4-indueced.xlsx"))                                          
.output_tidy_res_go(tidy_res=tidy_S1_5UTR, universe_dds=dds_cds_by_gene, return.DEInCat=TRUE,
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "5UTR_v2", "S1-5UTR-by_tx-exclude-DUX4-indueced.xlsx")) 
.output_tidy_res_go(tidy_res=tidy_S1_3UTR, universe_dds=dds_cds_by_gene, return.DEInCat=TRUE,
                    output_path=file.path(pkg_dir, "stats", "translation_changes", 
                                          "3UTR_v2", "S1-3UTR-by_tx-exclude-DUX4-indueced.xlsx"))     


#
# major figure for manuscript
#             
tmp_S1 <- tmp_S1 %>% add_column(scenario="S1")     
hline_S1 <- tmp_S1 %>% dplyr::filter(feature == "CDS") %>% summarize(median=median(log2FoldChange)) %>% pull(median)          
tmp_S2 <- tmp_S2 %>% add_column(scenario="S2")
hline_S2 <- tmp_S2 %>% dplyr::filter(feature == "CDS") %>% summarize(median=median(log2FoldChange)) %>% pull(median)
tmp_S6 <- tmp_S6 %>% add_column(scenario="S6")
hline_S6 <- tmp_S6 %>% dplyr::filter(feature == "CDS") %>% summarize(median=median(log2FoldChange)) %>% pull(median)
tmp <- tmp_S1 %>% add_row(tmp_S2) %>% add_row(tmp_S6)
df_hline <- data.frame(scenario=c("S1", "S2", "S6"), grp_median=c(hline_S1, hline_S2, hline_S6))
data_for_S1_S2_S6_TE=tmp

save(data_for_S1_S2_S6_TE, file=file.path(pkg_dir, "data", "data_for_S1_S2_S6_TE.rda"))

# correction made on 5/13/2022
# load(file.path(pkg_dir, "data", "data_for_S1_S2_S6_TE.rda"))
 data_for_S1_S2_S6_TE = data_for_S1_S2_S6_TE %>% dplyr::filter(scenario %in% c("S1", "S6")) %>% 
  add_row(tmp_S2 %>% add_column(scenario="S2"))
# save(data_for_S1_S2_S6_TE, file=file.path(pkg_dir, "data", "data_for_S1_S2_S6_TE.rda"))


#
#  use data_for_S1_S2_S56_TE: boxplot 
#
load(file.path(pkg_dir, "data", "data_for_S1_S2_S6_TE.rda"))
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "translation-changes")
tmp <- data_for_S1_S2_S6_TE 
tmp <- tmp %>% dplyr::mutate(scenario=factor(scenario, levels=c("S2", "S1", "S6")))


# (1) viline plot
ggplot(tmp, aes(x=feature, y=log2FoldChange, fill=feature)) +
  geom_violin(trim=TRUE, scale="width") +
  #geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  theme_bw() +
  facet_wrap(~scenario, nrow=1) +
  labs(y="Translation changes (log2FC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "none", panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank()) +
  #coord_cartesian(ylim=c(-5, 5)) +
  stat_summary(fun=median, geom="point", aes(stoke=0.3)) +
  scale_color_manual(values= wes_palette("Darjeeling1", n = 5)) +
  scale_fill_manual(values=wes_palette("Darjeeling1", n = 5)) 
  #stat_summary(fun=median, geom="line", size=1, aes(group=1), color="red") +
  #stat_summary(fun=median, geom="point", aes(stoke=0.3)) 
#  stat_summary(fun.data="mean_sdl", mult=1, 
#                 geom="crossbar", width=0.2 )
ggsave(file.path(fig_dir, "violine-log2FC-translation-changes-all-features-S1-S2-S6-v2.pdf"), width=6, height=4)

# (2) boxplot 1
tmp <- tmp %>% dplyr::mutate(scenario=factor(scenario, levels=c("S2", "S1", "S6")))
ggplot(tmp, aes(x=feature, y=log2FoldChange)) +
  geom_hline(data=df_hline, aes(yintercept=grp_median), size=1, color="gray25", alpha=0.5) +
  geom_boxplot(outlier.shape=NA, fill="gray25", color="gray25", alpha=0.3) +
  theme_bw() +
  #scale_fill_manual(values=my_color) +
  facet_wrap(~scenario, nrow=1) +
  labs(y="Translation changes (log2FC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        #panel.background = element_rect(fill = "#5e6472"),
        #plot.background = element_rect(fill = "#5e6472"),
        legend.position = "none", panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  coord_cartesian(ylim=c(-4, 3)) +
  #scale_color_manual(values= wes_palette("Royal2", n = 5)) +
  #scale_fill_manual(values=wes_palette("Royal2", n = 5)) +
  stat_summary(fun=median, geom="line", size=1, aes(group=1), color="red") +
  stat_summary(fun=median, geom="point", aes(stroke=1))
ggsave(file.path(fig_dir, "barplot-log2FC-translation-changes-all-features-S1-S2-S6.pdf"), width=5.5, height=4)


# (2) boxplot 2
# how about put them all together
tmp <- tmp %>% dplyr::mutate(scenario=factor(scenario, levels=c("S2", "S1", "S6")))
royal <- wes_palette("Royal2", n = 5)[c(1, 3,5)]
ggplot(tmp, aes(x=feature, y=log2FoldChange, color=scenario, fill=scenario)) +
  #geom_hline(data=df_hline, aes(yintercept=grp_median), size=1, color="gray25", alpha=0.5) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  theme_bw() +
  #scale_fill_manual(values=my_color) +
  #facet_wrap(~scenario, nrow=1) +
  labs(y="Translation changes (log2FC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        #panel.background = element_rect(fill = "#5e6472"),
        #plot.background = element_rect(fill = "#5e6472"),
        legend.position = "none", panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  coord_cartesian(ylim=c(-4, 3)) +
  scale_color_manual(values = royal) +
  scale_fill_manual(values = royal) +
  # scale_color_manual(values= my_color[1:3]) +
  # scale_fill_manual(values=my_color[1:3]) +
  stat_summary(fun=median, geom="line", size=0.5, aes(group=scenario), position = position_dodge(width = 0.9))
  #stat_summary(fun=median, geom="point", aes(stroke=1, group=scenario))
ggsave(file.path(fig_dir, "barplot-log2FC-translation-changes-all-features-S1-S2-S6-v3.pdf"), width=4, height=4)

# (3) p-value using mann-whitney tests
load(file.path(pkg_dir, "data", "data_for_S1_S2_S6_TE.rda"))
tmp <- data_for_S1_S2_S6_TE 
.summarize_wilcox <- function(df) {
  df <- df %>% drop_na(log2FoldChange)
  lfc_5utr <- df %>% dplyr::filter(feature=="5' UTR") %>% pull(log2FoldChange)
  lfc_tss <- df %>% dplyr::filter(feature=="TSS [-13, 13]") %>% pull(log2FoldChange)
  lfc_1st <- df %>% dplyr::filter(feature=="1st coding exon") %>% pull(log2FoldChange)
  lfc_cds <- df %>% dplyr::filter(feature=="CDS") %>% pull(log2FoldChange)
  lfc_3utr <- df %>% dplyr::filter(feature=="3' UTR") %>% pull(log2FoldChange)
  cat("Wilcox test (a.k.a. two-sample Mann-Whitney), alternative=two-sided")
  c(utr5_vs_cds = wilcox.test(x=lfc_cds, y=lfc_5utr, alternaive="two.sided")$p.value, 
    tss_vs_cds = wilcox.test(x=lfc_cds, y=lfc_tss, alternative="two.sided")$p.value,
    firstExon_vs_cds = wilcox.test(x=lfc_cds, y=lfc_1st, alternative="two.sided")$p.value,
    utr3_vs_cds = wilcox.test(x=lfc_cds, y=lfc_3utr, alternative="two.sided")$p.value)
}

s1 <- tmp %>% dplyr::filter(scenario=="S1")
s2 <- tmp %>% dplyr::filter(scenario=="S2")
s6 <- tmp %>% dplyr::filter(scenario=="S6")
.summarize_wilcox(s1)
.summarize_wilcox(s2)
.summarize_wilcox(s6)

.summarize_anova <- function(df) {
  res <- anova
  df <- df %>% drop_na(log2FoldChange)
  lfc_5utr <- df %>% dplyr::filter(feature=="5' UTR") %>% pull(log2FoldChange)
  lfc_tss <- df %>% dplyr::filter(feature=="TSS [-13, 13]") %>% pull(log2FoldChange)
  lfc_1st <- df %>% dplyr::filter(feature=="1st coding exon") %>% pull(log2FoldChange)
  lfc_cds <- df %>% dplyr::filter(feature=="CDS") %>% pull(log2FoldChange)
  lfc_3utr <- df %>% dplyr::filter(feature=="3' UTR") %>% pull(log2FoldChange)
  cat("Wilcox test (a.k.a. two-sample Mann-Whitney), alternative=two-sided \n")
  c(utr5_vs_cds = t.test(x=lfc_cds, y=lfc_5utr, alternaive="two.sided")$p.value, 
    tss_vs_cds = t.test(x=lfc_cds, y=lfc_tss, alternative="two.sided")$p.value,
    firstExon_vs_cds = t.test(x=lfc_cds, y=lfc_1st, alternative="two.sided")$p.value,
    utr3_vs_cds = t.test(x=lfc_cds, y=lfc_3utr, alternative="two.sided")$p.value)
}

# ok, let's plot the density of CDS and 3'UTR
s6 %>% dplyr::filter(feature %in% c("CDS", "3' UTR")) %>%
  ggplot(aes(x=log2FoldChange, group=feature, fill=feature)) +
    geom_density(alpha=0.3) +
    theme_bw()
ggsave(file.path(pkg_dir, "test.pdf"), width=4, height=2.5)    

aov <- lapply(levels(s6$feature)[-4], function(x) {
  data <- s6 %>% dplyr::filter(feature %in% c(x, "CDS"))
  anova(lm(log2FoldChange ~ feature, data=data))
})
names(aov) = levels(s6$feature)[-4]


aov2 <- lapply(levels(s2$feature)[-4], function(x) {
  data <- s2 %>% dplyr::filter(feature %in% c(x, "CDS"))
  anova(lm(log2FoldChange ~ feature, data=data))
})
names(aov2) = levels(s2$feature)[-4]

aov1 <- lapply(levels(s1$feature)[-4], function(x) {
  data <- s1 %>% dplyr::filter(feature %in% c(x, "CDS"))
  anova(lm(log2FoldChange ~ feature, data=data))
})
names(aov1) = levels(s2$feature)[-4]