# This script uses DESeq2 to determine the treatment effects on translation efficicancy.
# The iteraction team can be formulated in two ways, and both yield the same consequences:
#
# (1) Interaction term interprets whether the change in Ribo-seq are due to RNA-seq - 
#     estimate the difference between the treatment effects in Ribo-seq and in RNA-seq - 
#     𝑙𝑜g〖𝑡𝑟𝑒𝑎𝑡𝑒𝑑/𝑢𝑛𝑡𝑟𝑒𝑎𝑡𝑒𝑑|〗_(𝑅𝑖𝑏𝑜−𝑠𝑒𝑞)  − 𝑙𝑜𝑔 〖𝑡𝑟𝑒𝑎𝑡𝑒𝑑/𝑢𝑛𝑡𝑟𝑒𝑎𝑡𝑒𝑑|〗_(𝑅𝑁𝐴−𝑠𝑒𝑞)
#     ~ 𝑝𝑟𝑜𝑡𝑜𝑐𝑜𝑙+𝑡𝑟𝑒𝑎𝑡𝑚𝑒𝑛𝑡+𝑝𝑟𝑜𝑡𝑜𝑐𝑜𝑙:𝑡𝑟𝑒𝑎𝑡𝑚𝑒𝑛𝑡
# (2) conventional definition of translation efficiency - estimate the differences of 
#     the translational efficiency between treatment conditions - 
#     {𝑙𝑜𝑔〖(𝑅𝑖𝑏𝑜−𝑠𝑒𝑞)/(𝑅𝑁𝐴−𝑠𝑒𝑞)|〗_𝑡𝑟𝑒𝑎𝑡𝑒𝑑  −𝑙𝑜𝑔〖(𝑅𝑖𝑏𝑜−𝑠𝑒𝑞)/(𝑅𝑁𝐴−𝑠𝑒𝑞)|〗_𝑢𝑛𝑡𝑟𝑒𝑎𝑡𝑒𝑑 } 
#     ~ 𝑡𝑟𝑒𝑎𝑡𝑚𝑒𝑛𝑡+𝑝𝑟𝑜𝑡𝑜𝑐𝑜𝑙+𝑡𝑟𝑒𝑎𝑡𝑚𝑒𝑛𝑡:𝑝𝑟𝑜𝑡𝑜𝑐𝑜𝑙
#
# Specs:
# (1) combine rse_cds_by_gene and dds_cds_mRNA and subset by treatments of interest
# (2) do sizeFactor and dispersion separately
# (3) design = ~ treatment + protocol + treatment:protocol or protocol + treatment + protocol:treatment
# (4) get the results of the interaction term

library(DESeq2)
library(readxl)
library(writexl)
library(tidyverse)
library(corrr)
library(plyranges)
library(wesanderson)
library(goseq)
library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35

library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "CDS")
source(file.path(pkg_dir, "scripts", "06A-tools_TE.R"))
source(file.path(pkg_dir, "scripts", "tools.R"))

my_color <- wesanderson::wes_palette("Darjeeling1", n=5)[2:5]
names(my_color) <- levels(rse_cds_by_gene$treatment)

list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))

#
# load rse by genes
#
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda")) # cds by genes

#
# load DUX4_altered and INFg_altered
#
load(file.path(pkg_dir, "data", "DUX4_induced.rda")) #LFC > 2
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
.scatter_plot_TE <- function(res, label_x_pos=3, label_y_pos=-2, gap=0.5) {
  msg_up   <- sprintf("log2(Ribo / RNA) > 1: %4.0f", sum(res$status == "up"))
  msg_down <- sprintf("log2(Ribo / RNA) < -1: %4.0f", sum(res$status == "down"))  
  pearson <- cor(res$rna_lfc, res$ribo_lfc)             

  ggplot(res, aes(x=rna_lfc, y=ribo_lfc)) +
    geom_point(size=1.5, alpha=0.5, aes(color=status, shape=IFNg_induced)) +
    theme_bw() +
    #geom_vline(xintercept=c(-2, 2), linetype="dashed", alpha=0.5, color="gray50") +
    #geom_hline(yintercept=c(-2, 2), linetype="dashed", alpha=0.5, color="gray50") +
    annotate("text", x=label_x_pos, y=label_y_pos, label=paste0("Pearson = ", format(pearson, digit=2)),
             hjust = 0, vjust=1) +
    annotate("text", x=label_x_pos, y=label_y_pos-gap, label=msg_up, color="red", hjust = 0, vjust=1) +
    annotate("text", x=label_x_pos, y=label_y_pos-2*gap, label=msg_down, color="blue", hjust=0, vjust=1) +
    scale_color_manual(values=c("gray75", "red", "blue"), guide=FALSE) + #guide=FALSE
    scale_shape_manual(values=c(19, 1)) +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", legend.box="vertical")
}

#
# (1) S1 interaction term: null hypothesis: |LFC| = 0
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S1")
dds_S1 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, treatments=list_comp[[1]],
                         rna_mean_filter=10, ribo_mean_filter=5)
dds_S1 <- dds_S1[!rownames(dds_S1) %in% exclude_gene]
dds_S1 <- DESeq(dds_S1)
S1_ribo_over_rna <- results(dds_S1, name="protocolRibo.treatmentDOX_pulse", alpha=0.05) # translation efficiency changes
S1_rna <- results(dds_S1, name="treatment_DOX_pulse_vs_untreated", alpha=0.05)  
S1_ribo <- results(dds_S1, contrast=list(c("treatment_DOX_pulse_vs_untreated", "protocolRibo.treatmentDOX_pulse")), alpha=0.05)
tidy_S1 <- .tidy_res_ribo_over_rna(dds=dds_S1, inter_res=S1_ribo_over_rna, 
                                   ribo_res=S1_ribo, rna_res=S1_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((untreated_RNA_avg + DOX_pulse_RNA_avg)/2 >= 75) %>%
             dplyr::filter((untreated_Ribo_avg + DOX_pulse_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after=status)
tidy_S1_down <- tidy_S1 %>% dplyr::filter(status=="down")
tidy_S1_up <- tidy_S1 %>% dplyr::filter(status=="up")

# go analysis 
source(file.path(pkg_dir, "scripts", "tools.R"))
universal <- tidy_S1 %>% pull(ensembl) %>% str_replace("\\..*", "")
selected <- tidy_S1_down %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_genes=selected, 
                         return.DEInCat=TRUE, dds=dds_S1, p_value = 0.05)
selected <- tidy_S1_up %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_genes=selected, 
                         return.DEInCat=TRUE, dds=dds_S1, p_value = 0.05)
write_xlsx(list(`CDS` = tidy_S1,
                `translation_upregulation` = tidy_S1_up,
                `enriched_GO_term_up` = enriched_go_up,
                `translation_downregulation` = tidy_S1_down,
                `enriched_GO_term_down` = enriched_go_down), 
           path=file.path(stats_dir, "S1_DOX_pulse_vs_untreated.xlsx")) 

#How many translation down-regulated are associated with these GO term
tmp_S1 <- enriched_go_down %>%
  pull(DEInCat) %>% map(str_split, ",") %>% unlist() %>% unique()
length(tmp_S1)

# volcano
.plot_volcano(tidy_S1, title="DOX_pulse / untreated")  
ggsave(file.path(fig_dir, "volcano_DOX_pulse_vs_untreated_TEchanges.pdf"), width=4.5, height=3)

# scater plot: Ribo vs RNA
gg <- .scatter_plot_TE(res=tidy_S1, label_x_pos=-1, label_y_pos=-3) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DUX4 / untreated)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DUX4 / untreated)"), 
       title="S1: Ribo-seq vs. RNA-seq", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "TE_scatter_Ribo_vs_RNA_on_S1.pdf"), width=6, height=5)       

#
# (3) S2 - IFNg vs. untereated interaction term: null hypothesis: |LFC| = 0
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S2")
dds_S2 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, treatments=list_comp[["S2"]],
                         rna_mean_filter = 10, ribo_mean_filter = 5)
dds_S2 <- dds_S2[!rownames(dds_S2) %in% c(exclude_gene, IFNg_induced_v2$ensembl)]                         
dds_S2 <- DESeq(dds_S2)
S2_ribo_over_rna <- results(dds_S2, name="protocolRibo.treatmentIFNg", alpha=0.05) # translation efficiency changes
S2_rna <- results(dds_S2, name="treatment_IFNg_vs_untreated", alpha=0.05)  
S2_ribo <- results(dds_S2, contrast=list(c("treatment_IFNg_vs_untreated", "protocolRibo.treatmentIFNg")), alpha=0.05)
tidy_S2 <- .tidy_res_ribo_over_rna(dds=dds_S2, inter_res=S2_ribo_over_rna, 
                                   ribo_res=S2_ribo, rna_res=S2_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((untreated_RNA_avg + IFNg_RNA_avg)/2 >= 75) %>%
             dplyr::filter((untreated_Ribo_avg + IFNg_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after = status)
tidy_S2_down <- tidy_S2 %>% dplyr::filter(status=="down")
tidy_S2_up <- tidy_S2 %>% dplyr::filter(status=="up")

# volcano
.plot_volcano(tidy_S2, title="IFNg / untreated", numerator="IFNg")  
ggsave(file.path(fig_dir, "volcano_IFNg_vs_untreated_TEchanges.pdf"), width=4.5, height=3)

# scater plot: Ribo vs RNA
gg <- .scatter_plot_TE(res=tidy_S2, label_x_pos=-1, label_y_pos=-3.8, gap=0.7) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(IFNg / untreated)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(IFNg / untreated)"), 
       title="S2: Ribo-seq vs. RNA-seq", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "TE_scatter_Ribo_vs_RNA_on_S2.pdf"), width=6, height=5)     

# up-regulated translation changes
a <- tidy_S2 %>% dplyr::filter(rna_lfc >1, ribo_lfc >1, status=="up")
interested <- tidy_S2 %>% dplyr::filter(rna_lfc <1, ribo_lfc > -1, status=="up") 

# GO term for 82 up-regulated (interested)
universal <- tidy_S2 %>% pull(ensembl) %>% str_replace("\\..*", "")
selected <- interested %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go <- .do_goseq(universal=universal, selected_genes=selected, 
                         return.DEInCat=TRUE, dds=dds_S2, p_value = 0.1)

write_xlsx(list(`CDS`=tidy_S2, `translation_upregulation` = tidy_S2_up, 
                `translation_downregulation` = tidy_S2_down,
                `enriched_GO_term_up (FDR 0.1)` = enriched_go), 
           path=file.path(stats_dir, "S2_IFNg_vs_untreated.xlsx")) 

#
# (4) S3 - DOX_pulse_INFg vs. untereated interaction term: null hypothesis: |LFC| = 0
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S3")
dds_S3 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, treatments=list_comp[["S3"]],
                         rna_mean_filter = 10, ribo_mean_filter = 5)
dds_S3 <- dds_S3[!rownames(dds_S3) %in% exclude_gene]                         
dds_S3 <- DESeq(dds_S3)
S3_ribo_over_rna <- results(dds_S3, name="protocolRibo.treatmentDOX_pulse_IFNg", alpha=0.05) # translation efficiency changes
S3_rna <- results(dds_S3, name="treatment_DOX_pulse_IFNg_vs_untreated", alpha=0.05)  
S3_ribo <- results(dds_S3, contrast=list(c("treatment_DOX_pulse_IFNg_vs_untreated", "protocolRibo.treatmentDOX_pulse_IFNg")), alpha=0.05)
tidy_S3 <- .tidy_res_ribo_over_rna(dds=dds_S3, inter_res=S3_ribo_over_rna, 
                                   ribo_res=S3_ribo, rna_res=S3_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((untreated_RNA_avg + DOX_pulse_IFNg_RNA_avg)/2 >= 75) %>%
             dplyr::filter((untreated_Ribo_avg + DOX_pulse_IFNg_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after=status)
tidy_S3_down <- tidy_S3 %>% dplyr::filter(status=="down")
tidy_S3_up <- tidy_S3 %>% dplyr::filter(status=="up")
write_xlsx(tidy_S3, path=file.path(stats_dir, "S3_DOX_pulse_IFNg_vs_untreated.xlsx"))

#
# (4) S4 - DOX_pulse vs. IFNg 
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S4")
dds_S4 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, treatments=list_comp[["S4"]],
                         rna_mean_filter = 10, ribo_mean_filter = 5)
dds_S4 <- dds_S4[!rownames(dds_S4) %in% exclude_gene]    
dds_S4 <- DESeq(dds_S4)
S4_ribo_over_rna <- results(dds_S4, name="protocolRibo.treatmentDOX_pulse", alpha=0.05) # translation efficiency changes
S4_rna <- results(dds_S4, name="treatment_DOX_pulse_vs_IFNg", alpha=0.05)  
S4_ribo <- results(dds_S4, contrast=list(c("treatment_DOX_pulse_vs_IFNg", "protocolRibo.treatmentDOX_pulse")), alpha=0.05)
tidy_S4 <- .tidy_res_ribo_over_rna(dds=dds_S4, inter_res=S4_ribo_over_rna, 
                                   ribo_res=S4_ribo, rna_res=S4_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((IFNg_RNA_avg + DOX_pulse_RNA_avg)/2 >= 75) %>%
             dplyr::filter((IFNg_Ribo_avg + DOX_pulse_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after=status)
tidy_S4_down <- tidy_S4 %>% dplyr::filter(status=="down")
tidy_S4_up<- tidy_S4 %>% dplyr::filter(status=="up")

cor(tidy_S4$rna_lfc, tidy_S4$ribo_lfc) # 
write_xlsx(tidy_S4, path=file.path(stats_dir, "S4_DOX_pulse_vs_IFNg.xlsx"))

#
# (5) S5 - DOX_pulse_IFNg vs DOX_pulse 
#

fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S5")
dds_S5 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, treatments=list_comp[["S5"]],
                         rna_mean_filter = 10, ribo_mean_filter = 5)
dds_S5 <- dds_S5[!rownames(dds_S5) %in% exclude_gene]    
dds_S5 <- DESeq(dds_S5)
resultsNames(dds_S5) # sanity check
S5_ribo_over_rna <- results(dds_S5, name="protocolRibo.treatmentDOX_pulse_IFNg", alpha=0.05) # translation efficiency changes
S5_rna <- results(dds_S5, name="treatment_DOX_pulse_IFNg_vs_DOX_pulse", alpha=0.05)  
S5_ribo <- results(dds_S5, contrast=list(c("treatment_DOX_pulse_IFNg_vs_DOX_pulse", "protocolRibo.treatmentDOX_pulse_IFNg")), alpha=0.05)
tidy_S5 <- .tidy_res_ribo_over_rna(dds=dds_S5, inter_res=S5_ribo_over_rna, 
                                   ribo_res=S5_ribo, rna_res=S5_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((DOX_pulse_IFNg_RNA_avg + DOX_pulse_RNA_avg)/2 >= 75) %>%
             dplyr::filter((DOX_pulse_IFNg_Ribo_avg + DOX_pulse_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after=status) 
table(tidy_S5$status)
cor(tidy_S5$rna_lfc, tidy_S5$ribo_lfc) 
write_xlsx(tidy_S5, path=file.path(stats_dir, "S5_DOX_pulse_IFNg_vs_DOX_pulse.xlsx"))

#
# (6) S6 - DOX_pulse_INFg vs. INFg interaction term: null hypothesis: |LFC| = 0
#
fig_dir <- file.path(pkg_dir, "figures", "translation_changes", "S6")
dds_S6 <- .comb_RNA_Ribo(rse_rna = rse_cds_mRNA, rse_ribo = rse_cds_by_gene, 
                         treatments=list_comp[["S6"]],
                         rna_mean_filter=10, ribo_mean_filter=5)
dds_S6 <- dds_S6[!rownames(dds_S6) %in% exclude_gene]        
dds_S6 <- DESeq(dds_S6)
resultsNames(dds_S6)
S6_ribo_over_rna <- results(dds_S6, name="protocolRibo.treatmentDOX_pulse_IFNg", alpha=0.05) # translation efficiency changes
S6_rna <- results(dds_S6, name="treatment_DOX_pulse_IFNg_vs_IFNg", alpha=0.05)  
S6_ribo <- results(dds_S6, contrast=list(c("treatment_DOX_pulse_IFNg_vs_IFNg", "protocolRibo.treatmentDOX_pulse_IFNg")), alpha=0.05)

# tidy the results; add flag and DUX4_induced and IFNg_induced
tidy_S6 <- .tidy_res_ribo_over_rna(dds=dds_S6, inter_res=S6_ribo_over_rna, 
                                   ribo_res=S6_ribo, rna_res=S6_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter((IFNg_RNA_avg + DOX_pulse_IFNg_RNA_avg)/2 >= 75) %>%
             dplyr::filter((IFNg_Ribo_avg + DOX_pulse_IFNg_Ribo_avg)/2 >= 15) %>%
             dplyr::mutate(IFNg_induced = ensembl %in% IFNg_induced_v2$ensembl, .after="status")

tidy_S6_down <- tidy_S6 %>% dplyr::filter(status=="down")
tidy_S6_up <- tidy_S6 %>% dplyr::filter(status=="up")

# volcano
.plot_volcano(tidy_S6, title="DOX_pulse_IFNg / IFNg")  
ggsave(file.path(fig_dir, "volcano_DOX_pulse_IFNg_vs_IFNg_TEchanges.pdf"), width=4.5, height=3)

# scatter plot: 
gg <- .scatter_plot_TE(res=tidy_S6, label_x_pos=0, label_y_pos=-2, gap=0.3) +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       title="S6: Ribo-seq vs. RNA-seq", shape="IFNg_induced") 
ggsave(file.path(fig_dir, "TE_scatter_Ribo_vs_RNA_on_S6.pdf"), width=6, height=5)    

# go analysis on the `interested: down-regulated translation changes`
universal <- tidy_S6 %>% pull(ensembl) %>% str_replace("\\..*", "")
selected <- tidy_S6_down %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.2,
                         return.DEInCat=TRUE, dds=dds_S6)
selected <- tidy_S6_up %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go_up <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.2,
                            return.DEInCat=TRUE, dds=dds_S6)                         
write_xlsx(list(`CDS`=tidy_S6, `translation_upregulation` = tidy_S6_up,
                `enriched_GO_term_up` = enriched_go_up,
                `translation_downregulation` = tidy_S6_down,
                `enriched_GO_term_down` = enriched_go_down), 
           path=file.path(stats_dir, "S6_DOX_pulse_IFNg_vs_IFNg.xlsx"))                         


#
# save resutls
#
tidy_TE_CDS <- list(S1=tidy_S1, S2=tidy_S2, S3=tidy_S3, S4=tidy_S4, S5=tidy_S5, S6=tidy_S6)
save(tidy_TE_CDS, file=file.path(pkg_dir, "data", "tidy_TE_CDS.rda"))

# how many down-regualted translation changes genes are associated with the eniched GO terms?
tmp_S6 <- enriched_go_down %>%
  pull(DEInCat) %>% map(str_split, ",") %>% unlist() %>% unique()
length(tmp_S6)

# sanity check
# PRAMEF7 should have ~300 ribo-seq reads! but we only have single digit.

# per gene
.plot_counts(gene_id="ENSG00000204264.12", gene_name="PSMB8", dds=dds_S6)
ggsave(file.path(fig_dir, "PSMB8-S6.pdf"), width=4, height=3)

.plot_counts(gene_id="ENSG00000240065.8", gene_name="PSMB9", dds=dds_S6)
ggsave(file.path(fig_dir, "PSMB9-S6.pdf"), width=4, height=3)

.plot_counts(gene_id="ENSG00000205220.12", gene_name="PSMB10", dds=dds_S6)
ggsave(file.path(fig_dir, "PSMB10-S6.pdf"), width=4, height=3)

.plot_counts(gene_id="ENSG00000187837.4", gene_name="H1-2", dds=dds_S6)
ggsave(file.path(fig_dir, "H1-2-S6.pdf"), width=4, height=3)

#
# (7) S6 but with different threshold for DUX4_induced and INFg_induced
#

# exclude DUX4_incdued_v2
load(file.path(pkg_dir, "data", "IFNg_induced_v2.rda"))
load(file.path(pkg_dir, "data", "DUX4_induced_v2.rda"))
tidy_S6_v2 <- .tidy_res_ribo_over_rna(dds=dds_S6, inter_res=S6_ribo_over_rna, 
                                   ribo_res=S6_ribo, rna_res=S6_rna,
                                   alpha=0.05, lfc_threshold=1) %>%
             dplyr::filter(!ensembl %in% DUX4_induced_v2$ensembl) %>%
             dplyr::mutate(IFNg_induced_v2 = ensembl %in% IFNg_induced_v2$ensembl)

# scater plot: Ribo vs RNA
res <- tidy_S6_v2
msg_up   <- sprintf("log2(Ribo / RNA) > 1: %4.0f", sum(res$status == "up"))
msg_down <- sprintf("log2(Ribo / RNA) < -1: %4.0f", sum(res$status == "down"))  
pearson <- cor(res$rna_lfc, res$ribo_lfc)             
label_x_pos <- 0 
label_y_pos <- -2.5

ggplot(res, aes(x=rna_lfc, y=ribo_lfc)) +
  geom_point(size=2.5, alpha=0.5, aes(color=status, shape=IFNg_induced_v2)) +
  theme_bw() +
  #geom_vline(xintercept=c(-2, 2), linetype="dashed", alpha=0.5, color="gray50") +
  #geom_hline(yintercept=c(-2, 2), linetype="dashed", alpha=0.5, color="gray50") +
  annotate("text", x=label_x_pos, y=label_y_pos, label=paste0("Pearson = ", format(pearson, digit=2)),
           hjust = 0, vjust=1) +
  annotate("text", x=label_x_pos, y=label_y_pos-0.4, label=msg_up, color="red", hjust = 0, vjust=1) +
  annotate("text", x=label_x_pos, y=label_y_pos-0.8, label=msg_down, color="blue", hjust=0, vjust=1) +
  scale_color_manual(values=c("gray75", "red", "blue"), guide=FALSE) + #guide=FALSE
  scale_shape_manual(values=c(19, 1)) +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position="bottom", legend.box="vertical") +
  labs(x=bquote(~"RNA-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       y=bquote(~"Ribo-seq:" ~log[2]~"(DOX_pulse_IFNg / IFNg)"), 
       title="Ribo-seq vs. RNA-seq (exclude DUX4-induced)", shape="IFNg-induced") 
ggsave(file.path(fig_dir, "TE_scatter_Ribo_vs_RNA_on_S6-exclude-DUX4-induced.pdf"), width=6, height=5)     

# GO term
interested <- tidy_S6_v2 %>% dplyr::filter(status=="down")
library(goseq)
source(file.path(pkg_dir, "scripts", "tools.R"))
universal <- tidy_S6_v2 %>% pull(ensembl) %>% str_replace("\\..*", "")
selected <- interested %>% pull(ensembl) %>% str_replace("\\..*", "")
enriched_go <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.3,
                         return.DEInCat=TRUE, dds=dds_S6)
write_xlsx(list(`translation_downregulation` = interested,
                `Top_10_GO_term` = enriched_go[1:10, ]), 
           path=file.path(stats_dir, "S6_translation_downregulation-exclude-DUX4-indueced.xlsx"))                         

# how many down-regualted translation changes genes are associated with the eniched GO terms?
tmp_S6 <- readxl::read_xlsx(file.path(stats_dir, "S6_translation_downregulation.xlsx"), sheet=2) %>%
  pull(DEInCat) %>% map(str_split, ",") %>% unlist() %>% unique()
length(tmp_S6)

# are these RNA metabolic process associated translation changes also overlap with the DEU in Q1 or 1st exon?
metabolic <- unlist(str_split(enriched_go$DEInCat[5], ","))
metabolic_deu <- deu_Q1_S6 %>% dplyr::filter(gene_name %in% metabolic) # %>% dplyr::filter(log2fold_DOX_pulse_IFNg_IFNg.ribo < -0.5)

metabolic_deu <- deu_Q1_S6 %>% dplyr::filter(gene_name %in% metabolic, exonBaseMean > 5, log2fold_DOX_pulse_IFNg_IFNg.ribo > 0) %>%
  drop_na(padj.ribo, log2fold_DOX_pulse_IFNg_IFNg.ribo, log2fold_DOX_pulse_IFNg_IFNg.rna) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down, name) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down, -name) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) 


ggplot(metabolic_deu, aes(x=platform, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_point(aes(color=deu_down)) +
    geom_line(aes(group=gene_name, color=deu_down), show.legend=FALSE) +
    theme_bw() +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    #facet_wrap(~quantile, nrow=1) +
    theme(legend.position="bottom") + 
    scale_color_manual(values=c('#999999','#E69F00')) +
    labs(y="reletive exons usage LFC", title="metaboic: DOX_pulse_IFNg vs. IFNg", x="")  
ggsave(file.path(fig_dir, "test.pdf"))              

