#' Use DEXSeq to detect the differential exon usages
#' Use results from 07E-translation-specific-first-quarter-usage.R (S1_relative_exons_usage.rda, 
#' S2_relative_exons_usage.rda, S6_relative_exons_usage.rda)
#' to determine the elongation defect specificity in IFNg-stimulated and PSM
#' groups
#'

library(DEXSeq)
library(hg38.HomoSapiens.Gencode.v35)
library(tidyverse)
data(gene.anno)
gene.anno <- as.data.frame(gene.anno) 
txdb <- hg38.HomoSapiens.Gencode.v35

library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# parameters and load datasets
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "dexseq")
stat_dir <- file.path(pkg_dir, "stats", "differential_relative_exon_usage")

source(file.path(pkg_dir, "scripts", "07A-tools_dexseq.R"))
load(file.path(pkg_dir, "data", "S1_relative_exons_usage.rda"))
load(file.path(pkg_dir, "data", "S2_relative_exons_usage.rda"))
load(file.path(pkg_dir, "data", "S6_relative_exons_usage.rda"))
load(file.path(pkg_dir, "data", "DUX4_induced.rda"))
load(file.path(pkg_dir, "data", "IFNg_induced.rda"))

# there are about 90 PSMs; use the excel sheel from Danielle
psm_name <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "Gene Lists_Proteasome_MHC-I.xlsx"),
                              sheet=1, range="A4:B49") 
PSM <- gene.anno %>% dplyr::select(gene_name, gene_id) %>%
  dplyr::filter(gene_name %in% psm_name$Symbol)

.select_lowest_lfc <- function(.x)  {
 idx <- .x %>% dplyr::select(starts_with("log2fold") & ends_with(".ribo")) %>% 
   as.matrix() %>% which.min(.)
 .x[idx, ]
}

.select_best_deu <- function(.x)  {
  # select the best representative in the quarter
  # (1) ExonBaseMean (Ribo) > 5; (2) drop na in log2fold in ribo
  # (3) take the one with lowest LFC usage
 .x <- .x %>%  dplyr::filter(exonBaseMean >= 5) %>%
   tidyr::drop_na(starts_with("log2fold") & ends_with(".ribo")) %>%
   tidyr::drop_na(padj.ribo) 
 if (nrow(.x) > 1) {
   idx <- .x %>%   dplyr::select(starts_with("log2fold") & ends_with(".ribo")) %>% 
     as.matrix() %>% which.min(.)
   .x <- .x[idx, ]
 }
 return(.x)
}

#
# (1) S6: Q1
#
.tidy_Q6_per_quantile <- function(ribo, rna) {
  ribo %>% 
    left_join(rna, by="name", suffix=c(".ribo", ".rna")) %>%
    dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.ribo < -1,
                  rna_deu_down = padj.rna < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.rna < -1,
                  ribo_deu_up = padj.ribo < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.ribo > 1,
                  rna_deu_up = padj.rna < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.rna >1) %>%
    dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
    dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%                
    dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
    dplyr::mutate(IFNg_induced = groupID %in% IFNg_induced$ensembl) %>%
    rename(gene_id="groupID") %>%
    left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
    dplyr::relocate(gene_name, .after="gene_id") %>%
    arrange(padj.ribo)
}

Q1_S6 <- S6_relative_exons_usage$Q1
Q2_S6 <- S6_relative_exons_usage$Q2
Q3_S6 <- S6_relative_exons_usage$Q3
Q4_S6 <- S6_relative_exons_usage$Q4

deu_Q1_S6 <- .tidy_Q6_per_quantile(Q1_S6$ribo, Q1_S6$rna)
deu_Q2_S6 <- .tidy_Q6_per_quantile(Q2_S6$ribo, Q2_S6$rna)
deu_Q3_S6 <- .tidy_Q6_per_quantile(Q3_S6$ribo, Q3_S6$rna)
deu_Q4_S6 <- .tidy_Q6_per_quantile(Q4_S6$ribo, Q4_S6$rna)

#
# (1A) PSM
# 
# drop NA: padj.ribo and lfc on either) 18 remains and pick the one with lowest lfc

# testing: PSMB9 ENSG00000240065.8
.x <- deu_Q1_S6 %>% dplyr::filter(gene_id=="ENSG00000240065.8")
.x <- deu_Q1_S6 %>% dplyr::filter(gene_id=="ENSG00000125818.18")

deu_psm_Q1 <- deu_Q1_S6 %>% dplyr::filter(gene_id %in% PSM$gene_id)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

write_xlsx(deu_psm_Q1, path=file.path(stat_dir, "S6_Q1_PSM_best_deu.xlsx"))

deu_psm_Q2 <- deu_Q2_S6 %>% dplyr::filter(gene_id %in% PSM$gene_id)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_psm_Q3 <- deu_Q3_S6 %>% dplyr::filter(gene_id %in% PSM$gene_id)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_lowest_lfc(.x)) %>% ungroup() 

deu_psm_Q4 <- deu_Q4_S6 %>% dplyr::filter(gene_id %in% PSM$gene_id)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

# tidy up the PSM for print: visualize by trajectory plot in Q1 between ribo and rna (18)
q1 <- deu_psm_Q1 %>% drop_na(padj.rna) %>%
  #dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down, name) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down, -name) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q1")

q2 <- deu_psm_Q2 %>% drop_na(padj.rna) %>%
  #dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down, name) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down, -name) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q2")
q3 <- deu_psm_Q3 %>% drop_na(padj.rna) %>%
  #dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down, name) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down, -name) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q3")  
q4 <- deu_psm_Q4 %>% drop_na(padj.rna) %>%
  #dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down, name) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down, -name) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q4")

# viz usign the tidy dataset q1, q2, q3, q4
library(ggrepel)
q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>%
  ggplot(aes(x=platform, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_point(aes(color=deu_down)) +
    ggrepel::geom_text_repel(aes(label=gene_name), size=1.2, show.legend=FALSE) +
    geom_line(aes(group=gene_name, color=deu_down), alpha=0.5, show.legend=FALSE) +
    theme_bw() +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    facet_wrap(~quantile, nrow=1) +
    theme(legend.position="bottom") + 
    scale_color_manual(values=c('#999999','#E69F00')) +
    labs(y="reletive exons usage LFC", title="PSM: DOX_pulse_IFNg vs. IFNg", x="")            
ggsave(file.path(fig_dir, "S6-PSM-trajectory.pdf"), width=6, height=4)               

q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>% 
  dplyr::filter(platform=="Ribo") %>%
  ggplot(aes(x=quantile, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    theme_bw() + 
    labs(y="Ribo: reletive exons usage LFC", title="PSM: DOX_pulse_IFNg vs. IFNg", x="") 
ggsave(file.path(fig_dir, "S6-PSM-boxplot-all-quantiles.pdf"), width=4, height=3)                     

#
# (1B) IFNg
# 
# drop NA: padj.ribo and lfc on either) 18 remains and pick the one with lowest lfc
deu_IFNg_Q1 <- deu_Q1_S6 %>% dplyr::filter(IFNg_induced)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 


deu_IFNg_Q2 <- deu_Q2_S6 %>% dplyr::filter(IFNg_induced)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_IFNg_Q3 <- deu_Q3_S6 %>% dplyr::filter(IFNg_induced)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_IFNg_Q4 <- deu_Q4_S6 %>% dplyr::filter(IFNg_induced)  %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

q1 <- deu_IFNg_Q1 %>% drop_na(padj.rna) %>%
  dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q1")
q2 <- deu_IFNg_Q2 %>% drop_na(padj.rna) %>%
  dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q2")
q3 <- deu_IFNg_Q3 %>% drop_na(padj.rna) %>%
  dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q3")  
q4 <- deu_IFNg_Q4 %>% drop_na(padj.rna) %>%
  dplyr::filter((IFNg.rna + DOX_pulse_IFNg.rna) / 2 >= 5) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q4")

q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>%
  ggplot(aes(x=platform, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_point(aes(color=deu_down)) +
    geom_line(aes(group=gene_name, color=deu_down), alpha=0.3, show.legend=FALSE) +
    theme_bw() +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    facet_wrap(~quantile, nrow=1) +
    theme(legend.position="bottom") + 
    scale_color_manual(values=c('#999999','#E69F00')) +
    labs(y="reletive exons usage LFC", title="IFNg-stimulated: DOX_pulse_IFNg vs. IFNg", x="")            
ggsave(file.path(fig_dir, "S6-IFNg-stimulated-trajectory.pdf"), width=6, height=4)     

q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>% 
  dplyr::filter(platform=="Ribo") %>%
  ggplot(aes(x=quantile, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    theme_bw() + 
    labs(y="Ribo: reletive exons usage LFC", title="IFNg-stimulated: DOX_pulse_IFNg vs. IFNg", x="") 
ggsave(file.path(fig_dir, "S6-IFNg-stimulated-boxplot-all-quantiles.pdf"), width=4, height=3)  

#
# (1C) others
#
deu_others_Q1 <- deu_Q1_S6 %>% dplyr::filter(!IFNg_induced)  %>%
  dplyr::filter(!gene_id %in% PSM$gene_id) %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_others_Q2 <- deu_Q2_S6 %>% dplyr::filter(!IFNg_induced)  %>%
  dplyr::filter(!gene_id %in% PSM$gene_id) %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_others_Q3 <- deu_Q3_S6 %>% dplyr::filter(!IFNg_induced)  %>%
  dplyr::filter(!gene_id %in% PSM$gene_id) %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 

deu_others_Q4 <- deu_Q4_S6 %>% dplyr::filter(!IFNg_induced)  %>%
  dplyr::filter(!gene_id %in% PSM$gene_id) %>%
  group_by(gene_id) %>%
  group_modify(~.select_best_deu(.x)) %>% ungroup() 


q1 <- deu_others_Q1 %>%   drop_na(log2fold_DOX_pulse_IFNg_IFNg.rna)  %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q1")
q2 <- deu_others_Q2 %>% drop_na(log2fold_DOX_pulse_IFNg_IFNg.rna)  %>% 
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q2")
q3 <- deu_others_Q3 %>% drop_na(log2fold_DOX_pulse_IFNg_IFNg.rna)  %>% %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q3")  
q4 <- deu_others_Q4 %>% drop_na(log2fold_DOX_pulse_IFNg_IFNg.rna)  %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  dplyr::select(starts_with("log2fold"), gene_id, gene_name, deu_down) %>%
  gather(key=platform, value=log2fold_DOX_pulse_IFNg_IFNg, -gene_id, -gene_name, -deu_down) %>%
  dplyr::mutate(platform=if_else(grepl(".rna", platform), "RNA", "Ribo"))  %>%
  dplyr::mutate(platform = factor(platform, levels=c("RNA", "Ribo"))) %>%
  add_column(quantile="Q4")

q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>%
  ggplot(aes(x=platform, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    theme_bw() +
    geom_boxplot(aes(color=platform), width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    facet_wrap(~quantile, nrow=1) +
    theme(legend.position="bottom") + 
    labs(y="reletive exons usage LFC", title="Others: DOX_pulse_IFNg vs. IFNg", x="")  +
    ylim(c(-2.5, 2.5))          
ggsave(file.path(fig_dir, "S6-others-trajectory.pdf"), width=6, height=4)     

q1 %>% bind_rows(q2) %>% bind_rows(q3) %>% bind_rows(q4) %>% 
  dplyr::filter(platform=="Ribo") %>%
  ggplot(aes(x=quantile, y=log2fold_DOX_pulse_IFNg_IFNg)) +
    geom_boxplot(width = 0.3, fill="transparent", outlier.shape=NA, 
                 show.legend=FALSE) +
    theme_bw() + 
    labs(y="Ribo: reletive exons usage LFC", title="IFNg-stimulated: DOX_pulse_IFNg vs. IFNg", x="") 
ggsave(file.path(fig_dir, "S6-IFNg-stimulated-boxplot-all-quantiles.pdf"), width=4, height=3)  
#
# (1D) Q1 - PSM, IFNg-stimulated, and others
#


g1 <- dplyr::select(deu_others_Q1,  S6: PSM relative exon usage) %>% 
  add_column(group="Others")
g2 <-   dplyr::select(deu_psm_Q1,  log2fold_DOX_pulse_IFNg_IFNg.ribo) %>% 
  add_column(group="PSM")
g3 <- dplyr::select(deu_IFNg_Q1,  log2fold_DOX_pulse_IFNg_IFNg.ribo) %>% 
  add_column(group="IFNg-stimulated")
g1 %>% bind_rows(g2) %>% bind_rows(g3) %>%
  dplyr::mutate(group=factor(group, levels=c("PSM", "IFNg-stimulated", "Others"))) %>%
  ggplot(aes(x=group, y=log2fold_DOX_pulse_IFNg_IFNg.ribo)) +
    geom_boxplot(width=0.3, fill="transparent", outlier.shape=NA, show.legend=FALSE) +
    theme_bw() +
    labs(y="Ribo: reletive exons usage LFC", title="DOX_pulse_IFNg vs. IFNg in Q1") +
    ylim(c(-5, 2.5))
ggsave(file.path(fig_dir, "S6-Q1-LFC-boxplot.pdf"), width=4, height=3)    
wilcox.test(g1$log2fold_DOX_pulse_IFNg_IFNg.ribo, g2$log2fold_DOX_pulse_IFNg_IFNg.ribo)
wilcox.test(g1$log2fold_DOX_pulse_IFNg_IFNg.ribo, g3$log2fold_DOX_pulse_IFNg_IFNg.ribo)