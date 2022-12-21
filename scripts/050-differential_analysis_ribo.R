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
library(goseq)
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
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))
load(file.path(pkg_dir, "data", "rse_5p_3p_cds_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_around_TSS.rda"))
load(file.path(pkg_dir, "data", "DUX4_altered.rda"))
load(file.path(pkg_dir, "data", "IFNg_altered.rda"))


# list of comparisons
list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))

# get histone H1 and H2
histone_variants <- .histone_variants(gene.anno)
#
# load DUX4_altered and INFg_altered
#
load(file.path(pkg_dir, "data", "DUX4_induced.rda"))
load(file.path(pkg_dir, "data", "IFNg_induced.rda"))

#
# tools
#
# remove histone_variants (H1/H2)
.remove_histone_variants <- function(dds, histone_variants) {
  dds[!rowData(dds)$gene_id %in% histone_variants$gene_id]
}

.subset_by_treatment <- function(treatments, dds, threshold=12) {
    # treatments must include to treatments for comparison
    sub <- dds[, dds$treatment %in% treatments]
    sub <- sub[rowSums(counts(sub)) >= threshold, ]
    sub$treatment <- factor(sub$treatment, levels=treatments)
    sub <- DESeq(sub)
}

.subset_by_treatment_sizeFactor_byCDS <- function(treatments, res_tx, row_mean=5, rse_cds_by_gene) {
    # the dds_cds is used for sizeFactor
    # treatments must include to treatments for comparison
    dds_cds <- DESeqDataSet(rse_cds_by_gene[, rse_cds_by_gene$treatment %in% treatments])
    dds_cds <- estimateSizeFactors(dds_cds)

    dds_tx <- DESeqDataSet(rse_tx[, rse_tx$treatmetn %in% treatments])
    dds_tx$treatment <- factor(dds_tx$treatment, levels=treatments)
    dds_tx$sizeFactor <- dds_cds$sizeFactor
    dds_tx <- dds_tx[rowMeans(counts(dds_tx, normalized=TRUE)) >= row_mean, ]
    dds_tx <- DESeq(dds_tx)
}

#
# (1) CDS by genes: Comparison between groups of treatments
#
fig_dir <- file.path(pkg_dir, "figures", "cds_deseq2")

# (1a) check some gene counts

#ISG20, IFIH1, CD74, CXCL9, HLA-A, HLA-B, HLA-C are also good IFNg stimulated genes to check.
#PSMB5, PSMB6, PSMB7 are constitutive proteasome subunits that should be present in all the samples.
#RPL27, RPL13A, GAPDH are housekeeping protein coding genes to check.
psmb <- c("PSMB5", "PSMB6", "PSMB7" ,"PSMB8", "PSMB9", "PSMB10")
psmb_cnt <- .get_gene_counts(dds_cds_by_gene, psmb)
ggplot(psmb_cnt, aes(x=treatment, y=count)) +
  geom_point() +
  facet_wrap(~ gene_name, scale="free_y", nrow=3) +
  theme_bw() +
  labs(title="proteasome subunits", y="normalized count")
ggsave(file.path(fig_dir, "PSMB_norm_counts_R1.pdf"))

ifg <- c("ISG20", "IFIH1", "CD74", "CXCL9", "HLA-A", "HLA-B", "HLA-C")
ifg_cnt <- .get_gene_counts(dds_cds_by_gene, ifg)
ggplot(ifg_cnt, aes(x=treatment, y=count)) +
  geom_point() +
  facet_wrap(~ gene_name, scale="free_y", nrow=4) +
  theme_bw() +
  labs(title="IFNg stimuated genes", y="normalized count")
ggsave(file.path(fig_dir, "IFNg_norm_counts_R1.pdf"))

house_keeping <- c("RPL27", "RPL13A", "GAPDH")
cnt <- .get_gene_counts(dds_cds_by_gene, house_keeping)
ggplot(cnt, aes(x=treatment, y=count)) +
  geom_point() +
  facet_wrap(~ gene_name, scale="free", nrow=2) +
  theme_bw() +
  labs(title="House keeping genes", y="normalized count")
ggsave(file.path(fig_dir, "house_keeping_norm_counts_R1.pdf"))  

#
# (1b) CDS: comparison between groups of treatment: DESeq1
#      We made six comparisons which include all the factorial combination of the treatments. 
#      The hypothesis tests are based on $H_0: |logFC| \le 0$ and $H_1: |logFC| > 0$. The 
#      threshold of adjusted p-value is $0.05$ and logFC cutoff is 1. 
#      subject to CDS with 12 reads more or more RPF across six samples.
#

tmp_cds <- DESeqDataSet(rse_cds_by_gene, design = ~treatment)  
deseq2_cds <- lapply(list_comp, .subset_by_treatment, dds=tmp_cds)    
save(deseq2_cds, file=file.path(pkg_dir, "data", "deseq2_cds.rda"))

# Tidy resuls for all senarios: lfcThreshold=1; alpha=0.05   
load(file.path(pkg_dir, "data", "deseq2_cds.rda")) 

cds_summarize_1 <- map_dfr(deseq2_cds, function(x) {
  res <- results(x, lfcThreshold=0, alpha=0.05)
  as.data.frame(res) %>% dplyr::filter(padj < 0.05) %>%
    summarise(`LFC > 1 (up)` = sum(log2FoldChange > 1, na.rm=TRUE),
              `LFC < -1 (down)` = sum(log2FoldChange < -1, na.rm=TRUE))
}, .id="name") %>%
  add_column(comparison=map_chr(list_comp, function(x) paste0(x[2], "-vs-", x[1])), .after="name")

cds_summarize_1

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
file_name <- file.path(pkg_dir, "stats", "differential_analysis", "deseq2_cds_riboseq.xlsx")
writexl::write_xlsx(
  list(DOX_pulse_vs_untreated      = tidy_results[["S1"]] ,
       IFNg_vs_untreated           = tidy_results[["S2"]],
       DOX_pulse_IFNg_vs_untreated = tidy_results[["S3"]],
       DOX_pulse_vs_IFNg           = tidy_results[["S4"]],
       DOX_pulse_IFNg_vs_DOX_pulse = tidy_results[["S5"]],
       DOX_pulse_IFNg_vs_IFNg      = tidy_results[["S6"]]),
  path=file_name)

#
# S1 MA plot
# 

# (a) how many of DUX4_targets are also have translational changes?
tidy_results[["S1"]] %>% dplyr::filter(significant) %>% summarise(total=n(), dux4=sum(DUX4_targets))

# (b) MAplot: log10(baseMean) vs. -log10(FDR)
tidy_results[["S1"]] %>% 
  dplyr::mutate(log10_baseMean=log10(baseMean)) %>%
  ggplot(aes(x=log10_baseMean, y=log2FoldChange, color=significant)) +
    geom_point(alpha=0.5, show.legend=FALSE) +
    theme(legend.position="none") +
    theme_bw() +
    scale_color_manual(values=c("grey", "#F21A00")) +
    labs(title="DOX_pulse / untreated", 
         x="log10 (base mean)", y="Ribo-seq \nlog2FC (DOX_pulse / untreated)") +
         annotate("text", x=3, y=-5, label="up: 495 \ndown: 356", hjust = 0, vjust=1) 
ggsave(file.path(fig_dir, "S1-MAPlot.pdf"), width=3, height=3)  



##################################
#
# (3) around_TSS differential analysis: rowMeans filter: 6 across 6 samples
#
##################################
load(file.path(pkg_dir, "data", "rse_around_TSS.rda"))
load(file.path(pkg_dir, "data", "deseq2_cds.rda"))
dds_tss <- DESeqDataSet(rse_around_TSS, design = ~treatment)  
dds_tss <- .remove_histone_variants(dds_tss, histone_variants)

#
# (3a) S1
#
senario <- "S1" # DUX4 vs. untreated
sub <- dds[, dds_tss$treatment %in% list_comp[[senario]]]
sub <- sub[rowMeans(counts(sub)) >= 6, ]
sub$treatment <- factor(sub$treatment, levels=list_comp[[senario]])
sizeFactors(sub) <- sizeFactors(deseq2_cds[[senario]])
sub <- DESeq(sub)
res_tss  <- results(sub, alpha=0.05)
summary(res_tss)

# summarize down-regulated ribosome occupancy around TSS [-10, 10]
tss_down_S1 <- results(sub, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange < -1) %>%
    rename(tx_name="row") %>%
    mutate(gene_id   = rowData(dds_tss)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_tss)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(DUX4_induced = gene_id %in% DUX4_induced$ensembl) %>%
    dplyr::mutate(IFNg_induced = gene_id %in% IFNg_induced$ensembl, .after = DUX4_induced)
norm_cnt <- counts(sub[tss_down_S1$tx_name], normalized=TRUE)
tss_down_S1 <- tss_down_S1 %>% add_column(as.data.frame(norm_cnt))

# go analysis on the `interested`
source(file.path(pkg_dir, "scripts", "tools.R"))
universal <- rownames(deseq2_cds[[senario]]) %>% str_replace("\\..*", "")
selected <- tss_down_S1 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
length(selected) # 550
enriched_go <- .do_goseq(universal=universal, selected_genes=selected, 
                         return.DEInCat=TRUE, dds=deseq2_cds$S1, p_value = 0.01)
write_xlsx(list(`TSS_downreg_occupancy` = tss_down_S1, 
                `enriched_GO_term` = enriched_go), 
           path=file.path(pkg_dir, "stats", "differential_analysis", "TSS_S1_down_occupancy.xlsx")) 

# output transcripts that do not loss occupancy
tss_no_suppress <- results(sub, tidy=TRUE) %>% 
  dplyr::filter(abs(log2FoldChange) < 0.1) %>% 
  rename(tx_name="row") %>%
  mutate(as.data.frame(rowData(sub[tx_name])[, c("gene_id", "gene_name")])) %>%
  mutate(as.data.frame(counts(sub[tx_name], normalized=TRUE))) 
writexl::write_xlsx(tss_no_suppress, 
                    path=file.path(pkg_dir, "stats", "differential_analysis", "TSS_S1_not_suppressed_in_DUX4.xlsx"))  

#
# (3b) S6
#
senario <- "S6" # DOX-pulse-IFNg vs. IFNg
sub <- dds_tss[, dds_tss$treatment %in% list_comp[[senario]]]
sub <- sub[rowMeans(counts(sub)) >= 6, ]
sub$treatment <- factor(sub$treatment, levels=list_comp[[senario]])
sizeFactors(sub) <- sizeFactors(deseq2_cds[[senario]])
sub <- DESeq(sub)
res_tss  <- results(sub, alpha=0.05)
summary(res_tss)

# summarize down-regulated ribosome occupancy around TSS [-10, 10]
tss_S6_results <- results(sub, alpha=0.05, tidy=TRUE) %>%
    rename(tx_name="row") 

tss_down_S6 <- results(sub, alpha=0.05, tidy=TRUE) %>%
  dplyr::filter(padj < 0.05, log2FoldChange < -1) %>%
    rename(tx_name="row") %>%
    mutate(gene_id   = rowData(dds_tss)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_tss)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(DUX4_induced = gene_id %in% DUX4_induced$ensembl) %>%
    dplyr::mutate(IFNg_induced = gene_id %in% IFNg_induced$ensembl, .after = DUX4_induced) 

norm_cnt <- counts(sub[tss_down_S6$tx_name], normalized=TRUE)
tss_down_S6 <- tss_down_S6 %>% add_column(as.data.frame(norm_cnt))

tss_down_S6 %>% distinct(gene_name, .keep_all=TRUE) %>% filter(IFNg_induced) %>% nrow() #47

# go analysis on the `interested`
source(file.path(pkg_dir, "scripts", "tools.R"))
universal <- rownames(deseq2_cds[[senario]]) %>% str_replace("\\..*", "")
selected <- tss_down_S6 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
length(selected) # 744
enriched_go_S6 <- .do_goseq(universal=universal, selected_genes=selected, 
                         return.DEInCat=TRUE, dds=deseq2_cds$S6, p_value = 0.01)
write_xlsx(list(`TSS_downreg_occupancy` = tss_down_S6, 
                `enriched_GO_term` = enriched_go_S6), 
           path=file.path(pkg_dir, "stats", "differential_analysis", "TSS_S6_down_occupancy.xlsx")) 

# PSM genes?
psm_name <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "Gene Lists_Proteasome_MHC-I.xlsx"),
                              sheet=1, range="A4:B49") 
PSM <- as.data.frame(gene.anno) %>% dplyr::select(gene_name, gene_id) %>%
  dplyr::filter(gene_name %in% psm_name$Symbol)

tss_PSM_S6 <- results(sub, alpha=0.05, tidy=TRUE) %>%
    rename(tx_name="row") %>%
    dplyr::mutate(gene_id   = rowData(dds)[tx_name, "gene_id"], .before=baseMean) %>%
    dplyr::mutate(gene_name = rowData(dds)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj) %>%
    dplyr::filter(gene_name %in% PSM$gene_name) 
norm_cnt <- counts(sub[tss_PSM_S6$tx_name], normalized=TRUE)
tss_PSM_S6 <- tss_PSM_S6 %>% add_column(as.data.frame(norm_cnt))    

write_xlsx(list(`TSS_downreg_occupancy` = tss_down_S6, 
                `enriched_GO_term` = enriched_go_S6,
                `PSM` = tss_PSM_S6), 
           path=file.path(pkg_dir, "stats", "differential_analysis", "TSS_S6_down_occupancy.xlsx")) 

# IFNg genes?

##################################
#
# (3) 5'UTR differential analysis
#
##################################

load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))
load(file.path(pkg_dir, "data", "rse_5p_3p_cds_by_tx.rda"))
library(wesanderson)
my_color <- wes_palette("Darjeeling1", n=5)[2:5]
names(my_color) <- levels(dds_cds_by_gene$treatment)

# note: rse_5p, rse_3p and rse_cds are using the same sizeFactor inherited from dss_cds
rse_5p  <- rse_5p_3p_cds_by_tx$rse_5p
rse_3p  <- rse_5p_3p_cds_by_tx$rse_3p
rse_cds <- rse_5p_3p_cds_by_tx$rse_cds
dds_5p  <- DESeqDataSet(rse_5p, design = ~treatment)              
dds_3p  <- DESeqDataSet(rse_3p, design = ~treatment)  
dds_cds <- DESeqDataSet(rse_cds, design = ~treatment)  
# remove histone_variants (H1/H2)
.remove_histone_variants <- function(dds, histone_variants) {
  dds[!rowData(dds)$gene_id %in% histone_variants$gene_id]
}
dds_5p <- .remove_histone_variants(dds_5p, histone_variants)
dds_3p <- .remove_histone_variants(dds_3p, histone_variants)
dds_cds <- .remove_histone_variants(dds_cds, histone_variants)

#
# (3a) DOX_pulse_IFNg vs. IFNg (S6)
#
fig_dir <- file.path(pkg_dir, "figures", "utr_deseq2", "DOX_pulse_IFNg_vs_IFNg")
sub <- dds_5p[, dds_5p$treatment %in% c("IFNg", "DOX_pulse_IFNg")]
sub$treatment <- factor(sub$treatment, levels=c("IFNg", "DOX_pulse_IFNg"))
sub <- sub[rowSums(counts(sub)) >= 12]
#sub <- estimateSizeFactors(sub)
sub <- DESeq(sub, BPPARAM=bp_param)
res_5p <- results(sub, alpha=0.05, parallel=TRUE)
summary(res_5p)

# boxplot of down-regulated 5p (5'p/CDS/3'p)
sig_down <- as.data.frame(res_5p) %>%
    dplyr::filter(padj < 0.05, log2FoldChange < 0) %>%
    rownames_to_column(var="tx_name") %>%
    mutate(gene_id   = rowData(dds_5p)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_5p)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj)
sig_down_S6 <- sig_down
write_xlsx(sig_down_S6, path=file.path(pkg_dir, "stats", "differential_analysis", "5UTR_S6_down_occupancy.xlsx"))

#
# compare with tss_down_S6: overlapping genes?
#
length(unique(tss_down_S6$gene_name))
length(unique(sig_down_S6$gene_name))
sum(unique(tss_down_S6$gene_name) %in% unique(sig_down_S6$gene_name))
tmp <- sig_down_S6 %>% left_join(tss_S6_results, by="tx_name", suffix=c(".5p", ".TSS")) %>%
  drop_na(log2FoldChange.TSS) %>%
  dplyr::mutate(PSM=gene_name %in% PSM$gene_name)
ggplot(tmp, aes(x=log2FoldChange.5p, y=log2FoldChange.TSS)) +
  geom_point(aes(color=PSM), alpha=0.5) +
  geom_smooth(method="lm") +
  theme_minimal()
ggsave(file.path(fig_dir, "test.pdf"))  


#
# VIZ 1: boxplot of down-regulated occupancy in 5'UTR; combine 5', 3' and CDS reads  and TSS (tss_S6_results)
#
df <- .tidy_combo_2(sig=sig_down, dds_5p=dds_5p, dds_3p=dds_3p, dds_cds=dds_cds, dds_tss=dds_tss, 
                    samples=sub$sample_name) %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

df %>%
  dplyr::mutate(treatment = factor(treatment)) %>%
  group_by(region)%>%
  spread(key=treatment, log2_cnt) %>%
  summarise(DOX_pulse_IFNg_mean=mean(DOX_pulse_IFNg, na.rm=TRUE),
            IFNg_mean=mean(IFNg, na.rm=TRUE),
            p_value=t.test(DOX_pulse_IFNg, IFNg)$p.value,
            statistic=t.test(DOX_pulse_IFNg, IFNg)$statistic)
# boxplot            
ggplot(df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.8)  +
    facet_wrap(~ region, nrow=1) +
    theme_bw() +
    scale_fill_manual(values=my_color[c("IFNg", "DOX_pulse_IFNg")]) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="Significant 5' UTR (3092 tx; down)", y="log2(norm count +1)", x="samples")
ggsave(file.path(fig_dir, "five_prime_UTR_sig_down_occupancy_S6.pdf"), width=4, height=3)

#
# VIZ 2: boxplot of down-regulated occupacy in TSS; combine all four regions
#
df <- .tidy_combo_2(sig=tss_down_S6, dds_5p=dds_5p, dds_3p=dds_3p, dds_cds=dds_cds, dds_tss=dds_tss, 
                  samples=sub$sample_name) %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

df %>%
  dplyr::mutate(treatment = factor(treatment)) %>%
  group_by(region)%>%
  spread(key=treatment, log2_cnt) %>%
  summarise(DOX_pulse_IFNg_mean=mean(DOX_pulse_IFNg, na.rm=TRUE),
            IFNg_mean=mean(IFNg, na.rm=TRUE),
            p_value=t.test(DOX_pulse_IFNg, IFNg)$p.value,
            statistic=t.test(DOX_pulse_IFNg, IFNg)$statistic)

ggplot(df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.8)  +
    facet_wrap(~ region, nrow=1) +
    theme_bw() +
    scale_fill_manual(values=my_color[c("IFNg", "DOX_pulse_IFNg")]) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="Significant TSS (2640 tx; down)", y="log2(norm count +1)", x="samples")
ggsave(file.path(fig_dir, "TSS_sig_down_occupancy_S6.pdf"), width=4, height=3)

#
# VIS 3: plot all - size factor by CDS occupancy
#
tmp <- as.data.frame(res_5p) %>%
    rownames_to_column(var="tx_name") %>%
    mutate(gene_id   = rowData(dds_5p)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_5p)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj)
df <- .tidy_combo_2(sig=tmp, dds_5p=dds_5p, dds_3p=dds_3p, dds_cds=dds_cds, dds_tss=dds_tss, samples=sub$sample_name) %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

# t-test
t_test <- df %>%
  dplyr::mutate(treatment = factor(treatment)) %>%
  group_by(region)%>%
  spread(key=treatment, log2_cnt) %>%
  summarise(DOX_pulse_IFNg_mean=mean(DOX_pulse_IFNg, na.rm=TRUE),
            IFNg_mean=mean(IFNg, na.rm=TRUE),
            p_value=t.test(DOX_pulse_IFNg, IFNg)$p.value,
            statistic=t.test(DOX_pulse_IFNg, IFNg)$statistic)

ggplot(df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.8)  +
    facet_wrap(~ region, nrow=1) +
    theme_bw() +
    scale_fill_manual(values=my_color[c("IFNg", "DOX_pulse_IFNg")]) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="all 5' UTR", y="log2(norm count +1)", x="samples")
ggsave(file.path(fig_dir, "five_prime_UTR_occupancy_sizeFactorByCDS.pdf"), width=4, height=3)

#
# (2b) DOX_pulse vs. untreated (S1)
#
fig_dir <- file.path(pkg_dir, "figures", "utr_deseq2", "DOX_pulse_vs_untreated")
sub <- dds_5p[, dds_5p$treatment %in% c("DOX_pulse", "untreated")]
sub$treatment <- factor(sub$treatment, levels=c("untreated", "DOX_pulse"))
sub <- sub[rowSums(counts(sub)) >= 12]
sub <- DESeq(sub)
res_5p <- results(sub, alpha=0.05)
summary(res_5p)
res_5p_S1 <- res_5p

sig_down <- as.data.frame(res_5p) %>%
    dplyr::filter(padj < 0.05, log2FoldChange < 0) %>%
    rownames_to_column(var="tx_name") %>%
    mutate(gene_id   = rowData(dds_5p)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_5p)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj)
sig_down_S1 <- sig_down
write_xlsx(sig_down_S1, path=file.path(pkg_dir, "stats", "differential_analysis", "5UTR_S1_down_occupancy.xlsx"))

df <- .tidy_combo(sig=sig_down, dds_5p, dds_3p, dds_cds, samples=sub$sample_name) %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

# t-test
df %>%
  dplyr::mutate(treatment = factor(treatment)) %>%
  group_by(region)%>%
  spread(key=treatment, log2_cnt) %>%
  summarise(DOX_pulse_mean=mean(DOX_pulse, na.rm=TRUE),
            untreated_mean=mean(untreated, na.rm=TRUE),
            p_value=t.test(DOX_pulse, untreated)$p.value,
            statistic=t.test(untreated, untreated)$statistic)

ggplot(df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.7)  +
    facet_wrap(~ region, nrow=1) +
    theme_bw() +
    scale_fill_manual(values=my_color[c("DOX_pulse", "untreated")]) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="Significant 5' UTR (4408 tx; down)", y="log2(norm count +1)", x="samples")
ggsave(file.path(fig_dir, "five_prime_UTR_sig_down_occupancy.pdf"), width=4, height=3)

# all tx
tmp <- as.data.frame(res_5p) %>%
    rownames_to_column(var="tx_name") %>%
    mutate(gene_id   = rowData(dds_5p)[tx_name, "gene_id"], .before=baseMean) %>%
    mutate(gene_name = rowData(dds_5p)[tx_name, "gene_name"], .before=baseMean) %>%
    dplyr::arrange(padj)
df <- .tidy_combo(sig=tmp, dds_5p=dds_5p, dds_3p=dds_3p, dds_cds=dds_cds, samples=sub$sample_name) %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

# t-test
df %>%
  dplyr::mutate(treatment = factor(treatment)) %>%
  group_by(region)%>%
  spread(key=treatment, log2_cnt) %>%
  summarise(DOX_pulse_mean=mean(DOX_pulse, na.rm=TRUE),
            untreated_mean=mean(untreated, na.rm=TRUE),
            p_value=t.test(DOX_pulse, untreated)$p.value,
            statistic=t.test(untreated, untreated)$statistic)

ggplot(df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.7)  +
    facet_wrap(~ region, nrow=1) +
    theme_bw() +
    scale_fill_manual(values=my_color[c("DOX_pulse", "untreated")]) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="All 5' UTR", y="log2(norm count +1)", x="samples")
ggsave(file.path(fig_dir, "five_prime_UTR_all_occupancy.pdf"), width=4, height=3)            


#
# 3. all samples together
#
sample_info <- as.data.frame(colData(dds_5p))
# select transcripts that has some occupancy on CDS
tmp <- dds_cds[rowSums(counts(dds_cds)) >= 12]
tx_name <- rownames(tmp)
mat_5p <- counts(dds_5p[rownames(dds_5p) %in% tx_name],  normalized=TRUE) %>%
  as.data.frame() %>% 
  rownames_to_column(var="tx_name") %>%
  gather(key=sample_name, value=cnt, -tx_name) %>%
  add_column(region="5' UTR")
mat_3p <- counts(dds_3p[rownames(dds_3p) %in% tx_name], normalized=TRUE) %>%
  as.data.frame() %>% 
  rownames_to_column(var="tx_name") %>%
  gather(key=sample_name, value=cnt, -tx_name) %>%
  add_column(region="3' UTR")
mat_cds <- counts(dds_cds[rownames(dds_cds) %in% tx_name], normalized=TRUE) %>%
  as.data.frame() %>% 
  rownames_to_column(var="tx_name") %>%
  gather(key=sample_name, value=cnt, -tx_name) %>%
  add_column(region="CDS") 

comb_df <- add_row(mat_5p, mat_3p) %>%
  add_row(mat_cds) %>%
  mutate(region=factor(region, levels=c("5' UTR", "CDS", "3' UTR")))  %>%
  mutate(log2_cnt = log2(cnt+1)) %>%
  left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  arrange(treatment) %>%
  dplyr::mutate(sample_name=factor(sample_name, levels=unique(sample_name)))

ggplot(comb_df, aes(x=sample_name, y=log2_cnt, fill=treatment)) +
    geom_boxplot(outlier.colour="black", outlier.shape=NA, alpha=0.7)  +
    facet_wrap(~ region, nrow=1, scale="free_y") +
    theme_bw() +
    scale_fill_manual(values=my_color) +
    theme(axis.text.x = element_blank(), legend.position="bottom") +
    labs(title="All 5' UTR", y="log2(norm count +1)", x="samples") +
    coord_cartesian(ylim=c(0, 10))
ggsave(file.path(pkg_dir, "figures", "utr_deseq2", "five_prime_UTR_all_treatments.pdf"), 
       width=6, height=3)   


comb_df %>% dplyr::filter(region == "5' UTR") %>%
ggplot(aes(x=log2_cnt, fill=treatment)) +
  geom_density(alpha=0.3) +
  theme_minimal() +
  facet_wrap(~ treatment, nrow=4) +
  scale_fill_manual(values=my_color) + 
  coord_cartesian(xlim=c(0, 5))
ggsave(file.path(pkg_dir, "figures", "utr_deseq2", "test.pdf"), width=6, height=3)   

#
# go analysis: I am not sure this is nesseary or provide any meaningful information because
# the relatively low occupancy is global  
#

# find transcripts whose occupancy is not affected by DUX4