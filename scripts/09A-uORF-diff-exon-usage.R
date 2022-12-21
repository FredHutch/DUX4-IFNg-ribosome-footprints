# differential exon usage analysis for uORFs in untreated vs. IFNg or vs. DUX4
# use the excel sheets

library(DESeq2)
library(DEXSeq)
library(plyranges)
library(GenomicAlignments)
library(writexl)
library(readxl)
library(tidyverse)
library(corrr)
library(wesanderson)
library(BSgenome.Hsapiens.UCSC.hg38)
bs_genome <- BSgenome.Hsapiens.UCSC.hg38
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
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "uORFs")
ribotish_dir <- file.path(pkg_dir, "ribotish")
load(file.path(pkg_dir, "data", "p_sites.rda"))
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda")) # use for sizeFactor and colData

sample_info <- as(colData(dds_cds_by_gene)[, c("sample_name", "treatment", "bam_files")], "data.frame") %>%
  dplyr::rename(condition = "treatment") %>%
  mutate(sample_name = str_replace_all(sample_name, "-", "_"))

#
# DEXSeq using exonic parts (that will include 5' UTR)
#

# (A) Ribo-seq
exonic_parts <- exonicParts(txdb, linked.to.single.gene.only = TRUE) # GenomicFeatures
rse <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=exonic_parts, 
                    reads=pc, 
                    singelEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=FALSE)
})
rse_exon <- do.call(cbind, rse)
colData(rse_exon) <- as(sample_info, "DataFrame")
dxd_exonic_parts <- DEXSeqDataSetFromSE(rse_exon, design= ~ sample + exon + condition:exon)
save(dxd_exonic_parts, file=file.path(pkg_dir, "data", "dxd_exonic_parts.rda"))

# (B) RNA-seq
load(file.path(pkg_dir, "data", "rse_cds_mRNA.rda"))
sample_info <- as(colData(rse_cds_mRNA)[, c("sample_name", "treatment", "bam_files")], "data.frame") %>%
  dplyr::rename(condition = "treatment") %>%
  mutate(sample_name = str_replace_all(sample_name, "-", "_"))

rse <- summarizeOverlaps(features=exonic_parts, reads=sample_info$bam_files, 
                         singleEnd=TRUE, inter.feature=FALSE, ignore.strand=TRUE,
                         mode="IntersectionStrict", BPPARAM=bp_param)  
colData(rse) <- as(sample_info, "DataFrame")                         
dxd_exonic_parts_mRNA <- DEXSeqDataSetFromSE(rse, design= ~ sample + exon + condition:exon)
save(dxd_exonic_parts_mRNA, file=file.path(pkg_dir, "data", "dxd_exnoic_parts_mRNA.rda"))

# (C) do dexseq for Ribo-seq exonic parts
load(file.path(pkg_dir, "data", "DUX4_induced_v2.rda")) # LFC > 1
source(file.path(pkg_dir, "scripts", "07A-tools_dexseq.R"))
histone_variants <- .histone_variants(gene.anno)
exclude_gene <- c(histone_variants$gene_id, DUX4_induced_v2$ensembl)
list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))

dxd_exonic_parts <- dxd_exonic_parts[!rowData(dxd_exonic_parts)$groupID %in% c(histone_variants$gene_id, DUX4_induced_v2$ensembl)]
dxd_exonic_parts_mRNA <- dxd_exonic_parts_mRNA[!rowData(dxd_exonic_parts_mRNA)$groupID %in% c(histone_variants$gene_id, DUX4_induced_v2$ensembl)]

.filter_by_rowMeans <- function(rse, mean_filter = 15, treatments) {
     rse <- rse[, rse$treatment %in% treatments]
     rse <- rse[rowMeans(assays(rse)[["counts"]]) >= mean_filter]
}

rse_ribo <- .filter_by_rowMeans(rse=rse_cds_by_gene, treatments = list_comp[["S1"]])
keep_groupID <- rownames(rse_ribo)
dxd_sub <- dxd_exonic_parts[rowData(dxd_exonic_parts)$groupID %in% keep_groupID]
dxd_S2_ribo <- .do_dexseq(dxd_sub, treatments = list_comp[["S2"]])
dxd_S1_ribo <- .do_dexseq(dxd_sub, treatments = list_comp[["S1"]])
dxd_S6_ribo <- .do_dexseq(dxd_sub, treatments = list_comp[["S6"]])
res_S1_ribo <- DEXSeqResults(dxd_S1_ribo)
res_S2_ribo <- DEXSeqResults(dxd_S2_ribo)
res_S6_ribo <- DEXSeqResults(dxd_S6_ribo)

#
# uORFs usages in different treatment
#

# (A) load the ORF results
orf_files <- list.files(ribotish_dir, pattern="*_pred.txt", full.name=TRUE)
orf_list <- bplapply(orf_files, function(x) {
  tb <- read_delim(x, delim = "\t")
  position <- str_split(tb$GenomePos, ":")
  start_end <- str_split(sapply(position, "[[", 2), "-")
  gr <- data.frame(seqnames=sapply(position, "[[", 1),
                   strand=sapply(position, "[[", 3),
                   start=as.numeric(sapply(start_end, "[[", 1)),
                   end=as.numeric(sapply(start_end, "[[", 2))) %>%
    plyranges::as_granges()            
  mcols(gr) <- tb   
  gr %>% plyranges::filter(!duplicated(GenomePos))
})
names(orf_list) <- str_replace(basename(orf_files), "_pred.txt", "")

# (B) what's the dxd_exonic_parts say about uORFs (in 5'UTR) ?
features <- c("5'UTR", "5'UTR:CDSFrameOverlap", "5'UTR:Known")
uORF <- bplapply(names(orf_list), function(x) {
  orf_list[[x]] %>% plyranges::filter(TisType %in% features) %>%
    plyranges::mutate(source=x, overlap_status = "None")
})
names(uORF) <- names(orf_list)


#
# uORF: IFNg vs. untreated; check their exon usage!
#
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
# ? which uORFs are in untreated but not in IFNg?
ol <- findOverlaps(uORF[["untreated"]], uORF[["IFNg"]])
tmp <- uORF[["untreated"]][-queryHits(ol)] # untreated-specific

# ? does those untreated-specific uORFs have different exon usage?
ol <- findOverlaps(dxd_S2_ribo, tmp)
sel <- res_S2_ribo[queryHits(ol), ]
sel <- sel[!is.na(sel$padj) & sel$padj < 0.05, ]
pdf(file.path(fig_dir, "test.pdf"))
for (i in rownames(sel)) {
  plotDEXSeq(res_S2_ribo, sel[i, "groupID"], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
}
dev.off()

#
# Conclusion: ok, this doesn't really work, but get ideas why it is better to just look
# at the translation changes for the transcripts that have missing uORFs comparied
# to untreated
# 
load(file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))

ol <- findOverlaps(uORF[["untreated"]], uORF[["IFNg"]])
uORF[["untreated"]]$overlap_status[unique(queryHits(ol))] <- "overlap"
uORF[["IFNg"]]$overlap_status[unique(subjectHits(ol))] <- "overlap"

# find out their TE in 5'UTR and CDS? ; put on scatter plot; how about the counts?
path <- file.path(pkg_dir, "stats", "translation_changes", 
                  "5UTR_v2", "S2-5UTR-by_tx-exclude-DUX4-indueced.xlsx")
tb <- read_xlsx(path=path, sheet="all") 
rse_sub <- rse_5UTR_by_tx[tb$tx_name]
which_tx <- subsetByOverlaps(rse_sub, uORF[["untreated"]] %>% plyranges::filter(overlap_status == "None"),
                             ignore.strand=FALSE)
tb %>% dplyr::filter(tx_name %in% rownames(which_tx)) %>% 
ggplot(aes(x=rna_lfc, y=ribo_lfc)) +
    geom_point(size=1.5, alpha=0.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5))
ggsave(file.path(pkg_dir, "manuscript", "figures", "uORFs", "S2-5UTR_uORFs_in_untreated_only.pdf"))

# CDS
path <- file.path(pkg_dir, "stats", "translation_changes", 
                  "CDS", "S2_IFNg_vs_untreated.xlsx")
tb <- read_xlsx(path=path, sheet="CDS") 
sub_tb <- tb %>% dplyr::filter(ensembl %in% rowData(which_tx)$gene_id) 

cor(sub_tb$ribo_lfc, sub_tb$rna_lfc)
ggplot(sub_tb, aes(x=rna_lfc, y=ribo_lfc)) +
    geom_point(size=1.5, alpha=0.5) +
    theme_bw() +
    geom_abline(slope=1, intercept=0, linetype="dashed", alpha=0.5, color="gray50") +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5))
ggsave(file.path(pkg_dir, "manuscript", "figures", "uORFs", "S2-CDS_uORFs_in_untreated_only.pdf"))

# 
# IFNg vs. untreated
#
# ok, next (1) merge IFNg and untreated, check the p-sites count on scatter plot to see
# if the untreated really have more counts
# (2) check uORF and mORF activities, and check their TE In CDS! 
#
merge <- c(uORF[["untreated"]], uORF[["IFNg"]]) %>%
  plyranges::filter(!duplicated(GenomePos))

sample_idx <- c(1, 5, 6, 7:9)
uORF_IFNg_untreated <- bplapply(p_sites[sample_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_IFNg_untreated <- do.call(cbind, uORF_IFNg_untreated)
colData(uORF_IFNg_untreated) <- colData(dds_cds_by_gene)[names(p_sites)[sample_idx], ]
dds <- DESeqDataSet(uORF_IFNg_untreated, design=~treatment)
dds <- DESeq(dds)
res <- results(dds, alpha=0.05) 

cnt_df <- as.data.frame(counts(dds, normalized=TRUE)) %>%
  dplyr::mutate(untreated = log2(rowMeans(.[, 1:3])+1), IFNg=log2(rowMeans(.[, 4:6])+1)) %>%
  add_column(overlap_status=rowData(dds)$overlap_status, source=rowData(dds)$source) %>%
  dplyr::mutate(overlap_status = ifelse(overlap_status != "None", overlap_status, source))

ggplot(cnt_df, aes(x=IFNg, y=untreated, color=overlap_status, group=overlap_status, fill=overlap_status)) +
  geom_point(alpha=0.4) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  labs(x="IFNgs mean log2(count+1)", y="untreated mean log2(count+1)") +
  theme(legend.position = "bottom") 

ggsave(file.path(fig_dir, "scatter-untreated-IFNgs-uORF-counts.pdf"))    

# check the translation changes (in CDS) of those up-regulaed uORF in IFNg
up_reg <- as.data.frame(res) %>% 
  add_column(tx_name=rowData(dds)$Tid, symbol=rowData(dds)$Symbol) %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0)

path <- file.path(pkg_dir, "stats", "translation_changes", 
                  "CDS", "S2_IFNg_vs_untreated.xlsx")
tb <- read_xlsx(path=path, sheet="CDS") 
sub_tb <- tb %>% dplyr::filter(gene_name %in% up_reg$symbol) 

cor(sub_tb$ribo_lfc, sub_tb$rna_lfc)
ggplot(sub_tb, aes(x=rna_lfc, y=ribo_lfc)) +
    geom_point(size=1.5, alpha=0.5) +
    theme_bw() +
    geom_abline(slope=1, intercept=0, linetype="dashed", alpha=0.5, color="gray50") +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5))
ggsave(file.path(pkg_dir, "manuscript", "figures", "uORFs", "test.pdf"))

# 
# DUX4+IFNg vs. IFNg
#
merge <- c(uORF[["DOX-pulse_IFNg"]], uORF[["IFNg"]]) %>%
  plyranges::filter(!duplicated(GenomePos))

S6_idx <- c(2:4, 7:9)
uORF_S6 <- bplapply(p_sites[S6_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_S6 <- do.call(cbind, uORF_S6)
colData(uORF_S6) <- colData(dds_cds_by_gene)[names(p_sites)[S6_idx], ]
dds_S6 <- DESeqDataSet(uORF_S6, design=~treatment)
dds_S6 <- DESeq(dds_S6)
res_S6 <- results(dds_S6, alpha=0.05)

cnt_df <- as.data.frame(counts(dds_S6, normalized=TRUE)) %>%
  dplyr::mutate(DUX4_IFNg = log2(rowMeans(.[, 1:3])+1), IFNg=log2(rowMeans(.[, 4:6])+1)) %>%
  add_column(overlap_status=rowData(dds_S6)$overlap_status, source=rowData(dds_S6)$source) %>%
  dplyr::mutate(overlap_status = ifelse(overlap_status != "None", overlap_status, source))

ggplot(cnt_df, aes(x=DUX4_IFNg, y=IFNg, color=overlap_status, group=overlap_status, fill=overlap_status)) +
  geom_point(alpha=0.4) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  labs(x="DUX4+IFNg mean log2(count+1)", y="IFNg mean log2(count+1)") +
  theme(legend.position = "bottom") 

ggsave(file.path(fig_dir, "scatter-S6-DUX4_IFNg-IFNgs-uORF-counts.pdf"))    


# 
# DUX4 vs. untreated
#
merge <- c(uORF[["DOX-pulse"]], uORF[["untreated"]]) %>%
  plyranges::filter(!duplicated(GenomePos))

S1_idx <- c(2:4, 1, 5, 6)
uORF_S1 <- bplapply(p_sites[S1_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_S1 <- do.call(cbind, uORF_S1)
colData(uORF_S1) <- colData(dds_cds_by_gene)[names(p_sites)[S1_idx], ]
dds_S1 <- DESeqDataSet(uORF_S1, design=~treatment)
dds_S1 <- DESeq(dds_S1)
res_S1 <- results(dds_S1, alpha=0.05)

cnt_df <- as.data.frame(counts(dds_S1, normalized=TRUE)) %>%
  dplyr::mutate(DUX4 = log2(rowMeans(.[, 1:3])+1), untreated=log2(rowMeans(.[, 4:6])+1)) %>%
  add_column(overlap_status=rowData(dds_S1)$overlap_status, source=rowData(dds_S1)$source) %>%
  dplyr::mutate(overlap_status = ifelse(overlap_status != "None", overlap_status, source))

ggplot(cnt_df, aes(x=DUX4, y=untreated, color=overlap_status, group=overlap_status, fill=overlap_status)) +
  geom_point(alpha=0.4) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  labs(x="DUX4 mean log2(count+1)", y="untreated mean log2(count+1)") +
  theme(legend.position = "bottom") 

ggsave(file.path(fig_dir, "scatter-S1-DUX4-untreated-uORF-counts.pdf"))    
