# use the results determined by 07B-use_dexseq.R (dxd_S1 and dxd_S6) to 
# Aims: identify translation-specific differentail relative exon usage for S1 and S6 in 
# first exons, Q1, Q2, Q3, and Q4


library(DEXSeq)
library(hg38.HomoSapiens.Gencode.v35)
library(tidyverse)
library(writexl)
library(readxl)
data(gene.anno)
gene.anno <- as.data.frame(gene.anno) 
txdb <- hg38.HomoSapiens.Gencode.v35
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)
library(goseq)
library(plyranges)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "dexseq")
stat_dir <- file.path(pkg_dir, "stats", "differential_relative_exon_usage")

source(file.path(pkg_dir, "scripts", "tools.R")) # .do_goseq()
source(file.path(pkg_dir, "scripts", "07A-tools_dexseq.R"))
load(file.path(pkg_dir, "data", "dxd_S1.rda"))
load(file.path(pkg_dir, "data", "dxd_S6.rda"))
load(file.path(pkg_dir, "data", "DUX4_induced.rda"))
load(file.path(pkg_dir, "data", "IFNg_induced.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

histone_variants <- .histone_variants(gene_anno)
HLA <- gene.anno %>% 
  dplyr::filter(grepl("^HLA-(A|B|C|E|F|G|H)|^MR1|^B2M", gene_name)) %>%
  dplyr::select(gene_name, gene_id)
# there are about 90 PSMs; use the excel sheel from Danielle
psm_name <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "Gene Lists_Proteasome_MHC-I.xlsx"),
                              sheet=1, range="A4:B49") 
PSM <- gene.anno %>% dplyr::select(gene_name, gene_id) %>%
  dplyr::filter(gene_name %in% psm_name$Symbol)

#
# tools
#
.get_Q1_exons <- function(.x) {
    strand <- as.character(.x$genomicData.strand)[1]
    total_length <- sum(.x$genomicData.width)
    #num_q1_exons <- ifelse(strand == "+", sum(cumsum(.x$genomicData.width) / total_length < 0.25),
    #                                      sum(cumsum(rev(.x$genomicData.width)) / total_length < 0.25))

    if (strand %in% c("+", "*"))  {
        num_q1_exons <- sum(cumsum(.x$genomicData.width) / total_length < 0.25)
        q1_exons <- .x$featureID[1:num_q1_exons]
    }                  

    if (strand == "-") {
        n_exons <- nrow(.x)
        num_q1_exons <- sum(cumsum(rev(.x$genomicData.width)) / total_length < 0.25)
        q1_exons <- .x$featureID[seq(nrow(.x), by=-1, length.out=num_q1_exons)]
    }

    .x %>% dplyr::filter(featureID %in% q1_exons)
}

#
# S1
#
res_S1_ribo <- DEXSeqResults(dxd_S1$ribo)
res_S1_rna <- DEXSeqResults(dxd_S1$rna)

# (1) get Q1 exons
Q1_exons_ribo_S1 <- as.data.frame(res_S1_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() 
  # group_modify(~.get_Q1_exons(.x)) %>% ungroup() 

Q1_exons_rna_S1 <- as.data.frame(res_S1_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() %>%
  # group_modify(~.get_Q1_exons(.x)) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, DOX_pulse, log2fold_DOX_pulse_untreated)

# (2) candidates: for each groupID, determine 
deu_Q1_S1 <- Q1_exons_ribo_S1 %>% 
  left_join(Q1_exons_rna_S1, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_DOX_pulse_untreated.ribo < -1,
                rna_deu_down = padj.rna < 0.05 & log2fold_DOX_pulse_untreated.rna < -1,
                ribo_deu_up = padj.ribo < 0.05 & log2fold_DOX_pulse_untreated.ribo > 1, 
                rna_deu_up = padj.rna < 0.05 & log2fold_DOX_pulse_untreated.rna > 1) %>%
  dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
  dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)


deu_Q1_S1 %>% distinct(gene_name) %>% count()
deu_Q1_S1 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q1_S1 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

candidate_Q1_S1 <- deu_Q1_S1 %>% dplyr::filter(deu_down) 
candidate_Q1_S1_up <- deu_Q1_S1 %>% dplyr::filter(deu_up)
writexl::write_xlsx(list(down=candidate_Q1_S1, up=candidate_Q1_S1_up), 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                         "S1_Diff_Q1Exon_usage_candidates.xlsx"))

# (C) GO term analysis
universal <- deu_Q1_S1 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
selected <- candidate_Q1_S1 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_S1 <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene)  %>%
                    arrange(padj)
writexl::write_xlsx(enriched_go_S1, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                   "S1_Diff_Q1Exon_usage_candidates_GO.xlsx"))                         

enriched_go_S1 %>% dplyr::mutate(deu_to_category_ratio=numDEInCat/numInCat, log10adjpval=-10*log10(padj)) %>%
  arrange(desc(padj)) %>% dplyr::mutate(term=factor(term, levels=term)) %>%
  ggplot(aes(y=term, x=log10adjpval)) +
    geom_bar(stat="identity", width=0.2) +
    geom_point(aes(size=deu_to_category_ratio)) +
    theme_bw() + theme(legend.position="bottom") + labs(x="-10Log10(FDR)")
ggsave(file.path(fig_dir, "S1_Diff_Q1Exon_usage_candidates_GO.pdf"))    


# viz: sanity check
pdf(file.path(fig_dir, "S1-Q1-SLC6A6-ENSG00000131389.17-ribo.pdf"))
plotDEXSeq(res_S1_ribo, "ENSG00000131389.17", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf(file.path(fig_dir, "S1-Q1-SLC6A6-ENSG00000131389.17-rna.pdf"))
plotDEXSeq(res_S1_rna, "ENSG00000131389.17", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

#
# S6
#
histone_variants <- .histone_variants(gene_anno)
res_S6_ribo <- DEXSeqResults(dxd_S6$ribo)
res_S6_rna <- DEXSeqResults(dxd_S6$rna)

# (1) get Q1 exons
Q1_exons_ribo_S6 <- as.data.frame(res_S6_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() 
  # group_modify(~.get_Q1_exons(.x)) %>% ungroup() 

Q1_exons_rna_S6 <- as.data.frame(res_S6_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  #group_modify(~.get_Q1_exons(.x)) %>% ungroup() %>%
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() %>%
  dplyr::select(name, padj, IFNg, DOX_pulse_IFNg, log2fold_DOX_pulse_IFNg_IFNg)

deu_Q1_S6 <- Q1_exons_ribo_S6 %>% 
  left_join(Q1_exons_rna_S6, by="name", suffix=c(".ribo", ".rna")) %>%
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

deu_Q1_S6 %>% distinct(gene_name) %>% count()
deu_Q1_S6 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q1_S6 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

# get candidate 
candidate_Q1_S6 <- deu_Q1_S6 %>% dplyr::filter(deu_down)
candidate_Q1_S6_up <- deu_Q1_S6 %>% dplyr::filter(deu_up)
writexl::write_xlsx(list(down=candidate_Q1_S6, up=candidate_Q1_S6_up),
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                         "S6_Diff_Q1Exon_usage_candidates.xlsx"))
sum(unique(candidate_Q1_S6$gene_name) %in% unique(candidate_Q1_S1$gene_name))                         

# (C) GO term analysis
universal <- deu_Q1_S6 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
selected <- candidate_Q1_S6 %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_S6 <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene) %>% arrange(padj)
writexl::write_xlsx(enriched_go_S6, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                   "S6_Diff_Q1Exon_usage_candidates_GO.xlsx"))

enriched_go_S6 %>% dplyr::mutate(deu_to_category_ratio=numDEInCat/numInCat, log10adjpval=-10*log10(padj)) %>%
  arrange(desc(padj)) %>%
  dplyr::mutate(term=factor(term, levels=term)) %>%
  ggplot(aes(y=term, x=log10adjpval)) +
    geom_bar(stat="identity", width=0.2) +
    geom_point(aes(size=deu_to_category_ratio)) +
    theme_bw() + theme(legend.position="bottom") + labs(x="-10Log10(FDR)")
ggsave(file.path(fig_dir, "S6_Diff_Q1Exon_usage_candidates_GO.pdf"), height=8)  

# (D) Specificity in HLA, PSM, IFNg-stimulated?
Q1_HLA <- deu_Q1_S6 %>% dplyr::filter(gene_id %in% HLA$gene_id) %>%
  add_column(disc = "HLA")
# HLA-E and HLA-B (lfc =-1.15 and FDR=0.05008)
Q1_HLA %>% distinct(gene_name) %>% count()
Q1_HLA %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()

Q1_PSM <- deu_Q1_S6 %>% dplyr::filter(gene_id %in% PSM$gene_id) %>%
  add_column(disc = "PSM")
Q1_PSM %>% distinct(gene_name)
Q1_PSM %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
# 8/40 = 20% show DEU in Q1

Q1_IFNg <- deu_Q1_S6 %>% dplyr::filter(IFNg_induced) %>%
  add_column(disc = "IFNg-stimulated")
Q1_IFNg %>% distinct(gene_name) %>% count()
Q1_IFNg %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
# 30/247 = 12%

writexl::write_xlsx(list(Q1_HLA=Q1_HLA, Q1_PSM=Q1_PSM, Q1_IFNg=Q1_IFNg),
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                   "S6_Diff_Q1Exon_HLA_PSM_IFNg.xlsx"))

# (E) misc. statistics testings
half_quantile <- quantile(deu_Q1_S6$log2fold_DOX_pulse_IFNg_IFNg.ribo, na.rm=TRUE)[3]
Q1_HLA %>% bind_rows(Q1_PSM) %>% bind_rows(Q1_IFNg) %>%
  ggplot(aes(x=disc, y=log2fold_DOX_pulse_IFNg_IFNg.ribo)) +
    geom_violin() +
    geom_dotplot(binaxis='y', binwidth=0.2,
                 stackdir='center', dotsize=0.5, 
                 aes(color=flag, fill=flag), 
                 alpha=0.7) +
    geom_hline(yintercept=half_quantile, linetype="dashed", colour="red") +
    theme_minimal()
ggsave(file.path(fig_dir, "S6-Q1-violin-LFC-PSM-HLA-IFNg.pdf"))

# (E) misc stuff
cnt <- candidate_Q1_S6 %>% count(gene_id)
table(cnt$n)
a=cnt %>% dplyr::filter(n > 1) %>% pull(gene_id)
tmp <- candidate_Q1_S6 %>% dplyr::filter(gene_id %in% a)


writexl::write_xlsx(candidate_Q1_S6, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                         "S6_Diff_Q1Exon_usage_candidates.xlsx"))

# (F) viz: sanity check

# PSMB9 (lfc < -0.7)
pdf(file.path(fig_dir, "S6-Q1-PSMB9-ENSG00000240065.8-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000240065.8", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-PSMB9-ENSG00000240065.8-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000240065.8", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()                         

# ACADVL: exon 3 and 5
pdf(file.path(fig_dir, "S6-Q1-ACADVL-ENSG00000072778.20-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000072778.20", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-ACADVL-ENSG00000072778.20-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000072778.20", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

# PSME1: ENSG00000092010.15
pdf(file.path(fig_dir, "S6-Q1-PSME1-ENSG00000092010.15-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000092010.15", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-PSME1-ENSG00000092010.15-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000092010.15", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

# PSME2: ENSG00000100911.16
pdf(file.path(fig_dir, "S6-Q1-PSME2-ENSG00000100911.16-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000100911.16", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-PSME2-ENSG00000100911.16-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000100911.16", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

# PSMC2 ENSG00000161057.12
pdf(file.path(fig_dir, "S6-Q1-PSMC2-ENSG00000161057.12-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000161057.12", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-PSMC2-ENSG00000161057.12-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000161057.12", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

# PSMF1: dose not have DEU
pdf(file.path(fig_dir, "S6-Q1-PSMF1-ENSG00000125818.18-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000125818.18", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()
pdf(file.path(fig_dir, "S6-Q1-PSMF1-ENSG00000125818.18-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000125818.18", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

# ENSG00000143190.23    POU2F1

# sanity check:
.x <- as.data.frame(res_S6_ribo) %>% dplyr::filter(groupID=="ENSG00000169248.13")

#######################
#
# S6: Q2, Q3, and Q4 
#
#######################

# test "+": ENSG00000131389.17
.x <- as.data.frame(res_S6_ribo) %>% dplyr::filter(groupID == "ENSG00000131389.17")
q1x <- .get_exons_by_quantile(.x, n_quantile = 1, type="within")
q2x <- .get_exons_by_quantile(.x, n_quantile = 2, type="any")

# test: ENSG00000240065.8
.x <- as.data.frame(res_S6_ribo) %>% dplyr::filter(groupID == "ENSG00000240065.8")
q1x <- .get_exons_by_quantile(.x, n_quantile = 1, type="within")
q2x <- .get_exons_by_quantile(.x, n_quantile = 2, type="any")

# test "-": ENSG00000205220.12
.x <- as.data.frame(res_S6_ribo) %>% dplyr::filter(groupID == "ENSG00000205220.12")
q1x <- .get_exons_by_quantile(.x, n_quantile = 1, type="within")
q2x <- .get_exons_by_quantile(.x, n_quantile = 2, type="any")


#
# S6:Q2
#
Q2_exons_ribo_S6 <- as.data.frame(res_S6_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=2, type="any")) %>% ungroup() 

Q2_exons_rna_S6 <- as.data.frame(res_S6_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=2, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, IFNg, DOX_pulse_IFNg, log2fold_DOX_pulse_IFNg_IFNg)

deu_Q2_S6 <- Q2_exons_ribo_S6 %>% 
  left_join(Q2_exons_rna_S6, by="name", suffix=c(".ribo", ".rna")) %>%
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


deu_Q2_S6 %>% dplyr::filter(deu_down) %>% distinct(gene_id) %>% count()
deu_Q2_S6 %>% distinct(gene_id) %>% count() 

deu_Q2_S6 %>% dplyr::filter(deu_up) %>% distinct(gene_id) %>% count()
deu_Q2_S6 %>% distinct(gene_id) %>% count() 

#
# S6:Q3
#
Q3_exons_ribo_S6 <- as.data.frame(res_S6_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=3, type="any")) %>% ungroup() 

Q3_exons_rna_S6 <- as.data.frame(res_S6_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=3, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, IFNg, DOX_pulse_IFNg, log2fold_DOX_pulse_IFNg_IFNg)

deu_Q3_S6 <- Q3_exons_ribo_S6 %>% 
  left_join(Q3_exons_rna_S6, by="name", suffix=c(".ribo", ".rna")) %>%
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


deu_Q3_S6 %>% dplyr::filter(deu_down) %>% distinct(gene_id) %>% count()
deu_Q3_S6 %>% distinct(gene_id) %>% count() 

deu_Q3_S6 %>% dplyr::filter(deu_up) %>% distinct(gene_id) %>% count()
deu_Q3_S6 %>% distinct(gene_id) %>% count() 

#
# S6:Q4
#
Q4_exons_ribo_S6 <- as.data.frame(res_S6_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=4, type="any")) %>% ungroup() 

Q4_exons_rna_S6 <- as.data.frame(res_S6_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=4, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, IFNg, DOX_pulse_IFNg, log2fold_DOX_pulse_IFNg_IFNg)

deu_Q4_S6 <- Q4_exons_ribo_S6 %>% 
  left_join(Q4_exons_rna_S6, by="name", suffix=c(".ribo", ".rna")) %>%
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


deu_Q4_S6 %>% dplyr::filter(deu_down) %>% distinct(gene_id) %>% count()
deu_Q4_S6 %>% distinct(gene_id) %>% count() 

deu_Q4_S6 %>% dplyr::filter(deu_up) %>% distinct(gene_id) %>% count()
deu_Q4_S6 %>% distinct(gene_id) %>% count() 



#
# save S1 and S6
#
Q1_exons_S1 <- list(ribo=Q1_exons_ribo_S1, rna=Q1_exons_rna_S1)
Q2_exons_S1 <- list(ribo=Q2_exons_ribo_S1, rna=Q2_exons_rna_S1)
Q3_exons_S1 <- list(ribo=Q3_exons_ribo_S1, rna=Q3_exons_rna_S1)
Q4_exons_S1 <- list(ribo=Q4_exons_ribo_S1, rna=Q4_exons_rna_S1)
S1_relative_exons_usage <- list(Q1=Q1_exons_S1, Q2=Q2_exons_S1, Q3=Q3_exons_S1, Q4=Q4_exons_S1)
save(S1_relative_exons_usage, file=file.path(pkg_dir, "data", "S1_relative_exons_usage.rda"))


Q1_exons_S6 <- list(ribo=Q1_exons_ribo_S6, rna=Q1_exons_rna_S6)
Q2_exons_S6 <- list(ribo=Q2_exons_ribo_S6, rna=Q2_exons_rna_S6)
Q3_exons_S6 <- list(ribo=Q3_exons_ribo_S6, rna=Q3_exons_rna_S6)
Q4_exons_S6 <- list(ribo=Q4_exons_ribo_S6, rna=Q4_exons_rna_S6)
S6_relative_exons_usage <- list(Q1=Q1_exons_S6, Q2=Q2_exons_S6, Q3=Q3_exons_S6, Q4=Q4_exons_S6)
save(S6_relative_exons_usage, file=file.path(pkg_dir, "data", "S6_relative_exons_usage.rda"))


###############
#
# S2: Q2/Q3/Q4
#
###############
res_S2_ribo <- DEXSeqResults(dxd_S2$ribo)
res_S2_rna <- DEXSeqResults(dxd_S2$rna)

#
# S2: Q1
#
Q1_exons_ribo_S2 <- as.data.frame(res_S2_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() 

Q1_exons_rna_S2 <- as.data.frame(res_S2_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=1, type="within")) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, IFNg, log2fold_IFNg_untreated)

deu_Q1_S2 <- Q1_exons_ribo_S2 %>% 
  left_join(Q1_exons_rna_S2, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo < -1,
                rna_deu_down = padj.rna < 0.05 & log2fold_IFNg_untreated.rna < -1,
                ribo_deu_up = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo > 1, 
                rna_deu_up = padj.rna < 0.05 & log2fold_IFNg_untreated.rna > 1) %>%
  dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
  dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)

deu_Q1_S2 %>% distinct(gene_name) %>% count()
deu_Q1_S2 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q1_S2 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

#
# S2: Q2
#
Q2_exons_ribo_S2 <- as.data.frame(res_S2_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=2, type="any")) %>% ungroup() 

Q2_exons_rna_S2 <- as.data.frame(res_S2_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=2, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, IFNg, log2fold_IFNg_untreated)

deu_Q2_S2 <- Q2_exons_ribo_S2 %>% 
  left_join(Q2_exons_rna_S2, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo < -1,
                rna_deu_down = padj.rna < 0.05 & log2fold_IFNg_untreated.rna < -1,
                ribo_deu_up = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo > 1, 
                rna_deu_up = padj.rna < 0.05 & log2fold_IFNg_untreated.rna > 1) %>%
  dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
  dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)

deu_Q2_S2 %>% distinct(gene_name) %>% count()
deu_Q2_S2 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q2_S2 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

#
# S2: Q3
#
Q3_exons_ribo_S2 <- as.data.frame(res_S2_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=3, type="any")) %>% ungroup() 

Q3_exons_rna_S2 <- as.data.frame(res_S2_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=3, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, IFNg, log2fold_IFNg_untreated)

deu_Q3_S2 <- Q3_exons_ribo_S2 %>% 
  left_join(Q3_exons_rna_S2, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo < -1,
                rna_deu_down = padj.rna < 0.05 & log2fold_IFNg_untreated.rna < -1,
                ribo_deu_up = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo > 1, 
                rna_deu_up = padj.rna < 0.05 & log2fold_IFNg_untreated.rna > 1) %>%
  dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
  dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)

deu_Q3_S2 %>% distinct(gene_name) %>% count()
deu_Q3_S2 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q3_S2 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

#
# S2: Q4
#
Q4_exons_ribo_S2 <- as.data.frame(res_S2_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=4, type="any")) %>% ungroup() 

Q4_exons_rna_S2 <- as.data.frame(res_S2_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_exons_by_quantile(.x, n_quantile=4, type="any")) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, IFNg, log2fold_IFNg_untreated)

deu_Q4_S2 <- Q4_exons_ribo_S2 %>% 
  left_join(Q4_exons_rna_S2, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu_down = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo < -1,
                rna_deu_down = padj.rna < 0.05 & log2fold_IFNg_untreated.rna < -1,
                ribo_deu_up = padj.ribo < 0.05 & log2fold_IFNg_untreated.ribo > 1, 
                rna_deu_up = padj.rna < 0.05 & log2fold_IFNg_untreated.rna > 1) %>%
  dplyr::mutate(deu_down = if_else(ribo_deu_down & (!rna_deu_down), TRUE, FALSE)) %>%
  dplyr::mutate(deu_up = if_else(ribo_deu_up & (!rna_deu_up), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)

deu_Q4_S2 %>% distinct(gene_name) %>% count()
deu_Q4_S2 %>% dplyr::filter(deu_down) %>% distinct(gene_name) %>% count()
deu_Q4_S2 %>% dplyr::filter(deu_up) %>% distinct(gene_name) %>% count()

#
# save the results
#
Q1_exons_S2 <- list(ribo=Q1_exons_ribo_S2, rna=Q1_exons_rna_S2)
Q2_exons_S2 <- list(ribo=Q2_exons_ribo_S2, rna=Q2_exons_rna_S2)
Q3_exons_S2 <- list(ribo=Q3_exons_ribo_S2, rna=Q3_exons_rna_S2)
Q4_exons_S2 <- list(ribo=Q4_exons_ribo_S2, rna=Q4_exons_rna_S2)
S2_relative_exons_usage <- list(Q1=Q1_exons_S2, Q2=Q2_exons_S2, Q3=Q3_exons_S2, Q4=Q4_exons_S2)
save(S2_relative_exons_usage, file=file.path(pkg_dir, "data", "S2_relative_exons_usage.rda"))