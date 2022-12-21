# Differentail exon usages (DEXSeq)
# use the results determined by 07B-use_dexseq.R (dxd_S1 and dxd_S6) to 
# (1) identify translation-specific differentail relative 1st exon usage for S1 and S6 
# (2) infer whether these 1st exon usage are not random: compare results of S1 and S6

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

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "dexseq")
stat_dir <- file.path(pkg_dir, "stats", "dexseq")

source(file.path(pkg_dir, "scripts", "tools.R")) # .do_goseq(), .histone_variants()
load(file.path(pkg_dir, "data", "dxd_S1.rda"))
load(file.path(pkg_dir, "data", "dxd_S6.rda"))
load(file.path(pkg_dir, "data", "DUX4_induced.rda"))
load(file.path(pkg_dir, "data", "IFNg_induced.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

sig_down_S1 <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "5UTR_S1_down_occupancy.xlsx"))
sig_down_S6 <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "5UTR_S6_down_occupancy.xlsx"))

tss_down_S1 <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "S1_TSS_down_occupancy.xlsx"))
tss_no_suppress <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "S1_TSS_not_suppressed_in_DUX4.xlsx"))

#
# tools
# 
.get_1st_exon <- function(.x){
  strand <- as.character(.x$genomicData.strand)[1]
  exon1 <- ifelse(strand=="+",  .x$featureID[1], .x$featureID[nrow(.x)])
  .x %>% dplyr::filter(featureID %in% exon1)
}


#
# S1: differential relative 1st exon usage and compare to that of RNA-seq
#
histone_variants <- .histone_variants(gene.anno)
res_S1_ribo <- DEXSeqResults(dxd_S1$ribo)
res_S1_rna <- DEXSeqResults(dxd_S1$rna)

first_exon_ribo_S1 <- as.data.frame(res_S1_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_1st_exon(.x)) %>% ungroup() 

first_exon_rna_S1 <- as.data.frame(res_S1_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_1st_exon(.x)) %>% ungroup() %>%
  dplyr::select(name, padj, untreated, DOX_pulse, log2fold_DOX_pulse_untreated)

deu_1st_S1 <- first_exon_ribo_S1 %>% 
  left_join(first_exon_rna_S1, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu = padj.ribo < 0.05 & log2fold_DOX_pulse_untreated.ribo < -1,
                rna_deu = padj.rna < 0.05 & log2fold_DOX_pulse_untreated.rna < -1) %>%
  dplyr::mutate(flag = if_else(ribo_deu & (!rna_deu), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)

# extract candidates
candidate_S1 <- deu_1st_S1 %>% dplyr::filter(flag) 
writexl::write_xlsx(candidate_S1, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                  "S1_Diff_1stExon_usage_candidates.xlsx"))

# misc. analysis
# how many candidate show loss of occupancy in 5' UTR (differentially down-regulated)
sig_down_S1 <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "5UTR_S1_down_occupancy.xlsx"))
tss_down_S1 <- read_xlsx(file.path(pkg_dir, "stats", "differential_analysis", "S1_TSS_down_occupancy.xlsx"))

sum(candidate_S1$gene_id %in% sig_down_S1$gene_id)
sum(candidate_S1$gene_id %in% tss_down_S1$gene_id)
sum(candidate_S1$gene_id %in% tss_no_suppress$gene_id)
sum(sig_down_S1$tx_name %in% tss_down_S1$tx_name)


# GO terms analysis
library(goseq)
source(file.path(pkg_dir, "scripts", "tools.R")) # .do_goseq()
universal <- deu_1st_S1 %>% pull(gene_id) %>% str_replace("\\..*", "")
selected <- candidate_S1 %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_S1 <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene)

writexl::write_xlsx(enriched_go_S1, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                   "S1_Diff_1stExon_usage_candidates_GO.xlsx"))

#
# S6
#
load(file.path(pkg_dir, "data", "dxd_S6.rda"))
res_S6_ribo <- DEXSeqResults(dxd_S6$ribo)
res_S6_rna <- DEXSeqResults(dxd_S6$rna)

# sanity check: PSMB9, check third exon
as.data.frame(res_S6_ribo) %>% dplyr::filter(groupID=="ENSG00000240065.8")
as.data.frame(res_S6_rna) %>% dplyr::filter(groupID=="ENSG00000240065.8")

first_exon_ribo_S6 <- as.data.frame(res_S6_ribo) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_1st_exon(.x)) %>% ungroup() 

first_exon_rna_S6 <- as.data.frame(res_S6_rna) %>% rownames_to_column(var="name") %>%
  dplyr::filter(!groupID %in% histone_variants$gene_id) %>%
  group_by(groupID) %>% 
  group_modify(~.get_1st_exon(.x)) %>% ungroup() %>%
  dplyr::select(name, padj, IFNg, DOX_pulse_IFNg, log2fold_DOX_pulse_IFNg_IFNg)

deu_1st_S6 <- first_exon_ribo_S6 %>% 
  left_join(first_exon_rna_S6, by="name", suffix=c(".ribo", ".rna")) %>%
  dplyr::mutate(ribo_deu = padj.ribo < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.ribo < -1,
                rna_deu = padj.rna < 0.05 & log2fold_DOX_pulse_IFNg_IFNg.rna < -1) %>%
  dplyr::mutate(flag = if_else(ribo_deu & (!rna_deu), TRUE, FALSE)) %>%
  dplyr::mutate(DUX4_induced = groupID %in% DUX4_induced$ensembl) %>%
  dplyr::mutate(IFNg_induced = groupID %in% IFNg_induced$ensembl) %>%
  rename(gene_id="groupID") %>%
  left_join(dplyr::select(gene.anno, gene_id, gene_name), by="gene_id") %>%
  dplyr::relocate(gene_name, .after="gene_id") %>%
  arrange(padj.ribo)


# get candidate
candidate_S6 <- deu_1st_S6 %>% dplyr::filter(flag)

write_xlsx(candidate_S6, 
           path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                         "S6_Diff_1stExon_usage_candidates.xlsx"))

deu_1st <- list(S1=deu_1st_S1, S6=deu_1st_S6)
save(deu_1st, file=file.path(pkg_dir, "data", "deu_1st.rda"))

# how about IFNg-induced
candidate_S6 %>% dplyr::filter(IFNg_induced)
# viz: CXCL11

pdf(file.path(fig_dir, "S6-CXCL11-ENSG00000169248.13-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000169248.13", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf(file.path(fig_dir, "S6-CXCL11-ENSG00000169248.13-rna.pdf"))
plotDEXSeq(res_S6_rna, "ENSG00000169248.13", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf(file.path(fig_dir, "S6-PRAMEF10-ENSG00000187545.5-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000187545.5", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()
pdf(file.path(fig_dir, "S6-HLA-B-ENSG00000234745.11-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000234745.11", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()

pdf(file.path(fig_dir, "S6-HLA-A-ENSG00000206503.13-ribo.pdf"))
plotDEXSeq(res_S6_ribo, "ENSG00000206503.13", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()

# check with TSS/5'UTR suppressed tx:
sum(candidate_S6$gene_id %in% sig_down_S6$gene_id) # 108/411

# GO term analysis
universal <- deu_1st_S6 %>% pull(gene_id) %>% str_replace("\\..*", "")
selected <- candidate_S6 %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_S6 <- .do_goseq(universal=universal, selected_genes=selected,  p_value=0.05,
                         return.DEInCat=TRUE, dds=dds_cds_by_gene)

writexl::write_xlsx(enriched_go_S6, 
                    path=file.path(pkg_dir, "stats", "differential_relative_exon_usage", 
                                   "S6_Diff_1stExon_usage_candidates_GO.xlsx"))


#
# (3) compare S1 and S6; Diff_1stExon_usage_candidates candidate_S6 vs candidate_S1; scatter plot
#
deu_1st_S1 %>% dplyr::filter(gene_name=="HLA-B") %>% as.data.frame()
deu_1st_S6 %>% dplyr::filter(gene_name=="HLA-B") %>% as.data.frame()

deu_1st_S1 %>% dplyr::filter(gene_name=="COL1A1") %>% as.data.frame()
deu_1st_S6 %>% dplyr::filter(gene_name=="COL1A1") %>% as.data.frame()

#
# (3a) any overlaps / venn diagram
#
sum(candidate_S1$name %in% candidate_S6$name) #140 
library(VennDiagram)
venn.diagram(
  x = list(candidate_S1$name, candidate_S6$name),
  category.names = c("S1", "S6"),
  filename = file.path(fig_dir, "S1_vs_S6_candidates_venn_diagrame.png"),
  output=TRUE
)
#library(ggVennDiagram)

#
# (3b) scatter plot and correlation 
#
comb_candidate_S1 <- candidate_S1 %>% 
  dplyr::select(gene_id, gene_name, name, featureID, padj.ribo, log2fold_DOX_pulse_untreated.ribo, flag) %>%
  dplyr::left_join(select(deu_1st_S6, name, padj.ribo, log2fold_DOX_pulse_IFNg_IFNg.ribo, IFNg_induced, flag), by="name", suffix=c(".S1", ".S6")) %>%
  dplyr::mutate(overlap= flag.S1 & flag.S6)

cor(comb_candidate_S1$log2fold_DOX_pulse_IFNg_IFNg.ribo, comb_candidate_S1$log2fold_DOX_pulse_untreated.ribo)

ggplot(comb_candidate_S1, aes(x=log2fold_DOX_pulse_untreated.ribo, log2fold_DOX_pulse_IFNg_IFNg.ribo)) +
  geom_point(aes(color=overlap)) + 
  geom_abline(slope=1, intercept=0, color="gray50", linetype="dashed") +
  geom_smooth(method='lm', se = FALSE) +
  theme_minimal() +
  coord_cartesian(xlim = c(-5, 1), ylim=c(-5, 1.2)) +
  annotate("text", x=-1, y=-4, hjust=0, vjust=-1, label="cor =  0.3") +
  labs(title="S1 candidates: 1st exon usage")
ggsave(file.path(fig_dir, "S1_candidate_scatter_S6.pdf"), height=4, width=5)  

# what are the overlapped 

# what's in S1 GO term and what's in S6 GO term?
tmp1 <- enriched_go_S1 %>% pull(DEInCat) %>% map(str_split, ",") %>% unlist() %>% unique()
length(tmp1)

tmp6 <- enriched_go_S6 %>% pull(DEInCat) %>% map(str_split, ",") %>% unlist() %>% unique()
length(tmp6)