# 10B-manuscript-figures-others.R
# 1) GO terms in three genomic features: 5' UTR, TSS and 1st coding exon in S6
# 2) select "super-parent" GO terms
# 3) try different ways to present the data: dot plots
# what's the idea we try to convey here?

library(tidyverse)
library(GO.db)
library(readxl)
library(BiocParallel)
bp_param <- MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "top-GO-translation-changes-for-down-reg")
source(file.path(pkg_dir, "scripts", "tools.R"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

ls("package:GO.db")

select(GO.db, keys="GO:0016070", column="TERM", keytype="GOID")

#
# tools
#
BP_parent <- as.list(GOBPPARENTS)
BP_children <- as.list(GOBPCHILDREN)
BP_offspring <- as.list(GOBPOFFSPRING)
ntop <- 50

.get_genealogy <- function(category, ntop=50) {
  genealogy <- lapply(category, function(x) {
    parent <- BP_parent[[x]]
    children <- BP_children[[x]] 
    offspring <- BP_offspring[[x]]
    term <- mapIds(GO.db, keys=x, column="TERM", keytype="GOID", multiVals="first")
    data.frame(term=term, 
               parent_in_category = paste0(parent[parent %in% category], collapse=","),
               parent_in_top = any(parent %in% category[1:ntop]),
               number_of_parent = length(parent[parent %in% category]),
               children_in_category = paste0(children[children %in% category], collapse=","),
               children_in_top = any(children %in% category[1:ntop]),
               number_of_children = length(children[children %in% category]),
               offspring_in_category = paste0(offspring[offspring %in% category], collapse=","),
               number_of_offspring = length(offspring[offspring %in% category]))
   
  })
  names(genealogy) <- category
  genealogy <- do.call(rbind, genealogy) %>%
    rownames_to_column(var="category")
}

.parent_GO_barplot <- function(parent_go, title) {
    parent_go %>% dplyr::mutate(neglog10BH = -log10(padj)) %>%
      dplyr::arrange(neglog10BH) %>%
      dplyr::mutate(term=factor(term, levels=term)) %>%
      ggplot(aes(x=term, y=neglog10BH)) +
        geom_bar(stat="identity", width=0.7) +
        theme_bw() +
        coord_flip() +
        labs(y=TeX(r'(-$\log_{10}(padj))'), title=title) +
        theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=9), 
              panel.grid.major.y=element_blank(),
              axis.title.x = element_text(size=9),
              panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
}

#
# (A) acquire the top 50 terms and only keep the children
# 


# TSS
ntop <- 50
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
S6_TSS_tb <- read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
genealogy_tss <- .get_genealogy(category=S6_TSS_tb$category, ntop=ntop) 
sub_tss <- tss_genealogy[1:ntop, ]%>% dplyr::filter(!children_in_top)
childgo_tss <- S6_TSS_tb %>% left_join(sub_tss, by="category") 


# 5UTR
ntop <- 50
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
S6_5UTR_tb <- read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
genealogy_5utr <- .get_genealogy(category=S6_5UTR_tb$category, ntop=ntop) 
sub_5utr <- genealogy_5utr[1:ntop, ]%>% dplyr::filter(!children_in_top)

# 1st coding exon
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
S6_1stE_tb <- read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
ntop <- min(50, nrow(S6_1stE_tb))
genealogy_1stExon <- .get_genealogy(category=S6_1stE_tb$category, ntop=ntop)      
sub_1stExon <- genealogy_1stExon[1:ntop, ]%>% dplyr::filter(!children_in_top)


childgo_tss <- S6_TSS_tb %>% right_join(sub_tss, by="category") 
childgo_5utr <- S6_5UTR_tb %>% right_join(sub_5utr, by="category")
childgo_1stExon <- S6_1stE_tb %>% right_join(sub_1stExon, by="category") 
union_children <- union(childgo_5utr$category, childgo_tss$category) %>% union(childgo_1stExon$category) %>% unique()

# do the go analysis again to obtain the pval
.get_all_go_term <- function(path, feature) {
   tb <- read_xlsx(path=path, sheet="translation_downregulation")
   selected <- tb %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
   go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=1,
                        return.DEInCat=FALSE, dds=dds_cds_by_gene) %>%
                #dplyr::filter(category %in% union_children) %>%
                dplyr::mutate(feature = feature) 
}

universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
S6_1stE_go <- .get_all_go_term(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"), "1st coding exon")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
S6_5UTR_go <- .get_all_go_term(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"), "5' UTR")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
S6_TSS_go <- .get_all_go_term(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"), "TSS")

# prepare for the figure
# (1) flter out the processes that have values in 1st Exon (2) arrange by p-value of 5' UTR
keep_lv <- S6_1stE_go %>% dplyr::filter(category %in% union_children) %>%
  left_join(S6_5UTR_go, by="category") %>% 
  dplyr::filter(category != "GO:0002479") %>% # take out MHC class I, TAP-dependent
  arrange(desc(padj.x)) %>% dplyr::select(category, term.x)

# remove some terms: "cellular protein metabolic process: and "cellular component organization or biogenesis" "MHC class I, TAP-dependent"
keep_lv2 <- S6_1stE_go %>% dplyr::filter(category %in% union_children) %>%
  left_join(S6_5UTR_go, by="category") %>% 
  left_join(S6_TSS_go, by="category") %>% 
  dplyr::filter(!category %in% c("GO:0002479", "GO:0044267", "GO:0071840")) %>% # take out MHC class I, TAP-dependent
  dplyr::mutate(avg_logp = -(log10(padj.x)+ log10(padj.y) + log10(padj) / 3)) %>%
  arrange(avg_logp) %>% dplyr::select(category, term.x)

#tmp_term <- select(GO.db, keys=keep_lv$category, column="TERM", keytype="GOID") 

df <- S6_1stE_go %>%
  add_row(S6_5UTR_go) %>% add_row(S6_TSS_go) %>%
  dplyr::filter(category %in% keep_lv2$category) %>%
  dplyr::mutate(term = factor(term, levels=keep_lv2$term.x)) %>%
  dplyr::mutate(log10padj = -log10(padj), 
                overrepresented=ifelse(padj < 0.05, "YES", "NO")) %>%
  dplyr::mutate(feature=factor(feature, levels=c("5' UTR", "TSS", "1st coding exon")))

library(latex2exp)
ggplot(df, aes(y=term, x=feature)) +
  geom_point(aes(size=log10padj/2, color=overrepresented, fill=overrepresented)) +
  theme_bw() +
  labs(size=TeX("$-\\log_{10}(padj)$")) +
  scale_fill_manual(name ="over-represented",  values=c("ivory4", "firebrick4")) +
  scale_color_manual(name="over-represented",  values=c("ivory4", "firebrick4")) +
  theme(legend.position="bottom", legend.box = "vertical", 
        axis.title.y=element_blank(), axis.title.x=element_blank())
#ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-children-go.pdf"), width=7.5, height=7)
ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-children-go-avgp.pdf"), width=7.5, height=6)

#
# (B) Venn diagram for children GO terms
#

# TSS
ntop <- nrow(S6_TSS_tb)
genealogy_tss <- .get_genealogy(category=S6_TSS_tb$category, ntop=ntop) 
sub_tss <- tss_genealogy%>% dplyr::filter(!children_in_top)

# 5UTR
ntop <- nrow(S6_5UTR_tb)
genealogy_5utr <- .get_genealogy(category=S6_5UTR_tb$category, ntop=ntop) 
sub_5utr <- genealogy_5utr %>% dplyr::filter(!children_in_top)

# 1st coding exon
ntop <- nrow(S6_1stE_tb)
genealogy_1stExon <- .get_genealogy(category=S6_1stE_tb$category, ntop=ntop)      
sub_1stExon <- genealogy_1stExon %>% dplyr::filter(!children_in_top)

childgo_tss <- S6_TSS_tb %>% right_join(sub_tss, by="category") 
childgo_5utr <- S6_5UTR_tb %>% right_join(sub_5utr, by="category")
childgo_1stExon <- S6_1stE_tb %>% right_join(sub_1stExon, by="category")

library(ggvenn)
# children terms
x <- list("5' UTR" = childgo_5utr$category, 
          "TSS"=childgo_tss$category,
          "1st codind exon" = childgo_1stExon$category)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))
ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-venndiagramm-go-children.pdf"), width=5, height=5)          


# all terms
x <- list("5' UTR" = S6_5UTR_tb$category, "TSS"=S6_TSS_tb$category, "1st codind exon" = S6_1stE_tb$category)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))
ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-venndiagramm-go-all.pdf"), width=5, height=5)


#
# (C) visualization for down-regualted TE in TSS, 5' UTR, and 1st coding exons
#

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
TE_TSS <- read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                      sheet="all_TSS") # status gives "down/up/none"
down_TSS <- TE_TSS %>% dplyr::filter(status == "down")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
TE_5UTR <- read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                        sheet="all_5UTR")
down_5UTR <- TE_5UTR %>% dplyr::filter(status == "down")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
TE_1stEx <- read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                        sheet="all_1st_exon")
down_1stEx <- TE_1stEx %>% dplyr::filter(status == "down")                        

# 1. venn diagrame
library(ggvenn)        
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "translation-changes")

x <- list("5' UTR" = down_5UTR$tx_name, "TSS"=down_TSS$tx_name, "1st coding exon"=down_1stEx$tx_name)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))
ggsave(file.path(fig_dir, "venndiagram-S6-5UTR-TSS-1stE-downregulated.pdf"), width=5, height=5)

# 2. heatmap and parallel plots
.parallel <- function(selected_tx_name) {
  df_tss <- TE_TSS %>% dplyr::filter(tx_name %in% selected_tx_name) %>% add_column(feature="TSS")
  df_5utr <- TE_5UTR %>% dplyr::filter(tx_name %in% selected_tx_name) %>% add_column(feature="5' UTR")
  df_1stEx <- TE_1stEx %>% dplyr::filter(tx_name %in% selected_tx_name) %>% add_column(feature="1st coding Exon")
  df <- df_tss %>% add_row(df_5utr) %>% add_row(df_1stEx) %>%
    dplyr::mutate(feature = factor(feature, levels=c("5' UTR", "TSS", "1st coding Exon")))
  ggplot(df, aes(x=feature, y=log2FoldChange, group=tx_name)) +
    geom_line(alpha=0.1, color="gray50") +
    theme_bw()
}

gg <- .parallel(selected_tx_name=down_TSS$tx_name)
ggsave(file.path(fig_dir, "parrellal-S6-5UTR-TSS-1stE-on-TSS-downregulated.pdf"), width=6, height=3.5)

gg <- .parallel(selected_tx_name=down_5UTR$tx_name)
ggsave(file.path(fig_dir, "parrellal-S6-5UTR-TSS-1stE-on-5UTR-downregulated.pdf"), width=6, height=3.5)

gg <- .parallel(selected_tx_name=down_1stEx$tx_name)
ggsave(file.path(fig_dir, "parrellal-S6-5UTR-TSS-1stE-on-1stExon-downregulated.pdf"), width=6, height=3.5)

###################### ignore the rest ########################################
#
# TSS parent go
#
ntop <- 50
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
S6_TSS_tb <- read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
tss_genealogy <- .get_genealogy(category=S6_TSS_tb$category, ntop=ntop) 

super_TSS <- S6_TSS_tb %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)
gg <- .parent_GO_barplot(super_TSS[1:10,], title="Around translation start sites [-13, 13]") 
plot(gg)
ggsave(file.path(fig_dir, "S6-TSS-translation-changes-parent-go.pdf"), width=6, height=3)

#
# 5'UTR parent go
#
ntop <- 25
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
S6_5UTR_tb <- read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
utr_genealogy <- .get_genealogy(category=S6_5UTR_tb$category, ntop=ntop) 

super_5UTR <- S6_5UTR_tb[1:30, ] %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)

gg <- .parent_GO_barplot(super_5UTR[1:10, ], title="5' UTR") 
plot(gg)
ggsave(file.path(fig_dir, "S6-5UTR-translation-changes-parent-go.pdf"), width=6, height=3)

#
# 1'st coding exons parent go
#
ntop <- 25
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
S6_1stE_tb <- read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet="Enriched_GO_for_down")
genealogy <- .get_genealogy(category=S6_1stE_tb$category, ntop=ntop)               
super_1stE <- S6_1stE_tb[1:30, ] %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)

gg <- .parent_GO_barplot(super_1stE[1:10, ], title="First coding exons") 
plot(gg)
ggsave(file.path(fig_dir, "S6-1stCodingExons-translation-changes-parent-go.pdf"), width=6, height=3)

#
# CDS parent go
#
ntop <- 25
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "CDS")
S6_CDS_tb <- read_xlsx(path=file.path(stats_dir, "S6_DOX_pulse_IFNg_vs_IFNg.xlsx"),
                sheet="enriched_GO_term_down")
genealogy <- .get_genealogy(category=S6_CDS_tb$category, ntop=ntop)               
super_CDS <- S6_CDS_tb[1:ntop, ] %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)

gg <- .parent_GO_barplot(super_CDS[1:10, ], title="CDS") 
plot(gg)
ggsave(file.path(fig_dir, "S6-CDS-translation-changes-parent-go.pdf"), width=8, height=3)

#
# top 25 "super-parents" among 5' UTR, TSS, and 1st coding exons
#
ntop = 25
inner <- S6_TSS_tb %>% inner_join(S6_5UTR_tb, by="category") %>% inner_join(S6_1stE_tb, by="category")
genealogy <- .get_genealogy(category=inner$category, ntop=ntop)               
super_inner <- inner[1:ntop, ] %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)

# Venn diagram
library(ggvenn)
x <- list("5' UTR" = S6_5UTR_tb$category, "TSS"=S6_TSS_tb$category, "1st codind exon" = S6_1stE_tb$category)
ggvenn(x, stroke_size = 0.5, set_name_size = 4, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))
ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-venndiagramm-go.pdf"), width=5, height=5)


#
# subset: get offspring union, size dots by adjusted p-value
#

#
# take the super-parent union and print out their p-value
#
union_super_parents <- super_5UTR %>% full_join(super_TSS, by="category") %>% full_join(super_1stE, by="category") %>%
  pull(category)

.sub_tb <- function(tb, union_category, feature) {
    tb %>% dplyr::filter(category %in% union_category) %>%
      dplyr::mutate(feature = feature) 
      #dplyr::select(-DEInCat)
}

# get all the GO term statistics
.get_all_go_term <- function(path, feature) {
   tb <- read_xlsx(path=path, sheet="translation_downregulation")
   selected <- tb %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
   go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=1,
                        return.DEInCat=FALSE, dds=dds_cds_by_gene) %>%
                dplyr::filter(category %in% union_super_parents) %>%
                dplyr::mutate(feature = feature) 
}

universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
S6_1stE_go <- .get_all_go_term(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"), "1st coding exon")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
S6_5UTR_go <- .get_all_go_term(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"), "5' UTR")

stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
S6_TSS_go <- .get_all_go_term(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"), "TSS")



df <- S6_1stE_go %>%
  add_row(S6_5UTR_go) %>% add_row(S6_TSS_go) %>%
  dplyr::mutate(log10padj = -log10(padj), 
                enriched=ifelse(padj < 0.05, "YES", "NO")) %>%
  dplyr::mutate(feature=factor(feature, levels=c("5' UTR", "TSS", "1st coding exon")))


ggplot(df, aes(y=term, x=feature)) +
  geom_point(aes(size=log10padj, color=enriched, fill=enriched)) +
  theme_bw() +
  labs(size=TeX("$-\\log_{10}(padj)$")) +
  scale_fill_manual(name ="enriched",  values=c("ivory4", "firebrick4")) +
  scale_color_manual(name="enriched",  values=c("ivory4", "firebrick4")) +
  theme(legend.position="bottom", legend.box = "vertical", 
        axis.title.y=element_blank(), axis.title.x=element_blank())
ggsave(file.path(fig_dir, "S6-5UTR-TSS-1stE-parent-go.pdf"), width=5.5, height=7)

# caption: take the union of the top super-parent GO terms among three features        

#
# intersection of transcripts down-regulated in 5' UTR, TSS and 1st exons and then GO analysis
#
sheet <- "translation_downregulation"
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
S6_TSS_tb <- read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet=sheet)
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
S6_5UTR_tb <- read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet=sheet)
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
S6_1stE_tb <- read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                sheet=sheet)     

library(ggVenn)
tmp1 <- S6_TSS_tb %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "") # 1063
tmp2 <- S6_5UTR_tb %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "") # 1745
tmp3 <- S6_1stE_tb %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "") # 592
# intersect is only 115 -> not much enrichment in GO terms

inner <- S6_TSS_tb %>% inner_join(S6_5UTR_tb, by="gene_id") %>% inner_join(S6_1stE_tb, by="gene_id")
universal <- rownames(dds_cds_by_gene) %>% str_replace("\\..*", "")
selected <- inner %>% distinct(gene_id) %>% pull(gene_id) %>% str_replace("\\..*", "")
enriched_go_down <- .do_goseq(universal=universal, selected_gene=selected,  p_value=1,
                              return.DEInCat=FALSE, dds=dds_cds_by_gene) # dds must be "by gene"

a <- .do_goseq(universal=universal, selected_gene=selected,  p_value=0.4,
                              return.DEInCat=TRUE, dds=dds_cds_by_gene)

tmp <- enriched_go_down %>% filter(over_represented_pvalue < 0.005)
genealogy <- .get_genealogy(category=tmp$category, ntop=86)               
super_inner <- tmp %>% left_join(genealogy, by="category") %>%
  dplyr::filter(!parent_in_top)
