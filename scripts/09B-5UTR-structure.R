# 09B-5UTR-structure.R
# This script calculate the RNA folding structure and free energy of 5' UTR (MEF) using Vienna.
# Note: adjutsted MEF is MEF per 100 nucleotides
# (1) Take a selected set of genes and compare their 5'UTR free energy against that of all 5' UTR
# The collections of genes include down-regulated translation changes in CDS, 1st exons, etc.s

library(Biostrings)
library(ShortRead)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyranges)
library(readxl)
library(BiocParallel)
library(ggplot2)
library(stringr)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)
bs_genome <- BSgenome.Hsapiens.UCSC.hg38
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)

#
# load feature 5' UTR and get sequences
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fold_dir <- file.path(pkg_dir, "RNAfold")

load(file.path(pkg_dir, "data", "tx_based_features.rda")) # feature_5p
feature_5p <- tx_based_features$feature_5p # note that those UTR are unique
setwd(fold_dir)

#
# tool: parse the output file
#
.parse_RNAfold <- function(path) {
    r <- readLines(path)
    name_idx <- seq(1,length(r), by=3)
    df <- data.frame(name = stringr::str_replace(r[name_idx], ">", ""),
                     width = nchar(r[name_idx + 1]),
                     MEF = sapply(stringr::str_split(r[name_idx + 2], " \\("), "[[", 2)) %>%
            dplyr::mutate(MEF = as.numeric(stringr::str_trim(stringr::str_replace(MEF, "\\)", ""))))
}

#
# all 5'UTR sequence
#

# task 1
library(tictoc)
tic("task 1")
unlist_feature <- unlist(feature_5p) %>% 
    plyranges::mutate(tx_name = names(.))
seq_5p <- getSeq(BSgenome.Hsapiens.UCSC.hg38, unlist_feature)
tmp <- splitAsList(seq_5p, names(seq_5p))
cat_exons = sapply(tmp, function(x) {
  comb <- c(unlist(x))
})
cat_exons <- as(cat_exons, "DNAStringSet")
toc()

#ShortRead::writeFasta(seq_5p, file=file.path(tmp_dir, "all-5UTR.fa"))
ShortRead::writeFasta(cat_exons, file=file.path(fold_dir, "all-5UTR.fa"))
system("rm all-5UTR.txt")
cmt <- sprintf("RNAfold -i %s --outfile=%s --noPS --job=12", 
               file.path(fold_dir, "all-5UTR.fa"), "all-5UTR.txt")
system(cmt)

df_MEF <- .parse_RNAfold(path=file.path(tmp_dir, "all-5UTR.txt")) %>%
  dplyr::mutate(adj_MEF = 100 * MEF / width)

df_MEF %>% ggplot(aes(x=width, y=MEF)) +
  geom_hex(bins=70) +
   theme_bw() +
   scale_fill_continuous(type = "viridis") 
ggsave(file.path(tmp_dir, "hex-width-vs-MEF.pdf"), width=5, height=4)   

df_MEF %>% 
  ggplot(aes(x=width, y=adj_MEF)) +
    geom_hex(bins=70) +
     theme_bw() +
     scale_fill_continuous(type = "viridis") 
ggsave(file.path(tmp_dir, "hex-width-vs-adjMEF.pdf"), width=5, height=4)

#
# RNA folding: 5'UTR + 1st coding exons; polysome (high / sub)
#

# (1) 5' UTR + 1st coding exons
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", 
                       "5UTR-plus-1st-exon")
S6_5UTR_1stexon_down <- read_xlsx(path=file.path(stats_dir, 
                                  "S6-5UTR-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                                  sheet="translation_downregulation")

# (2) polysome (high / total) - normalized by total
polysome <- read_xlsx(path=file.path(fold_dir, "Poly_lists_5'UTR MFE analysis.xlsx"), 
                      sheet = 1, skip=7)[, 1] %>%
  rename(`high/total`=`down...1`) %>%                      
  tidyr::drop_na(`high/total`)       
gene_id <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% polysome$`high/total`) %>% pull(gene_id)
df_polysome <- AnnotationDbi::select(txdb, keys=gene_id, columns="TXNAME", keytype="GENEID")

# (3) polysome - new list with different noramlization suggested by Andrew
polysome_new <- read_xlsx(path=file.path(fold_dir, "Poly-seq (new list).xlsx"), 
                      sheet = 1) 
gene_id <- as.data.frame(gene.anno) %>% 
  dplyr::filter(gene_name %in% polysome_new$`TE down Poly-seq (High/Sub)`) %>% 
  pull(gene_id)
df_polysome_new <- AnnotationDbi::select(txdb, keys=gene_id, columns="TXNAME", keytype="GENEID")


# (4) get MEF of all 5'UTR, and the rest of the group
df_MEF <- .parse_RNAfold(path=file.path(fold_dir, "all-5UTR.txt")) %>%
  dplyr::mutate(adj_MEF = 100 * MEF / width)

# get MEF of DUX4 targets 84
df_MEF_DUX4 <- .parse_RNAfold(path=file.path(fold_dir, "DUX4-targets-polysome-84.txt")) %>%
  dplyr::mutate(adj_MEF = 100 * MEF / width)

df_MEF_poly <- df_MEF %>% dplyr::filter(name %in% df_polysome$TXNAME) 
df_MEF_poly_new <- df_MEF %>% dplyr::filter(name %in% df_polysome_new$TXNAME) 

df_MEF_5p_1stExon <- df_MEF %>% dplyr::filter(name %in%  S6_5UTR_1stexon_down$tx_name)

# (5) combine all the categroies of gene lists
comb_MEF <- df_MEF %>% add_column(group = "all 5' UTR") %>%
  bind_rows(df_MEF_5p_1stExon %>% add_column(group = "Down in 5' UTR regions")) %>%
  bind_rows(df_MEF_DUX4 %>% add_column(group="DUX4 targets")) %>%
  bind_rows(df_MEF_poly %>% add_column(group = "Down polysome")) %>%
  dplyr::mutate(group=factor(group, 
                             level=c("all 5' UTR", "Down in 5' UTR regions", "Down polysome", "DUX4 targets")))

comb_MEF_new <- df_MEF %>% add_column(group = "all 5' UTR") %>%
  bind_rows(df_MEF_5p_1stExon %>% add_column(group = "Down in 5' UTR regions")) %>%
  bind_rows(df_MEF_DUX4 %>% add_column(group="DUX4 targets")) %>%
  bind_rows(df_MEF_poly_new %>% add_column(group = "Down polysome")) %>%
  dplyr::mutate(group=factor(group, 
                             level=c("all 5' UTR", "Down in 5' UTR regions", "Down polysome", "DUX4 targets")))

# (5) combine all the categroies of gene lists
segment_df <- data.frame(x=c(1, 1, 1), y=c(-70, -78, 0), xend=c(2, 3, 4), yend=c(-70, -78, 0),
                         group=c(1, 2, 3))
all_DUX4 <- format(wilcox.test(df_MEF$adj_MEF, df_MEF_DUX4$adj_MEF)$p.value, format="3", digit=2)
all_down5p <- format(wilcox.test(df_MEF$adj_MEF, df_MEF_5p_1stExon$adj_MEF)$p.value, format="e", digit=2)
all_polysome <- format(wilcox.test(df_MEF$adj_MEF, df_MEF_poly$adj_MEF)$p.value, format="e", digit=2)
all_polysome_new <- format(wilcox.test(df_MEF$adj_MEF, df_MEF_poly_new$adj_MEF)$p.value, format="e", digit=2)


comb_MEF %>% 
  ggplot(aes(x=group, y=adj_MEF)) +
    geom_boxplot(width=0.7, aes(fill=group), outlier.shape=NA, show.legend=FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          axis.title.x = element_blank()) + 
    labs(y=" Minimal Free Energy (MEF)") +
    scale_fill_manual(values=c("transparent", '#E69F00', '#999999', '#56B4E9')) 
    #geom_segment(data=segment_df, aes(x=x, y=y, xend=xend, yend=yend)) +
    #geom_text(aes(family="serif"), x=1.5, y=-70, hjust = 0.5, vjust=-0.7, label=paste0("*** p=", all_down5p), size=2.5)  +
    #geom_text(aes(family="serif"), x=2, y=-78, hjust = 0.5, vjust=-0.7, label="* p=6.3e-4", size=2.5) +
    #geom_text(aes(family="serif"), x=2.5, y=0, hjust = 0.5, vjust=-0.7, label=paste0("*** p=", all_DUX4), size=2.5) +
ggsave(file.path(fold_dir, "boxplot-AdjMEF-5UTR-2.pdf"), width=2, height=4)  


comb_MEF_new %>% 
  ggplot(aes(x=group, y=adj_MEF)) +
    geom_boxplot(width=0.7, aes(fill=group), outlier.shape=NA, show.legend=FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          axis.title.x = element_blank()) + 
    labs(y=" Minimal Free Energy (MEF)") +
    scale_fill_manual(values=c("transparent", '#E69F00', '#999999', '#56B4E9')) 
ggsave(file.path(fold_dir, "boxplot-AdjMEF-5UTR-new-polysome.pdf"), width=2, height=4)  

comb_MEF_new %>% group_by(group) %>%
 summarize(mu = mean(adj_MEF))

comb_MEF %>% group_by(group) %>%
 summarize(mu = mean(adj_MEF))
 
#
# list fro polysome-seq
#
df_MEF <- .parse_RNAfold(path=file.path(tmp_dir, "all-5UTR.txt")) %>%
  dplyr::mutate(adj_MEF = 100 * MEF / width) %>%
  tibble::add_column(group="all")


polysome_list <- read_xlsx(path=file.path(pkg_dir, "extdata", "Sub-up_poly-down_overlap gene list.xlsx")) 
sub <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% polysome_list$`Gene List`) 
df <- AnnotationDbi::select(txdb, keys=sub$gene_id, columns="TXNAME", keytype="GENEID")
poly_MEF <- df_MEF %>% dplyr::filter(name %in% df$TXNAME) %>%
  dplyr::mutate(group="polysome-selected") 

df_MEF %>% bind_rows(poly_MEF) %>%
    ggplot(aes(x=group, y=adj_MEF)) +
      geom_boxplot(width=0.5, outlier.shape=NA) +
      theme_bw() +
      labs(y="adjusted MEF (MEF/#nucleotide/100)")
ggsave(file.path(tmp_dir, "boxplot-AdjMEF-polysome-selected.pdf"), width=2.5, height=4)  

#
# another list of polysome-seq: use poly-seq(high/total)
#
polysome_list <- read_xlsx(path=file.path(fold_dir, "Poly_lists_5'UTR MFE analysis.xlsx"), 
                           sheet = 1, skip=7)[, -c(3, 6, 9)] 
names(polysome_list) <- paste(rep(c("high/total", "low/total", "sub/total", "high/low", "high/sub"), each=2), 
                                rep(c("down", "up"), 5), sep="-")
polysome_list <- polysome_list %>%
  tidyr::gather(key=group, value=gene_name) %>%
  tidyr::drop_na(gene_name) %>%
  dplyr::filter(!group %in% c("low/total-down", "low/total-up", "high/total-up", "high/low-up")) %>%
  dplyr::mutate(group=factor(group))
sub <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% polysome_list$gene_name) 
df <- AnnotationDbi::select(txdb, keys=sub$gene_id, columns="TXNAME", keytype="GENEID")
poly_MEF <- map_dfr(levels(polysome_list$group), function(x) {
  sub_poly_list <- polysome_list %>% dplyr::filter(group == x)
  sub_anno <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% sub_poly_list$gene_name ) 
  df <- AnnotationDbi::select(txdb, keys=sub_anno$gene_id, columns="TXNAME", keytype="GENEID")
  df_MEF %>% dplyr::filter(name %in%  df$TXNAME) %>%
    dplyr::mutate(group=x)
}) %>% dplyr::mutate(ratio = sapply(str_split(as.character(group), "-"), "[[", 1)) %>%
  dplyr::mutate(group = factor(group, levels=c("high/total-down", "sub/total-down", "sub/total-up",
                                               "high/sub-down", "high/sub-up", "high/low-down")))

df_MEF %>% bind_rows(poly_MEF) %>%
    ggplot(aes(x=group, y=adj_MEF, fill=ratio)) +
      geom_boxplot(width=0.5, outlier.shape=NA) +
      theme_bw() +
      labs(y="adjusted MEF (MEF/#nucleotide/100)") +
      theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), 
           axis.title.x=element_blank(), legend.position="bottom")
ggsave(file.path(tmp_dir, "boxplot-AdjMEF-all-poly-lists.pdf"), width=5, height=4)  


