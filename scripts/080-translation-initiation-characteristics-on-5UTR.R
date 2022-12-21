# Ideas from "Profiling of Small Ribosomal Subunits Reveals Modes and Regulation of Translation Initiation"
#  On 5'UTR can we find 
# 1. frequence of the first nucleotide - A/C/G/T (U)? For those affected by DUX4, which nucleotide do they start with (frequency)?
# 2. Do those 5'UTR (affected by DUX4) have TOP motif (CTTTT or more T)?
# 3. Are those 5'UTR affected by DUX4 short < 30 nt?

#
# load library
#
library(DESeq2)
library(readxl)
library(writexl)
library(readxl)
library(tidyverse)
library(corrr)
library(plyranges)
library(wesanderson)
library(goseq)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
bs_genome <- BSgenome.Hsapiens.UCSC.hg38
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# load dataset / sheets
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "UTR5-translation-initiation-featrues")
load(file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_TSS_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_5p_3p_cds_by_tx.rda"))
load(file.path(pkg_dir, "data", "rse_5UTR_1stExon_by_tx.rda"))

grl_5UTR <- rowRanges(rse_5p_3p_cds_by_tx[["rse_5p"]])
top_motif <- "CTTTT"

#
# tools
#

.get_1st_exon <- function(.x){
  exon1_name <- ifelse(.x$strand[1] %in% c("+", "*"),  .x$exon_name[1], .x$exon_name[nrow(.x)])
  .x %>% dplyr::filter(exon_name %in% exon1_name)
}

.first_nucleotide_by_tx <- function(tx_name, grl) {
    # sanitized: (does not filter duplicated exons)
    #   subset grl by tx_name input -> 
    #     convert to data.frme ->
    #       get the first exon ->
    #         re-arrange to tx_name order
  
  gr <- unlist(grl)
  sanitized <- gr %>% plyranges::filter(names(.) %in% tx_name) %>%
    plyranges::mutate(tx_name=names(.), coordinates=as.character(.))

  df_gr <- as.data.frame(mcols(sanitized)) %>%
    tibble::add_column(start=start(sanitized), 
                       end=end(sanitized), strand=as.character(strand(sanitized)),
                       seqnames=as.character(seqnames(sanitized))) %>%
    group_by(tx_name) %>%
    group_modify(~.get_1st_exon(.)) %>% ungroup() %>% plyranges::as_granges()

  #names(df_gr) <- df_gr$tx_name
  #df_gr <- df_gr[tx_name] # rearrange to tx_name

  df_sn <- promoters(df_gr, upstream=0, downstream=1 )
  df_gr$first_nucleotide  <- Biostrings::getSeq(bs_genome, df_sn)
  return(df_gr)
}

.first_nucleotide_frequency <- function(tx_name, grl) {
    # sanitized: 
    #   subset grl by tx_name input -> 
    #     convert to data.frme ->
    #       get the first exon ->
    #         remove duplicated exons by exon_name
  gr <- unlist(grl)
  sanitized <- gr %>% plyranges::filter(names(.) %in% tx_name) %>%
    plyranges::mutate(tx_name=names(.), coordinates=as.character(.))
  
  df <- as.data.frame(mcols(sanitized)) %>%
    tibble::add_column(start=start(sanitized), 
                       end=end(sanitized), strand=as.character(strand(sanitized)),
                       seqnames=as.character(seqnames(sanitized))) %>%
    group_by(tx_name) %>%
    group_modify(~.get_1st_exon(.)) %>% ungroup() %>%
    dplyr::distinct(exon_name, .keep_all=TRUE) %>% plyranges::as_granges()

  df_sn <- promoters(df, upstream=0, downstream=1 )
  seq <- getSeq(bs_genome, df_sn)

  table(as.character(seq))
}

.grl_convert_to_df <- function(tx_name, grl) {
    gr <- unlist(grl)
    sanitized <- gr %>% plyranges::filter(names(.) %in% tx_name) %>%
      plyranges::mutate(tx_name=names(.), coordinates=as.character(.))
  df <- as.data.frame(mcols(sanitized)) %>%
    tibble::add_column(start=start(sanitized), 
                       end=end(sanitized), strand=as.character(strand(sanitized)),
                       exon_width=width(sanitized),
                       seqnames=as.character(seqnames(sanitized))) 
}

.count_top_motif <- function(tx_name, grl) {
  .first_nucleotide_by_tx(tx_name=tx_name, grl) %>%
    plyranges::filter(!duplicated(exon_name)) %>%
    plyranges::mutate(seq=getSeq(bs_genome, .)) %>%
    plyranges::mutate(count_top_motif=Biostrings::vcountPattern("CTTTT", seq)) %>%
    plyranges::mutate(exon_width = width(.))
}


#
# whole genome
#

# ATCG frequecy
whole_genome_freq <- .first_nucleotide_frequency(tx_name = names(grl_5UTR), grl=grl_5UTR)
whole_genome_freq /sum(whole_genome_freq)

# frequency of CTTTT in 5'UTR; 
whole_genome <- .first_nucleotide_by_tx(tx_name=names(grl_5UTR), grl=grl_5UTR) %>%
  plyranges::filter(!duplicated(exon_name)) %>%
  plyranges::mutate(seq=getSeq(bs_genome, .)) %>%
  plyranges::mutate(count_top_motif=Biostrings::vcountPattern("CTTTT", seq)) %>%
  plyranges::mutate(exon_width = width(.))

summary_DF <- whole_genome %>% group_by(tx_name) %>%    
  summarize(count=sum(count_top_motif), 
            total_length=sum(exon_width))

sum(summary_DF$count > 0) / nrow(summary_DF) 
quantile(summary_DF$total_length)
# lengh distribution of CTTTT motif (use histogram?)


#
# subject: 5'UTR translation changes
#
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR_v2")
utr5_down <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="translation_downregulation")
utr5_up <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                             sheet="translation_upregulation")
utr5_tec <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-by_tx-exclude-DUX4-indueced.xlsx"),
                             sheet="all")

# utr5_tec (individual transripts)
df_gr <- .first_nucleotide_by_tx(utr5_tec$tx_name, grl=grl_5UTR)
df <- as.data.frame(mcols(df_gr)) %>% select(tx_name, first_nucleotide)
utr5_tec <- utr5_tec %>% left_join(df, by="tx_name")
pdf(file.path(fig_dir, "5UTR-first-nt-vs-translation-changes-S6.pdf" ), width=4, height=3.5)
ggplot(utr5_tec, aes(log2FoldChange, color=first_nucleotide)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.1, 0.76), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="5' UTR", x="Translation changes (LFC, S6)", y="Cumulative fraction")
dev.off()  

# (1)  up/down-specific                          
freq_down <- .first_nucleotide_frequency(tx_name = utr5_down$tx_name, grl=grl_5UTR)
freq_down / sum(freq_down)

freq_up <- .first_nucleotide_frequency(tx_name = utr5_up$tx_name, grl=grl_5UTR)
freq_up / sum(freq_up)

# (2) TOP motif
sanitized_gr_up <- .count_top_motif(tx_name=utr5_up$tx_name, grl=grl_5UTR)
sanitized_gr_down <- .count_top_motif(tx_name=utr5_down$tx_name, grl=grl_5UTR)

summary_df_down <- sanitized_gr_down %>% group_by(tx_name) %>%    
  summarize(count=sum(count_top_motif), 
            total_length=sum(exon_width))
sum(summary_df_down$count > 0) / nrow(summary_df_down) 

summary_df_up <- sanitized_gr_up %>% group_by(tx_name) %>%    
  summarize(count=sum(count_top_motif), 
            total_length=sum(exon_width))
sum(summary_df_up$count > 0) / nrow(summary_df_up) 

# TOP cumulative fraction
sanitized_5utr <- .count_top_motif(tx_name=utr5_tec$tx_name, grl=grl_5UTR)
sanitized_5utr <- sanitized_5utr %>%
  plyranges::mutate(top_motif = if_else(count_top_motif > 0, "TOP motif", "Others"))
df <- as.data.frame(mcols(sanitized_5utr))  
tmp <- utr5_tec %>% inner_join(df, by="tx_name") %>%
  dplyr::mutate(top_motif = factor(top_motif, levels=c("TOP motif", "Others")))

pdf(file.path(fig_dir, "5UTR-TOPmotif-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(tmp, aes(log2FoldChange, color=top_motif)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.18, 0.85), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="5' UTR", x="Translation changes", y="Cumulative fraction")
dev.off() 

#
# subject: 1st exons
#
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "1st_exon_v2")
first_exon_tec <- readxl::read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                                    sheet="all")
first_exon_down <- readxl::read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="translation_downregulation")
first_exon_up <- readxl::read_xlsx(path=file.path(stats_dir, "S6-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                             sheet="translation_upregulation")
df_gr <- .first_nucleotide_by_tx(first_exon_tec$tx_name, grl=grl_5UTR)
df <- as.data.frame(mcols(df_gr)) %>% select(tx_name, first_nucleotide)
first_exon_tec <- first_exon_tec %>% 
  left_join(df, by="tx_name") %>%
  drop_na(first_nucleotide)


pdf(file.path(fig_dir, "first-exon-5UTR-first-nucleotide-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(first_exon_tec, aes(log2FoldChange, color=first_nucleotide)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.1, 0.76), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="1st exon", x="Translation changes (LFC; S6)", y="Cumulative fraction")
dev.off() 

# (1)  up/down-specific                          
freq_down <- .first_nucleotide_frequency(tx_name = first_exon_down$tx_name, grl=grl_5UTR)
freq_down / sum(freq_down)

freq_up <- .first_nucleotide_frequency(tx_name = first_exon_up$tx_name, grl=grl_5UTR)
freq_up / sum(freq_up)

# (2) TOP motif
sanitized_1st_down <- .count_top_motif(tx_name=first_exon_down$tx_name, grl=grl_5UTR)
summary_df_down <- sanitized_1st_down %>% group_by(tx_name) %>%    
  summarize(count=sum(count_top_motif), 
            total_length=sum(exon_width))

sum(summary_df_down$count > 0) / nrow(summary_df_down) 
quantile(summary_df_down$total_length)

sanitized_1st_exons <- .count_top_motif(tx_name=first_exon_tec$tx_name, grl=grl_5UTR)
sanitized_1st_exons <- sanitized_1st_exons %>%
  plyranges::mutate(top_motif = if_else(count_top_motif > 0, "TOP motif", "Others"))
df <- as.data.frame(mcols(sanitized_1st_exons))  
tmp <- first_exon_tec %>% inner_join(df, by="tx_name") %>%
  dplyr::mutate(top_motif = factor(top_motif, levels=c("TOP motif", "Others")))

pdf(file.path(fig_dir, "first-exon-5UTR-TOPmotif-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(tmp, aes(log2FoldChange, color=top_motif)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.18, 0.85), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="1st coding exons", x="Translation changes", y="Cumulative fraction")
dev.off() 

# up
sanitized_1st_up <- .count_top_motif(tx_name=first_exon_up$tx_name, grl=grl_5UTR)
summary_df_up <- sanitized_1st_up %>% group_by(tx_name) %>%    
  summarize(count=sum(count_top_motif), 
            total_length=sum(exon_width))

sum(summary_df_up$count > 0) / nrow(summary_df_up) 
quantile(summary_df_up$total_length)

#
# subject CDS (S6)
#

#
# Subject: TSS that loss peaks (S6: ribo-seq differential analysis only)
#
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "around_TSS_v2")
tss_down <- readxl::read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="translation_downregulation") 
tss_up <- readxl::read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="translation_upregulation")  
tss_tec <- readxl::read_xlsx(path=file.path(stats_dir, "S6-TSS-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="all")                                                                 

# cumulated frequency
df_gr <- .first_nucleotide_by_tx(tss_tec$tx_name, grl=grl_5UTR)

df <- as.data.frame(mcols(df_gr)) %>% dplyr::select(tx_name, first_nucleotide)
first_exon_tec <- first_exon_tec %>% 
  left_join(df, by="tx_name") %>%
  drop_na(first_nucleotide)


pdf(file.path(fig_dir, "TSS-first-nucleotide-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(first_exon_tec, aes(log2FoldChange, color=first_nucleotide)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.1, 0.76), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="Translation start sites [-13, 13]", x="Translation changes (LFC; S6)", y="Cumulative fraction")
dev.off() 

freq_down <- .first_nucleotide_frequency(tx_name = tss_down$tx_name, gr=grl)
freq_down / sum(freq_down) 
#          A          C          G          T
# 0.36051159 0.16946443 0.39008793 0.07993605
freq_up <- .first_nucleotide_frequency(tx_name = tss_up$tx_name, gr=grl)
freq_up / sum(freq_up)
#          A          C          G          T
# 0.31707317 0.17073171 0.41951220 0.09268293

# TOP cumulative fraction
sanitized_tss <- .count_top_motif(tx_name=tss_tec$tx_name, grl=grl_5UTR)
sanitized_tss <- sanitized_tss %>%
  plyranges::mutate(top_motif = if_else(count_top_motif > 0, "TOP motif", "Others"))
df <- as.data.frame(mcols(sanitized_tss))  
tmp <- tss_tec %>% inner_join(df, by="tx_name") %>%
  dplyr::mutate(top_motif = factor(top_motif, levels=c("TOP motif", "Others")))

pdf(file.path(fig_dir, "TSS-TOPmotif-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(tmp, aes(log2FoldChange, color=top_motif)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.18, 0.85), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="Translation start sites [-13, 13]", x="Translation changes", y="Cumulative fraction")
dev.off() 


#
# subject: 5'UTR + 1st exons
#
stats_dir <- file.path(pkg_dir, "stats", "translation_changes", "5UTR-plus-1st-exon")
first_exon_tec <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                                    sheet="all")
first_exon_down <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                              sheet="translation_downregulation")
first_exon_up <- readxl::read_xlsx(path=file.path(stats_dir, "S6-5UTR-1st-exon-by_tx-exclude-DUX4-indueced.xlsx"),
                             sheet="translation_upregulation")
df_gr <- .first_nucleotide_by_tx(first_exon_tec$tx_name, grl=grl_5UTR)

df <- as.data.frame(mcols(df_gr)) %>% dplyr::select(tx_name, first_nucleotide)
first_exon_tec <- first_exon_tec %>% 
  left_join(df, by="tx_name") %>%
  drop_na(first_nucleotide)


pdf(file.path(fig_dir, "5UTR-plus-first-exon-first-nucleotide-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(first_exon_tec, aes(log2FoldChange, color=first_nucleotide)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.1, 0.76), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="5'UTR+1st exon", x="Translation changes (LFC; S6)", y="Cumulative fraction")
dev.off() 

# first nucleotides freqquency
freq_down <- .first_nucleotide_frequency(tx_name = first_exon_down$tx_name, gr=grl_5UTR)
freq_down / sum(freq_down) 
#         A         C         G         T
# 0.3004559 0.1939494 0.4053046 0.1002901
freq_up <- .first_nucleotide_frequency(tx_name = first_exon_up$tx_name, gr=grl_5UTR)
freq_up / sum(freq_up)

# TOP cumulative fraction
sanitized_5utr_1exon <- .count_top_motif(tx_name=first_exon_tec$tx_name, grl=grl_5UTR)
sanitized_5utr_1exon <- sanitized_5utr_1exon %>%
  plyranges::mutate(top_motif = if_else(count_top_motif > 0, "TOP motif", "Others"))
df <- as.data.frame(mcols(sanitized_5utr_1exon))  
tmp <- tss_tec %>% inner_join(df, by="tx_name") %>%
  dplyr::mutate(top_motif = factor(top_motif, levels=c("TOP motif", "Others")))

pdf(file.path(fig_dir, "5UTR-plus-1st-exon-TOPmotif-vs-translation-changes.pdf"), width=4, height=3.5)
ggplot(tmp, aes(log2FoldChange, color=top_motif)) +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position=c(0.18, 0.85), legend.title = element_blank()) +
  coord_cartesian(xlim = c(-4, 2)) + 
  labs(title="5' UTR + 1st coding exon", x="Translation changes", y="Cumulative fraction")
dev.off() 