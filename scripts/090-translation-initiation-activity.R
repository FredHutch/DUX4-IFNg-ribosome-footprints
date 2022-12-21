# 090-translation-initiation-activity.R 
# TIA: translation initiation activity
# This script processes Ribo-TISH results.
# previously I used Ribo-TISH to predict the ORFs in four conditions and 
# intent to use the `difftis` function but keeps getting errors. So I 
# decided to write a function to do differential analysis on initiation activity in R.
#
# 1. features: open reading frame on 5' UTR (predicted by RiboTish)
# 2. merge the features detected from two conditions
# 3. get p-sites counts for the features (output to xlsx sheets for Danielle)
# 4. (Maybe) find translation initiatin activity by DESeq2; size factor determined by p-site 
#    counts on CDS `dds_cds_by_gene.rda`.

# NOTE: Tis refers to Translation Initiation Site; TIA - translation initiation activity

#
# load libraries
#

library(DESeq2)
library(plyranges)
library(GenomicAlignments)
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
# define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
ribotish_dir <- file.path(pkg_dir, "ribotish")
load(file.path(pkg_dir, "data", "p_sites.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda")) # use for sizeFactor and colData
ignore.strand <- FALSE


#
# TIA flow begins
#

# 1. load the ORF results
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
  gr <- gr %>% plyranges::filter(!duplicated(GenomePos))
  gr                
})
names(orf_list) <- str_replace(basename(orf_pred), "_pred.txt", "")

# 2. make bed files
bplapply(names(orf_list), function(x) {
    rtracklayer::export(orf_list[[x]], 
                        con=file.path(dirname(orf_files)[1], paste0(x, "_pred.bed")), format="BED")
})

# 3. get p-site count and print out the excell sheets
lapply(names(orf_list), function(x) {
    message(x)
    se <- bplapply(p_sites, function(pc) {
            summarizeOverlaps(features=orf_list[[x]],
                              reads=pc,
                              inter.feature=FALSE,
                              ignore.strand=FALSE)
    })
    se <- do.call(cbind, se)
    col_data <- colData(dds_cds_by_gene)
    colData(se) <- col_data 
    dds <- DESeq2::DESeqDataSet(se, design=~treatment)
    # 3. output the results; normalization use dds_cds_by_gene sizeFactors
    cnt_norm <- counts(dds, normalized=TRUE) %>% as.data.frame()
    colnames(cnt_norm) <- paste0(colnames(cnt_norm), "_norm")
    output <- as.data.frame(rowData(dds)) %>%
      add_column(cnt_norm)
    write_xlsx(x=output, path=file.path(dirname(orf_files)[1], paste0(x, "_pred.xlsx")))  
})



#
# pull out TisType == 5' UTR
#

# define features
features <- c("5'UTR", "5'UTR:CDSFrameOverlap", "5'UTR:Known")

#
# IFNg vs. untreated adn 5'UTR only?
#
orf_IFNg <- orf_list[["IFNg"]] %>% plyranges::filter(TisType %in% features) %>%
  plyranges::mutate(source="IFNg", overlap_status = "None")
orf_untreated <- orf_list[["untreated"]] %>% plyranges::filter(TisType %in% features) %>%
  plyranges::mutate(source="untreated", overlap_status = "None")

ov <- findOverlaps(orf_IFNg, orf_untreated) # 320 / 604 are overlapped
orf_IFNg$overlap_status[unique(queryHits(ov))] <- "overlap_untreated"
orf_untreated$overlap_status[unique(subjectHits(ov))] <- "overlap_IFNg"

# 1. merge; distict ORF based on GenomePos
merge <- c(orf_IFNg, orf_untreated) %>%
  plyranges::filter(!duplicated(GenomePos))

# 2. get the p-site counts (untreated and IFNg only)
sample_idx <- c(1, 5, 6, 7:9)
uORF_IFNg_untreated <- bplapply(p_sites[sample_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_IFNg_untreated <- do.call(cbind, uORF_IFNg_untreated)
col_data <- colData(dds_cds_by_gene)[names(p_sites)[sample_idx], ]
colData(uORF_IFNg_untreated) <- col_data
uORF_IFNg_untreated <- DESeq2::DESeqDataSet(uORF_IFNg_untreated, design=~treatment)

# 3. output the results
cnt_norm <- counts(uORF_IFNg_untreated, normalized=TRUE) %>% as.data.frame()
colnames(cnt_norm) <- paste0(colnames(cnt_norm), "_norm")
S2_output <- as.data.frame(rowData(uORF_IFNg_untreated)) %>%
  add_column(cnt_norm)

#
# DOX_pulse vs. untreated adn 5'UTR only?
#
orf_DOX <- orf_list[["DOX-pulse"]] %>% plyranges::filter(TisType %in% features) %>%
  plyranges::mutate(source="DOX-pulse", overlap_status = "None") %>%
  add_column(overlap_status=rwoData(dds)$overlap_status)


ov <- findOverlaps(orf_DOX, orf_untreated) 
length(unique(queryHits(ov))) # 190 / 283 overlap 
orf_DOX$overlap_status[unique(queryHits(ov))] <- "overlap_untreated"
orf_untreated$overlap_status[unique(subjectHits(ov))] <- "overlap_DOX"

# 1. merge; distict ORF based on GenomePos
merge <- c(orf_DOX, orf_untreated) %>%
  plyranges::filter(!duplicated(GenomePos))

# 2. get the p-site counts (untreated and IFNg only)
sample_idx <- c(1, 5, 6, 10:12)
uORF_DOX_untreated <- bplapply(p_sites[sample_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_DOX_untreated <- do.call(cbind, uORF_DOX_untreated)
col_data <- colData(dds_cds_by_gene)[names(p_sites)[sample_idx], ]
colData(uORF_DOX_untreated) <- col_data
uORF_DOX_untreated <- DESeq2::DESeqDataSet(uORF_DOX_untreated, design=~treatment)
# 3. output the results
cnt_norm <- counts(uORF_DOX_untreated, normalized=TRUE) %>% as.data.frame()
colnames(cnt_norm) <- paste0(colnames(cnt_norm), "_norm")
S1_output <- as.data.frame(rowData(uORF_DOX_untreated)) %>%
  add_column(cnt_norm)

#
# DOX_pulse vs. IFNg adn 5'UTR only?
#
orf_DOX_IFNg <- orf_list[["DOX-pulse_IFNg"]] %>% plyranges::filter(TisType %in% features) %>%
  plyranges::mutate(source="DOX_IFNg", overlap_status = "None")

ov <- findOverlaps(orf_DOX_IFNg, orf_IFNg) 
length(unique(queryHits(ov))) # 170/ 232 overlap 
orf_DOX_IFNg$overlap_status[unique(queryHits(ov))] <- "overlap_IFNg"
orf_IFNg$overlap_status[unique(subjectHits(ov))] <- "overlap_DOX_IFNg"

# 1. merge; distict ORF based on GenomePos
merge <- c(orf_DOX_IFNg, orf_IFNg) %>%
  plyranges::filter(!duplicated(GenomePos))

# 2. get the p-site counts (untreated and IFNg only)
sample_idx <- c(2:4, 7:9 )
uORF_DOX_IFNg_IFNg <- bplapply(p_sites[sample_idx], function(pc) {
    summarizeOverlaps(features=merge,
                      reads=pc,
                      inter.feature=FALSE,
                      ignore.strand=FALSE)
})
uORF_DOX_IFNg_IFNg <- do.call(cbind, uORF_DOX_IFNg_IFNg)
col_data <- colData(dds_cds_by_gene)[names(p_sites)[sample_idx], ]
colData(uORF_DOX_IFNg_IFNg) <- col_data
uORF_DOX_IFNg_IFNg <- DESeq2::DESeqDataSet(uORF_DOX_IFNg_IFNg, design=~treatment)
# 3. output the results
cnt_norm <- counts(uORF_DOX_IFNg_IFNg, normalized=TRUE) %>% as.data.frame()
colnames(cnt_norm) <- paste0(colnames(cnt_norm), "_norm")
S6_output <- as.data.frame(rowData(uORF_DOX_IFNg_IFNg)) %>%
  add_column(cnt_norm)


#
# output merged uORFs of IFNg+untreated, DOX+untreated, and DOX-IFNg+IFNg
#
write_xlsx(x=list(`IFNg and untreated`=S2_output, `DOX and untreated`=S1_output,
                  `DOX-pulse_IFNg and IFNg`=S6_output),
           path=file.path(pkg_dir, "ribotish", "merged-5primeUTR-ORFs.xlsx"))

dds_merged_uORFs <- list(`S1`=uORF_DOX_untreated, `S2`=uORF_IFNg_untreated, `S6`=uORF_DOX_IFNg_IFNg)
save(dds_merged_uORFs, file=file.path(pkg_dir, "data", "dds_merged_uORFs.rda"))

# bed files
library(rtracklayer)
rtracklayer::export(orf_DOX_IFNg, con=file.path(pkg_dir, "ribotish", "uORF_DOX_IFNg.bed"), format="BED")
rtracklayer::export(orf_IFNg, con=file.path(pkg_dir, "ribotish", "uORF_IFNg.bed"), format="BED")
rtracklayer::export(orf_DOX, con=file.path(pkg_dir, "ribotish", "uORF_DOX.bed"), format="BED")
rtracklayer::export(orf_untreated, con=file.path(pkg_dir, "ribotish", "uORF_untreated.bed"), format="BED")

#
# Venn diagrames of uORFs between two sources
#
library(ggvenn)
features <- c("5'UTR", "5'UTR:CDSFrameOverlap", "5'UTR:Known")
orf_IFNg <- orf_list[["IFNg"]] %>% plyranges::filter(TisType %in% features) 
orf_untreated <- orf_list[["untreated"]] %>% plyranges::filter(TisType %in% features)
orf_DOX_IFNg <- orf_list[["DOX-pulse_IFNg"]] %>% plyranges::filter(TisType %in% features) 
orf_DOX <- orf_list[["DOX-pulse"]] %>% plyranges::filter(TisType %in% features) 

genome_pos <- list(`DOX_IFNg`=orf_DOX_IFNg$GenomePos,
                   `DOX`=orf_DOX$GenomePos,
                   `IFNg`=orf_IFNg$GenomePos,
                   `untreated`=orf_untreated$GenomePos)
fill_color <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")                   

ggvenn(genome_pos[c(2, 4)], fill_color = fill_color[c(2,4)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-uORFs-DOX-untreated.pdf"), width=5, height=5)

ggvenn(genome_pos[c(3, 4)], fill_color = fill_color[c(3,4)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-uORFs-IFNg-untreated.pdf"), width=5, height=5)

ggvenn(genome_pos[c(1, 3)], fill_color = fill_color[c(1,3)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-uORFs-DOX+IFNg-IFNg.pdf"), width=5, height=5)

ggvenn(genome_pos, fill_color = fill_color,
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-uORFs-four-treatements.pdf"))


#
# how about internal?
#
features <- c("Internal", "Internal:CDSFrameOverlap", "Internal:Known")
orf_IFNg <- orf_list[["IFNg"]] %>% plyranges::filter(TisType %in% features) 
orf_untreated <- orf_list[["untreated"]] %>% plyranges::filter(TisType %in% features)
orf_DOX_IFNg <- orf_list[["DOX-pulse_IFNg"]] %>% plyranges::filter(TisType %in% features) 
orf_DOX <- orf_list[["DOX-pulse"]] %>% plyranges::filter(TisType %in% features) 

genome_pos <- list(`DOX_IFNg`=orf_DOX_IFNg$GenomePos,
                   `DOX`=orf_DOX$GenomePos,
                   `IFNg`=orf_IFNg$GenomePos,
                   `untreated`=orf_untreated$GenomePos)
fill_color <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")                   

ggvenn(genome_pos[c(2, 4)], fill_color = fill_color[c(2,4)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-internals-DOX-untreated.pdf"), width=5, height=5)

ggvenn(genome_pos[c(3, 4)], fill_color = fill_color[c(3,4)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-internals-IFNg-untreated.pdf"), width=5, height=5)

ggvenn(genome_pos[c(1, 3)], fill_color = fill_color[c(1,3)],
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-internals-DOX+IFNg-IFNg.pdf"), width=5, height=5)

ggvenn(genome_pos, fill_color = fill_color,
       stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(pkg_dir, "ribotish", "venn-internals-four-treatements.pdf"))
