# applye DEXseq to both Ribo-seq and RNA-seq

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
# define parameters and load datasets
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
dxd_ribo <- get(load(file.path(pkg_dir, "data", "dxd_cds_parts.rda")))
dxd_rna <- get(load(file.path(pkg_dir, "data", "dxd_cds_parts_mRNA.rda")))
fig_dir <- file.path(pkg_dir, "figures", "dexseq")
source(file.path(pkg_dir, "scripts", "07A-tools_dexseq.R"))
load(file.path(pkg_dir, "data", "rse_cds_by_gene.rda"))

# list of comparisons
list_comp <- list(S1 = c("untreated", "DOX_pulse"),
                  S2 = c("untreated", "IFNg"),
                  S3 = c("untreated", "DOX_pulse_IFNg"),
                  S4 = c("IFNg", "DOX_pulse"),
                  S5 = c("DOX_pulse", "DOX_pulse_IFNg"),
                  S6 = c("IFNg", "DOX_pulse_IFNg" ))
#
# load DUX4_altered and INFg_altered
#
load(file.path(pkg_dir, "data", "DUX4_induced.rda"))
load(file.path(pkg_dir, "data", "IFNg_induced.rda"))

#
# define histone
#
gene_name <- gene.anno$gene_name
h2 <- grep("\\<H2A|\\<H2B", gene_name)
h1 <- grep("\\<H1-", gene_name)
histone_variants <- gene.anno[c(h1, h2), ]
dim(histone_variants) #96

#
# exclude histone and DUX4-induced
#
dxd_ribo <- dxd_ribo[!rowData(dxd_ribo)$groupID %in% c(histone_variants$gene_id, DUX4_induced$ensembl)]
dxd_rna <- dxd_rna[!rowData(dxd_rna)$groupID %in% c(histone_variants$gene_id, DUX4_induced$ensembl)]

# 
# tools
#
.filter_by_rowMeans <- function(rse, mean_filter = 15, treatments) {
     rse <- rse[, rse$treatment %in% treatments]
     rse <- rse[rowMeans(assays(rse)[["counts"]]) >= mean_filter]
}

#
# (1) S1 / Ribo-seq and RNA-seq
#

# (1a) filter
rse_ribo <- .filter_by_rowMeans(rse=rse_cds_by_gene, treatments = list_comp[["S1"]])
keep_groupID <- rownames(rse_ribo)

# (1B) Ribo-seq
dxd_S1_ribo <- dxd_ribo[rowData(dxd_ribo)$groupID %in% keep_groupID]
dxd_S1_ribo <- .do_dexseq(dxd_S1_ribo, treatments = list_comp[["S1"]])

# (1C) RNA-seq
dxd_S1_rna <- dxd_rna[rowData(dxd_rna)$groupID %in% keep_groupID]
dxd_S1_rna <- .do_dexseq(dxd_S1_rna, treatments=list_comp[["S1"]])

dxd_S1 <- list(ribo=dxd_S1_ribo, rna=dxd_S1_rna)
save(dxd_S1, file=file.path(pkg_dir, "data", "dxd_S1.rda"))

#
# (2) S2 / Ribo-seq and RNA-seq
#

# (2A) filter
rse_ribo <- .filter_by_rowMeans(rse=rse_cds_by_gene, treatments = list_comp[["S2"]])
keep_groupID <- rownames(rse_ribo)

# (2B) Ribo-seq
dxd_S2_ribo <- dxd_ribo[rowData(dxd_ribo)$groupID %in% keep_groupID]
dxd_S2_ribo <- .do_dexseq(dxd_S2_ribo, treatments = list_comp[["S2"]])

# (2C) RNA-seq
dxd_S2_rna <- dxd_rna[rowData(dxd_rna)$groupID %in% keep_groupID]
dxd_S2_rna <- .do_dexseq(dxd_S_rna, treatments=list_comp[["S2"]])

dxd_S2 <- list(ribo=dxd_S2_ribo, rna=dxd_S2_rna)
save(dxd_S2, file=file.path(pkg_dir, "data", "dxd_S2.rda"))


#
# (6) S6 / Ribo-seq and RNA-seq
#
rse_ribo <- .filter_by_rowMeans(rse=rse_cds_by_gene, treatments = list_comp[["S6"]])
keep_groupID <- rownames(rse_ribo)

dxd_S6_ribo <- dxd_ribo[rowData(dxd_ribo)$groupID %in% keep_groupID]
dxd_S6_ribo <- .do_dexseq(dxd_S6_ribo, treatments=list_comp[["S6"]])

dxd_S6_rna <- dxd_rna[rowData(dxd_ribo)$groupID %in% keep_groupID]
dxd_S6_rna <- .do_dexseq(dxd_S6_rna, treatments=list_comp[["S6"]])
dxd_S6 <- list(ribo=dxd_S6_ribo, rna=dxd_S6_rna)
save(dxd_S6, file=file.path(pkg_dir, "data", "dxd_S6.rda"))

# sanity check: PSMB9
assays(dxd_S6_rna["ENSG00000240065.8:E003"])[[1]]
idx <- dxd_ribo$condition %in% list_comp[["S6"]]
assays(dxd_ribo["ENSG00000240065.8:E003"])[[1]][, idx]
assays(dxd_S6_ribo["ENSG00000240065.8:E003"])[[1]]

# sanity check
as.data.frame(rowData(dxd_S6_rna)) %>% dplyr::filter(groupID=="ENSG00000240065.8")
as.data.frame(rowData(dxd_S6_ribo)) %>% dplyr::filter(groupID=="ENSG00000240065.8")
