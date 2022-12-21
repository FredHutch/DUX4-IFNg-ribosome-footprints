# This script uses Suja's MB135 HDUX4 RNA-seq to determine the promoter
# of DUX4 targets using the SGSeq package. Need to cite the reference:
# Goldstein LD, Cao Y, Pau G, Lawrence M, Wu TD, Seshagiri S, Gentleman R (2016) 
# Prediction and Quantification of Splice Events from RNA-Seq Data. PLoS ONE 11(5): 
# e0156132. doi:10.1371/journal.pone.0156132
#
# Here we use the DUX4 targets list is prepared by Danielle: 'What DUX4 target genes to use_.xlsx'

# 
library(Biostrings)
library(ShortRead)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyranges)
library(readxl)
library(BiocParallel)
library(ggplot2)
library(stringr)
library(SGSeq)
library(org.Hs.eg.db)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)
bs_genome <- BSgenome.Hsapiens.UCSC.hg38
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
data(gene.anno)

# 
# define parameters 
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
fig_dir <- file.path(pkg_dir, "figures", "DUX4-targets-splice-graphs")
hdux4_dir <- file.path("/fh/fast/tapscott_s/CompBio", "RNA-Seq", "hg38.MB135.HDUX4")
targets <- readxl::read_xlsx(file.path(pkg_dir, "extdata", "What DUX4 target genes to use_.xlsx"),
                             sheet=1) %>%
  dplyr::select(`High Polysome  (basemean>50)         (86 genes)`, 
                `44 genes present in all 3 lists:`)  %>%
  rename(high_polysome=`High Polysome  (basemean>50)         (86 genes)`,
         combo=`44 genes present in all 3 lists:`)   

# get the SYMBOL and GENEID: 
targets_polysome <-  targets %>% dplyr::distinct(high_polysome) %>% 
  tidyr::drop_na(high_polysome) %>%
  dplyr::mutate(gene_id = mapIds(org.Hs.eg.db, keys=high_polysome, keytype="SYMBOL", 
                column="ENTREZID", multiVals="first")) %>%
  rename(gene_name=high_polysome)                
targets_combo <- targets %>% dplyr::distinct(combo) %>% 
  tidyr::drop_na(combo) %>%
  dplyr::mutate(gene_id = mapIds(org.Hs.eg.db, keys=combo, keytype="SYMBOL", 
                column="ENTREZID", multiVals="first")) %>%
  rename(gene_name = combo)                


load(file.path(hdux4_dir, "data", "sgfc.rda")) #
load(file.path(hdux4_dir, "data", "sgvc.rda")) #
load(file.path(hdux4_dir, "data","sgvc.AS.rda")) # findAltFirstExon.LTR() from `RNA-Seq/R_scripts/findAltFirstExon.R`




#
# A. combo: find alternative TSS using Splicing variants (AFE or AS)
#    According to sgvc, there are three types of gene models and give the nomenclature M1, M2, S3
#      S1: according to sgvc, it has aternative start sites (AS/AEF)
#      S2: do not have alternative AS/AEF
#      S3: do not have alternative AS/AEF, but they do, by observation
#

# which sgvc.AS contains the genes in targets_combo?
AS_DUX4_flag = bplapply(rowRanges(sgvc.AS), function(x) {
  flag = any(unlist(geneName(x)) %in% targets_combo$gene_id )
}, BPPARAM=bp_param) 

# extract the SG gene ID from sgvc.AS: data.frame with SG gene ID, gene name
# this results 13 AS gene model with 19 genes
M1 <- map_dfr(which(unlist(AS_DUX4_flag)), function(idx) {
  rng <- rowRanges(sgvc.AS)[[idx]]
  gene_name <- targets_combo %>% dplyr::filter(gene_id %in% unique(unlist(geneName(rng)))) %>% 
    pull(gene_name) %>% paste0(., collapse="/")
  c(SG_geneID = geneID(rng)[1], gene_name=gene_name)
}) %>% dplyr::distinct(SG_geneID, .keep_all=TRUE)

# 1. splice events for DUX4-targets that have alternative start sites; print by SG gene ID
pdf(file.path(fig_dir, "SpliceEvents-AS-combo44.pdf")) # 13 model, 19 genes
for (i in 1:nrow(M1))
  plotFeatures(sgfc, geneID=M1$SG_geneID[i], 
               main=M1$gene_name[i], color_novel="red", cex=0.5)
dev.off()

# 2. splice events for DUX4-targets that do not have aternative start sites
M2 <- targets_combo %>%
  dplyr::filter(!gene_name %in% unlist(str_split(M1$gene_name, "/"))) 

M2$SG_geneID <- sapply(S2$gene_id, function(x){
  model <- SGSeq:::restrictFeatures(sgfc, geneName = x) 
  if (nrow(model) == 0) {
    return()
  }
  unique(geneID(rowRanges(model)))
})
 
# plot splicign graphs
pdf(file.path(fig_dir, "SpliceEvents-non-AS-combo44.pdf"))  
for (i in 1:nrow(M2)) {
  plotFeatures(sgfc, geneName=M2$gene_id[i], 
               main=M2$gene_name[i], color_novel="red", cexRow=0.6)
}
dev.off() # note: some don't have splicing events (e.g. KDM4E)

# 3. observe and see which one might have alternative promoter?
# MBD3L2, MBD3L4, TPRX1, TRIM49B, and ZNF705G might have alternative promoter
# # these have very strange altered model (e.g. MBD3L2/MBD3L4/MBD3L5 are all link together)
# - get their SG geneID and find their gene model
M3 <- targets_combo %>% dplyr::filter(gene_name %in% c("MBD3L2", "MBD3L4", "TPRX1", "TRIM49B", "ZNF705G"))
M3$SG_geneID<- sapply(M3$gene_id, function(x){
  gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneName = x) )
  unique(geneID(gr))
})
pdf(file.path(fig_dir, "SpliceEvents-non-AS-combo44-subset.pdf"))
for (i in 1: nrow(M3))
 plotFeatures(sgfc, geneID=M3$SG_geneID[i], main=M3$gene_name[i], color_novel="red", cexRow=0.6)
dev.off()

#
# B. combo - find the 5'UTR GRanges
#

#
# (1) M1: probably have to hand pick some ; need to re-order the GRanges
#
alt_5UTR <- list()

# 1: SG_geneID 71 PRAMEF1/PRAMEF12/PRAMEF2/PRAMEF7 (+): also use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[1]))
alt_5UTR[["71-PRAMEF1/PRAMEF12/PRAMEF2/PRAMEF7-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["71-PRAMEF1/PRAMEF12/PRAMEF2/PRAMEF7-1"]] <- gr[[3]][type(gr[[3]]) == "E"]

# 2: SG geneID 853: PRAMEF11 (-): also use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[2]))
alt_5UTR[["853-PRAMEF11"]] <- gr[[1]][type(gr[[1]]) == "E"]

# 3: SG geneID 854: HNRNPCL1/PRAMEF6 (-): also use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[3]))
alt_5UTR[["854-HNRNPCL1/PRAMEF6-1"]] <- gr[[1]][type(gr[[1]]) == "E"] # need to switch the order by start
alt_5UTR[["854-HNRNPCL1/PRAMEF6-1"]] <- gr[[2]][type(gr[[2]]) == "E"] # need to switch the order by start

# 4: 859-PRAMEF13/PRAMEF8: use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[4]))
alt_5UTR[["859-PRAMEF13/PRAMEF8-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-2"]] <- gr[[2]][type(gr[[2]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-3"]] <- gr[[3]][type(gr[[3]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-4"]] <- gr[[4]][type(gr[[4]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-5"]] <- gr[[5]][type(gr[[5]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-6"]] <- gr[[6]][type(gr[[6]]) == "E"]
alt_5UTR[["859-PRAMEF13/PRAMEF8-7"]] <- gr[[7]][type(gr[[7]]) == "E"]


# 5: 1788      TRIM43: use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[5]))
alt_5UTR[["1788-TRIM43"]] <- gr[[1]][type(gr[[1]]) == "E"] # all three have the same exon, but different junction


# 6:  2338      TRIM43B: use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[6]))
alt_5UTR[["2338-TRIM43B-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["2338-TRIM43B-2"]] <- gr[[2]][type(gr[[2]]) == "E"]

# 7: 3532      SLC34A2: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[7]))
alt_5UTR[["3532-SLC34A2"]] <- gr[[1]][type(gr[[1]]) == "E"]

# 8:8354      TRIM49C/TRIM53AP: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[8]))
alt_5UTR[["8354-TRIM49C/TRIM53AP-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["8354-TRIM49C/TRIM53AP-2"]] <- gr[[2]][type(gr[[2]]) == "E"]
alt_5UTR[["8354-TRIM49C/TRIM53AP-3"]] <- gr[[3]][type(gr[[3]]) == "E"]
alt_5UTR[["8354-TRIM49C/TRIM53AP-4"]] <- gr[[4]][type(gr[[4]]) == "E"]
alt_5UTR[["8354-TRIM49C/TRIM53AP-5"]] <- gr[[5]][type(gr[[5]]) == "E"]

# 9: 8775      TRIM49/TRIM53AP: TRIM49 don't use annotated
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[9]))
alt_5UTR[["8775-TRIM49/TRIM53AP-1"]] <- gr[[1]][type(gr[[1]])=="E"]
alt_5UTR[["8775-TRIM49/TRIM53AP-3"]] <- gr[[3]][type(gr[[3]])=="E"]
alt_5UTR[["8775-TRIM49/TRIM53AP-5"]] <- gr[[5]][type(gr[[5]])=="E"]



# 10: 9714      CCNA1: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[10]))
alt_5UTR[["9714-CCNA1-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["9714-CCNA1-2"]] <- gr[[2]][type(gr[[2]]) == "E"]

# 11 : 12948     LEUTX: dont' use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[11]))
alt_5UTR[["12948-LEUTX"]] <- gr[[1]][type(gr[[1]]) == "E"]

# 12: 13613     DUXA: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[12]))
alt_5UTR[["13613-DUXA-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["13613-DUXA-2"]] <- gr[[2]][type(gr[[2]]) == "E"]


# 13: 14427     RFPL2: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=M1$SG_geneID[13]))
alt_5UTR[["14427-RFPL2-1"]] <- gr[[1]][type(gr[[1]]) == "E"]
alt_5UTR[["14427-RFPL2-2"]] <- gr[[2]][type(gr[[2]]) == "E"]

# M3: TPRX1(13527), TRIM49B (8169) and ZNF705G (6567)
# 1. M3: MBD3L2: this involved MDB3L2/3/4/5 ... so just use the annotated ones
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[1]))


# 3. TPRX1     284355         13527: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[3]))
alt_5UTR[["13527-TPRX1"]] <- gr[9]

# 4 TRIM49B: E1/E2/E3 # don't use annotation  TRIM49B   283116          8169
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[4]))
alt_5UTR[["8169-TRIM49B"]] <- gr[type(gr)=="E"][1:3]

# 5  ZNF705G   100131980       6567: dont' use annotatoin
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[5]))
alt_5UTR[["6567-ZNF705G"]] <- gr[c(37, 39, 41)]


#
# which targets use alternative and which use annotated 5' UTR?
# - M1 and M3[3:4] genes use alternative 5' UTR and the rest uses 
#

#
# get alternative 5' UTR sequence (targets_combo)
#
seq_alt <- lapply(alt_5UTR, function(gr) {
  # re-order by coordinates
  gr <- gr[order(start(gr))]
  gr <- as(gr, "GRanges")
  seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  c(unlist(seq))
})
seq_alt <- as(seq_alt, "DNAStringSet")

#
# get annotated 5' UTR sequence (targets_combo) and MEF
#
load(file.path(pkg_dir, "data", "tx_based_features.rda")) # feature_5p
feature_5p <- tx_based_features$feature_5p # note that those UTR are unique
genes <- targets_combo %>%
  dplyr::filter(!gene_name %in% unlist(str_split(M1$gene_name, "/"))) %>%
  dplyr::filter(!gene_name %in% M3$gene_name[3:5]) %>% pull(gene_name)
gene_id <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% genes) %>% pull(gene_id)
tx_name <- AnnotationDbi::select(hg38.HomoSapiens.Gencode.v35, 
                                 keys=gene_id, column="TXNAME", 
                                 keytype="GENEID") %>% pull(TXNAME)
sub <- feature_5p[names(feature_5p) %in% tx_name]
unlist_sub <- unlist(sub) %>% 
    plyranges::mutate(tx_name = names(.))
tmp <- getSeq(BSgenome.Hsapiens.UCSC.hg38, unlist_sub)
tmp <- splitAsList(tmp, names(tmp))
seq_anno = sapply(tmp, function(x) {
  comb <- c(unlist(x))
})
seq_anno <- as(seq_anno, "DNAStringSet")
seq_5p <- c(seq_alt, seq_anno)

fold_dir <- file.path(pkg_dir, "RNAfold")
ShortRead::writeFasta(seq_5p, file=file.path(fold_dir, "DUX4-targets-comb-44.fa"))
cmt <- sprintf("RNAfold -i %s --outfile=%s --noPS --job=4", 
               file.path(fold_dir, "DUX4-targets-comb-44.fa"), "DUX4-targets-comb-44.txt")
system(cmt)

#############################
## part II: targets_polysome
#############################

#
# M1 model
#

# 1) which genes have alternative transcript start site (AS/AFE)
AS_DUX4_flag = bplapply(rowRanges(sgvc.AS), function(x) {
  flag = any(unlist(geneName(x)) %in% targets_polysome$gene_id )
}, BPPARAM=bp_param) 

# 2) M1: extract the SG gene ID - edited gene model involving the gene_id of interest(ENTRENZID)

M1 <- map_dfr(which(unlist(AS_DUX4_flag)), function(idx) {
  rng <- rowRanges(sgvc.AS)[[idx]]
  gene_name <- targets_polysome %>% dplyr::filter(gene_id %in% unique(unlist(geneName(rng)))) %>% 
    pull(gene_name) %>% paste0(., collapse="/")
  c(SG_geneID = geneID(rng)[1], gene_name=gene_name)
}) %>% dplyr::distinct(SG_geneID, .keep_all=TRUE)

# 3) print out splice events for DUX4-targets that have alternative start sites; print by SG gene ID
pdf(file.path(fig_dir, "SpliceEvents-AS-polysome84.pdf")) # 18 model, 61 genes
for (i in 1:nrow(M1))
  plotFeatures(sgfc, geneID=M1$SG_geneID[i], 
               main=M1$gene_name[i], color_novel="red", cex=0.5)
dev.off()

# 4) get the GRagnes of AS/AEF
alt_5UTR <- bplapply(M1$SG_geneID, function(id) {
  grl <- rowRanges(SGSeq:::restrictFeatures(sgvc.AS, geneID=id))
  # extracting alternative 5' UTR
  AS <- lapply(grl, function(gr) {
    gr[type(gr) == "E"]
  })
}, BPPARAM=bp_param)
alt_5UTR <- unlist(alt_5UTR)

# 5) keep unique 5'UTR model: 49 (model_name will be used to generate the output file)
model_name <- map_dfr(alt_5UTR, function(x) {
  # model_name consists of geneID (SGSeq) + featureID
  c(SG_geneID=geneID(x)[1], gene_id = paste0("ENTREZID:", paste0(unique(unlist(geneName(a))), collapse="/")),
    model_name=paste0("g", geneID(x)[1], "-f", paste0(featureID(x), collapse="-")))
}) %>% tibble::rownames_to_column(var="idx_keep") %>%
  dplyr::mutate(idx_keep = as.numeric(idx_keep)) %>%
  dplyr::distinct(model_name, .keep_all=TRUE)
alt_5UTR <- alt_5UTR[model_name$idx_keep]  
names(alt_5UTR) <- model_name$model_name

#
# M2 and M3 model
# 

# 1) which genes are not in M1
M2 <- targets_polysome %>%
  dplyr::filter(!gene_name %in% unlist(str_split(M1$gene_name, "/"))) 

M2$SG_geneID <- sapply(S2$gene_id, function(x){
  model <- SGSeq:::restrictFeatures(sgfc, geneName = x) 
  if (nrow(model) == 0) {
    return()
  }
  unique(geneID(rowRanges(model)))
})
 
# 2) plot splicign graphs
pdf(file.path(fig_dir, "SpliceEvents-non-AS-polysome84.pdf"))  
for (i in 1:nrow(M2)) {
  plotFeatures(sgfc, geneName=M2$gene_id[i], 
               main=M2$gene_name[i], color_novel="red", cexRow=0.6)
}
dev.off() # note: some don't have splicing events (e.g. KDM4E)

# 3) M3: one that might have alternative 5'TUR; plot by SG gene model
candidates <- c("TPRX1", "TRIM49B", "ZNF705G", 
                "CLEC17A", "FAM90A27P", "ZNF705A",
                "ZNF679")
# NOTE: ZIM3 has a very complated model and intertwined with DUXA and I am not going to consider it.

M3 <- targets_polysome %>% dplyr::filter(gene_name %in% candidates)
M3$SG_geneID<- sapply(M3$gene_id, function(x){
  gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneName = x) )
  unique(geneID(gr))
})
pdf(file.path(fig_dir, "SpliceEvents-non-AS-polysome84-subset.pdf"))
for (i in 1: nrow(M3))
 plotFeatures(sgfc, geneID=M3$SG_geneID[i], main=M3$gene_name[i], color_novel="red", cexRow=0.6)
dev.off()

# 1: CLEC17A   388512        12796
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[1]))
alt_5UTR[["12796-CLEC17A"]] <- gr[type(gr)=="E"][1:2]

# 2. FAM90A27P 646508        13087
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[2]))
alt_5UTR[["13087-FAM90A27P"]] <- gr[type(gr)=="E"][1]

# 3. TPRX1     284355         13527: don't use annotation
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[3]))
alt_5UTR[["13527-TPRX1"]] <- gr[9]

# 4. TRIM49B: E1/E2/E3 # don't use annotation  TRIM49B   283116          8169
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[4]))
alt_5UTR[["8169-TRIM49B"]] <- gr[type(gr)=="E"][1:3]

# 5. ZNF679    168417         5635
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[5]))
alt_5UTR[["5635-ZNF679"]] <- gr[type(gr)=="E"][1]

# 6. ZNF705A   440077         8916; 5'UTR end at 8,172,626
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[6]))
alt_5UTR[["8916-ZNF705A"]] <- gr[type(gr)=="E"][1:4]
end(alt_5UTR[["8916-ZNF705A"]])[4] <- 8172626

# ZNF705G   100131980      6567
gr <- rowRanges(SGSeq:::restrictFeatures(sgfc, geneID=M3$SG_geneID[7]))
alt_5UTR[["6567-ZNF705G"]] <- gr[c(37, 39, 41)]

#
# get alternative 5' UTR sequence (targets_combo)
#
seq_alt <- lapply(alt_5UTR, function(gr) {
  # re-order by coordinates
  gr <- gr[order(start(gr))]
  gr <- as(gr, "GRanges")
  seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  c(unlist(seq))
})
seq_alt <- as(seq_alt, "DNAStringSet")

#
# Get sequence: alternative, M1 and M3, and non-alternative, M2
#
load(file.path(pkg_dir, "data", "tx_based_features.rda")) # feature_5p
feature_5p <- tx_based_features$feature_5p # note that those UTR are unique
genes <- targets_polysome %>%
  dplyr::filter(!gene_name %in% unlist(str_split(M1$gene_name, "/"))) %>%
  dplyr::filter(!gene_name %in% M3$gene_name) %>% pull(gene_name)

gene_id <- as.data.frame(gene.anno) %>% dplyr::filter(gene_name %in% genes) %>% pull(gene_id)
tx_name <- AnnotationDbi::select(hg38.HomoSapiens.Gencode.v35, 
                                 keys=gene_id, column="TXNAME", 
                                 keytype="GENEID") %>% pull(TXNAME)
sub <- feature_5p[names(feature_5p) %in% tx_name]
unlist_sub <- unlist(sub) %>% 
    plyranges::mutate(tx_name = names(.))
tmp <- getSeq(BSgenome.Hsapiens.UCSC.hg38, unlist_sub)
tmp <- splitAsList(tmp, names(tmp))
seq_anno = sapply(tmp, function(x) {
  comb <- c(unlist(x))
})
seq_anno <- as(seq_anno, "DNAStringSet")
seq_5p <- c(seq_alt, seq_anno)

fold_dir <- file.path(pkg_dir, "RNAfold")
ShortRead::writeFasta(seq_5p, file=file.path(fold_dir, "DUX4-targets-polysome-84.fa"))
cmt <- sprintf("RNAfold -i %s --outfile=%s --noPS --job=4", 
               file.path(fold_dir, "DUX4-targets-polysome-84.fa"), "DUX4-targets-polysome-84.txt")
system(cmt)


#
# out put the 5'UTR of DUX4 - DUX4-induced alternative and annotated
#


alt_ls <- map_dfr(names(alt_5UTR), function(x) {
  SG_ID <- str_replace(str_split(x, "-")[[1]][1], "g", "")
  if (SG_ID %in% M1$SG_geneID)
    assosciated_genes <- M1 %>% dplyr::filter(SG_geneID == SG_ID) %>% pull(gene_name)
  if (SG_ID %in% M3$SG_geneID)
    assosciated_genes <- M3 %>% dplyr::filter(SG_geneID == SG_ID) %>% pull(gene_name)  
  gr <- alt_5UTR[[x]]
  gr <- gr[order(start(gr))]  
  as.data.frame(alt_5UTR[[x]]) %>% 
    tibble::add_column(novel_5UTR=x, assosciated_genes = assosciated_genes, .before="seqnames")     
})

df <- AnnotationDbi::select(hg38.HomoSapiens.Gencode.v35, 
                            keys=names(sub), column="GENEID", 
                            keytype="TXNAME", multiVals="first") %>% 
  rename(gene_id="GENEID") %>%
  dplyr::left_join(dplyr::select(as.data.frame(gene.anno), gene_id, gene_name), by="gene_id")
                                 
anno_ls <- map_dfr(names(sub), function(x) {
  gr <- sub[[x]]
  as.data.frame(gr) %>% 
    tibble::add_column(tx_name=x, gene_id=dplyr::filter(df, TXNAME==x)$gene_id,
                        gene_name=dplyr::filter(df, TXNAME==x)$gene_name, .before="seqnames")
})

write_xlsx(list("novel-5UTR-exons" = alt_ls, "annotated-5UTR-exons" = anno_ls),
           path=file.path(pkg_dir, "manuscript/suppl-tables", "DUX4-targets-novel-and-annotated-5'UTR.xlsx"))