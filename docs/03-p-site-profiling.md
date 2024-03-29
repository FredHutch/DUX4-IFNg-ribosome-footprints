# P-sites profiling on genomic features {#profiling}

In this chapter, we used `ribosomeProfilingQC::getPsiteCoordinates()` to get p-sites coordinates. Due to its large size, the `p-sites` dataset is not included in the repo. However, we do include the p-sites profiling on gene-based CDS and transcript-based genomic features such as 5' UTR, 13 nts up/downstream from translation start site, first coding exon, and 3' UTR, saved as `DESeqDatSet` instances.

Code chunk below loads libraries and defines local tools:

```r
library(DESeq2)
library(GenomicAlignments)
library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
txdb <- hg38.HomoSapiens.Gencode.v35


library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

cds_by_gene <- cdsBy(txdb, by="gene")
ignore.strand <- FALSE

#
# tools
#
.get_col_row_data <- function(rse, txdb, dds_cds_by_gene) {
  # colnames and append sample_info
  colnames(rse) <- colnames(dds_cds_by_gene)
  colData(rse) <- colData(dds_cds_by_gene)
  # rowData: tx_name -> gene_id, gene_name, gene_type
  tx_name <- rownames(rse)
  df <- AnnotationDbi::select(txdb, keys=tx_name, columns="GENEID", 
                              keytype="TXNAME",
                              multiVals="first") %>% 
    as.data.frame() %>%
    dplyr::distinct(TXNAME, .keep_all=TRUE) %>%
    dplyr::rename(tx_name=TXNAME, gene_id=GENEID) %>%
    dplyr::left_join(as.data.frame(gene.anno), by="gene_id") %>%
    dplyr::select(tx_name, gene_id, gene_type, gene_name, hgnc_id)                       
  rownames(df) <- df$tx_name
  rowData(rse) <- df[rownames(rse), ]
  
  return(rse)
}
```

## Define genomic features
The code chunk presented below defines the annotated transcript-based genomic features, including the 5' UTR, 13 nucleotides up/downstream from translation start sites, first coding exons, and 3' UTR. To prevent duplication resulting from isoforms sharing the same 5' UTR, exons, and 3' UTR, we retained only unique features. As a result, we obtained a list of `GRanges` instances (named `tx_based_features`) represent the coordinates of these genomic features. We didn't include this object in our repository.


```r
#
# define features
#
feature_5p  <- fiveUTRsByTranscript(txdb, use.name=TRUE)
feature_3p  <- threeUTRsByTranscript(txdb, use.name=TRUE)
feature_cds <- cdsBy(txdb, by="tx", use.name=TRUE)

# ensure to just include unique UTRs
.unique_UTRs <- function(utrs) {
    # exclude UTRs that are not unique
    exons_names_by_tx <- bplapply(utrs, function(gr)
      paste(gr$exon_name, collapse=","))
    keep_tx <- as.data.frame(unlist(exons_names_by_tx)) %>%
      rownames_to_column(var="tx_name") %>%
      dplyr::rename(exons_names=`unlist(exons_names_by_tx)`) %>%
      dplyr::distinct(exons_names, .keep_all=TRUE)
    utrs[keep_tx$tx_name]                 
}

unique_feature_5p <- .unique_UTRs(utrs=feature_5p)
unique_feature_3p <- .unique_UTRs(utrs=feature_3p)

# first coding exon and TSS [-13, 13]
first_exon_per_tx <- bplapply(feature_cds, function(gr) {
  st <- as.character(strand(gr[1]))
  if (st %in% c("+", "*")) x <- gr[1]
  if (st == "-") x <- gr[length(gr)]
  return(x)
})
first_exon_per_tx <- GRangesList(first_exon_per_tx)

.distinct_ranges <- function(gr) { # keeping distinct 1st coding exons only
    mcols(gr)$rng <- as.character(gr)
    gr_mcols <- as.data.frame(mcols(gr)) %>% dplyr::distinct(rng)
    gr <- gr[names(gr) %in% rownames(gr_mcols)]
}


first_exon <- .distinct_ranges(unlist(first_exon_per_tx))
#around_TSS <- promoters(first_exon, upstream=2, downstream=3)
around_TSS <- promoters(first_exon, upstream=13, downstream=14) # 27 nucleotides, 13 up/downstream of TSS

tx_based_features <- list(feature_5p=unique_feature_5p, 
                          feature_3p = unique_feature_3p,
                          first_exon = first_exon, 
                          around_TSS = around_TSS)
save(tx_based_features, file=file.path(pkg_dir, "data", "tx_based_features.rda"))     
```

## Get p-sites coordinates
Here we focused on the footprints with dominant lengths 26 to 29 and obtained their p-site coordinates using the offset of 13 nucleotides, as previously defined. To accomplish this, we employed the `ribosomeProfilingQC::getPsiteCoordinates()` function. While we did not store the resulting p_sites object in the repository, we did include the p-site profiling in the genomic features of interest.


```r
ignore.strand <- FALSE
yieldSize <- 1e+07
best_p_site <- 13
reads_len <- c(26:29)
anchor <- "5end"

# (a) get p site coordinates (we load from )
reads <- bplapply(sample_info$bam_files, function(f) {
  bam_file <- BamFile(file = f)
  pc <- getPsiteCoordinates(bam_file, bestpsite = best_p_site,
                            anchor = anchor)
  pc.sub <- pc[pc$qwidth %in% reads_len]
})
names(reads) <- sample_info$sample_name
p_sites <- reads
save(p_sites, file=file.path(pkg_dir, "data", "p_sites.rda"))
```

### Profiling on CDS
Use `GenomicAlignments::summarizeOverlaps()` to count p-sites on CDS and save as `DESeqDataSet` with size factor estimated by _DESeq2_. This size factor will also be used to normalize the counts on other transcript-based genomic features.


```r
# profiling
rse_cds_by_gene <- bplapply(reads, function(pc) {
  summarizeOverlaps(features=cds_by_gene, 
                    reads=pc, 
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_cds_by_gene <- do.call(cbind, rse_cds_by_gene)

# tidy colData and rowData
colnames(rse_cds_by_gene) <- sample_info$sample_name
colData(rse_cds_by_gene) <- append(colData(rse_cds_by_gene), 
                                   as(sample_info, "DataFrame"))
rowData(rse_cds_by_gene) <- gene.anno[rownames(rse_cds_by_gene), 
                                      c("gene_id", "gene_type",
                                        "gene_name", "hgnc_id")]
# dds with loose filtering with row sum > 12
dds_cds_by_gene <- 
  DESeqDataSet(rse_cds_by_gene, design = ~treatment) # loose filtering             
dds_cds_by_gene <- 
  dds_cds_by_gene[rowSums(counts(dds_cds_by_gene)) > 12]
dds_cds_by_gene <- estimateSizeFactors(dds_cds_by_gene)
save(dds_cds_by_gene, file=file.path(pkg_dir, "data",
                                    "dds_cds_by_gene.rda"))
```

### Profiling on other genomic features 

Profiling near translation start site [-13, 13]:

```r
load(file.path(pkg_dir, "data", "p_sites.rda")) # same as "reads" made previously
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda")) # need the sizeFactor and column data
load(file.path(pkg_dir, "data", "tx_based_features.rda")) # 5'UTR/TSS/1stExon/3'UTR
around_TSS <- tx_based_features$around_TSS
feature_5p <- tx_based_features$feature_5p
feature_3p <- tx_based_features$feature_3p
first_exon <- tx_based_features$first_exon
feature_cds <- cdsBy(txdb, by="tx", use.name=TRUE)
ignore.strand <- FALSE
inter.feature <- FALSE

rse_TSS_by_tx <- bplapply(p_sites, function(pc) {
    summarizeOverlaps(features=around_TSS, 
                      reads=pc, 
                      inter.feature=FALSE,
                      ignore.strand=ignore.strand)
})
rse_TSS_by_tx <- do.call(cbind, rse_TSS_by_tx)
rse_TSS_by_tx <- .get_col_row_data(rse_TSS_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)

save(rse_TSS_by_tx, file=file.path(pkg_dir, "data", "rse_TSS_by_tx.rda")) 
```


```r
# profile 5' UTR
rse_5UTR_by_tx <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=feature_5p, 
                    reads=pc, 
                    singleEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_5UTR_by_tx <- do.call(cbind, rse_5UTR_by_tx)
rse_5UTR_by_tx <- .get_col_row_data(rse_5UTR_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)
save(rse_5UTR_by_tx, file=file.path(pkg_dir, "data", "rse_5UTR_by_tx.rda"))
```


```r
# profile on 1st exon
rse_1st_exon_by_tx <- bplapply(p_sites, function(pc) {
    summarizeOverlaps(features=first_exon, 
                    reads=pc, 
                    inter.feature=FALSE,
                    ignore.strand=FALSE)
})
rse_1st_exon_by_tx <- do.call(cbind, rse_1st_exon_by_tx)
rse_1st_exon_by_tx <- .get_col_row_data(rse=rse_1st_exon_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)      
save(rse_1st_exon_by_tx, file=file.path(pkg_dir, "data", "rse_1st_exon_by_tx.rda")) 
```


```r
# profile by 3UTR
rse_3UTR_by_tx <- bplapply(p_sites, function(pc) {
  summarizeOverlaps(features=feature_3p, 
                    reads=pc, 
                    singleEnd=TRUE,
                    inter.feature=FALSE,
                    ignore.strand=ignore.strand)
})
rse_3UTR_by_tx <- do.call(cbind, rse_3UTR_by_tx)
rse_3UTR_by_tx <- .get_col_row_data(rse_3UTR_by_tx, dds_cds_by_gene=dds_cds_by_gene, txdb=txdb)
save(rse_3UTR_by_tx, file=file.path(pkg_dir, "data", "rse_3UTR_by_tx.rda"))
```



