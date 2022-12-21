
.make_pca_by_treatment <- function(dds, sample_labels=FALSE) {
    require(ggplot2)
    require(ggrepel)
    rlg <- rlog(dds, blind=TRUE)
    data <- plotPCA(rlg,
                    intgroup=c("treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))                
    gg <- ggplot(data, aes(PC1, PC2, color=treatment)) +
      geom_point(size=1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      labs(title="PC1 and PC2 ") +
      theme_minimal() 

    if(sample_labels) {
      gg <- gg + 
        geom_text_repel(aes(label=rownames(data)), size=3, show.legend=FALSE)
    }
    
    return(gg)
}

.tidy_combo <- function(sig, dds_5p, dds_3p, dds_cds, samples) {
  # combine counts on four regions: 5p, TSS, CDS, and 3p

  sample_info <- as.data.frame(colData(dds_5p))
  mat_5p <- counts(dds_5p[rownames(dds_5p) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>% 
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="5' UTR")
  mat_3p <- counts(dds_3p[rownames(dds_3p) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>% 
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="3' UTR")
  mat_cds <- counts(dds_cds[rownames(dds_cds) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>% 
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="CDS") 

  comb_df <- add_row(mat_5p, mat_3p) %>%
    add_row(mat_cds) %>%
    mutate(region=factor(region, levels=c("5' UTR", "CDS", "3' UTR")))  %>%
    dplyr::filter(sample_name %in% samples) %>%
    mutate(log2_cnt = log2(cnt+1)) %>%
    left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") 
}

.tidy_combo_2 <- function(sig, dds_5p, dds_3p, dds_cds, dds_tss, samples) {
  # combine counts on four regions: 5p, TSS, CDS, and 3p
  # at this point, dds_tss have not been assigne size factor; use the same one with dds_5p
  sizeFactors(dds_tss) <- sizeFactors(dds_5p)
  sample_info <- as.data.frame(colData(dds_5p))

  mat_5p <- counts(dds_5p[rownames(dds_5p) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="5' UTR")
  mat_3p <- counts(dds_3p[rownames(dds_3p) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="3' UTR")
  mat_cds <- counts(dds_cds[rownames(dds_cds) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="CDS")
  mat_tss <- counts(dds_tss[rownames(dds_tss) %in% sig$tx_name], normalized=TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var="tx_name") %>%
    gather(key=sample_name, value=cnt, -tx_name) %>%
    add_column(region="TSS")  

  comb_df <- add_row(mat_5p, mat_3p) %>%
    add_row(mat_cds) %>%
    add_row(mat_tss) %>%
    mutate(region=factor(region, levels=c("5' UTR", "TSS", "CDS", "3' UTR")))  %>%
    dplyr::filter(sample_name %in% samples) %>%
    mutate(log2_cnt = log2(cnt+1)) %>%
    left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name")
}

.cat2DEgenes <- function(GOID, pwf) {
    gene2cat <- getgo(rownames(pwf), "hg38", "ensGene", fetch.cats = "GO:BP")
    names(gene2cat) <- rownames(pwf)
    cat2gene <- goseq:::reversemapping(gene2cat)
    #' sanity check
    doesIDexist <- GOID %in% names(cat2gene)
    if (!all(doesIDexist)) stop("GOID is not found")
    sig_gene <- rownames(pwf)[pwf$DEgenes==1]
    geneInGOID <- cat2gene[[GOID]]
}

.do_goseq <- function(universal, selected_genes, p_value=0.01,
                      return.DEInCat=FALSE, dds=NULL) {
    require(goseq)
    library(org.Hs.eg.db)
    library(GO.db)
    library(hg38.HomoSapiens.Gencode.v35)

    txsByGene <- transcriptsBy(hg38.HomoSapiens.Gencode.v35, by="gene")
    names(txsByGene) <- str_replace(names(txsByGene), "\\..*", "")
    a <- duplicated( names(txsByGene))
    txsByGene <- txsByGene[!a]
    lengthData <- median(width(txsByGene))

    universal <- universal[!duplicated(universal)]
    isDEGs <- as.integer(universal %in% selected_genes)
    names(isDEGs) <- universal
    bias.data <- lengthData[names(isDEGs)]
    pwf <- nullp(isDEGs, bias.data=bias.data, plot.fit=FALSE)
    GO.BP <- goseq(pwf, "hg38", "ensGene", test.cats=c("GO:BP"))

    enriched.BP <- GO.BP %>%
      mutate(padj = p.adjust(over_represented_pvalue, method="BH")) %>%
      dplyr::filter(padj < p_value) %>%
      dplyr::select(-under_represented_pvalue)

    if (return.DEInCat & !is.null(dds)) {
        DEInCat <- purrr::map_chr(enriched.BP$category, function(GOID) {
          cat_genes <- .cat2DEgenes(GOID, pwf=pwf)
          cat_genes <- cat_genes[cat_genes %in% selected_genes]
          rownames(dds) <- str_replace(rownames(dds), "\\..*", "") # clean name
          cat_genename <- rowData(dds[cat_genes])$gene_name
          paste(cat_genename, collapse=",")
        }) %>% 
          as.data.frame() %>% 
          dplyr::rename(DEInCat=".") %>% 
          add_column(category=enriched.BP$category)
      enriched.BP <- enriched.BP %>% left_join(DEInCat, by="category")
    }
    return(enriched.BP)
}    

.get_gene_counts <- function(dds, targets) {
  # targets by gene names  
  col_data <- as.data.frame(colData(dds))
  targets_anno <- as.data.frame(rowData(dds)) %>%
    dplyr::filter(gene_name %in% targets)
  #targets_anno <- as(gene.anno, "data.frame") %>%
  #  dplyr::filter(gene_name %in% targets)

  targets_cnt <- as.data.frame(counts(dds[targets_anno$gene_id], normalized=TRUE)) %>%
    rownames_to_column(var="gene_id") %>%
    add_column(gene_name=targets_anno$gene_name) %>%
    dplyr::mutate(gene_name=factor(gene_name, levels=targets)) %>%
    gather(key=sample_name, value=count, -gene_name, -gene_id) %>%
    left_join(col_data %>% dplyr::select(sample_name, treatment), by="sample_name")
}


.metaPlot <- function (UTR5coverage, CDScoverage, UTR3coverage, sample, xaxis = c("RPFs",
    "mRNA"), bins = c(UTR5 = 100, CDS = 500, UTR3 = 100), ...)
{
    xaxis <- match.arg(xaxis)
    if (!is.list(UTR5coverage) || !is.list(UTR3coverage) || !is.list(CDScoverage)) {
        stop("UTR5coverage, CDScoverage and UTR3coverage must be\n         output of coverageDepth")
    }
    if (!xaxis %in% names(UTR5coverage) || !xaxis %in% names(UTR3coverage) ||
        !xaxis %in% names(CDScoverage)) {
        stop("UTR5coverage, CDScoverage and UTR3coverage must be\n         output of coverageDepth.",
            "And RPFs or mRNA must be available.")
    }
    if (!is.numeric(sample)) {
        sample <- which(names(CDScoverage[[xaxis]]$coverage) %in% sample)
    }
    if (length(sample) > 1) {
        sample <- sample[1]
        message("Only first sample will be plotted.")
    }
    if (!all(c("UTR5", "CDS", "UTR3") %in% names(bins))) {
        stop("bins should be a integer vector with names UTR5, CDS and UTR3.")
    }

    # coverage in CompressedRleList
    regions <- list(UTR5 = UTR5coverage[[xaxis]][["coverage"]][[sample]],
        CDS = CDScoverage[[xaxis]][["coverage"]][[sample]], UTR3 = UTR3coverage[[xaxis]][["coverage"]][[sample]])
    bins <- bins[names(regions)]
    regions <- mapply(regions, bins, FUN = function(.ele, .len) {
        .ele[lengths(.ele) >= .len]
    }, SIMPLIFY = FALSE)
    ids <- Reduce(intersect, lapply(regions, names))
    # convert to RleList
    regions <- lapply(regions, `[`, i = ids)
    metagene <- mapply(regions, bins, FUN = function(.ele, .len) {
        l <- lengths(.ele)
        ir <- IRanges(1, width = l)
        tile <- tile(ir, n = .len)
        names(tile) <- names(.ele)
        vws <- Views(.ele, tile)
        vms <- viewMeans(vws)
        vms <- do.call(rbind, as.list(vms))
        vms <- colMeans(vms)
    }, SIMPLIFY = FALSE)
    v <- unlist(metagene, use.names = TRUE)
}

.histone_variants <- function(gene.anoo) {
  gene.anno <- as.data.frame(gene.anno)
  gene_name <- gene.anno$gene_name
  h2 <- grep("\\<H2A|\\<H2B", gene_name)
  h1 <- grep("\\<H1-", gene_name)
  histone_variants <- gene.anno[c(h1, h2), ]
}