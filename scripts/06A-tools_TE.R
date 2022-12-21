# tools_TE.R
#
# tools:
#
.comb_RNA_Ribo <- function(rse_rna, rse_ribo, treatments, 
                           rna_mean_filter=80, ribo_mean_filter=15) {
  # input: rse from RNA-seq and Ribo-seq; output: dds
  # filter (12) and merge rse_cds_mRNA and rse_cds_by_gene and get subset of specific treatments 
  design = ~ protocol + treatment + protocol:treatment #
  col_data_column <- c("bam_files", "sample_name", "treatment", "protocol")

  # subsetting RNA-seq; get size factors
  rse_rna <- rse_rna[, rse_rna$treatment %in% treatments]
  rse_rna$treatment <- factor(rse_rna$treatment, levels=treatments)
  rse_rna <- rse_rna[rowMeans(assays(rse_rna)[["counts"]]) >= rna_mean_filter]
  rse_rna$treatment <- factor(rse_rna$treatment)
  rse_rna$protocol <- "RNA"
  colData(rse_rna) <- colData(rse_rna)[, col_data_column]
  dds_rna <- DESeqDataSet(rse_rna, design = ~ treatment)
  dds_rna <- estimateSizeFactors(dds_rna)

  # subsetting Ribo-seq; get size factors
  rse_ribo <- rse_ribo[, rse_ribo$treatment %in% treatments]
  rse_ribo$treatment <- factor(rse_ribo$treatment, levels=treatments)
  rse_ribo <- rse_ribo[rowMeans(assays(rse_ribo)[["counts"]]) >= ribo_mean_filter]
  rse_ribo$protocol <- "Ribo"
  colData(rse_ribo) <- colData(rse_ribo)[, col_data_column]
  dds_ribo <- DESeqDataSet(rse_ribo, design = ~ treatment)
  dds_ribo <- estimateSizeFactors(dds_ribo)

  # combine
  colnames(dds_rna) <- paste0(colnames(dds_rna), "_RNA")
  colnames(dds_ribo) <- paste0(colnames(dds_ribo), "_Ribo")
  keep <- intersect(rownames(dds_ribo), rownames(dds_rna))
  dds <- cbind(dds_rna[keep], dds_ribo[keep])
  dds$protocol <- factor(dds$protocol, levels=c("RNA", "Ribo"))
  design(dds) <- design
  return(dds)
}

.comb_RNA_Ribo_sizeFactor_by_cds <- function(rse_rna_by_tx, rse_ribo_by_tx, 
                                             rse_rna_cds_by_gene, rse_ribo_cds_by_gene, treatments, 
                                             rna_mean_filter=10, ribo_mean_filter=5) {
  # input: rse_rna_by_tx and rse_ribo_by_tx of small features by tx (TSS, 5'UTR, 3'UTR, and 1st exons) 
  # for both RNA-seq and Ribo-seq 

  design = ~ protocol + treatment + protocol:treatment #
  col_data_column <- c("bam_files", "sample_name", "treatment", "protocol", "sizeFactor")
  .subset_convert2dds <- function(rse, treatments) {
    rse <- rse[, rse$treatment %in% treatments]
    rse$treatment <- factor(rse$treatment, levels=treatments)
    DESeqDataSet(rse, design = ~ treatment)
  }

  # (1) subsetting RNA-seq and assign sizeFactor of  rse_rna_cds_by_gene to rse_rna_by_tx
  dds_rna_cds_by_gene <- .subset_convert2dds(rse_rna_cds_by_gene, treatments=treatments)
  dds_rna_cds_by_gene <- estimateSizeFactors(dds_rna_cds_by_gene)

  dds_rna_by_tx <- .subset_convert2dds(rse_rna_by_tx, treatments=treatments)
  sizeFactors(dds_rna_by_tx) <- sizeFactors(dds_rna_cds_by_gene) 
  dds_rna_by_tx <- dds_rna_by_tx[rowMeans(counts(dds_rna_by_tx, normalized=TRUE)) >= rna_mean_filter]
  dds_rna_by_tx$protocol <- "RNA"
  colData(dds_rna_by_tx) <- colData(dds_rna_by_tx)[, col_data_column]
  
  # (2) subsetting Ribo-seq and assinge sizeFactor by rse_ribo_cds_by_gene
  dds_ribo_cds_by_gene <- .subset_convert2dds(rse_ribo_cds_by_gene, treatments)
  dds_ribo_cds_by_gene <- estimateSizeFactors(dds_ribo_cds_by_gene)

  dds_ribo_by_tx <- .subset_convert2dds(rse_ribo_by_tx, treatments)
  sizeFactors(dds_ribo_by_tx) <- sizeFactors(dds_ribo_cds_by_gene)
  dds_ribo_by_tx <- dds_ribo_by_tx[rowMeans(counts(dds_ribo_by_tx, normalized=TRUE)) >= ribo_mean_filter]
  dds_ribo_by_tx$protocol <- "Ribo"
  colData(dds_ribo_by_tx) <- colData(dds_ribo_by_tx)[, col_data_column]


  # combine
  colnames(dds_rna_by_tx) <- paste0(colnames(dds_rna_by_tx), "_RNA")
  colnames(dds_ribo_by_tx) <- paste0(colnames(dds_ribo_by_tx), "_Ribo")
  keep <- intersect(rownames(dds_ribo_by_tx), rownames(dds_rna_by_tx))
  rowData(dds_ribo_by_tx) <- rowData(dds_ribo_by_tx)[, names(rowData(dds_ribo_by_tx))]
  dds <- cbind(dds_rna_by_tx[keep], dds_ribo_by_tx[keep])
  dds$protocol <- factor(dds$protocol, levels=c("RNA", "Ribo"))
  design(dds) <- design
  return(dds)
}

.tidy_res_ribo_over_rna <- function(dds, inter_res, ribo_res, rna_res,
                              alpha=0.05, lfc_threshold=0) {
  # sanity check: lfc_threshold and baseMean_threshold must be numerical
  ribo_res <- as.data.frame(ribo_res) %>% rownames_to_column(var="ensembl") %>%
    dplyr::rename(ribo_lfc = log2FoldChange, ribo_padj = padj) %>%
    dplyr::select(ensembl, ribo_lfc, ribo_padj)
  rna_res  <- as.data.frame(rna_res)  %>% rownames_to_column(var="ensembl") %>%
    dplyr::rename(rna_lfc = log2FoldChange, rna_padj = padj) %>%
    dplyr::select(ensembl, rna_lfc, rna_padj)

  cnt <- as.data.frame(counts(dds[rownames(inter_res)], normalized=TRUE)) %>%
    rownames_to_column(var="ensembl")
  
  cnt_avg <- cnt %>% gather(sample_name, cnt, -ensembl) %>% 
    dplyr::mutate(treatment=colData(dds)[sample_name, "treatment"]) %>%
    dplyr::mutate(protocol=colData(dds)[sample_name, "protocol"]) %>%
    group_by(ensembl, treatment, protocol) %>% summarise(avg=mean(cnt)) %>%
    dplyr::mutate(group=paste0(treatment, "_", protocol)) %>% 
    dplyr::mutate(group=factor(group, 
                               levels=c(paste0(levels(treatment), "_", levels(protocol)[1]), 
                                        paste0(levels(treatment), "_", levels(protocol)[2])))) %>%
    ungroup() %>% 
    dplyr::select(-treatment, -protocol) %>% 
    tidyr::spread(group, avg) %>%
    dplyr::rename_with(~paste0(.x, "_avg")) %>% rename(ensembl = "ensembl_avg") %>%
    dplyr::left_join(cnt, by="ensembl")

    
  df <- as.data.frame(inter_res) %>% 
    rownames_to_column(var="ensembl") %>%
    dplyr::mutate(gene_name=rowData(dds[ensembl])$gene_name,
                  gene_type=rowData(dds[ensembl])$gene_type, .before=2)  %>%  
    left_join(ribo_res, by="ensembl") %>%
    left_join(rna_res, by="ensembl") %>%                 
    dplyr::filter(!is.na(padj)) %>%               
    dplyr::mutate(status = "non_sig") %>%
    dplyr::mutate(status = ifelse(padj < alpha & log2FoldChange > lfc_threshold, "up", status)) %>%
    dplyr::mutate(status = ifelse(padj < alpha & log2FoldChange < -lfc_threshold, "down", status)) %>%
    dplyr::mutate(status = factor(status, levels=c("non_sig", "up", "down"))) %>%
    left_join(cnt_avg, by="ensembl") %>%
    arrange(padj)
}

.plot_volcano <- function(res, title=NULL, numerator="DOX_pulse") {
  #res <- as.data.frame(res)
  res %>% dplyr::mutate(log10_padj = -log10(padj)) %>%
   ggplot(aes(x=log2FoldChange, y=log10_padj, color=status)) +
     geom_point(size=0.8, alpha=0.5) +
     geom_vline(xintercept=1, linetype="dashed", color="gray50") +
     geom_vline(xintercept=-1, linetype="dashed", color="gray50") +
     geom_hline(yintercept = -log10(0.05), linetype="dashed", color="gray50") +
     theme_classic() +
     labs(title=title, x="LFC", y=bquote(~-log[10]~(padj))) +
     theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
     scale_color_manual(values=c("gray", my_color[numerator], my_color[numerator])) +
     scale_x_continuous(minor_breaks = c(-1, 1))
}

.plot_scatter_ribo_vs_rna <- function(res, title=NULL, label_x_pos = 4, label_y_pos = -2) {
    msg_up   <- sprintf("log2(Ribo / RNA) > 1: %4.0f", sum(res$status == "up"))
    msg_down <- sprintf("log2(Ribo / RNA) < -1: %4.0f", sum(res$status == "down"))  
    pearson <- cor(res$rna_lfc, res$ribo_lfc)             

    ggplot(res, aes(x=rna_lfc, y=ribo_lfc)) +
      geom_point(size=2.5, alpha=0.5, aes(color=status, shape=DUX4_induced)) +
      theme_bw() +
      geom_vline(xintercept=c(-1, 1), linetype="dashed", alpha=0.5, color="gray50") +
      geom_hline(yintercept=c(-1, 1), linetype="dashed", alpha=0.5, color="gray50") +
      annotate("text", x=label_x_pos, y=label_y_pos, label=paste0("Pearson = ", format(pearson, digit=2)),
               hjust = 0, vjust=1) +
      annotate("text", x=label_x_pos, y=label_y_pos-0.7, label=msg_up, color="red", hjust = 0, vjust=1) +
      annotate("text", x=label_x_pos, y=label_y_pos-1.4, label=msg_down, color="blue", hjust=0, vjust=1) +
      scale_color_manual(values=c("gray75", "red", "blue")) +
      theme(legend.position="none", panel.grid.minor = element_blank(), 
            plot.title = element_text(hjust = 0.5)) +
      scale_shape_manual(values=c(19, 1))            
}

.plot_counts <- function(gene_id, gene_name, dds) {
  df <- data.frame(cnt=log2(counts(dds[gene_id], normalized=TRUE)[1, ] +1),
                   treatment=colData(dds)$treatment,
                   protocol=colData(dds)$protocol) 
  df_means <- df %>% group_by(treatment, protocol) %>%
    summarise(mean=mean(cnt))

  gg <- ggplot(df, aes(x=treatment, y=cnt)) +
    geom_point(size=1) +
    #stat_summary(fun.y=mean, colour="red", geom="line")
    facet_wrap( ~ protocol) + 
    geom_line(data=df_means, mapping=aes(x=treatment, y=mean, group=1), colour="red") +
    theme_bw() +
    labs(title=gene_name, y="log2(nrom_counts + 1)", x="")                
}
