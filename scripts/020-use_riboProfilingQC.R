#####################
# ribosome footpring quality control: use the ribosomeProfilingQC package
#####################
# ml R/4.1.2-foss-2021b

library(ribosomeProfilingQC)
library(tidyverse)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
fig_dir <- file.path(pkg_dir, "figures", "QC")
source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "fork_readsEndPlot.R"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

#
# sample information and CDS
#
bam_dir <- file.path(scratch_dir, "bam", "merged_bam_runs")
bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)
sample_info <- data.frame(
  bam_files = bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)) %>%
    dplyr::mutate(sample_name = str_replace(basename(bam_files), ".bam", ""),
                  treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
    dplyr::mutate(treatment = factor(treatment, levels=c("untreated", "DOX-pulse", "IFNg", "DOX-pulse_IFNg")))

CDS <- ribosomeProfilingQC::prepareCDS(txdb)

#
# tools
#
.plot_read_end_hist <- function(vec) {
  tidy_df <- as.data.frame(vec) %>% 
    rownames_to_column(var="dist") %>%
    dplyr::mutate(dist=factor(dist, levels=dist)) %>%
    dplyr::mutate(perc_vec = 100 * (vec / sum(vec)) )
  gg <- ggplot(tidy_df, aes(x=dist, y=vec)) +
      geom_bar(stat="identity", width=0.7) +
      theme_minimal() +
      labs(x="Distance from 5' end of reads to start codon", y="counts") +
      geom_vline(xintercept = which(tidy_df$dist == -1) + 0.5, 
                 linetype="dashed", alpha=0.3, show.legend=FALSE) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, size=6))
}

.plot_read_end_line <- function(vec) {
  tidy_df <- as.data.frame(vec) %>% 
    rownames_to_column(var="dist") %>%
    dplyr::mutate(dist=factor(dist, levels=dist)) %>%
    dplyr::mutate(perc_vec = 100 * (vec / sum(vec)) )
  gg <- ggplot(tidy_df, aes(x=dist, y=vec, group=1)) +
      geom_line(color="blue", show.legend=FALSE) +
      theme_minimal() +
      labs(x="Distance from 5' end of reads to start codon", y="counts") +
      geom_vline(xintercept = which(tidy_df$dist == -1) + 0.5, linetype="dashed", 
                 alpha=0.3, show.legend=FALSE) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, size=6))
}

.plot_distance2codon <- function(distance) {
  tmp <- distance %>% as.data.frame(stringsAsFactors=FALSE) %>%
    dplyr::rename(index="Var1", Frequency="Freq") %>%
    dplyr::mutate(index=as.numeric(index), Frequency=as.numeric(Frequency)) %>%
    dplyr::mutate(frame = as.factor(index %% 3)) 
  gg <- ggplot(tmp, aes(x=index, y=Frequency, fill=frame)) +
    geom_bar(stat="identity", width=0.7) +
    theme_minimal() +
    theme(legend.position=c(0.9, 0.75), legend.key.size = unit(0.4, 'cm')) +
    labs(x="P site relative to start codon") +
    scale_x_continuous(breaks=seq(0, 50, 3)) +
    scale_fill_brewer(palette="Dark2")
}

.tidy_dist_data <- function(dist_list) {
  dist <- map_dfr(names(dist_list), function(x) {
    as.data.frame(dist_list[[x]]) %>%
    dplyr::rename(counts = `dist_list[[x]]`) %>%
    tibble::rownames_to_column(var = "dist") %>%
    tibble::add_column(sample_name = x) %>%
    dplyr::mutate(dist = factor(dist, levels=dist)) 
  }) %>%
    dplyr::left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
    dplyr::arrange(treatment) %>%
    dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
}

#
# (1) esitmate the optimal read lengths
#

# (a) read length frequency
read_length_freq <- bplapply(sample_info$bam_files, function(x) {
  bam_file <- BamFile(x)
  p_site <- estimatePsite(bam_file, CDS, genome)
  pc <- getPsiteCoordinates(bam_file, bestpsite = p_site)
  read_length_freq <- summaryReadsLength(pc, widthRange = c(25:39), plot=FALSE)
})
names(read_length_freq) <- sample_info$sample_name

# (b) tidy data
length_freq <- map_dfr(names(read_length_freq), function(x) {
  as.data.frame(read_length_freq[[x]]) %>%
    dplyr::rename(length="Var1") %>%
    add_column(sample_name=x) %>%
    dplyr::mutate(order = as.numeric(length)) %>%
    dplyr::mutate(length = as.numeric(as.character(length)))
}) %>%
  left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))

# (c) plot
ggplot(length_freq, aes(x=length, y=Freq)) +
  geom_bar(stat="identity", width=0.7) +
  theme_bw() +
  facet_wrap( ~ sample_name, nrow=4) +
  labs(x="Read length", y="Frequency")
ggsave(file.path(fig_dir, "freqment_size_frequency.pdf"))  

tidy_freq <- map_dfr(levels(length_freq$treatment), function(x) {
  length_freq %>% dplyr::filter(treatment == x) %>%
    dplyr::select(-order) %>%
          spread(sample_name, Freq) %>%
          dplyr::mutate(avg_freq = rowMeans(.[, 3:5]))
})

ggplot(tidy_freq, aes(x=length, y=avg_freq)) +
  geom_bar(stat="identity", width=0.7) +
  theme_bw() +
  facet_wrap(~treatment, nrow=1) +
  labs(x="Fragment size", y="Frequency") +
  theme(panel.grid.minor.y=element_blank())
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "QC")
ggsave(file.path(fig_dir, "fragment_size_frequency_per_treatment.pdf"), width=6, height=2)  

# (1b) % of [26:29] or [27:29]?
length_freq %>% filter(length %in% c(26:29)) %>%
  group_by(sample_name) %>%
  summarize(sum=sum(Freq))

#
# 2. distance from 5'end to start codon [-29, 30]
#
read_length <- c(26:29)
CDS_neg <- CDS[strand(CDS) == "-"]
CDS_pos <- CDS[strand(CDS) == "+"]
# (a) distance to start codon [-29, 30]
start_codon_30 <- bplapply(sample_info$bam_files, function(x) {
  bam_file <- BamFile(x)
  fork_readsEndPlot(bam_file, CDS, toStartCodon=TRUE, readLen=read_length, window=c(-29, 30))
  #ribosomeProfilingQC::readsEndPlot(bam_file, CDS_pos, toStartCodon=TRUE, readLen=read_length,
  #             window= c(-29, 30)) 
})
names(start_codon_30) <- sample_info$sample_name

# tidy data and hist (bar)
dist <- .tidy_dist_data(start_codon_30)
ggplot(dist, aes(x=dist, y=counts)) +
  geom_bar(stat="identity", width=0.7) +
  theme_bw() +
  labs(x="Distance from 5' end of reads to start codon", y="counts") +
  facet_wrap( ~ sample_name, nrow=4, scale="free") +
  geom_vline(xintercept = which(dist$dist == 1), 
             linetype="dashed", alpha=0.3, show.legend=FALSE) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, size=4),
        panel.grid.major = element_blank(), #panel.grid.minor = element_blank(), 
        panel.background = element_blank())
ggsave(file.path(fig_dir, "distance_from_5end_to_start_codon_30-fork-readsEndPlot.pdf"))#, width=8, height=6)

# (b) distance to start codon [-149, 150]
start_codon_150 <- bplapply(sample_info$bam_files, function(x) {
  bam_file <- BamFile(x)
  fork_readsEndPlot(bam_file, CDS, toStartCodon=TRUE, readLen=read_length, window=c(-149, 150))

})
names(start_codon_150) <- sample_info$sample_name

# tidy data
dist <- .tidy_dist_data(start_codon_150)
vline <- which(unique(dist$dist) == 1) 
ggplot(dist, aes(x=dist, y=counts, group=sample_name)) +
  geom_line(color="blue", show.legend=FALSE, size=0.3) +
  facet_wrap( ~ sample_name, nrow=6, scale="free") +
  theme_minimal() +
  labs(x="Distance from 5' end of reads to start codon [-149, 150]") +
  geom_vline(xintercept = vline, linetype="dashed", 
             alpha=0.3, show.legend=FALSE) +
  theme(axis.text.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank())            
ggsave(file.path(fig_dir, "distance_from_5end_start_codon_150-fork-readsEndPlot.pdf"), width=8, height=6)

# scale the counts by size factor; set equal y-axis
start_codon_150_scaled <- lapply(names(start_codon_150), function(x) {
  start_codon_150[[x]] / sizeFactors(dds_cds_by_gene)[x]
})
names(start_codon_150_scaled) <- names(start_codon_150)

dist <- .tidy_dist_data(start_codon_150_scaled)
vline <- which(unique(dist$dist) == 1) 
ggplot(dist, aes(x=dist, y=counts, group=sample_name)) +
  geom_line(color="blue", show.legend=FALSE, size=0.3) +
  facet_wrap( ~ sample_name, nrow=6) +
  theme_minimal() +
  labs(y="scaled counts; limited to [0, 1e6]", 
       x="Distance from 5' end of reads to start codon [-149, 150]") +
  geom_vline(xintercept = vline, linetype="dashed", 
             alpha=0.3, show.legend=FALSE) +
  theme(axis.text.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank()) 
  # coord_cartesian(ylim = c(0, 1e6))          
ggsave(file.path(fig_dir, "distance_from_5end_start_codon_150_scaled-fork-readsEndPlot.pdf"), width=8, height=6)

# (d) distance to end codon [-29, 30]
end_codon_30 <- bplapply(sample_info$bam_file, function(x) {
  bam_file <- BamFile(x)
  fork_readsEndPlot(bam_file, CDS, toStartCodon=FALSE, readLen=read_length, window=c(-29, 30))
  #readsEndPlot(bam_file, CDS, toStartCodon=FALSE, readLen=read_length, window= c(-29, 30))    
})   
names(end_codon_30) <- sample_info$sample_name
# tidy data and hist plot
dist <- .tidy_dist_data(end_codon_30)
vline <- which(unique(dist$dist) == 1) 
gg <- ggplot(dist, aes(x=dist, y=counts)) +
      geom_bar(stat="identity", width=0.7) +
      theme_bw() +
      labs(x="Distance from 5' end of reads to end codon [-29, 30]", y="counts") +
      facet_wrap( ~ sample_name, nrow=4, scale="free") +
      geom_vline(xintercept = vline, 
                 linetype="dashed", alpha=0.3, show.legend=FALSE) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, size=4),
            panel.grid.major = element_blank(), #panel.grid.minor = element_blank(), 
            panel.background = element_blank())
pdf(file.path(fig_dir, "distance_from_5end_end_codon_30-fork-readsENdPlot.pdf"), width=8, height=6)
plot(gg)
dev.off()

#
# (4) reading frame
#
reading_frame <- bplapply(sample_info$bam_files, function(x) {
  bam_file <- BamFile(x)
  p_site <- estimatePsite(bam_file, CDS, genome)
  pc <- getPsiteCoordinates(bam_file, bestpsite = p_site)
  pc_sub <- pc[pc$qwidth %in% read_length]
  pc_sub <- assignReadingFrame(pc_sub, CDS)
  distance <- plotDistance2Codon(pc_sub)
})

names(reading_frame) <- sample_info$sample_name

# tidy reading_frame tool
.tidy_reading_frame <- function(rf_list) {
  rf <- map_dfr(names(rf_list), function(x) {
    rf_list[[x]] %>% as.data.frame(stringsAsFactors=FALSE) %>%
      dplyr::rename(index="Var1", Frequency="Freq") %>%
      dplyr::mutate(index=as.numeric(index), Frequency=as.numeric(Frequency)) %>%
      dplyr::mutate(frame = as.factor(index %% 3)) %>%
      add_column(sample_name = x)
  })
}  

# ggplot
tidy_rf <- .tidy_reading_frame(reading_frame) %>%
  left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))

ggplot(tidy_rf, aes(x=index, y=Frequency, fill=frame)) +
    geom_bar(stat="identity", width=0.7) +
    theme_minimal() +
    facet_wrap( ~ sample_name, nrow=4, scale="free") +
    theme(legend.position=c(0.25, 0.93), legend.key.size = unit(0.3, 'cm'), 
          axis.text.x=element_text(size=5), panel.grid.major = element_blank()) +
    labs(x="P site relative to start codon", y="counts") +
    scale_x_continuous(breaks=seq(0, 50, 3)) +
    scale_fill_brewer(palette="Dark2")
ggsave(file.path(fig_dir, "reading_frame_psite_to_start_codon.pdf"), width=8, height=6)   


#
# (5) metacoverage: 5'UTR -> CDS -> 3'UTR
#

cov_utr5 <- coverageDepth(RPFs = sample_info$bam_files, gtf = txdb, region="utr5", 
                          readLen = read_length, bestpsite=13)
cov_cds  <- coverageDepth(RPFs = sample_info$bam_files, gtf = txdb, region="cds", 
                          readLen = read_length, bestpsite=13)
cov_utr3 <- coverageDepth(RPFs = sample_info$bam_files, gtf = txdb, region="utr3", 
                          readLen = read_length, bestpsite=13)

cov_utr5_cds_utr3 <- list(cov_utr5=cov_utr5, cov_cds=cov_cds, cov_utr3=cov_utr3)
save(cov_utr5_cds_utr3, file=file.path(pkg_dir, "data", "cov_utr5_cds_utr3.rda"))

#
# plot
#
load(file.path(pkg_dir, "data", "cov_utr5_cds_utr3.rda"))
cov_utr5 <- cov_utr5_cds_utr3$cov_utr5
cov_cds <- cov_utr5_cds_utr3$cov_cds
cov_utr3 <- cov_utr5_cds_utr3$cov_utr3

bins <- c(UTR5 = 100, CDS = 1000, UTR3 = 100)
meta_cov <- bplapply(sample_info$bam_files, function(x) {
  x <- basename(x)
  .metaPlot(UTR5coverage=cov_utr5, CDScoverage=cov_cds, UTR3coverage=cov_utr3, 
            sample=x, xaxis="RPFs", bins=bins) %>%
    as_tibble() %>%
    rownames_to_column(var="index") %>%
    dplyr::mutate(index = as.numeric(index))  %>%
    add_column(sample_name = str_replace(x, ".bam", "")) %>%
    add_column(regions=rep(names(bins), bins))
})

# tidy metaPlot coverages
tidy_meta_cov <- do.call(rbind, meta_cov) %>%
  dplyr::mutate(regions = factor(regions, levels=names(bins))) %>%
  dplyr::left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
 
bins_sum <- cumsum(bins)
ggplot(tidy_meta_cov, aes(group=sample_name)) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=1, xmax=bins_sum[1], ymin=-Inf, ymax=Inf), 
            fill="#66C2A5", alpha=0.1) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=bins_sum[1]+1, xmax=bins_sum[2], ymin=-Inf, ymax=Inf), 
            fill="#FC8D62", alpha=0.1) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=bins_sum[2]+1, xmax=bins_sum[3], ymin=-Inf, ymax=Inf), 
            fill="#8DA0CB", alpha=0.1) +            
  geom_line(aes(x=index, y=value, group=sample_name), color="blue", size=0.3) +
  theme_minimal() +
  facet_wrap( .~ sample_name, ncol=2, scale="free_y") +
  labs(y="mean of coverage", x="") +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank()) +
  geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
  geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
  scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                     labels=c("5'UTR", "CDS", "3'UTR")) 
ggsave(file.path(fig_dir, "meta_coverage.pdf"), width=10, height=8)  

#
# scaled meta_coverage; change bins
#
meta_cov_scale <- bplapply(meta_cov, function(x) {
  x %>% dplyr::mutate(scaled_value = value / sizeFactors(dds_cds_by_gene)[sample_name[1]]) 
})

tidy_meta_cov <- do.call(rbind, meta_cov_scale) %>%
  dplyr::mutate(regions = factor(regions, levels=names(bins))) %>%
  dplyr::left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
bins_sum <- cumsum(bins)

ggplot(tidy_meta_cov, aes(group=sample_name)) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=1, xmax=bins_sum[1], ymin=-Inf, ymax=Inf), 
            fill="#66C2A5", alpha=0.1) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=bins_sum[1]+1, xmax=bins_sum[2], ymin=-Inf, ymax=Inf), 
            fill="#FC8D62", alpha=0.1) +
  geom_rect(data=data.frame(sample_name=sample_info$sample_name),
            aes(xmin=bins_sum[2]+1, xmax=bins_sum[3], ymin=-Inf, ymax=Inf), 
            fill="#8DA0CB", alpha=0.1) +            
  geom_line(aes(x=index, y=scaled_value, group=sample_name), color="blue", size=0.3) +
  theme_minimal() +
  facet_wrap( .~ sample_name, ncol=2) +
  labs(y="mean of coverage (scaled, y-axis limited to 1)", x="") +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank()) +
  geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
  geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
  scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                     labels=c("5'UTR", "CDS", "3'UTR")) +
  coord_cartesian(ylim = c(0, 1))      

ggsave(file.path(fig_dir, "meta_coverage_scaled.pdf"), width=10, height=8)  

#
# combine triplicates + facet_wrap metacoverages
#
tidy_avg <- map_dfr(levels(tidy_meta_cov$treatment), function(x) {
  df <- tidy_meta_cov %>% dplyr::filter(treatment == x) %>%
    dplyr::select(-value) %>%
    spread(sample_name, scaled_value)
  df <- df %>% dplyr::mutate(avg = rowMeans(df[, c(4,5,6)]))  %>%
    dplyr::select(index, regions, treatment, avg)
}) %>%
  dplyr::mutate(log_avg=log2(avg+1))

as.data.frame(colData(dds_cds_by_gene)) %>% 
  group_by(treatment) %>% 
  summarise(mean_records = mean(records),
            mean_scaled = mean(records/ sizeFactor))

# range [0, 0.3]
ggplot(tidy_avg, aes(x=index, y=avg, fill=treatment, group=treatment)) +
    geom_area(stat="identity", alpha=0.7) + 
    coord_cartesian(ylim = c(0, 0.3)) + #2.26
    theme_minimal() +
    labs(y="mean coverage (scaled)") +
    facet_wrap(~treatment) +
    theme(legend.position="none")+
          #panel.grid.major = element_blank()) +
    #geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
    #geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                       labels=c("5'UTR", "CDS", "3'UTR"))                       
ggsave(file.path(fig_dir, "meta_coverage_scale_facet_wrap_0.3.pdf"), width=6.5, height=4)    

# range [0, 0.5]
ggplot(tidy_avg, aes(x=index, y=avg, fill=treatment, group=treatment)) +
    geom_area(stat="identity", alpha=0.8) + 
    coord_cartesian(ylim = c(0, 0.5)) + #2.26
    theme_minimal() +
    labs(y="mean coverage (scaled)") +
    facet_wrap(~treatment) +
    theme(legend.position="none")+
          #panel.grid.major = element_blank()) +
    #geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
    #geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                       labels=c("5'UTR", "CDS", "3'UTR"))                       
ggsave(file.path(fig_dir, "meta_coverage_scale_facet_warp_0.5.pdf"), width=6.5, height=4) 

# range [0, 2.26]
ggplot(tidy_avg, aes(x=index, y=avg, fill=treatment, group=treatment)) +
    geom_area(stat="identity", alpha=0.8) + 
    coord_cartesian(ylim = c(0, 2.26)) + #2.26
    theme_minimal() +
    labs(y="mean coverage (scaled)") +
    facet_wrap(~treatment) +
    theme(legend.position="none")+
          #panel.grid.major = element_blank()) +
    #geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
    #geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                       labels=c("5'UTR", "CDS", "3'UTR"))                       
ggsave(file.path(fig_dir, "meta_coverage_scale_facet_warp_2.pdf"), width=6.5, height=4) 

## log scale
ggplot(tidy_avg, aes(x=index, y=log_avg, fill=treatment, group=treatment)) +
    geom_area(stat="identity", alpha=0.6) + 
    coord_cartesian(ylim = c(0, 1.7)) + #2.26
    theme_minimal() +
    labs(y="log2 (norm mean coverage +1)") +
    facet_wrap(~treatment)
    theme(legend.position="top",
          panel.grid.major = element_blank()) +
    geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
    geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                       labels=c("5'UTR", "CDS", "3'UTR")) 
ggsave(file.path(fig_dir, "meta_coverage_scale_superpose_log.pdf"), width=6.5, height=4)  

tmp <- tidy_meta_cov %>% 
  dplyr::filter(sample_name %in% c("DCH1_untreated_1", "DCH7_DOX-pulse_1", "DCH4_IFNg_1", "DCH10_DOX-pulse_IFNg_1"))
ggplot(tmp, aes(x=index, y=scaled_value, fill=treatment)) +
    geom_area(stat="identity", alpha=0.6) + 
    coord_cartesian(ylim = c(0, 0.5)) + #2.26
    theme_minimal() +
    facet_wrap(~treatment) +
    labs(y="mean coverage (scaled)") +
    theme(legend.position="top",
          panel.grid.major = element_blank()) +
    geom_vline(aes(xintercept=bins_sum[1]+1), alpha=0.4, linetype="dashed", color="grey50") +
    geom_vline(aes(xintercept=bins_sum[2]+1), linetype="dashed", color="grey50") +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                       labels=c("5'UTR", "CDS", "3'UTR")) 
ggsave(file.path(fig_dir, "test.pdf"), width=6.5, height=4)  

#
# (6) sense and antisense strand
#
reads <- bplapply(sample_info$bam_files, function(f) {
  bam_file <- BamFile(file = f, yieldSize = yieldSize)
  pc <- getPsiteCoordinates(bam_file, bestpsite = best_p_site,
                            anchor = anchor)
  pc.sub <- pc[pc$qwidth %in% reads_len]
})
names(reads) <- sample_info$sample_name

# strand plot
strand <- bplapply(reads, function(pc, CDS) {
    if (length(intersect(seqlevelsStyle(pc), seqlevels(CDS))) ==
        0) {
        seqlevelsStyle(pc) <- seqlevelsStyle(CDS)[1]
    }
    ol <- findOverlaps(pc, CDS, ignore.strand = FALSE)
    a <- unique(queryHits(ol))
    reads.rev <- ribosomeProfilingQC:::switch.strand(pc)
    ol.anti <- findOverlaps(reads.rev, CDS, ignore.strand = FALSE)
    b <- unique(queryHits(ol.anti))
    b <- b[!b %in% a]
    per <- c(sense = length(a), antisense = length(b))/length(pc) * 100
}, CDS)
strand <- do.call(rbind, strand) 

strand %>% as.data.frame() %>%
  tibble::rownames_to_column(var="sample_name") %>%
  gather(key=direction, value=frequency, -sample_name) %>%
  dplyr::left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name))) %>%
  dplyr::mutate(format=format(frequency, digits=1)) %>%
  ggplot(aes(x=direction, y=frequency)) +
    geom_bar(stat="identity", width=0.3) +
    facet_wrap( ~ sample_name, nrow=4) +
    theme_bw() +
    geom_text(aes(label=format), vjust=0.3, size=3) +
    labs(y="frequency (%)", x="")
ggsave(file.path(fig_dir, "strand_plot.pdf"))    

#
# (7) viz the supreseed near-start-codon enrichment
#
load(file.path(pkg_dir, "data", "cov_utr5_cds_utr3.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))

cov_utr5 <- cov_utr5_cds_utr3$cov_utr5
cov_cds  <- cov_utr5_cds_utr3$cov_cds
cov_utr3  <- cov_utr5_cds_utr3$cov_utr3
bins <- c(UTR5 = 100, CDS = 1000, UTR3=100)
meta_cov <- bplapply(sample_info$bam_files, function(x) {
  x <- basename(x)
  .metaPlot(UTR5coverage=cov_utr5, CDScoverage=cov_cds, UTR3coverage=cov_utr3,
            sample=x, xaxis="RPFs", bins=bins) %>%
    as_tibble() %>%
    rownames_to_column(var="index") %>%
    dplyr::mutate(index = as.numeric(index))  %>%
    add_column(sample_name = str_replace(x, ".bam", "")) %>%
    add_column(regions=rep(names(bins), bins))
})

# ggplot(df, aes(x=x, y=y, fill=status)) + geom_area(stat="identity")
