# 10A-manuscript-figures-QC.R
# This script makes metacoverage for each treatment (merge p-sites of the triplicates) 
# 


#
# load library
#
library(ribosomeProfilingQC)
library(tidyverse)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(latex2exp)
library(ggthemes)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# define parameters
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
source(file.path(pkg_dir, "scripts", "tools.R"))
load(file.path(pkg_dir, "data", "p_sites.rda")) # p-sites of read length 26:29
load(file.path(pkg_dir, "data", "cov_utr5_cds_utr3.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))
sample_info <- as.data.frame(colData(dds_cds_by_gene))
CDS <- ribosomeProfilingQC::prepareCDS(txdb)
# assigne reading frame to p-sites: assignReadingFrame(pc_sub, CDS)


#
# meta-plot
# 
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

# (a) clos-up meta-plot (one track per treatment)
meta_cov_scale <- bplapply(meta_cov, function(x) {
  x %>% dplyr::mutate(scaled_value = value / sizeFactors(dds_cds_by_gene)[sample_name[1]]) 
})

tidy_meta_cov <- do.call(rbind, meta_cov_scale) %>%
  dplyr::mutate(regions = factor(regions, levels=names(bins))) %>%
  dplyr::left_join(dplyr::select(sample_info, sample_name, treatment), by="sample_name") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
bins_sum <- cumsum(bins)

tidy_avg <- map_dfr(levels(tidy_meta_cov$treatment), function(x) {
  df <- tidy_meta_cov %>% dplyr::filter(treatment == x) %>%
    dplyr::select(-value) %>%
    spread(sample_name, scaled_value)
  df <- df %>% dplyr::mutate(avg = rowMeans(df[, c(4,5,6)]))  %>%
    dplyr::select(index, regions, treatment, avg)
})

treatment_labs <- c(untreated="untreated", DOX_pulse="DOX_pulse", 
                    IFNg=TeX(r'(IFN{$\gamma})'), DOX_pulse_IFNg=TeX(r'(DOX\_pulse+IFN{$\gamma})'))
levels(tidy_avg$treatment) <- treatment_labs

close_up_cov <- tidy_avg %>% dplyr::filter(index < 131)
ggplot(close_up_cov, aes(x=index, y=avg, group=treatment)) +
    geom_area(stat="identity", alpha=0.5, fill="gray50") + 
    geom_line(size=1, color="gray25") +
    theme_minimal() +
    labs(y="mean coverage") +
    facet_wrap(~treatment, nrow=4, labeller = label_parsed) +
    theme(legend.position="none", axis.title.x = element_blank(), axis.title.y=element_blank()) +
    scale_x_continuous(breaks=c(101), labels=c("TSS")) +
    scale_y_sqrt(breaks=c(0, 0.5, 2))                                     
ggsave(file.path(fig_dir, "meta-coverage", "meta-coverage-close-up-scale_y_sqrt.pdf"), width=1.8, height=4)   

ggplot(tidy_avg, aes(x=index, y=avg, group=treatment)) +
    geom_area(stat="identity", alpha=0.5, fill="gray50") + 
    geom_line(size=0.4, color="gray25") +
    theme_minimal() +
    labs(y="mean coverage") +
    facet_wrap(~treatment, nrow=4, labeller = label_parsed) +
    theme(legend.position="none", axis.title.x = element_blank(), 
          panel.grid.minor.x=element_blank()) +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                        labels=c("5' UTR", "CDS", "3' UTR")) +  
    scale_y_sqrt(breaks=c(0, 0.5, 2))                                     
ggsave(file.path(fig_dir, "meta-coverage", "meta-coverage-all-scale_y_sqrt.pdf"), width=5, height=5)   

ggplot(tidy_avg, aes(x=index, y=avg, group=treatment)) +
    geom_area(stat="identity", alpha=0.8, fill="gray10") + 
    theme_minimal() +
    labs(y="mean coverage") +
    facet_wrap(~treatment, nrow=4, labeller=label_parsed) +
    theme(legend.position="none", axis.title.x = element_blank(), 
          panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank()) +
    scale_x_continuous(breaks=c(bins_sum[1]/2, mean(bins_sum[c(1,2)]), mean(bins_sum[c(2,3)])), 
                        labels=c("5' UTR", "CDS", "3' UTR")) +
    scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2))                        
ggsave(file.path(fig_dir, "meta-coverage", "meta-coverage-all.pdf"), width=5, height=5)   


#
# p-site reading frame periodicity on CDS (merge tripliecates)
#

