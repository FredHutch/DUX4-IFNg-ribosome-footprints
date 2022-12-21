# 10C-manuscript-figures-reading-frame-periodicity.R

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
library(plyranges)
library(hg38.HomoSapiens.Gencode.v35)
txdb <- hg38.HomoSapiens.Gencode.v35
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(BiocParallel)
bp_param=MulticoreParam(workers = 4L)
register(bp_param, default=TRUE)

#
# define parameters and load datasets
#
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
fig_dir <- file.path(pkg_dir, "manuscript", "figures", "QC")
source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "tools_assign_utr_frame.R"))
load(file.path(pkg_dir, "data", "p_sites.rda")) # p-sites of read length 26:29
load(file.path(pkg_dir, "data", "cov_utr5_cds_utr3.rda"))
load(file.path(pkg_dir, "data", "dds_cds_by_gene.rda"))
sample_info <- as(colData(dds_cds_by_gene), "data.frame")
read_length <- c(26:29)
best_offsite <- 13
sample_info <- as.data.frame(colData(dds_cds_by_gene))
CDS <- ribosomeProfilingQC::prepareCDS(txdb)

#
# assigne reading frame to p-sites
#
rf_psites <- bplapply(p_sites, ribosomeProfilingQC::assignReadingFrame, CDS)
pc_by_regions <- bplapply(rf_psites, .assign_utr_frame, CDS, txdb)

frame_freq_by_regions <- map_dfr(pc_by_regions, function(pc) {
  tb <- map_dfr(pc, function(x) {
    table(x$readingFrame) / length(x)
  }, .id="region") %>%
    gather(key=Frame, value=Percent, -region)
}, .id="sample_name") %>%
  dplyr::mutate(Percent=as.numeric(Percent)) %>%
  dplyr::mutate(region = case_when(region == "cds" ~ "CDS",
                                   region == "utr5" ~ "5' UTR",
                                   region == "utr3" ~ "3' UTR")) %>%
  dplyr::mutate(region = factor(region, levels=c("5' UTR", "CDS", "3' UTR"))) %>%
  dplyr::mutate(treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
  dplyr::mutate(treatment = factor(treatment, levels=c("untreated", "DOX-pulse", "IFNg", "DOX-pulse_IFNg"))) %>%
  dplyr::arrange(treatment) %>% 
  dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))

#
# viz by p-site distribution on frame
#
library(viridis)
frame_freq_by_regions %>% #dplyr::filter(region == "five_utr") %>%
  ggplot(aes(x=Frame, y=sample_name)) +
    geom_tile(aes(fill=Percent)) +
    facet_wrap(~region) +
    theme_classic() +
    scale_fill_viridis() +
    labs(y="")
ggsave(file.path(fig_dir, "frame_frequency_by_regions.pdf"), height=3, width=5)  

# take average over treatment
frame_freq_by_regions  %>% 
  spread(key=Frame, value=Percent) %>% 
  group_by(treatment, region) %>%
  summarise(`0` = mean(`0`), `1` = mean(`1`), `2`=mean(`2`)) %>%
  gather(key=Frame, value=Percent, -treatment, -region) %>% 
  ggplot(aes(x=Frame, y=treatment)) +
    geom_tile(aes(fill=Percent)) +
    facet_wrap(~region) +
    theme_classic() +
    scale_fill_viridis() +
    labs(y="") +
    theme(legend.key.size = unit(0.3, "cm"), legend.title=element_blank())

ggsave(file.path(fig_dir, "frame_frequency_by_regions_average_over_treatment.pdf"), height=1.7, width=4)  
   
#
# viz frame frequency by meta-gene plot
#

# get p-sites by frames
pcs <- pc_by_regions[[2]] 
frame <- 0
anchor = "5end"
ext = 5000
level = "tx"
regions <- c("utr5", "cds", "utr3")
# frame 0, utr5/cds/utr3
pcs_frame0 <- lapply(pcs, function(x) 
     x %>% plyranges::filter(readingFrame == 0))
regions <- names(pcs_frame0)     
region = regions[1]
reads <- pcs_frame0$utr5
reads_utr5_frame0 <- lapply(pc_by_regions, function(x) {
  pc <- x$utr5
  pc %>% plyranges::filter(readingFrame == 0)
})

.coverageDepth_per_region_by_frame <- function(pc_by_regions, txdb, anchor="5end", ext=5000, level="tx", frame=0) {
    # this function mimics coverageDepth() from ribosomeProfilingQC
    regions <- names(pc_by_regions[[1]])
    cvgs_per_region <- bplapply(regions, function(region) {
        cvgs <- list()
        reads <- lapply(pc_by_regions, function(x) {
          x[[region]] %>% plyranges::filter(readingFrame == frame)
        })
        # if utr5 or utr3, exclude reads that overlaps with CDS
        cd <- .getCvgs(reads, txdb, level, anchor, region, ext)
        cvgs[["RPFs"]] <- cd
        return(cvgs)
    })
    names(cvgs_per_region) <- regions
    cvgs_per_region
}

cov_frame0 <- .coverageDepth_per_region_by_frame(pc_by_regions, txdb, frame=0)
cov_frame1 <- .coverageDepth_per_region_by_frame(pc_by_regions, txdb, frame=1)
cov_frame2 <- .coverageDepth_per_region_by_frame(pc_by_regions, txdb, frame=2)

# build meta coverage
# use .metaplot() from tools.R
.get_meta_cov <- function(cov_frame, bins, reading_frame) {
  meta_cov <- bplapply(sample_info$sample_name, function(x) {
    .metaPlot(UTR5coverage=cov_frame[["utr5"]], 
              CDScoverage=cov_frame[["cds"]], UTR3coverage=cov_frame[["utr3"]],  
              sample=x, xaxis="RPFs", bins=bins) %>%
      as_tibble() %>%
      rownames_to_column(var="index") %>%
      dplyr::mutate(index = as.numeric(index))  %>%
      add_column(regions=rep(names(bins), bins), sample_name=x, reading_frame=reading_frame) 
  })
  # tidy meta_cov
  do.call(rbind, meta_cov) %>%
    dplyr::mutate(regions = factor(regions, levels=names(bins))) %>%
    dplyr::left_join(dplyr::select(sample_info, sample_name, treatment, sizeFactor), by="sample_name") %>%
    dplyr::arrange(treatment) %>%
    dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
}

bins <- c(UTR5 = 15, CDS = 100, UTR3 = 15)
meta_cov_frame0 <- .get_meta_cov(cov_frame0, bins, reading_frame=0)
meta_cov_frame1 <- .get_meta_cov(cov_frame1, bins, reading_frame=1) 
meta_cov_frame2 <- .get_meta_cov(cov_frame2, bins, reading_frame=2) 
meta_cov_by_frame <- list(frame0=meta_cov_frame0, frame1=meta_cov_frame1, frame2=meta_cov_frame2)

save(meta_cov_by_frame, file=file.path(pkg_dir, "data", "meta_cov_by_frame.rda"))

# viz by frames
do.call(rbind, meta_cov_by_frame) %>% # normalized by sizefactor
  #dplyr::filter(sample_name == sample_info$sample_name[1]) %>% 
  dplyr::mutate(index = index + 0.3*reading_frame) %>% 
  dplyr::mutate(norm_value = value / sizeFactor) %>%
  dplyr::mutate(reading_frame = factor(as.character(reading_frame), levels=c("0", "1", "2"))) %>%
  ggplot(aes(x=index, y=norm_value, fill=reading_frame)) +
    geom_bar(stat="identity", alpha=0.7, width=0.3) +
    theme_bw() +
    scale_fill_brewer(palette="Dark2", name="frame") +
    scale_y_sqrt(breaks=c(0, 0.5, 2)) +
    labs(y="mean coverage") +
    scale_x_continuous(breaks=c(1, bins[1]+1, bins[2]+bins[1]+1), 
                      labels=c("5' UTR", "CDS", "3' UTR")) +
    facet_wrap(~sample_name, nrow=4, scale="free_x") +    
    theme(legend.position=c(0.3, 0.95), legend.title = element_text(size=10),
          panel.grid.major.x=element_blank(),
          axis.title.x = element_blank(), legend.key.size = unit(0.3, 'cm'),
          panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank()) 

ggsave(file=file.path(fig_dir, "reading_frame_periodicity_all_samples_norm.pdf"), width=10, height=8)

#
# take average of the triplicates of each treatment; generate four meta-plots
#
treatments <- levels(sample_info$treatment)

cov_per_treatment <- bplapply(meta_cov_by_frame, function(frame_cov) {
    map_dfr(treatments, function(x) {
        frame_cov %>% dplyr::filter(treatment == x) %>%
          dplyr::mutate(norm_value = value / sizeFactor) %>%
          dplyr::select(-value, -sizeFactor) %>%
          spread(sample_name, norm_value) %>%
          dplyr::mutate(avg_norm_value = rowMeans(.[, 5:7]))
    })
})  

do.call(rbind, cov_per_treatment) %>% # normalized by sizefactor
  dplyr::mutate(index = index + 0.3*reading_frame) %>% 
  dplyr::mutate(reading_frame = factor(as.character(reading_frame), levels=c("0", "1", "2"))) %>%
  ggplot(aes(x=index, y=avg_norm_value, fill=reading_frame)) +
    geom_bar(stat="identity", alpha=0.8, width=0.3) +
    theme_bw() +
    scale_fill_brewer(palette="Dark2", name="frame") +
    scale_y_sqrt(breaks=c(0, 0.5, 2)) +
    labs(y="p-site mean coverage") +
    scale_x_continuous(breaks=c(1, bins[1]+1, bins[2]+bins[1]+1), 
                      labels=c("5' UTR", "TSS", "3' UTR")) +
    facet_wrap(~treatment, nrow=2, scale="free_x") +    
    theme(legend.position=c(0.95, 0.88), legend.title = element_text(size=10),
          legend.background=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x = element_text(size=7),
          axis.title.x = element_blank(), legend.key.size = unit(0.25, 'cm'),
          panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank()) 

ggsave(file=file.path(fig_dir, "reading_frame_periodicity_by_treatment_norm.pdf"), width=6, height=3)




################## the old stuff that might be wrong #####################
# 
# get p_sites from each reading frame
#
.extract_p_sites_by_frame <- function(rf_psites, reading_frame=0)  {
    bplapply(rf_psites, function(pc) {
        pc %>% plyranges::filter(readingFrame == reading_frame)
    })
}


frame_0 <- .extract_p_sites_by_frame(rf_psites, reading_frame=0)
frame_1 <- .extract_p_sites_by_frame(rf_psites, reading_frame=1)
frame_2 <- .extract_p_sites_by_frame(rf_psites, reading_frame=2)
source(file.path(pkg_dir, "scripts", "tools_getCvg.R"))

# paramters
anchor = "5end"
ext = 5000
level = "tx"
regions <- c("utr5", "cds", "utr3")

# get coverage for each frame
.coverageDepth_per_region <- function(reads, txdb, anchor="5end", ext=5000, level="tx") {
    # this function mimics coverageDepth() from ribosomeProfilingQC
    regions <- c("utr5", "cds", "utr3")
    cvgs_per_region <- bplapply(regions, function(region) {
        cvgs <- list()
        # if utr5 or utr3, exclude reads that overlaps with CDS
        cd <- .getCvgs(reads, txdb, level, anchor, region, ext)
        cvgs[["RPFs"]] <- cd
        return(cvgs)
    })
    names(cvgs_per_region) <- regions
    cvgs_per_region
}

cov_frame0 <- .coverageDepth_per_region(frame_0, txdb)
cov_frame1 <- .coverageDepth_per_region(frame_1, txdb)
cov_frame2 <- .coverageDepth_per_region(frame_2, txdb)

# use .metaplot() from tools.R
.get_meta_cov <- function(cov_frame, bins, reading_frame) {
  meta_cov <- bplapply(sample_info$sample_name, function(x) {
    .metaPlot(UTR5coverage=cov_frame[["utr5"]], 
              CDScoverage=cov_frame[["cds"]], UTR3coverage=cov_frame[["utr3"]],  
              sample=x, xaxis="RPFs", bins=bins) %>%
      as_tibble() %>%
      rownames_to_column(var="index") %>%
      dplyr::mutate(index = as.numeric(index))  %>%
      add_column(regions=rep(names(bins), bins), sample_name=x, reading_frame=reading_frame) 
  })
  # tidy meta_cov
  do.call(rbind, meta_cov) %>%
    dplyr::mutate(regions = factor(regions, levels=names(bins))) %>%
    dplyr::left_join(dplyr::select(sample_info, sample_name, treatment, sizeFactor), by="sample_name") %>%
    dplyr::arrange(treatment) %>%
    dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name)))
}

bins <- c(UTR5 = 15, CDS = 100, UTR3 = 15)
meta_cov_frame0 <- .get_meta_cov(cov_frame0, bins, reading_frame=0)
meta_cov_frame1 <- .get_meta_cov(cov_frame1, bins, reading_frame=1) 
meta_cov_frame2 <- .get_meta_cov(cov_frame2, bins, reading_frame=2) 
meta_cov_by_frame <- list(frame0=meta_cov_frame0, frame1=meta_cov_frame1, frame2=meta_cov_frame2)
save(meta_cov_by_frame, file=file.path(pkg_dir, "data", "meta_cov_by_frame.rda"))



# 
# viz meta-coverage by reading frames (0, 1, 2) 
#

# a.  just one sample, not normalized
do.call(rbind, meta_cov_by_frame) %>% # normalized by sizefactor
  #dplyr::filter(sample_name == sample_info$sample_name[1]) %>% 
  dplyr::mutate(index = index + 0.3*reading_frame) %>% 
  dplyr::mutate(norm_value = value / sizeFactor) %>%
  dplyr::mutate(reading_frame = factor(as.character(reading_frame), levels=c("0", "1", "2"))) %>%
  ggplot(aes(x=index, y=norm_value, fill=reading_frame)) +
    geom_bar(stat="identity", alpha=0.7, width=0.3) +
    theme_bw() +
    scale_fill_brewer(palette="Dark2", name="frame") +
    scale_y_sqrt(breaks=c(0, 0.5, 2)) +
    labs(y="mean coverage") +
    scale_x_continuous(breaks=c(1, bins[1]+1, bins[2]+bins[1]+1), 
                      labels=c("5' UTR", "TSS", "3' UTR")) +
    facet_wrap(~sample_name, nrow=4, scale="free_x") +    
    theme(legend.position=c(0.3, 0.95), legend.title = element_text(size=10),
          panel.grid.major.x=element_blank(),
          axis.title.x = element_blank(), legend.key.size = unit(0.3, 'cm'),
          panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank()) 


ggsave(file=file.path(fig_dir, "reading_frame_periodicity_all_samples_norm.pdf"), width=10, height=8)

# b. take average of the triplicates of each treatment; generate four meta-plots
treatments <- levels(sample_info$treatment)

cov_per_treatment <- bplapply(meta_cov_by_frame, function(frame_cov) {
    map_dfr(treatments, function(x) {
        frame_cov %>% dplyr::filter(treatment == x) %>%
          dplyr::mutate(norm_value = value / sizeFactor) %>%
          dplyr::select(-value, -sizeFactor) %>%
          spread(sample_name, norm_value) %>%
          dplyr::mutate(avg_norm_value = rowMeans(.[, 5:7]))
    })
})  

do.call(rbind, cov_per_treatment) %>% # normalized by sizefactor
  dplyr::mutate(index = index + 0.3*reading_frame) %>% 
  dplyr::mutate(reading_frame = factor(as.character(reading_frame), levels=c("0", "1", "2"))) %>%
  ggplot(aes(x=index, y=avg_norm_value, fill=reading_frame)) +
    geom_bar(stat="identity", alpha=0.8, width=0.3) +
    theme_bw() +
    scale_fill_brewer(palette="Dark2", name="frame") +
    scale_y_sqrt(breaks=c(0, 0.5, 2)) +
    labs(y="p-site mean coverage") +
    scale_x_continuous(breaks=c(1, bins[1]+1, bins[2]+bins[1]+1), 
                      labels=c("5' UTR", "TSS", "3' UTR")) +
    facet_wrap(~treatment, nrow=2, scale="free_x") +    
    theme(legend.position=c(0.95, 0.88), legend.title = element_text(size=10),
          legend.background=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.title.x = element_blank(), legend.key.size = unit(0.3, 'cm'),
          panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank()) 

ggsave(file=file.path(fig_dir, "reading_frame_periodicity_by_treatment_norm.pdf"), width=8, height=4)

