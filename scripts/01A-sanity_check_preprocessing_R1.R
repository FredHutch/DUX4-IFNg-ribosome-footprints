# sanity_check_preprocessing_R1.sanity_check_preprocessing_R1
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

# use riboProfilingQC to check the fragment length
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint"
scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"


# 
CDS <- prepareCDS(txdb)
bam_file <- list.files(file.path(scratch_dir, "run_2", "bam"), pattern=".sortedByCoord.out.bam$",
                       full.names=TRUE)

x <- BamFile(bam_file[1])
p_site <- estimatePsite(x, CDS, genome)
pc <- getPsiteCoordinates(x, bestpsite = p_site)
read_length_freq <- summaryReadsLength(pc, widthRange = c(25:39), plot=FALSE)

read_length_freq <- bplapply(bam_file, function(x) {
  bam_file <- BamFile(x)
  p_site <- estimatePsite(bam_file, CDS, genome)
  pc <- getPsiteCoordinates(bam_file, bestpsite = p_site)
  read_length_freq <- summaryReadsLength(pc, widthRange = c(25:39), plot=FALSE)
})
names(read_length_freq) <- str_replace(basename(bam_file), "_Aligned.sortedByCoord.out.bam", "")

# tidy
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
gg <- ggplot(length_freq, aes(x=length, y=Freq)) +
  geom_bar(stat="identity", width=0.7) +
  theme_bw() +
  facet_wrap( ~ sample_name, nrow=4) +
  labs(x="Read length", y="Frequency")

ggsave(file.path(pkg_dir, "figures", "QC", "freqment_size_frequency.2.pdf"))  

