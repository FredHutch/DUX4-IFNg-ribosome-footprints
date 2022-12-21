# call_preprocessing_R1.R / remember to ml SAMtools/1.10-GCCcore-8.3.0
# (1) runs preprocessing_R1.sh
# (2) merge bam files run1 + run2
# (3) merge bam files by treatment (combine all triplicates)


# Run FASTQ folders
pkg_dir <- "/fh/fast/tapscott_s/CompBio/Ribo-seq/hg38.DUX4.IFN.ribofootprint.2"
run_1 <- "/shared/ngs/illumina/dhamm/201107_VH00319_7_AAAHJGCM5"
run_2 <- "/shared/ngs/illumina/dhamm/201110_VH00319_8_AAAGV7HM5"
library(tidyverse)
setwd(file.path(pkg_dir, "scripts"))

# RUN1
run <- "run_1"
fastq_files <- list.files(file.path(run_1, "Unaligned/Project_dhamm"),
                          pattern=".fastq.gz", full.names = TRUE )
sample_name <- unique(sapply(str_split(basename(fastq_files), "_S"), "[[", 1))
map(sample_name, function(x) {
  fqs <- str_subset(fastq_files, x)
  fq1 <- str_subset(fqs, "R1")
  cmt <- sprintf("sbatch -n4 ./000-preprocessing_R1.sh %s %s", fq1, run)
  system(cmt)
})  

# RUN2
run <- "run_2"
fastq_files <- list.files(file.path(run_2, "Unaligned/Project_dhamm"),
                          pattern=".fastq.gz", full.names = TRUE )
sample_name <- unique(sapply(str_split(basename(fastq_files), "_S"), "[[", 1))
map(sample_name, function(x) {
  fqs <- str_subset(fastq_files, x)
  fq1 <- str_subset(fqs, "R1")
  cmt <- sprintf("sbatch -n4 ./000-preprocessing_R1.sh %s %s", fq1, run)
  system(cmt)
})  

#
# merge run 1 and run 2
#
library(tidyverse)
library(BiocParallel)
bp_param=MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

scratch_dir <- "/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
des_dir <- file.path(scratch_dir, "bam", "merged_bam_runs")
bam_file_1 <- data.frame(
  bam_file = list.files(file.path(scratch_dir, "run_1", "bam"), pattern=".sortedByCoord.out.bam$",
                        full.names=TRUE)) %>%
  dplyr::mutate(sample_name =  str_replace(basename(bam_file), "_Aligned.sortedByCoord.out.bam", ""))

bam_file_2 <- data.frame(
  bam_file = list.files(file.path(scratch_dir, "run_2", "bam"), pattern=".sortedByCoord.out.bam$",
                        full.names=TRUE)) %>%
  dplyr::mutate(sample_name =  str_replace(basename(bam_file), "_Aligned.sortedByCoord.out.bam", ""))

merge_df <- left_join(bam_file_1, bam_file_2, by="sample_name") %>%
  dplyr::mutate(out_bam = file.path(des_dir, paste0(sample_name, ".bam")))
  
bplapply(1:12, function(i){
  cmt <- sprintf("samtools merge -f %s %s %s", merge_df$out_bam[i], merge_df$bam_file.x[i], merge_df$bam_file.y[i])
  system(cmt)
  # SAMtools::sort and index
  tmp <- file.path(des_dir, str_replace(basename(merge_df$out_bam[i]), ".bam", "_sort.bam"))
  system(sprintf("samtools sort -@ 2 %s > %s", merge_df$out_bam[i], tmp))
  system(sprintf("mv %s %s", tmp, merge_df$out_bam[i]))
  system(sprintf("samtools index %s", merge_df$out_bam[i]))
})                  

#
# merged bam by triplicates (total four treatments)
#
bam_dir <- file.path(scratch_dir, "bam", "merged_bam_runs")
des_dir <- file.path(scratch_dir, "bam", "merged_bam_triplicates")
bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)
merge_df <- data.frame(
  bam_files = bam_files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)) %>%
    dplyr::mutate(sample_name = str_replace(basename(bam_files), ".bam", ""),
                  treatment = str_replace(str_sub(sample_name, start=1L, end=-3L), "[^_]+_", "")) %>%
    group_by(treatment) %>%
    summarise(in_bam = paste(bam_files, collapse=" ")) %>%
    dplyr::mutate(out_bam = file.path(des_dir, paste0(treatment, ".bam")))    

bplapply(1:4, function(i) {
  cmt <- sprintf("samtools merge %s %s", merge_df$out_bam[i], merge_df$in_bam[i])
  system(cmt)
  # SAMtools::sort and index
  tmp <- file.path(des_dir, str_replace(basename(merge_df$out_bam[i]), ".bam", "_sort.bam"))
  system(sprintf("samtools sort -@ 2 %s > %s", merge_df$out_bam[i], tmp))
  system(sprintf("mv %s %s", tmp, merge_df$out_bam[i]))
  system(sprintf("samtools index %s", merge_df$out_bam[i]))
} ) 