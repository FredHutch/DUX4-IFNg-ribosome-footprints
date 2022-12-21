# make_rnaseq_ucsc_track.R
# ml deepTools
# ml SAMtools/1.11-GCC-10.2.0
library(BiocParallel)
library(tidyverse)
bp_param <- MulticoreParam(workers = 12L)
register(bp_param, default=TRUE)

bam_dir <- "/fh/scratch/delete90/tapscott_s/DCH_Paired_RiboSeq_RNA-Seq_BAM_files"
out_dir <- "/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/hg38_MB135iDUX4_IFNg_RNAseq"
file.exists(out_dir)

Bam_Files <- list.files(bam_dir, pattern=".bam$", full.names=TRUE)
sample_info <- data.frame(sample_name = str_replace(basename(Bam_Files), ".bam", ""),
                          bam_files = Bam_Files,
                          bw_files = file.path(out_dir, str_replace(basename(Bam_Files), ".bam", ".bw")))

#
# (1) make sort and index
#
setwd(bam_dir)
bplapply(sample_info$bam_files, function(x) {
    sample_name <- str_replace(basename(x), ".bam", "")
    sorted <- paste0(x, ".sorted")
    cmt_1 <- sprintf("samtools sort -@ 4 -T %s %s -o %s", sample_name, x, sorted)
    system(cmt_1)
    cmt_2 <- sprintf("mv %s %s", sorted, x)
    system(cmt_2)
    cmt_3 <- sprintf("samtools index %s", x)
    system(cmt_3)
})

#
#  (2) convert BAM to bigWig
#

bplapply(1:nrow(sample_info), function(i) {
  # deepTools::bamCoverage
  cmt <- sprintf("bamCoverage --binSize 10 -b %s -o %s",
                 sample_info$bam_files[i], sample_info$bw_files[i])
  system(cmt)
})

#
# make track file: use xfile URL then change to S3
#
track_name <- "MB135iDUX4_IFNg_RNAseq"
track_file_name <- file.path(out_dir, paste0(track_name, "_HubTrack_s3.txt"))
createHubTrackLine_bigwig(bwDir=out_dir, trackName=track_name,
                          TrackFileName=track_file_name,
                          shortLabel=track_name,
                          longLabel = "MB135 iDUX4 pulse and IFNg treated RNA-seq", 
                          url="https://fh-pi-tapscott-s-eco-public.s3.us-west-2.amazonaws.com/ucsc-tracks/bigWig")
createHubTrackLine_bigwig(bwDir=out_dir, trackName=track_name,
                          TrackFileName=file.path(out_dir, paste0(track_name, "_HubTrack_xfile.txt")),
                          shortLabel=track_name,
                          longLabel = "MB135 iDUX4 pulse and IFNg treated RNA-seq",
                          url="http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/bigWig")

# in case you need to switch URL
# convert from old S3 to new S3 URL
tmp <- "https://fh-pi-tapscott-s-eco-public.s3.us-west-2.amazonaws.com"
new <- "https://fh-pi-tapscott-s-eco-public.s3.us-west-2.amazonaws.com/ucsc-tracks"
x <- "hg38.trackDb.txt"
tx  <- readLines(x)
tx2  <- gsub(pattern = tmp, replace = new, x = tx)
writeLines(tx2, con=x)

