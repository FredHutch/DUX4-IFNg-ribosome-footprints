#!/bin/bash
#./preprocess_R1.sh
#SBATCH -n4 -t 1-0 -p campus-new --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

# pipeline for ribo-seq forward read only
# called by call_preprocessing_R1.R

fq1=$1
run=$2

rRNA_ref="/fh/fast/tapscott_s/CompBio/genome_reference/rRNA.homo_sapiens/Sequence/Bowtie2Index/rRNA.homo_sapiens"
tsoRNA_ref="/fh/fast/tapscott_s/CompBio/genome_reference/tsoRNA.homo_sapiens/Sequence/Bowtie2Index/tsoRNA.homo_sapiens"
genome_fasta="/fh/fast/tapscott_s/CompBio/genome_reference/GRCh38/Sequence/WholeGenomeFasta/GRCh38.p13.genome.fa"
gtf_file="/fh/fast/tapscott_s/CompBio/genome_reference/GRCh38/Annotation/gencode.v35.annotation.gtf"
star_index="/fh/fast/tapscott_s/CompBio/genome_reference/GRCh38/Sequence/STARIndex"
# define folders
scratch_dir="/fh/scratch/delete90/tapscott_s/hg38.DUX4.IFN.ribofootprint.R1"
trim_dir=$scratch_dir/$run/trimmed
filtered_dir=$scratch_dir/$run/filtered
fastqc_dir=$scratch_dir/$run/fastqc
bam_dir=$scratch_dir/$run/bam

mkdir -p $trim_dir
mkdir -p $fastqc_dir
mkdir -p $filtered_dir
mkdir -p $bam_dir

# module load
ml purge
ml FastQC/0.11.9-Java-11
ml cutadapt/2.7-GCCcore-8.3.0-Python-3.7.4
ml Bowtie2/2.4.1-GCCcore-8.3.0
ml SAMtools/1.10-GCCcore-8.3.0
ml STAR/2.7.6a-foss-2019b

#
# (1) cutadpt
#
base_fq1="$(basename -- $fq1)"
sample_name="${base_fq1//_S*.fastq.gz}"
echo $sample_name
cd $trim_dir

echo trimming
cutadapt --cores=4 -O 13 -a N{4}TGGAATTCTCGGGTGCCAAGG -u 4 -m 25 --untrimmed-output ${sample_name}_R1_untrimmed.fastq.gz -o ${sample_name}_R1_trimmed.fastq.gz $fq1 > ${sample_name}_cutadapt_log.txt

fq=$trim_dir/${sample_name}_R1_trimmed.fastq.gz

# (2) rRNA filter: input $fq
cd $filtered_dir
echo filtering

bowtie2 --seed=123 --threads=32 -x $rRNA_ref -U $fq -S ${sample_name}_R1.rRNA.sam \
  --un-gz=${sample_name}_R1_rFiltered.fq.gz \
  > ${sample_name}_R1_rRNA.bowtie2.log 2> ${sample_name}_R1_rRNA.bowtie2.log2

bowtie2 --seed=123 --threads=32 -x $tsoRNA_ref  \
  -U ${sample_name}_R1_rFiltered.fq.gz \
  -S ${sample_name}_R1.tsoRNA.sam \
  --un-gz=${sample_name}_R1_rtsFiltered.fq.gz \
  > ${sample_name}_R1_tsoRNA.bowtie2.log 2> ${sample_name}_R1_tsoRNA.bowtie2.log2

fq=$filtered_dir/${sample_name}_R1_rtsFiltered.fq.gz 

# samtool view and index ${sample_name}_R*.rRNA.sam
samtools view -bS ${sample_name}_R1.rRNA.sam > ${sample_name}_R1.rRNA.bam
samtools sort -@ 4 ${sample_name}_R1.rRNA.bam -o ${sample_name}_R1.rRNA.sorted.bam
mv ${sample_name}_R1.rRNA.sorted.bam ${sample_name}_R1.rRNA.bam
samtools index ${sample_name}_R1.rRNA.bam
rm ${sample_name}_R1.rRNA.sam

# fastqc
fastqc -t 4 ${sample_name}_R1_rtsFiltered.fq.gz -o $fastqc_dir

#
# (3) STAR
#
echo aligning

STAR --runThreadN 4 --outSAMattributes NH HI NM MD AS nM jM jI XS \
  --genomeDir $star_index \
  --readFilesIn $filtered_dir/${sample_name}_R1_rtsFiltered.fq.gz \
  --outFileNamePrefix $bam_dir/${sample_name}_ \
  --sjdbOverhang 25 \
  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 20 \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate

samtools index $bam_dir/${sample_name}_Aligned.sortedByCoord.out.bam