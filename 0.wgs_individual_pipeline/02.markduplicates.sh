#!/bin/bash
#SBATCH --job-name=192picard      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                # memory per node
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send mail when process begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Identify duplicated reads with sorted bam or sam files
# High coverage ones includes barcode as
# 368, 345, 192, 169, 281, 381, 294, 234, 323, 165, 338, 334, 370,
# 308, 332, 211, 258, 203, 218, 249, 245, 242, 202, 304, 235, 309,
# 330, 205
# define file head
barcode=Sample_Barcode-192
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R2.fastq.gz
in_picard=/Genomics/grid/users/yushit/.local/bin/picard.jar
# define output
out_fastq1=/scratch/tmp/yushi/$barcode.trim.R1.fastq.gz
out_fastq2=/scratch/tmp/yushi/$barcode.trim.R2.fastq.gz
out_sam=/scratch/tmp/yushi/$barcode.hg38.sam
out_bam=/scratch/tmp/yushi/$barcode.hg38.bam
out_bam_sorted=/scratch/tmp/yushi/$barcode.hg38.sorted.bam
out_bam_duplicates=/scratch/tmp/yushi/$barcode.hg38.dup.bam
out_txt_dupmetrics=/scratch/tmp/yushi/$barcode.hg38.dup.txt

echo 'samtools sorting...' # sort bam file

module load java
module load samtools

samtools sort -m 3G -o $out_bam_sorted $out_bam

echo 'picard duplications...' # run picard to mark duplicates

java -Xmx60g -jar $in_picard MarkDuplicates I=$out_bam_sorted O=$out_bam_duplicates M=$out_txt_dupmetrics
