#!/bin/bash
#SBATCH --job-name=497pathobwa             # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=8                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=66G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Pick unmapped human genome sequences after bwa
# convert sorted sam files into fastq
barcode=Sample_Barcode-497
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R2.fastq.gz
# define output
out_sam=/scratch/tmp/yushi/$barcode.hg38.sam
out_bam=/scratch/tmp/yushi/$barcode.hg38.bam
out_bam_sorted=/scratch/tmp/yushi/$barcode.hg38.sort.bam
out_fastq=/scratch/tmp/yushi/$barcode.unmapped.hg38.fastq

# align to the human genome
echo 'mapping bwa...'
originalpath=$PATH
export PATH=/Genomics/grid/users/yushit/.local/bin/bwa-0.7.17/:$PATH
bwa mem -t 8 /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa $in_fastq1 $in_fastq2 > $out_sam
PATH=$originalpath

# pick the unmapped sequences
echo 'picking unmapped reads...'
module load samtools
samtools view -U $out_bam -o /dev/null $out_sam

echo 'converting to fastq...'
samtools sort -T /scratch/tmp/yushi/temp -m 36G $out_bam | samtools bam2fq > $out_fastq
