#!/bin/bash
#SBATCH --job-name=22tri_bwa              # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=4                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                          # memory per node
#SBATCH --time=24:00:00                    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Trimming adapter and bwa mapping
# High coverage ones includes barcode as
# 192, 169, 368, 345, 281, 381, 294, 234, 323, 165, 338, 334, 370,
# 308, 332, 211, 258, 203, 218, 249, 245, 242, 202, 304, 235, 309,
# 330, 205
# define file head
barcode=Sample_Barcode-22
#fastcode1=981-read-1
#fastcode2=981-read-4
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13970_B01_NAN_Lane.2019-04-19/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13970_B01_NAN_Lane.2019-04-19/$barcode/fastq/*R2.fastq.gz
in_picard=/Genomics/grid/users/yushit/.local/bin/picard.jar
in_gatk=/Genomics/grid/users/yushit/.local/bin/gatk-4.1.3.0
in_sambamba=/Genomics/grid/users/yushit/.local/bin/sambamba-0.7.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_variants=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
# define output
out_fastq1=/scratch/tmp/yushi/$barcode.trim.R1.fastq.gz
out_fastq2=/scratch/tmp/yushi/$barcode.trim.R2.fastq.gz
out_sam=/scratch/tmp/yushi/$barcode.hg38.sam
out_bam=/scratch/tmp/yushi/$barcode.hg38.bam
out_bam_sorted=/scratch/tmp/yushi/$barcode.hg38.sorted.bam
out_bam_duplicates=/scratch/tmp/yushi/$barcode.hg38.dup.bam
out_txt_dupmetrics=/scratch/tmp/yushi/$barcode.hg38.dup.txt
out_bam_readgroups=/scratch/tmp/yushi/$barcode.hg38.reg.bam
out_table_recalibr=/scratch/tmp/yushi/$barcode.rbqs.table
out_bam_recalibrat=/scratch/tmp/yushi/$barcode.rbqs.bam
out_g_vcf=/scratch/tmp/yushi/turkana_wgs_vcf/$barcode.hg38.g.vcf

module load java
module load samtools

echo 'trimming...'
originalpath=$PATH
export PATH=/Genomics/grid/users/yushit/.local/bin/:$PATH
cutadapt --nextseq-trim 20 -e 0.05 --overlap 2 --minimum-lengt=20 --trim-n -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $out_fastq1 -p $out_fastq2 $in_fastq1 $in_fastq2
PATH=$originalpath


echo 'mapping...'
originalpath=$PATH
export PATH=/Genomics/grid/users/yushit/.local/bin/bwa-0.7.17/:$PATH
bwa mem -t 4 /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa $out_fastq1 $out_fastq2 > $out_sam
PATH=$originalpath

echo 'bam counting...'
samtools view -Sbq 1 $out_sam > $out_bam
echo $barcode
samtools view $out_bam | wc -l

rm -f $out_sam
rm -f $out_fastq1
rm -f $out_fastq2
