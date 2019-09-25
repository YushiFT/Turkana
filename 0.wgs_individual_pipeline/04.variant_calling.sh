#!/bin/bash
#SBATCH --job-name=169vcfcal     # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=66G                # memory per node
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send mail when process begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Calling GATK variants and export VCF
# High coverage ones includes barcode as
# 192, 169, 368, 345, 281, 381, 294, 234, 323, 165, 338, 334, 370,
# 308, 332, 211, 258, 203, 218, 249, 245, 242, 202, 304, 235, 309,
# 330, 205
# define file head
barcode=Sample_Barcode-169
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R2.fastq.gz
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
out_g_vcf=/scratch/tmp/yushi/$barcode.hg38.g.vcf

module load java

echo 'calling GATK variants as VCF and zip...' # sort bam file

$in_gatk/gatk --java-options "-Xmx60g" HaplotypeCaller -R $in_genome -I $out_bam_recalibrat -O $out_g_vcf -ERC GVCF --max-alternate-alleles 2
