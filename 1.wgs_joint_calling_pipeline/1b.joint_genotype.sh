#!/bin/bash
#SBATCH --job-name=1jointcall              # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=8                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Perform joint genotyping on samples pre-called with HaplotypeCaller
# define input
chrom=chr1
in_gatk=/Genomics/grid/users/yushit/.local/bin/gatk-4.1.3.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_database=/scratch/tmp/yushi/cohort/high_28/$chrom
# define output
out_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.vcf.gz

module load java
module load samtools


echo 'calling joint genotype...' # sort bam file

$in_gatk/gatk --java-options "-Xmx65g" GenotypeGVCFs \
            -R $in_genome \
            -V gendb://$in_database \
            -O $out_g_vcf_gz \
            --tmp-dir=/scratch/tmp/yushi
