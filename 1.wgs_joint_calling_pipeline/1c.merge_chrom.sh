#!/bin/bash
#SBATCH --job-name=vqsr                    # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=7                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=50G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Perform joint genotyping on samples pre-called with HaplotypeCaller
# define input
in_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vcf/chr*.high_28.vcf.gz
# define output
out_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vcf/all_chr.high_28.vcf.gz

echo 'merging variants in all chromosomes...'

bcftools concat -o high.cov.allchr.vcf.gz $in_g_vcf_gz
