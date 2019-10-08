#!/bin/bash
#SBATCH --job-name=1phasing               # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=30G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

#define input
chrom=chr1
chromnum=1
in_bcftools=/Genomics/grid/users/yushit/.local/bin/bcftools
in_eagle=/Genomics/grid/users/yushit/.local/bin/Eagle_v2.4.1
in_ref_vcf=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
in_raw_vcf=/scratch/tmp/yushi/cohort/high_28_vqsr_vcf/$chrom.high_28_vqsr.vcf.gz
in_unzip_vcf=/scratch/tmp/yushi/cohort/high_28_vqsr_vcf/$chrom.high_28_vqsr.vcf
in_ref_m3vcf=/Genomics/ayroleslab2/yushi/ref/public_datasets/1000GP_VCF_P3/$chromnum.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
#define output
out_split_vcf=/scratch/tmp/yushi/cohort/high_28_vqsr_vcf/$chrom.high_28_split.vcf
out_phased_vcf=/scratch/tmp/yushi/cohort/high_28_phased_vcf/$chrom.high_28_vqsr



echo 'spliting multi-alleles...'

module load samtools

#gunzip $in_raw_vcf

$in_bcftools/bcftools norm -m - $in_unzip_vcf > $out_split_vcf

echo 'phasing...'

$in_eagle/eagle --vcf $out_split_vcf \
                --geneticMapFile $in_eagle/tables/genetic_map_hg38_withX.txt.gz \
                --outPrefix $out_phased_vcf
