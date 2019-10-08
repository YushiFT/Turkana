#!/bin/bash
#SBATCH --job-name=2impute1k                # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=30G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

#define input
chrom=chr2
chromnum=2
in_minimac3=/Genomics/grid/users/yushit/.local/bin/Minimac3/bin
in_minimac4=/Genomics/grid/users/yushit/.local/bin/Minimac4/build
in_ref_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28_nochr.vcf
in_raw_vcf=/scratch/tmp/yushi/cohort/high_28_phased_vcf/$chrom.high_28_vqsr.vcf
in_raw_nonchr_vcf=/scratch/tmp/yushi/cohort/high_28_phased_vcf/$chrom.high_28_vqsr_nochr.vcf
in_ref_m3vcf=/Genomics/ayroleslab2/yushi/ref/public_datasets/1000GP_VCF_P3/$chromnum.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
#define output
out_prefix=/scratch/tmp/yushi/cohort/high_28_impute_vcf/$chrom.high_28_imputed

gunzip $in_raw_vcf.gz

perl -pe 's/^chr//' $in_raw_vcf > $in_raw_nonchr_vcf

$in_minimac4/minimac4 --refHaps $in_ref_m3vcf \
                      --haps $in_raw_nonchr_vcf \
                      --prefix $out_prefix \
                      --format GT,DS,GP \
                      --cpus 15
