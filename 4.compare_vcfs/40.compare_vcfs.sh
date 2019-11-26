#!/bin/bash
#SBATCH --job-name=vcftoref                # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=30G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

#define input
in_bcftools=/Genomics/grid/users/yushit/.local/bin/bcftools
in_test_vcf_gz=/scratch/tmp/yushi/cohort/high_28_impute_vcf/chr10.high_28_imputed.dose.vcf.gz
in_test_vcf=/scratch/tmp/yushi/cohort/high_28_impute_vcf/chr10.high_28_imputed.dose.vcf
in_array_vcf=/scratch/tmp/yushi/liftover_array/array_3batches_hg38_nonchr.vcf
#define output
out_insec=/scratch/tmp/yushi/compare_vcf/insec_high_28_array

gunzip $in_test_vcf_gz

$in_bcftools/bcftools isec -p $out_insec -Oz $in_test_vcf_gz $in_array_vcf
