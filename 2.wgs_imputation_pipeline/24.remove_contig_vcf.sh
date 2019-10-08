#!/bin/bash
#SBATCH --job-name=20mvconti               # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=40G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu


# using shapeit to phase genotype data

# define input
chrom=chr20
in_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.vcf
# define output
out_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28_nochr.vcf

perl -pe 's/^chr//' $in_vcf > $out_vcf
