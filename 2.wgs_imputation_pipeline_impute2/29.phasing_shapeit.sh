#!/bin/bash
#SBATCH --job-name=20shapeit               # create a short name for your job
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
in_shapeit=/Genomics/grid/users/yushit/.local/bin/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin
in_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.vcf
in_map=/Genomics/ayroleslab2/yushi/ref/public_datasets/1000GP_Phase3/genetic_map_chr20_combined_b37.txt
# define output
out_gen_phased=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.phased

$in_shapeit/shapeit --input-vcf $in_vcf \
            -M $in_map \
            -O $out_gen_phased \
            -T 16
