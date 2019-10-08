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
chrom=chr20
in_minimac3=/Genomics/grid/users/yushit/.local/bin/Minimac3/bin
in_ref_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28_nochr.vcf
#define output

$in_minimac3/Minimac3 --refHaps $in_ref_vcf \
                      --processReference \
                      --prefix refPanel \
                      --allTypedSites
