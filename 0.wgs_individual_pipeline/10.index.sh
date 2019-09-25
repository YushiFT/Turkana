#!/bin/bash
#SBATCH --job-name=vcfindex      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                # memory per node
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send mail when process begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

in_variants=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
in_gatk=/Genomics/grid/users/yushit/.local/bin/gatk-4.1.3.0

out_variants=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

module load java

$in_gatk/gatk IndexFeatureFile -F $in_variants 
