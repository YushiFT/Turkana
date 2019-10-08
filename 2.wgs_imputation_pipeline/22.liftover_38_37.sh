#!/bin/bash
#SBATCH --job-name=20liftove               # create a short name for your job
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
in_picard=/Genomics/grid/users/yushit/.local/bin/picard.jar
in_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28_nochr.vcf
in_chain=/Genomics/ayroleslab2/yushi/ref/public_datasets/liftover/hg38ToHg19.over.chain
in_genome=/Genomics/ayroleslab2/yushi/ref/hg37_10kb/hg19_v0_Homo_sapiens_assembly19.fasta
# define output
out_reject=/scratch/tmp/yushi/cohort/high_28_vcf/liftover/$chrom.high_28_38to37_reject.vcf
out_lifted_vcf=/scratch/tmp/yushi/cohort/high_28_vcf/liftover/$chrom.high_28_hg37_lifted.vcf

module load java

java -jar $in_picard LiftoverVcf \
   I=$in_vcf \
   O=$out_lifted_vcf \
   CHAIN=$in_chain \
   REJECT=$out_reject \
   R=$in_genome \
   WARN_ON_MISSING_CONTIG=true
