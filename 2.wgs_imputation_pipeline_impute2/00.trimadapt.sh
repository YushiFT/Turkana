#!/bin/bash
#SBATCH --job-name=refpanel           # create a short name for your job
#SBATCH --nodes=1                     # node count
#SBATCH --ntasks-per-node=1           # total number of tasks across all nodes
#SBATCH --cpus-per-task=8             # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                     # memory per node
#SBATCH --time=24:00:00 --quos=1wk    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin             # send mail when process begins
#SBATCH --mail-type=end               # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Filter 1000G reference panel with only AFRs with IMPUTE2

# define file head
barcode=Sample_Barcode-192
# define input
in_impute2=//Genomics/grid/users/yushit/.local/bin/impute_v2.3.2_x86_64_static
# define output


$in_impute2/impute2 \
 -m ./Example/example.chr22.map \
 -h ./Example/example.chr22.1kG.haps \
 -l ./Example/example.chr22.1kG.legend \
 -g ./Example/example.chr22.study.gens \
 -strand_g ./Example/example.chr22.study.strand \
 -int 1 50000 \
 -int 100000 150000 \
 -Ne 20000 \
 -o ./Example/example.chr22.one.phased.impute2
