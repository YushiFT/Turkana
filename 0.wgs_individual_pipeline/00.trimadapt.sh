#!/bin/bash
#SBATCH --job-name=bwamapping      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=1000                # memory per node
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send mail when process begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=yushifelixtang@gmail.com

# Remove adapter sequences, primers, poly-A tails and other unwanted sequence
# from WGS data.
# High coverage ones includes barcode as
# 368, 345, 192, 169, 281, 381, 294, 234, 323, 165, 338, 334, 370,
# 308, 332, 211, 258, 203, 218, 249, 245, 242, 202, 304, 235, 309,
# 330, 205
# define file head
barcode=Sample_Barcode-192
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R2.fastq.gz
# define output
out_fastq1=/scratch/tmp/yushi/$barcode.trim.R1.fastq.gz
out_fastq2=/scratch/tmp/yushi/$barcode.trim.R2.fastq.gz

originalpath=$PATH
export PATH=/Genomics/grid/users/yushit/.local/bin/:$PATH
cutadapt --nextseq-trim 20 -e 0.05 --overlap 2 --minimum-lengt=20 --trim-n -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $out_fastq1 -p $out_fastq2 $in_fastq1 $in_fastq2
PATH=$originalpath
