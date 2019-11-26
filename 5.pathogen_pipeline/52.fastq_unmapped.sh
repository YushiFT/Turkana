#!/bin/bash
#SBATCH --job-name=312unmap_pa             # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=3                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=66G                          # memory per node
#SBATCH --time=24:00:00                    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Pick unmapped human genome sequences after bwa
# convert sorted sam files into fastq
barcode=Sample_Barcode-312
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Project_AYR_13560_B01_NAN_Lane.2018-09-28/$barcode/fastq/*R2.fastq.gz
in_magic_blast=/Genomics/grid/users/yushit/.local/bin/ncbi-magicblast-1.5.0/bin
in_ref_all=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/all/my_ref
# define output
out_sam=/scratch/tmp/yushi/$barcode.hg38.sam
out_bam=/scratch/tmp/yushi/$barcode.hg38.bam
out_bam_sorted=/scratch/tmp/yushi/$barcode.hg38.sort.bam
out_fastq=/scratch/tmp/yushi/$barcode.unmapped.hg38.fastq
out_all_table=/scratch/tmp/yushi/pathogen_blast/all/$barcode.txt

# pick the unmapped sequences
echo 'picking unmapped reads...'
module load samtools
samtools view -S -b -f 4 $out_sam > $out_bam

echo 'sorting bam...'
samtools sort -T /scratch/tmp/yushi/temp -m 36G -o $out_bam_sorted $out_bam

echo 'converting to fastq...'
samtools bam2fq $out_bam_sorted > $out_fastq

rm -f $out_bam
rm -f $out_bam_sorted

# pick the unmapped sequences
echo 'blast against all pathogen...'
$in_magic_blast/magicblast -query $out_fastq -db $in_ref_all -paired -infmt fastq -outfmt tabular -out $out_all_table
