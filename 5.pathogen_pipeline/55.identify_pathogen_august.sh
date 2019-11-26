#!/bin/bash
#SBATCH --job-name=1106pathogen             # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=4                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=60G                          # memory per node
#SBATCH --time=24:00:00                    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Pick unmapped human genome sequences after bwa
# convert sorted sam files into fastq
barcode=1106
fastcode1=1106-read-1
fastcode2=1106-read-4
# define input
in_fastq1=/Genomics/ayroleslab2/alea/archive_raw_fastq/Turkana_NovaSeq_25Aug19/$fastcode1.fastq.gz
in_fastq2=/Genomics/ayroleslab2/alea/archive_raw_fastq/Turkana_NovaSeq_25Aug19/$fastcode2.fastq.gz
in_magic_blast=/Genomics/grid/users/yushit/.local/bin/ncbi-magicblast-1.5.0/bin
in_ref_all=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/all/my_ref
# define output
out_sam=/scratch/tmp/yushi/$barcode.hg38.sam
out_bam=/scratch/tmp/yushi/$barcode.hg38.bam
out_bam_sorted=/scratch/tmp/yushi/$barcode.hg38.sort.bam
out_fastq=/scratch/tmp/yushi/$barcode.unmapped.hg38.fastq
out_all_table=/scratch/tmp/yushi/pathogen_blast/all/$barcode.txt

# align to the human genome
echo 'mapping bwa...'
originalpath=$PATH
export PATH=/Genomics/grid/users/yushit/.local/bin/bwa-0.7.17/:$PATH
bwa mem -t 4 /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa $in_fastq1 $in_fastq2 > $out_sam
PATH=$originalpath

# pick the unmapped sequences
echo 'picking unmapped reads...'
module load samtools
samtools view -S -b -f 4 $out_sam > $out_bam

rm -f $out_sam

echo 'sorting bam...'
samtools sort -T /scratch/tmp/yushi/temp -m 36G -o $out_bam_sorted $out_bam

echo 'converting to fastq...'
samtools bam2fq $out_bam_sorted > $out_fastq

rm -f $out_bam
rm -f $out_bam_sorted

# pick the unmapped sequences
echo 'blast against all pathogen...'
$in_magic_blast/magicblast -query $out_fastq -db $in_ref_all -paired -infmt fastq -outfmt tabular -out $out_all_table

rm -f $out_fastq
