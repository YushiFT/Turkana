#!/bin/bash
#SBATCH --job-name=192magic_bl             # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=50G                          # memory per node
#SBATCH --time=24:00:00                    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Pick unmapped human genome sequences after bwa
# convert sorted sam files into fastq
barcode=Sample_Barcode-192
# define input
in_magic_blast=/Genomics/grid/users/yushit/.local/bin/ncbi-magicblast-1.5.0/bin
in_fastq=/scratch/tmp/yushi/$barcode.unmapped.hg38.fastq
in_ref_all=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/all/my_ref
in_ref_cholera=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/cholera/my_ref
in_ref_leishmania=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/leishmania/my_ref
in_ref_malaria=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/malaria/my_ref
in_ref_mycobacterium=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/mycobacterium/my_ref
in_ref_plasmodium=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/plasmodium/my_ref
in_ref_trichuris=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/trichuris/my_ref
in_ref_trypanosoma=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/trypanosoma/my_ref
# define output
out_all=/scratch/tmp/yushi/pathogen_blast/all/$barcode.txt
out_cholera=/scratch/tmp/yushi/pathogen_blast/cholera/$barcode.txt
out_leishmania=/scratch/tmp/yushi/pathogen_blast/leishmania/$barcode.txt
out_malaria=/scratch/tmp/yushi/pathogen_blast/malaria/$barcode.txt
out_mycobacterium=/scratch/tmp/yushi/pathogen_blast/mycobacterium/$barcode.txt
out_plasmodium=/scratch/tmp/yushi/pathogen_blast/plasmodium/$barcode.txt
out_trichuris=/scratch/tmp/yushi/pathogen_blast/trichuris/$barcode.txt
out_trypanosoma=/scratch/tmp/yushi/pathogen_blast/trypanosoma/$barcode.txt


# pick the unmapped sequences
echo 'blast against all pathogen...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_all -paired -infmt fastq -outfmt tabular -out $out_all

echo 'blast against cholera...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_cholera -paired -infmt fastq -outfmt tabular -out $out_cholera

echo 'blast against leishmania...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_leishmania -paired -infmt fastq -outfmt tabular -out $out_leishmania

echo 'blast against malaria...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_malaria -paired -infmt fastq -outfmt tabular -out $out_malaria

echo 'blast against mycobacterium...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_mycobacterium -paired -infmt fastq -outfmt tabular -out $out_mycobacterium

echo 'blast against plasmodium...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_plasmodium -paired -infmt fastq -outfmt tabular -out $out_plasmodium

echo 'blast against trichuris...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_trichuris -paired -infmt fastq -outfmt tabular -out $out_trichuris

echo 'blast against trypanosoma...'
$in_magic_blast/magicblast -query $in_fastq -db $in_ref_trypanosoma -paired -infmt fastq -outfmt tabular -out $out_trypanosoma
