#!/bin/bash
#SBATCH --job-name=RefDatBLAST             # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=4                  # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=66G                          # memory per node
#SBATCH --time=24:00:00                    # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Pick unmapped human genome sequences after bwa
# convert sorted sam files into fastq
barcode=Sample_Barcode-192
# define input
in_magic_blast=/Genomics/grid/users/yushit/.local/bin/ncbi-magicblast-1.5.0/bin
in_all=/Genomics/ayroleslab2/yushi/data/pathogen_ref/pathogen_combined_ref_ncbi.fasta
in_cholera=/Genomics/ayroleslab2/yushi/data/pathogen_ref/cholera_pasteurella_multocida_genomic.fna
in_leishmania=/Genomics/ayroleslab2/yushi/data/pathogen_ref/leishmania_donovani_genomic.fna
in_malaria=/Genomics/ayroleslab2/yushi/data/pathogen_ref/malaria_vivax_genomic.fna
in_mycobacterium=/Genomics/ayroleslab2/yushi/data/pathogen_ref/mycobacterium_tuberculosis.fna
in_plasmodium=/Genomics/ayroleslab2/yushi/data/pathogen_ref/plasmodium_genomic.fna
in_trichuris=/Genomics/ayroleslab2/yushi/data/pathogen_ref/trichuris_trichiura_genomic.fna
in_trypanosoma=/Genomics/ayroleslab2/yushi/data/pathogen_ref/trypanosoma_brucei.fna
# define output
out_all=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/all/my_ref
out_cholera=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/cholera/my_ref
out_leishmania=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/leishmania/my_ref
out_malaria=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/malaria/my_ref
out_mycobacterium=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/mycobacterium/my_ref
out_plasmodium=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/plasmodium/my_ref
out_trichuris=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/trichuris/my_ref
out_trypanosoma=/Genomics/ayroleslab2/yushi/ref/pathogen_ref/trypanosoma/my_ref


echo 'construct reference database for all...'
$in_magic_blast/makeblastdb -in $in_all -out $out_all -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_cholera -out $out_cholera -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_leishmania -out $out_leishmania -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_malaria -out $out_malaria -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_mycobacterium -out $out_mycobacterium -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_plasmodium -out $out_plasmodium -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_trichuris -out $out_trichuris -parse_seqids -dbtype nucl
$in_magic_blast/makeblastdb -in $in_trypanosoma -out $out_trypanosoma -parse_seqids -dbtype nucl
