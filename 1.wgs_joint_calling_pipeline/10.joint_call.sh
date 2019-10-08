#!/bin/bash
#SBATCH --job-name=20joincal               # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                          # memory per node
#SBATCH --time=165:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Perform joint genotyping on samples pre-called with HaplotypeCaller
# define input
chrom=chr20
in_gatk=/Genomics/grid/users/yushit/.local/bin/gatk-4.1.3.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_database=/scratch/tmp/yushi/cohort/high_28/$chrom

in_variant1=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-165.hg38.g.vcf
in_variant2=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-169.hg38.g.vcf
in_variant3=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-192.hg38.g.vcf
in_variant4=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-202.hg38.g.vcf
in_variant5=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-203.hg38.g.vcf
in_variant6=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-205.hg38.g.vcf
in_variant7=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-211.hg38.g.vcf
in_variant8=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-218.hg38.g.vcf
in_variant9=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-234.hg38.g.vcf
in_variant10=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-235.hg38.g.vcf
in_variant11=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-242.hg38.g.vcf
in_variant12=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-245.hg38.g.vcf
in_variant13=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-249.hg38.g.vcf
in_variant14=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-258.hg38.g.vcf
in_variant15=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-281.hg38.g.vcf
in_variant16=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-294.hg38.g.vcf
in_variant17=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-304.hg38.g.vcf
in_variant18=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-308.hg38.g.vcf
in_variant19=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-309.hg38.g.vcf
in_variant20=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-323.hg38.g.vcf
in_variant21=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-330.hg38.g.vcf
in_variant22=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-332.hg38.g.vcf
in_variant23=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-334.hg38.g.vcf
in_variant24=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-338.hg38.g.vcf
in_variant25=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-345.hg38.g.vcf
in_variant26=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-368.hg38.g.vcf
in_variant27=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-370.hg38.g.vcf
in_variant28=/Genomics/ayroleslab2/yushi/data/turkana_wgs_vcf/Sample_Barcode-381.hg38.g.vcf

# define output
out_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.vcf.gz


module load java
module load samtools


echo 'merge all vcfs...' # sort bam file

$in_gatk/gatk --java-options "-Xmx65g -Xms65g" GenomicsDBImport \
            -V $in_variant1 \
            -V $in_variant2 \
            -V $in_variant3 \
            -V $in_variant4 \
            -V $in_variant5 \
            -V $in_variant6 \
            -V $in_variant7 \
            -V $in_variant8 \
            -V $in_variant9 \
            -V $in_variant10 \
            -V $in_variant11 \
            -V $in_variant12 \
            -V $in_variant13 \
            -V $in_variant14 \
            -V $in_variant15 \
            -V $in_variant16 \
            -V $in_variant17 \
            -V $in_variant18 \
            -V $in_variant19 \
            -V $in_variant20 \
            -V $in_variant21 \
            -V $in_variant22 \
            -V $in_variant23 \
            -V $in_variant24 \
            -V $in_variant25 \
            -V $in_variant26 \
            -V $in_variant27 \
            -V $in_variant28 \
            --genomicsdb-workspace-path /scratch/tmp/yushi/cohort/high_28/$chrom \
            --tmp-dir=/scratch/tmp/yushi \
            -L $chrom \
            --reader-threads 16

echo 'calling joint genotype...'


$in_gatk/gatk --java-options "-Xmx65g" GenotypeGVCFs \
            -R $in_genome \
            -V gendb://$in_database \
            -O $out_g_vcf_gz \
            --tmp-dir=/scratch/tmp/yushi
