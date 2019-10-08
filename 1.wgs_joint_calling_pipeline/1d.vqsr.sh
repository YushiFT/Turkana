#!/bin/bash
#SBATCH --job-name=1vqsr                   # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Perform joint genotyping on samples pre-called with HaplotypeCaller
# define input
chrom=chr1
in_gatk=/Genomics/grid/users/yushit/.local/bin/gatk-4.1.3.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.high_28.vcf.gz
in_hapmap=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf
in_omni=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf
in_1000G=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
in_dbsnp=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
in_resourceDIR=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
# define output
out_r_plots=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.vqsr_plots.R
out_snps_recal=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.snps.recal
out_snps_tranc=/scratch/tmp/yushi/cohort/high_28_vcf/$chrom.snps.tranches
out_g_vcf_gz=/scratch/tmp/yushi/cohort/high_28_vqsr_vcf/$chrom.high_28_vqsr.vcf.gz

module load java
module load samtools

echo 'index vcfs...'

$in_gatk/gatk IndexFeatureFile -F $in_g_vcf_gz

echo 'calculating VQSLOD tranches for SNPs...'

$in_gatk/gatk --java-options "-Xmx65g -Xms65g" VariantRecalibrator \
            -V $in_g_vcf_gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $in_resourceDIR/hapmap_3.3.hg38.vcf.gz \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 $in_resourceDIR/1000G_omni2.5.hg38.vcf.gz \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 $in_resourceDIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            -O $out_snps_recal \
            --tranches-file $out_snps_tranc \
            -R $in_genome \
            --rscript-file $out_r_plots

echo 'filtering SNPs with VQSLOD...'

$in_gatk/gatk --java-options "-Xmx65g -Xms65g" ApplyVQSR \
            -V $in_g_vcf_gz \
            --recal-file $out_snps_recal \
            --tranches-file $out_snps_tranc \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O $out_g_vcf_gz \
            -R $in_genome

rm -f $out_snps_recal
rm -f $out_snps_tranc
