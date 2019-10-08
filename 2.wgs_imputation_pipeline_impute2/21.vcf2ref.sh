#!/bin/bash
#SBATCH --job-name=vcftoref                # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=70G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu


# command line variables
#my $vcf_file = $arg{-vcf};                    # VCF file to process; can be phased or unphased
#my $gen_file = $arg{-gen};                    # name of new file to print; suffix and gzipping are auto
#my $chrom = $arg{-chr};                       # chromosome to include in output files
#my $start_pos = $arg{-start};                 # first position to include in output files
#my $end_pos = $arg{-end};                     # last position to include in output files
#my $samp_include_file = $arg{-samp_include};  # file containing a list of samples to include in gen file
#my $samp_exclude_file = $arg{-samp_exclude};  # file containing a list of samples to exclude from gen file
#my $is_samp_snptest = $arg{-samp_snptest};    # boolean: should we print the sample file in SNPTEST format?
#my $is_special_gt = $arg{-special_gt};        # boolean: allow 'special' genotype coding in VCF?
#my $is_gt_like = $arg{-gt_like};              # boolean: print genotype likelihoods rather than hard calls?
#my $no_mono = $arg{-no_mono};                 # boolean: omit monomorphic sites from processed files?
#my $snps_only = $arg{-snps_only};             # boolean: omit any sites that are not biallelic SNPs?
#my $indels_only = $arg{-indels_only};         # boolean: omit any sites that are not biallelic INDELs?
#my $svs_only = $arg{-svs_only};               # boolean: omit any sites that are not biallelic SVs?



perl ./20.vcf2impute_gen.pl \
 -vcf /scratch/tmp/yushi/cohort/high_28_vcf/chr20.high_28.vcf \
 -gen ./test.gen \
 -chr chr20 \
 -start 60000 \
 -end 64334159
