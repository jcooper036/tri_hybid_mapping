#!/bin/bash
#SBATCH -t 2:00:00                  # max time job will take
#SBATCH -n 64                       # number of cores (usally 64 max per node)
#SBATCH -N 1                        # number of nodes
#SBATCH -A owner-guest              # account
#SBATCH -p kingspeak-guest          # partition (guest is free) other option: ember won't boot, takes forever to run
#SBATCH --mail-user=jcooper036@gmail.com    # 
#SBATCH --requeue       ## put this back in the queue if it does not finish

##### load modules
module load gatk
module load samtools

set -e

###### variables, mostly file paths
dir=/scratch/general/lustre/u0885126/rec_map
ref=${dir}/ref_genomes/Dmel_r6_180905
gatk36Path=gatk/gatk.jar
javaVersion=1.8
trimoPath=trimmomatic/trimmomatic-0.36.jar
scripts=genomics/scripts
check=${dir}/check
mer=${dir}/var_calling/merged
CORES=64
gatk_command=gatk

## wrapper for the make files for geneome alignment
samples=(
    s_J10075
    s_J10224
    s_J20039
    s_J40022
    s_J40023
    s_J40024
    s_rep1_female
    s_rep1_male
    s_rep2_female
    s_rep2_male
    s_rep3_female
    s_rep3_male
)

## merge the bam files
make -f genomics/scripts/merge_to_vcf.mk

## get only the SNPs that passed
cat ${mer}_Q30_SNPs.vcf | grep 'PASS\|^#' > ${mer}_SNPs_cleaned.vcf
cat ${mer}_Q30_InDel.vcf | grep 'PASS\|^#' > ${mer}_InDels_cleaned.vcf

## make a table of the SNPs
${gatk_command} \
    -R ${ref}.fasta \
    -T VariantsToTable \
    -V ${mer}_SNPs_cleaned.vcf \
    -F CHROM -F POS -F ID -F QUAL -F AC \
    -o ${mer}_SNPs_cleaned.table
