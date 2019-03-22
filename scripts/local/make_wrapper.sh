#! /bin/bash
set -e

###### variables, mostly file paths
dir=genomics
ref=${dir}/ref_genomes/Dmel_r6_180905
merged=${dir}/temp/parents_merge
gatk36Path=gatk/gatk.jar
javaVersion=1.8
trimoPath=trimmomatic/trimmomatic-0.36.jar
scripts=genomics/scripts
check=genomics/check
# java_path=/usr/libexec/java_home -v ${javaVersion} --exec java
java_path=java
trim_command="${java_path} -Xmx8G -jar ${trimoPath}"
gatk_command="${java_path} -Xmx8G -jar ${gatk36Path}"
CORES=3
mer=${dir}/var_calling/merged

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

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
    echo $sample
    export samp=$sample
    make -f genomics/scripts/fastq_to_read_groups.mk || exit
done

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
