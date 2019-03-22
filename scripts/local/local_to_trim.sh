#! /bin/bash
set -e

###### variables, mostly file paths
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
ref=${dir}/ref_genomes/mel_sim_sec_concat
merged=${dir}/temp/parents_merge
gatk36Path=gatk/gatk.jar
javaVersion=1.8
trimoPath=trimmomatic/trimmomatic-0.36.jar
scripts=${dir}/scripts
check=${dir}/check
# java_path=/usr/libexec/java_home -v ${javaVersion} --exec java
java_path=java
trim_command="${java_path} -Xmx8G -jar ${trimoPath}"
gatk_command="${java_path} -Xmx8G -jar ${gatk36Path}"
CORES=3
mer=${dir}/var_calling/merged

## wrapper for the make files for geneome alignment
samples=(
    J10075_15364X1
    J10224_15364X2
    J20039_15364X3
    J40022_15364X4
    J40023_15364X5
    J40024_15364X6
    rep1_female_15367X1
    rep1_male_15367X2
    rep2_female_15367X4
    rep2_male_15367X5
    rep3_female_15367X3
    rep3_male_15367X6
)

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
    echo $sample
    export samp=$sample
    make -f genomics/scripts/local/fastq_to_trim_groups.mk || exit
done