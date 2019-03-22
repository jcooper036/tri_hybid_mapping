#!/bin/bash
set -e

###### variables, mostly file paths
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
ref=${dir}/ref_genomes/mel_sim_sec_concat
scripts=/Users/Jacob/genomics/scripts/local
check=${dir}/check
CORES=3

## wrapper for the make files for geneome alignment
samples=(
    rep1_male_15367X2
    rep2_male_15367X5
    rep3_male_15367X6
    J40022_15364X4
    J20039_15364X3
    rep1_female_15367X1
    rep2_female_15367X4
    rep3_female_15367X3    
    J10075_15364X1
    J10224_15364X2
    J40023_15364X5
    J40024_15364X6
)

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
        echo $sample
        export samp=$sample
        make -f ${scripts}/trim_to_bam_local.mk || exit
done