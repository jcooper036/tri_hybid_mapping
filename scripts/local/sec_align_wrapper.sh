#! /bin/bash
set -e

###### variables, mostly file paths
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping/

## wrapper for the make files for geneome alignment
samples=(
    J40022_15364X4
    J40023_15364X5
    J40024_15364X6
)

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
    echo $sample
    export samp=$sample
    make -f genomics/scripts/local/sec_align.mk || exit
done

