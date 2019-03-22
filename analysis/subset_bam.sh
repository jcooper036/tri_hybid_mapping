#!/bin/bash
set -e

## set variables, make the D.sec consensus
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
ref=${dir}/ref_genomes/Dmel_r6_180905.fasta
bam=${dir}/subbams/dsec_3Lsub.bam
out=genomics/analysis/consensus_building/sec_3Lsub.fastq
chr=2R
min=17100000 #Lhr = 17429032
max=17700000

## sim sample
sample=J20039_15364X3
samtools index ${dir}/merge/${sample}_sorted_rg.bam
samtools view -b ${dir}/merge/${sample}_sorted_rg.bam "${chr}:${min}-${max}" > ${dir}/subbams/${sample}_sim_2Rsub.bam

## sec samples
sample=J40022_15364X4
samtools index ${dir}/merge/${sample}_sorted_rg.bam
samtools view -b ${dir}/merge/${sample}_sorted_rg.bam "${chr}:${min}-${max}" > ${dir}/subbams/${sample}_sec_2Rsub.bam

sample=J40023_15364X5
samtools index ${dir}/merge/${sample}_sorted_rg.bam
samtools view -b ${dir}/merge/${sample}_sorted_rg.bam "${chr}:${min}-${max}" > ${dir}/subbams/${sample}_sec_2Rsub.bam

sample=J40024_15364X6
samtools index ${dir}/merge/${sample}_sorted_rg.bam
samtools view -b ${dir}/merge/${sample}_sorted_rg.bam "${chr}:${min}-${max}" > ${dir}/subbams/${sample}_sec_2Rsub.bam

samples=(
    s_J40022
    s_J40023
    s_J40024
)

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
    echo $sample
    export samp=$sample
    make -f genomics/scripts/local/sec_align.mk || exit
done