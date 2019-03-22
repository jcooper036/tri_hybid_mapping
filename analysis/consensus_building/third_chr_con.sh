#! /bin/bash


## set variables, make the D.sec consensus
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
ref=${dir}/ref_genomes/Dmel_r6_180905.fasta
chr=2R
bam=${dir}/subbams/dsec_${chr}sub.bam
out=genomics/analysis/consensus_building/sec_${chr}sub.fastq

samtools mpileup -uf ${ref} ${bam} | bcftools call -c | vcfutils.pl vcf2fq >${out}

## change some variables, make the D.sim consensus
bam=${dir}/subbams/J20039_15364X3_${chr}sub.bam
out=genomics/analysis/consensus_building/sim_${chr}sub.fastq

samtools mpileup -uf ${ref} ${bam} | bcftools call -c | vcfutils.pl vcf2fq >${out}