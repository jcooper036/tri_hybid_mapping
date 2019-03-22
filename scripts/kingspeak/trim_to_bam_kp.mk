SHELL := /bin/bash

.PHONY: # not sure what this is

.DELETE_ON_ERROR: # either files to delete, or if left blank means they all get deleted

## sample from the environment
sample = ${samp}

###### variables, mostly file paths
dir=/scratch/general/lustre/u0885126/three_genomes
ref=${dir}/ref_genomes/mel_sim_sec_concat
check=genomics/check
CORES=12
bam=${dir}/bam/${sample}_sorted

all: ${check}/${sample}_index.txt # should be the files that need to be complete


## align the reads to the reference genome
${check}/${sample}_align.txt: ${dir}/${sample}_trim_pair_R1.fastq ${dir}/${sample}_trim_pair_R2.fastq
	bwa mem -t ${CORES} \
	${ref}.fasta \
	${dir}/fastq/${sample}_R1.fastq ${dir}/fastq/${sample}_R2.fastq \
	| samtools view -b > ${dir}/bam/${sample}.bam
	touch ${check}/${sample}_align.txt

## sort the reads
${check}/${sample}_sort.txt: ${check}/${sample}_align.txt
	samtools sort ${dir}/bam/${sample}.bam > ${bam}.bam
	touch ${check}/${sample}_sort.txt

## index the bam
${check}/${sample}_index.txt: ${check}/${sample}_sort.txt
	samtools index ${bam}.bam
	touch ${check}/${sample}_index.txt

