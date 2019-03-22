SHELL := /bin/bash

.PHONY: # not sure what this is

.DELETE_ON_ERROR: # either files to delete, or if left blank means they all get deleted

## sample from the environment
sample = ${samp}

###### variables, mostly file paths
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
ref=${dir}/ref_genomes/Dsec13
merged=${dir}/temp/parents_merge
gatk36Path=gatk/gatk.jar
javaVersion=1.8
trimoPath=trimmomatic/trimmomatic-0.36.jar
scripts=genomics/scripts
check=genomics/check
# java_path=/usr/libexec/java_home -v ${javaVersion} --exec java
java_path=java
trim_command=${java_path} -Xmx8G -jar ${trimoPath}
gatk_command=${java_path} -Xmx8G -jar ${gatk36Path}
CORES=3
bam=${dir}/sec_bams/${sample}_sorted

all: ${bam}_rg.bam # should be the files that need to be complete

## make an index for the reference genome
${check}/ref_index.txt:
	samtools faidx ${ref}.fasta
	touch ${check}/ref_index.txt

## bwa index
${check}/ref_bwaindex.txt: ${check}/ref_index.txt
	bwa index ${ref}.fasta
	touch ${check}/ref_bwaindex.txt

## using picard tools, make a reference dictionary from the reference fasta file
${check}/ref_dict.txt: ${check}/ref_bwaindex.txt
	picard CreateSequenceDictionary R=${ref}.fasta O=${ref}.dict
	touch ${check}/ref_dict.txt

## trim the reads target: pre-req
${check}/${sample}_trim.txt: ${check}/ref_dict.txt
	${trim_command} \
	PE -threads ${CORES} -phred33 \
	${dir}/fastq/${sample}_R1.fastq ${dir}/fastq/${sample}_R2.fastq \
	${dir}/fastq/${sample}_trim_pair_R1.fastq ${dir}/fastq/${sample}_trim_unpair_R1.fastq \
	${dir}/fastq/${t}_trim_pair_R2.fastq ${dir}/fastq/${sample}_trim_unpair_R2.fastq \
	ILLUMINACLIP:/path/to/adapter/sequecnes.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:30 # par for trimming
	touch ${check}/${sample}_trim.txt

## align the reads to the reference genome
${check}/${sample}_align.txt: ${check}/${sample}_trim.txt
	bwa mem -t 3 \
	${ref}.fasta \
	${dir}/fastq/${sample}_R1.fastq ${dir}/fastq/${sample}_R2.fastq \
	| samtools view -b > ${dir}/sec_bams/${sample}.bam
	touch ${check}/${sample}_align.txt

## sort the reads
${check}/${sample}_sort.txt: ${check}/${sample}_align.txt
	samtools sort ${dir}/sec_bams/${sample}.bam > ${bam}.bam
	touch ${check}/${sample}_sort.txt

## index the bam
${check}/${sample}_index.txt: ${check}/${sample}_sort.txt
	samtools index ${bam}.bam
	touch ${check}/${sample}_index.txt

## add read groups
${bam}_rg.bam: ${check}/${sample}_index.txt
	picard AddOrReplaceReadGroups \
	I=${bam}.bam \
	O=${bam}_rg.bam \
	SORT_ORDER=coordinate \
	RGPL=illumina \
	RGPU=unit1 \
	RGLB=lib1 \
	RGID=${sample} \
	RGSM=${sample}
