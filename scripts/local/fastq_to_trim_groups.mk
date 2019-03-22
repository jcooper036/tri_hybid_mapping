SHELL := /bin/bash

.PHONY: # not sure what this is

.DELETE_ON_ERROR: # either files to delete, or if left blank means they all get deleted

## sample from the environment
sample = ${samp}

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
trim_command=${java_path} -Xmx8G -jar ${trimoPath}
gatk_command=${java_path} -Xmx8G -jar ${gatk36Path}
CORES=3
bam=${dir}/bam/${sample}_sorted

all: ${check}/${sample}_trim.txt # should be the files that need to be complete

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
	${dir}/fastq/${sample}_trim_pair_R2.fastq ${dir}/fastq/${sample}_trim_unpair_R2.fastq \
	ILLUMINACLIP:/path/to/adapter/sequecnes.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:30 # par for trimming
	touch ${check}/${sample}_trim.txt