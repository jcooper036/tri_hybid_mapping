SHELL := /bin/bash

.PHONY: # not sure what this is

.DELETE_ON_ERROR: # either files to delete, or if left blank means they all get deleted

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
trim_command=${java_path} -Xmx8G -jar ${trimoPath}
gatk_command=${java_path} -Xmx8G -jar ${gatk36Path}
CORES=3
mer=${dir}/var_calling/merged

all: ${mer}_Q30_SNPs.vcf ${mer}_Q30_InDel.vcf  # should be the files that need to be complete

## merge the bam files
${check}/merge_test.txt:
	samtools merge ${mer}.bam ${dir}/merge/*_sorted_rg.bam
	touch ${check}/merge_test.txt

## index the merged bam file
${check}/merge_index.txt: ${check}/merge_test.txt
	samtools index ${mer}.bam
	touch ${check}/merge_index.txt

## re-align the reads with GATK
${mer}.intervals: ${check}/merge_test.txt
	${gatk_command} \
	-T RealignerTargetCreator \
	-nt ${CORES} \
	-R ${ref}.fasta \
	-I ${mer}.bam \
	-o ${mer}.intervals

## make an alignment of the indels (long process)
${mer}_post_interval.bam: ${mer}.intervals
	${gatk_command} \
	-T IndelRealigner \
	-R ${ref}.fasta \
	-I ${mer}.bam \
	-targetIntervals ${mer}.intervals \
	-LOD 3.0 \
	-o ${mer}_post_interval.bam

## find all the raw SNPs
${mer}_raw_snp_Q30.vcf: ${mer}_post_interval.bam
	${gatk_command} \
	-T UnifiedGenotyper \
	-nt ${CORES} \
	-R ${ref}.fasta \
	-I ${mer}_post_interval.bam \
	-gt_mode DISCOVERY \
	-stand_call_conf 30 \
	-o ${mer}_raw_snp_Q30.vcf

## annotate those SNPs
${mer}_raw_snp_Q30_annotated.vcf: ${mer}_raw_snp_Q30.vcf
	${gatk_command} \
	-T VariantAnnotator \
	-nt ${CORES} \
	-R ${ref}.fasta \
	-I ${mer}_post_interval.bam \
	-G StandardAnnotation \
	-V:variant,VCF ${mer}_raw_snp_Q30.vcf \
	-XA SnpEff \
	-o ${mer}_raw_snp_Q30_annotated.vcf

## finds all raw indels (long step)
${mer}_inDels_Q30.vcf: ${mer}_raw_snp_Q30_annotated.vcf
	${gatk_command} \
	-T UnifiedGenotyper \
	-nt ${CORES} \
	-R ${ref}.fasta \
	-I ${mer}_post_interval.bam \
	-gt_mode DISCOVERY \
	-glm INDEL \
	-stand_call_conf 30 \
	-o ${mer}_inDels_Q30.vcf

## quality score filtration for SNPs
${mer}_Q30_SNPs.vcf: ${mer}_inDels_Q30.vcf
	${gatk_command} \
	-T VariantFiltration \
	-R ${ref}.fasta \
	-V ${mer}_raw_snp_Q30_annotated.vcf \
	--mask ${mer}_inDels_Q30.vcf \
	--maskExtension 5 \
	--maskName InDel \
	--clusterWindowSize 10 \
	--filterExpression "MQ0 >= 4.0 && ((MQ0 / (1.0 * DP)) > 0.1)" \
	--filterName "BadValidation" \
	--filterExpression "QUAL < 30.0" \
	--filterName "LowQual" \
	--filterExpression "QD < 5.0" \
	--filterName "LowVQCBD" \
	--filterExpression "FS > 60.0" \
	--filterName "FisherStrand" \
	-o ${mer}_Q30_SNPs.vcf

## quality filters for the indels
${mer}_Q30_InDel.vcf: ${mer}_Q30_SNPs.vcf
	${gatk_command} \
	-T VariantFiltration \
	-R ${ref}.fasta \
	-V ${mer}_Q30_SNPs.vcf \
	--clusterWindowSize 10 \
	--filterExpression "MQ0 >= 4.0 && ((MQ0 / (1.0 * DP)) > 0.1)" \
	--filterName "BadValidation" \
	--filterExpression "QUAL < 30.0" \
	--filterName "LowQual" \
	--filterExpression "QD < 5.0" \
	--filterName "LowVQCBD" \
	--filterExpression "FS > 60.0" \
	--filterName "FisherStrand" \
	-o ${mer}_Q30_InDel.vcf