#! /bin/bash

## point of this is to coordinate the different analysis steps

# variables
# dir=genomics
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
vcf=${dir}/vcfs/parents_merge_181009.vcf
vcfTable=${dir}/vcfs/parents_merge_anyDiff_181010.tsv
bamfile=${dir}/subbams/rep3_male_15367X6_3Lsub.bam
# bamfile=${dir}/var_calling/merged_post_interval.bam

# parse the vcf - this should be a make file really
# genomics/analysis/vcf_parse.py ${vcf} ${vcfTable}

# read a bam and write an output table of the read types
genomics/analysis/read_sort.py <(samtools view -h ${bamfile}) ${vcfTable} Desktop/test_out.txt

