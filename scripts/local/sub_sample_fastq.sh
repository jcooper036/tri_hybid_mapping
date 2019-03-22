#!/bin/bash

dir1=/Volumes/Jacob_2TB_storage/180905_rec_mapping_parents/fastq
dir2=/Volumes/Jacob_2TB_storage/180925_rec_samples/fastq
dir3=genomics/fastq

reads=400000 # 100,000 reads, since every 4 lines is one read

## parents
head -${reads} ${dir1}/J10075_15364X1_R1.fastq > ${dir3}/s_J10075_R1.fastq
head -${reads} ${dir1}/J10075_15364X1_R2.fastq > ${dir3}/s_J10075_R2.fastq

head -${reads} ${dir1}/J10224_15364X2_R1.fastq > ${dir3}/s_J10224_R1.fastq
head -${reads} ${dir1}/J10224_15364X2_R2.fastq > ${dir3}/s_J10224_R2.fastq

head -${reads} ${dir1}/J20039_15364X3_R1.fastq > ${dir3}/s_J20039_R1.fastq
head -${reads} ${dir1}/J20039_15364X3_R2.fastq > ${dir3}/s_J20039_R2.fastq

head -${reads} ${dir1}/J40022_15364X4_R1.fastq > ${dir3}/s_J40022_R1.fastq
head -${reads} ${dir1}/J40022_15364X4_R2.fastq > ${dir3}/s_J40022_R2.fastq

head -${reads} ${dir1}/J40023_15364X5_R1.fastq > ${dir3}/s_J40023_R1.fastq
head -${reads} ${dir1}/J40023_15364X5_R2.fastq > ${dir3}/s_J40023_R2.fastq

head -${reads} ${dir1}/J40024_15364X6_R1.fastq > ${dir3}/s_J40024_R1.fastq
head -${reads} ${dir1}/J40024_15364X6_R2.fastq > ${dir3}/s_J40024_R2.fastq

## samples
head -${reads} ${dir2}/rep1_female_15367X1_R1.fastq > ${dir3}/s_rep1_female_R1.fastq
head -${reads} ${dir2}/rep1_female_15367X1_R2.fastq > ${dir3}/s_rep1_female_R2.fastq
head -${reads} ${dir2}/rep1_male_15367X2_R1.fastq > ${dir3}/s_rep1_male_R1.fastq
head -${reads} ${dir2}/rep1_male_15367X2_R2.fastq > ${dir3}/s_rep1_male_R2.fastq

head -${reads} ${dir2}/rep2_female_15367X4_R1.fastq > ${dir3}/s_rep2_female_R1.fastq
head -${reads} ${dir2}/rep2_female_15367X4_R2.fastq > ${dir3}/s_rep2_female_R2.fastq
head -${reads} ${dir2}/rep2_male_15367X5_R1.fastq > ${dir3}/s_rep2_male_R1.fastq
head -${reads} ${dir2}/rep2_male_15367X5_R2.fastq > ${dir3}/s_rep2_male_R2.fastq

head -${reads} ${dir2}/rep3_female_15367X3_R1.fastq > ${dir3}/s_rep3_female_R1.fastq
head -${reads} ${dir2}/rep3_female_15367X3_R2.fastq > ${dir3}/s_rep3_female_R2.fastq
head -${reads} ${dir2}/rep3_male_15367X6_R1.fastq > ${dir3}/s_rep3_male_R1.fastq
head -${reads} ${dir2}/rep3_male_15367X6_R2.fastq > ${dir3}/s_rep3_male_R2.fastq