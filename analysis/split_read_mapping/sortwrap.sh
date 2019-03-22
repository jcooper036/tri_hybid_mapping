## point of this is to coordinate the different analysis steps

# variables
dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
threeGenomes=${dir}/subbams/rep1_male_threespec_3L.bam
crossRef=${dir}/subbams/rep1_male_melOnly_3Lsub.bam
vcfTable=${dir}/vcfs/parents_merge_anyDiff_181010.tsv
outFile=Desktop/test.txt

# read a bam and write an output table of the read types
genomics/analysis/split_read_mapping/three_genome_sort.py <(samtools view -h ${threeGenomes}) <(samtools view -h ${threeGenomes}) <(samtools view -h ${crossRef}) ${vcfTable} ${outFile}


# genomics/analysis/split_read_mapping/three_genome_sort.py < (samtools view -h ${threeGenomes}) (samtools view -h ${crossRef}) Desktop/test_out.txt
