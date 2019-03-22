#! /bin/bash

dir=/Volumes/Jacob_2TB_storage/rec_resc_mapping
out=${dir}/subbams/dsec_2Rsub

samtools merge ${out}.bam ${dir}/subbams/*_sec_2Rsub.bam