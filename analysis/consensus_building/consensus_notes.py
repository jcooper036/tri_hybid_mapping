#!/usr/bin/env python3

samtools mpileup -uf chr5B_leaf_rust.fasta M27454_5B_reads_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq
