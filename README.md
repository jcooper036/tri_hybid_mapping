# Tri-hybrid cross reveals new major hybrid incompatbility locus between *D. melanogaster* and *D. simulans*.  
These are the scripts and make files that I used for processing genomic data for our publication on mapping hybrid incompatibility genes. When we have a link to the paper, I will include that here.

## Requirements
* trimmomatic
* java version 1.8
* GATK 3.6
* samtools
* picard tools
* bwa

## Pipline
* trim reads
* align to the reference genome
* sort the reads
* index
* merge bam files
* call variants with GATK
* build VCF file from SNP table

## Use
The "make_wrapper.sh" in /scripts/local/ should contain all the parameters that are required for going from fastq files to a vcf, given the reference genomes. The scripts in /analysis/ are for use once there is a vcf. There is an R script for making the plots from correctly parsed SNP tables.  
  
There are also a series of helpful commands for running these programs on a SLURM managed cluster. This was run on the CHPC at the University of Utah, so your expereinece may differ. They are located in /scripts/kingspeak/