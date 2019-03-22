#!/bin/bash
#SBATCH -t 60:00:00                 # max time job will take
#SBATCH -n 12                       # number of cores (usally 64 max per node)
#SBATCH -N 1                        # number of nodes
#SBATCH -A owner-guest              # account
#SBATCH -p ember-guest          # partition (guest is free) other option: ember won't boot, takes forever to run
#SBATCH --mail-user=jcooper036@gmail.com    # 
#SBATCH --requeue       ## put this back in the queue if it does not finish

##### load modules
module load samtools
module load bwa

set -e

###### variables, mostly file paths
dir=/scratch/general/lustre/u0885126/three_genomes
ref=${dir}/ref_genomes/Dmel_r6_180905
scripts=/uufs/chpc.utah.edu/common/home/u0885126/scripts
check=${dir}/check
CORES=12

## wrapper for the make files for geneome alignment
samples=(
    J10075_15364X1
    J10224_15364X2
    J20039_15364X3
    J40022_15364X4
    J40023_15364X5
    J40024_15364X6
    rep1_female_15367X1
    rep1_male_15367X2
    rep2_female_15367X4
    rep2_male_15367X5
    rep3_female_15367X3
    rep3_male_15367X6
)

## turn the fastq files into bam files with read groups
for sample in ${samples[@]}; do
        echo $sample
        export samp=$sample
        make -f ${scripts}/three_genomes/trim_to_bam_kp.mk || exit
done