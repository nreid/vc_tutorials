#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o ../log_files/%x_%j.out
#SBATCH -e ../log_files/%x_%j.err

hostname
date

# load software
module load samtools

# set directory for input and output files
DIR=../align_stepwise

#SAM to BAM CONVERSION
module load samtools/1.7
# son
samtools view -bhS $DIR/son.sam >$DIR/son.bam
# mom
samtools view -bhS $DIR/mom.sam >$DIR/mom.bam
# dad
samtools view -bhS $DIR/dad.sam >$DIR/dad.bam

