#!/bin/bash 
#SBATCH --job-name=quality_control
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load fastqc

# run fastqc
fastqc -t 4 -o ../fastqc ../rawdata/*fq
