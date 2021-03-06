#!/bin/bash 
#SBATCH --job-name=tabix
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load htslib


# move to variant
cd ../variants_bcftools

tabix -p vcf chinesetrio.vcf.gz

date