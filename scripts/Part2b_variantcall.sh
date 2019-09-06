#!/bin/bash 
#SBATCH --job-name=variantcall
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

module load bcftools


# move to variant
cd ../variants_bcftools

bcftools call -m -v -Oz -o chinesetrio.vcf.gz chinesetrio.pileup

date