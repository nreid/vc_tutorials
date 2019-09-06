#!/bin/bash 
#SBATCH --job-name=pileup
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

# make directory if it doesn't exist. 
mkdir -p ../variants_bcftools

# move to alignment directory
cd ../align_stepwise

# make a list of bam files
ls *mkdup.bam >list.bam

# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

bcftools mpileup \
	-f $GEN \
	-b list.bam \
	-q 20 -Q 30 \
	-r chr20:29400000-34400000 >../variants_bcftools/chinesetrio.pileup

date

