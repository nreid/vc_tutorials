#!/bin/bash 
#SBATCH --job-name=gatk_HC
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


# load software
module load GATK/4.0

gcc

module list

# change directory to the bam file location
cd ../align_pipe

# make a list of bam files
ls *bam >bam.list

GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk HaplotypeCaller \
     -R $GEN \
     -I son.bam \
     -ERC GVCF \
     -L ../coverage_stats/targets.bed \
     --output ../variants_gatk/test.$SLURM_JOB_ID.g.vcf
