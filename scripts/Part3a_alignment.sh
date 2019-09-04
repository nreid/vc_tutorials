#!/bin/bash 
#SBATCH --job-name=alignment
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load samtools
module load samblaster
module load bwa

# make sure alignment directory exists
mkdir -p ../align_pipe

# set a variable 'GEN' that gives the location and base name of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38

# execute the pipe for the son:
bwa mem -t 4 -R '@RG\tID:son\tSM:son' $GEN ../rawdata/son.1.fq ../rawdata/son.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/son - >../align_pipe/son.bam

# execute the pipe for the mom:
bwa mem -t 4 -R '@RG\tID:mom\tSM:mom' $GEN ../rawdata/mom.1.fq ../rawdata/mom.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/mom - >../align_pipe/mom.bam

# execute the pipe for the dad:
bwa mem -t 4 -R '@RG\tID:dad\tSM:dad' $GEN ../rawdata/dad.1.fq ../rawdata/dad.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/dad - >../align_pipe/dad.bam


