#!/bin/bash 
#SBATCH --job-name=fb_mono
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



hostname
date

# load required software

module load bedtools
module load bamtools
module load htslib
module load freebayes


# make a directory for results if it doesn't exist
mkdir -p ../variants_freebayes 

# change directory to the bam file location
cd ../align_pipe

# make a list of bam files
ls *bam >bam.list

# set a variable for the reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

OUTLIERWINDOWS=../coverage_stats/coverage_outliers.bed.gz

# note that bamtools region specification uses ".." instead of "-"
bamtools merge -list bam.list -region chr20:29400000..34400000 | \
bamtools filter -in stdin -mapQuality ">30" -isProperPair true | \
bedtools intersect -v -a stdin -b $OUTLIERWINDOWS | \
freebayes -f $GEN --stdin --report-monomorphic | \
bgzip -c >../variants_freebayes/chinesetrio_fb_mono.vcf.gz

tabix -p vcf ../variants_freebayes/chinesetrio_fb_mono.vcf.gz

date
