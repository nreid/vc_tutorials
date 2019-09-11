#!/bin/bash
#SBATCH --job-name=bcf_vcf
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools/1.6
module load htslib

mkdir -p ../filtered_vcfs

cd ../filtered_vcfs

bcftools filter -s LowQual -e '%QUAL<20' ../variants_freebayes/chinesetrio_fb.vcf.gz | bgzip -c > fb_filter.vcf.gz
bcftools filter -s LowQual -e '%QUAL<20' ../variants_gatk/chinesetrio.vcf | bgzip -c > gatk_filter.vcf.gz
bcftools filter -s LowQual -e '%QUAL<20' ../variants_bcftools/chinesetrio.vcf.gz | bgzip -c > bcf_filter.vcf.gz











module load htslib


# get the dbsnp set for chromosome 20
tabix -h ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz 20 | \
sed 's/^20/chr20/' | \
bgzip -c >chr20.dbsnp.vcf.gz