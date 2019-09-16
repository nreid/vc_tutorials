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

cd ../filtered_vcfs

# get the dbsnp set for chromosome 20
	# do a bunch of reformatting to make sure it works with the following steps
	
	# download only a section of chr20 from dbsnp
	tabix -h ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz 20:28000000-35000000 | \
	sed 's/^20/chr20/' | \
	bgzip -c >chr20.dbsnp.vcf.gz
	# update the sequence dictionary
	gatk UpdateVCFSequenceDictionary -V chr20.dbsnp.vcf.gz --source-dictionary ../align_pipe/son.bam --output chr20.dbsnp.contig.vcf.gz --replace=true
	# make sure all the filter columns say PASS so we can filter below and not lose all the variants
	bcftools filter -s LowQual chr20.dbsnp.contig.vcf.gz | bgzip -c >chr20.dbsnp.filter.vcf.gz
	# index
	tabix -p vcf chr20.dbsnp.filter.vcf.gz
	# remove intermediate files
	rm chr20.dbsnp.contig.vcf.gz chr20.dbsnp.vcf.gz

# set variables for software (need to get these updated/installed for module)
VT=~/bin/vt/vt
VAP=~/bin/vcflib/bin/vcfallelicprimitives

# create a pedigree file:
echo giabct	son	dad	mom	1 -9 >ct.ped
echo giabct	son	0	0	1 -9 >>ct.ped
echo giabct	son	0	0	1 -9 >>ct.ped


$VT partition -f PASS <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz)
$VT partition -f PASS <($VAP fb_filter.vcf.gz) bcf_filter.vcf.gz

$VT multi_partition -f PASS <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz) bcf_filter.vcf.gz
$VT multi_partition -f PASS <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz) bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz
$VT multi_partition <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz) bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz




bcftools isec -p out fb_filter.vcf.gz gatk_filter.vcf.gz bcf_filter.vcf.gz

# all variants
bcftools isec -p gatk_bcf gatk_filter.vcf.gz bcf_filter.vcf.gz
bcftools isec -p fb_gatk fb_filter.vcf.gz gatk_filter.vcf.gz
bcftools isec -p fb_bcf fb_filter.vcf.gz bcf_filter.vcf.gz
bcftools isec -p gatk_db gatk_filter.vcf.gz chr20.filter.dbsnp.vcf.gz
bcftools isec -p fb_db fb_filter.vcf.gz chr20.dbsnp.filter.vcf.gz
bcftools isec -p bcf_db bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz


# only those that PASS (or have no filter flag set)
bcftools isec -f PASS,. -p gatk_bcf_p gatk_filter.vcf.gz bcf_filter.vcf.gz
bcftools isec -f PASS,. -p fb_gatk_p fb_filter.vcf.gz gatk_filter.vcf.gz
bcftools isec -f PASS,. -p fb_bcf_p fb_filter.vcf.gz bcf_filter.vcf.gz
bcftools isec -f PASS,. -p gatk_db_p gatk_filter.vcf.gz chr20.dbsnp.filter.vcf.gz
bcftools isec -f PASS,. -p fb_db_p fb_filter.vcf.gz chr20.dbsnp.filter.vcf.gz
bcftools isec -f PASS,. -p bcf_db_p bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz

# this leads to obvious problems where the variant caller everyone thinks doesn't work as well hits the most snps. 

# must decompose complex variants. use vt

VT=~/bin/vt/vt
VAP=~/bin/vcflib/bin/vcfallelicprimitives
$VT multi_partition -f PASS <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz) bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz
$VT multi_partition <($VAP fb_filter.vcf.gz) <($VAP gatk_filter.vcf.gz) bcf_filter.vcf.gz chr20.dbsnp.filter.vcf.gz

