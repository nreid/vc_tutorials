#!/bin/bash 
#SBATCH --job-name=gatk_HC
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



hostname
date

# make sure partition is specified as `xeon` to prevent slowdowns on amd processors. 

# load required software

module load GATK/4.0
module load htslib
module load bedtools

mkdir -p ../variants_gatk

# make the targets file
cd ../coverage_stats
# select the regions we wish to call variants on, merge them using bedtools
tabix ../coverage_stats/coverage_1kb.bed.gz chr20:29400000-34400000  | awk '$6 > 850 && $6 < 2550' | bedtools merge >../coverage_stats/targets.bed
 # it turns out we could use the outlier bed with -XL instead of making a target file and using -L

# change directory to the bam file location
cd ../align_pipe

# make a list of bam files
ls *bam >bam.list

# set a variable for the reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk HaplotypeCaller \
     -R $GEN \
     -I son.bam \
     -ERC GVCF \
     -L ../coverage_stats/targets.bed \
     --output ../variants_gatk/son.g.vcf

date 

gatk HaplotypeCaller \
     -R $GEN \
     -I mom.bam \
     -ERC GVCF \
     -L ../coverage_stats/targets.bed \
     --output ../variants_gatk/mom.g.vcf

date 

gatk HaplotypeCaller \
     -R $GEN \
     -I dad.bam \
     -ERC GVCF \
     -L ../coverage_stats/targets.bed \
     --output ../variants_gatk/dad.g.vcf

date

