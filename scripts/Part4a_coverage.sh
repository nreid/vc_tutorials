#!/bin/bash 
#SBATCH --job-name=coverage
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

module load bedtools
module load bamtools
module load htslib

mkdir -p ../coverage_stats
cd ../coverage_stats

# genome index file from samtools faidx
FAI=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta.fai

# make a "genome" file, required by bedtools makewindows command, set variable for location
cut -f 1-2 $FAI > Homo_sapiens_assembly38.fasta.genome
GFILE=../coverage_stats/Homo_sapiens_assembly38.fasta.genome

# make 1kb window bed file, set variable for location
bedtools makewindows -g $GFILE -w 1000 >hg38_1kb.bed
WIN1KB=../coverage_stats/hg38_1kb.bed

# go to sequence alignment folder
cd ../align_pipe

# make a list of bam files
ls *bam >bam.list

# pipe:
	# 1) merge bam files
	# 2) filter by quality and proper pairing
	# 3) convert alignments to bed format
	# 4) map alignments to 1kb windows, counting (but also getting the mean and median of the mapping quality score from column 5)

bamtools merge -list bam.list | \
bamtools filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
bedtools bamtobed -i stdin | \
bedtools map \
-a $WIN1KB \
-b stdin \
-c 5 -o mean,median,count \
-g $GFILE \
>../coverage_stats/coverage_1kb.bed

# bgzip compress and tabix index the resulting file
bgzip ../coverage_stats/coverage_1kb.bed
tabix -p bed ../coverage_stats/coverage_1kb.bed.gz

# select and merge outlier windows
zcat ../coverage_stats/coverage_1kb.bed.gz | awk '$6 < 850 || $6 > 2250' | bedtools merge | bgzip >../coverage_stats/coverage_outliers.bed.gz 

date

