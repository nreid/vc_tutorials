#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load bwa/0.7.17
module load samtools

mkdir -p ../align_stepwise

# set input and output directories
indir=../rawdata
outdir=../align_stepwise

# current location of indexed HG38
# may need to be changed. 
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38

# note that read group info is added during alignment. 

# each line aligns one family member's sequences
# son
bwa mem -t 12 -R '@RG\tID:son\tSM:son' $GEN $indir/son.1.fq $indir/son.2.fq -o $outdir/son.sam 
date
# mom
bwa mem -t 12 -R '@RG\tID:mom\tSM:mom' $GEN $indir/mom.1.fq $indir/mom.2.fq -o $outdir/mom.sam 
date
# dad
bwa mem -t 12 -R '@RG\tID:dad\tSM:dad' $GEN $indir/dad.1.fq $indir/dad.2.fq -o $outdir/dad.sam 
date

