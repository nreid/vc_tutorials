# Stepwise QC, alignment, post-alignment processing. #

## Setting up working directory ##
```bash
mkdir -p vc_workshop/rawdata vc_workshop/fastqc vc_workshop/align_stepwise vc_workshop/scripts
cd vc_workshop
```

## download data ##

NIST GIAB data, chinese trio. 
son: 45x coverage bam files for chr20:10000000-15000000
mother & father: 100x coverage bam files for chr20:10000000-15000000

Data from:
https://www.nist.gov/programs-projects/genome-bottle
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/

samtools view -bh <file> <region> outputs the given region from the file, includes the header, outputs as bam. 
bedtools bamtofastq converts bam format back to fastq (so we can turn it back into a bam!)


son:

```bash
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'

samtools view -uh $SON chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq son.1.fq -fq2 son.2.fq
```

mother:

```bash
MOM=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam

samtools view -uh $MOM chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq mom.1.fq -fq2 mom.2.fq
```

father:

```bash
DAD='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam'

samtools view -uh $DAD chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq dad.1.fq -fq2 dad.2.fq
```

## Assessing read quality ##



## Quality trimming ##


## Alignment ##


## Marking duplicates ##


## Sorting ##


## Exploring SAM files ##

Explain SAM flags:
https://broadinstitute.github.io/picard/explain-flags.html
