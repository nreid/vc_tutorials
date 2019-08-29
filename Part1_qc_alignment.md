# Stepwise QC, alignment, post-alignment processing. #

## Setting up working directory ##
```bash
mkdir -p vc_workshop/rawdata vc_workshop/fastqc vc_workshop/align_stepwise vc_workshop/scripts
cd vc_workshop
```

## Download data ##

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

scripts:
- scripts/Part1a_datadownload.sh

## Assess read quality ##

FastQC is used to evaluate the quality of the raw sequencing data. 

scripts: 
- scripts/Part1b_fastqc.sh

## Quality trim ##

not generally necessary, but sickle, trimmomatic

## Alignment and compression ##

bwa mem, samtools

scripts:	
- scripts/Part1c_align.sh<br>
- scripts/Part1d_compress.sh

## Sort reads by genome position ##

samtools

scripts:	
- scripts/Part1e_sort.sh

## Marking duplicates ##

picard

scripts:
- scripts/Part1f_markduplicates.sh

## Index alignment files ##

samtools

scripts:
- scripts/Part1g_indexbams.sh

## Explore SAM files ##

Now that we have processed alignment files we can learn about sam files and do some exploring. The ability to inspect and summarize various aspects of an alignment file can be especially helpful when something has gone wrong with sequencing or bioinformatic analysis. 

There are no scripts for the section below. Run the code interactively. 

Add discussion of SAM format here. 

_Useful links:_
- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

- [Explain SAM flags (column 2)](https://broadinstitute.github.io/picard/explain-flags.html)

We can get a lot of basic stats on the SAM file using samtools stats:

```bash
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
samtools stats -r $GEN son.mkdup.bam >son.samstat.txt
```
This file is rather messy, but we can remove specific parts of it using grep. 

These are some basic stats about the alignment file:

```bash
grep ^SN son.samstat.txt | cut -f 2-
```

This is a histogram of per base coverage:

```bash
grep ^COV son.samstat.txt | cut -f 2-
```
And there is much more information. 

We can look at the per base coverage of individual regions easily using "samtools depth"

```bash
ls *mkdup.bam >list.bam
samtools depth -f list.bam -r chr20:13934265-13934558
```

We expect the parents to have 100x coverage and the son to have 50x coverage. Is there anything unusual about this region?

Samtools has implemented a simple alignment viewer in 'tview'. IGV is better, but tview is right in the terminal. Let's look at this region. 

```bash
# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools tview
samtools tview --reference $GEN dad.mkdup.bam
```

After opening tview, type 'g' then enter 'chr20:13934265' to visit a particular location. 
"?" will bring up a help dialog
"q" exits


igv
