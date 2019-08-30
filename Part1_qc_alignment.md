# Stepwise QC, alignment, post-alignment processing. #

## Introduction

This section of the tutorial deals with quality control, alignment to a reference genome, and post-processing of high-throughput short-read sequencing data in preparation for variant calling. In this section, steps are executed individually. For an example using pipes, see Part 3. 

Steps here will use the following software packages:

- [samtools](http://www.htslib.org/doc/samtools.html)
- [picard tools](https://broadinstitute.github.io/picard/)
- [bwa](http://bio-bwa.sourceforge.net/)
- along with a variety of utilities available on unix-like operating systems. 

Each step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) scheduler. The code can easily be modified to run interactively, or in other contexts. 


## Contents
  
1.    [ Set up a working directory ](#Set-up-a-working-directory)  
2.    [ Prepare reference genome ](#Prepare-Reference-genome)
3.    [ Download data ](#Download-data)
4.    [ Assess read quality ](#Assess-read-quality)
5.    [ Quality trim ](#Quality-trim)
6.    [ Align and compress ](#Align-and-compress)
7.    [ Sort reads by genome position ](#Sort-reads-by-genome-position)
8.    [ Mark duplicates ](#Mark-duplicates)
9.    [ Index alignment files ](#Index-alignment-files)
10.   [ Exploring SAM files ](#Exploring-SAM-files)


## Set up a working directory ##

We will begin the tutorial by setting up a working directory to organize the files we'll generate. From your home directory, enter the code below on the command line. 

```bash
mkdir -p vc_workshop/rawdata vc_workshop/fastqc vc_workshop/align_stepwise vc_workshop/scripts
cd vc_workshop
```

## Prepare Reference genome

Most software packages that align short-read sequencing data or otherwise manipulate a reference genome require that genome to be indexed in some way. We will generate indexes using both `bwa` and `samtools`. For the workshop, a pre-indexed human genome is provided, but the code below shows how it's done:

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
bwa index $GEN
samtools faidx $GEN
```

## Download data ##

For all following steps of the workshop, we'll use data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project, hosted by the [National Institute of Standards and Technology](https://www.nist.gov). The data were collected as part of an effort to create reference standards to compare genomic workflows. Large quantities of sequencing data have been generated on multiple platforms (Illumina, PacBio, etc) for seven individuals. 

We are going to use data from a trio (mother, father, and son) of Chinese ancestry. The data consist of 250bp paired end reads sequenced on an Illumina HiSeq 2500. To ensure the analyses run quickly, we'll only use data from 5 megabases of chromosome 20 (10mb - 15mb). The expected sequencing coverage is 100x for the mother and father, and 45x for the son. 

More information about the data can be found at the links below:
https://www.nist.gov/programs-projects/genome-bottle
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/

samtools view -bh <file> <region> outputs the given region from the file, includes the header, outputs as bam. 
'bedtools bamtofastq' converts bam format back to fastq (so we can practice turning it back into a bam later!)

For example, this code would download data for the son:

```bash
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'

samtools view -uh $SON chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq son.1.fq -fq2 son.2.fq
```


scripts:
- scripts/Part1a_datadownload.sh

## Assess read quality ##

FastQC is used to evaluate the quality of the raw sequencing data. 

scripts: 
- scripts/Part1b_fastqc.sh

## Quality trim ##

not generally necessary for variant calling, but sickle and trimmomatic can be used. 

## Align and compress ##

bwa mem, samtools

scripts:	
- scripts/Part1c_align.sh    
- scripts/Part1d_compress.sh

## Sort reads by genome position ##

samtools

scripts:	
- scripts/Part1e_sort.sh

## Mark duplicates ##

picard

scripts:
- scripts/Part1f_markduplicates.sh

## Index alignment files ##

samtools

scripts:
- scripts/Part1g_indexbams.sh

___

# Exploring SAM files #

_perhaps this should be broken out into a separate sub-section_

Now that we have processed alignment files we can learn about sam files and do some exploring. The ability to inspect and summarize various aspects of an alignment file can be especially helpful when something has gone wrong with sequencing or bioinformatic analysis. 

There are no scripts for the section below. Run the code interactively. 

Add discussion of SAM format here. 

_Useful links:_
- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

- [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)

___

We can get a lot of basic stats on the SAM file using "samtools stats":

```bash
# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools stats
samtools stats -r $GEN son.mkdup.bam >son.samstat.txt
```
This file is rather messy, but we can pull out specific parts of it using grep. 

These are some basic stats about the alignment file:

```bash
grep ^SN son.samstat.txt | cut -f 2-
```

This is a histogram of per base coverage:

```bash
grep ^COV son.samstat.txt | cut -f 2-
```
And there is much more information. 

___

We can look at the per base coverage of individual regions easily using "samtools depth"

```bash
# make a list of bam files
ls *mkdup.bam >list.bam
# run samtools depth using that list
samtools depth -f list.bam -r chr20:13934265-13934558
```

We expect the parents to have 100x coverage and the son to have 50x coverage. Is there anything unusual about this region?

___

Samtools has implemented a simple alignment viewer in 'tview'. IGV is better, but tview is right in the terminal. Let's look at the same region as above. 

```bash
# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools tview
samtools tview --reference $GEN dad.mkdup.bam
```

After opening tview, type 'g' then enter 'chr20:13934265' to visit that location. 
"?" will bring up a help dialog
"q" exits

___

IGV is a much better way of visualizing alignment files...
