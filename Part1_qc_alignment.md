# Stepwise QC, alignment, post-alignment processing. #

## Introduction

This section of the tutorial deals with quality control, alignment to a reference genome, and post-processing of high-throughput short-read sequencing data in preparation for variant calling. In this section, steps are executed individually. For an example using pipes, see Part 3. 

Steps here will use the following software packages:

- [samtools](http://www.htslib.org/doc/samtools.html)
- [picard tools](https://broadinstitute.github.io/picard/)
- [bwa](http://bio-bwa.sourceforge.net/)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- along with a variety of utilities available on unix-like operating systems. 
`
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

We're going to use data from a trio (mother, father, and son) of Chinese ancestry. The data consist of 250bp paired end reads sequenced on an Illumina HiSeq 2500. To ensure the analyses run quickly, we'll only use data from 5 megabases of chromosome 20 (10mb - 15mb). The expected sequencing coverage is 100x for the mother and father, and 45x for the son. 

More information about the data can be found at the links below:
https://www.nist.gov/programs-projects/genome-bottle
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/

To download the data, we'll use `samtools`. The data have already been aligned to a reference genome, and the resulting BAM file has been compressed and indexed. This will allow us to get reads only from the region we're interested in. Conveniently, `samtools` can read BAM files from an ftp server, provided the index is present, so we won't need to download the whole dataset. We'll then convert the data back to the unaligned fastq format using `bedtools` so we can continue with the tutorial.

We'll accomplish this with a unix pipeline, where the symbol `|` indicates that the output of the command to the left should be redirected as input to the command to the right. This will be discussed in more detail in Part 3 of the tutorial. 

The commands to be chained together are as follows:

`samtools view -bh <file> <region>` outputs the given `region` from the `file`, includes the header, and outputs as a compressed BAM file.  
`samtools sort -n -` sorts the reads by name so that read pairs will be found together in the file. The `-` indicates the data should be read from the pipe.    
`bedtools bamtofastq -i /dev/stdin/ -fq <read 1> -fq2 <read 2>` converts bam format back to fastq. In this case `/dev/stdin/` indicates the data should be read from the pipe. 

Below is an example of this code put together to download data for the son:

```bash
# set a variable 'SON' that gives the location of the file containing the data from the son:
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'
# download the data, sort it, reformat to fastq
samtools view -uh $SON chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq son.1.fq -fq2 son.2.fq
```

scripts:
- [scripts/Part1a_datadownload.sh](scripts/Part1a_datadownload.sh)

### Inspecting fastq files ###

This set of commands will allow you to open the FASTQ file and inspect the data.  Each read is represented by 4 lines.  Details on the FASTQ file format can be found <a href="https://bio Informatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>.   
 

Change the directory to the `rawdata/` folder using:   
```bash
cd rawdata/  
```  

To inspect the first few lines in the FASTQ file can use the `head` command, as follows which will print the first few lines into the terminal window.     

```bash
head son.1.fq
```

>@61CC3AAXX100125:7:118:2538:5577/1
>GACACCTTTAATGTCTGAAAAGAGACATTCACCATCTATTCTCTTGGAGGGCTACCACCTAAGAGCCTTCATCCCC
>+
>?>CADFEEEBEDIEHHIDGGGEEEEHFFGIGIIFFIIEFHIIIIHIIFFIIIDEIIGIIIEHFFFIIEHIFA@?==
>@61CC3AAXX100125:7:1:17320:13701/1
>CTCAGAAGACCCTGAGAACATGTGCCCAAGGTGGTCACAGTGCATCTTAGTTTTGTACATTTTAGGGAGATATGAG
>+
>?BCAAADBBGGHGIDDDGHFEIFIIIIFGEIFIIFIGIGEFIIGGIIHEFFHHHIHEIFGHHIEFIIEECE?>@89
>@61CC3AAXX100125:7:93:5100:14497/1
>CTCAACTGGCTGAAAGTATTATCAATAGAAAGGAATGTTCAGGTTCTTCAATTTTAGAGTGCCCTGGCCTAGAAGA
 
The command below will let you inspect more part of the file.
```bash
less son.1.fq
```   

 
The command below will help you count number of reads in the file
``` 
-bash-4.2$ grep -c "^@61CC" NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
```   


## Assess read quality ##

In order to evaluate the general quality of reads in the file we will be using `FASTQC` package.  The command can take multiple files as input and outputs a quality report in html format. Once the file is generated you have to transfer the file to your local computer to open it and examine the results carefully. To copy the file from the Xanadu cluster please use the `transfer.cam.uchc.edu` node.  

scripts: 
- [scripts/Part1b_fastqc.sh](scripts/Part1b_fastqc.sh)

## Quality trim ##

not generally necessary for variant calling, but sickle and trimmomatic can be used. 

## Align and compress ##

bwa mem, samtools

scripts:	
- [scripts/Part1c_align.sh](scripts/Part1c_align.sh)    
- [scripts/Part1d_compress.sh](scripts/Part1d_compress.sh)

## Sort reads by genome position ##

samtools

scripts:	
- [scripts/Part1e_sort.sh](scripts/Part1e_sort.sh)

## Mark duplicates ##

picard

scripts:
- [scripts/Part1f_markduplicates.sh](scripts/Part1f_markduplicates.sh)

## Index alignment files ##

samtools

scripts:
- [scripts/Part1g_indexbams.sh](scripts/Part1g_indexbams.sh)

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
# set a variable 'GEN' that gives the location of the reference genome:
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
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools tview
samtools tview --reference $GEN dad.mkdup.bam
```

After opening tview, type 'g' then enter 'chr20:13934265' to visit that location. 
"?" will bring up a help dialog
"q" exits

___

IGV is a much better way of visualizing alignment files...
