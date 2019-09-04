# Alignment, post-alignment processing using a pipe

## Introduction

In Part 1 of the tutorial, we downloaded data QC'ed it, aligned it to a reference genome, and did post-alignment processing. Each step was conducted using a discrete script. Each script required that all the data be read, modified, and then re-written. This creates a couple problems. For non-trivial datasets, all of this reading and writing takes time. It also generates many copies of the data. If our first approach had been applied to a 300 gigbyte dataset, by the time we finished we'd be using somewhere in the neighborhood of 1.8 terabytes of disk space. 

In this section of the tutorial we'll see how we can streamline this approach. To do this we'll make use of a core feature of unix-like operating systems: the pipe. 

The pipe, `|`, allows you to [redirect the output of one command to the input of another](https://en.wikipedia.org/wiki/Standard_streams). In unix-speak, you would say that you are redirecting the "standard output" stream of command 1 to the "standard input" of command 2. 

In this way, you can chain together many specialized programs to achieve some goal without having to repeatedly read and write files. This mitigates both of the problems described above. 

Many programs in bioinformatics are designed with pipes in mind, including core utilities such as `samtools` and the variant caller `Freebayes` we'll use in Part 4. 

Steps here will use the following software packages:

- [ samtools ](http://www.htslib.org/doc/samtools.html)
- [ bwa ](http://bio-bwa.sourceforge.net/)
- [ samblaster ](https://github.com/GregoryFaust/samblaster)

We'll assume we have already downloaded and QC'ed fastq files as in steps 1-5 from [Part 1](/Part1_qc_alignment.md). 

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Piped alignment and post-processing ](#Piped-alignment-and-post-processing)
3.    [ Index alignment files ](#Index-alignment-files)

## Update your working directory


## Piped alignment and post-processing


## Index alignment files


assume we start with QC'ed fastq files. 

```bash
bwa mem | samblaster | samtools sort >sample.bam.
```

using bwa mem we first align and tag sequences with read groups. 

next we stream the output using a pipe to samblaster, which marks duplicates. it can also identify split and discordantly mapping reads and write them to a separate file. useful for structural variation identification.

finally we stream that output to samtools sort and write a compressed bam file. 

instead of reading and writing our sequences 5 times, resulting in 6 copies of our sequence data, we now only have two copies, the original fastq files, and processed bam files. 

if we have multiple sets of fastq files for each sample, we can then merge them into a single file. this is convenient, although not strictly necessary: both GATK and Freebayes identify which reads belong to which samples by the read group tags, and not the bam file from which they originated. 

last, we index the files. this allows accessing reads from specific genomic regions without scanning the entire file. 