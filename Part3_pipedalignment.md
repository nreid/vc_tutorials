# Alignment, post-alignment processing using a pipe

## Introduction

In Part 1 of the tutorial, we downloaded our data QC'ed it, aligned it to a reference genome, and did post-alignment processing. Each step was conducted using a discrete script. Each script required that all the data be read, modified, and then re-written. This creates a couple problems. For non-trivial datasets, all of this reading and writing takes time. It also generates many copies of the data. If our first approach had been applied to a 300 gigbyte dataset, by the time we finished we'd be using somewhere in the neighborhood of 1.8 terabytes of disk space. 

In this section of the tutorial we'll see how we can streamline this approach. To do this we'll make use of a core feature of unix-like operating systems: the pipe. 

The pipe, `|` allows you to redirect the output of one command to the input of another. In unix-speak, you would say that you are redirecting the "standard output" stream of program 1 to the "standard input" of progam 2. 

in this way, you can chain together many specialized programs to achieve some goal without having to repeatedly read and write files. this is helpful because you can avoid making a big mess of intermediate files that you don't need and, in the case of high throughput sequencing data, take up quite a lot of space. 

many programs in bioinformatics are designed with pipes in mind, including core u

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