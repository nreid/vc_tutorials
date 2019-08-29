# Alignment, post-alignment processing using a pipe


Here we will see how we can get most of our big alignment and processing tasks from part 1 done in one fell swoop without creating many intermediate files. 

to do this we make use of a core feature of unix-like operating systems: the pipe. 

the pipe "|" allows you to redirect the output of one command to the input of another. in unix-speak, you would say that you are redirecting the "standard output" stream of program 1 to the "standard input" of progam 2. 

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