# Variant discovery tutorials

In this repository, re-working variant detection tutorials for UConn CBC workshop. 

## Introduction

This repository is a usable, publicly available tutorial for introduction to the basics of variant calling from high-throughput short-read sequencing data. Most steps have associated bash scripts tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) scheduler that can be modified to run interactively or with another job scheduler.  

Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should submit each script to the scheduler as `sbatch scriptname.sh` after modifying it as appropriate.  

Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this tutorial, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, you can get one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

__Proposed structure:__

1. [ Stepwise QC, alignment, post-alignment processing ](/Part1_qc_alignment.md)

2. [ Variant Calling: bcftools ](/Part2_bcftools.md)

3. Part 1, but a piped example

4. 
	a. Variant calling: Freebayes, post-filtering

	b. Variant calling: GATK, joint calling using gvcf. post-filtering

	c. Beyond variant calling: genotype likelihoods. (depending on audience?)

5. Variant annotation. 

__Proposed data:__

NIST Genome in a Bottle asian trio. chr20:10000000-15000000

an arbitrary 5mb region of the genome. 
- takes only a few minutes to align. 
- 50-100x coverage for each of 3 individuals. 

source:
- https://www.nist.gov/programs-projects/genome-bottle
- ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
