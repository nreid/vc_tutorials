# Variant discovery using freebayes

## Introduction

This section of the tutorial introduces variant calling using `freebayes`. It assumes you have completed Part 1 and Part 3 and have QC'ed, aligned, post-processed bam files in the directory structure first created in Part 1. 

Steps here will use the following software packages:

- [ freebayes ](https://github.com/ekg/freebayes)
- [ bedtools ](https://bedtools.readthedocs.io/en/latest/)
- [ bamtools ](https://github.com/pezmaster31/bamtools)
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Decide where to call variants ]()
3.    [ Call variants ]()


## Motivation


## Update your working directory