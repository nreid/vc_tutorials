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

In [Part 2](/Part2_bcftools.md) we used `bcftools` to call variants. If you remember, the two step protocol in `bcftools` is to first summarize the read pileup, base by base, across the genome, and then evaluate the evidence for variation at each site. This means `bcftools` takes the alignments in the read pileup literally. Two factors combine to make this problematic for variant calling. First, indels and complex variants (indels + SNPs) may not have a single unambiguous representation. Second, we align reads (or read pairs) independently. This means that a single underlying biological sequence that differs from the reference may be represented by multiple different alignments in the read pileup. This can lead to both false negative and false positive variant calls. 

`bcftools` deals with this problem by downgrading base qualities to account for uncertainty in the sequence alignment. This reduces the likelihood of false positive calls, but it also reduces the sensitivity to indels and nearby SNP variants. 



## Update your working directory