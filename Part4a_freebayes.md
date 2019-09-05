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

In [Part 2](/Part2_bcftools.md) we used `bcftools` to call variants. If you remember, the two step protocol in `bcftools` is to first summarize the read pileup, base by base, across the genome, and then evaluate the evidence for variation at each site. This means `bcftools` takes the alignments in the read pileup literally. 

Two factors combine to make this problematic for variant calling. First, indels and complex variants (indels + SNPs) may not have a single unambiguous representation. Second, we align reads (or read pairs) independently. This means that a single underlying biological sequence that differs from the reference may be represented by multiple different alignments in the read pileup. This can lead to both false negative and false positive variant calls. 

`bcftools` deals with this problem by downgrading base qualities to account for uncertainty in the sequence alignment (using "base alignment qualities"). This reduces the probability of false positive calls, but it also reduces the sensitivity to indels and nearby SNP variants. 

The variant caller we will use here, `freebayes` and another, the GATK's `HaplotypeCaller` (discussed in [Part 4b](Part4b_gatk.md)) take different approaches to this problem. Instead of calling variants at single sites, `freebayes` calls __haplotypes__. Haplotypes are combinations of alleles at different sites found on a single chromosome. To do this, `freebayes` takes in the aligned reads overlapping a short region of the genome, standardizes the representation of indels, identifies candidate haplotypes from the reads (which must be shorter than the read length), and then evaluates the evidence for variation at the site using a Bayesian model. 

INSERT FIGURE

The key innovation here is the use of haplotype alleles. These allow the evidence for indels and complex variants with many possible alignments to be evaluated more rigorously, improving the sensitivity and specificity for these calls. 

`freebayes` has a few other advantages as well. The biggest among them is that it is designed with piping in mind. It can accept a single stream of BAM input merged from multiple samples. That means that filtering or other alignment manipulations can be performed on the fly, delivering data to the caller without generating intermediate alignment files. We will see how this can work below. 


## Update your working directory

First we'll create a new directory to hold the results we'll generate. Make sure you're in the directory vc_workshop and type:

```bash
mkdir -p variants_freebayes coverage_stats
```

## Decide where to call variants

The most common approach to discovering genomic variation, short-read sequencing followed by reference mapping presents a few problems:
- Most genomes are drafts. This means that many regions present in the true genome of the reference individual are absent from the reference sequence. 
- Many genomes contain regions of repetitive or highly similar sequence. 
- The genomic composition of individuals within a species can be highly variable, with some regions copied, sometimes many times, or deleted. 

These issues mean that when you map 







