# Variant discovery using GATK

## Introduction

This section of the tutorial introduces variant calling using `GATK`. It assumes you have completed Part 1 and Part 3 and have QC'ed, aligned, post-processed bam files in the directory structure first created in Part 1. 

Steps here will use the following software packages:

- [ GATK ](https://software.broadinstitute.org/gatk/)
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)
- [ tabix ](http://www.htslib.org/doc/tabix.html)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 


## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  
3.    [ Generate single sample GVCF files ]()
4.    [ Consolidate GVCF files]()
5.    [ Do joint variant calling on the GVCFs ]()

## Motivation

In [Part 4a](/Part4a_freebayes.md) we used `freebayes` to call variants. We learned that a key advantage of `freebayes` over `bcftools` is the use of haplotype alleles to overcome uncertainty and stochasticity in the representation of complex variants. Here we're going to use a variant caller implemented in the [Broad Institute's Genome Analysis Toolkit, a.k.a. GATK](https://software.broadinstitute.org/gatk/) called `HaplotypeCaller`. This tool also calls haplotype alleles, but its approach varies slightly from that of `freebayes`. Where freebayes identifies candidate alleles from the reference alignments, `HaplotypeCaller` does full _de novo_ assembly from reads within a small region to find candidate alleles, then aligns the reads from individual samples to those alleles to call genotypes. This is more computationally intensive, and as such, `HaplotypeCaller` is generally a bit slower than `freebayes`. The end results are similar, and both approaches improve sensitivity and specificity for indels and nearby variants. 

`GATK` is a much broader set of tools than just the variant calling algorithm. It actively developed and supported by a team at the Broad and they regularly update a set of best practices for using their tools. The approach they recommend for calling variants on multiple samples differs in a couple ways from others:

1. For model systems they have tools that recalibrate the quality scores for
	a. bases in aligned sequence data (base quality score recalibration)
	b. called variants (variant quality score recalibration)
 Recalibrating the scores simply means that the phred-scaled error probabilities they represent become more accurate. These tools require training data in the form of independently ascertained variants that are likely to occur in your samples. These are easily obtainable for model systems, but not so much in others. 
2. They recommend a two-step procedure for multi-sample variant calling. This procedure makes scaling to calling variants on very large cohorts more efficient. 

Here we will not demonstrate the quality-score recalibration, but we will use the two-step multi-sample variant calling procedure. 

## Update your working directory


## Call variants on single samples to generate GVCF files

## Consolidate GVCF files

## Jointly genotype the GVCFs to generate a VCF file

