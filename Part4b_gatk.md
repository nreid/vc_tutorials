# Variant discovery using GATK

## Introduction

This section of the tutorial introduces variant calling using `GATK`. It assumes you have completed Part 1, Part 3 and Part4a and have QC'ed, aligned, post-processed bam files in the directory structure first created in Part 1. 

Steps here will use the following software packages:

- [ GATK ](https://software.broadinstitute.org/gatk/)
- [ bedtools ](https://bedtools.readthedocs.io/en/latest/)
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)
- [ tabix ](http://www.htslib.org/doc/tabix.html)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 


## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  
3.    [ Call variants on single samples to generate GVCF files ](#Call-variants-on-single-samples-to-generate-GVCF-files)
4.    [ Consolidate GVCF files ](#Consolidate-GVCF-files)
5.    [ Jointly genotype the GVCFs to generate a VCF file ](#Jointly-genotype-the-GVCFs-to-generate-a-VCF-file)

## Motivation

In [Part 4a](/Part4a_freebayes.md) we used `freebayes` to call variants. We learned that a key advantage of `freebayes` over `bcftools` is the use of haplotype alleles to overcome uncertainty and stochasticity in the representation of complex variants. Here we're going to use a variant caller implemented in the [Broad Institute's Genome Analysis Toolkit, a.k.a. GATK](https://software.broadinstitute.org/gatk/) called `HaplotypeCaller`. This tool also calls haplotype alleles, but its approach varies slightly from that of `freebayes`. Where freebayes identifies candidate alleles from the reference alignments, `HaplotypeCaller` does full _de novo_ assembly from reads within a small region to find candidate alleles, then aligns the reads from individual samples to those alleles to call genotypes. This is more computationally intensive, and as such, `HaplotypeCaller` is generally a bit slower than `freebayes`. The end results are similar, and both approaches improve sensitivity and specificity for indels and nearby variants. 

`GATK` is a much broader set of tools than just the variant calling algorithm. It actively developed and supported by a team at the Broad and they regularly update a set of [best practices](https://software.broadinstitute.org/gatk/best-practices/) for using their tools. The approach they recommend for calling variants on multiple samples differs in a couple ways from others:

1. For model systems they have tools that recalibrate the quality scores for:  
    a. Bases in aligned sequence data (base quality score recalibration).  
    b. Called variants (variant quality score recalibration).  

 Recalibrating the scores simply means that the phred-scaled error probabilities they represent become more accurate. These tools require training data in the form of independently ascertained variants that are likely to occur in your samples. These are easily obtainable for model systems, but not so much in others. 

2. They have developed a three-step procedure for multi-sample variant calling. According to the GATK team, this procedure makes scaling to calling variants on very large cohorts more efficient. The biggest advantage is that in no part of the pipeline is all the read data for all the individuals read into memory simultaneously. There are, however, no published tests (that we're aware of) demonstrating that it is faster, and if so, on what scale, or whether it exhibits better statistical properties than joint calling all samples at once, e.g. in `freebayes`. 

Here we will not demonstrate the quality-score recalibration, but we will use the three-step multi-sample variant calling procedure. In [Part 5](/Part5_variant_annotation.md), we will annotate and compare variant calls from all three approaches in the tutorial. 

## Update your working directory

First we'll create a new directory to hold the results we'll generate. Make sure you're in the directory vc_workshop and type:

```bash
mkdir -p variants_gatk
```

## Call variants on single samples to generate GVCF files

You can jointly call variants across multiple samples using HaplotypeCaller in a single step, as we did with `freebayes`, but for this tutorial, we'll demonstrate the multi-step approach. The steps are: 

1. Summarize sequence variation across the genome for each sample individually in a GVCF file (`HaplotypeCaller -ERC GVCF`).
2. Combine the GVCF files into a database (`GenomicsDBImport`).
3. Jointly genotype the samples (`GenotypeGVCFs`)

The first step must be run independently for each sample, e.g.:

```bash
# set a variable for the reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk HaplotypeCaller \
     -R $GEN \
     -I son.bam \
     -ERC GVCF \
     -L ../coverage_stats/targets.bed \
     --output ../variants_gatk/son.g.vcf
```

Here we've provided a single BAM file and called `HaplotypeCaller` with the flag `-ERC GVCF` to produce the GVCF file. We've also provided a set of intervals to analyze that are the inverse of the outlier windows we identified in Part 4a. The script we submit to the scheduler runs this command on each BAM file separately. 

___
scripts:
[scripts/Part4c_gatk_gvcf.sh](scripts/Part4c_gatk_gvcf.sh)

## Consolidate GVCF files

In the next step, we consolidate the GVCF files into a single database. 

```bash
gatk --java-options "-Xmx10g -Xms4g" GenomicsDBImport \
  -V ../variants_gatk/mom.g.vcf \
  -V ../variants_gatk/dad.g.vcf \
  -V ../variants_gatk/son.g.vcf \
  --genomicsdb-workspace-path ../variants_genomicsdb \
  --overwrite-existing-genomicsdb-workspace true \
  -L chr20:29000000-35000000
```

With many samples, you can, instead of calling `-V` for each one, provide a table of sample names and file paths. 

___
scripts:
[scripts/Part4d_gatk_genomicsDBimport.sh](scripts/Part4d_gatk_genomicsDBimport.sh)

## Jointly genotype the GVCFs to generate a VCF file

Finally, we jointly genotype each sample, producing a VCF file. 

```bash
# set a variable for the reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk GenotypeGVCFs \
    -R $GEN \
    -V gendb://../variants_genomicsdb \
    -O ../variants_gatk/chinesetrio.vcf 
```

With all variant calling approaches, high sequencing depth and/or many samples

___
scripts:
[scripts/Part4d_gatk_genotype.sh](scripts/Part4d_gatk_genotype.sh)

