# Variant discovery using samtools. 

## Introduction

This section of the tutorial introduced variant calling using the methods implemented in `bcftools`. 

We will assume that you have completed Part 1 and have available QC'ed and processed sequence alignment files in bam format and the directory structure first created there. 

Steps here will use the following software packages:

- [bcftools](http://www.htslib.org/doc/bcftools.html)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Generate a pileup file ](#Generate-a-pileup-file)  
3.    [ Call variants ](#Call-variants)  
4.    [ The VCF format ](#The-VCF-format)  
5.    [ Filter variants ](#Filter-variants)  
6.    [ Compress and index the VCF file ](#Compress-and-index-the-VCF-file)  


## Update your working directory

First we'll make a new directory to hold variant calls. 

```bash
mkdir variants_bcftools
```

## Generate a pileup file

`bcftools` uses a two step procedure to call variants. First it generates a pileup file using `bcftools mpileup`. The pileup file summarizes the per base coverage at each site in the genome. Each row represents a genomic position, and each position that has sequencing coverage is present in the file. The second step, using `bcftools call` actually applies a statistical model to evaluate the evidence for variation represented in that summary and generates output in the variant call format (VCF). 

Because pileup files for whole genome sequencing contain every site in the genome, they can be huge. We also don't typically need to look at them, or do anything with them, so if you use this approach on your own data, you will usually simply pipe the output of `bcftools mpileup` directly to `bcftools call` (for further discussion of the use of pipes, see section 3). 

For demonstration purposes here, we'll keep the steps separate. 

Here is a call to `bcftools mpileup`:

```bash
# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

bcftools mpileup \
	-f $GEN \
	-b list.bam \
	-q 20 -Q 30 \
	-r chr20:10000000-15000000 >../variants_bcftools/chinesetrio.pileup
```

We give `bcftools` the reference genome with `-f`, a list of bam files with `-b`, tell it to exclude bases with quality lower than 30 and reads with mapping quality lower than 20 with `-q` and `-Q` and ask it to generate the pileup only for the region we're focusing on here with `-r`. 

___
scripts:
- [scripts/Part2a_mpileup.sh](scripts/Part2a_mpileup.sh)

## Call variants

The second step is to evaluate the evidence (summarized in the pileup file) that the sequence variation we observe is true biological variation, and not errors introduced during library preparation and sequencing. Here we use `bcftools call`. 

```bash
bcftools call -m -v -Oz -o chinesetrio.vcf.gz chinesetrio.pileup
```

We use the `-o` flag to indicate the output file name. The flag `-m` specifies one of two possible variant calling routines, `-v` says that only variable sites should be output, and `-Oz` indicates the output should be compressed in a version of the gzip compression format. 

___
scripts:
- [scripts/Part2b_variantcall.sh](scripts/Part2b_variantcall.sh)


## The VCF format


## Filter variants


## Compress and index the VCF file

