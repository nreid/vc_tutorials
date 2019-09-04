# Variant annotation

## Introduction


## Introduction


Steps here will use the following software packages:

-
- 

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  

## Motivation


## Update your working directory


## Variant normalization and decomposition

vcflib, vt, vcfeval (https://www.realtimegenomics.com/products/rtg-tools)
vgraph (https://github.com/bioinformed/vgraph)
bcftools normalize
https://github.com/Illumina/hap.py

?? how do we discuss annotation ?? 
?? how do we compare variant call sets ??

discussion of normalization:
https://genome.sph.umich.edu/wiki/Variant_Normalization


bcftools annotate does simple matching. does not account for complex alleles or incompatible representations. 