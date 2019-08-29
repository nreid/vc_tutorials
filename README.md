# Variant discovery tutorials

In this repository, re-working variant detection tutorials for UConn CBC workshop. 

__Proposed structure:__

1. Stepwise QC, alignment, post-alignment processing. 

2. Variant Calling: samtools mpileup. 

3. Part 1, but a piped example. 

4. 
	a. Variant calling: Freebayes, post-filtering. 

	b. Variant calling: GATK, joint calling using gvcf. post-filtering. 

	c. Beyond variant calling: genotype likelihoods. (depending on audience?)

5. Variant annotation. 

__Proposed data:__

NIST Genome in a Bottle asian trio. chr20:10000000-15000000

an arbitrary 5mb region of the genome. 
- takes only a few minutes to align. 
- 50-100x coverage for each of 3 individuals. 

source:
https://www.nist.gov/programs-projects/genome-bottle

ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
