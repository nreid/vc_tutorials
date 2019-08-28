# Stepwise QC, alignment, post-alignment processing. 


## download data

NIST GIAB data, chinese trio. 
45x coverage bam files for chr20:10000000-15000000

samtools view -bh <file> <region> outputs the given region from the file, includes the header, outputs as bam. 
bedtools bamtofastq converts bam format back to fastq (so we can turn it back into a bam!)


son:

```bash
SON=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-221080150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam

samtools view -bh $SON chr20:10000000-15000000 | \
bedtools bamtofastq -i stdin - -fq son.1.fq -fq2 son.2.fq
```

mother:


father:

