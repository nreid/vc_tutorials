# Stepwise QC, alignment, post-alignment processing. #

## Introduction

This section of the tutorial deals with quality control, alignment to a reference genome, and post-processing of high-throughput short-read sequencing data in preparation for variant calling. In this section, steps are executed individually. For an example using pipes, see Part 3. 

Steps here will use the following software packages:

- [samtools](http://www.htslib.org/doc/samtools.html)
- [picard tools](https://broadinstitute.github.io/picard/)
- [bwa](http://bio-bwa.sourceforge.net/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [sickle](https://github.com/najoshi/sickle)
- along with a variety of utilities available on unix-like operating systems.    

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 


## Contents
  
1.    [ Set up a working directory ](#Set-up-a-working-directory)  
2.    [ Prepare reference genome ](#Prepare-Reference-genome)
3.    [ Download data ](#Download-data)
4.    [ Assess read quality ](#Assess-read-quality)
5.    [ Quality trim ](#Quality-trim)
6.    [ Align and compress ](#Align-and-compress)
7.    [ Sort reads by genome position ](#Sort-reads-by-genome-position)
8.    [ Mark duplicates ](#Mark-duplicates)
9.    [ Index alignment files ](#Index-alignment-files)
10.   [ Exploring SAM files ](#Exploring-SAM-files)


## Set up a working directory ##

We will begin the tutorial by setting up a working directory to organize the files we'll generate. From your home directory, enter the code below on the command line. 

```bash
mkdir -p vc_workshop/rawdata vc_workshop/fastqc vc_workshop/align_stepwise vc_workshop/scripts
cd vc_workshop
```

## Prepare Reference genome

Most software packages that align short-read sequencing data or otherwise manipulate a reference genome require that genome to be indexed in some way. We will generate indexes using both `bwa` and `samtools`. For the workshop, a pre-indexed human genome is provided, but the code below shows how it's done:

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
bwa index $GEN
samtools faidx $GEN
```

## Download data ##

For all following steps of the workshop, we'll use data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project, hosted by the [National Institute of Standards and Technology](https://www.nist.gov). The data were collected as part of an effort to create reference standards to compare genomic workflows. Large quantities of sequencing data have been generated on multiple platforms (Illumina, PacBio, etc) for seven individuals. 

We're going to use data from a trio (mother, father, and son) of Chinese ancestry. The data consist of 250bp paired end reads sequenced on an Illumina HiSeq 2500. To ensure the analyses run quickly, we'll only use data from 5 megabases of chromosome 20 (10mb - 15mb). The expected sequencing coverage is 100x for the mother and father, and 45x for the son. 

More information about the data can be found at the links below:  
- https://www.nist.gov/programs-projects/genome-bottle
- ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
___

To download the data, we'll use `samtools`. The data have already been aligned to a reference genome, and the resulting BAM file has been compressed and indexed. This will allow us to get reads only from the region we're interested in. Conveniently, `samtools` can read BAM files from an ftp server, provided the index is present, so we won't need to download the whole dataset. We'll then convert the data back to the unaligned fastq format using `bedtools` so we can continue with the tutorial.

We'll accomplish this with a unix pipeline, where the symbol `|` indicates that the output of the command to the left should be redirected as input to the command to the right. This will be discussed in more detail in Part 3 of the tutorial. 

The commands to be chained together are as follows:

- `samtools view -bh <file> <region>` outputs the given `region` from the `file`, includes the header, and outputs as a compressed BAM file.  
- `samtools sort -n -` sorts the reads by name so that read pairs will be found together in the file. The `-` indicates the data should be read from the pipe.    
- `bedtools bamtofastq -i /dev/stdin/ -fq <forward reads> -fq2 <reverse reads>` converts bam format back to fastq. In this case `/dev/stdin/` indicates the data should be read from the pipe. 

The result will be two fastq files, one containing the forward reads and one containing reverse reads with members of each pair on the same lines of their corresponding files. 

Below is an example of this code put together to download data for the son:

```bash
# set a variable 'SON' that gives the location of the file containing the data from the son:
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'
# download the data, sort it, reformat to fastq
samtools view -uh $SON chr20:10000000-15000000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq son.1.fq -fq2 son.2.fq
```
___
scripts:
- [scripts/Part1a_datadownload.sh](scripts/Part1a_datadownload.sh)

### Inspecting fastq files ###

This set of commands will allow you to open the FASTQ file and inspect the data.  Each read is represented by 4 lines.  Details on the FASTQ file format can be found <a href="https://bio Informatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>.   
 

Change the directory to the `rawdata/` folder using:   
```bash
cd rawdata/  
```  

To inspect the first few lines in the FASTQ file can use the `head` command, as follows which will print the first few lines into the terminal window.     

```bash
head son.1.fq
```
A single fastq sequence entry looks like this:
```
@H2YVCBCXX:1:1:383:0/1
GCCCTGCCACACTTAGTGTATTGGCTTTTGTCTTCAGGGTTTTCACCTCGTGTTTTTGATATGGCTACCACGACTCCTAGACTACCTCATTACACAGCATTTACAGGCCAGGATAGGCTCCTCTTCATGTCATTATATGACAATAATTTTCATGTCCTTTTTAAAGCATAAAATAAAAATTCCCACACACCTCCAATAAAACTTTCTGTTCTGTCCCATAAGCTAGAACTGCATCTCTTGCCCATTCCT
+
DDDDDIIIIIHIHIIIEHHIIIIIIIGHHIIIIIIIIIGGHHIIHHIIIIIIHIGIHHHIHIFIHIHHHHIGHGHIIHHIHHIHIIIIHIIIIIIGIIEFHDCHHHIIIHIHHIIIHIHIIIIIIIHIGIIGHHIIIGEHHGGIGGHIECEHHIIIIIIHIGHH@H@GHHHEHHHHHIIF@@F?CHIH/CGGEHH@G/FHHHHHHFHHHHE@.DCHHHE..:GHE@.8.ABEHHEAHEG.BG.9B@BG.
```

The first line is the sequence name. The second line is the DNA sequence. The third line, which always begins with a "+" can contain other optional information, but usually does not. The fourth line encodes the base quality scores in ASCII characters. 

The command below will let you inspect more part of the file.
```bash
less son.1.fq
```   

The command below will count the number of reads in the file. 
```bash
grep -c "^+" son.1.fq
```   

It's worth noting that while we are using uncompressed fastq files in this example, all of the programs we are working with will accept files compressed using gzip (generally denoted with the file extension '.gz'), and it's good practice to keep fastq files compressed so they use less storage space. 


## Assess read quality ##

In order to evaluate the general quality of reads in the file we will use the `FastQC` package.  The command can take multiple files as input and outputs a quality report in html format. 

Our script runs fastqc as follows:

```bash
fastqc -t 4 -o ../fastqc ../rawdata/*fq
```

The `*fq` indicates that fastqc should run on all files ending in "fq" in directory `rawdata/`. 

Once the files are generated you have to transfer the file to your local computer to open it and examine the results carefully. To copy the file from the Xanadu cluster please use the `transfer.cam.uchc.edu` node.  

You can use an ftp program, or the unix utility `scp` as:
```bash
scp user_name@transfer.cam.uchc.edu:/FULL_PATH_to_FILES/*.html . 
```

Again, `*html` will indicate that `scp` should copy all files ending in "html"

INSERT FIGURES HERE? OR LEAVE THEM TO WORKSHOP? 
___
scripts: 
- [scripts/Part1b_fastqc.sh](scripts/Part1b_fastqc.sh)

## Quality trim ##

Current variant callers account for uncertainties in mapping (conditional on the quality of the reference genome) and in base calling, so quality trimming is not generally necessary for this application (the worrisome sources of error in variant calling are "known unknowns", like the incompleteness of the reference genome, or systematic error arising from library prep). However, if you have a dataset plagued by adapter contamination or poor quality reads, you may want to try trimming to salvage it and/or remove some of the noise. 

`sickle` is a commonly used tool for this task. 

We could trim the reads for the son as follows:

```bash
SEQ=son
sickle pe -t sanger \
	-l 100 \
	-f ../rawdata/$SEQ.1.fq \
	-r ../rawdata/$SEQ.2.fq \
	-o ../rawdata/$SEQ.trim.1.fq \
	-p ../rawdata/$SEQ.trim.2.fq \
	-s ../rawdata/$SEQ.trim.0.fq
```

This would discard any read trimmed shorter than 100bp, and if its pair was longer than 100bp, it would be placed in the file given by `-s`, which would be read as `../rawdata/son.trim.0.fq`. 

We can run our script, which also runs `FastQC` on the trimmed data, but as we will see, it will have little impact for this particular dataset. The next steps will use the untrimmed data. 
___
scripts:	
- [scripts/Part1b2_sickle_fastqc.sh](scripts/Part1b2_sickle_fastqc.sh)    


## Align and compress ##

Now that we have QC-ed our sequence data, it's time to align it to a reference genome. For that we'll use `bwa`, one of the most widely used short-read aligners. `bwa` implements several alignment methods, but `mem` is best for our application. We previously indexed our reference genome, so we're ready to go here. 

```bash
# set a variable 'GEN' that gives the location and base name of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38
# execute bwa mem
bwa mem -t 4 -R '@RG\tID:son\tSM:son' $GEN ../rawdata/son.1.fq ../rawdata/son.2.fq -o ../align_stepwise/son.sam
```
Here `-t 4` indicates that the program should use four processors, and we feed `bwa mem` both the location of the reference genome and both files of paired end reads.   

The `-R` flag specifies the read group information. Adding read group information is critical, though it does not have to be done at this stage (picard has a tool for adding read groups to alignment files). The read group can specify the source of the reads, including the library, sequencing machine, run, lane (for more see the SAM specification) but critically, for multisample variant calling it should include the sample ID.  

When simultaneously calling variants on many samples, variant callers do not track reads by their alignment file of origin, but using read group information. If this information is absent, all reads will be treated as if they came from a single sample, and a single genotype call will result. Maintaining detailed read group information can also be helpful if some sequencing runs turn out to be problematic. In that case, even if all reads are pooled in a single bam file, there are tools you can use to filter out undesirable read groups on the fly.  

Finally, we'll compress the resulting alignment file:

```bash
samtools view -bhS ../rawdata/son.sam >../align_stepwise/son.bam
```

Because we're doing all steps individually, you may note that at this point (if we've done the quality trimming) we now have 4 copies of our sequence data. 
___
scripts:	
- [scripts/Part1c_align.sh](scripts/Part1c_align.sh)    
- [scripts/Part1d_compress.sh](scripts/Part1d_compress.sh)

## Sort reads by genome position ##

To call variants at a given position in the reference genome, we need to look at all the reads that overlap that position. In order to do this efficiently, we need to sort the reads in the alignment files by their positions in the reference genome. We'll use `picard tools` for this. For example:

```bash
IN=../align_stepwise/son.bam
OUT=../align_stepwise/son.sort.bam
java -jar $PICARD SortSam \
        INPUT=$IN \
        OUTPUT=$OUT \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True
```

scripts:	
- [scripts/Part1e_sort.sh](scripts/Part1e_sort.sh)

## Mark duplicates ##

Duplicate sequences are those which originate from the same molecule after extracting and shearing genomic DNA. There are two types: _optical_ and _polymerase chain reaction (PCR)_ duplicates. Optical duplicates are an error introduced by the sequencer. PCR duplicates are introduced by library prepartion protocols that use PCR. Duplicates cause 2 types of artifacts that mislead variant callers. 
- __First__, errors introduced by the polymerase can be propagated to multiple copies of a given fragment. Because these errors are actually part of the DNA sequence, they are likely to have high base qualities. If many sequences from the fragment containing the error are present, the variant caller can be deceived into identifying it as biological variation. 
- __Second__, when variant callers call genotypes, they assume that heterozygous sites will have equal representation of both alleles in the sequence pool (as they should for germ-line mutations). Dramatically unbalanced coverage of an allele can be a signal that variation is spurious. Because of its exponential reproduction of fragments, PCR can randomly alter allele balance, causing a variant caller to incorrectly call genotypes as homozygotes, or a whole site as invariant. 

For these reasons we need to exclude duplicate sequences from variant calling. They can be identified most easily from paired-end data as those sequences for which both reads have identical start sites. This may eliminate some sequences which are in fact derived from unique fragments in the original library, but if fragmentation is actually random, identical fragments should be rare. Once identified, duplicate sequences can be marked and ignored during variant calling (or other types of analyses) downstream. 

Here is some example code:

```bash
# assign input and output file names to variables
IN=../align_stepwise/son.sort.bam
OUT=../align_stepwise/son.mkdup.bam
# run MarkDuplicates
java -jar $PICARD MarkDuplicates \
        INPUT=$IN \
        OUTPUT=$OUT \
        METRICS_FILE=$IN.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=True
```

In this example, duplicate sequences remain in the file, but they are flagged as such. 

___

scripts:
- [scripts/Part1f_markduplicates.sh](scripts/Part1f_markduplicates.sh)

## Index alignment files ##

The last step in preparing the reads is to index the bam files. This needs to be done to enable fast access to reads overlapping a given position of the genome. Without the index, if you wanted to access reads at the beginning of chromosome 20, you'd need to read through chromosomes 1-19 until you got there. In this case, `picard` has done that for us as part of MarkDuplicates, but if you used a different approach, you'd need to use `samtools` to do this. For example:   

```bash
samtools index ../align_stepwise/*mkdup.bam
```

scripts:
- [scripts/Part1g_indexbams.sh](scripts/Part1g_indexbams.sh)

___

# Exploring SAM files #

_perhaps this should be broken out into a separate sub-section_

Now that we have processed alignment files we can learn about sam files and do some exploring. The ability to inspect and summarize various aspects of an alignment file can be especially helpful when something has gone wrong with sequencing or bioinformatic analysis. 

There are no scripts for the section below. Run the code interactively. 

Add discussion of SAM format here. 

_Useful links:_
- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

- [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)

___

We can get a lot of basic stats on the SAM file using "samtools stats":

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools stats
samtools stats -r $GEN son.mkdup.bam >son.samstat.txt
```
This file is rather messy, but we can pull out specific parts of it using grep. 

These are some basic stats about the alignment file:

```bash
grep ^SN son.samstat.txt | cut -f 2-
```

This is a histogram of per base coverage:

```bash
grep ^COV son.samstat.txt | cut -f 2-
```
And there is much more information. 

___

We can look at the per base coverage of individual regions easily using "samtools depth"

```bash
# make a list of bam files
ls *mkdup.bam >list.bam
# run samtools depth using that list
samtools depth -f list.bam -r chr20:13934265-13934558
```

We expect the parents to have 100x coverage and the son to have 50x coverage. Is there anything unusual about this region?

___

Samtools has implemented a simple alignment viewer in 'tview'. IGV is better, but tview is right in the terminal. Let's look at the same region as above. 

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools tview
samtools tview --reference $GEN dad.mkdup.bam
```

After opening tview, type 'g' then enter 'chr20:13934265' to visit that location. 
"?" will bring up a help dialog
"q" exits

___

IGV is a much better way of visualizing alignment files...
