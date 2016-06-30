#Mapping and counting of metagenome data
Authored by Jin Choi for EMAMAME2016 based on a previous tutorial by Adina Howe

## Overarching Goal
* This tutorial will contribute towards an understanding of **metagenome data**

##Learning Objectives
* Understanding mapping file
* Use mapping program(Bowtie2, samtools, bcftools)
* Counting

## Mapping and counting useful when,,
* If you want to find point mutation when you sequence pure culture
* If you assemble your multiple sample by yourself when you sequence metagenomic sample
* If you want to generate gene expression table when you sequence transcriptomic/metatranscriptomic sample

## updating the operating system
```
apt-get update
apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat
```

## Install software
```
cd 
wget wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
mv bowtie2-2.2.9 BT2_HOME
```
Set up path(you need to set up path again if you are logged in later):
```
PATH=$PATH:~/BT2_HOME
export PATH
```

install samtools
```
apt-get -y install samtools
```

## Download data
First, download referene sequence
```
cd ~/
mkdir mapping_counting
cd mapping_counting
curl -O http://athyra.idyll.org/~t/REL606.fa.gz
gunzip REL606.fa.gz
```
Then, download sequencing file
```
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098038/SRR098038.fastq.gz
```

## Do the mapping
Now let’s map all of the reads to the reference. Start by indexing the reference genome:
```
cd ~/mapping_counting

bowtie2-build REL606.fa REL_reference

```
Now, do the mapping of the raw reads to the reference genome (this would be done with -1 and -2 if these were paired-end reads):
```
bowtie2 -x REL_reference -U SRR098038.fastq.gz -S SRR098038.sam
```

This file contains all of the information about where each read hits on the reference.

Next, index the reference genome with samtools:

```
samtools faidx REL606.fa
```

Convert the SAM into a BAM file:
```

samtools import REL606.fa.fai SRR098038.sam SRR098038.bam
```

Sort the BAM file:
```
samtools sort SRR098038.bam SRR098038.sorted
```

And index the sorted BAM file:
```
samtools index SRR098038.sorted.bam

```

## Visualizing alignment
At this point you can visualize with samtools tview or Tablet.

‘samtools tview’ is a text interface that you use from the command line; run it like so:
```
samtools tview SRR098038.sorted.bam REL606.fa
```
The ‘.’s are places where the reads align perfectly in the forward direction, and the ‘,’s are places where the reads align perfectly in the reverse direction. Mismatches are indicated as A, T, C, G, etc.

You can scroll around using left and right arrows; to go to a specific coordinate, use ‘g’ and then type in the contig name and the position. For example, type ‘g’ and then ‘rel606:553093<ENTER>’ to go to position 553093 in the BAM file.

Use ‘q’ to quit.

For the Tablet viewer, click on the link and get it installed on your local computer. Then, start it up as an application. To open your alignments in Tablet, you’ll need three files on your local computer: REL606.fa, SRR098042.sorted.bam, and SRR098042.sorted.bam.bai. You can copy them over using Dropbox, for example.
## counting alignments
This command:
```
samtools view -c -f 4 SRR098038.bam
```
will count how many reads DID NOT align to the reference (207029).

This command:

```
samtools view -c -F 4 SRR098038.bam
```
will count how many reads DID align to the reference (6839602).

And this command:

```
gunzip -c SRR098038.fastq.gz | wc
```

will tell you how many lines there are in the FASTQ file (28186524). Reminder: there are four lines for each sequence.
## Calling SNPs
You can use samtools to call SNPs like so:
```
samtools mpileup -uD -f REL606.fa SRR098038.sorted.bam | bcftools view -bvcg - > SRR098038.raw.bcf
```
(See the ‘mpileup’ docs here.)

Now convert the BCF into VCF:
```
bcftools view SRR098038.raw.bcf > SRR098038.vcf
```
You can check out the VCF file by using ‘tail’ to look at the bottom:
```
tail *.vcf
```
Each variant call line consists of the chromosome name (for E. coli REL606, there’s only one chromosome - rel606); the position within the reference; an ID (here always ‘.’); the reference call; the variant call; and a bunch of additional information about

Again, you can use ‘samtools tview’ and then type (for example) ‘g’ ‘rel606:4616538’ to go visit one of the positions. The format for the address to go to with ‘g’ is ‘chr:position’.

You can read more about the VCF file format here.

If you are interested in mapping metagenome file or transcriptome file, visit [here](https://github.com/metajinomics/tutorials_en/blob/gh-pages/metagenome/mapping_counting.md) for more tutorials