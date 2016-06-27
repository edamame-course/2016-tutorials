#Mapping and counting of metagenome data
Authored by Jin Choi

## Overarching Goal
* This tutorial will contribute towards an understanding of **metagenome data**

##Learning Objectives
* Understanding mapping file
* Use mapping program
* Counting

## updating the operating system
```
apt-get update
apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat
```

## Install software
```
cd /root
wget -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download

tar xvfj bwa-0.7.10.tar.bz2
cd bwa-0.7.10
make

cp bwa /usr/local/bin
```

install samtools
```
apt-get -y install samtools
```

## Download data
```
cd /mnt

curl -O http://athyra.idyll.org/~t/REL606.fa.gz
gunzip REL606.fa.gz

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098038/SRR098038.fastq.gz
```

## Do the mapping
Now let’s map all of the reads to the reference. Start by indexing the reference genome:
```
cd /mnt

bwa index REL606.fa

```
Now, do the mapping of the raw reads to the reference genome:
```
bwa aln REL606.fa SRR098038.fastq.gz  > SRR098038.sai
```
Make a SAM file (this would be done with ‘sampe’ if these were paired-end reads):
```
bwa samse REL606.fa SRR098038.sai SRR098038.fastq.gz > SRR098038.sam
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

```
samtools tview SRR098038.sorted.bam REL606.fa
```

## counting alignments

```
samtools view -c -f 4 SRR098038.bam
```

## Calling SNPs
```
samtools mpileup -uD -f REL606.fa SRR098038.sorted.bam | bcftools view -bvcg - > SRR098038.raw.bcf
```

## Make count table

