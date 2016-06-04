#Xander

Authored by Taylor Dunivin for EDAMAME2016 based on a previous tutorial by Qiong Wang for EDAMAME2015 

[EDAMAME-2016 wiki] (https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

##Learning Objectives
* Prepare gene references for assembly
* Understand Xander assembly parameters 
* Assemble megatenome based on specified gene
* Examine the abundance and diversity of genes of interest from metagenome data

---

###Citations
Wang, Q., J. A. Fish, M. Gilman, Y. Sun, C. T. Brown, J. M. Tiedje and J. R. Cole. 2015. Xander: Employing a Novel Method for Efficient Gene-Targeted Metagenomic Assembly. Microbiome. 3:32. DOI: [10.1186/s40168-015-0093-6] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6). 

Fish, J. A., B. Chai, Q. Wang, Y. Sun, C. T. Brown, J. M. Tiedje, and J. R. Cole. 2013. FunGene: the Functional Gene Pipeline and Repository. Front. Microbiol. 4: 291. DOI: [10.3389/fmicb.2013.00291] (http://www.ncbi.nlm.nih.gov/pubmed/24101916).

###Required tools
* [RDPTools] (https://github.com/rdpstaff/RDPTools)
* Python 2.7+
* Java 2.6+
* [HMMER 3.1] (http://hmmer.org) (If using HMMER 3.0 add --allcol to bin/run_xander_skel.sh )
* UCHIME 

###1 Connect to Xander AMI
For this tutorial, we will use an existing AMI that contains all the necessary tools to run Xander. Search for the AMI name "RDP-Edamame-2015" or ID "ami-e973b782"

Click launch.


###2 Prepare gene reference
As Xander is a gene-targeted metagenome assembler, the first steps involve preparing a profile hidden Markov model of your gene of interest. This will ultimately be used to guide the assembly. 

The Xander assembler, which is included in the RDPTools, can be accessed and downloaded from the [RDP staff GiHhub] (https://github.com/rdpstaff/RDPTools). Since this is already included in the AMI,  we can skip this download in the tutorial and navigate to the Xander assembler directory.

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler
```

In the Xander assembler, you should see several directories including ```gene_resource```. If you navigate to this directory, you will see a list of genes that have pre-made hidden Markov models that can be used for assembly of your metagenome. We will now go through how to make your own profile hidden Markov model in case your gene of interest is not pre-processed. 

```gene_resource``` contains reference sequence files and models. The analysis pipeline is preconfigured with the _rplB_ phylogenetic marker gene, and nitrogen cycling genes including _nirK_, _nirS_, _nifH_, _nosZ_ and _amoA_.

In this tutorial, we will use the preexisting _rplB_ files and models, but we'll briefly go over how to set up this directory. 

For each individual gene of interest, you should make a separate directory in ```gene_resource``` with a subdirectory ```originaldata```. 

```originaldata``` will hold the four required files for Xander. This step requires biological insight!
* **gene.seeds**: a small set of **full length, high quality** protein sequences in FASTA format, used to build gene.hmm, forward and reverse HMMs. Can be downloaded from RDP's FunGene database. 
* **gene.hmm**: this is the HMM built from gene.seeds using original HMMER3. This is used to build for_enone.hmm and align contigs after assembly. Can be downloaded from RDP's FunGene database. 
* **framebot.fa**: a large near full length known protein set for identifying start kmers and FrameBot nearest matching. More diversity is better, more sequences means more starting points (more computational time) but less susceptible to noise than model creation. Prefer near full-length and well-annotated sequences. Filter with Minimum HMM Coverage at least 80 (%).
* **nucl.fa**: a large near full length known set used by UCHIME chimera check.

These files can be downloaded from the RDP's FunGene database. In your web browser, open the RDP's FunGene database [http://fungene.cme.msu.edu] (http://fungene.cme.msu.edu). Here you will find a list of gene families based on searches of the NCBI non-redundant protein database using "training sequences" or **gene.seeds**. Select _rplB_ in the Phylogenetic markers heading. 

To download the gene.seeds and gene.hmm files, click on the links at the top left of a page. Note this may not work with Safari.
(add image of fungene with links boxed in red)

To obtain the framebot.fa and nucl.fa files, click on the show/hide filter options link at the top right of the page. Here you can limit the the score, minimum size aa, and minimum HMM coverage. 
* The score will vary by gene. You may choose to leave this area blank. You can also manually look though the sequences to see if there is a large drop off in score at a specific number and use this number as a cutoff. 
* The minimum size aa will depend on the actual protein size. You'll want to find a balance between near-full length sequences and diversity. 
* The minimum HMM coverage must be greater than 80%. What percent cutoff you choose will require biological insight. Lower % HMM coverage will increase diversity but may lower the quality of your search. Once you have set your parameters, click filter. Then you will need to select sequences. It is not recommended to blindly select all sequences. Instead, manually go through and make sure no oddballs are included. For reference, the Xander paper used over 700 near full-length sequences for their .fa files. Once your sequences are selected, click begin analysis at the top right of the page. This will take you to a page where you can download your framebot.fa (protein) and nucl.fa (nucleotide) sequences. Then change the names to framebot.fa and nucl.fa respectively. Move all of these files to ```originaldata```. 

From these files, we need to make three new files
* **for_enone.hmm** and **rev_enone.hmm** for the forward and reverse HMMs respectively. This is used to assemble gene contigs.
* **ref_aligned.faa** file containing a set of protein reference sequences aligned with for_enone.hmm. This is used to identify starting kmers. Need to manually examine the alignment using Jalview or alignment viewing tool to spot any badly aligned sequences. If found, it is likely there are no sequences in gene.seeds close these sequences. You need to valide these problem sequences to see they are from the gene you are interested, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps.

This can be accomplished using the ```prepare_gene_ref.sh``` script, but first, we need to edit this file. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/bin
nano prepare_gene_ref.sh
```

We need to change this file to reflect our gene of interest and our paths.
* line 4: change gene to ```rplB```
* line 9: change the jar directory to ```/home/ubuntu/tools/RDPTools```
* line 10: change the reference directory to ```/home/ubuntu/tools/RDPTools/Xander_assembler```
* line 14: change the hmmer xanderpatch location to ```hmmer_xanderpatch=/home/ubuntu/tools/third_party_tools/hmmer-3.0-xanderpatch```

Now save the changes you've made to the shell script. 

Excellent! You now have a script ready to prepare the gene reference for Xander. Let's run it. All you need to do is excecute this script and specify your gene of interest. In our case, this is _rplB_. 

```
./prepare_gene_ref.sh rplB
```

Check and see that all of your output files are in your gene directory. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/rplB
```

You should now see three new files: **ref_aligned.faa**, **for_enone.hmm**, and **rev_enone.hmm**

Remember that when running Xander on your own gene of interest, you will need to manually check ```ref_aligned.faa``` for any poorly aligned sequences. This can be done using Jalview or another alignment viewing tool. 

Now that we have our gene reference files, we can begin to set up the assembler. 

###3 Set up metagenomic data

Since you may want to run Xander multiple times, it can be useful to make a directory for each project that includes the gene and dataset used. We will do this now. 

In our case, we will use data provided by the creators of Xander. This demo_reads file here contains a subset of reads from one of seven corn rhizosphere replicates used in the original Xander publication. This subset is enriched in reads matching rplB, nirK and nifH genes. **These paired-end reads have already been quality trimmed and merged using RDP's read assembler**.

Let's navigate to the ```Xander_assembler``` directory, make a new directory in Xander assembler, and then move the practice data there. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler
mkdir rplB_demo
scp /home/ubuntu/demo/data/demo_reads.fa /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
```

Now we're ready to adjust paramters for analysis!

##4 Set up environmental variables

We will need to edit one shell script to run Xander. You may want to keep the original file for future reference, so we will copy it to into ```rplB_demo``` instead of editing them directly. 

```
scp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
scp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
```

We need to change ```xander_setenv.sh``` to reflect our directories and gene of interest. This is also where we can adjust Xander parameters. First let's think about our parameter options. 

####How to choose the FILTER_SIZE for your dataset?
For count 1 bloom, the size of the bloom filter is approximately 2^(FILTER_SIZE-3) bytes. The memory needed for Java is less than double the size of the bloom filter. Increase the FILTER_SIZE if the predicted false positive rate (in output file *_bloom_stat.txt) is greater than 1%. Based on our experience with soil metagenome data, FILTER_SIZE 32 (1 GB memory) for data file size of 2GB, 35 (8 GB) for file size of 6GB, 38 (64 GB ) for file size of 70GB, 40 (256 GB) for file size of 350GB were appropriate. For count 2 bloom filter, double the memory sizes.

####Analysis Parameters
* SEQFILE -- Absolute path to the sequence files. Can use wildcards to point to multiple files (fasta, fataq or gz format)
genes -- Genes to assemble (supported out of the box: rplB, nirK, nirS, nifH, nosZ, amoA)
* SAMPLE_SHORTNAME -- a short name for your sample, prefix of contig IDs (needed when pool contigs from multiple samples)

####DBG Parameters
* MAX JVM HEAP -- Maximum amount of memory DBG processes can use (must be larger than FILTER_SIZE below)
* K_SIZE -- K-mer size to assemble at, must be divisible by 3 (recommend 45, maximum 63)
* FILTER SIZE -- size of the bloom filter, 2**FILTER_SIZE, 38 = 32 GB, 37 = 16 GB, 36 = 8 GB, 35 = 4 GB, increase 
* FILTER_SIZE if the bloom filter predicted false positive rate is greater than 1%
* MIN_COUNT=1 -- minimum kmer occurrence in SEQFILE to be included in the final bloom filter

####Contig Search Parameters
* PRUNE=20 -- prune the search if the score does not improve after n_nodes (default 20, set to 0 to disable pruning)
* PATHS=1 -- number of paths to search for each starting kmer, default 1 returns the shortest path
* LIMIT IN SECS=100 -- number of seconds a search allowed for each kmer, recommend 100 secs for 1 shortest path, need to increase if PATHS is large

####Contig Merge Parameters
* MIN BITS=50 --mimimum assembled contigs bit score
* MIN LENGTH=150 -- minimum assembled protein contigs

Now that we understand the paramters a bit more, lets edit our environment. 

```nano run_xander_setenv.sh```

Directories must match the absolute path that we are using. 

```
SEQFILE=/home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo/demo_reads.fa
WORKDIR=/home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign
```

Now we want to select our kmer size. For this experiment, we will use ```45```, but numbers divisible by three between 45 and 63 are permissible. Higher numbers will yield more stringent results. 

Next we select the filter size. This dataset is relatively small (~2MB), so we will specify ```32```. 

We also need to change the security of this file so that we can excecute it. 

```
chmod 777 run_xander_setenv.sh
```

Now we are ready to roll!

###4 Run Xander
To run Xander, we use one simple command. Note that if you want to run multiple genes at once, simply say all of them in the command instead of one.

```
./run_xander_skel.sh xander_setenv.sh “build find search” “rplB”
```
