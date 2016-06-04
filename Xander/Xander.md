#Xander

Authored by Taylor Dunivin for EDAMAME2016 based on a previous tutorial by Qiong Wang for EDAMAME2015 

[EDAMAME-2016 wiki] (https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!
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


###2 Prepare gene reference
As Xander is a gene-targeted metagenome assembler, the first steps involve preparing a profile hidden Markov model of your gene of interest. This will ultimately be used to guide the assembly. 

The Xander assembler, which is included in the RDPTools, can be accessed and downloaded from the [RDP staff GiHhub] (https://github.com/rdpstaff/RDPTools). Since this is already included in the AMI,  we can skip this download in the tutorial and navigate to the Xander assembler directory.

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler
```

In the Xander assembler, you should see several directories including ```gene_resource```. If you navigate to this directory, you will see a list of genes that have pre-made hidden Markov models that can be used for assembly of your metagenome. We will now go through how to make your own profile hidden Markov model in case your gene of interest is not pre-processed. 

For each individual gene of interest, you should make a separate directory in ```gene_resource```. Let's try it, and make a new directory for the gene _rpoB_. 

```
mkdir rpoB
```

Now navigate to this directory, and make a new directory called ```original_data```

```
cd rpoB
mkdir originaldata
```

In this directory, we will the four required files for Xander. This step requires biological insight!
* **gene.seeds**: a small set of **full length, high quality** protein sequences in FASTA format, used to build gene.hmm, forward and reverse HMMs. 
* **gene.hmm**: this is the HMM built from gene.seeds using original HMMER3. This is used to build for_enone.hmm and align contigs after assembly. 
* **framebot.fa**: a large near full length known protein set for identifying start kmers and FrameBot nearest matching. More diversity is better, more sequences means more starting points (more computational time) but less susceptible to noise than model creation. Prefer near full-length and well-annotated sequences. Filter with Minimum HMM Coverage at least 80 (%).
* **nucl.fa**: a large near full length known set used by UCHIME chimera check.

These files can be downloaded from the RDP's FunGene database. In your web browser, open the RDP's FunGene database [http://fungene.cme.msu.edu] (http://fungene.cme.msu.edu). Here you will find a list of gene families based on searches of the NCBI non-redundant protein database using "training sequences" or **gene.seeds**. Select _rpoB_ in the Phylogenetic markers heading. 

To download the gene.seeds and gene.hmm files, click on the links at the top left of a page. Note this may not work with Safari.
(add image of fungene with links boxed in red)

To obtain the framebot.fa and nucl.fa files, click on the show/hide filter options link at the top right of the page. Here you can limit the the score, minimum size aa, and minimum HMM coverage. The score will vary by gene. You can leave this area blank, or look though the sequences to see if there is a large drop off in score at a specific number. The minimum size aa will depend on the actual protein size. The minimum HMM coverage must be greater than 80%. What percent cutoff you choose will require biological insight. Lower % HMM coverage will increase diversity but may lower the quality of your search. Once you have set your parameters, you will select sequences. Then click begin analysis. This will take you to a page where you can download your framebot.fa (protein) and nucl.fa (nucleotide) sequences. Then change the names to framebot.fa and nucl.fa respectively. Move all of these files to ```originaldata```. 

From these files, we need to make three new files
* **for_enone.hmm** and **rev_enone.hmm** for the forward and reverse HMMs respectively. This is used to assemble gene contigs.
* **ref_aligned.faa** file containing a set of protein reference sequences aligned with for_enone.hmm. This is used to identify starting kmers. Need to manually examine the alignment using Jalview or alignment viewing tool to spot any badly aligned sequences. If found, it is likely there are no sequences in gene.seeds close these sequences. You need to valide these problem sequences to see they are from the gene you are interested, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps.

This can be accomplished using the ```prepare_gene_ref.sh``` script, but first, we need to edit this file. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/bin
nano prepare_gene_ref.sh
```

We need to change four parts of this file. 
* line 4: change gene to ```rplB```
* line 9: change the jar directory to ```/home/ubuntu/tools/RDPTools```
* line 10: change the reference directory to ```/home/ubuntu/tools/RDPTools/Xander_assembler```
* line 14: change the hmmer xanderpatch location to ```hmmer_xanderpatch=/home/ubuntu/tools/third_party_tools/hmmer-3.0-xanderpatch```

Now save the changes you've made to the shell script. 

Excellent! You now have a script ready to prepare the gene reference for Xander. Let's run it. All you need to do is excecute this script and specify your gene of interest. In our case, this is _rpoB_. 

```
./prepare_gene_ref.sh rpoB
```

Check and see that all of your output files are in your gene directory. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/rpoB
```

Now that we know the script ran successfully, we need to check ```ref_aligned.faa``` for any poorly aligned sequences. 
(should I include this part?)



###3 Set up environmental variables 

Since you may want to run Xander multiple times, it can be useful to make a directory for each project. In our case, we will use the pre-existing directory ```testdata``` as this contains practice data provided by the creators of Xander. 

We will need to edit two shell scripts to run Xander. You may want to keep the original files for future reference, so we will copy them to into ```testdata```. 

```
scp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/tools/RDPTools/Xander_assembler/testdata
scp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/tools/RDPTools/Xander_assembler/testdata
```

Now that we have these scripts copied, we can edit them. We first need to change the base directory in ```run_xander_skel.sh```. 

```
nano run_xander_skel.sh
```

In line 4, change the base directory to reflect our base directory ```/home/ubuntu/tools/RDPTools/Xander_assembler/bin```

Awesome! That's all we need to change here, so let's save it and move on. 

Now we need to change ```xander_setenv.sh``` to reflect our directories and gene of interest. This is also where we can adjust Xander parameters. 

```nano run_xander_setenv.sh```

* SEQFILE= ```/home/ubuntu/tools/RDPTools/Xander_assembler/testdata/???.fastq```
* WORKDIR= ```/home/ubuntu/tools/RDPTools/Xander_assembler/testdata```
* REFDIR= ```/home/ubuntu/tools/RDPTools/Xander_assembler```
* JARDIR= ```/home/ubuntu/tools/RDPTools```







