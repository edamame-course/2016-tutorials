#Xander

Authored by Taylor Dunivin for EDAMAME2016 based on a previous tutorial by Qiong Wang for EDAMAME2015 

[EDAMAME-2016 wiki] (https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **gene-targeted metagenome assembly**

##Learning Objectives
* Prepare gene references for assembly
* Understand Xander assembly parameters 
* Assemble megatenome based on specified gene
* Examine the abundance and diversity of genes of interest from metagenome data

---

###Citation
Wang, Q., J. A. Fish, M. Gilman, Y. Sun, C. T. Brown, J. M. Tiedje and J. R. Cole. Xander: Employing a Novel Method for Efficient Gene-Targeted Metagenomic Assembly. Microbiome.2015, 3:32. DOI: [10.1186/s40168-015-0093-6] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6). 

###Required tools
* [RDPTools] (https://github.com/rdpstaff/RDPTools)
* Python 2.7+
* Java 2.6+
* [HMMER 3.1] (http://hmmer.org) (If using HMMER 3.0 add --allcol to bin/run_xander_skel.sh )
* UCHIME 

###1 Connect to Xander AMI
For this tutorial, we will use an existing AMI that contains all the necessary tools to run Xander. Search for the AMI name "RDP-Edamame-2015" or ID "ami-e973b782"

###2 
