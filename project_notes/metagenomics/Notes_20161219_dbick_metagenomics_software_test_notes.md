# Testing Metagenomics datasets with new methods
---
*12/19/2016**

These are my notes on tinkering around on previously generated metagenomics and ribotyping experiments. I will try to slowly advance towards the use of new techniques on the data so that I can garner some publications in the field.

## Table of Contents
* [The Austrian Acidosis dataset summary](#austriansummary)


<a name="austriansummary"></a>
## The Austrian Acidosis dataset summary

I found this dataset on [SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB9353) and thought it would be a good test case for the following reasons:

* The Austrians did not label the data well
* They tested the epimural layer
* They tested different animal timepoint stages
* They considered animals in acidosis

Here is their project summary:

> This study aimed to investigate the impact of diet induced subacute rumen acidosis (SARA) in a transient feeding model on the ruminal bacterial epithelial microbiome in dairy cattle. Eight dry dairy cows were fed high-grain diet, mainly consisting of barley grain, wheat, corn, and rapeseed meal. Baseline diet was a forage-only diet consisting of hay and grass silage Feeding model was conducted as follow: one week adaptation (10% daily from 0 to 60%, on dry matter basis), one week SARA challenge, one week only forage-feeding and two weeks SARA challenge again (this time using a two-days adaptation phase). DNA was isolated from the rumen mucosa and used for Illumina MiSeq 16S rRNA gene amplicon sequencing.

My goal is to: 

1. Perform a cursory ribotyping experiment using their data
2. Use principal component analysis to isolate putative timepoint layers

It should be a straightforward exercise, but the differentiation of the different samples will be difficult as they have two SARA challenges, and I don't know how well they performed their sample collection.

Let's collect the data from a list of accessions that I've gathered from NCBI.

> pwd: /home/dbickhart/share/metagenomics/austrian_metagenomics_test_resources

```bash
# Test download
~/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --split-files ERR920940

# There was only one fastq file... the Austrians  must have combined their reads from the MiSeq into a long, overlapping read
# Getting all the rest
for i in `cat SraAccList.txt`; do echo $i; ~/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --split-files $i; done

# Moving them all to a separate folder
mkdir fastqs
mv ./*.fastq ./fastqs/
```

