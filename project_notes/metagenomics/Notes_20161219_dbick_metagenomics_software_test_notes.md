# Testing Metagenomics datasets with new methods
---
*12/19/2016**

These are my notes on tinkering around on previously generated metagenomics and ribotyping experiments. I will try to slowly advance towards the use of new techniques on the data so that I can garner some publications in the field.

## Table of Contents
* [The Austrian Acidosis dataset summary](#austriansummary)
* [Meta-analysis of metagenomics WGS data](#metanalysis)


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

<a name="metanalysis"></a>
## Meta-analysis of metagenomics WGS data

I think that there are a number of studies on SRA that consist of WGS experiments on cattle rumen. I am hoping to download just the cattle data, process it through [MASH](https://github.com/marbl/Mash) and then identify inter-individual variability/WGS coverage bias using the dataset.

I have selected a large group of WGS experiments from SRA, and I have removed mischaracterized datasets that point to chicken, pig or other ruminants. The list of files is in SharedFolders/metagenomics/metanalysis/sra_file_accession_list.csv

I am now going to download the list of fastqs from SRA to perform the downstream analysis.

> Fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/datasources

```bash
# sending off the list of accession entries to fastq-dump
cat accession_list.txt | xargs -I {} sbatch --nodes 1 --tasks-per-node 1 --mem 500 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/sratoolkit.2.8.1-centos_linux64/bin/fastq-dump.2 -I --split-files {}"
```

Now I need to try to generate Mash sketches. First, I'd like to see how kmer size impacts MASH sketch attributes. Let's test this out on a range of sizes.

```bash
mkdir sketch_test
#noise filter = 2, 15mer
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 15 -r -m 2 -o sketch_test/DRR017219_15 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"

# std output follows
Sketching datasources/DRR017219_1.fastq...
Sketching datasources/DRR017219_2.fastq...
Estimated genome size: 2.1573e+08
Estimated coverage:    7.202
Estimated genome size: 2.29383e+08
Estimated coverage:    7.105
Writing to sketch_test/DRR017219_15.msh...

mv slurm-653995.out ./sketch_test/DRR017219_15.out

sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 18 -r -m 2 -o sketch_test/DRR017219_18 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -r -m 2 -o sketch_test/DRR017219_21 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 24 -r -m 2 -o sketch_test/DRR017219_24 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 32 -r -m 2 -o sketch_test/DRR017219_32 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"


# In order to estimate the distance between sketches, they must be the same kmer size
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 15 -r -m 2 -o sketch_test/DRR019503_15 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 18 -r -m 2 -o sketch_test/DRR019503_18 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -r -m 2 -o sketch_test/DRR019503_21 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 24 -r -m 2 -o sketch_test/DRR019503_24 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 32 -r -m 2 -o sketch_test/DRR019503_32 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
