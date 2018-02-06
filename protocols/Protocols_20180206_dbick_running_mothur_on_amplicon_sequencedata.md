# Mothur analysis of 16S amplicons
---
*2/6/2018*

These are my notes on how to perform a Mothur analysis on 16S amplicons using Madison's best practices workbook.



## Preparing for the analysis

```bash
mkdir buccal_sop
cp FASTQ_Generation_2017_12_17_10_37_25Z_67404143/*.fastq ./buccal_sop/

wget https://www.mothur.org/w/images/3/32/Silva.nr_v132.tgz
wget http://www.mothur.org/w/images/1/19/Gg_13_8_99.refalign.tgz
wget http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz

tar -xvf Silva.nr_v132.tgz
tar -xvf Gg_13_8_99.refalign.tgz; tar -xvf Gg_13_8_99.taxonomy.tgz
```

## Filtering and preparing sequence files for assembly/redundancy

> /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/buccal_samples/pilot_project

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/mothur/mothur

make.file(inputdir=buccal_sop, type=fastq)
make.contigs(file=buccal_sop/stability.files, processors=60)
