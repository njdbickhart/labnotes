# Polled TF investigation
---
*3/27/2019*

These are my notes on trying to identify transcription factor binding sites in the regions of the polled gene for future gene editing experiments.

## Table of contents


## Recap of conversation

We want to try to identify any functional element within the polled locus that may be causal for the phenotype. One possibility is to identify TF factors in the region. I will try to do a liftover from my previous analysis. Worst case -- I may need to rerun TFLOC or find human sequence conservation sites to do this.

I just reviewed the Nguyen et al. manuscript and it appears that their analysis overlapped fairly well with mine. I think that their results are going to be a good approximation of what I did albeit without all of the tuning that I would need to do on my analysis pipeline! I will use their data as a stepping off point. 

| Locus | UMD3 coords | Btau4 coords |
|:--- | :--- | :--- |
|Celtic | chr1:1705834-1706045 | chr1:1,517,223-1,517,434 |

## Preparing files

I want to gather all of the materials I need. Unfortunately, there's a paucity of chromatin accessibility data in cattle. I want to use this data to validate associations among TF's found in the Australian publication, but the only tissue that the FAANG group uploaded was for [CD-4 cells!](https://www.ebi.ac.uk/ena/data/view/ERX2628476) Oh well... maybe something sticks out? 

I'm going to download a bunch of data for analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/unmapped_scraping/side_test

```bash
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/000/ERR2611830/ERR2611830_1.fastq.gz"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/000/ERR2611830/ERR2611830_2.fastq.gz"
```

I ran some ContraTF analysis on the entire region and on the Celtic locus. I'm looking for any TF that could be in the region. The 212 bp duplicated in the Celtic locus has some multi-species conservation sites.

The mongolian locus has an enhancer site located in this region: chr1:1974983-1978220 