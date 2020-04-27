# Project notes for CRISPR gRNA design
---
*3/25/2020*

These are my notes for the project.

## Table of Contents

## Design tests

### Environmental setup

I'm going to be testing out [FlashFry](https://github.com/mckennalab/FlashFry) as the hard point program for use in this test. I'll be designing a wrapper for ease of use so that it can be used by other members of the team. Eventually, I could fork the repo and add on a rudimentary GUI. 

> Ceres: /project/forage_assemblies/r(tab)

```bash
module load java/1.8.0_121

wget https://github.com/mckennalab/FlashFry/releases/download/1.10/FlashFry-assembly-1.10.jar

# Testing on human chromosome 22 as per the vignette
wget https://raw.githubusercontent.com/aaronmck/FlashFry/master/test_data/quickstart_data.tar.gz
tar -xf quickstart_data.tar.gz

# indexing
mkdir tmp
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
 index \
 --tmpLocation ./tmp \
 --database chr22_cas9ngg_database \
 --reference chr22.fa.gz \
 --enzyme spcas9ngg

# Discovery
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
 discover \
 --database chr22_cas9ngg_database \
 --fasta EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta \
 --output EMX1.output

# Scoring
java -Xmx4g -jar FlashFry-assembly-1.10.jar \
  score \
  --input EMX1.output \
  --output EMX1.output.scored \
  --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
  --database chr22_cas9ngg_database
```