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

#### Scoring existing metrics

They gave me several sites to check. Let's migrate the data and generate the calls. Hopefully I can add the data via a tab delimited file and use scoring to estimate the likelihood of the calls.

> Ceres: /lustre/project/forage_assemblies/r(tab)

```bash
module load java/1.8.0_121

# Indexing
sbatch --nodes=1 --mem=50000 --ntasks-per-node=10 -p priority -q msn --wrap="java -Xmx50g -jar FlashFry-assembly-1.10.jar index --tmpLocation ./tmp --database ARS-UCD1.2_Btau5.0.1Y_cas9_database --reference /lustre/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa --enzyme spcas9ngg"

mkdir cpf1tmp; sbatch --nodes=1 --mem=50000 --ntasks-per-node=10 -p priority -q msn --wrap="java -Xmx50g -jar FlashFry-assembly-1.10.jar index --tmpLocation ./cpf1tmp --database ARS-UCD1.2_Btau5.0.1Y_cpf1_database --reference /lustre/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa --enzyme cpf1"

# Now running the discovery and scoring
java -Xmx14g -jar FlashFry-assembly-1.10.jar discover --database ARS-UCD1.2_Btau5.0.1Y_cas9_database --fasta s_hdr.fa --output s_hdr.output

# It "worked" but didn't actually get the sites I wanted. Let's try the exact sequence from the probe taken from UCSC genome browser
java -Xmx14g -jar FlashFry-assembly-1.10.jar discover --database ARS-UCD1.2_Btau5.0.1Y_cas9_database --fasta s_site.fa --output s_site.output

# It worked! I need to include the CAS9 binding site though to get it to run. Now to score
java -Xmx4g -jar FlashFry-assembly-1.10.jar score --input s_site.output --output s_site.output.scored --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot --database ARS-UCD1.2_Btau5.0.1Y_cas9_database


# Now to queue up the rest
for i in hp1_site hp2_site; do echo $i; java -Xmx14g -jar FlashFry-assembly-1.10.jar discover --database ARS-UCD1.2_Btau5.0.1Y_cpf1_database --fasta ${i}.fa --output ${i}.output; java -Xmx14g -jar FlashFry-assembly-1.10.jar score --input ${i}.output --output ${i}.output.scored --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot --database ARS-UCD1.2_Btau5.0.1Y_cpf1_database; done
```

Now let's queue this up and run it on all assemblies. I copied over a select listing of data points from the excel sheet that was sent to me and then used a script to format it all into a json config file.

```bash
module load java/1.8.0_121 samtools blat/36

# The formatting
perl -ne 'chomp; @s = split(/\t/); $s[0] =~ s/\s+/_/g; $s[0] =~ s/\./_/g; $s[4] =~ s/\%//g; print "  \"$s[0]\" : {\n"; print "   \"enzyme\" : \"$s[2]\",\n"; print "   \"seq\" : \"$s[1]\",\n"; print "   \"edit\" : \"$s[3]\",\n"; print "   \"editperc\" : \"$s[4]\"\n  },\n";' < template.tab >> config.json

# Testing
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake -np -s snakeFile

# OK, looks good. Let's test this bugger out
mkdir logs; sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake -s snakeFile --cluster "sbatch --nodes=1 --ntasks-per-node=3 --mem=15000 --partition=priority -q msn -o logs/{rule}.{wildcards.names}.out" -p --jobs 250

# Alright, there were issues with flashfry finding a PAM for a sequence with an "N" base, so I requested the actual PAM sequence from the cooperator.
# I think that I can work with this now
perl -ne 'chomp; @s = split(/\t/); $s[0] =~ s/\s+/_/g; $s[0] =~ s/\./_/g; $s[4] =~ s/\%//g; print "  \"$s[0]\" : {\n"; print "   \"pam\" : \"$s[2]\",\n"; print "   \"enzyme\" : \"$s[3]\",\n"; print "   \"seq\" : \"$s[1]\",\n"; print "   \"edit\" : \"$s[4]\",\n"; print "   \"editperc\" : \"$s[5]\"\n  },\n";' < template2.tab >> config.json

# I copied over the snakefile and renamed it. Let's try this again!
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake -s fryPamCheck --cluster "sbatch --nodes=1 --ntasks-per-node=3 --mem=15000 --partition=priority -q msn -o logs/{rule}.{wildcards.names}.out" -p --jobs 250
```

## Primer design pipeline test

Testing this on taurine sequence with 10 added designs to test

> Ceres: /lustre/project/forage_assemblies/recombinetics/fryTaurus

```bash
module load exonerate/2.2.0 java/1.8.0_121 samtools blat/36

mkdir logs
sbatch --nodes=1 --mem=5000 -c 2 -p priority -q msn  snakemake -s fryPamCheck --cluster "sbatch --nodes=1 -c 3 --mem=15000 --partition=priority -q msn  -o logs/{rule}.{wildcards.names}.out" -p --jobs 250

# Printing the DAG for showing to the group
snakemake --forceall --dag -s fryPamCheck > fryPamCheck.dag
```

## Diagnosing Eui-soo's offtarget site analysis

I suspect that Eui-soo's offtargets are in highly repetitive kmers. I will check this using Jellyfish k-count queries.

> Ceres: /lustre/project/cattle_genome_assemblies/dominette/ARS-UCD1.2_Btau5.0.1Y

```bash
module load jellyfish2/2.2.9

sbatch --nodes=1 --mem=200000 --ntasks-per-node=30 -p priority -q msn --wrap="jellyfish count -m 23 -C -s 500M -o ARS-UCD1.2_Btau5.0.1Y.23mer.jf ARS-UCD1.2_Btau5.0.1Y.fa"

jellyfish query ARS-UCD1.2_Btau5.0.1Y.23mer.jf -s /lustre/project/forage_assemblies/recombinetics/fryTaurus/test_offt_sites2.fasta
# Nothing

jellyfish query ARS-UCD1.2_Btau5.0.1Y.23mer.jf -s /lustre/project/forage_assemblies/recombinetics/fryTaurus/test_offt_sites.fasta
# Everything had at least one kmer
```