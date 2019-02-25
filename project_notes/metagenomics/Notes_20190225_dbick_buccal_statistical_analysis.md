# Buccal data analysis
---
*2/25/2019*

These are my notes on the statistical analysis of the Buccal swab OTU data with the overall goal of identifying which OTU are most likely part of the "oral-only" microbial community. 

## Table of Contents

## Preparing the data for analysis

I ran into some trouble getting Mothur to consistently analyze the samples I wanted to process. The problem was that the distance matrix I tried to calculate was absolutely huge! I am going to automate this via a script I've written and then use a few metrics to reduce the number of comparisons. 

I also want to include all of the samples taken thus far. The more the better. I want to include rumen contents and buccal samples to better train my models.

First, let's make sure all of the data is in the right place and is properly accounted for:

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/buccal

```bash
mkdir combined_fastqs
mv buccal_samples/evergreen_view/combined_fastq/*.fastq ./combined_fastqs/
mv buccal_samples/amelie_data/4_26_18/*.fastq ./combined_fastqs/
mv buccal_samples/pilot_project/buccal_sop/*.fastq ./combined_fastqs/

ls combined_fastqs/ | wc -l
248

mkdir tax_files
cp buccal_samples/pilot_project/silva.v132.v4.align ./tax_files/
cp buccal_samples/pilot_project/silva.nr_v132.tax ./tax_files/
cp buccal_samples/pilot_project/gg_13_8_99.fasta ./tax_files/
cp buccal_samples/pilot_project/gg_13_8_99.gg.tax ./tax_files/

# Now to make sure all of the fastq files are represented in the input file before queueing this up
perl -ne '$_ =~ s/4_26_18\///g; print $_;' < buccal_samples/amelie_data/true_listings.files > combined_fastqs/combined_fq_total.files
cat buccal_samples/evergreen_view/combined_fastq/db1.files >> combined_fastqs/combined_fq_total.files
perl -ne '$_ =~ s/FASTQ_Generation_2017-12-17_10_37_25Z-67404143\///g; print $_;' < buccal_samples/pilot_project/FASTQ_Generation_2017_12_17_10_37_25Z_67404143/pilot.true.files >> combined_fastqs/combined_fq_total.files

# Now to queue it up. The mothur batch file is: mothur_script.batch; The mothur shell wrapper is mothur_combined_run.sh

# Argh! The "input" directory setting appears to be useless! I need to put paths for all of the fastq files in the "files" file
perl -lane 'print "$F[0]\tcombined_fastqs/$F[1]\tcombined_fastqs/$F[2]";' < combined_fq_total.files > temp
mv temp combined_fq_total.files

sbatch mother_combined_run.sh
```