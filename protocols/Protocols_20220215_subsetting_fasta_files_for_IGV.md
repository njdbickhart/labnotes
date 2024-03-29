# Notes on how to subset MAG BAM files for analysis

Location:

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye


## Setup

Here is the location on Ceres:

```bash
pwd
	/lustre/project/forage_assemblies/sheep_project/complete_flye
```

And here are the relevant files:

* flye4.contigs.fasta	The HiFi assembly
* flye4.contigs.fasta.ccs.bam 	HiFi read alignments to the HiFi assembly contigs
* flye4.prodigal.shortform.tab	The ORF predictions on the HiFi assembly


Our targets are the two largest contigs in the assembly. Here are those contigs:

* contig_65409
* contig_652

```bash
grep contig_65409 flye4.contigs.fasta.fai
	contig_65409    5546585 2723763400      60      61

grep -P 'contig_652\t' flye4.contigs.fasta.fai
	contig_652      4500737 2713110351      60      61
```

## Generating the relevant files for IGV viewing

```bash
# Making things tidy by making a new directory for output
mkdir test_subset_largest


## READ ALIGNMENT FILES (BAM)
# First, let's subset the BAM file to select only the reads that map to the two contigs:
module load samtools
#                Output file name                              Input BAM                   [Regions, these can be contig names or subsets of contigs]
samtools view -o test_subset_largest/flye4.largest_two.ccs.bam flye4.contigs.fasta.ccs.bam contig_65409 contig_652


# In order to use the BAM file, we must index it
samtools index test_subset_largest/flye4.largest_two.ccs.bam

## FASTA SEQUENCE OF REFERENCE
# We can also use samtools to generate a subset of our two largest contigs for plotting
samtools faidx flye4.contigs.fasta contig_65409 contig_652 > test_subset_largest/flye4.largest_two.fasta

# This file only contains sequence from the largest two contigs

## GENE BED FILE
# This one is trickier and required allot of other programs and custom scripts to generate. We need a BED file for all of the predicted ORFs in the contigs
# We will pull them from Prodigal ORF predictions
perl -lane '%check = ("contig_65409" => 1, "contig_652" => 1); @j = split(/_/, $F[0]); $name = "$j[0]_$j[1]"; if($check{$name}){print "$name\t$F[1]\t$F[2]\t$F[0]";}' < flye4.prodigal.shortform.tab > test_subset_largest/flye4.largest_two.orfpreds.bed

# Packaging everything up for easy transfer
tar -czvf test_subset_largest.tar.gz test_subset_largest/
```

## Output files and descriptions

Here are the final files in the subset folder:

* flye4.largest_two.ccs.bam	Read alignments to only the two largest contigs
* flye4.largest_two.ccs.bam.bai	The read alignment file index (important to have in the same directory -- allows fast access to alignment data)
* flye4.largest_two.fasta	A subset fasta file containing sequence from only the two largest contigs in the assembly
* flye4.largest_two.orfpreds.bed	A bed file containing the coordinates of all of the predicted ORFs in the two contigs

## Generating Fastq files from BAM files

There are specialized tools to do this, but I always forget their arguments and I prefer to use Perl in any case here! Also, note that the alignments do not contain the usual quality scores because I aligned fasta files (not fastq!) to the contigs. In order to generate a suitable FASTQ, I need to generate the quality scores de novo.

```bash
module load samtools

samtools view test_subset_largest/flye4.largest_two.ccs.bam contig_65409 | perl -lane '$q = "I" x length($F[9]); print "\@$F[0]\n$F[9]\n+\n$q";' > test_subset_largest/contig_65409.fastq
samtools view test_subset_largest/flye4.largest_two.ccs.bam contig_652 | perl -lane '$q = "I" x length($F[9]); print "\@$F[0]\n$F[9]\n+\n$q";' > test_subset_largest/contig_652.fastq

```