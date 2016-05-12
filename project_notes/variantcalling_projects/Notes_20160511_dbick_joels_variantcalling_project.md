# Identifying variants within Joel's bulls
---
*5/11/2016*

These are my notes on the alignment and variant calling performed on Joel's bulls as part of our BARD funded grant.

## Table of Contents
* [Organizing the data](#organizing)
* [Generating 1000 bulls SNP and INDEL annotations](#onethousand)
	* [Output tab file columns](#outheads)

<a name="organizing"></a>
## Organizing the data

Several of Joel's bulls have been sequenced already and I would like to perform the analysis right away to give him some data to work with. I will make a list of all previously sequenced animals that we have and then run my analysis pipeline when ready.

#### Here is the pipeline:
* SNP and INDEL calling with Samtools
* RD CNV calling with JaRMS
* PE+SR CNV calling with RAPTR-SV
* Variant annotation with my java program and SNPEff

#### Here are the animals and their attributed contributors
* Genome Canada
	* HOLCANM000005279989.bam      
	* HOLDEUM000000253642.bam      
	* HOLUSAM000002265005.bam      
	* HOLUSAM000017349617.bam
	* HOLCANM000006026421.bam      
	* HOLGBRM000000598172.bam      
	* HOLUSAM000002297473.bam      
	* HOLUSAM000123066734.bam
	* HOLCANM000100745543.bam      
	* HOLITAM006001001962.bam      
	* HOLUSAM000017129288.bam      
	* HOLUSAM000132973942.bam
* Our lab
	* HODEU000000253642
	* HOGBR000000598172
	* HOUSA000002040728
	* HOUSA000002290977
	* HOUSA000017349617
	* HOUSA000122358313
* 1000 bulls
	* HOUSA000129800008
	* HOITA006001001962
	* HOUSA000131823833
	* HOUSA000130588960
	* HOUSA000123066734
	* HOUSA000001697572
	* HOCAN000010705608
	* HOUSA000017129288
	* HOAUS000H00930377
	* HOUSA000002297473
	* HOCAN000006820564
	* HOUSA000060372887
	* HOUSA000002103297
	* HOCAN000006947936
	* HOCAN000006026421
	* HOCAN000005279989
	* HOUSA000002205082
	* HODEU000578448776
	* HOUSA000002147486
	* HOUSA000002183007
	* HOCAN000006961162

That makes **39** animals already available for some form of processing. I will explain the caveats about using the 1000 bulls data to Joel later. Here are a list of these animals on the 3850 for future reference.

> 3850: /seq1/bickhart/side_projects/joels_bulls

* **final_file_joel_bulls.list** -> Contains bam files that we currently maintain from Genome Canada and our project
* **1000_bulls_sequenced_joels_bulls.list** -> Contains the names (reformatted) of bulls that were in the 1000 bulls data files but we do not have the bams for.

<a name="onethousand"></a>
## Generating 1000 bulls SNP and INDEL annotations

I am going to create tabular excel output for Joel by splitting the 1000 bull genomes files into the tabular output that I've generated previously for Tad. Joel would like the entire chromosome from the 1kbulls VCF file. While this makes things easier for me, it may take longer to process. I'll get started immediately.

I just rewrote my [vcf subsectioning](https://github.com/njdbickhart/perl_toolchain/blob/master/vcf_utils/filterAndSubsectionVCFfile.pl) script to handle the processing of the vcf.

> 3850: /seq1/1kbulls_annotatedvcf

```bash
bgzip Chr5-Beagle-Run5.eff.vcf
bcftools index Chr5-Beagle-Run5.eff.vcf.gz

perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f Chr5-Beagle-Run5.eff.vcf.gz -o Chr5_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list

```

Here are the output data columns:

<a name="outheads"></a>
#### Output tab file columns

1. Chromosome
2. Base pair position
3. Reference allele
4. Alternate allele
5. Quality score (not present in 1k bulls data)
6. Type (can be either a "SNP" or an "INDEL")
7. Mutation (different classifications based on the mutation's presence in or outside of a gene)
8. Priority (can be "HIGH," "MODERATE," "LOW," or "MODIFIER" in decreasing order of interest)
9. Gene (gene name or ensembl accession)
10. AA or Amino Acid (only present in cases of nonsynonymous mutations)
11. Columns 11 through the last column are animal genotypes