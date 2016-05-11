# Identifying variants within Joel's bulls
---
*5/11/2016*

These are my notes on the alignment and variant calling performed on Joel's bulls as part of our BARD funded grant.

## Table of Contents
* [Organizing the data](#organizing)

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