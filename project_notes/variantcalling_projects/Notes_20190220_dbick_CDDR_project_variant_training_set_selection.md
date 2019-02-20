# CDDR project variant selection list
---
*2/20/2019*

These are my notes on the generation of a list of training set variant sites for the CDDR analysis.

## Table of Contents

## Summary of what's been done

Let's first go through what's been done according to Kiranmayee's notes and tables.

#### Preprocessing of data
* Alignment of reads to ARS-UCDv1.4 and variant calling using samtools mpileup
* [Filtration of SNPs based on expected statistics and training set size](https://github.com/bkiranmayee/CDDR_Variants_filtering/blob/master/notes.md)
* [Conversion of ARS-UCDv1.4 SNP sites to ARS-UCDv1.2 sites](https://github.com/bkiranmayee/CDDR-Assembly-liftover-Project/blob/master/labnotes.md)
* Realignment of sequence reads to ARS-UCDv1.2 reference directly
* [Download and transfer of additional files for cross-validation](https://github.com/bkiranmayee/My_Labnotes/blob/master/Additional_Holsteins_WGS_alignment.md)

#### Creating test and training sets for models
* Filtering variant sites based on hard filters
* [Using hard filtered variant sites to id holstein specific variants on the same reference dataset](https://github.com/bkiranmayee/My_Labnotes/blob/master/vqsr_results_analysis.md)
* [Progressively filtering SNPs based on assumptions of homozygosity in haplotype regions](https://github.com/bkiranmayee/My_Labnotes/blob/master/progressiveSelectionKBv3.pl)


We've tested some assumptions of the data, and it looks like there are very big problems with determining what is a reliable, true variant site from the noise. We need to do a better job of identifying (a) the problems that cause bias and (b) variant sites that provide improved prediction accuracy. Using statistical estimates will be key here as well as getting a good set of training data.

## What we need moving forward

Now, let's summarize what's been going on in the literature and what we can do to improve the state of the field. 

Recently, work by [Hubert Pausch's group](https://www.biorxiv.org/content/10.1101/460345v2) has tested the utility of graph-based genome alignment for SNP/INDEL discovery. The basis for the improvement is in the removal of reference-allele-alignment bias, so if such regions in the cattle genome could be discovered apart from this study, they could be of use in our modeling scheme. Work by [Ben Langmeade's group](file:///C:/Users/dbickhart/Zotero/storage/5G2MYUU2/s13059-018-1595-x.html) provides a means for us to interrogate this. 


#### Datasets

##### SNP calls for model training
* SNP calls for first filter (172 bulls only) (consensus of three callers):
	* Transfer the ARS-UCDv14 aligned bams to the rumenmicrobial_project folder (**Kiranmayee**)
		* Realign to ARS-UCDv1.2
		* Mark duplicates and immediately delete the previous realigned bam and the ARS-UCDv14 aligned bam
		* Realign INDELs using GATK INDEL target creator. 
	* Samtools mpileup calls on ARS-UCDv1.2 (**Kiranmayee**). Use the script but in < 1000 segments. 
	* GATK BQSR(?; Test on one BAM and assess the results)-modified raw variant calls (**Kiranmayee**)
	* FreeBayes SNP calls on BAMS  (**Derek**)
* SNP calls for second filter (Haplotype prior information):
	* Filtered dataset based on haplotype assignment (using [perl script](https://github.com/bkiranmayee/My_Labnotes/blob/master/progressiveSelectionKBv3.pl)) (**Kiranmayee**)
	* Assess Variant quality and VCF metrics before deciding if a [third filter is needed](https://software.broadinstitute.org/gatk/documentation/article?id=11069). (**Both**)
* SNP calls for third filter (hard filter -- if needed!):
	* Use GATK raw calls as basis for consensus SNPs
	* Use recommended quality metrics in the GATK article to remove problematic SNPs
	* Remove BovineHD and 1kBulls snps first. (**Kiranmayee**)
* Products:
	* A list of SNPs for model training (truth=true, training=true) in GATK VQSR.
	* A list of SNPs from all three callers
	* A list of metrics and plots from the variants from all three callers. 

##### Model evaluation
* BovineHD variants for 172 animals (genotype concordance estimates) (**Derek**)
* Parent-progeny information for 172 animals (opposing homozygotes) (**Derek**)
* Variant datasets for GATK training:
	* BovineHD; known=false, training=true, truth=true; prior=15.0
	* 1kbull_HQF; known=false, training=true, truth=false, prior=10.0
	* **NovelHolstein**; known=false, training=true, truth=false, prior=8.0
	* dbsnp; known=true, training=false, truth=false, prior=2.0
* Cross-validation
	* Check model performance in original dataset vs original calls without breed-specific data. (already done)
	* Check model performance (with breed-specific data vs non-breed-specific data) in SRA datasets aligned to ARS-UCDv1.2 ("leave-one-out" bootstrapping with VQSR model training (?)) (**Kiranmayee**)
	* Check variant calls against those made by graphtyper using the consensus SNPs fed into the index from the original consensus calls (see first bullet point above in the "SNP calls for model training" section). (**Derek**)


## SNP calling for consensus

Let's start by generating SNP calls for the three callers.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/kiranmayee/CDDR

```bash
# Generate the windows for the samtools mpileup segment calling
perl -lane 'for($x = 1; $x < $F[1] + 10000000; $x += 10000000){$e = $x + 10000000; if($x > $F[1]){$e = $F[1];} print "$F[0]\:$x\-$e";}' < ARSUCD1.2.current_ref.fa.fai > ARSUCD1.2.current_ref.fa.samtools.mpileup.wins
```