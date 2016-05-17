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

I am missing 2 bulls from Joel's list. Here are my notes on how to reconcile the lists.

> 3850: /seq1/bickhart/side_projects/joels_bulls

```bash
perl -lane 'print $F[1];' < /work1/grw/grw1/Joel/With_Seq_41.txt > georges_list_41.txt
cat final_file_joel_bulls.list 1000_bulls_sequenced_joels_bulls.list > my_list_bulls.txt

# I used vim to format the my_list_bulls.txt file so that the names were uniform
# Still quite a few were missing
```

I need to go back to the source and use the last list of bulls to reconcile the differences. 

```bash
for i in ../../../1000_bulls_bams/HO*.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; done > 100_bulls_holsteins.list

# ID bulls in 1000 bulls data
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt
mv group_1_2.txt 1000_bulls_joel_with_canada.list
mv group_2.txt bulls_minus_1000_data.list

# I realized that I can generate the association for all files at once
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: canadian_bulls.list
File Number 4: 100_bulls_holsteins.list
Set     Count
1       364
1;2     19
1;2;3   9
1;2;3;4 3
1;2;4   5
1;4     26
2       34
2;4     1
4       2

mv group_1_2.txt 1000_bulls_only_joel.list
cat group_1_2_4.txt group_1_2_3.txt group_1_2_3_4.txt group_2_4.txt > joel_already_controlled.list

# George has a list of bulls that have sequence, supposedly. Let's use that list instead
perl -lane 'print $F[1];' < With_Seq_41.txt > georges_list_sequenced.list
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt georges_list_sequenced.list canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: georges_list_sequenced.list
File Number 4: canadian_bulls.list
File Number 5: 100_bulls_holsteins.list
Set     Count
1       364
1;2;3   19
1;2;3;4 9
1;2;3;4;5       3
1;2;3;5 5
1;5     26
2       30
2;3     4
2;3;5   1
5       2

# OK, let's find out who is in the 2;3 group and where they're located
head group_2_3.txt
HODEU000000254210
HOUSA000002247437
HONLD000839380546
HODEU000000830287

for i in HODEU000000254210 HOUSA000002247437 HONLD000839380546 HODEU000000830287; do getids $i | perl -e '<>; <>; while(<>){chomp; @s = split(/\s+/); print "$s[1]\n";}' > $i.names; done
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl 1000_bulls_sequenced_reformatted.list HODEU000000254210.names
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: HODEU000000254210.names
Set     Count
1       425
1;2     1
2       5
# I did this for the rest, and only found discrepencies with HOUSA000002247437 and HONLD000839380546

# Testing George Liu's list for preferred ids
for i in `cat george_lius_list.txt`; do getids $i | perl -e '<>; <>; while(<>){chomp; @s = split(/\s+/); print "$s[1]\n";}'; done > george_lius_list_altnames.txt
grep HOUSA000002247437 george_lius_list_altnames.txt
grep HONLD000839380546 george_lius_list_altnames.txt

# Nothing
# HONLD000839380546 was in the 1000 bulls list, but I was unable to find it because of my permissive substitution script
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt georges_list_sequenced.list canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: georges_list_sequenced.list
File Number 4: canadian_bulls.list
File Number 5: 100_bulls_holsteins.list
Set     Count
1       363
1;2;3   20
1;2;3;4 9
1;2;3;4;5       3
1;2;3;5 5
1;5     26
2       30
2;3     3
2;3;5   1
5       2

# So we remove HOUSA000002247437 but process the rest
grep -v 'HOUSA000002247437' group_2_3.txt > alt_1000_bulls_ids.txt
mv group_1_2_3.txt 1000_bulls_presumptive_list.list
cat group_2_3_5.txt group_1_2_3_5.txt group_1_2_3_4_5.txt group_1_2_3_4.txt > joels_bulls_we_already_have.list

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list HODEU000000254210.names
# This is HODEU000578194407
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list HODEU000000830287.names
# This is HODEU002261530135
```

The files that contain the proper sequenced bulls are:
* /seq1/bickhart/side_projects/joels_bulls/1000_bulls_presumptive_list.list
* /seq1/bickhart/side_projects/joels_bulls/joels_bulls_we_already_have.list

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

I need to process the following, remaining chromosomes for Joel: 15, 3, 20, 6, 14, 16, 1, 18, 19, 7, 13, 28, 10.

I've created a short bash script to thread this to speed things up.

#### /seq1/1kbulls_annotatedvcf/process_script.sh
```bash
# $1 = chromosome

beaglevcf=${1}-Beagle-Run5.eff.vcf.gz
beagleuncomp=${1}-Beagle-Run5.eff.vcf
progout=${1}_joels_holstein_subsection.tab

# Uncompress and recompress files
gunzip $beaglevcf
bgzip $beagleuncomp
bcftools index $beaglevcf

perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f $beaglevcf -o $progout -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list
```

And here is the code that I'm using to farm out the jobs.

> 3850: /seq1/1kbulls_annotatedvcf/

```bash
for i in Chr15 Chr3 Chr20 Chr6 Chr14 Chr16 Chr1 Chr18 Chr7 Chr13 Chr28 Chr10; do echo $i; sh process_script.sh $i & done

# Chr19 was already bgzipped, so I will process that directly
perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f Chr19-Beagle-Run5.eff.vcf.gz -o Chr19_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list
```

#### I just redid the 1000 bulls list, so I need to reprocess the data so that Joel has the proper animals

```bash
for i in Chr10 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr28 Chr29 Chr6 Chr20 Chr1 Chr3; do echo $i; perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f ${i}-Beagle-Run5.eff.vcf.gz -o ${i}_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_presumptive_list_reformatted.list; done
```