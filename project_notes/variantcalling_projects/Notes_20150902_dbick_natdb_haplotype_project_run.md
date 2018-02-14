# Natdb haplotype sequencing project
---
*9/2/2015*

These are my notes and cmdline commands for running the variant analysis pipeline on sequence data for the haplotype sequencing project.

## Table of Contents
* [Testing the setup on NextSeq500 data](#testone)
* [Restarting the pipeline](#restart)
* [Pilot phase 1 UMD3 align](#umd1)
* [Alignment to modified ARSUCD-v14](#arsucd)

<a name="testone"></a>
## Testing the setup on NextSeq500 data

```bash
# Save fastq list
ls /mnt/cifs/bickhart-qnap/NextSeq/150818_NS500432_0013_AH1513BGXX/H1513BGXX/DBickhart-NatDB-2015-08-18/*/*.fastq.gz > H1513BGXX.fastq.list

```

Now, I used VIM to process the spreadsheet. It needs the following ordered columns for each line (delimited by tabs):

1. Fastq for read1
2. Fastq for read2
3. Library name
4. Sample name

I used the Animal's common name as the library name, and the Holstein ID as the sample name. I'm going to run the program in a screen so that I can safely leave the process running the entire time.

> Blade14: /mnt/iscsi/vnx_gliu_7/haplotype_project

```bash
# Attaching a "screen"
screen -a

# Now to invoke the pipeline
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs H1513BGXX.spreadsheet.tab --output H1513BGXX --reference /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa --coords /mnt/iscsi/vnx_gliu_7/reference/samtools_chr_segs.txt --threads 10
```

Let's dissect the input:
* **--fastqs** H1513BGXX.spreadsheet.tab <- This the spreadsheet I just manually created with all of the fastq locations
* **--output** H1513BGXX <- This is a folder that will be created that will contain all of the files
* **--reference** /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa <- This is a previously created reference genome fasta file that I'm using for this run. It's indexed with BWA and Samtools
* **--coords** /mnt/iscsi/vnx_gliu_7/reference/samtools_chr_segs.txt <- This is a file that splits the chromosomes into 50 megabase chunks for faster Samtools processing
* **--threads** 10 <- I want the program to use 10 threads

-- 
*9/8/2015*

I had a chance to check the output, and it looks like the pipeline terminated prematurely, let's see why:

> Blade14: /mnt/iscsi/vnx_gliu_7/haplotype_project

```bash
# Let's see what's in the directory
ls H1513BGXX/*
	H1513BGXX/H1513BGXX.spreadsheet.tab.fastqc  H1513BGXX/MergedBamPipeline.1441245029.log
	
	H1513BGXX/fastqc:
	HOCAN000006392464  HOCAN000007383072  HOUSA000055181279  HOUSA000130498623
	
	H1513BGXX/HOCAN000006392464:
	HOCAN000006392464.1.nodup.bam.metrics  HOCAN000006392464.3.nodup.bam.metrics  HOCAN000006392464.merged.bam
	HOCAN000006392464.2.nodup.bam.metrics  HOCAN000006392464.4.nodup.bam.metrics  HOCAN000006392464.merged.bam.bai
	
	H1513BGXX/HOCAN000007383072:
	HOCAN000007383072.1.nodup.bam.metrics  HOCAN000007383072.3.nodup.bam.metrics  HOCAN000007383072.merged.bam
	HOCAN000007383072.2.nodup.bam.metrics  HOCAN000007383072.4.nodup.bam.metrics  HOCAN000007383072.merged.bam.bai
	
	H1513BGXX/HOUSA000055181279:
	HOUSA000055181279.1.nodup.bam          HOUSA000055181279.2.nodup.bam.bai      HOUSA000055181279.4.nodup.bam
	HOUSA000055181279.1.nodup.bam.bai      HOUSA000055181279.2.nodup.bam.metrics  HOUSA000055181279.4.nodup.bam.bai
	HOUSA000055181279.1.nodup.bam.metrics  HOUSA000055181279.3.bam                HOUSA000055181279.4.nodup.bam.metrics
	HOUSA000055181279.2.nodup.bam          HOUSA000055181279.3.bam.bai
	
	H1513BGXX/HOUSA000130498623:
	HOUSA000130498623.1.nodup.bam.metrics  HOUSA000130498623.3.nodup.bam.metrics  HOUSA000130498623.merged.bam
	HOUSA000130498623.2.nodup.bam.metrics  HOUSA000130498623.4.nodup.bam.metrics  HOUSA000130498623.merged.bam.bai
	
	H1513BGXX/vcfs:

# No SNP vcfs, and no HOUSA000055181279.merged.bam. Apart from that, everything else completed fine.

# Let's look at the log file:
tail H1513BGXX/MergedBamPipeline.1441245029.log
	...
	05/08/115 - 00:50:57 | [Main] - Submit: runSamtoolsBCFCaller, args => /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa, H1513BGXX/vcfs/combined.samtools.chr29:1-51505224.bcf, H1513BGXX/HOCAN000007383072/HOCAN000007383072.merged.bam,H1513BGXX/HOUSA000130498623/HOUSA000130498623.merged.bam,H1513BGXX/HOUSA000055181279/HOUSA000055181279.merged.bam,H1513BGXX/HOCAN000006392464/HOCAN000006392464.merged.bam, chr29:1-51505224, 1
	05/08/115 - 00:50:57 | [Main] - Submit: runSamtoolsBCFCaller, args => /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa, H1513BGXX/vcfs/combined.samtools.chrX:1-50000000.bcf, H1513BGXX/HOCAN000007383072/HOCAN000007383072.merged.bam,H1513BGXX/HOUSA000130498623/HOUSA000130498623.merged.bam,H1513BGXX/HOUSA000055181279/HOUSA000055181279.merged.bam,H1513BGXX/HOCAN000006392464/HOCAN000006392464.merged.bam, chrX:1-50000000, 1
	05/08/115 - 00:50:57 | [Main] - Submit: runSamtoolsBCFCaller, args => /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa, H1513BGXX/vcfs/combined.samtools.chrX:50000001-100000000.bcf, H1513BGXX/HOCAN000007383072/HOCAN000007383072.merged.bam,H1513BGXX/HOUSA000130498623/HOUSA000130498623.merged.bam,H1513BGXX/HOUSA000055181279/HOUSA000055181279.merged.bam,H1513BGXX/HOCAN000006392464/HOCAN000006392464.merged.bam, chrX:50000001-100000000, 1
	05/08/115 - 00:50:57 | [Main] - Submit: runSamtoolsBCFCaller, args => /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa, H1513BGXX/vcfs/combined.samtools.chrX:100000001-148823899.bcf, H1513BGXX/HOCAN000007383072/HOCAN000007383072.merged.bam,H1513BGXX/HOUSA000130498623/HOUSA000130498623.merged.bam,H1513BGXX/HOUSA000055181279/HOUSA000055181279.merged.bam,H1513BGXX/HOCAN000006392464/HOCAN000006392464.merged.bam, chrX:100000001-148823899, 1
	05/08/115 - 00:50:57 | End - Run completed successfully

# I don't have error messages saved to the log file, unfortunately. The time on the log entries is far too short though. They all failed one after another.

# the Java Virtual Machine (JVM) runtime created an error log
less hs_err_pid26982.log
	# There is insufficient memory for the Java Runtime Environment to continue.
	# Native memory allocation (mmap) failed to map 5,632,425,984 bytes for committing reserved memory.
	# Possible reasons:
	#   The system is out of physical RAM or swap space
	#   In 32 bit mode, the process size limit was hit
	# Possible solutions:
	#   Reduce memory load on the system
	#   Increase physical memory or swap space
	#   Check if swap backing store is full
	#   Use 64 bit Java on a 64 bit OS
	#   Decrease Java heap size (-Xmx/-Xms)
	#   Decrease number of Java threads
	#   Decrease Java thread stack sizes (-Xss)
	#   Set larger code cache with -XX:ReservedCodeCacheSize=
	# This output file may be truncated or incomplete.
	#
	#  Out of Memory Error (os_linux.cpp:2718), pid=26982, tid=140373394138880

```

Based on the files I have above, and the JVM error, I suspect that the server ran out of allocatable memory during the "mark duplicates" step, as that is written in Java. Apparently, the memory on the server was completely used up by competing processes. Once that job failed, the rest of the pipeline stopped. 

Because I only delete files if the end result file exists, I can pick up where the pipeline left off and complete the dataset. I want to capture all error messages in a separate log as well, so I need to implement that as soon as possible.

<a name="restart"></a>
## Restarting the pipeline

In case something goes wrong, there are several ways to correct stages in the pipeline that do not require rerunning the whole thing from scratch! Let's work with the above error and walk it through to the end.

> Blade14: /mnt/iscsi/vnx_gliu_7/haplotype_project/H1513BGXX

```bash
# We need to mark optical/pcr duplicates first
java -jar ~/picard-tools-1.85/MarkDuplicates.jar I=HOUSA000055181279/HOUSA000055181279.3.bam O=HOUSA000055181279/HOUSA000055181279.3.nodup.bam M=HOUSA000055181279/HOUSA000055181279.3.nodup.bam.metrics VALIDATION_STRINGENCY=SILENT
	# the VALIDATION_STRINGENCY option is a requirement for BWA aligned bams. 
	# The author of picard tools decided to make BAMs fail validation if they had non-zero mapping scores for soft-clip alignments
	# The author of BWA decided to give these alignments scores in the first place.
	# Hence, the conflict

# now to index the bam
samtools index HOUSA000055181279/HOUSA000055181279.3.nodup.bam

# now to merge all of the HOUSA000055181279 bams
# But first, I need to create a header file because the normal bam merge process destroys the readgroups
for i in HOUSA000055181279/*.bam; do samtools view -H $i; done | sort | uniq > HOUSA000055181279/HOUSA000055181279.header.sam
# Now I'm going to manually edit the file to change the sort order and remove redundant header info

# Ok, merger step
samtools merge -h HOUSA000055181279/HOUSA000055181279.header.sam -@ 5 HOUSA000055181279/HOUSA000055181279.merge.bam HOUSA000055181279/HOUSA000055181279.1.nodup.bam HOUSA000055181279/HOUSA000055181279.2.nodup.bam HOUSA000055181279/HOUSA000055181279.3.nodup.bam HOUSA000055181279/HOUSA000055181279.4.nodup.bam

# Indexing the merged bam
samtools index HOUSA000055181279/HOUSA000055181279.merge.bam
```

<a name="umd1"></a>
## Pilot phase 1 UMD3 align

I am going to set up an initial alignment of the reads against the UMD3 reference genome while I wait for the cattle assembly to go public on NCBI. First, I need to tag all bam sequence sent by WUSTL in my spreadsheet tab file.

The relational data for library name and bam ID are found in this file in my sequencing projects folder:
* final_UNALIGNED_samplemap_june_2017.xlsx

> Assembler 2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing

```bash
# Generating a file list
ls /mnt/nfs/SequenceData/wustl/set*/*/*/*/*.bam > wustl_bam_file_locs.list

# Copying library and animal name info to the server
vim wustl_relational_info.tab
dos2unix wustl_relational_info.tab

# Now, removing the "BTAN-" identifier and formatting the spreadsheet information
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; $_ =~ s/BTAN-//g; @s = split(/\t/); $data{$s[2]} = [$s[0], $s[1]];} close IN; open(IN, "< $ARGV[1]"); while($l = <IN>){chomp $l; @s = split(/\//, $l); $d = $data{$s[-1]}; print "$l\t$d->[0]\t$d->[1]\n";} close IN;' wustl_relational_info.tab wustl_bam_file_locs.list > wustl_formatted_bam_data.tab
```

Scratch the UMD3 align -- Bob will take care of that faster than I can. I will instead work on alignment to a modified new reference containing the IGC haplotypes.

<a name="arsucd"></a>
## Alignment to modified ARSUCD-v14

I've copied Ben's ARS-UCDv14 to the following directory:

> /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114/ARS-UCD1.0.14.clean.fasta

The fasta has 73 missing supplementary scaffolds (contaminants or vectors) and 5 chromosomes/scaffolds with changes compared to v13.

```R
v13 <- read.delim("../ars_ucd_113/ARS-UCD1.0.13.fasta.fai", header=FALSE)
v14 <- read.delim("ARS-UCD1.0.14.clean.fasta.fai", header=FALSE)
library(dplyr)

joined_summary <- left_join(v13[,c(1,2)], v14[,c(1,2)], by="V1")
filter(joined_summary, is.na(V2.y)) %>% count
# 73 missing
filter(joined_summary, V2.x != V2.y)
#                     V1      V2.x      V2.y
#1                    10 103152713 103152712
#2 Leftover_ScbfJmS_1601     18130     18110
#3 Leftover_ScbfJmS_1857     26100     26078
#4 Leftover_ScbfJmS_1867    214469    214448
#5 Leftover_ScbfJmS_1961     37852     37820
```

I am going to add several scaffolds to this assembly and make a separate (v14+igc) assembly version for alignment. Based on Doro's email, I should keep/modify the following fastas:

> Contigs for SNP discovery
I had a look through your list of BAC clones and compared them with the final BAC overview list I created last year (attached to this email). The good news is that a lot of the BACs do not need to be included in the alignment/mapping exercise (see table below).

>1)	NKC – all clones that contain NKC can be taken out. The reference assembly is good for this region and we have not found haplotypes that differed in gene content. 

>2)	MHC – most of the BAC clones that contain MHC aligned to regions that are not part of the classical class I polymorphic region. The one clone to keep is the LIB14427_MHC which contains a different class I haplotype to the reference genome. 

>3)	LRC – most LRC clones are actually not LRC or contain polymorphic LILR or KIR genes. The one to keep in the analysis is the LIB14413_LRC, as it contains LILR genes.

>We do have a few more additional files that can be added to the haplotype list, which are attached as fasta files to this email. The summary below includes all files to be included.

>MHC
>1)	TPI_4222_MHC classI haplotype: This was sequenced and assembled from BAC clones at Pirbright. This haplotype contains 4 classical MHC genes (gene 1, gene 2, gene 4, gene 6).
>2)	Dominette_MHC class I haplotype: This is essentially the classical class I region from the reference genome UMD3.1. This haplotype contains 2 classical MHC genes (gene 5 and gene 2).
>3)	LIB14427_MHC class I haplotype (from your list): This is a partial haplotype containing 2 classical MHC genes (gene 3 and gene2).

>LRC
>1)	TPI_4222_KIR haplotype 1: This was sequenced and assembled from BAC clones at Pirbright. This haplotype contains the KIR region.
>2)	LIB14413_LRC (from your list), contains LILR genes
>3)	LIB14604_CH240_391K10: Another partial KIR haplotype.
>4)	LIB14602_CH240-370M3: contains LILR genes

>We would say that before you decide to include any of the above please take a quick look at how similar they are to the ARS reference you are using so you get an idea of what to expect.

So, I need to do nucmer comparisons with the novel haplotype fragments, and then keep only two of the previously assembled BAC haplotypes: LIB14413_LRC and LIB14427_MHC. First the nucmer aligns.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
java -jar ../../binaries/CombineFasta/store/CombineFasta.jar standardize -f /mnt/nfs/nfs2/dbickhart/igc/TPI4222_A14_MHCclassI.fas -o TPI4222_A14_MHCclassI.fasta
java -jar ../../binaries/CombineFasta/store/CombineFasta.jar standardize -f /mnt/nfs/nfs2/dbickhart/igc/CH240_391K10_polished.fasta -o CH240_391K10_polished.fasta
java -jar ../../binaries/CombineFasta/store/CombineFasta.jar standardize -f /mnt/nfs/nfs2/dbickhart/igc/Domino_MHCclassI_gene2-5hapl.fasta -o Domino_MHCclassI_gene2-5hapl.fasta
java -jar ../../binaries/CombineFasta/store/CombineFasta.jar standardize -f /mnt/nfs/nfs2/dbickhart/igc/LRC_CH240_370M3.fas -o LRC_CH240_370M3.fasta
java -jar ../../binaries/CombineFasta/store/CombineFasta.jar standardize -f /mnt/nfs/nfs2/dbickhart/igc/TPI_4222_LRC_hap1.fas -o TPI_4222_LRC_hap1.fasta


sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta /mnt/nfs/nfs2/dbickhart/igc/CH240_391K10_polished.fasta
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta /mnt/nfs/nfs2/dbickhart/igc/Domino_MHCclassI_gene2-5hapl.fasta
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta /mnt/nfs/nfs2/dbickhart/igc/LRC_CH240_370M3.fas
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta /mnt/nfs/nfs2/dbickhart/igc/TPI4222_A14_MHCclassI.fas
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta /mnt/nfs/nfs2/dbickhart/igc/TPI_4222_LRC_hap1.fas

# Damn, the reference file size was too large. Going to use the whole chromosome as the reference instead
# chr23 for MHC and chr18 for LRC. NKC is apparently well covered in this assembly!
samtools faidx ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta 23 > ver14_chr23.fa
samtools faidx ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta 18 > ver14_chr18.fa

sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ver14_chr23.fa TPI4222_A14_MHCclassI.fasta
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ver14_chr23.fa Domino_MHCclassI_gene2-5hapl.fasta

sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ver14_chr18.fa TPI_4222_LRC_hap1.fasta
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ver14_chr18.fa CH240_391K10_polished.fasta
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ver14_chr18.fa LRC_CH240_370M3.fasta

# Rewriting fasta headers according to the filenames and associations
head *.fasta
==> CH240_391K10_polished.fasta <==
>CH240_391K10_Vector_Trim_polished
==> Domino_MHCclassI_gene2-5hapl.fasta <==
>CH103G5_and_CH463M1_and_CH252J17
==> LRC_CH240_370M3.fasta <==
>CH240_370M3
==> TPI4222_A14_MHCclassI.fasta <==
>TPI4222_A14_MHCclassI
==> TPI_4222_LRC_hap1.fasta <==
>HF_LRC_hap1


cat *.fasta > pirbright_combined.fa
vim pirbright_combined.fa

# Now selecting the USDA haplotypes that Doro suggested were not redundant
samtools faidx /mnt/nfs/nfs2/dbickhart/igc/usda_cumulative_igc_haplotypes.fa LIB14427_MHC LIB14413_LRC > usda_nonredundant_haplotypes.fa

# And finally, sending them all to the final fasta
cat ../ars_ucd_114/ARS-UCD1.0.14.clean.fasta pirbright_combined.fa usda_nonredundant_haplotypes.fa > ARS-UCD1.0.14.clean.wIGCHaps.fasta
# Parity check to ensure that it is all present and clear
samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta

sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 --wrap="module load bwa; bwa index ARS-UCD1.0.14.clean.wIGCHaps.fasta"

# OK, now we're ready to run. Just want to make sure that I'm generating a proper script, so I will test one bam file for alignment now
head -n1 ../../../bickhart-users/natdb_sequencing/wustl_formatted_bam_data.tab > test_spreadsheet_subset.tab
```

I need to switch folders to generate the scripts since my pipeline was not designed to run outside of the working directory!

> Assembler2: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project

```bash
perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsSlurm.pl -b test -t /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/test_spreadsheet_subset.tab -f /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta
sbatch test/001HO05072/scripts/bwaAlign_555751224267203.sh

# I learned that I added error codes about the fasta reference being in the right directory. Since the NFS2 mount
# is pretty hammered, it might make sense to copy the reference files to this directory too.
cp /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta* ./

# Steve asked me to submit all of the jobs so that they're ready to go when he restarts the cluster
perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsSlurm.pl -b aligns -t /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/wustl_formatted_bam_data.tab -f ARS-UCD1.0.14.clean.wIGCHaps.fasta -m
#	Generated 2629 alignment scripts for 173 samples!
```

## Round 2 alignment

We just received the remainder of the bam files and I need to queue them up for alignment. I am going to align them in a separate folder and then queue up a separate merger later. 

> Assembler2: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project 

```bash
ls /mnt/nfs/SequenceData/wustl/2017-12/xfer.genome.wustl.edu/gxfer1/90727643070659/*.bam > wustl_round2_bam_file.list

# the samplemap csv file was carriage return delimited
perl -e 'while(<>){$_ =~ s/\r/\n/g; print $_;}' < /mnt/nfs/SequenceData/wustl/2017-12/xfer.genome.wustl.edu/gxfer1/90727643070659/Holstein.samplemap.csv > round2_sample_map.csv

# Formatting the alignment input tab file
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %h; while(<IN>){chomp; @s = split(/,/); @bsegs = split(/-/, $s[6]); $s[6] = $bsegs[1]; $h{$s[1]}->{$s[2]}->{$s[3]} = $s[6] . "";} close IN; open(IN, "< $ARGV[1]"); while($y = <IN>){chomp $y; @s = split(/\//, $y); @j = split(/[_\.]/, $s[-1]); $t = $h{$j[1]}->{$j[2]}->{$j[3]}; print "$y\t$t\-lib2\t$t\n";} close IN;' round2_sample_map.csv wustl_round2_bam_file.list > wustl_round2_bam_remaining_files.tab

# Running the pipeline
perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b round2 -t wustl_round2_bam_remaining_files.tab -f ARS-UCD1.0.14.clean.wIGCHaps.fasta -m
```

## Preparing for variant calling

I still need to mark dups, realign indels and then merge the 35 improved samples to the list. I will do this with a custom script.

#### indel_realign_markdups.sh

```bash
#!/usr/bin/sh
#SBATCH --nodes=1
#SBATCH --mem=25000
#SBATCH --ntasks-per-node=5

# This is a oneshot script designed to realign indels, mark dups and otherwise prepare BAMs for variant calling
# $1 is the sorted merged bam to realign/modify
# $2 is the reference fasta
realign_targs=${1}.targets
realigned_bam=${1}.realn.bam
duplicates_metrics=${1}.dup.metrics
dedup_bam=${1}.dedup.bam

START=$(date +%s.%N)
module load bwa samtools gatk/3.7 picard/2.9.2

echo $GenomeAnalysisTK -T RealignerTargetCreator -R $2 -I $1 -o $realign_targs

$GenomeAnalysisTK -T RealignerTargetCreator -R $2 -I $1 -o $realign_targs

echo $GenomeAnalysisTK -T IndelRealigner -R $2 -I $1 -targetIntervals $realign_targs -o $realigned_bam

$GenomeAnalysisTK -T IndelRealigner -R $2 -I $1 -targetIntervals $realign_targs -o $realigned_bam

samtools index $realigned_bam

echo java -jar $PICARD MarkDuplicates I=$realigned_bam O=$dedup_bam M=$duplicates_metrics

java -jar $PICARD MarkDuplicates I=$realigned_bam O=$dedup_bam M=$duplicates_metrics
samtools index $dedup_bam
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "Executed in $DIFF time from $START and $END"
```

And now to queue up all the tasks to get things going.

> Assembler2: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project

```bash
ls ./*/*/*.sorted.merged.bam > sorted.merged.bam.list

cat sorted.merged.bam.list | xargs -I {} sbatch indel_realign_markdups.sh {} ARS-UCD1.0.14.clean.wIGCHaps.fasta

# The GATK strikes again! This time it complained about the fact that I named the indel regions as a ".targets" file instead of a ".intervals" file.
```

## Variant calling

I need to select suitable BAM files for the analysis. It looks like there are some bams with very low coverage that may impact the variant calling. I would want to remove them to avoid biasing our variant calls.

> Assembler2: 

```bash
ls ./*/*/*.dedup.bam > dedup.final.bam.list
module load samtools; perl ~/sperl/sequence_data_scripts/getBamStats.pl -n dedup.final.bam.list -o summary_stats_dedup_bams.2018.02.13.tab

# Let's see how many bams are under 5X coverage
perl -lane 'if($F[4] < 6){print $_;}' < summary_stats_dedup_bams.2018.02.13.tab
BamName TotalReads      MappedReads     UnmappedReads   RawXCov MapXCov AvgRawChrcov    AvgMapChrcov
./aligns/001HO11393/001HO11393.sorted.merged.bam.dedup.bam      164363326       162629074       1734252 2.70031330886414        2.67182140698742 19.4390582126377 19.3689079864248
./aligns/007HO11737/007HO11737.sorted.merged.bam.dedup.bam      172415432       169853345       2562087 3.58796098994332        3.53464402114182 29.5061149049186 29.3991967532581
./aligns/200HO00405/200HO00405.sorted.merged.bam.dedup.bam      172835797       170806510       2029287 2.46090600427256        2.43201219494965 22.6916124034632 22.6223531838276
./round2/023HO06635/023HO06635.sorted.merged.bam.dedup.bam      64439402        62251035        2188367 3.55242707180079        3.43178637786891 23.115716386796  23.052625147377

# Not too many. Now let's check high deviation from the mapping percentages
perl -lane 'if($F[5] == 0){next;} if($F[2]/$F[1] < 0.97){print $_;}' < summary_stats_dedup_bams.2018.02.13.tab
./aligns/007HO10676-0676_qia_test/007HO10676-0676_qia_test.sorted.merged.bam.dedup.bam  449513407       433637102       15876305        24.7808568453848  23.905625016229 132.590328508592        132.130627379593
./round2/014HO07519/014HO07519.sorted.merged.bam.dedup.bam      492766949       472660610       20106339        27.165345974444 26.0569200617019 200.678179098475 199.926318938763
./round2/023HO06635/023HO06635.sorted.merged.bam.dedup.bam      64439402        62251035        2188367 3.55242707180079        3.43178637786891 23.115716386796  23.052625147377	

# The round2 version of the 014HO07519 bam has a 95% mapping rate! That's pretty suspicious, considering the round1 bam has a 98% mapping rate

# So I will drop all bams under 5X coverage and the < 97% mapping rate bams (except the qia_test bam, as I do not have a replicate of that sample!
./aligns/001HO11393/001HO11393.sorted.merged.bam.dedup.bam
./aligns/007HO11737/007HO11737.sorted.merged.bam.dedup.bam
./aligns/200HO00405/200HO00405.sorted.merged.bam.dedup.bam
./round2/023HO06635/023HO06635.sorted.merged.bam.dedup.bam
./round2/014HO07519/014HO07519.sorted.merged.bam.dedup.bam

# The filtered list of bams is dedup.filtered.final.bam.list
# Now to queue up the samtools mpileup scripts
perl -ne '$_ =~ s/^\./\/mnt\/nfs\/nfs1\/derek.bickhart\/CDDR-Project/; print $_;' < dedup.filtered.final.bam.list > dedup.filtered.final.bam.fullpath.list
perl ~/sperl/sequence_data_pipeline/samtoolsMpileupBamsSlurm.pl -b vcfs -t dedup.filtered.final.bam.fullpath.list -f ARS-UCD1.0.14.clean.wIGCHaps.fasta