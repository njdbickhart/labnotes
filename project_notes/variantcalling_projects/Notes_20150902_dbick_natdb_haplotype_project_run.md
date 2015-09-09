# Natdb haplotype sequencing project
---
*9/2/2015*

These are my notes and cmdline commands for running the variant analysis pipeline on sequence data for the haplotype sequencing project.

## Table of Contents
* [Testing the setup on NextSeq500 data](#testone)
* [Restarting the pipeline](#restart)

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
