# Natdb haplotype sequencing project
---
*9/2/2015*

These are my notes and cmdline commands for running the variant analysis pipeline on sequence data for the haplotype sequencing project.

## Table of Contents
* [Testing the setup on NextSeq500 data](#testone)

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

> Blade14: 

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

