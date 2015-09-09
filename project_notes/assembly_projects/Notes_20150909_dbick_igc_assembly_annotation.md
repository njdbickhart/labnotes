# IGC assembly and annotation work
---
*9/9/2015*

These are my notes on the assembly/association of BAC clones/contigs to their most likely physical map locations.

## Table of Contents
* [First batch unitig assignments](#firstbatch)
* [Initial BAC end mapping information table](#bacend)

<a name="firstbatch"></a>
## First batch unitig assignments

I just received an email from Tim with six assembled unitigs from the first round of sequencing. Here's Tim's email explaining the situation:

> The first two BAC pools were sequenced over this past holiday (in U.S.) weekend.  Went pretty well, there were two pools of four BACs, and six appear to have assembled, an excellent outcome in my view.  However, since we don’t know what sequence to expect exactly for each BAC, I cannot directly assign them without further analysis which I am looking for help on since I have a torrent of other stuff to get done before end of September.

...

> Note that the pools contain one BAC each from the four genomic regions Derek identified in his spreadsheet; chr18, LRC, MHC, and NKC.  In  LIB14363, the clones included were 389p14, 170i19, 503o23, and 49i22 respectively.  In LIB14370, the BACs were 280L6, 315a1, 203k5, and 60g5.  If someone could assign each unitig to the clone, I will rename the files with the clone name and run quiver on the adjusted contig to get the low quality areas to higher quality, and those 6 should be finished !  Let me know if it isn’t clear what I am asking.

Let's first associate the BAC plate IDs with the proper originating locations and libraries.

<a name="bacend"></a>
#### Initial BAC end mapping information

| Pool | Plate ID | Library | Mapped Location (UMD3.1) | IGC | Notes |
| :--- | :--- | :---| :--- | :--- | :--- |
LIB14363 | 389p14 | CH240 | chr18:56,954,654-57,129,335 | chr18| Btau4 coords |
LIB14363 | 170i19 | RP42 | chr4:86,685,432-? |LRC | | 
LIB14363 | 503o23 | CH240 | chr23:28,268,254-28,468,033 | MHC | |
LIB14363 | 49i22 | CH240 | chr5:106,206,907-106,445,457 | NKC |Btau4 coords |
LIB14370 | 280L6 | CH240 | chr18:57,092,237-57,268,067 | chr18 | Btau4 coords |
LIB14370 | 315a1 | CH240 | chrX:?-? | LRC | John Hammond's BAC. |
LIB14370 | 203k5 | CH240 | chr23:28,672,885-28,900,932 | MHC | |
LIB14370 | 60g5 | CH240| chr5:106,309,704-106,537,019 | NKC | Btau4 coords |

Now that I have the initial assignments accounted for, let's align the unitigs to the UMD3 genome. First, let's grab some information on the unitigs using samtools, then I think that I need to write a script to automatically subsection the assembled contig, grep out the sam alignments and try to draw consensus placement across the genome.

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_9

```bash
# Simple line counts
wc -l *.fastq
      4 LIB14363_unitig_0_oriented_vector_trimmed.fastq
      4 LIB14363_unitig_1_oriented_vector_trimmed.fastq
      4 LIB14363_unitig_2_oriented_vector_trimmed.fastq
      4 LIB14363_unitig_3_oriented_vector_trimmed.fastq
      4 LIB14370_unitig_0_oriented_vector_trimmed.fastq
      4 LIB14370_unitig_1_oriented_vector_trimmed.fastq
      4 LIB14370_unitig_2_oriented_vector_trimmed.fastq
     28 total

# So, it's fastq format with one unitig in each file.
# The quality scores are in a range that doesn't correspond to normal fastq phred format. It's Sanger (+33) but extends from 0 to at least 51

# Let's get estimates of the length for each assembled unitig
for i in *.fastq; do echo -ne "$i\t"; perl -e '<>; $s = <>; chomp $s; print length($s); print "\n";' < $i; done
	LIB14363_unitig_0_oriented_vector_trimmed.fastq 205,556
	LIB14363_unitig_1_oriented_vector_trimmed.fastq 175,603
	LIB14363_unitig_2_oriented_vector_trimmed.fastq 207,270
	LIB14363_unitig_3_oriented_vector_trimmed.fastq 190,538
	LIB14370_unitig_0_oriented_vector_trimmed.fastq 240,874
	LIB14370_unitig_1_oriented_vector_trimmed.fastq 230,348
	LIB14370_unitig_2_oriented_vector_trimmed.fastq 213,051

# These unitigs just about correspond to the predicted lengths of the clones based on their end sequence alignments
```

Now, Tim mentioned several regions of low quality within the assembled fragments. In order to hedge my bets, I'm going to take 1kb subsections of the fastq files, align them with BWA mem, then get the consensus position based on the reads.


