# IGC assembly and annotation work
---
*9/9/2015*

These are my notes on the assembly/association of BAC clones/contigs to their most likely physical map locations.

## Table of Contents
* [First batch unitig assignments](#firstbatch)
* [Initial BAC end mapping information table](#bacend)
* [Unitig clone association and statistics](#unitigassoc)

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

I've written a script that automates this process and produces condensed alignment information output.

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_9

```bash
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14363_unitig_0_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14363_unitig_0_alignment.tab

# The script worked beautifully!
# Parsing the output for table form:
perl -lane '$d = shift(@F); $d = substr($d, 0, 17); unshift(@F, $d); print join("|", @F);' < LIB14363_unitig_0_alignment.tab
```

|unitig | unitigstart | unitigend | chr | chrstart | chrend | quality scores|
| :--- | ---: | ---: | :--- | ---: | ---: | :--- |
LIB14363_unitig_0|0|53000|chr5|99323438|99376651|57;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|52000|53000|chr3|23605239|23605314|0
LIB14363_unitig_0|53000|68000|chr5|99376650|99391666|60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|68000|69000|chr1|95672763|95673763|0
LIB14363_unitig_0|68000|69000|chr5|99391666|99391887|0
LIB14363_unitig_0|69000|70000|chr24|24578139|24579139|0
LIB14363_unitig_0|69000|70000|chr16|5352863|5353286|0
LIB14363_unitig_0|69000|70000|chr20|46114211|46114284|0
LIB14363_unitig_0|69000|70000|chr15|48866568|48866668|0
LIB14363_unitig_0|70000|83000|chr5|99391973|99405090|60;60;60;58;60;60;60;60;60;60;60;60;60;4;60;60
LIB14363_unitig_0|82000|83000|chr10|76749422|76749542|2
LIB14363_unitig_0|82000|83000|chr16|62059739|62059839|0
LIB14363_unitig_0|82000|83000|chr13|59875910|59875962|3
LIB14363_unitig_0|83000|93000|chr5|99405090|99414502|60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|92000|93000|chr7|12251476|12251624|8
LIB14363_unitig_0|93000|177000|chr5|99428994|99499162|60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|176000|177000|chr12|25972252|25972331|0
LIB14363_unitig_0|177000|178000|chr7|34670721|34671721|4
LIB14363_unitig_0|177000|192000|chr5|99499068|99428858|4;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|191000|192000|chr16|65102305|65102421|0
LIB14363_unitig_0|191000|192000|chr12|28233552|28233627|27
LIB14363_unitig_0|191000|192000|chr28|36448542|36448582|0
LIB14363_unitig_0|192000|206000|chr5|99427566|99415503|60;60;60;60;60;60;60;60;60;60;60;60;60;60
LIB14363_unitig_0|206000|207000|*|0|0|0

Looking good! It's clear that LIB14363, unitig 0 belongs to the NKC cluster. Also, the alignment length is 175,630 on the UMD3.1 genome, but the unitig's length is 207,000. Looks like a 25-30kb increase in length.

Let's process the rest and continue.

```bash
# 14363, unitig1
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14363_unitig_1_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14363_unitig_1_alignment.tab

# 14363, unitig2
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14363_unitig_2_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14363_unitig_2_alignment.tab

# 14363, unitig3
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14363_unitig_3_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14363_unitig_3_alignment.tab

# 14370, unitig0
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14370_unitig_0_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14370_unitig_0_alignment.tab

# 14370, unitig1
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14370_unitig_1_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14370_unitig_1_alignment.tab

# 14370, unitig2
perl ../../../../programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14370_unitig_2_oriented_vector_trimmed.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14370_unitig_2_alignment.tab
```

Now, let's associate the regions with the table so that Tim knows which unitig belongs to which BAC clone.

<a name="unitigassoc"></a>
## Unitig clone association and statistics

I have to do this by eye, as sub-alignments and misassemblies are going to be difficult to predict and account for in a script. So I will take a look at each alignment.tab file I generated with my script above and then select the best aligned region to incorporate into the following table.

Let's first generate the full, non-redundant aligned region for the predominant chromosome mapping for each unitig. This way, I can estimate the likely "true" length of the region that we've just assembled to show how much of the reference assembly is missing.

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_9

```bash
# 14363, unitig0
# Predominant chromosome is chr5, or the NKC cluster
perl -e '$v = "chr5"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14363_unitig_0_alignment.tab
	chr5    99323438        99499162

# 14363, unitig1
# Predominant chromosome is chr4, or the LRC cluster
perl -e '$v = "chr4"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14363_unitig_1_alignment.tab
	chr4    86526740        86686045

# 14363, unitig2
# Predominant chromosome is chr23, or the MHC cluster
perl -e '$v = "chr23"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14363_unitig_2_alignment.tab
	chr23   28253056        28445603

# 14363, unitig3
# Predominant chromosome is chr18, or the chr18 misassembled region
perl -e '$v = "chr18"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14363_unitig_3_alignment.tab
	chr18   57436394        57586436

# 14370, unitig0
# Predominant chromosome is chr5, or the NKC cluster
perl -e '$v = "chr5"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14370_unitig_0_alignment.tab
	chr5    99414367        99598079

# 14370, unitig1
# Predominant chromosome is chr23, or the MHC cluster
perl -e '$v = "chr23"; $min = 500000000; $max = 0; while(<>){chomp; @s = split(/\t/); if($s[3] eq $v){if($s[4] < $min){$min = $s[4];} if($s[4] > $max){$max = $s[4];} if($s[5] < $min){$min = $s[5];} if($s[5] > $max){$max = $s[5];}}} print "$v\t$min\t$max\n";' < LIB14370_unitig_1_alignment.tab
	chr23   28658242        28879843

# 14370, unitig2
# Lots of alignments to chr18, but it looks misassembled. The middle region aligns to a big chunk of 18, but distributed across the following regions:
chr18   57566101        40333977 
```

I have enough information. Let's put together the final table. I also want to align the NKC pacbio unitig ends to see if the overlapping regions match up well together.

#### Unitig association table.

| Pool | Unitig |Plate ID | Library | BAC Mapped Location | IGC | Notes |
| :--- | :--- | :--- | :---| :--- | :--- | :--- |
LIB14363 | unitig3| 389p14 | CH240 | chr18:56,954,654-57,129,335 | chr18| Btau4 coords |
LIB14363 | unitig1| 170i19 | RP42 | chr4:86,685,432-? |LRC | | 
LIB14363 | unitig2| 503o23 | CH240 | chr23:28,268,254-28,468,033 | MHC | |
LIB14363 | unitig0| 49i22 | CH240 | chr5:106,206,907-106,445,457 | NKC |Btau4 coords |
LIB14370 | part unitig2| 280L6 | CH240 | chr18:57,092,237-57,268,067 | chr18 | Btau4 coords |
LIB14370 | N/A |315a1 | CH240 | chrX:?-? | LRC | John Hammond's BAC. |
LIB14370 | unitig1 |203k5 | CH240 | chr23:28,672,885-28,900,932 | MHC | |
LIB14370 | unitig0 |60g5 | CH240| chr5:106,309,704-106,537,019 | NKC | Btau4 coords |

#### Unitig stats table.


| Pool | Unitig |Plate ID | Library | IGC | Unitig len | Map len | size difference |
| :--- | :--- | :--- | :---| :--- | ---: | ---:| ---:|
LIB14363 | unitig3| 389p14 | CH240 | chr18 |  190,538 | 150,042 | 40,496
LIB14363 | unitig1| 170i19 | RP42 | LRC | 175,603 | 159,305 | 16,298
LIB14363 | unitig2| 503o23 | CH240 | MHC | 207,270 | 192,547 | 14,723
LIB14363 | unitig0| 49i22 | CH240 | NKC | 205,556 | 175,724 | 29,832
LIB14370 | part unitig2| 280L6 | CH240 | chr18 | 213,051 | N/A | N/A
LIB14370 | N/A |315a1 | CH240 | LRC | N/A | N/A | N/A
LIB14370 | unitig1 |203k5 | CH240 | MHC | 230,348 | 221,601 | 8,747
LIB14370 | unitig0 |60g5 | CH240| NKC | 240,874 | 183,712 | 57,162

So, the NKC cluster stands out as the largest gain over the currently assembled region. Interestingly, it should also have two overlapping bac clones from the same library. Let's try to align the regions that should overlap to see what falls out.

