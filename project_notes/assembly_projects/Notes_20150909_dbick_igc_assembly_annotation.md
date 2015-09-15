# IGC assembly and annotation work
---
*9/9/2015*

These are my notes on the assembly/association of BAC clones/contigs to their most likely physical map locations.

## Table of Contents
* [First batch unitig assignments](#firstbatch)
	* [Initial BAC end mapping information table](#bacend)
* [Unitig clone association and statistics](#unitigassoc)
	* [First batch unitig association data tables](#unitigtables)
* [Sub alignments to interrogate individual regions](#subaligns)
* [Reassembly without circularization](#reassembly)
	* [New Assembly alignments table](#newassembletable)
* [Subsequent unitig alignments](#subsequent)
	* [Subsequent unitig assignment table](#subsequentassign)
	* [9/15/2015 updated information table](#915summaries)

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

--
*9/10/2015*

I've written [a script](https://github.com/njdbickhart/perl_toolchain/blob/master/assembly_scripts/alignUnitigSectionsToRef.pl) that automates this process and produces condensed alignment information output.

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

<a name="unitigtables"></a>
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

<a name="subaligns"></a>
## Sub alignments to interrogate individual regions

I want to take a look at individual regions to try to understand what's going on. First, let's look at the repetitive region in the chr18 BAC: chr18: 57,566,725-57,571,376 

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_9

```bash
bwa mem ../../../../umd3_data/umd3_kary_unmask_ngap.fa temp.fq > chr18_unitig3_alignment.sam
samtools view -bS chr18_unitig3_alignment.sam | samtools sort -T chr18temp -o chr18_unitig3_alignment.bam -
samtools index chr18_unitig3_alignment.bam

# Only pulling relevant bam information:
samtools view chr18_unitig3_alignment.bam chr18:57566725-57571376 | perl -lane '@nsegs = split(/\./, $F[0]); print "$nsegs[1]\t$nsegs[2]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[11]\t$F[12]";' | sort -k1,1

# Looks like 4 periods of a 6-7kb repetitive element are condensed in the assembly in this region.
```

Let's now pull the overlapping portions of the two NKC unitigs to see how well they align to each other. The shared unitig coordinates should be:  

chr5:99,414,367-99,499,162

I'm going to grep out the coordinates of each unitig (14363, unitig0; 14370, unitig0) that correspond to the shared alignable distance.

```bash
perl -lane 'if($F[3] eq "chr5" && (($F[4] >= 99414367 && $F[4] <= 99499162) || ($F[5] >= 99414367 && $F[5] <= 99499162))){print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < LIB14363_unitig_0_alignment.tab
	83000   93000   chr5    99405090        99414502
	93000   177000  chr5    99428994        99499162
	177000  192000  chr5    99499068        99428858
	192000  206000  chr5    99427566        99415503

perl -lane 'if($F[3] eq "chr5" && (($F[4] >= 99414367 && $F[4] <= 99499162) || ($F[5] >= 99414367 && $F[5] <= 99499162))){print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < LIB14370_unitig_0_alignment.tab
	155000  156000  chr5    99414367        99415367
	156000  171000  chr5    99414973        99512938
	171000  241000  chr5    99510938        99442936

# There appears to be an inversion there between the two unitigs. 
# I'm going to pull the largest amount of sequence and let BWA do the heavy lifting
perl -e '$h = <>; chomp $h; $h =~ s/[\r\@]//g; $s = <>; chomp $s; $s =~ s/\r//g; print "$h\n"; $b = substr($s, 83000, 123000); $b =~ s/(.{60})/$1\n/g; print "$b\n";' < LIB14363_unitig_0_oriented_vector_trimmed.fastq > LIB14363_unitig_0_83000_206000.fa

perl -e '$h = <>; chomp $h; $h =~ s/[\r\@]//g; $s = <>; chomp $s; $s =~ s/\r//g; print ">$h\n"; $b = substr($s, 155000, 86000); $b =~ s/(.{60})/$1\n/g; print "$b\n";' < LIB14370_unitig_0_oriented_vector_trimmed.fastq > LIB14370_unitig_0_155000_241000.fa

# bwa indexing and alignment
bwa index LIB14363_unitig_0_83000_206000.fa
bwa mem LIB14363_unitig_0_83000_206000.fa LIB14370_unitig_0_155000_241000.fa > overlapping_nkc_align.sam

# There was alignment overlap, but the edit distance was large!
	NM:i:3703
	1542M1I597M3I18550M1D3636M2D5143M7D820M5D110M3D15M4D409M3I747M2I11M12I102M1D270M12D295M1D946M6D95M6D7M4I26M1
D1269M1D194M1I101M1I319M6I5M3I112M3D1785M3I5160M4I1623M2I2228M1D1761M37D1305M3D2475M116I24M1D862M1D716M36D431M4I26
1M99D4M11D6M54D10M1I11M3D11M3D2M4D4M3D2M3I18M1D12M3I4M1I9M18D5M7D2M1D4M7I2M5D7M1D10M10I9M4I7M15D7M3D11M1D3M3I3M1D2
M4I10M8D15M3D8M8I10M10D5M15I7M3I8M8I8M23D7M2D3M4I5M2I4M5I10M2D4M9D12M24I7M13I6M4I4M1D5M1I7M2D6M2D8M6I3M8I4M2I5M2I6
M2I10M66I3M9I5M1I3M6D4M4D4M3I3M2D7M6D4M1D13M7I11M9D13M1D9M3I20M13D6M6D7M1I13M2I3M1D8M1D4M2I15M1D12M3D22M28D2M6D8M4
D14M6D2M4D7M19D4M12I4M4I13M2I1M1D5M4D12M5D2M10D10M4D7M5I4M9I5M18I5M4D11M10I3M12D15M6D6M9I2M3D5M4I5M5I2M1I3M17D14M8
D4M8D7M2D17M1I4M2D5M20I7M6D2M24I7M1I3M14I4M4I6M6D5M2I7M8D10M25I3M5D12M2D9M2D2M2D6M1D4M3D11M4D6M6D4M16D4M2D10M8I7M2
4D5M20D14M6I2M4D6M10I4M2D4M4I4M1D4M9D7M6I5M5D4M20I4M8D9M20I11M9I7M1I10M5D9M1D3M8I4M3I6M2I7M1I4M2D6M8D10M1I9M1I10M7
D6M5I11M16D9M23D4M2I2M1D2M4D6M5D16M4D18M4D6M6D7M20I11M5I4M1D14M37I2M7I5M14I10M36I6M4D3M1D6M9D6M2D12M25I14M4D2M2I4M
3I6M3I6M2I7M2D6M17D6M6I4M10D7M4D4M1I7M8D12M7I11M14I3M12D6M2D8M10I3M4D10M2D6M7D4M10D8M4D7M13I6M1D10M2D4M3I4M18D3M14
D3M3D5M10I8M13I6M9I6M11D6M6I21M5I7M14I15M8D3M6I6M4D6M21I8M2D7M8D4M5D8M15I7M19D10M3I6M22I9M11D2M1I3M6D7M1D5M23I10M1
D5M4I6M3D3M1D14M9D6M12I11M1D11M6I9M4I7M6I5M5I4M5D14M7I4M5I5M14D7M4D8M1D27M9D3M21D9M11D10M6I6M1I3M1I2M25D7M19I5M11I
2M3I2M1D10M4I6M2I8M3I9M13D3M2D7M12D7M12D5M17D8M2D16M29D7M4D4M4D10M1I12M7I9M26I247M1D102M3D78M1D17M1D86M14D73M2I115
M2I84M1I599M1D206M2I149M6I1149M1D260M4I30M2I404M1D1323M1I314M1D335M2D600M3I1123M11I50M5D4M1D144M7D1517M1D84M3I330M
1D1688M13I130M1D23M1I668M10I126M1D131M16D219M2D145M5I8M4D274M7D243M3I238M2I553M11D1636M2I37M7I4148M1D651M6I464M24I
5M10I347M2I268M2I433M2I924M4I1213M5I1240M1I2523M975S
```

So, we have overlap, but it's difficult to determine if these deletions/insertions are due to quiver mistakes or actual genetic diversity between heterozygous strands. We'll see after Tim re-quivers!

<a name="reassembly"></a>
## Reassembly without circularization
*9/14/2015*

Tim identified a huge drop in read depth in the previously assembled contigs. He determined that it was due to a circularization step in his automated assembly script that he's previously used to assemble bacterial genomes. Let's check the alignment of these fastqs and then compare them to the previous assemblies to see if there are any noticeable differences. I added a new "longest alignment" printout feature to my [align unitig script](https://github.com/njdbickhart/perl_toolchain/blob/master/assembly_scripts/alignUnitigSectionsToRef.pl) so that I can track problem regions more easily.

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_11

```bash
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f CH240-389p14.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o CH240-389p14.alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr18   57436092        57586436        150344
                        chr11   48330220        48330506        286
                        chr2    95374529        95374606        77

head CH240-389p14.alignment.tab ../tim_assemblies_2015_9_9/LIB14363_unitig_3_alignment.tab
# there are some slight differences, but they're difficult to tell without comparisons. Let's do a filtered diff

cut -f2,3,4,5,6,7 CH240-389p14.alignment.tab > new.temp.tab
cut -f2,3,4,5,6,7 ../tim_assemblies_2015_9_9/LIB14363_unitig_3_alignment.tab > old.temp.tab
diff new.temp.tab old.temp.tab
	7,12c7,18
	< 29000 36000   chr18   57571316        57566725        53;60;60;60;60;60;60
	< 35000 36000   chr4    70667239        70667584        0
	< 35000 36000   chr28   6264630 6264677 0
	< 36000 39000   chr18   57571254        57570254        60;60;60
	< 39000 40000   chr7    24121627        24122627        11
	< 39000 176234  chr18   57565725        57436092        11;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
	---
	> 29000 36000   chr18   57571316        57566716        53;60;60;60;60;60;60
	> 35000 36000   chr11   50634599        50634962        2
	> 35000 36000   chr2    95374529        95374605        0
	> 36000 43000   chr18   57571236        57566716        60;60;38;60;60;60;60;60
	> 43000 44000   chr5    78743182        78744182        11
	> 44000 50000   chr18   57570618        57566725        60;60;60;60;60;60
	> 49000 50000   chr28   6264630 6264677 0
	> 49000 50000   chr11   48330461        48330506        0
	> 50000 54000   chr18   57571556        57566725        16;60;60;60
	> 53000 54000   chr17   73910126        73910424        0
	> 54000 191000  chr18   57571301        57436394        60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
	> 191000        192000  *       0       0       0

# OK, so the newer assembly is ~ 12kb smaller, and has fewer repeats in the center. 
# Also, the shuffling in this region is far more pronounced in the new assembly.
```

Interesting, so the unitig was smaller and the aligned distance is smaller as well. I suspect that this was a consequence of the removal of the circularized region from Tim's script. Let's try aligning 49I22 as Tim sent an obvious RD coverage drop for this clone.

```bash
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f CH240-49i22.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o CH240-49i22.alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr16   5352863 65102421        59749558
                        chr5    99323438        99512899        189461

# The chr5 alignment is the correct alignment -- the chr16 alignments result from two sub-alignments that were further apart on chr16
# Let's do the same diff comparison as before

cut -f2,3,4,5,6,7 CH240-49i22.alignment.tab > new.temp.tab
cut -f2,3,4,5,6,7 ../tim_assemblies_2015_9_9/LIB14363_unitig_0_alignment.tab > old.temp.tab
diff new.temp.tab old.temp.tab
	16,24c16,24
	< 93000 162000  chr5    99428994        99498704        60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;58;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
	< 162000        163000  chr9    42703888        42704888        3
	< 162000        164000  chr5    99498705        99500068        3;54
	< 163000        164000  chr17   12297629        12297806        0
	< 164000        177000  chr5    99499882        99512899        60;60;60;60;60;60;60;60;60;60;60;60;60
	< 176000        177000  chr16   65102305        65102421        0
	< 176000        177000  chr28   36448542        36448582        0
	< 177000        191000  chr5    99428024        99415961        60;60;60;60;60;60;60;60;60;60;60;60;60;60
	< 191000        191012  *       0       12      0
	---
	> 93000 177000  chr5    99428994        99499162        60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
	> 176000        177000  chr12   25972252        25972331        0
	> 177000        178000  chr7    34670721        34671721        4
	> 177000        192000  chr5    99499068        99428858        4;60;60;60;60;60;60;60;60;60;60;60;60;60;60;60
	> 191000        192000  chr16   65102305        65102421        0
	> 191000        192000  chr12   28233552        28233627        27
	> 191000        192000  chr28   36448542        36448582        0
	> 192000        206000  chr5    99427566        99415503        60;60;60;60;60;60;60;60;60;60;60;60;60;60
	> 206000        207000  *       0       0       0

# Looks like the large middle alignment reduced in size in the new unitig. 
# I'm going to have to align WGS reads to this to see which one is better.
# Let's finalize the alignments on the two remaining unitigs. 
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f CH240-503o23.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o CH240-503o23.alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr23   28253676        28445603        191927

perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f RP42-170i19.fastq -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o RP42-170i19.alignment.tab
...
Longest aligments:      chr     start   end     length
                        chrX    96671056        121322701       24651645
                        chr4    86526740        86686044        159304
                        chr3    23081974        23082974        1000
                        chr5    30209424        30210424        1000
                        chr2    93761877        93762160        283
                        chr1    119137314       119137406       92
                        chr14   84089633        84089712        79

# OK, let's create a comparison table for the results
# Unitig lengths first:
for i in *.fastq; do echo -ne "$i\t"; perl -e '<>; $s = <>; chomp $s; print length($s); print "\n";' < $i; done
	CH240-389p14.fastq      176236
	CH240-49i22.fastq       191014
	CH240-503o23.fastq      193890
	RP42-170i19.fastq       161302

```

<a name="newassembletable"></a>
#### New Assembly alignments

| Clone | 9/9 align loc | 9/9 align len | 9/9 unitig len | 9/14 align loc | 9/14 align len | 9/14 unitig len | 
| :--- | :--- | ---: | ---: | :--- | ---: | ---: |
CH240-389p14 | chr18:57,436,394-57,586,436 | 150,042 | 190,538 | chr18:57,436,092-57,586,436 | 150,344 | 176,236
CH240-49i22 | chr5:99,323,438-99,499,162 | 175,724 | 205,556 |  chr5:99,323,438-99,512,899 | 189,461 | 191,014
CH240-503o23 | chr23:28,253,056-28,445,603 | 192,547 | 207,270 | chr23:28,253,676-28,445,603 | 191,927 | 193,890
RP42-170i19 | chr4:86,526,740-86,686,045 | 159,305 | 175,603 | chr4:86,526,740-86,686,044 | 159,304 | 161,302

<a name="subsequent"></a>
## Subsequent unitig alignments
*9/15/2015*

Tim sent me three more unitigs from a batch of four pooled BACs. Here is the information on the BACs from [my previous notes file](https://github.com/njdbickhart/labnotes/blob/master/project_notes/assembly_projects/Notes_20150503_dbick_igc_bac_selection.md):


|Library | Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :--- | :---: | :--- |
| RPCI-42 | chr18 | RPCI42_154M1 | chr18:57,371,161-57,531,510 | Flanking 5' region|
| RPCI-42 | LRC | RPCI42_145P10 | chr18:63,005,919-? | The end sequence might not currently be on this assembly |
| RPCI-42 | MHC | RPCI42_133J13 | chr23:28,321,536-28,448,607 | The 3' end is a split alignment on this assembly |
| RPCI-42 | NKC | RPCI42_154D6 | chr5:99,735,207-99,968,671 | |

Let's get started on assigning the unitig ids to the BAC clones so that Tim can start polishing them. One problem: Tim sent me fasta files this time. I'll update my script to determine file format and process entries appropriately.

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/assemblies/tim_assemblies_2015_9_15

```bash
# unitig 170
# RPCI42_133J13
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14397_unitig_170_Vector_Trim.fasta -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14397_unitig_170_alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr23   28,321,453        28,518,608        197,155
                        chr2    134262461       134263461       1000
                        chr4    112685209       112686209       1000
                        chr13   17090722        17090855        133
                        chr5    38979253        38979313        60

# Interesting, so this is the MHC BAC. The 3' split end of my initial alignment suggested otherwise!

# unitig 174
# RPCI41_154D6
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14397_unitig_174_quiver_Vector_Trim.fasta -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14397_unitig_174_alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr5    99,736,625        99,969,406        232,781
                        chr17   50062891        50063009        118
                        chr2    34746550        34746649        99
                        chr7    98140412        98140479        67

# unitig 1
# RPCI42_154M1
perl ~/share/programs_source/Perl/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f LIB14397_unitig_1_vector_trimmed.fasta -r ../../../../umd3_data/umd3_kary_unmask_ngap.fa -o LIB14397_unitig_1_alignment.tab
...
Longest aligments:      chr     start   end     length
                        chr18   57371117        57532624        161507
                        *       0       1000    1000
                        chr27   22286043        22286293        250
                        chr1    78760548        78760659        111
                        chrX    44952057        44952154        97
                        chr4    3045190 3045260 70

# NOTE: from the chr18 alignment file:
	LIB14397_unitig_1_vector_trimmed        2000    10000   *       0       1000    0;0;0;0;0;0;0;0
# that's a 8kb region that is not on the reference genome. Also the whole thing is rearranged on the genome 

# Let's get the fasta lengths
for i in *.fasta; do echo $i; samtools faidx $i; done
for i in *.fai; do head $i; done
	LIB_14397_unitig_170_Vector_Trim        167,891  35      80      82
	LIB13497_unitig_174|quiver_Vector_Trim  155,313  41      80      82
	LIB14397_unitig_1_vector_trimmed        160,202  35      80      82
```

Let's assign the unitigs to BAC clone IDs.

<a name="subsequentassign"></a>
#### BAC clone assignments

| Fasta file | BAC ID | IGC | Map coordinates | Map Len | Unitig len | Difference| 
| :--- | :--- | :--- | :--- | ---: | ---: | ---: |
LIB14397_unitig_170_Vector_Trim | RPCI42_133J13 | MHC | chr23:28,321,453-28,518,608 | 197,155 | 167,891 | -29,264
LIB14397_unitig_174_Vector_Trim | RPCI41_154D6 | NKC | chr5:99,736,625-99,969,406 | 232,781 | 155,313 | -77,468
LIB14397_unitig_1_vector_trim | RPCI42_154M1 | chr18 | chr18:57,371,117-57,532,624 | 161,507 | 160,202 | -1.305


OK, now I'm going to tabulate all finished unitigs so far just so that I have a picture of which regions of the genome have been covered.

<a name="915summaries"></a>
#### 9/15/2015 current completed clones

| BAC ID | IGC | Map Chr | Map Start | Map End | Map Len | Unitig Len | Difference|
| :--- | :--- | :--- | ---: | ---: | ---: |---: | ---: |
RPCI42_154M1 | chr18 | chr18 | 57,371,117 | 57,532,624 | 161,507 | 160,202 | -1.305
CH240-389p14 | chr18 | chr18 | 57,436,092 | 57,586,436 | 150,344 | 176,236 | 25,892
CH240-49i22 | NKC |chr5 | 99,323,438 | 99,512,899 | 189,461 | 191,014 | 1,553
RPCI41_154D6 | NKC | chr5 | 99,736,625 | 99,969,406 | 232,781 | 155,313 | -77,468
CH240-503o23 | MHC | chr23 | 28,253,676 | 28,445,603 | 191,927 | 193,890 | 1,963
RPCI42_133J13 | MHC | chr23 | 28,321,453 | 28,518,608 | 197,155 | 167,891 | -29,264
RPCI42-170i19 | LRC | chr4 | 86,526,740 | 86,686,044 | 159,304 | 161,302 | 1,998

#### 9/15/2015 current clone fastas/fastqs

| fasta file | BAC ID | IGC |
| :--- | :--- | :---
LIB14397_unitig_170_Vector_Trim.fasta | RPCI42_133J13 | MHC 
CH240-503o23.fastq | CH240-503o23 | MHC
LIB14397_unitig_174_quiver_Vector_Trim.fasta | RPCI41_154D6 | NKC
CH240-49i22.fastq | CH240-49i22 | NKC 
LIB14397_unitig_1_vector_trimmed.fasta | RPCI42_154M1 | chr18 
CH240-389p14.fastq | CH240-389p14 | chr18
RP42-170i19.fastq | RPCI42-170i19 | LRC


