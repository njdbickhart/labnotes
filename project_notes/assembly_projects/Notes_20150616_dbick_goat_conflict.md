# Goat assembly conflict resolution
---
*6/16/2015*

These are my notes on resolving conflicts between the different technologies used to create the Goat reference assembly.

## Table of Contents
* [Overall strategy](#Strategy)
* [Initial example](#Example)
* [Automated resolution](#Automated)

<a name="Strategy"></a>
## Overall strategy

My goals are to get a consensus from four different types of evidence here:

* PacBio sequence
* RH map conversion
* BNG Irys data
* Short read sequence data

I need a really good first example to work with first, so that I can demonstrate the utility of this consensus approach.

<a name="Example"></a>
## Initial example piece

One of the major errors in the Irys-PacBio scaffolding was with utg3, the PacBio contig. It maps to two different Irys contigs. 

Brian found that the region also has different RH mappings. From Brian's email:

> SS3 (SuperScaffold3)
> - The initial ~2.1 Mb of contig utg3 maps to Chr24, while rest (~15 Mb) maps
> to Chr 15. However, this scaffold returns to Chr15 after the mapping error.
> The RH map continues without error on either side of the error section.

Let's pull the divets from the sequence realignment from this region so that we can see if there are transchr reads that map here.

> Blade 14: /mnt/nfs/nfs2/GoatData/testanimal

```bash
# First, let's see how many times utg3 pops up in the search, and what is the disc read profile:
perl -lane 'if($F[1] eq "utg3"){print $_;}' < ERR405776.preprocess.D.divet | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 9 -m
```

|Entry     | Count|
|:---------|-----:|
|deletion  |   814|
|eversion  |  1614|
|insertion | 27841|
|insinv    |   152|
|inversion |    65|
|maxdist   |   182|
|transchr  | 29326|

A decent amount of transchr alignments. Let's try pulling out regions that align only to the breakpoints. I'll try a blunt search first just to see if it works sufficiently.

```bash
perl -lane 'if(($F[1] eq "utg3" && $F[2] < 22000000 && $F[3] > 20000000) || ($F[5] eq "utg3" && $F[6] < 22000000 && $F[7] > 20000000)){print $_;}' < ERR405776.preprocess.D.divet > utg3_problem_region.divet

wc -l utg3_problem_region.divet
	0 utg3_problem_region.divet

# Well, that didn't work. Let's try another brute force check
perl -lane 'if($F[1] eq "utg3" && $F[9] eq "transchr"){print "$F[1]\t$F[2]\t$F[3]\t$F[5]:$F[6]-$F[7]";}elsif($F[5] eq "utg3" && $F[9] eq "transchr"){print "$F[5]\t$F[6]\t$F[7]\t$F[1]:$F[2]-$F[3]";}' < ERR405776.preprocess.D.divet | sortBedFileSTDIN.pl > utg3_problem_region.divet

wc -l utg3_problem_region.divet
	58646 utg3_problem_region.divet
# OK, that generated some meaningful results. Let's visually scan it first.

	utg3    2118268 2118369 utg29351:21251945-21252046
	utg3    2118791 2118892 utg30677:3103738-3103839 <- most obvious pattern from the scan
	utg3    2118791 2118892 utg30677:3103735-3103836
	utg3    2118791 2118892 utg30677:3103738-3103839
	utg3    2119868 2119969 utg42471:5343618-5343719
	utg3    2119913 2120014 utg28094:1744152-1744253
	utg3    2119913 2120014 utg41555:4857215-4857316

# Let's try to generate an RD map in the region.
echo -e 'utg3\t2117768\t2118268\nutg3\t2118268\t2118768\nutg3\t2118768\t2119268\nutg3\t2119268\t2119768' > utg3_problem_region.window.bed

samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam utg3:2116768-2120768 > temp.bam
bedtools coverage -abam temp.bam -b utg3_problem_region.window.bed > utg3_problem_region.RD.bed

# The output had complete coverage everywhere.
# OK, so that's not the region. Let's expand the window map abit.
perl -e 'for($x = 2100000; $x < 2200000; $x += 500){$e = $x + 500; print "utg3\t$x\t$e\n";}' > utg3_problem_region.window.bed
samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam utg3:2100000-2200000 > temp.bam
bedtools coverage -abam temp.bam -b utg3_problem_region.window.bed > utg3_problem_region.RD.bed

perl -lane 'if($F[6] < 0.15){print $_;}' < utg3_problem_region.RD.bed
	utg3    2172500 2173000 2       22      500     0.0440000
	utg3    2173000 2173500 2       12      500     0.0240000
	utg3    2173500 2174000 1       59      500     0.1180000

# Let's expand out a tiny bit...
perl -e 'for($x = 2000000; $x < 2200000; $x += 500){$e = $x + 500; print "utg3\t$x\t$e\n";}' > utg3_problem_region.window.bed
samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam utg3:2000000-2200000 > temp.bam
bedtools coverage -abam temp.bam -b utg3_problem_region.window.bed > utg3_problem_region.RD.bed
perl -lane 'if($F[6] < 0.15){print $_;}' < utg3_problem_region.RD.bed                                             	utg3    2172500 2173000 2       22      500     0.0440000
	utg3    2173000 2173500 2       12      500     0.0240000
	utg3    2173500 2174000 1       59      500     0.1180000

# OK, so the same regions, despite expanding 1 Mb upstream. Let's check the divet file again
	utg3    2172008 2172109 utg8631:10179-10280		<- upstream
	utg3    2172485 2172586 utg49918:11751620-11751721
	utg3    2172590 2172691 utg49918:11838423-11838524
	utg3    2172590 2172691 utg658:184665-184766
	utg3    2172590 2172691 utg30877:4142566-4142667
	utg3    2172590 2172691 utg963:5957169-5957270
	utg3    2172590 2172691 utg49918:11751525-11751626
	utg3    2172820 2172921 utg31066:159815-159916
	utg3    2172820 2172921 utg22:682717-682818
	utg3    2172820 2172921 utg31351:10075109-10075210
	utg3    2172820 2172921 utg51744:1967429-1967530
	utg3    2172820 2172921 utg23407:4394306-4394407
	utg3    2172820 2172921 utg31333:3826784-3826885
	utg3    2172946 2173047 utg50398:2818879-2818980	<- after this, no coverage?
	utg3    2173532 2173633 utg30751:13525611-13525712	<- downstream

# Now to check the bam more closely. It looks like there are just two reads across the entire 1000 bp region
samtools view ERR405776.1.USDA_V3_noheader_sorted.bam utg3:2172300-2174000 | less
	ERR405776.22726001      369     utg3    2172354 0       63H38M  utg87   12548280        0       ATCCATTATGATAAAAACTCTCCAGAAAGCAGGAATAG  ...  NM:i:0  MD:Z:38 AS:i:38 XS:i:38 SA:Z:utg45965,2433800,+,30S71M,0,0;
	ERR405776.109178896     433     utg3    2172354 0       55H46M  utg42030        10757819        0       ATCCATTATGATAAAAACTCTCCAGAAAGCAGGAATAGAAGGAACA  ...  NM:i:0  MD:Z:46 AS:i:46 XS:i:46 SA:Z:utg2841,61362,-,66M35S,0,0;
	ERR405776.205827251     321     utg3    2172354 0       64H37M  utg1311 13962   0       ATCCATTATGATAAAAACTCTCCAGAAAGCAGGAATA   ...   NM:i:0  MD:Z:37 AS:i:37 XS:i:37 SA:Z:utg1704,1031023,+,72M29S,0,0;
	ERR405776.50657621      321     utg3    2172979 0       34M67H  utg1927 3398406 0       TACAAATTCAATGCAATCCCTATCAAGCTACCAT      ...      NM:i:0  MD:Z:34 AS:i:34 XS:i:34 SA:Z:utg50398,2818768,+,31S70M,40,1;
	ERR405776.84342956      433     utg3    2172982 0       31M70H  utg9342 7126    0       AAATTCAATGCAATCCCTATCAAGCTACCAT ... NM:i:0  MD:Z:31 AS:i:31 XS:i:31 SA:Z:utg50398,2818768,-,28S73M,42,1;
	ERR405776.168441381     65      utg3    2173532 0       59M42S  utg30751        13525611        0       GATCTAATTAAAATTAAAAGCTTCTGCACAACAAAGGAAACTATAAGCAAGGTAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTC   ...   NM:i:0  MD:Z:59 AS:i:59 XS:i:58

# they're all split reads
# Last aligning read maps to 2172089 on that contig
```

OK, so that's confirmed. This segment has RH, SEQ and BNG information that suggest that the PacBio contig should be split. Let's see if I can convey that with a figure.

Let's start with a simple histogram for read depth in the region.

```R
library(ggplot2)
data <- read.table("utg3_problem_region.RD.bed", header=FALSE)
histdata <- data.frame(region=paste0(data$V1, ":", data$V2, "-", data$V3), count=data$V4)
pdf(file="rdplot.pdf", useDingbats=FALSE)
plot <- ggplot(histdata, aes(histdata$count)) + geom_histogram(binwidth=2) + theme_bw()
dev.off()

library(GenomicRanges)
gr <- GRanges(seqnames=data$V1, IRanges(start = data$V2, end = data$V3), strand = "+", value= data$V4)
library(ggbio)
pdf(file="rdplot.pdf", useDingbats=FALSE)
autoplot(gr, stat = "coverage", geom = "area") + theme_bw()
dev.off()

library(ggbio)
pdf(file="rdplot.pdf", useDingbats=FALSE)
 #Damn, that didn't work at all
```

*6/17/2015*

--

In order to create a coverage plot, I'm going to have to use R samtools to read the alignments. 

```R
library(Rsamtools)
library(ggbio)
bamfile <- system.file("extdata", "temp.bam", package="Rsamtools")
pdf(file = "rdplot.pdf", useDingbats = FALSE)
bf <- BamFile("temp.bam")
library(GenomicAlignments)
gr <- readGAlignments("temp.bam", "temp.bam.bai", format="BAM", ScanBamParam(what=c("qname")), use.names=FALSE)
qplot(histdata$region, histdata$count, data=histdata, geom="area") + theme_bw()

# After a series of terrible mistakes, this actually worked:
library(ggplot2)
pdf(file = "rdplot.pdf", useDingbats = FALSE)
qplot(histdata$region, histdata$count, data=histdata, geom="density", fill=as.factor(histdata$count)) + theme_bw()
dev.off()
```

Now to make the composite figure, I'm going to have to take the entire contig's RD. I'll reduce it down to 5kb window counts which will make about 3000 windows. It's going to be rough, but I'll deal.

```bash
samtools view -H ERR405776.1.USDA_V3_noheader_sorted.bam | grep 'utg3' | head
	@SQ     SN:utg3 LN:17184534

perl -e 'for($x = 0; $x < 17184534; $x += 5000){$e = $x + 5000; print "utg3\t$x\t$e\n";}' > utg3_windows.bed
samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam utg3:1-17184534 > utg3.bam
samtools index utg3.bam

bedtools coverage -abam utg3.bam -b utg3_windows.bed > utg3_full_cov.RD.bed
cat utg3_full_cov.RD.bed | sortBedFileSTDIN.pl > utg3_full_cov.RD.sorted.bed

cat utg3_full_cov.RD.bed | cut -f4 | statStd.pl                  total   3437
	Minimum 160
	Maximum 1275
	Average 735.458248
	Median  730
	Standard Deviation      105.202653
	Mode(Highest Distributed Value) 769

# Hmmm... 5kb windows were too broad, let's try 3kb instead
perl -e 'for($x = 0; $x < 17184534; $x += 3000){$e = $x + 3000; print "utg3\t$x\t$e\n";}' > utg3_windows.bed
bedtools coverage -abam utg3.bam -b utg3_windows.bed | sortBedFileSTDIN.pl > utg3_full_cov.RD.bed

cat utg3_full_cov.RD.bed | cut -f4 | statStd.pl                  total   5729
	total   5729
	Minimum 16
	Maximum 897
	Average 446.912725
	Median  444
	Standard Deviation      73.443072
	Mode(Highest Distributed Value) 456

# OK, it's R time now
```

I will start off with the simple ggplot figure now in R.

```R
library(ggplot2)
data <- read.table("utg3_full_cov.RD.bed", header=FALSE)
pdf(file = "rdplot.pdf", useDingbats = FALSE)
ggplot(data=data, aes(x=row.names(data), y=V4, group=1)) + geom_line()
dev.off()

# Hmm... the plot is way too large, sadly. Let's focus in on a smaller window
data <- read.table("utg3_problem_region.RD.bed", header=FALSE)
histdata <- data.frame(region=paste0(data$V1, ":", data$V2, "-", data$V3), count=data$V4)
ggplot(data=histdata, aes(x=histdata$region, y=histdata$count, group=1)) + geom_line()
```
<a name="Automated"></a>
## Automated resolution

Ben and Brian have highlighted some target regions that I need to focus on. I want to pull read depth near the breakpoints of these conflicts to try to identify the exact BP regions for Alex to split.

Here are the contigs that Ben found matched up with Brian's issues:

| contig | issue start | issue end | comments |
| :--- | ---: | ---: | :--- |
utg3	| 2,172,189 | 2,173,532 | already fixed
utg6	| 12627251 | 12693378 | single split
utg93	| 468307 | 994847 | single split
utg171	| 505855 | 728596 | single and scaffold split
utg56637 | 487143 | 1386675 | single split 
utg60845 | 1897236 | 2953761 | single and scaffold split
utg512	| 2621300 | 2702353 | single and scaffold split
utg512 | 2927342 | 2994794 | second single and scaffold split
utg523	| 258245 | 829685 | single split
utg1358	| 3286011 | 3359657 | single split
utg23289	| 4374669 | 4687863 | single split
utg23300	| 2382137 | 2802969 | single split
utg29133	| 25654509 | 25678509 | single split
utg29351	| 7521722	| 7572019 | single and scaffold split
utg30737	| 261887	| 1463206	| scaffold split
utg31351	| 645649	| 865449	| scaffold split
utg32260	| 4621055	| 4673465	| single split
utg32281	| 19013441	| 19056981	| single and scaffold split
utg32439	| 2075075	| 2163333	| single split
utg40937	| 5330270 | 5507704	| single split

I'm going to rip them all and push them to a bed file. Then I'll create a shell script designed to automate the task that I did above (albeit without the frustration of working with crappy R libraries!).


> Blade14: /mnt/nfs/nfs2/GoatData/testanimal

```bash
# this is to paste the above table and process it to a bed file
vim test_regions.file
perl -ne '@F = split(/\s+\|\s+/); print "$F[0]\t$F[1]\t$F[2]\n";' < test_regions.file > test_regions.bed

# OK, let's automate this
perl -lane 'my $name = "$F[0]_$F[1]_$F[2]"; 
	print STDERR $name; 
	my $ucsc = "$F[0]:$F[1]-$F[2]"; 
	open(OUT, "> $name.regions.bed"); 
	for($x = $F[1]; $x < $F[2]; $x += 500){
		$e = $x + 500; 
		print OUT "$F[0]\t$x\t$e";
	} 
	close OUT; 
	system("samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam $ucsc > $name.temp.bam"); 
	system("bedtools coverage -abam $name.temp.bam -b $name.regions.bed > $name.coverage.bed");' < test_regions.bed

# Now it's just a matter of grepping out the low coverage regions in these target sites
# First, let's tidy up
rm utg*.temp.bam
rm utg*.regions.bed
mkdir coverage_regions
mv *.coverage.bed ./coverage_regions/

# OK, now to pull regions that have less than 20% coverage of bases within the 500bp windows
cd coverage_regions/
for i in *.bed; do echo $i; perl -lane 'if($F[6] < 0.2){print $_;}' < $i; done
```

#### Obvious examples:

| contig | bp start | bp end | reads |
| :--- | ---: | ---: | :--- |
utg3	| 2,172,189 | 2,173,532 | already fixed
utg23289 | 4473169 | 4478669 | 3
utg30737 | 642387 | 643387 | 0
utg31351 | 673649 | 674149 | 0
utg32260 | 4642055 | 4643055 | 5 (at end)
utg32281 | 19021441 | 19023441 | 2
utg32439 | 2159075 | 2160075 | 1
utg40937 | 5456770 | 5457270 | 0
utg512 | 2667300 | 2668800 | 1
utg56637 | 1183643 | 1184143 | 2
utg56637 | 1321143 | 1322643 | 6
utg60845 | 2864736 | 2865236 | 1
utg93 | 633807 | 634807 | 1

##### Needed closer inspection:

| contig | bp start | bp end | reads | comments |
| :--- | ---: | ---: | :--- | :--- |
utg1358 | 3336011 | 3336511 | 10 | no hard-clipping
utg171 | 562355 | 562855 | 11	| hard-clipping near 562764-562840; utg389,-3281397,25S76M, 
utg171 | 564355 | 565355 | 10 | 10 bp deletion near 564694
utg23300 | 2389137 | 2389637 | 9	| no hard-clipping
utg23300 | 2801137 | 2802137 | 19 | no hard-clipping
utg29133 | 25659509	| 25662009 | 42 | no hard-clipping
utg29351 | 7558222 | 7560222 | 17 | split read mapping: 7559460 with pair mapping to utg29163
utg512 | 2968842 | 2969342 | 25 | disc pair mapping: 2969027 with pair mapping to utg609
utg6 | 12655751 | 12656751 | 24 | no hard-clipping


I am going to err on the side of caution with the "needed closer inspection" contigs and try to split them anyways. Let's create a bed file with the split locations so that I can generate the split fastas for Alex to rescaffold.

I'm going to go to town with a perl one-time script designed to pull the regions from the fasta file using the fai information. Here is the script: [splitFastaWBreakpointBed.pl](https://github.com/njdbickhart/perl_toolchain/blob/master/sequence_data_scripts/splitFastaWBreakpointBed.pl)

```bash
# Creating the bed file for the program
vim breakpoint_regions.file
perl -ne '@F = split(/\s+\|\s+/); print "$F[0]\t$F[1]\t$F[2]\n";' < breakpoint_regions.file > breakpoint_regions.bed
perl ~/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa -b breakpoint_regions.bed -o goat_split_19ctg_assembly.fa

samtools faidx goat_split_19ctg_assembly.fa
gzip goat_split_19ctg_assembly.fa
```

*7/2/2015*

--

Ben and Brian have finalized the list of regions to split. I am just going to double check their work, create a full, final list of split regions, and then create the split fasta for Alex.

Here are the regions from Ben's excel file (GOAT_confilcts_break_contig.xlsx).

| Scaff | start | end | type |
| :--- | ---: | ---: | :--- |
utg443 | 1676022 | 1676027 | yellow small
utg443 | 1730200 | 1730205 | yellow small
utg443 | 1950600 | 1950605 | yellow small
utg751 | 4353500 | 4353505 | red small
utg1746 | 110150 | 112000 | red
utg23187 | 1317750 | 1317755 | blue
utg23285 | 5797500 | 5797505 | blue
utg23407 | 7526750 | 7526755 | yellow small
utg23417 | 119545 | 119550 | yellow small
utg2076 | 688021 | 690011 | red
utg28063 | 129100 | 132700 | red
utg32563 | 5565000 | 5566300 | red
utg3647 | 2006550 | 2006950 | yellow (maybe?)
utg41967 | 231500 | 231505 | blue
utg49095 | 334820 | 336000 | blue

> Blade14: /mnt/nfs/nfs2/GoatData/testanimal

```bash
# Creating a new breakpoints file so that I can keep this analysis separate
vim breakpoint_regions_new.file

# Repeating the automation steps listed previously
perl -ne '@F = split(/\s+\|\s+/); print "$F[0]\t$F[1]\t$F[2]\n";' < breakpoint_regions_new.file > breakpoint_regions_new.bed

perl -lane 'my $name = "$F[0]_$F[1]_$F[2]"; 
	print STDERR $name; 
	my $ucsc = "$F[0]:$F[1]-$F[2]"; 
	open(OUT, "> $name.regions.bed"); 
	for($x = $F[1]; $x < $F[2]; $x += 500){
		$e = $x + 500; 
		print OUT "$F[0]\t$x\t$e";
	} 
	close OUT; 
	system("samtools view -b ERR405776.1.USDA_V3_noheader_sorted.bam $ucsc > $name.temp.bam"); 
	system("bedtools coverage -abam $name.temp.bam -b $name.regions.bed > $name.coverage.bed");' < breakpoint_regions_new.bed

rm utg*.temp.bam
rm utg*.regions.bed

# Now to pull the regions that have less than 20% coverage
for i in *coverage.bed; do echo $i; perl -lane 'if($F[6] < 0.2){print $_;}' < $i; done
```

**Refined coordinates**

| Scaff | start | end | type | confirmed |
| :--- | ---: | ---: | :--- | :--- | 
utg443 | 1676022 | 1676027 | yellow small | REMOVE
utg443 | 1730200 | 1730205 | yellow small | confirm
utg443 | 1950600 | 1950605 | yellow small | confirm
utg751 | 4353500 | 4354000 | red small | confirm
utg1746 | 111736 | 112118 | red	| confirm
utg23187 | 1317750 | 1317755 | blue | confirm
utg23285 | 5797000 | 5798380 | blue  | confirm
utg23407 | 7526750 | 7527250 | yellow small | confirm
utg23417 | 119400 | 119550 | yellow small | confirm
utg2076 | 688521 | 690011 | red	| confirm
utg28063 | 129600 | 133100 | red | confirm
utg32563 | 5565000 | 5566300 | red | confirm
utg3647 | 2006266 | 2007232 | yellow (maybe?) | confirm
utg41967 | 231290 | 231592 | blue | confirm	
utg49095 | 334820 | 336000 | blue | confirm

```bash
# I resaved the breakpoint_regions_new.file file with the updated coordinates above
perl -ne '@F = split(/\s+\|\s+/); print "$F[0]\t$F[1]\t$F[2]\n";' < breakpoint_regions_new.file > breakpoint_regions_saved.bed

# Now to combine all of the regions, rerun my script, and then gzip the resulting fasta file
cat breakpoint_regions_saved.bed breakpoint_regions_saved.bed > final_ctg_breakpoints.bed
perl ~/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa -b final_ctg_breakpoints.bed -o goat_split_28ctg_assembly.fa

# Just making sure that everything is in order...
samtools faidx goat_split_28ctg_assembly.fa
# Nope! It's not! There was a problem with the two breakpoint files that I think have to do with carriage returns
dos2unix breakpoint_regions_saved.bed
dos2unix breakpoint_regions.bed
cat breakpoint_regions_saved.bed breakpoint_regions.bed > final_ctg_breakpoints.bed

# Let's try it again
perl ~/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f ../Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa -b final_ctg_breakpoints.bed -o goat_split_36ctg_assembly.fa
samtools faidx goat_split_36ctg_assembly.fa
gzip goat_split_36ctg_assembly.fa

wc -l goat_split_36ctg_assembly.fa.fai
	3110 goat_split_36ctg_assembly.fa.fai
# 36 split events + 3074 = 3110, so we're good.
```

