# Buffalo manuscript analysis
---
*2/23/2017*

These are my notes on the generation of SV, MEI and other variant calls on the buffalo sequence data compared to the cattle reference assembly.

## Table of Contents
* [Setting the stage](#stage)
* [Tangram run](#tangram)
* [Melt run](#melt)


## Setting the stage
I need to generate alignments and to start the data analysis. I plan on aligning to UMD3 because of the gene annotation present on that assembly. That will be critical in the identification of MEI and selective sweep variants that are present in the 5'UTR of genes.

Actually, I have everything aligned to UMD3 and can work with the BAMs that are already in existence. Much of this work was in the previous (Notes_20151102_dbick_pan_ruminant_mge_project.md) note file. My new contribution is the generation of CNmops data on buffalo to include with the MEI and SV data.

## Running CN.MOPS to determine I/NI windows

My thoughts are: if we can determine shared windows that are not-informative, then these likely indicate regions of "misassembly" or haplotype divergence in Buffalo vs cattle.

> Blade 14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```R
# Getting a list of BAM files
library("cn.mops")
bamfiles <- list.files(recursive=TRUE, include.dirs=TRUE, pattern="ITWB.+merged.bam$")
bamfiles <- append(bamfiles, c("PC1/PC1.merged.bam"))

# Generating windows and getting read counts
chrnames <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chrX")

bamDataRanges <- lapply(chrnames, function(chr){var <- getReadCountsFromBAM(bamfiles, mode="paired", WL=1000, refSeqName=chr); return(c(var, chr))})

# It crashed
```

OK, I need to generate the windows myself, load them with values, and then load them into R to run cn.mops.

```bash
# Making blind, naive, non overlapping windows
perl -lane 'for(my $x = 1; $x < $F[1]; $x += 1000){my $e = ($x + 1000 > $F[1])? $F[1] : $x + 1000; if($e - $x < 500){next;} print "$F[0]\t$x\t$e";}' < ../../reference/umd3_kary_nmask_hgap.fa.fai > umd3_naive_1kb_nonovlp_wins.bed

# Doing RD counts in those windows
for i in ITWB*/ITWB*.merged.bam; do folder=`echo $i | cut -d'/' -f1`; echo $folder; samtools bedcov umd3_naive_1kb_nonovlp_wins.bed $i > ${folder}/${folder}.1kb.rd.bed; done
``` 

Now to load them into cn.mops! Not sure how to begin here as the data entry is a specific type of GRanges object. I can attempt to load them as a raw data matrix instead.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```R
bedfiles <- list.files(recursive=TRUE, indlue.dirs=TRUE, pattern="ITWB.+.1kb.rd.bed$")
fileIds <- sapply(strsplit(bedfiles, "/"), "[[", 1)

# reading all files in as a list for subsetting later
rdtables <- lapply(bedfiles, function(bed){var <- read.delim(file=bed, sep="\t", header=FALSE, col.names=c("chr", "start", "end", "value")); fname <- sapply(strsplit(bed, "/"), "[[", 1); return(c(var, fname))})

coords <- data.frame(chr=rdtables[[1]]$chr, start=rdtables[[1]]$start, end=rdtables[[1]]$end)
library("GenomicRanges")
# Taking advantage of built-in coercion here
GRangeObj <- as(coords, "GRanges")

rdVals <- data.frame(sapply(rdtables, "[[", 4))
colnames(rdVals) <- sapply(strsplit(bedfiles, "/"), "[[", 1)
for(i in colnames(rdVals)){mcols(GRangeObj)[i] <- rdVals[,i]}

library(cn.mops)
res <- cn.mops(GRangeObj)
res <- calcIntegerCopyNumbers(res)

# Now we want to (a) extract the data into tabular format and (b) get the ini calls
# I can replicate this later for overlapping cnvs to generate read depth plots
# First, the ini calls to tab
# In the future: set "quote=FALSE" in write.table to remove the double quotes
write.table(as.data.frame(iniCall(res)), file="buffalo_inicalls.tab", sep="\t")
write.table(as.data.frame(cnvs(res)), file="buffalo_cnmops_cnvs.tab", sep="\t")
write.table(as.data.frame(cnvr(res)), file="buffalo_cnmops_cnvrs.tab", sep="\t")

# I saved the workspace so I can load this dataset again if needed for plotting
```

A successful run! Now to try to sort out the results. I need to see if any one sample is screwing up the rest and remove it from consideration. My goal is to find genomic regions that are in the following categories:

1. Overlap with previous CNV calls and are INI+
2. Overlap with previous CNV calls and are INI- (likely regions of large sequence divergence)
3. Do not overlap with previous CNV calls and are INI- (likely normal regions in the genome)

```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f buffalo_cnmops_cnvs.tab -c 6
Entry   Count
"ITWB1" 1365
"ITWB10"        1004
"ITWB11"        948
"ITWB12"        2066
"ITWB13"        1024
"ITWB14"        2556
"ITWB15"        899
"ITWB2" 1456
"ITWB3" 1297
"ITWB4" 1243
"ITWB5" 1493
"ITWB6" 2942
"ITWB7" 1273
"ITWB9" 25335		<- This is a problem! It may mean that I need to rerun cn.mops without ITWB9
"median"        1

```
OK, that's an unfortunate setback. Let's reprocess the data quickly.

```R
bedfiles <- bedfiles[! bedfiles %in% "ITWB9/ITWB9.1kb.rd.bed"]
fileIds <- sapply(strsplit(bedfiles, "/"), "[[", 1)

rdtables <- lapply(bedfiles, function(bed){var <- read.delim(file=bed, sep="\t", header=FALSE, col.names=c("chr", "start", "end", "value")); fname <- sapply(strsplit(bed, "/"), "[[", 1); return(c(var, fname))})
library("GenomicRanges")
coords <- data.frame(chr=rdtables[[1]]$chr, start=rdtables[[1]]$start, end=rdtables[[1]]$end)
GRangeObj <- as(coords, "GRanges")

rdVals <- data.frame(sapply(rdtables, "[[", 4)); colnames(rdVals) <- sapply(strsplit(bedfiles, "/"), "[[", 1); for(i in colnames(rdVals)){mcols(GRangeObj)[i] <- rdVals[,i]}

library(cn.mops)
res <- cn.mops(GRangeObj); res <- calcIntegerCopyNumbers(res)

write.table(as.data.frame(iniCall(res)), file="buffalo_inicalls.tab", sep="\t", quote=FALSE)
write.table(as.data.frame(cnvs(res)), file="buffalo_cnmops_cnvs.tab", sep="\t", quote=FALSE)
write.table(as.data.frame(cnvr(res)), file="buffalo_cnmops_cnvrs.tab", sep="\t", quote=FALSE)
# Again, saved the workspace
```

OK, let's begin again!

```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f buffalo_cnmops_cnvs.tab -c 6
Entry   Count
ITWB1   1304
ITWB10  1161
ITWB11  1250
ITWB12  2480
ITWB13  1294
ITWB14  3150
ITWB15  1041
ITWB2   1368
ITWB3   1170
ITWB4   1157
ITWB5   1496
ITWB6   3873
ITWB7   1189
median  1

# That's a big difference! Let's continue to try to pull out the INI call regions
perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\t$s[2]\t$s[3]\t$s[6]\n";}' < buffalo_inicalls.tab > buffalo_cnmops_inicalls.bed
perl -e '<>; while(<>){chomp; @s = split(/\t/); ($c) = $s[9] =~ /CN(\d+)/; print "$s[1]\t$s[2]\t$s[3]\t$s[6]\t$s[7]\t$s[8]\t$c\n";}' < buffalo_cnmops_cnvs.tab > buffalo_cnmops_cnvs.bed
perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\t$s[2]\t$s[3]\n";}' < buffalo_cnmops_cnvrs.tab > buffalo_cnmops_cnvrs.bed

# Now to start the association 
# Making CNVR files
cat jarms/buffalo_deletions_regions.dels.bed jarms/buffalo_fixed_segdups.bed | ~/bin/mergeBed -i stdin > jarms/buffalo_combined_cnvrs.bed

# Simple stats
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/table_bed_length_sum.pl jarms/buffalo_combined_cnvrs.bed buffalo_cnmops_cnvrs.bed
```

FName |  IntNum | TotLen | LenAvg | LenStdev  |      LenMedian   |    SmallestL  |     LargestL
:--- | ---: | ---: | ---: | ---: | ---: | ---: | ---:
jarms/buffalo_combined_cnvrs.bed    |    4549  |  44247951   |     9726.9621894922 |16309.2004859253     |   5999 |   999  |   394999
buffalo_cnmops_cnvrs.bed   |     8548  |  57708909   |     6751.15921853065  |      25260.8368258282  |      4000 |   3000 |   1540000

So cn.mops has smaller events by average but higher deviation. Also, 57 Mbp are copy number variable in the reference using that dataset.

```bash
intersectBed -a jarms/buffalo_combined_cnvrs.bed -b buffalo_cnmops_cnvrs.bed | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1670
        Total Length:           11290322
        Length Average:         6760.67185628742
        Length Median:          4000
        Length Stdev:           12902.5686437736
        Smallest Length:        499
        Largest Length:         253499
intersectBed -a jarms/buffalo_combined_cnvrs.bed -b buffalo_cnmops_cnvrs.bed -r -f 0.5 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       695
        Total Length:           6035712
        Length Average:         8684.47769784173
        Length Median:          4999
        Length Stdev:           16693.3149568773
        Smallest Length:        1500
        Largest Length:         253499

# So not the best agreement. Maybe 30%. Let's see how Jarms calls face up to INI metrics
# INI is a statistic that represents the likelihood that an area is copy number variable
# The higher the INI value, the more likely the region is a true CNV
# I want to test how this statistic is applied to the CNV calls
perl -lane 'if($F[3] >= 0.05){print "$F[0]\t$F[1]\t$F[2]";}' < buffalo_cnmops_inicalls.bed > buffalo_cnmops_inicalls.filtered.bed
mergeBed -i buffalo_cnmops_inicalls.filtered.bed > buffalo_cnmops_inicalls.filtered.merged.bed
# there were some cnmops CNVs missing when I intersected the data.
# Rechecking number distributions against the INI calls
intersectBed -v -a buffalo_cnmops_cnvrs.bed -b buffalo_cnmops_inicalls.filtered.merged.bed -wa | uniq | intersectBed -a stdin -b buffalo_cnmops_inicalls.bed -wb | cut -f7 | statStd.pl
total   67
Minimum 0.0400902829861389
Maximum 0.0498240870973898

# So it look like my first threshold wasn't too far off. I think that 0.04 is the minimum
perl -lane 'if($F[3] >= 0.04){print "$F[0]\t$F[1]\t$F[2]";}' < buffalo_cnmops_inicalls.bed > buffalo_cnmops_inicalls.filtered.bed
mergeBed -i buffalo_cnmops_inicalls.filtered.bed > buffalo_cnmops_inicalls.filtered.merged.bed

# cn.mops INI calls against cnvrs
intersectBed -b buffalo_cnmops_cnvrs.bed -a buffalo_cnmops_inicalls.filtered.merged.bed -wa | uniq | wc -l
11274 / 357,511 total INI regions

# JaRMs calls against significant INI regions
intersectBed -a jarms/buffalo_combined_cnvrs.bed -b buffalo_cnmops_inicalls.filtered.merged.bed -wa | uniq | wc -l
3372 / 4549 = 74.1%

# JaRMs calls in non-significant INI regions
intersectBed -v -a jarms/buffalo_combined_cnvrs.bed -b buffalo_cnmops_inicalls.filtered.merged.bed -wa | uniq > jarms/buffalo_combined_cnvrs.ini.filtered.bed
# Obviously 1177 such regions!
# Duplication intersections here are true seg dups
# Deletion intersections here are likely buffalo-specific deletions

# Segdup check
intersectBed -a jarms/buffalo_fixed_segdups.bed -b jarms/buffalo_combined_cnvrs.ini.filtered.bed | wc -l
22  <- not a reassuring sign. Theyre all small too

# Deletion check
intersectBed -a jarms/buffalo_deletions_regions.gt50perc.bed -b jarms/buffalo_combined_cnvrs.ini.filtered.bed | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       183
        Total Length:           1079817
        Length Average:         5900.63934426229
        Length Median:          4999
        Length Stdev:           6206.91906722537
        Smallest Length:        2499
        Largest Length:         55499
# That's better. We'll check against genes to see if my previous hypotheses hold up
intersectBed -a jarms/buffalo_deletions_regions.gt50perc.bed -b jarms/buffalo_combined_cnvrs.ini.filtered.bed > jarms/buffalo_deletions_regions.ini.filtered.bed
intersectBed -a jarms/buffalo_deletions_regions.ini.filtered.bed -b /mnt/iscsi/vnx_gliu_7/ruminant_project/gene_data/umd3_ensgene.bed -wb | wc -l
42 <- genes with deletions

# Interesting gene to follow up on: ENSBTAG00000035438.3
# Was part of my "5 reciprocally deleted genes" in all buffalo in the mge project
# Is a cn.mops predicted deletion and jarms predicted deletion
# Is apparently INI significant
echo -e "chrX\t143452001\t143490945" | intersectBed -a stdin -b buffalo_cnmops_inicalls.bed -wb
# CNVR: 8521 in cnmops
```

Trying to print this one out just to see.

```R
library(cn.mops)
pdf(file="buffalo_cnmops_shroom2_test_8521.pdf", useDingbats=FALSE)
plot(res, which=8521)
dev.off()
```
<a name="tangram"></a>
## Tangram run

I figure that more than one method is needed to confirm the MEI associations I made previously. 

I need to convert the bams from BWA format to Mosaik first. Then I'll run Tangram to try to ID the MEI.

> /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
for i in */*.merged.bam; do name=`echo $i | cut -d'/' -f1`; echo $name; ~/Tangram/bin/tangram_bam --input $i --ref /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa --output $name/${name}.tangram.bam; samtools index $name/${name}.tangram.bam; done
ITWB10
Segmentation fault (core dumped)
```

This is a known issue on the github repo. Since the repo is dead, other users are suggesting that another tool is used instead. Shortest software run, ever!

<a name="melt"></a>
## Melt run

