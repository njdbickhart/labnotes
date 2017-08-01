# Buffalo manuscript analysis
---
*2/23/2017*

These are my notes on the generation of SV, MEI and other variant calls on the buffalo sequence data compared to the cattle reference assembly.

## Table of Contents
* [Setting the stage](#stage)
* [Tangram run](#tangram)
* [Melt run](#melt)
* [Yak comparison](#yak)


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

I need to generate a list of mobile element consensus sequences. I'll target 3 to start. My fasta sampling sites are within the Repeatmasker repeat libraries. Specifically here:

> Blade14: ~/RepeatMasker/Libraries/20150807/bos_taurus/*

Let's start by pulling BovB, a good L1 LINE and maybe a SINE.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
mkdir melt
samtools faidx /home/dbickhart/RepeatMasker/Libraries/20150807/bos_taurus/longlib
samtools faidx /home/dbickhart/RepeatMasker/Libraries/20150807/bos_taurus/longlib 'BovB#LINE/RTE-BovB' > melt/BovB.fa

samtools faidx /home/dbickhart/RepeatMasker/Libraries/20150807/bos_taurus/longlib 'L1_BT#LINE/L1' > melt/L1_BT.fa

samtools faidx /home/dbickhart/RepeatMasker/Libraries/20150807/bos_taurus/shortcutlib
samtools faidx /home/dbickhart/RepeatMasker/Libraries/20150807/bos_taurus/shortcutlib 'MIR#SINE/MIR' > melt/mir.fa

# Now I need to generate the bed file for the locations on UMD3
cd melt/
~/RepeatMasker/RepeatMasker -pa 15 -species cow -no_is /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa

perl -e '<>; <>; <>; while(<>){$_ =~ s/^\s+//; chomp; @s = split(/\s+/); if($s[9] eq "L1_BT"){print "$s[4]\t$s[5]\t$s[6]\t$s[9]\n";}}' < /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa.out > umd3_kary.L1_BT.bed
perl -e '<>; <>; <>; while(<>){$_ =~ s/^\s+//; chomp; @s = split(/\s+/); if($s[9] eq "BovB"){print "$s[4]\t$s[5]\t$s[6]\t$s[9]\n";}}' < /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa.out > umd3_kary.BovB.bed
perl -e '<>; <>; <>; while(<>){$_ =~ s/^\s+//; chomp; @s = split(/\s+/); if($s[9] eq "MIR"){print "$s[4]\t$s[5]\t$s[6]\t$s[9]\n";}}' < /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa.out > umd3_kary.MIR.bed

# I also need Bowtie-2 indicies
bowtie2-build /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.bowtie2

# Now to make the zip files
# I'm going to use an error rate of 10 to start with each -- not sure if this is a good error rate!
# Except for the L1 which is an error rate of 3 in the human data they cite
~/jdk1.8.0_05/bin/java -jar ~/MELTv2.0.2/MELT.jar BuildTransposonZIP BovB.fa umd3_kary.BovB.bed BovB 10
~/jdk1.8.0_05/bin/java -jar ~/MELTv2.0.2/MELT.jar BuildTransposonZIP L1_BT.fa umd3_kary.L1_BT.bed L1BT 3
~/jdk1.8.0_05/bin/java -jar ~/MELTv2.0.2/MELT.jar BuildTransposonZIP mir.fa umd3_kary.MIR.bed MIR 10

# I need to make different working directories
mkdir BovB
mkdir L1BT
mkdir MIR

# Preprocessing
for i in ../ITWB*/*.bam; do echo $i; ~/jdk1.8.0_05/bin/java -Xmx2G -jar ~/MELTv2.0.2/MELT.jar Preprocess $i /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa & done

# BovB alignment
for i in `ls ../ITWB*/*merged.bam | grep -v 9`; do echo $i; ~/jdk1.8.0_05/bin/java -Xmx2G -jar ~/MELTv2.0.2/MELT.jar IndivAnalysis -w BovB/ -l $i -c 20 -h /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa -t BovB/BovB_MELT.zip & done

# I got allot of memory overhead errors. Increasing the JVM max memory size
for i in `ls ../ITWB*/*merged.bam | grep -v 9`; do echo $i; ~/jdk1.8.0_05/bin/java -Xmx6G -jar ~/MELTv2.0.2/MELT.jar IndivAnalysis -w BovB/ -l $i -c 20 -h /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa -t BovB/BovB_MELT.zip & done
```

<a name="yak"></a>
## Yak comparison

This is a list of the commands I used to process the Yak files similar to the buffalo. I need to run the following commands on each file individually:

* Jarms 500bp windows
* RAPTR-SV preprocessing
* MEIDivet on the RAPTR-SV divet file
* SAMTOOLS bedcov on raw read counts for cn.mops

And for the group collectively:

* Samtools mpileup
* cn.mops on the bedcov

I will make up two scripts to process these on the cluster so I can generate the information automatically.

#### individual_analysis.sh

```bash
#!/usr/bin/sh
# This is a collection script to automate the processing of each yak bam into useful data
# $1 = bam used in the analysis
# $2 = reference genome used
#SBATCH --nodes=1
#SBATCH --mem=35000
#SBATCH --ntasks-per-node=7

module load java/jdk1.8.0_121
module load samtools

# adding direct link to binary so that RAPTR-SV runs
BINARY=/mnt/nfs/nfs2/bickhart-users/binaries
export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/mrsfast:$PATH

BASE=`basename $1 | cut -d'.' -f1`
# JARMS
mkdir jarms
java -Xmx30g -jar $BINARY/JARMs/store/JaRMS.jar call -i $1 -f $2 -o jarms/${BASE}.jarms -w 500

# RAPTR-SV
mkdir raptr
mkdir temp
java -Xmx34g -jar $BINARY/RAPTR-SV/store/RAPTR-SV.jar preprocess -i $1 -o raptr/${BASE}.raptr -r $2 -g -m 100000000 -p temp -t 6

# MEI-DIVET
mkdir mei
java -Xmx30g -jar $BINARY/MEIDivetID/store/MEIDivetID.jar -i raptr/${BASE}.raptr.preprocess.D.divet -r /mnt/nfs/nfs2/dbickhart/buffalo/gene_data/umd3_ensgene_2kb_upstream.bed -o mei/${BASE}.MEIDivet

# Generating raw read windows
mkdir rdwins
samtools bedcov /mnt/nfs/nfs2/dbickhart/buffalo/gene_data/umd3_naive_1kb_nonovlp_wins.bed $1 > rdwins/${BASE}.rawrd.wins.bed
```

> fry: /mnt/nfs/nfs2/dbickhart/buffalo/yak/yak_bams

```bash
for i in `ls */*.sorted.bam`; do echo $i; sbatch ../../individual_analysis.sh $i /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa ; done

# Whoops! JARMS v9's jar is corrupt! Downgrading and trying again with a shortened script
for i in `ls */*.sorted.bam`; do echo $i; sbatch ../../jarms_sep_run.sh $i /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa ; done

# Running the SNP calling as a big block 
cd snpcalls
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 --wrap="module load samtools; module load bcftools; samtools mpileup -ugf /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa ../SRR3112415/SRR3112415.sorted.bam ../SRR3112417/SRR3112417.sorted.bam ../SRR3112418/SRR3112418.sorted.bam ../SRR3112421/SRR3112421.sorted.bam ../SRR3112425/SRR3112425.sorted.bam ../SRR3112426/SRR3112426.sorted.bam ../SRR3112428/SRR3112428.sorted.bam ../SRR3112431/SRR3112431.sorted.bam ../SRR3112432/SRR3112432.sorted.bam ../SRR3112433/SRR3112433.sorted.bam ../SRR3112434/SRR3112434.sorted.bam ../SRR3112439/SRR3112439.sorted.bam ../SRR3112440/SRR3112440.sorted.bam ../SRR3112441/SRR3112441.sorted.bam | bcftools call -vmO z -o yak_snp_calls.vcf.gz"
cd ..
```

Damn, the RAPTR-SV processing is running into a wall on Steve's cluster. I need to generate just the divets, so let me try that first.

TODO: write samToDivet.pl in the following directory

> fry: /mnt/nfs/nfs2/dbickhart/buffalo


## Generating figures and tables

Figure strategy:
1. Small schema of survey
2. Venn of genes with SVs, MEI, or SNPs upstream
3. Examples of identified SVs, MEI, etc

Table strategy:
1. Input dataset statistics (samples, mapped reads, etc)
2. Number of variants identified 

Supplemental strategy
1. Large SV tables
2. Additional SVPV plots
3. Heatmaps of certain genes for read depth profiles

I am going to try to run with just the buffalo data at this time. My goal is to ID real TFBS and real CNV locations using SVPV for a figure.

SVPV needs SV VCF files as input (annoying). I'll have to generate these in order to run the program. Hopefully I will generate a huge list of pdfs and then have the leisure of scanning them for interesting events.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
# I need to generate a list of chromosome IDs from the BAM headers
samtools view -H ../ITWB15/ITWB15.merged.bam | grep '@SQ' | perl -lane '$F[1] =~ s/SN://g; $F[2] =~ s/LN://g; print "##contig=<ID=$F[1],length=$F[2]>";'

# this is the command I will use to join the different files
cat ../jarms/ITWB10.filtered.dels.bed ../jarms/ITWB11.filtered.dels.bed | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4,5 -o collapse | head

chr1    94001   100000  ITWB10  0.2656999881467868
chr1    113501  121000  ITWB11  0.22925861608387316
chr1    875001  879500  ITWB11,ITWB10   0.1635907305519774,0.13096004668439276
chr1    3727501 3735500 ITWB10  0.4147666729709639
chr1    5256501 5261500 ITWB11,ITWB10   0.10947693270024103,0.1343105166689682

chr15   70302685        .     G       <DUP>   .       PASS    IMPRECISE;SVTYPE=DUP;SVMETHOD=EMBL.DELLYv0.7.3;CHR2=chr15;END=70303402;INSLEN=0;PE=85;MAPQ=60;CT=5to3;CIPOS=-9,9;CIEND=-9,9;RDRATIO=1.26233;AC=2;AN=6     GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-77.3372,0,-245.037:10000:PASS:141:438:177:3:44:17:0:0        0/1:-124.327,0,-244.227:10000:PASS:187:461:155:3:46:25:0:0

perl convert_file_to_vcf.pl -a ../jarms/ITWB10.filtered.dels.bed,../jarms/ITWB10.filtered.dups.bed,../jarms/ITWB11.filtered.dels.bed,../jarms/ITWB11.filtered.dups.bed,../jarms/ITWB12.filtered.dels.bed,../jarms/ITWB12.filtered.dups.bed,../jarms/ITWB13.filtered.dels.bed,../jarms/ITWB13.filtered.dups.bed,../jarms/ITWB14.filtered.dels.bed,../jarms/ITWB14.filtered.dups.bed,../jarms/ITWB15.filtered.dels.bed,../jarms/ITWB15.filtered.dups.bed,../jarms/ITWB1.filtered.dels.bed,../jarms/ITWB1.filtered.dups.bed,../jarms/ITWB2.filtered.dels.bed,../jarms/ITWB2.filtered.dups.bed,../jarms/ITWB3.filtered.dels.bed,../jarms/ITWB3.filtered.dups.bed,../jarms/ITWB4.filtered.dels.bed,../jarms/ITWB4.filtered.dups.bed,../jarms/ITWB5.filtered.dels.bed,../jarms/ITWB5.filtered.dups.bed,../jarms/ITWB6.filtered.dels.bed,../jarms/ITWB6.filtered.dups.bed,../jarms/ITWB7.filtered.dels.bed,../jarms/ITWB7.filtered.dups.bed -b ITWB10,ITWB10,ITWB11,ITWB11,ITWB12,ITWB12,ITWB13,ITWB13,ITWB14,ITWB14,ITWB15,ITWB15,ITWB1,ITWB1,ITWB2,ITWB2,ITWB3,ITWB3,ITWB4,ITWB4,ITWB5,ITWB5,ITWB6,ITWB6,ITWB7,ITWB7 -f /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa -t jarms -o jarms_calls.vcf

# That worked quite well! Now to do the same for cnmops
head ../buffalo_cnmops_cnvrs.tab
seqnames        start   end     width   strand  ITWB1   ITWB10  ITWB11  ITWB12  ITWB13  ITWB14  ITWB15  ITWB2   ITWB3   ITWB4   ITWB5   ITWB6   ITWB7
1       chr1    96001   100001  4001    *       CN2     CN2     CN2     CN2     CN2     CN3     CN2     CN2     CN2     CN2     CN2     CN2     CN2
2       chr1    511001  515001  4001    *       CN2     CN3     CN3     CN3     CN3     CN6     CN3     CN2     CN2     CN2     CN2     CN2     CN2

perl convert_file_to_vcf.pl -a ../buffalo_cnmops_cnvrs.tab -b cnmops -f /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa -t cnmops -o cnmops_calls.vcf

# OK, now I need to try to include the MEI events
# On second thought, it's too difficult and is likely to not coincide with my dup/del calls

wget http://hgdownload.soe.ucsc.edu/goldenPath/bosTau8/database/refGene.txt.gz
gunzip refGene.txt.gz

mkdir plots
for i in `ls ../ITWB*/ITWB*.merged.bam | grep -v 9`; do echo -n "${i},"; done; echo

~/SVPV/SVPV -vcf JaRMS:jarms_calls.vcf,cnMops:cnmops_calls.vcf -o plots -aln ../ITWB10/ITWB10.merged.bam,../ITWB11/ITWB11.merged.bam,../ITWB12/ITWB12.merged.bam,../ITWB13/ITWB13.merged.bam,../ITWB14/ITWB14.merged.bam,../ITWB15/ITWB15.merged.bam,../ITWB1/ITWB1.merged.bam,../ITWB2/ITWB2.merged.bam,../ITWB3/ITWB3.merged.bam,../ITWB4/ITWB4.merged.bam,../ITWB5/ITWB5.merged.bam,../ITWB6/ITWB6.merged.bam,../ITWB7/ITWB7.merged.bam -samples ITWB10,ITWB11,ITWB12,ITWB13,ITWB14,ITWB15,ITWB1,ITWB2,ITWB3,ITWB4,ITWB5,ITWB6,ITWB7 -ref_gene refGene.txt -rgi
# Note: I had to remove the blank line and "start" artifact in cnmops output to get this to work
```

Now I want to identify the number of genes that have each category of variant in the 2kb upstream regions.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
cat buffalo_umd3_comparative_snp.hom.2kb_upstream.bed | cut -f4 | sort | uniq > gene_venns/buffalo_upstream_snp_genes.list
intersectBed -a ../gene_data/umd3_ensgene_2kb_upstream.bed -b buffalo_cnmops_cnvrs.bed | cut -f4 | sort | uniq > gene_venns/buffalo_upstream_cnmops_genes.list
intersectBed -a ../gene_data/umd3_ensgene_2kb_upstream.bed -b jarms/buffalo_combined_cnvrs.bed | cut -f4 | sort | uniq > gene_venns/buffalo_upstream_jarms_genes.list
cp mei_calls/buffalo.group_13_plus.txt gene_venns/buffalo_upstream_mei_genes.list

wc -l *.list
  1163 buffalo_upstream_cnmops_genes.list
   409 buffalo_upstream_jarms_genes.list
   312 buffalo_upstream_mei_genes.list
 24497 buffalo_upstream_snp_genes.list
 26381 total

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl buffalo_upstream_cnmops_genes.list buffalo_upstream_jarms_genes.list buffalo_upstream_mei_genes.list buffalo_upstream_snp_genes.list
	File Number 1: buffalo_upstream_cnmops_genes.list
	File Number 2: buffalo_upstream_jarms_genes.list
	File Number 3: buffalo_upstream_mei_genes.list
	File Number 4: buffalo_upstream_snp_genes.list
	Set     Count
	1       3
	1;2     6
	1;2;4   115
	1;3;4   35
	1;4     1004
	2       32
	2;3;4   3
	2;4     253
	3;4     274
	4       22813

# Interesting that there are no 4-way entries!
# These genes are likely insertions
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 2_3_4 buffalo_upstream_cnmops_genes.list buffalo_upstream_jarms_genes.list buffalo_upstream_mei_genes.list buffalo_upstream_snp_genes.list
	File Number 2: buffalo_upstream_jarms_genes.list
	File Number 3: buffalo_upstream_mei_genes.list
	File Number 4: buffalo_upstream_snp_genes.list
ENSBTAG00000047603.1	<- unknown with 3 exons
ENSBTAG00000045654.1	<- known pseudogene
ENSBTAG00000036098.3	<- unknown with 2 exons

# now to make a venn!
```

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams/gene_venns

```R
# I was getting errors because the SNP list is so large, it could not redesign the edges without making them negative
venn.plot <-draw.quad.venn(area1 = 1163, area2 = 409, area3 = 312, area4 = 1684, n12 = 6, n13 = 0, n123 = 0, n124 = 115, n134 = 35, n14 = 1005, n23 = 0, n234 = 3, n24 = 253, n34 = 274, n1234 = 0, category = c("cn.mops", "JaRMS", "MEI", "SNP"), fill = c("orange", "red", "green", "blue"), lty = "dashed", cex = 2, cat.cex = 2, direct.area = TRUE);
pdf(file="upstream_gene_intersections_buffalo.pdf", useDingbats=FALSE)
grid.draw(venn.plot)
dev.off()

# Didn't work!
# I think I get it: it needs tallies of all of the pairwise comparisons instead of a direct count
venn.plot <-draw.quad.venn(area1 = 1163, area2 = 409, area3 = 312, area4 = 22813, n12 = 121, n13 = 35, n14 = 1154, n23 = 3, n24 = 368, n34 = 309, n123 = 0, n124 = 115, n134 = 35, n234 = 3, n1234 = 0, category = c("cn.mops", "JaRMS", "MEI", "SNP"), fill = c("orange", "red", "green", "blue"), lty = "dashed", cex = 2, cat.cex = 2, cat.col = c("orange", "red", "green", "blue"))
pdf(file="upstream_gene_intersections_buffalo.pdf", useDingbats=FALSE)
grid.draw(venn.plot)
dev.off()
```

#### Table generation

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams/gene_venns

```bash
# Calculating segdups
cd /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams/jarms

perl -lane 'print "$F[1]\t$F[2]\t$F[3]";' < ./buffalo_fixed_segdups.tab | bed_length_sum.pl
        Interval Numbers:       130
        Total Length:           349870
        Length Average:         2691.30769230769
        Length Median:          1249
        Min Length:             1499
        Max Length:             999
        Length Stdev:           7148.79888119586

intersectBed -a /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams/jarms/buffalo_jarms_common_del_regions.bed -b /mnt/iscsi/vnx_gliu_7/100_base_run/jarms/cattle_background_dels_cnvrs.bed -v | bed_length_sum.pl
        Interval Numbers:       5061
        Total Length:           40826500
        Length Average:         8066.8840150168
        Length Median:          2500
        Min Length:             17500
        Max Length:             4000
        Length Stdev:           24375.0320464921

# Seg dup intersections
perl -lane 'print "$F[1]\t$F[2]\t$F[3]";' < ./buffalo_fixed_segdups.tab | intersectBed -a ../../gene_data/umd3_ensgene.bed -b stdin | cut -f4 | sort | uniq | wc -l
30

# Now, the easy stuff where I just copy over the tables I've already generated
mv snpEff_genes.txt /mnt/nfs/nfs2/dbickhart/buffalo/tables/combined_buffalo_snpeff_genes_table.tab
cp buffalo_cnmops_cnvrs.tab /mnt/nfs/nfs2/dbickhart/buffalo/tables/combined_buffalo_cn.mops_detected_cnvrs.tab

# Now I want to make a full annotated Jarms CNVR dataset for deletions with the cn.mops INI significance
cp buffalo_cnmops_inicalls.filtered.merged.bed ./jarms/
cd jarms

~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d genedb.list -i ITWB_dup_file.list -o buffalo_jarms_deletion_combined_cn.mops -t
# The INI significant entries greatly outweigh the count of actual CNVRs. I'm going to trim them.

perl -e '$h = <>; print "$h"; while(<>){chomp; @s = split(/\t/); if($s[4] == 1 && $s[5] =~ /INI_SIG/){next;}else{print join("\t", @s); print "\n";}}' < buffalo_jarms_deletion_combined_cn.mops_regions.tab > buffalo_jarms_deletion_combined_cn.mops_regions.filtered.tab

wc -l buffalo_jarms_deletion_combined_cn.mops_regions.tab buffalo_jarms_deletion_combined_cn.mops_regions.filtered.tab
  360654 buffalo_jarms_deletion_combined_cn.mops_regions.tab
   13445 buffalo_jarms_deletion_combined_cn.mops_regions.filtered.tab
  374099 total

# Now to do the same with the MEI calls
cd mei_calls
# Condensing down my MEIDivet file into flat bed format
for i in ITWB*.filtered.putative.mei; do echo $i; perl -lane 'print "$F[0]\t$F[1]\t$F[2]";' < $i > $i.bed; done
for i in *.filtered.putative.mei.bed; do basename=`echo $i | cut -d'.' -f1`; echo -e "$i\t$basename"; done > mei_file_list.tab
# I need to add the upstream gene region annotations

~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d genedb.list -i mei_file_list.tab -o buffalo_combined_full_MEI_detections -t
cp buffalo_combined_full_MEI_detections_regions.tab /mnt/nfs/nfs2/dbickhart/buffalo/tables/
