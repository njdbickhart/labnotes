# Buffalo manuscript analysis
---
*2/23/2017*

These are my notes on the generation of SV, MEI and other variant calls on the buffalo sequence data compared to the cattle reference assembly.

## Table of Contents
* [Setting the stage](#stage)


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

```R
bedfiles <- list.files(recursive=TRUE, indlue.dirs=TRUE, pattern="ITWB.+.1kb.rd.bed$")
fileIds <- sapply(strsplit(bedfiles, "/"), "[[", 1)

# reading all files in as a list for subsetting later
rdtables <- lapply(bedfiles, function(bed){var <- read.delim(file=bed, sep="\t", header=FALSE, col.names=c("chr", "start", "end", "value")); fname <- sapply(strsplit(bed, "/"), "[[", 1); return(c(var, fname))})

coords <- data.frame(chr=rdtables[[1]]$chr, start=rdtables[[1]]$start, end=rdtables[[1]]$end)
library("GenomicRanges")
# Taking advantage of built-in coercion here
GRangeObj <- as(coords, "GRanges")

```
