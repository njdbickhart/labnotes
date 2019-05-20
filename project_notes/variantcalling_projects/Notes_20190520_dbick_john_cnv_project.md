# John CNV collaboration
---
*5/20/2019*

These are my notes for interrogating CNVs found in John's project.


## Table of Contents


## Looking at the stats

John prepared a lumpy SV callset at this location: 

> Ceres: ~/bostauruscnv/docker/annotated/btcnv.smoove.square.anno.vcf.gz

Let's take a look at the variant sets and see what we can get out of this analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/unmapped_scraping

```bash
# First, some quick stats
gunzip -c ~/bostauruscnv/docker/annotated/btcnv.smoove.square.anno.vcf.gz | grep -v BND | perl -ne 'if($_ =~ /^\#/){next;} @s = split(/\t/); ($end) = $_ =~ /SVLEN=(\-?\d{1,8})/; $end = abs($end); print "$end\n";' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   46498
Sum:    40735777978
Minimum 1
Maximum 99550579
Average 876075.916771
Median  1212
Standard Deviation      6278789.817535
Mode(Highest Distributed Value) 51

# Hmm... that's quite a bit! I suspect that the CNV calls overlap and that I'm not going to get a fair reading on the variability of the cattle genome!

# Let's try to parse the data into a table to find the trends easier.
# I'm going to ditch all of the BND calls because they are very hard to interpret and parse
gunzip -c ~/bostauruscnv/docker/annotated/btcnv.smoove.square.anno.vcf.gz | grep -v BND | perl -ne 'if($_ =~ /^\#/){next;} @s = split(/\t/); ($end) = $_ =~ /SVLEN=(\-?\d{1,8})/; $end = abs($end); ($type) = $_ =~ /SVTYPE=(.{2,4})\;/; print "$s[0]\t$s[1]\t$s[5]\t$type\t$end\n";' > john_cole_sv_regions.tab

# OK, now for a list of the total number of regions and the count of variants per chromosome
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f john_cole_sv_regions.tab -c 3 -d '\t' -m
```

|Entry | Value|
|:-----|-----:|
|DEL   | 23452|
|INV   | 14915|
|DUP   |  8131|

OK, I'm not a big fan of Lumpy's INV calls, but let's get some more stats and then see how much these calls add to the total dataset.

```bash
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f john_cole_sv_regions.tab -c 0 -d '\t' -m
```

|Entry          | Value|
|:--------------|-----:|
|1              |  2855|
|5              |  2403|
|4              |  2233|
|2              |  2188|
|6              |  2159|
|3              |  2145|
|8              |  1941|
|12             |  1908|
|10             |  1877|
|15             |  1875|
|9              |  1802|
|7              |  1747|
|X              |  1730|
|23             |  1700|
|11             |  1544|
|21             |  1530|
|18             |  1461|
|16             |  1326|
|20             |  1326|
|17             |  1195|
|14             |  1192|
|13             |  1107|
|29             |  1089|
|19             |   976|
|24             |   928|
|26             |   871|
|28             |   852|
|27             |   789|
|22             |   668|
|25             |   503|
|Y              |    75|

So the number of calls isn't proportional to the chromosome length, but that could make biological sense. Let's run it through R to generate some plots.

```R
library(dplyr)
library(ggplot2)

data <- read.delim("john_cole_sv_regions.tab", header=FALSE)
colnames(data) <- c("chr", "pos", "qual", "type", "len")
data$chr <- as.factor(data$chr)

# Checking SV regions larger than a megabase
summary(data[data$len > 1000000,])
      chr            pos                 qual          type
 1      : 179   Min.   :    62343   Min.   :    2.0   DEL: 563
 X      : 156   1st Qu.: 14421240   1st Qu.:   13.4   DUP:1219
 10     : 135   Median : 31330588   Median :   54.8   INV: 148
 2      : 118   Mean   : 39198962   Mean   :  196.6
 9      : 101   3rd Qu.: 57275747   3rd Qu.:  160.7
 6      :  99   Max.   :155952417   Max.   :56523.5
 (Other):1142
      len
 Min.   : 1001237
 1st Qu.: 2424221
 Median :10011280
 Mean   :20490794
 3rd Qu.:31913430
 Max.   :99550579

# Removing unplaced contigs
data.filt <- data %>% filter(!startsWith(as.character(chr), "NKL"))

# This suggests to me that 1/8th of the duplication calls are erroneous, and there are a handful of inversion and deletion calls that are pretty bad too
pdf(file="svlen_by_chr.pdf", useDingbats=FALSE)
ggplot(data.filt, aes(x=chr, y=len, fill=chr)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + theme_bw() + scale_y_log10() + labs(title="SV length (Log10) by Chr")
dev.off()

# Now for CNVs by type
pdf(file="svlen_by_type.pdf", useDingbats=FALSE)
ggplot(data.filt, aes(x=type, y=len, fill=type)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + theme_bw() + scale_y_log10() + labs(title="SV length (Log10) by Type")
dev.off()


# Finally, quality (which is going to be partly an analog of allele count) by type
pdf(file="svqual_by_type.pdf", useDingbats=FALSE)
ggplot(data.filt, aes(x=type, y=qual, fill=type)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + theme_bw() + scale_y_log10() + labs(title="SV quality (Log10) by Type")
dev.off()
```