# Rumen Metagenome manuscript Reviewer Responses
---
*4/11/2019*

My previous note files were getting too large, so I moved this to a new file given the number of calculations I had to make.

## Table of Contents

## Starting out

Here are my notes for generating statistics and data to respond to the reviewer comments on our manuscript.

> Ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool

```bash
# Getting the Bowers standard counts:
# High-quality draft (HQ)
perl -lane 'if($F[-2] > 90 && $F[-1] < 5){print $F_;}' < pacbio_final_dastool_DASTool_summary.txt | wc -l
10
# Medium-quality draft (MQ)
perl -lane 'if($F[-2] > 50 && $F[-1] < 10){print $F_;}' < pacbio_final_dastool_DASTool_summary.txt | wc -l
103

# High-quality (HQ)
perl -lane 'if($F[-2] > 90 && $F[-1] < 5){print $F_;}' < illumina_final_dastool_DASTool_summary.txt | wc -l
42
# Medium-quality (MQ)
perl -lane 'if($F[-2] > 50 && $F[-1] < 10){print $F_;}' < illumina_final_dastool_DASTool_summary.txt | wc -l
325

# Now, printing them out so that I can generate new master tables
# Illumina
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 5 || $s[-2] < 90){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa");while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_DASTool_HQD.tab
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 10 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa");while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_DASTool_MQD.tab

# Pacbio
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 5 || $s[-2] < 90){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_final_dastool_DASTool_HQD.tab
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 10 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_final_dastool_DASTool_MQD.tab
```

#### Generating new master tables and calculating stats

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Editting the json files
cp illumina_data_files_11_2018.json illumina_data_files_4_2019.json
cp pacbio_data_files_9_2018.json pacbio_data_files_4_2019.json

python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j illumina_data_files_4_2019.json -t illumina_preliminary_blobtools_table.tab -o illumina_megahit_master_table_4_1_2019.tab
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j pacbio_data_files_4_2019.json -t pacbio_preliminary_blobtools_table.tab -o pacbio_final_pilon_master_table_4_1_2019.tab

# Subsectioning them for ease of use in my other analysis pipelines
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < illumina_megahit_master_table_4_1_2019.tab > illumina_megahit_master_table_4_1_2019.MQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < illumina_megahit_master_table_4_1_2019.tab > illumina_megahit_master_table_4_1_2019.HQbins.short.tab

# Full dataset subsectioning
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,9,22,23,27,31,35,39,43,47,48,49,50,51,52,53,56,57,58]) . "\n";' < illumina_megahit_master_table_4_1_2019.tab > illumina_megahit_master_table_4_1_2019.short.tab
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,9,22,23,27,31,35,39,43,47,48,49,50,51,52,53,56,57,58]) . "\n";' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_master_table_4_1_2019.short.tab
```

#### Calculating new stats for the manuscript

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
# 16S stats
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab -c 0 -l mick_16s_pacbio_summary.contig.list | perl -ne 'chomp; @F = split(/\t/); print "$F[11]\t$F[0]\t$F[4]\t$F[10]\t$F[12]\n";' > mick_16s_pacbio_summary.reformat.MQbins.tab
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab -c 0 -l mick_16s_pacbio_summary.contig.list | perl -ne 'chomp; @F = split(/\t/); print "$F[11]\t$F[0]\t$F[4]\t$F[10]\t$F[12]\n";' > mick_16s_pacbio_summary.reformat.HQbins.tab

wc -l mick_16s_pacbio_summary.reformat.MQbins.tab mick_16s_pacbio_summary.reformat.HQbins.tab               64 mick_16s_pacbio_summary.reformat.MQbins.tab
  14 mick_16s_pacbio_summary.reformat.HQbins.tab
  78 total

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f mick_16s_pacbio_summary.reformat.MQbins.tab -c 2 -d '\t'
Entry   Value
Bacteria        59
Eukaryota       4
Archaea 1
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Dataset comparison counts
perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   22
Sum:    317319
Minimum 3552
Maximum 30215
Average 14423.590909
Median  13725
Standard Deviation      7221.517982

perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < illumina_megahit_master_table_4_1_2019.HQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3650
Sum:    25491578
Minimum 1000
Maximum 90246
Average 6983.993973
Median  4899.5
Standard Deviation      6453.684480

perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   1254
Sum:    16151563
Minimum 1508
Maximum 367662
Average 12880.034290
Median  9669
Standard Deviation      14777.483832

# The "does not align to any other dataset" stats
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[1]\n";}' < illumina_megahit_master_table_4_1_2019.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   27120
Sum:    87817564
Minimum 1000
Maximum 71737
Average 3238.110767
Median  2455
Standard Deviation      2838.136951

perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[1]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   20
Sum:    137081
Minimum 1762
Maximum 16216
Average 6854.050000
Median  7242.5
Standard Deviation      4053.889237

# Now to find out if we still have firmicutes in the PACBIO dataset
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[0]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_MQ_unique_ctgs.list
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[0]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_MQ_unique_ctgs.list

# Still Firmicutes and the statistic still stands
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab -c 0 -d '\t' -l pacbio_final_pilon_MQ_unique_ctgs.list | perl -lane 'print $F[2];' > pacbio_final_pilon_MQ_unique_ctgs.list.gc
perl -lane 'if($F[0] eq "name"){next;}else{print $F[2];}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.short.gc

##### Generating new supplementary table for HQ bins for Reviewer #2
# First, let's calculate the total coverage for all datasets from the short-reads
# The coverage values are averages, so I need to multiply by the lengths to get an estimate
perl -e '$sum = 0; $sa = 0; $n = 0; while(<>){chomp; @F = split(/\t/); if($F[0] eq "name"){next;} $sum += ($F[9] * $F[1]); $sa += $F[9]; $n++}; $avg = $sa / $n; print "$sum\t$sa\t$n\t$avg\n";' < pacbio_final_pilon_master_table_4_1_2019.tab
66888601387.3526        4652765.74600004        77670   59.9042840993954
perl -e '$sum = 0; $sa = 0; $n = 0; while(<>){chomp; @F = split(/\t/); if($F[0] eq "name"){next;} $sum += ($F[9] * $F[1]); $sa += $F[9]; $n++}; $avg = $sa / $n; print "$sum\t$sa\t$n\t$avg\n";' < illumina_megahit_master_table_4_1_2019.tab
129991919911.382        54111883.2900027        2182263 24.7962245109791

# OK, it's a slight overestimate, but still in the ballpark
# Let's generate a "proportion of aligned reads" statistic per HQ bin

perl consolidate_HQ_bin_table.pl illumina_megahit_master_table_4_1_2019.tab 129991919911 illumina_megahit_HQBIN_supplementary_summary.tab
perl consolidate_HQ_bin_table.pl pacbio_final_pilon_master_table_4_1_2019.tab 66888601387 pacbio_final_HQBIN_supplementary_summary.tab

#OK so most of the data worked out great, but the proportional coverage is bad. Damn. I'm going to have to calculate this by hand from the bams!

# Generating prodigal supplementary tables for the MQ bins
perl -lane 'if($F[0] eq "name"){next;}else{ print $F[0];}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.ctg.list
perl -lane 'if($F[0] eq "name"){next;}else{ print $F[0];}' < illumina_megahit_master_table_4_1_2019.MQbins.short.tab > illumina_megahit_master_table_4_1_2019.MQbins.ctg.list

# Generating the tables from the MQ dataset
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} close IN; print "ContigID\tStart\tEnd\tOrient\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\n"; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/(\_\d{1,4}$)//; if(exists($h{$s[0]})){$s[0] .= $1; print join("\t", @s) . "\n";}} close IN;' pacbio_final_pilon_master_table_4_1_2019.MQbins.ctg.list supplementary_table_10_long_read_prodigal.tab > supplementary_table_10_long_read_prodigal.revis.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} close IN; print "ContigID\tStart\tEnd\tOrient\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\n"; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/(\_\d{1,4}$)//; if(exists($h{$s[0]})){$s[0] .= $1; print join("\t", @s) . "\n";}} close IN;' illumina_megahit_master_table_4_1_2019.MQbins.ctg.list supplementary_table_9_short_read_prodigal.tab > supplementary_table_9_short_read_prodigal.revis.tab

# Now to regenerate the average counts and retest them for significance
perl -e 'print "Contig\tLength\tComplete\tPartial\tTech\n"; while(<>){chomp; @s = split(/\t/); if($s[0] eq "name"){next;}else{print "$s[0]\t$s[1]\t$s[14]\t$s[15]\tPacBio\n";}}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.orfs.tab
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "name"){next;}else{print "$s[0]\t$s[1]\t$s[14]\t$s[15]\tIllumina\n";}}' < illumina_megahit_master_table_4_1_2019.MQbins.short.tab > illumina_megahit_master_table_4_1_2019.MQbins.orfs.tab

cat pacbio_final_pilon_master_table_4_1_2019.MQbins.orfs.tab illumina_megahit_master_table_4_1_2019.MQbins.orfs.tab > combined_4_1_2019_MQbins.orfs.tab
```

```R
library(dplyr)
orfdata <- read.delim("combined_4_1_2019_MQbins.orfs.tab", header=TRUE)

summary(orfdata[orfdata$Tech == "PacBio", "Complete"])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   0.00    9.00   15.00   22.35   26.00  663.00
summary(orfdata[orfdata$Tech == "PacBio", "Partial"])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  0.000   1.000   1.000   1.269   2.000   2.000
summary(orfdata[orfdata$Tech == "Illumina", "Complete"])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  0.000   1.000   2.000   3.752   4.000 238.000
summary(orfdata[orfdata$Tech == "Illumina", "Partial"])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  0.000   1.000   2.000   1.499   2.000   2.000

# Now I want to test averages within contig length quantiles
orfdata <- within(orfdata, quantile <- as.integer(cut(Length, quantile(Length, probs=0:5/5), include.lowest = TRUE)))
for(x in c(1,2,3,4,5)){print(ks.test(orfdata[orfdata$quantile == x & orfdata$Tech == "Illumina","Complete"], orfdata[orfdata$quantile == x & orfdata$Tech == "PacBio","Complete"], alternative="greater"))}

# All except for the smallest contigs passed the test.
orfdata %>% group_by(quantile, Tech) %>% summarize(mean = mean(Complete), num = n(), sum=sum(Complete))
# A tibble: 10 x 5
# Groups:   quantile [5]
   quantile Tech       mean   num    sum
      <int> <fct>     <dbl> <int>  <int>
 1        1 Illumina  0.465 20555   9556
 2        1 PacBio    1.47     17     25

 3        2 Illumina  1.19  20434  24318
 4        2 PacBio    2.48    137    340

 5        3 Illumina  2.29  20243  46261
 6        3 PacBio    3.98    337   1342

 7        4 Illumina  4.35  19220  83654
 8        4 PacBio    7.26   1344   9751

 9        5 Illumina 13.2   14568 192679
10        5 PacBio   27.3    6001 163703


######
## Recalculating table 2
######

illumina <- read.delim("illumina_megahit_master_table_4_1_2019.short.tab", header=TRUE)
illumina %>% mutate(cat = ifelse(DASBin_HQ != "NOBIN", "HQ", ifelse(DASBin_MQ != "NOBIN", "MQ", "UNBIN"))) %>% group_by(cat) %>% summarize(mean = mean(CompleteORFs))
1 HQ     7.66
2 MQ     3.39
3 UNBIN  1.31

illumina %>% mutate(cat = ifelse(DASBin_HQ != "NOBIN", "HQ", ifelse(DASBin_MQ != "NOBIN", "MQ", "UNBIN"))) %>% group_by(cat, superkingdom.t.24) %>% summarize(count = sum(as.numeric(length)))

pacbio <- read.delim("pacbio_final_pilon_master_table_4_1_2019.short.tab", header=TRUE)
pacbio %>% mutate(cat = ifelse(DASBin_HQ != "NOBIN", "HQ", ifelse(DASBin_MQ != "NOBIN", "MQ", "UNBIN"))) %>% group_by(cat) %>% summarize(mean = mean(CompleteORFs))
1 HQ     48.2
2 MQ     20.8
3 UNBIN  14.6

pacbio %>% mutate(cat = ifelse(DASBin_HQ != "NOBIN", "HQ", ifelse(DASBin_MQ != "NOBIN", "MQ", "UNBIN"))) %>% group_by(cat, superkingdom.t.24) %>% summarize(count = sum(as.numeric(length)))
```

#### Checking partial ORF counts near contig ends

I need to redo the prodigal short form data. Then, I'm going to apply the quantile approach to segregate the datasets as per above. I also want to get the average contig coverage data integrated to see if contig coverage plays a role in the number of completed ORFs.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
# Generating the input
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[0] ne "name"){print "$F[0]\t$F[1]\t$F[9]\n";}' < ../master_tables/illumina_megahit_master_table_4_1_2019.tab > illumina_MQ_bins_contiglens.tab
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[0] ne "name"){print "$F[0]\t$F[1]\t$F[9]\n";}' < ../master_tables/pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_MQ_bins_contiglens.tab

perl count_partial_ORFs_near_contig_ends.pl -s illumina_megahit_prodigal_proteins.shortform.tab -f illumina_MQ_bins_contiglens.tab -o  illumina_MQ_bins_orfs
perl count_partial_ORFs_near_contig_ends.pl -s pacbio_final_prodigal_proteins.shortform.tab -f pacbio_MQ_bins_contiglens.tab -o pacbio_MQ_bins_orfs
```

Now to do the analysis...

```R
library(dplyr)
library(ggplot2)
illmq <- read.delim("illumina_MQ_bins_orfs.partial.cats", header=TRUE)

table(illmq$Location, illmq$Completion)

          Full MissBoth MissEnd MissStart
  End     3664        0   68705         0
  Mid   349240        0       0         0
  Start   3564     4238     289     69202

# Total 142434 partial orfs
# 96.8% occur within start or end o

chisq.test(table(illmq$Location, illmq$Completion))

        Pearson's Chi-squared test

data:  table(illmq$Location, illmq$Completion)
X-squared = 936070, df = 6, p-value < 2.2e-16

pbmq <- read.delim("pacbio_MQ_bins_orfs.partial.cats", header=TRUE)

table(pbmq$Location, pbmq$Completion)

          Full MissBoth MissEnd MissStart
  End      565        0    5007         0
  Mid   174036        0       0         0
  Start    560        3       0      4934

chisq.test(table(pbmq$Location, pbmq$Completion))

        Pearson's Chi-squared test

data:  table(pbmq$Location, pbmq$Completion)
X-squared = 331520, df = 6, p-value < 2.2e-16

illmq$Tech <- c("Illumina")
pbmq$Tech <- c("PacBio")

combmq <- bind_rows(illmq, pbmq)

summary(lm(illmq$OrfLen ~ illmq$CtgCov))

Call:
lm(formula = illmq$OrfLen ~ illmq$CtgCov)

Residuals:
   Min     1Q Median     3Q    Max
 -1127   -466   -140    287  36173

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  8.965e+02  9.914e-01 904.332  < 2e-16 ***
illmq$CtgCov 2.130e-02  7.594e-03   2.805  0.00503 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 679.9 on 498900 degrees of freedom
Multiple R-squared:  1.577e-05, Adjusted R-squared:  1.377e-05
F-statistic:  7.87 on 1 and 498900 DF,  p-value: 0.005026

summary(lm(pbmq$OrfLen ~ pbmq$CtgCov))

Call:
lm(formula = pbmq$OrfLen ~ pbmq$CtgCov)

Residuals:
    Min      1Q  Median      3Q     Max
-2836.2  -447.5  -165.1   264.5 12331.9

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 782.91845    2.12745  368.01   <2e-16 ***
pbmq$CtgCov   0.59307    0.02424   24.47   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 643.9 on 185103 degrees of freedom
Multiple R-squared:  0.003224,  Adjusted R-squared:  0.003219
F-statistic: 598.7 on 1 and 185103 DF,  p-value: < 2.2e-16

cor.test(pbmq$CtgCov, pbmq$OrfLen)

        Pearson's product-moment correlation

data:  pbmq$CtgCov and pbmq$OrfLen
t = 24.469, df = 185100, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.0522386 0.0613203
sample estimates:
       cor
0.05678063

# So there is an effect on ORFLen, but the OrfLen isn't correlated with coverage. Let's plot it though
pdf("pacbio_ctgcov_by_orflen.pdf", useDingbats=FALSE)
ggplot(pbmq, aes(x=CtgCov, y=OrfLen)) + stat_summary(fun.data = mean_se, geom = "errorbar") + geom_smooth(method='lm')
dev.off()

# A logit regression shows an influence based on location of the ORF
summary(glm(pbmq$Location ~ pbmq$OrfLen, family=binomial(link="logit")))

Call:
glm(formula = pbmq$Location ~ pbmq$OrfLen, family = binomial(link = "logit"))

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-4.0567   0.2190   0.2505   0.2729   0.2993

Coefficients:
             Estimate Std. Error z value Pr(>|z|)
(Intercept) 3.049e+00  2.322e-02  131.31   <2e-16 ***
pbmq$OrfLen 5.816e-04  2.875e-05   20.23   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 50014  on 185104  degrees of freedom
Residual deviance: 49521  on 185103  degrees of freedom
AIC: 49525

Number of Fisher Scoring iterations: 6

# It's probably the different classes of ORFs that are to blame
combmq %>% group_by(Location, Tech) %>% summarize(mean = mean(OrfLen), num = n(), sum=sum(OrfLen), median=median(OrfLen), stdev = sd(OrfLen))
# A tibble: 6 x 7
# Groups:   Location [3]
  Location Tech      mean    num       sum median stdev
  <fct>    <chr>    <dbl>  <int>     <int>  <dbl> <dbl>
1 End      Illumina  744.  72369  53868969    584  662.
2 End      PacBio    651.   5572   3625196    470  624.
3 Mid      Illumina  950. 349240 331923038    815  665.
4 Mid      PacBio    831. 174036 144650616    668  645.
5 Start    Illumina  800.  77293  61818935    626  728.
6 Start    PacBio    636.   5497   3494252    458  606.
```


## Creating read subsamples for Phase team

We're testing the hypothesis that PE links are more informative than Hi-C data in binning metagenome assembly contigs. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/pilot_project/illumina

```bash
# I think that the files need to be unzipped for this, unfortunately!
for i in YMPrepCannula_S1_*.fastq.gz; do sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p msn --wrap="gunzip $i"; done

sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq.gz -s YMPrepCannula_S1_L001_R2_001.fastq.gz -f YMPrepCannula_S1_L002_R1_001.fastq.gz -s YMPrepCannula_S1_L002_R2_001.fastq.gz -f YMPrepCannula_S1_L003_R1_001.fastq.gz -s YMPrepCannula_S1_L003_R2_001.fastq.gz -f YMPrepCannula_S1_L004_R1_001.fastq.gz -s YMPrepCannula_S1_L004_R2_001.fastq.gz -o illumina_pe_reads_m1_subsample -l 22487509"

# Ooops! Looks like gziping doesn't work for file line counts
sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq -s YMPrepCannula_S1_L001_R2_001.fastq -f YMPrepCannula_S1_L002_R1_001.fastq -s YMPrepCannula_S1_L002_R2_001.fastq -f YMPrepCannula_S1_L003_R1_001.fastq -s YMPrepCannula_S1_L003_R2_001.fastq -f YMPrepCannula_S1_L004_R1_001.fastq -s YMPrepCannula_S1_L004_R2_001.fastq -o illumina_pe_reads_m1_subsample -l 22487509 -i subsample_m1.log -b 2280780244"
sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq -s YMPrepCannula_S1_L001_R2_001.fastq -f YMPrepCannula_S1_L002_R1_001.fastq -s YMPrepCannula_S1_L002_R2_001.fastq -f YMPrepCannula_S1_L003_R1_001.fastq -s YMPrepCannula_S1_L003_R2_001.fastq -f YMPrepCannula_S1_L004_R1_001.fastq -s YMPrepCannula_S1_L004_R2_001.fastq -o illumina_pe_reads_s3_subsample -l 40889499 -i subsample_s3.log -b 2280780244"
```

I sent the files off to Max for analysis. We'll see what he finds!

## Generating log differential coverage plots

The first reviewer wants to see log differential coverage plots. I had generated these with Blobplots, but I think he wants us to use his pipeline. First, let's align the Hi-C reads and error corrected reads to the HQ bins I just prepared.

Given the fact that MMGenome2 can plot a grid of differential coverage, let's see if we can plot ALL datasets (Hi-C, PacbioEC, normal short-read, 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/diff_cov

```bash
module load bwa samtools metabat/2.12.1
# Generating the lists I need
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print join("\t", @F[(0,1,2,9,14,15,16,17,21,8,23,27,39,51,52,53)]) . "\n";}' < ../master_tables/illumina_megahit_master_table_4_1_2019.tab > illumina_HQBIN_only_shortcov.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print join("\t", @F[(0,1,2,9,14,15,16,17,21,8,23,27,39,51,52,53)]) . "\n";}' < ../master_tables/pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_HQBIN_only_shortcov.tab

# Generating the HQ bin fastas
perl -lane 'if($F[0] eq "name"){next;}else{print $F[0];}' < illumina_HQBIN_only_shortcov.tab > illumina_HQBIN_only_shortcov.ctg.list
samtools faidx -r illumina_HQBIN_only_shortcov.ctg.list ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa > illumina_megahit_final_onlyHQBIN.fa
samtools faidx illumina_megahit_final_onlyHQBIN.fa
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn --wrap="bwa index illumina_megahit_final_onlyHQBIN.fa"

perl -lane 'if($F[0] eq "name"){next;}else{print $F[0];}' < pacbio_HQBIN_only_shortcov.tab > pacbio_HQBIN_only_shortcov.ctg.list
samtools faidx -r pacbio_HQBIN_only_shortcov.ctg.list ../../assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > pacbio_final_pilon_onlyHQBIN.fa
samtools faidx pacbio_final_pilon_onlyHQBIN.fa
sbatch --nodes=1 --ntasks-per-node=1 --mem=10000 -p msn --wrap="bwa index pacbio_final_pilon_onlyHQBIN.fa"

# OK, I've got our illumina coverage down, but I haven't aligned the pacbio EC reads using bwa mem yet. Let's queue that up
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem illumina_megahit_final_onlyHQBIN.fa ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta | samtools sort -o illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam -T illumina_pbtemp - ; samtools index illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam"
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem pacbio_final_pilon_onlyHQBIN.fa ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta | samtools sort -o pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam -T pacbio_pbtemp - ; samtools index pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam"

# Let's also queue up the hi-c read pairs
ls /beegfs/project/rumen_longread_metagenome_assembly/sequence_data/pilot_project/hic/*.fastq > hic_reads_spreadsheet.tab
# Editted the file to fit my pipeline's preferred format
vim hic_reads_spreadsheet.tab

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b illumina_hic -t hic_reads_spreadsheet.tab -f /beegfs/project/rumen_longread_metagenome_assembly/analysis/diff_cov/illumina_megahit_final_onlyHQBIN.fa -m -p msn
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b pacbio_hic -t hic_reads_spreadsheet.tab -f /beegfs/project/rumen_longread_metagenome_assembly/analysis/diff_cov/pacbio_final_pilon_onlyHQBIN.fa -m -p msn

sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p msn --wrap="jgi_summarize_bam_contig_depths --outputDepth illumina_hic_jgi_depth.tab illumina_hic/M1HIC/M1HIC.sorted.merged.bam illumina_hic/S3HIC/S3HIC.sorted.merged.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p msn --wrap="jgi_summarize_bam_contig_depths --outputDepth pacbio_hic_jgi_depth.tab pacbio_hic/M1HIC/M1HIC.sorted.merged.bam pacbio_hic/S3HIC/S3HIC.sorted.merged.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="jgi_summarize_bam_contig_depths --outputDepth pacbio_final_pilon_onlyHQBIN.pbecreads.depth.tab pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="jgi_summarize_bam_contig_depths --outputDepth illumina_megahit_final_onlyHQBIN.pbecreads.depth.tab illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam"


# Now to start consolidating everything together
perl -lane 'print "$F[0]\t$F[3]";' < pacbio_HQBIN_only_shortcov.tab > pacbio_HQBIN_only_shortcov.temp.tab
perl -lane 'print "$F[0]\t$F[3]\t$F[5]";' < pacbio_hic_jgi_depth.tab > pacbio_hic_jgi_depth.temp
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' pacbio_hic_jgi_depth.temp pacbio_HQBIN_only_shortcov.temp.tab > pacbio_HQBIN_shortread_avgcov.tab

perl -lane 'print "$F[0]\t$F[3]";' < illumina_HQBIN_only_shortcov.tab > illumina_HQBIN_only_shortcov.temp.tab
perl -lane 'print "$F[0]\t$F[3]\t$F[5]";' < illumina_hic_jgi_depth.tab > illumina_hic_jgi_depth.temp.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' illumina_hic_jgi_depth.temp.tab illumina_HQBIN_only_shortcov.temp.tab > illumina_HQBIN_shortread_avgcov.tab

# Now for the PB reads vs Illumina reads on both assemblies
perl -lane 'print "$F[0]\t$F[3]";' < illumina_megahit_final_onlyHQBIN.pbecreads.depth.tab > illumina_megahit_final_onlyHQBIN.pbecreads.depth.temp
perl -lane 'print "$F[0]\t$F[3]";' < pacbio_final_pilon_onlyHQBIN.pbecreads.depth.tab > pacbio_final_pilon_onlyHQBIN.pbecreads.depth.temp


perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' illumina_megahit_final_onlyHQBIN.pbecreads.depth.temp illumina_HQBIN_only_shortcov.temp.tab > illumina_HQBIN_sr_vs_lr_avgcov.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' pacbio_final_pilon_onlyHQBIN.pbecreads.depth.temp pacbio_HQBIN_only_shortcov.temp.tab > pacbio_HQBIN_sr_vs_lr_avgcov.tab
```

##### New ideas: 
1. Mapping efficiency plays a big role in Pacbio read alignment, so subsection the reads?
2. Also check to see if there's a phylum bias in the bins. This will be easier. Not sure if the pacbio alignments will change much, but we can do them to address any questions about alignment bias


```bash
module load bwa samtools metabat/2.12.1
# Generating new PacBio short read fastqs
sbatch --nodes=1 --ntasks-per-node=6 --mem=20000 -p msn --wrap="python3 ~/python_toolchain/sequenceData/sampleLongReadToShortRead.py -f /beegfs/project/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta -s rumen_pacbio_corrected.sr.R1.fq -r rumen_pacbio_corrected.sr.R2.fq -t 6"

sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem illumina_megahit_final_onlyHQBIN.fa rumen_pacbio_corrected.sr.R1.fq rumen_pacbio_corrected.sr.R2.fq | samtools sort -o illumina_HQBIN.pbecreads.leapfrog.sorted.bam -T illumina_pbtemp -; samtools index illumina_HQBIN.pbecreads.leapfrog.sorted.bam; jgi_summarize_bam_contig_depths --outputDepth illumina_pb_ecreads_leapfrog_depth.tab illumina_HQBIN.pbecreads.leapfrog.sorted.bam"
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem pacbio_final_pilon_onlyHQBIN.fa rumen_pacbio_corrected.sr.R1.fq rumen_pacbio_corrected.sr.R2.fq | samtools sort -o pacbio_HQBIN.pbecreads.leapfrog.sorted.bam -T pacbio_pbtemp -; samtools index pacbio_HQBIN.pbecreads.leapfrog.sorted.bam; jgi_summarize_bam_contig_depths --outputDepth pacbio_pb_ecreads_leapfrog_depth.tab pacbio_HQBIN.pbecreads.leapfrog.sorted.bam"

wc -l rumen_pacbio_corrected.sr.R1.fq
399818640 rumen_pacbio_corrected.sr.R1.fq
# That's approximately 99 million read pairs

# Now to prep them the same way as before and generate the plots for comparison
perl -lane 'if($F[0] eq "contigName"){next;} print "$F[0]\t$F[3]";' < pacbio_pb_ecreads_leapfrog_depth.tab > pacbio_pb_ecreads_leapfrog_depth.temp
perl -lane 'if($F[0] eq "contigName"){next;} print "$F[0]\t$F[3]";' < illumina_pb_ecreads_leapfrog_depth.tab > illumina_pb_ecreads_leapfrog_depth.temp

perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq "contigName"){print join("\t", @s) . "\n";next;} $long = $h{$s[0]}; print "$s[0]\t$s[1]\t$long\n";} close IN;' illumina_pb_ecreads_leapfrog_depth.temp illumina_HQBIN_sr_vs_lr_avgcov.tab > illumina_HQBIN_sr_vs_lr_avgcov.leapfrog.tab
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq "contigName"){print join("\t", @s) . "\n";next;} if(!exists($h{$s[0]})){print STDERR "Error!\n";} $long = $h{$s[0]}; print "$s[0]\t$s[1]\t$long\n";} close IN;' pacbio_pb_ecreads_leapfrog_depth.temp pacbio_HQBIN_sr_vs_lr_avgcov.tab > pacbio_HQBIN_sr_vs_lr_avgcov.leapfrog.tab
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Generating the tax consensus for the contigs in the HQ bins
perl generate_tax_table.pl pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.HQbins.taxtable.tab
perl generate_tax_table.pl illumina_megahit_master_table_4_1_2019.HQbins.short.tab > illumina_megahit_master_table_4_1_2019.HQbins.taxtable.tab

# I am just pulling the BlobTools phyla information (most stable taxonomic category in the bins it appears)
perl -ne 'chomp; @F = split(/\t/); print "$F[0]\t$F[5]\n";' < pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab > ../diff_cov/pacbio_HQbins_taxassignment.tab
perl -ne 'chomp; @F = split(/\t/); print "$F[0]\t$F[5]\n";' < illumina_megahit_master_table_4_1_2019.HQbins.short.tab > ../diff_cov/illumina_HQbins_taxassignment.tab
```


Now to plot the coverage differential plots in mmgenome2.

> F:/SharedFolders/metagenomics/pilot_manuscript/diff_cov/

```R
# PacbioHIC
library(mmgenome2)
covdataframe <- read.delim("pacbio_HQBIN_shortread_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe)

pdf(file="pacbio_short_read_M1HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_M1HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

pdf(file="pacbio_short_read_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_S3HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

pdf(file="pacbio_M1HIC_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_M1HIC", y = "cov_S3HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

#IlluminaHIC
covdataframe <- read.delim("illumina_HQBIN_shortread_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe)

pdf(file="illumina_short_read_M1HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_M1HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_short_read_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_S3HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_M1HIC_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_M1HIC", y = "cov_S3HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()


# Illumina sr vs lr
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe)

pdf(file="illumina_sr_vs_lr_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

# Generating coverage by len plots to see if there are biases
mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="illumina_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="illumina_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()

# I just want to make sure that the short read and long read data doesn't have an obvious bias. Let's get the coefficients from a linear model to make sure
ilsrlm <- lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)
summary(ilsrlm)

Call:
lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)

Residuals:
   Min     1Q Median     3Q    Max 
-33.22 -14.55  -7.16   4.38 807.70 

Coefficients:
                  Estimate  Std. Error t value Pr(>|t|)    
(Intercept)    45.26095831  1.85465448  24.404  < 2e-16 ***
mm.trim$length -0.00019230  0.00002624  -7.328 2.61e-13 ***
mm.trim$gc     -0.28991344  0.03881071  -7.470 9.02e-14 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 28.81 on 6855 degrees of freedom
Multiple R-squared:  0.01357,	Adjusted R-squared:  0.01328 
F-statistic: 47.14 on 2 and 6855 DF,  p-value: < 2.2e-16

illrlm <- lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)
summary(illrlm)

Call:
lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
 -4.246  -2.629  -1.719   1.120 202.217 

Coefficients:
                   Estimate   Std. Error t value Pr(>|t|)    
(Intercept)    -0.356198367  0.339818384  -1.048   0.2946    
mm.trim$length -0.000010549  0.000004808  -2.194   0.0283 *  
mm.trim$gc      0.067749436  0.007111078   9.527   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 5.278 on 6855 degrees of freedom
Multiple R-squared:  0.0151,	Adjusted R-squared:  0.01481 
F-statistic: 52.54 on 2 and 6855 DF,  p-value: < 2.2e-1


# Pacbio sr vs lr
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe)

pdf(file = "pacbio_sr_vs_lr_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="pacbio_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="pacbio_lr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()

pbsrlm <- lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)
summary(pbsrlm)

Call:
lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
-88.386 -18.075  -5.624  11.008 174.108 

Coefficients:
                   Estimate   Std. Error t value    Pr(>|t|)    
(Intercept)    109.37309258  11.39537216   9.598     < 2e-16 ***
mm.trim$length   0.00024681   0.00002309  10.691     < 2e-16 ***
mm.trim$gc      -1.22724830   0.24173592  -5.077 0.000000573 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 29.84 on 430 degrees of freedom
Multiple R-squared:  0.2597,	Adjusted R-squared:  0.2562 
F-statistic: 75.42 on 2 and 430 DF,  p-value: < 2.2e-16

pblrlm <- lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)
summary(pblrlm)

Call:
lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
-13.517  -2.942  -0.721   2.172  41.477 

Coefficients:
                   Estimate   Std. Error t value Pr(>|t|)    
(Intercept)    -0.093548111  1.781133609  -0.053  0.95814    
mm.trim$length  0.000033695  0.000003608   9.338  < 2e-16 ***
mm.trim$gc      0.110666364  0.037784107   2.929  0.00358 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.664 on 430 degrees of freedom
Multiple R-squared:  0.1757,	Adjusted R-squared:  0.1719 
F-statistic: 45.83 on 2 and 430 DF,  p-value: < 2.2e-16


# Let's try a plot with the contigs colored by phyla
# Pacbio by phyla
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
phyladatafram <- read.delim("pacbio_HQbins_taxassignment.tab", header=TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe, taxonomy = phyladatafram)

# One of the contigs was misassigned to the Arthropoda, let's change it to no-hit
mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda", "no-hit"))

pdf(file = "pacbio_sr_vs_lr_covplot.byphylum.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Illumina by phyla
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header =TRUE)
phyladatafram <- read.delim("illumina_HQbins_taxassignment.tab", header = TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe, taxonomy = phyladatafram)

mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda" | phylum == "Chordata" | phylum == "Mollusca" | phylum == "Nematoda", "no-hit"))
mm.trim <- mm.trim %>% mutate(phylum=replace(phylum, phylum == "Acidobacteria" | phylum == "Candidatus Desantisbacteria" | phylum == "Candidatus Melainabacteria" | phylum == "Candidatus Moranbacteria" | phylum == "Candidatus Nomurabacteria" | phylum == "Candidatus Saccharibacteria" | phylum == "Chlorobi" | phylum == "Chlorophyta" | phylum == "Elusimicrobia" | phylum == "Gemmatimonadetes" | phylum == "Kiritimatiellaeota" | phylum == "Thermotogae", "other"))

pdf(file="illumina_sr_vs_lr_covplot.byphylum.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Now trying the Leapfrog read alignments
# Illumina first
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.leapfrog.tab", header=TRUE)
phyladatafram <- read.delim("illumina_HQbins_taxassignment.tab", header = TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe, taxonomy = phyladatafram)

mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda" | phylum == "Chordata" | phylum == "Mollusca" | phylum == "Nematoda", "no-hit"))
mm.trim <- mm.trim %>% mutate(phylum=replace(phylum, phylum == "Acidobacteria" | phylum == "Candidatus Desantisbacteria" | phylum == "Candidatus Melainabacteria" | phylum == "Candidatus Moranbacteria" | phylum == "Candidatus Nomurabacteria" | phylum == "Candidatus Saccharibacteria" | phylum == "Chlorobi" | phylum == "Chlorophyta" | phylum == "Elusimicrobia" | phylum == "Gemmatimonadetes" | phylum == "Kiritimatiellaeota" | phylum == "Thermotogae", "other"))

pdf(file="illumina_sr_vs_lr_covplot.bygc.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_sr_vs_lr_covplot.byphylum.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Testing zero coverage contigs based on mapping bias
prevcovdata <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
count(prevcovdata[prevcovdata$long == 0,])
1  2836
count(covdataframe[covdataframe$long == 0,])
1   178  <- less mapping bias, but it's still there!

# Now Illumina coverage by length plots
mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="illumina_sr_vs_length_binplot.leap.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="illumina_lr_vs_length_binplot.leap.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()

# pacbio
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.leapfrog.tab", header = TRUE)
phyladatafram <- read.delim("pacbio_HQbins_taxassignment.tab", header=TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe, taxonomy = phyladatafram)

pdf(file = "pacbio_sr_vs_lr_covplot.bygc.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file = "pacbio_sr_vs_lr_covplot.byphylum.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

prevcovdata <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
cor.test(covdataframe$short, covdataframe$long)

	Pearson's product-moment correlation

data:  covdataframe$short and covdataframe$long
t = 20.767, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6567436 0.7513733
sample estimates:
      cor 
0.7072118

count(prevcovdata[prevcovdata$long == 0,])
1     8
count(covdataframe[covdataframe$long == 0,])
1     0


# Now for the pacbio leapfrog coverage by length plots
mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="pacbio_sr_vs_length_binplot.leap.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="pacbio_lr_vs_length_binplot.leap.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
```

#### Calculating the total count of reads mapping to HQ bins

I need to do this from the BAMs, otherwise I will not calculate it correctly!

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/

```bash
module load samtools

# TOTAL mapped reads
samtools idxstats publicdb/USDA/USDA.sorted.merged.bam | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[2] + $s[3];} print "$c\n";'
1148721358

# Mapped to contigs
samtools idxstats publicdb/USDA/USDA.sorted.merged.bam | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[2];} print "$c\n";'
535457601 = 46.6%

samtools idxstats publicdb/USDA/USDA.sorted.merged.bam | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -c 0 -l ~/rumen_longread_metagenome_assembly/analysis/master_tables/pacbio_final_pilon_master_table_4_1_2019.HQbins.ctg.list -d '\t' | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[2] + $s[3];} print "$c\n";'
14774478 = 1.3% of short reads

gunzip -c usda_pacbio_secpilon.USDA.sorted.merged.bam.cov.gz | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -c 0 -l ~/rumen_longread_metagenome_assembly/analysis/master_tables/pacbio_final_pilon_master_table_4_1_2019.MQbins.ctg.list -d '\t' | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
85551359 = 7.4% of short reads
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
# It turns out that the blobtools coverage datasets have this all calculated as well!
gunzip -c usda_illum_megahit.USDA.sorted.merged.bam.cov.gz | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -c 0 -l ~/rumen_longread_metagenome_assembly/analysis/master_tables/illumina_megahit_master_table_4_1_2019.HQbins.ctg.list -d '\t' | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
14270518 = 1.3% of short reads! Incredible!

gunzip -c usda_illum_megahit.USDA.sorted.merged.bam.cov.gz | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -c 0 -l ~/rumen_longread_metagenome_assembly/analysis/master_tables/illumina_megahit_master_table_4_1_2019.MQbins.ctg.list -d '\t' | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
111230244 = MQbins accounted for 9.7% of short reads
```