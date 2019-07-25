# Yaklander genome assembly analysis
---
*3/25/2019*

These are my notes on running some analysis on the Yaklander (Yak x Highland cross) triobinned assembly.

## Table of Contents
* [Preparing the assembly fastas](#prep)
* [Running the repeat analysis](#repeat)
* [Plotting the data and generating summary stats](#plot)

<a name="prep"></a>
## Preparing the assembly fastas

I just need to gather the materials I need for the analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander

```bash
module load bwa
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.sire.gapfilled.arrow2.fasta"
```

<a name="repeat"></a>
## Running the repeat analysis

OK, now that everything is ready, let's queue up RepeatMasker and generate the files I need. I've already done repeatmasking on the ARS-UCDv1.2 assembly and the UMD3.1 assembly, so those comparative repeat lengths can be added to the plots if needed for comparison.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander
 
```bash
module unload perl/5.24.1
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.sire.gapfilled.arrow2.fasta"

# I've combined the gap and repeat analysis into one rolling script
sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.arrow2.fasta yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch generate_repeat_counts.sh yaklander.sire.gapfilled.arrow2.fasta yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa

# Grrr! The "|arrow" delimiters strike again! I need to remove them 
module load bwa samtools
perl -e 'while(<>){if($_ =~ /^>/){$_ =~ s/\|arrow//g;} print $_;}' < yaklander.dam.gapfilled.arrow2.fasta > yaklander.dam.gapfilled.rfmt.fa
perl -e 'while(<>){if($_ =~ /^>/){$_ =~ s/\|arrow//g;} print $_;}' < yaklander.sire.gapfilled.arrow2.fasta > yaklander.sire.gapfilled.rfmt.fa
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index yaklander.dam.gapfilled.rfmt.fa; samtools faidx yaklander.dam.gapfilled.rfmt.fa"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index yaklander.sire.gapfilled.rfmt.fa; samtools faidx yaklander.sire.gapfilled.rfmt.fa"

sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.arrow2.fasta yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch generate_repeat_counts.sh yaklander.sire.gapfilled.arrow2.fasta yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa

sbatch generate_repeat_counts.sh yaklander.dam.gapfilled.rfmt.fa yakdam ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa
sbatch --dependency=afterok:659870 generate_repeat_counts.sh yaklander.sire.gapfilled.rfmt.fa yaksire ../dominette/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa


# Ack! I screwed up my script and read in the fasta file instead of the repeatmasker output!
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[4] =~ s/\|arrow//g; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < yaklander.dam.gapfilled.arrow2.fasta.out > yakdam.repeat.bed
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f yakdam.repeat.bed -c 5 -d '\t' -m > yakdam.rep.count.md
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[5] =~ s/\?//g; $F[5] =~ tr/\//_/; print "$F[5]\t$F[6]\t$ARGV[1]\n";} close IN;' yakdam.repeat.bed yakdam > yakdam.repeat.lens

perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[4] =~ s/\|arrow//g; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < yaklander.sire.gapfilled.arrow2.fasta.out > yaksire.repeat.bed
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f yaksire.repeat.bed -c 5 -d '\t' -m > yaksire.rep.count.md
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[5] =~ s/\?//g; $F[5] =~ tr/\//_/; print "$F[5]\t$F[6]\t$ARGV[1]\n";} close IN;' yaksire.repeat.bed yaksire > yaksire.repeat.lens

# Getting the other data together...
cp ../dominette/repeatmasker/umd3_reference_genome.fasta.out.rep.lens ./
cp ../dominette/repeatmasker/ARS-UCD1.2_Btau5.0.1Y.fa.out.rep.lens ./

cat *.rep*.lens > combined_yaklander_repeat.lens
```

<a name="plot"></a>
## Plotting the data and generating summary stats

Now to plot the most divergent repeat class groups.

```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(scales)

data <- read.delim("combined_yaklander_repeat.lens", header=FALSE)
colnames(data) <- c("Class", "Len", "ASM")

nrow(data[data$Len == 0,])
[1] 291   <- This shouldn't be possible...
summary(data[data$Len == 0,])

              Class          Len         ASM
 LINE_L1         :169   Min.   :0   ARSUCD :73
 LTR_ERV1        : 29   1st Qu.:0   UMD3   :72
 LTR_ERVL-MaLR   : 28   Median :0   yakdam :68
 DNA_hAT-Charlie : 13   Mean   :0   yaksire:78
 LTR_ERVL        : 12   3rd Qu.:0
 DNA_TcMar-Tigger: 10   Max.   :0
 (Other)         : 30

# Looks like a handful of parsing errors on my part? Or perhaps an issue with RepeatMasker? 
# Still, they make up a very small proportion of repeats.
data.filt <- data[data$Len > 0,]

# Dam ttest for significant average repeat lenght differences with UMD3
ttest_yakdam <- group_by(data.filt, ASM, Class) %>% summarize(Len = list(Len)) %>% spread(ASM, Len) %>% group_by(Class) %>% mutate(p_value = t.test(unlist(yakdam), unlist(UMD3))$p.value, t_value = t.test(unlist(yakdam), unlist(UMD3))$statistic)

# Yak Dam significant size enrichments compared to UMD3
ttest_yakdam[ttest_yakdam$p_value < 0.05,]
# A tibble: 6 x 7
# Groups:   Class [6]
  Class        ARSUCD       UMD3        yakdam      yaksire      p_value t_value
  <fct>        <list>       <list>      <list>      <list>         <dbl>   <dbl>
1 Low_complex… <int [85,90… <int [82,7… <int [84,0… <int [83,9… 6.34e- 4    3.42
2 LTR_ERVK     <int [122,8… <int [113,… <int [113,… <int [119,… 3.58e- 7    5.09
3 Satellite_c… <int [11,98… <int [6,57… <int [8,06… <int [8,12… 2.56e-78   18.9
4 Simple_repe… <int [554,8… <int [536,… <int [535,… <int [546,… 2.87e-20    9.22
5 SINE_Core-R… <int [385,4… <int [378,… <int [379,… <int [379,… 4.93e- 7    5.03
6 SINE_tRNA-C… <int [1,038… <int [1,02… <int [1,02… <int [1,02… 7.97e- 6    4.47

# Sire ttest 
ttest_yaksire <- group_by(data.filt, ASM, Class) %>% summarize(Len = list(Len)) %>% spread(ASM, Len) %>% group_by(Class) %>% mutate(p_value = t.test(unlist(yaksire), unlist(UMD3))$p.value, t_value = t.test(unlist(yaksire), unlist(UMD3))$statistic)

# Yak Sire significant size enrichments vs UMD3
ttest_yaksire[ttest_yaksire$p_value < 0.05,]
# A tibble: 9 x 7
# Groups:   Class [9]
  Class        ARSUCD       UMD3        yakdam      yaksire      p_value t_value
  <fct>        <list>       <list>      <list>      <list>         <dbl>   <dbl>
1 LINE_L1      <int [871,2… <int [855,… <int [851,… <int [851,… 4.76e- 2   -1.98
2 LINE_RTE-Bo… <int [725,6… <int [711,… <int [707,… <int [707,… 2.69e- 2    2.21
3 Low_complex… <int [85,90… <int [82,7… <int [84,0… <int [83,9… 4.77e- 3    2.82
4 LTR_ERVK     <int [122,8… <int [113,… <int [113,… <int [119,… 1.78e-36   12.6
5 rRNA         <int [1,224… <int [1,31… <int [1,24… <int [1,48… 4.38e- 2   -2.02
6 Satellite_c… <int [11,98… <int [6,57… <int [8,06… <int [8,12… 3.53e-15    7.88
7 Simple_repe… <int [554,8… <int [536,… <int [535,… <int [546,… 3.51e-16    8.15
8 SINE_Core-R… <int [385,4… <int [378,… <int [379,… <int [379,… 1.37e- 5    4.35
9 SINE_tRNA-C… <int [1,038… <int [1,02… <int [1,02… <int [1,02… 8.18e- 7    4.93

# So the same class of characters, but the LINE and BOVB are new as are the rRNA
# Let's check against ARS-UCD -- should be fewer in both cases
# Yak dam
ttest_ayakdam <- group_by(data.filt, ASM, Class) %>% summarize(Len = list(Len)) %>% spread(ASM, Len) %>% group_by(Class) %>% mutate(p_value = t.test(unlist(yakdam), unlist(ARSUCD))$p.value, t_value = t.test(unlist(yakdam), unlist(ARSUCD))$statistic)

ttest_ayakdam[ttest_ayakdam$p_value < 0.05,]
# A tibble: 4 x 7
# Groups:   Class [4]
  Class       ARSUCD       UMD3       yakdam      yaksire        p_value t_value
  <fct>       <list>       <list>     <list>      <list>           <dbl>   <dbl>
1 LTR_ERVK    <int [122,8… <int [113… <int [113,… <int [119,…    2.55e-9   -5.96
2 LTR_ERVL-M… <int [145,5… <int [143… <int [143,… <int [143,…    5.55e-3   -2.77
3 Satellite_… <int [11,98… <int [6,5… <int [8,06… <int [8,12…    1.52e-3    3.17
4 Simple_rep… <int [554,8… <int [536… <int [535,… <int [546,…    2.51e-3    3.02

# Aha! Larger centromeres you say? 

# Yak Sire
ttest_ayaksire <- group_by(data.filt, ASM, Class) %>% summarize(Len = list(Len)) %>% spread(ASM, Len) %>% group_by(Class) %>% mutate(p_value = t.test(unlist(yaksire), unlist(ARSUCD))$p.value, t_value = t.test(unlist(yaksire), unlist(ARSUCD))$statistic)

ttest_ayaksire[ttest_ayaksire$p_value < 0.05,]
# A tibble: 5 x 7
# Groups:   Class [5]
  Class        ARSUCD       UMD3        yakdam      yaksire      p_value t_value
  <fct>        <list>       <list>      <list>      <list>         <dbl>   <dbl>
1 LINE_L1      <int [871,2… <int [855,… <int [851,… <int [851,… 3.51e- 2   -2.11
2 LTR_ERVL-Ma… <int [145,5… <int [143,… <int [143,… <int [143,… 4.03e- 2   -2.05
3 Satellite_c… <int [11,98… <int [6,57… <int [8,06… <int [8,12… 3.78e-31  -11.6
4 Simple_repe… <int [554,8… <int [536,… <int [535,… <int [546,… 4.33e- 3    2.85
5 SINE_tRNA-C… <int [1,038… <int [1,02… <int [1,02… <int [1,02… 3.32e- 2    2.13

# Hmm... and the sire's centromeres are smaller? Let's plot the centromeres first
data.cent <- data.filt[data.filt$Class == "Satellite_centr",]
summary(data.cent)
               Class            Len              ASM
 Satellite_centr  :34747   Min.   :     9   ARSUCD :11988
 DNA              :    0   1st Qu.:   432   UMD3   : 6574
 DNA_Harbinger    :    0   Median :   651   yakdam : 8061
 DNA_hAT          :    0   Mean   :  3322   yaksire: 8124
 DNA_hAT-Ac       :    0   3rd Qu.:  1409
 DNA_hAT-Blackjack:    0   Max.   :523277

# More centromeric repeats in ARSUCD
pdf(file="centromere_lengths_ridges.pdf", useDingbats=FALSE)
ggplot(data.cent, aes(x=Len, y=ASM, fill=ASM)) + geom_density_ridges(scale=3, rel_min_height = 0.01) + scale_fill_brewer(palette="Dark2") + labs(title = "Centromere/Satellite Region Lengths", xlab="Length (Log10 bp)", ylab="Assembly") + scale_x_log10(labels=comma)
dev.off()

# Now for a boxplot of the "usual suspect" ttest groups
# The Yak Sire dataset had the most groups. Let's plot that.
toplot <- ttest_yaksire[ttest_yaksire$p_value < 0.05,"Class"] %>% pull(Class)
toplot <- as.factor(as.character(toplot))

data.sig <- data.filt[data.filt$Class %in% toplot,]
data.sig$Class <- as.factor(as.character(data.sig$Class))

rtable <- data.sig %>% group_by(ASM, Class) %>% summarize(mean = mean(Len), sd = sd(Len))

pdf(file="sig_repeat_category_averages.pdf", useDingbats=FALSE)
ggplot(rtable, aes(x=Class, y=mean, fill=ASM)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position=position_dodge(0.9)) + labs(title = "Significantly different average Repeat Class Lengths")
dev.off()

# Since the sd was far larger than the mean in most cases, this made the plot very unwieldy!
# Let's just plot the means
pdf(file="sig_repeat_cat_avg_woerror.pdf", useDingbats=FALSE)
ggplot(rtable, aes(x=Class, y=mean, fill=ASM)) + geom_bar(stat="identity", color="black", position=position_dodge()) + labs(title="Significantly different average Repeat Class Lengths") + scale_fill_brewer(palette="Dark2") + theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

# Some of the repeats have different distributions that are not adequately matched by the means. Let's plot them separately
mixed <- as.factor(c("LINE_L1", "LINE_RTE-BovB", "Low_complexity", "Simple_repeat"))

data.mixed <- data.sig[data.sig$Class %in% mixed,]
data.mixed$Class <- as.factor(as.character(data.mixed$Class))

pdf(file="repeat_mixed_dist_plot.pdf", useDingbats=FALSE)
ggplot(data.mixed, aes(y=Len, x=Class, fill=ASM)) + geom_violin(position=position_dodge(1)) + scale_fill_brewer(palette = "Dark2") + labs(title = "Repeat classes with different length distributions") + scale_y_log10(labels=comma) + theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

# Finally, let's tabulate most of the results so that Ed can go through them himself
ttest_table <- group_by(data.filt, ASM, Class) %>% summarize(Len = list(Len)) %>% spread(ASM, Len) %>% group_by(Class) %>% mutate(ARSNum = length(unlist(ARSUCD)), ARSMean = mean(unlist(ARSUCD)), ARSSd = sd(unlist(ARSUCD)), UMD3Num = length(unlist(UMD3)), UMD3Mean = mean(unlist(UMD3)), UMD3Sd = sd(unlist(UMD3)), DamNum = length(unlist(yakdam)), DamMean = mean(unlist(yakdam)), DamSd = sd(unlist(yakdam)), SireNum = length(unlist(yaksire)), SireMean = mean(unlist(yaksire)), SireSd = sd(unlist(yaksire)), SireVsUMDPvalue = t.test(unlist(yaksire), unlist(UMD3))$p.value, SireVsARSPvalue = t.test(unlist(yaksire), unlist(ARSUCD))$p.value, DamVsUMDPvalue = t.test(unlist(yakdam), unlist(UMD3))$p.value, DamVsARSPvalue = t.test(unlist(yakdam), unlist(ARSUCD))$p.value)

write.table(ttest_table[,!(names(ttest_table) %in% c("ARSUCD", "UMD3", "yakdam", "yaksire"))], file="RepeatClassSummaries.tab", row.names=FALSE, quote=FALSE, sep="\t")
```

## Repeat gap intersection analysis

I want to see just how many gaps intersect with major repeat classes.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander

```bash
module load bedtools
python3 ~/python_toolchain/sequenceData/intersectGapFlanksWithRepeats.py -g yaksire.gapstatus.tab -r yaksire.repeat.bed -o yaksire.gaprepeat

python3 ~/python_toolchain/sequenceData/intersectGapFlanksWithRepeats.py -g yakdam.gapstatus.tab -r yakdam.repeat.bed -o yakdam.gaprepeat
```

#### Sire Stats

|Class   |RepeatConsistent |RepeatComplex  |None  |
|:-------|----------------:|--------------:|-----:|
|Closed  |35               |175            |9     |
|Trans   |25               |142            |12    |
|Total   |60               |317            |21    |

#### Dam Stats

|Class   |RepeatConsistent |RepeatComplex  |None  |
|:-------|----------------:|--------------:|-----:|
|Closed  |29               |172            |12    |
|Trans   |25               |152            |8     |
|Total   |54               |324            |20    |

I've noticed that some of the gap regions are different in terms of repeat classification between the assemblies. Let's test this.

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[2], $s[3], $s[4]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $other = $data{$s[0]}; $cons = ($other->[0] eq $s[2] && $other->[1] eq $s[3] && $other->[2] eq $s[4])? 1 : 0; print join("\t", @s) . "\t$cons\t" . join("\t", @{$other}) . "\n";} close IN;' yakdam.gaprepeat.comp yaksire.gaprepeat.comp | perl -lane 'print $F[5];' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   398
Sum:    212
Minimum 0
Maximum 1
Average 0.532663
Median  1
Standard Deviation      0.499560
Mode(Highest Distributed Value) 1

# Just above 50%? Wonder if it's due to gap region though
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], $s[2], $s[3], $s[4]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $other = $data{$s[0]}; $cons = ($other->[1] eq $s[2] && $other->[2] eq $s[3] && $other->[3] eq $s[4])? 1 : 0; print join("\t", @s) . "\t$cons\t" . join("\t", @{$other}) . "\n";} close IN;' yakdam.gaprepeat.comp yaksire.gaprepeat.comp > yaklander.gaprepeat.comp

perl -lane 'if($F[1] eq $F[6]){print "1\n";}else{print "0\n";}' < yaklander.gaprepeat.comp | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   398
Sum:    337
Minimum 0
Maximum 1
Average 0.846734
Median  1
Standard Deviation      0.360697
Mode(Highest Distributed Value) 1  # 80% overlap in gap closure status

perl -lane 'if($F[1] eq $F[6]){}else{@j = split(/-/, $F[0]); print $j[0]}' < yaklander.gaprepeat.comp | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m
```

#### gap discrepancy counts
|Entry          | Value|
|:--------------|-----:|
|X              |    13|
|15             |     4|
|13             |     3|
|27             |     3|
|Unplaced       |    19|
|18             |     2|
|23             |     2|
|5              |     2|
|Y              |     2|
|10             |     1|
|12             |     1|
|14             |     1|
|16             |     1|
|1              |     1|
|28             |     1|
|29             |     1|
|2              |     1|
|3              |     1|
|8              |     1|
|9              |     1|

```bash
# After removing the discrepancy regions...
perl -lane 'if($F[1] eq $F[6]){print $F[5];}' < yaklander.gaprepeat.comp | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   337
Sum:    211
Minimum 0
Maximum 1
Average 0.626113					<-- only improved by 10% 
Median  1
Standard Deviation      0.484554
Mode(Highest Distributed Value) 1

# Let's ID how many unique repeat classes are causing the gaps
perl -lane 'print "$F[2]-$F[3]-$F[4]\n$F[7]-$F[8]-$F[9]";' < yaklander.gaprepeat.comp | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m

perl -lane 'foreach my $c (2,3,4,7,8,9){print "$F[$c]";}' < yaklander.gaprepeat.comp | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m
```

#### Highlight combinations

|Entry                                                    | Value|
|:--------------------------------------------------------|-----:|
|None-None-None                                           |    41|
|LINE/L1-LINE/L1-LINE/L1                                  |    25|
|LINE/L1-None-LINE/L1                                     |    24|
|LINE/RTE-BovB-LINE/RTE-BovB-LINE/RTE-BovB                |    15|
|LINE/RTE-BovB-None-LINE/RTE-BovB                         |    11|
|None-Simple_repeat-None                                  |    11|
|Satellite/centr-Satellite/centr-Satellite/centr          |    10|
|Satellite/centr-None-Satellite/centr                     |    10|

#### Repeat flanking counts

|Entry              | Value|
|:------------------|-----:|
|None               |   801|
|LINE/L1            |   391|
|LINE/RTE-BovB      |   288|
|SINE/tRNA-Core-RTE |   167|
|Simple_repeat      |   118|
|SINE/MIR           |    77|
|SINE/tRNA          |    73|
|SINE/Core-RTE      |    73|
|LTR/ERV1           |    72|
|Satellite/centr    |    65|
|LINE/L2            |    53|
|LTR/ERVL-MaLR      |    41|
|LTR/ERVL           |    40|
|DNA/hAT-Charlie    |    30|
|LTR/ERVK           |    26|
|Low_complexity     |    21|
|LINE/CR1           |    10|
|DNA/hAT-Tip100     |    10|


Now to test this out in R and find other correlations.

```R
yakgap <- read.delim("yaklander.gaprepeat.comp", header=FALSE)

table(yakgap$V2, yakgap$V6)

           0   1
  Closed 109 110
  Trans   77 102

# So, it looks like it's not really biased towards any one phenomenon. Though I do recall that the closed and trans status differs between assemblies.
```