# Water Buffalo aCGH seg dup discovery
---
*10/27/2017*

These are my notes on converting the Buffalo aCGH results from 2014 over to the released UMD_CASPUR WB2.0 assembly and generating figures and tables for the manuscript.

## Table of Contents
* [dACGH results](#dacgh)
* [ACGH ratio gathering](#acgh)

## Coordinate conversions

I need to start first with the full coordinate conversion of UMD3 against UMD_CASPUR. Let's download that file and begin the conversion.

```bash


```
## Btau4 comparison

I need to compare marker coordinates against George's previous cattle survey to see how the data matches up.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/buffalo_acgh

```bash
# I downloaded the CNVRs from George's 2010 paper and need to reformat them
perl -lane '$F[2] =~ s/,//g; $F[3] =~ s/,//g; $F[4] =~ s/,//g; if($F[7] =~ /^\d+/){$F[7] =~ s/\%//g; $F[7] /= 100;} print join("\t", @F);' < liu_acgh_2010_btau4.tab > liu_acgh_2010_btau4.fixed.tab

perl -lane 'print "$F[1]\t$F[2]\t$F[3]\t$F[5]\t$F[6]";' < liu_acgh_2010_btau4.fixed.tab > liu_acgh_2010_btau4.fixed.bed

# I converted the excel datasheets into tab delimited files first
for i in ITWB*.txt; do name=`echo $i | cut -f1`; dos2unix $i; done
for i in ITWB*.txt; do name=`echo $i | cut -f1`; echo $name; perl -lane 'unless($F[0] =~ /\d+/){next;} $score = ($F[6] == 0)? $F[7]: $F[6]; $class = ($F[6] > 0)? "gain" : "loss"; print "$F[1]\t$F[2]\t$F[3]\t$score\t$class";' < $i > $name.bed; done

# Reformatting the bed data for my annotation program
for i in ITWB*.bed; do name=`echo $i | cut -d'.' -f1`; echo $name; perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[1]\t$s[2]\t$s[4]\t$s[3]\n";}' < $i > $name.format.bed; done
```


## UMD3 comparison

This is to check to see if the UMD3 coordinates of my JaRMS analysis match up with the aCGH results. aCGH coordinates are based on UMD3.1. Using the format bed files generated in the above Btau4 comparison section.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/buffalo_acgh

```bash
cp /mnt/nfs/nfs2/dbickhart/buffalo/tables/buffalo_jarms_fixed_segdups.tab ./

cp /mnt/nfs/nfs2/dbickhart/buffalo/gene_data/*.bed ./gene_data/
ls /mnt/nfs/nfs2/bickhart-users/buffalo_acgh/gene_data/*.bed > gene_data/db_list
ls /mnt/nfs/nfs2/bickhart-users/buffalo_acgh/*.format.bed | perl -lane '@bsegs = split(/\//, $F[0]); @dsegs = split(/\./, $bsegs[-1]); print "$F[0]\t$dsegs[0]";' > animal_format_bed.list

# Creating CNVRs
java -jar /mnt/nfs/nfs2/bickhart-users/binaries/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d gene_data/db_list -i animal_format_bed.list -o ITWB_acgh_combined -t

perl -lane 'if($F[0] =~ /^Region/){next;}else{print "$F[1]\t$F[2]\t$F[3]\t$F[4]";}' < ITWB_acgh_combined_regions.tab > ITWB_acgh_combined_regions.bed

wc -l ITWB_acgh_combined_regions.bed
805 ITWB_acgh_combined_regions.bed

perl -lane 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]";' < buffalo_jarms_fixed_segdups.tab > buffalo_jarms_fixed_segdups.bed
bedtools intersect -a ITWB_acgh_combined_regions.bed -b buffalo_jarms_fixed_segdups.bed -wa -wb | wc -l
30 <- out of 130 jarms segdups (23% overlap)

cp /mnt/nfs/nfs2/dbickhart/buffalo/tables/buffalo_jarms_fixed_del_regions.bed ./
bedtools intersect -a ITWB_acgh_combined_regions.bed -b buffalo_jarms_fixed_del_regions.bed -wa -wb | wc -l  
684 <- out of 5478 jarms fixed dels (12.5% overlap)

bedtools intersect -a ITWB_acgh_combined_regions.bed -b buffalo_jarms_fixed_del_regions.bed -wa | uniq | wc -l
# 122 uniq del regions
# plus 30 uniq dup regions
# 152 / 805 = 18% overlap


# Now for the full CNVR list
bedtools sort -i /mnt/nfs/nfs2/dbickhart/buffalo/tables/buffalo_combined_cnvrs.bed > buffalo_combined_cnvrs.bed
cat buffalo_combined_cnvrs.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       4549
        Total Length:           44247951
        Length Average:         9726.9621894922
        Length Median:          5999
        Length Stdev:           16309.2004859253
        Smallest Length:        999
        Largest Length:         394999

# Non-uniq overlap
bedtools intersect -a ITWB_acgh_combined_regions.bed -b buffalo_combined_cnvrs.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1192
        Total Length:           16381468
        Length Average:         13742.8422818792
        Length Median:          7709
        Length Stdev:           21549.6011346469
        Smallest Length:        259
        Largest Length:         332059

# Uniq overlap
bedtools intersect -a ITWB_acgh_combined_regions.bed -b buffalo_combined_cnvrs.bed -wa | uniq | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       221
        Total Length:           178052739
        Length Average:         805668.502262443
        Length Median:          21759
        Length Stdev:           9565400.98294736
        Smallest Length:        1459
        Largest Length:         142490159
```
<a name="dacgh"></a>
## "dACGH" results

I know that George wants to mirror our 2012 GR manuscript, but many of the methods that we used were highly subjective and required extensive reference vetting. I will use the statistical models to try to replicate the same methods we used, but with some added veracity in prediction.

Unfortunately, I also need to realign all of the damn reads to work with the files!

> Assembler2: /mnt/nfs/nfs2/bickhart-users/buffalo_acgh

```bash
# Generating a spreadsheet for my alignment script
for i in /mnt/cifs/bickhart-qnap/buffalo_sequence/Sample_*; do name=`basename $i | cut -d'_' -f2`; echo $name; perl -e 'chomp @ARGV; @f = `ls $ARGV[0]`; %h; foreach $b (grep(/.+gz$/, @f)){chomp $b; @bsegs = split(/[\._]/, $b); $h{$bsegs[2]}->{$bsegs[4]}->{$bsegs[3]} = $b;} foreach my $l (keys(%h)){foreach my $n (keys(%{$h{$l}})){print "$ARGV[0]\/$h{$l}->{$n}->{R1}\t$ARGV[0]\/$h{$l}->{$n}->{R2}\t$ARGV[1]\t$ARGV[1]\n";}}' $i $name > $name.seqreads.tab; done

# Running the buffalo reads
for i in *.seqreads.tab; do name=`echo $i | cut -d'.' -f1`; echo $name; perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b $name -t $i -f umd3_kary_unmask_ngap.fa -m -p assemble1; done

ln -s /mnt/nfs/nfs2/dbickhart/dominette_asm/umd3/dominette.merged.umd3.bam
samtools index dominette.merged.umd3.bam

# I need to replicate the aCGH findings that George had in the past. Let's use similar settings.
# In the 2010 GR manuscript, he used 0.5 log2 ratio and 3 consecutive windows (0.5_3) as well as 0.3_5 and 0.3_3. 
# I will start with 0.5_3 here and a window size of 1kb (from our 2012 GR manuscript)
module load samtools; sbatch -p assemble1 ../binaries/cnv-seq/cnv-seq.pl --test ITWB1/ITWB1/ITWB1.sorted.merged.bam --ref dominette.merged.umd3.bam --window-size 1000 --genome-size 2800000000 --annotate

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7 ITWB9 PC1; do echo $i; sbatch -p assemble1 ../binaries/cnv-seq/cnv-seq.pl --test $i/$i/$i.sorted.merged.bam --ref dominette.merged.umd3.bam --window-size 1000 --genome-size 2800000000 --annotate; done

# Dammit all!! The reference genome for that umd3 copy was the NCBI version!
perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b dominette -t dominette_sequence_files.tab -f umd3_kary_unmask_ngap.fa -p assemble1 -m

# Now trying it on the dominette-umd3 comparison
module load samtools; sbatch -p assemble1 ../binaries/cnv-seq/cnv-seq.pl --test ITWB1/ITWB1/ITWB1.sorted.merged.bam --ref dominette/dominette/dominette.sorted.merged.bam --window-size 1000 --genome-size 2800000000 --annotate

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7 ITWB9 PC1; do echo $i; sbatch -p assemble1 ../binaries/cnv-seq/cnv-seq.pl --test $i/$i/$i.sorted.merged.bam --ref dominette/dominette/dominette.sorted.merged.bam --window-size 1000 --genome-size 2800000000 --annotate; done

# OK, now I need to load the package in R to call CNVs individually
```

```R
library(cnv)
files <- list.files(pattern = "*.count")
data <- cnv.cal("ITWB1.sorted.merged.window-1000.minw-4.count", log2.threshold=0.5, minimum.window=3, annotate=TRUE)

# It worked, but the annotation takes too long to segment CNV calls
for(f in files){data <- cnv.cal(f, log2.threshold=0.5, minimum.window=3, annotate=FALSE); write.table(data, file=paste0(f, ".tab"), sep="\t", quote=FALSE);}
```

Hold the phone -- George said that PC1 was the reference set for the aCGH comparison! Let me recalculate the ratios based on the PC1 reads.

```bash
module load samtools
for i in ITWB1 ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7 ITWB9; do echo $i; sbatch -p assemble1 ../binaries/cnv-seq/cnv-seq.pl --test $i/$i/$i.sorted.merged.bam --ref PC1/PC1/PC1.sorted.merged.bam --window-size 1000 --genome-size 2800000000; done

```

<a name="acgh"></a>
## Gathering aCGH ratios and generating comparisons

George uploaded the raw data on the cluster. Unfortunately the raw data was incredibly, incredibly discursive and complex. I was able to find the logratio in the files and I will test out some stats to see if it is suitable for comparison.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/buffalo_acgh
 
```bash
perl -e 'for($x = 0; $x < 10; $x++){<>;} while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); printf("%.10g\n", $s[15]);}' < aCGH_buffalo/US83800207_255144910001_S01_CGH_1105_Oct12.txt |  perl ~/sperl/bed_cnv_fig_table_pipeline/statStd.pl
total   972478
Minimum -1.234553856
Maximum 2.196189466
Average -0.001129
Median  0.0032069637105
Standard Deviation      0.070513
Mode(Highest Distributed Value) 0

# Hmm... I can't find the probe coordinates. I will try to infer them from the probe sequence 
perl -lane 'if($F[6] =~ /[ACGT]+/){print ">$F[1]\n$F[6]";}' < aCGH_buffalo/US83800207_255144910001_S01_CGH_1105_Oct12.txt > acgh_probes.fa
bwa mem /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa acgh_probes.fa > acgh_probes.umd3.sam

perl -lane 'if($F[0] =~ /^@/){next;}else{$e = length($F[9]) + $F[3]; print "$F[2]\t$F[3]\t$e\t$F[0]";}' < acgh_probes.umd3.sam > acgh_probes.umd3.bed

# OK, so now I need to grep out all of the logr values, then associate them with specific coordinates
for i in US83800207_255144910049_S01_CGH_1105_Oct12.txt US83800207_255144910050_S01_CGH_1105_Oct12.txt US83800207_255144910090_S01_CGH_1105_Oct12.txt US83800207_255144910091_S01_CGH_1105_Oct12.txt US83800207_255144910092_S01_CGH_1105_Oct12.txt US83800207_255144910093_S01_CGH_1105_Oct12.txt US83800207_255144910094_S01_CGH_1105_Oct12.txt US83800207_255144910095_S01_CGH_1105_Oct12.txt US83800207_255144910096_S01_CGH_1105_Oct12.txt US83800207_255144910097_S01_CGH_1105_Oct12.txt US83800207_255144910098_S01_CGH_1105_Oct12.txt US83800207_255144910099_S01_CGH_1105_Oct12.txt US83800207_255144910103_S01_CGH_1105_Oct12.txt US83800207_255144910104_S01_CGH_1105_Oct12.txt; do name=`echo $i | cut -d'_' -f1,2`; echo $name; perl -e 'for($x = 0; $x < 10; $x++){<>;} while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); printf("%s\t%.10g\n", $s[1], $s[15]);}' < aCGH_buffalo/$i > $name.logratios; done

# I changed them all to their ITWB callsigns
# Now to swap probe numbers for the coordinates
for i in *.logratios; do echo $i; perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[3]} = [$s[0], $s[1], $s[2]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){print join("\t", @{$h{$s[0]}}) . "\t$s[1]\n";}}' acgh_probes.umd3.bed $i > $i.bed; done

for i in *logratios.bed; do name=`echo $i | cut -d'.' -f1`; echo $name; bedtools sort -i $i > $name.logratios.sorted.bed; done

# modifying the data so that I can compare it
for i in *-vs-PC1*.count; do name=`echo $i | cut -d'.' -f1`; echo $name; mv $i $name.sorted.merged.window-1000.minw-4.count; done

# Ok, now I'm going to do a direct intersection and then generate the data comparison
for i in ITWB*.count.tab; do echo $i; perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\t$s[2]\t$s[3]\t$s[7]\t$s[8]\n";}' < $i > $i.bed; done
for i in *.count.tab.bed; do name=`echo $i | cut -d'.' -f1`; echo $name; bedtools sort -i $i > $name.dacgh.sorted.bed; done

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do bedtools intersect -a $i.dacgh.sorted.bed -b $i.logratios.sorted.bed -wa -wb > $i.dacgh.acgh.combined.tab; done

# Now to generate the input data for the correlation
for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); print "$s[3]\t$s[8]\t$ARGV[1]\n";} close IN;' $i.dacgh.acgh.combined.tab $i >> complete.dacgh.acgh.values.tab; done

# And the CNV intersected regions, only
# NOTE: I'm only using the aCGH CNV regions here
for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do bedtools intersect -a ITWB_acgh_combined_regions.bed -b $i.dacgh.acgh.combined.tab -wb | perl -e 'chomp(@ARGV); while(<STDIN>){chomp; @s = split(/\t/); print "$s[7]\t$s[12]\t$ARGV[0]\n";}' $i >> cnv.dacgh.acgh.values.tab; done

wc -l cnv.dacgh.acgh.values.tab complete.dacgh.acgh.values.tab
   2237430 cnv.dacgh.acgh.values.tab
  25113675 complete.dacgh.acgh.values.tab

# There was a problem with Infinite values
for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[3] =~ /Inf/){next;} print "$s[3]\t$s[8]\t$ARGV[1]\n";} close IN;' $i.dacgh.acgh.combined.tab $i >> complete.dacgh.acgh.values.tab; done

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do bedtools intersect -a ITWB_acgh_combined_regions.bed -b $i.dacgh.acgh.combined.tab -wb | perl -e 'chomp(@ARGV); while(<STDIN>){chomp; @s = split(/\t/); if($s[7] =~ /Inf/){next;} print "$s[7]\t$s[12]\t$ARGV[0]\n";}' $i >> cnv.dacgh.acgh.values.tab; done
```

Now to attempt to plot the correlations.

```R
library(ggplot2)
complete <- read.delim("complete.dacgh.acgh.values.tab", header=FALSE)
colnames(complete) <- c("dacgh", "acgh", "Sample")

cor.test(complete$dacgh, complete$acgh, use = "complete.obs")

        Pearson's product-moment correlation

data:  complete$dacgh and complete$acgh
t = 1307.6, df = 22825000, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.2636033 0.2643666
sample estimates:
     cor
0.263985

# So the complete dataset has poor correlation

cnv <- read.delim("cnv.dacgh.acgh.values.tab", header=FALSE)
colnames(cnv) <- c("dacgh", "acgh", "Sample")

cor.test(cnv$dacgh, cnv$acgh, use = "complete.obs")

        Pearson's product-moment correlation

data:  cnv$dacgh and cnv$acgh
t = 1204.1, df = 1822100, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6648626 0.6664798
sample estimates:
     cor
0.665672

# Much better, but not amazing
cnv$Sample <- as.factor(cnv$Sample)
png(file="correlation_cnvregions_dacgh.png", height=2000, width=2000)
ggplot(cnv, aes(x=dacgh, y=acgh, color=Sample)) + geom_point(shape=19, alpha=0.3) + geom_smooth(method=lm, se=FALSE) + theme_set(theme_gray(base_size = 24))


```

#### Original MrsFAST CN values

The original MrsFAST computed CN values are located here:

> Assembler2: /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/

Let's try to generate the associations again so that we can reestimate the corelations between the aCGH and the original MrsFAST data.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/buffalo_acgh 

```bash
# First, let's calculate the log2 ratio of the test / reference copy number
for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do echo $i; perl calculate_log2_ratio.pl /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/PC1_umd3/hits_umd3_template_file3.bed.gc.depth.normalized.CN /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/Sample_${i}/hits_umd3_template_file3.bed.gc.depth.normalized.CN ${i}_mrsfast_log2.bed; done


module load bedtools/2.26.0; for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do bedtools intersect -a ${i}_mrsfast_log2.bed -b $i.logratios.sorted.bed -wa -wb > $i.mrsfast.dacgh.acgh.combined.bed; done

for i in *mrsfast.dacgh.acgh.combined.bed; do name=`echo $i | cut -d'.' -f1`; echo $name; perl -lane 'print "$F[3]\t$F[7]";' < $i > $name.mrsfast.dacgh.acgh.values.tab; done

# OK, I forgot to select only the CNV region areas! Let's do that quickly
cat /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/*/*.final* | bedtools sort -i stdin | bedtools merge -i stdin > mrsfast_buffalo_cnv_regions.bed

cat mrsfast_buffalo_cnv_regions.bed | perl ~/sperl/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       1371
        Total Length:           60,062,616
        Length Average:         43809.3479212254
        Length Median:          22815
        Length Stdev:           68334.5305998886
        Smallest Length:        989
        Largest Length:         885429

# Looks pretty similar to my previous counts. Let's intersect it with the cn values
mkdir mrsfast_cnv_cns
bedtools intersect -a mrsfast_buffalo_cnv_regions.bed -b /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/PC1_umd3/hits_umd3_template_file3.bed.gc.depth.normalized.CN -wb | perl -lane 'print "$F[3]\t$F[4]\t$F[5]\t$F[6]";' > mrsfast_cnv_cns/PC1_umd3_cnv_cns.bed

wc -l mrsfast_cnv_cns/PC1_umd3_cnv_cns.bed
24188 mrsfast_cnv_cns/PC1_umd3_cnv_cns.bed # Not as many windows! about 100 fold less

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do echo $i; bedtools intersect -a mrsfast_buffalo_cnv_regions.bed -b /mnt/nfs/nfs2/bickhart-users/dbickhart_blade2/Buffalo/umd3_calls/Sample_${i}/hits_umd3_template_file3.bed.gc.depth.normalized.CN -wb | perl -lane 'print "$F[3]\t$F[4]\t$F[5]\t$F[6]";' > mrsfast_cnv_cns/${i}_umd3_cnv_cns.bed; done

for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do echo $i; perl calculate_log2_ratio.pl mrsfast_cnv_cns/PC1_umd3_cnv_cns.bed mrsfast_cnv_cns/${i}_umd3_cnv_cns.bed mrsfast_cnv_cns/${i}_mrsfast_cnvcns_log2.bed; done

# Now to combine the acgh with the dacgh calls
mkdir mrsfast_cnv_acgh/
for i in ITWB10 ITWB11 ITWB12 ITWB13 ITWB14 ITWB15 ITWB2 ITWB1 ITWB3 ITWB4 ITWB5 ITWB6 ITWB7; do bedtools intersect -a mrsfast_cnv_cns/${i}._mrsfast_cnvcns_log2.bed -b $i.logratios.sorted.bed -wa -wb > mrsfast_cnv_acgh/$i.mrsfast.dacgh.acgh.combined.bed; done

# Stripping it down to just the values:
for i in mrsfast_cnv_acgh/*.bed; do name=`basename $i | cut -d'.' -f1`; echo $name; perl -lane 'print "$F[3]\t$F[7]";' < $i > mrsfast_cnv_acgh/$name.mrsfast.dacgh.acgh.cnvcn.values.tab; done
```

OK, now to plot them individually as George wants to see individual correlations and plots

```R
buff <- c("ITWB10", "ITWB11", "ITWB12", "ITWB13", "ITWB14", "ITWB15", "ITWB1", "ITWB2", "ITWB3", "ITWB4", "ITWB5", "ITWB6", "ITWB7")

library(ggplot2); for(b in buff){print(b); outpng <- paste0(b, ".mrsfast.corr.png"); infile <- paste0(b, ".mrsfast.dacgh.acgh.values.tab"); cnv <- read.delim(infile, header=FALSE); colnames(cnv) <- c("MrsFAST", "aCGH"); cor <- cor.test(cnv$MrsFAST, cnv$aCGH, use = "complete.obs"); print(cor$estimate);  png(file=outpng, height=2000, width=2000); print(ggplot(cnv, aes(x=MrsFAST, y=aCGH)) + geom_point(shape=19, alpha=0.5) + geom_smooth(method=lm, se=FALSE) + theme_set(theme_gray(base_size=24)) + ggtitle(paste0("dACGH correlation for ", b, " (r = ", cor$estimate, " )"))); dev.off();}
```