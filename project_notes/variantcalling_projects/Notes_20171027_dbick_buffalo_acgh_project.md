# Water Buffalo aCGH seg dup discovery
---
*10/27/2017*

These are my notes on converting the Buffalo aCGH results from 2014 over to the released UMD_CASPUR WB2.0 assembly and generating figures and tables for the manuscript.

## Coordinate conversions

I need to start first with the full coordinate conversion of UMD3 against UMD_CASPUR. Let's download that file and begin the conversion.

```bash


```
## Btau4 comparison

I need to compare marker coordinates against George's previous cattle survey to see how the data matches up.

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