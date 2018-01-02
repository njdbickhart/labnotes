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