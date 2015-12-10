# Ruminant variant genome-wide comparison project
---
*11/2/2015*

These are my notes on the pan ruminant comparative genomics project I'm hoping to present at PAG. 

Let's start by first generating the BAMs I need. Buffalo will take some time so let's get that started.

> Blade2: /mnt/iscsi/md3200i_4/dbickhart/Buffalo/new_umd3

```bash
# Generating the spreadsheet
ls ../../../schroeder/Project_WaterBuffalo/full/ | perl -e '%h; while(<>){chomp; @s = split(/_/); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > new_umd3_spreadsheet.tab
```

Damn, I can't install my version of Perl on the system! A prior local-lib installation is interfering with everything and I don't have time to correct it. I'm going to have to make due with the cow4 alignments I guess!

OH! Good news! Thanks to my checking my old notes, Steve had apparently transferred the old Buffalo data to my QNAP drive! Nice! Let's get the spreadsheets ready on Blade14 instead (where I know my pipeline can run).

> Blade14: /home/dbickhart/vnx_gliu_7/ruminant_project

```bash
# Getting the buffalo spreadsheet ready
ls /mnt/cifs/bickhart-qnap/buffalo_sequence/Sample_*/*.gz | perl -e 'use File::Basename; %h; while(<>){chomp; $b = basename($_); @s = split(/_/, $b); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > buffalo_fq_spreadsheet.tab

# Now the goat
ls /mnt/nfs/nfs2/GoatData/African-Goat-Resequencing/NextSeq500/150820_NS500432_0014_AH15AFBGXX/H15AFBGXX/African-Goat-ReSeq/G0*/*.gz | perl -e 'use File::Basename; %h; while(<>){chomp; $b = basename($_); @s = split(/_/, $b); $h{$s[0]}->{$s[4]}->{$s[3]} = $_;} foreach $an (keys(%h)){foreach $num (keys(%{$h{$an}})){ $fq1 = $h{$an}->{$num}->{"R1"}; $fq2 = $h{$an}->{$num}->{"R2"}; print "$fq1\t$fq2\t$an\t$an\n";}}' > goat_spreadsheet.tab

# Just some minor housecleaning -- I don't want to waste time with EFF and SNP calling!
cat buffalo_fq_spreadsheet.tab goat_spreadsheet.tab > combined_spreadsheet.tab
cp ~/.mergedpipeline.cnfg ./temp_pipeline.cnfg


perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs combined_spreadsheet.tab --output goat_buff_bams --reference ../reference/umd3_kary_unmask_ngap.fa --config temp_pipeline.cnfg --threads 20
```

OK, I had to remove the error message profiling (open3 wasn't working with BWA for some reason!) but I'm now aligning samples to UMD3. As far as the 100 bulls project, I have samples that are aligned by BWA (version 6.X, not 7.X) to the right reference genome, but with the wrong algorithm (SAMPE, not MEM). Let's lay out the overall plan:

#### Project/Manuscript workflow
* Identify MEI, misassembly and large scale genome divergence in Goat, Cattle and Buffalo
* Annotate genomic regions and identify gene families under expansion/contraction
* Focus on the transcriptional regulation regions with MEI discovery
* Additionally SNPs and INDELs in the regulatory regions
* Identify cis-regulatory elements that are conserved/altered between cattle/goat/buffalo
* Also focus on gene deletions
* Identify individual diversity profiles and remove them as potential CNV "noise"
* TARGET: PNAS if I can get enough data

#### 100 bulls data
* For the abstract: Run JaRMS on samples and identify misassembled regions
* Later: Run RAPTR-SV and pull trans-chr reads to identify MEI
* Use this data as the "background"

#### Buffalo and Goat
* Align to UMD3 and remove misassembled regions
	* Misassembled regions are defined as:
	* Regions that are deletions, but are not gaps and not repeats
* Run RAPTR-SV and identify trans-chr MEIs
* Estimate count of unmapped reads per genome (supplementary data)

## JaRMS test on 100 bulls data

OK I know that the BAM files for the 100 bulls samples follow this format within the 100_base_run folder:

> (sample_name)/(sample_name).rg.rd.full.sorted.merged.bam

I'm currently using 75% of the memory of the blade, so let's try looping through these files one at a time to get the preliminary data first.

> Blade14: /mnt/iscsi/vnx_gliu_7/100_base_run

```bash
# OK! Running! 
~/jdk1.8.0_05/bin/java -jar ~/JaRMS/store/JaRMS.jar call -i BIBR02/BIBR02.rg.rd.full.sorted.merged.bam -f ../reference/umd3_kary_unmask_ngap.fa -o BIBR02.jarms.calls -t 2

# NOTE: my equation for estimating the X coverage of the genome always underestimates it by a factor of 1 x 10^10. 
# I need to check what's wrong, adjust it, and add an option to manually set window size

# OK, I uploaded a new version of JaRMs, but I might just keep the 500bp windows (or lower them to 400?)
~/jdk1.8.0_05/bin/java -jar ~/JaRMS/store/JaRMS.jar call -i BIBR02/BIBR02.rg.rd.full.sorted.merged.bam -f ../reference/umd3_kary_unmask_ngap.fa -o BIBR02.jarms.400 -t 2 -w 400


for i in *cnvs.bed; do cat $i | ~/bin/bed_length_sum.pl; done
400bp wins
        Interval Numbers:       21171
        Total Length:           286758029
        Length Average:         13544.8504558122
        Length Median:          1199
        Min Length:             16399
        Max Length:             3999
        Length Stdev:           34721.7105164687
500bp wins
        Interval Numbers:       16882
        Total Length:           241222118
        Length Average:         14288.7168581922
        Length Median:          6749
        Min Length:             54499
        Max Length:             3999
        Length Stdev:           33876.4349035586

# That was the total, let's try just the duplications (I suspect the 400bp wins picked up more gaps)
for i in BIBR02.*cnvs.bed; do echo $i; grep 'duplication' $i | perl ~/bin/bed_length_sum.pl; done
BIBR02.jarms.400cnvs.bed
        Interval Numbers:       8848
        Total Length:           120026752
        Length Average:         13565.4104882459
        Length Median:          125599
        Min Length:             16399
        Max Length:             1999
        Length Stdev:           39448.5304808339
BIBR02.jarms.callscnvs.bed
        Interval Numbers:       6494
        Total Length:           87859006
        Length Average:         13529.2587003388
        Length Median:          2249
        Min Length:             54499
        Max Length:             1999
        Length Stdev:           31939.4778282759

# Hmmm, so window length really changes the proportion of calls. Let's see how many novel calls there are.

intersectBed -a BIBR02.jarms.callscnvs.bed -b BIBR02.jarms.400cnvs.bed -v | bed_length_sum.pl
        Interval Numbers:       5142
        Total Length:           29343858
        Length Average:         5706.70128354726
        Length Median:          2999
        Min Length:             2999
        Max Length:             2499
        Length Stdev:           8718.15676037973

# Interesting. At least 50% of the calls are 3kb length, but there are still some differences
intersectBed -b BIBR02.jarms.callscnvs.bed -a BIBR02.jarms.400cnvs.bed -v | bed_length_sum.pl
        Interval Numbers:       8862
        Total Length:           63116338
        Length Average:         7122.13247573911
        Length Median:          9199
        Min Length:             27999
        Max Length:             9199
        Length Stdev:           13124.7575031819

# In this case, all of the calls are smaller

# I'm going to go with a flat 500bp window for everything. This way, I can compare the low coverage results against the higher coverage ones.
# Let's automate this for now and then take care of the issues later.
for i in `ls -d B?????`; do echo $i; ~/jdk1.8.0_05/bin/java -jar ~/JaRMS/store/JaRMS.jar call -i $i/${i}.rg.rd.full.sorted.merged.bam -f ../reference/umd3_kary_unmask_ngap.fa -o ${i}.auto.jarms -t 2 -w 500; done
```

The deadline is approaching, so I need to get some very preliminary results in order to pad my abstract. I think that I'm going to take the JaRMS calls that have finished for the Indicus animals and compare them to the RAPTR-SV results I previously generated (March 2 2015). My hope is to tie tandem dups with predicted duplication regions -- non-duplicated tandem repeat signatures indicate assembly shuffling.

Let's talk about the signatures I'm looking for:

* Genomic region rearrangement - tandem dups without RD signal, or trans-chr alignments
* MEI - transchr or maxdist repetitive reads
* Genomic deletion  deletion present in all individuals

OK, some summary statistics before I begin:

> Blade14: /mnt/iscsi/vnx_gliu_7/100_base_run

```bash
wc -l raptr-sv/finalfilter/BI*.tand.bed
   166 raptr-sv/finalfilter/BIBR02.finalfilter.raptr.tand.bed
   616 raptr-sv/finalfilter/BIBR03.finalfilter.raptr.tand.bed
   132 raptr-sv/finalfilter/BIBR04.finalfilter.raptr.tand.bed
   615 raptr-sv/finalfilter/BIBR05.finalfilter.raptr.tand.bed
  1447 raptr-sv/finalfilter/BIBR07.finalfilter.raptr.tand.bed
   268 raptr-sv/finalfilter/BIBR08.finalfilter.raptr.tand.bed
   130 raptr-sv/finalfilter/BIBR09.finalfilter.raptr.tand.bed
    75 raptr-sv/finalfilter/BIGI01.finalfilter.raptr.tand.bed
   348 raptr-sv/finalfilter/BIGI02.finalfilter.raptr.tand.bed
   415 raptr-sv/finalfilter/BIGI05.finalfilter.raptr.tand.bed
   109 raptr-sv/finalfilter/BIGI06.finalfilter.raptr.tand.bed
  1085 raptr-sv/finalfilter/BIGI07.finalfilter.raptr.tand.bed
   347 raptr-sv/finalfilter/BIGI08.finalfilter.raptr.tand.bed
  6723 raptr-sv/finalfilter/BINE01.finalfilter.raptr.tand.bed
   376 raptr-sv/finalfilter/BINE04.finalfilter.raptr.tand.bed
   187 raptr-sv/finalfilter/BINE07.finalfilter.raptr.tand.bed
   203 raptr-sv/finalfilter/BINE09.finalfilter.raptr.tand.bed
   352 raptr-sv/finalfilter/BINE10.finalfilter.raptr.tand.bed
    61 raptr-sv/finalfilter/BINE12.finalfilter.raptr.tand.bed
   187 raptr-sv/finalfilter/BINE13.finalfilter.raptr.tand.bed
   250 raptr-sv/finalfilter/BINE23.finalfilter.raptr.tand.bed
 14092 total

# So, BINE1 is pretty bad. Let's focus on Brahman for the abstract
# Let's go through my analysis pipeline
intersectBed -a raptr-sv/finalfilter/BIBR02.finalfilter.raptr.tand.bed -b BIBR02.auto.jarmscnvs.bed | wc -l
	23
grep 'duplication' BIBR02.auto.jarmscnvs.bed | wc -l
	5075
wc -l raptr-sv/finalfilter/BIBR02.finalfilter.raptr.tand.bed
	166 raptr-sv/finalfilter/BIBR02.finalfilter.raptr.tand.bed

# That leaves about 143 tandem events unaccounted for
# Let's confirm by checking their supporting read counts
# I noticed that "3" was the mode for this sample
intersectBed -a raptr-sv/finalfilter/BIBR02.finalfilter.raptr.tand.bed -b BIBR02.auto.jarmscnvs.bed -v | perl -lane 'if($F[4] > 3){print $_;}' | wc -l
	36 <- 36 candidate rearrangements

# Let's get the numbers for the rest of the samples
for i in BIBR02 BIBR03 BIBR04 BIBR05 BIBR07 BIBR08 BIBR09; do echo -en "$i\t"; intersectBed -a raptr-sv/finalfilter/${i}.finalfilter.raptr.tand.bed -b ${i}.auto.jarmscnvs.bed -v | perl -lane 'if($F[4] > 3){print $_;}' | wc -l; done
BIBR02  36
BIBR03  Error: The requested bed file (BIBR03.auto.jarmscnvs.bed) could not be opened. Exiting!
0
BIBR04  32
BIBR05  95
BIBR07  Error: The requested bed file (BIBR07.auto.jarmscnvs.bed) could not be opened. Exiting!
0
BIBR08  76
BIBR09  32

# Weird! The BIBR03 and BIBR07 files have idx files... oh well. I'll debug later.
# Let's print out the data to text files so I can intersect them
for i in BIBR02 BIBR03 BIBR04 BIBR05 BIBR07 BIBR08 BIBR09; do echo -en "$i\n"; intersectBed -a raptr-sv/finalfilter/${i}.finalfilter.raptr.tand.bed -b ${i}.auto.jarmscnvs.bed -v | perl -lane 'if($F[4] > 3){print $_;}' > temp.${i}.unsupported.tand; done

# Ok, I'm testing by grepping out chromosomes and comparing by eye. Here are some I've found
# chr2
temp.BIBR02.unsupported.tand:chr2       65185285        65185604        TANDEM  4.0
temp.BIBR04.unsupported.tand:chr2       65185274        65185584        TANDEM  4.0
temp.BIBR08.unsupported.tand:chr2       65185312        65185642        TANDEM  8.0

# chr4
temp.BIBR02.unsupported.tand:chr4       6117380 6117469 TANDEM  11.0
temp.BIBR04.unsupported.tand:chr4       6117387 6117469 TANDEM  5.0
temp.BIBR08.unsupported.tand:chr4       6117381 6117467 TANDEM  9.0

# chr10
temp.BIBR02.unsupported.tand:chr10      26147295        26148209        TANDEM  5.0
temp.BIBR04.unsupported.tand:chr10      26147248        26148178        TANDEM  6.0
temp.BIBR08.unsupported.tand:chr10      26147232        26148206        TANDEM  4.0
```

Hmm, I just realized that many of these were smaller regions that would not be detected by JaRMS at a 500bp window size anyways...

## JaRMS on Goat and Buffalo

This is my first-pass test of JaRMS on the Goat and Buffalo samples. Let's see how it performs.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
for i in `ls */*.bam`; do name=`echo $i | cut -d'/' -f1`; echo $name; ~/jdk1.8.0_05/bin/java -jar ~/JaRMS/store/JaRMS.jar call -i $i -f ../../reference/umd3_kary_unmask_ngap.fa -o ${name}.jarms -t 5 -w 500; done

wc -l *.bed
    5941 AG280.jarmscnvs.bed
    5875 AG304.jarmscnvs.bed
    8814 AG306.jarmscnvs.bed
   25795 ITWB10.jarmscnvs.bed
   26677 ITWB11.jarmscnvs.bed
   26853 ITWB12.jarmscnvs.bed
   24230 ITWB13.jarmscnvs.bed
   21375 ITWB14.jarmscnvs.bed
   25063 ITWB15.jarmscnvs.bed
   21298 ITWB1.jarmscnvs.bed
   23390 ITWB2.jarmscnvs.bed
   23708 ITWB3.jarmscnvs.bed
   23884 ITWB4.jarmscnvs.bed
   22783 ITWB5.jarmscnvs.bed
   17310 ITWB6.jarmscnvs.bed
   23141 ITWB7.jarmscnvs.bed
   24781 PC1.jarmscnvs.bed
  350918 total

# Quite a few more calls in Buffalo! Might be a dataset bias... but we can't focus on that for now
```

## RAPTR-SV on new Goat and Buffalo bams

I need to generate the RAPTR-SV calls for Goat and Buffalo. Let's do this in an automated fashion for now.

> Blade14: /mnt/iscsi/vnx_gliu_7/ruminant_project/goat_buff_bams

```bash
for i in `ls */*.bam`; do name=`echo $i | cut -d'/' -f1`; echo $name; ~/jdk1.8.0_05/bin/java -Xmx50g -jar ~/RAPTR-SV/store/RAPTR-SV.jar preprocess -i $i -o ${name}.raptr.preprocess -r ../../reference/umd3_kary_nmask_hgap.fa -t 10 -p ../../tmp/; done

```

## TFBS preparation

I want to see the impact of variants on TFBSs, so I better do the liftover from btau4 to umd3. Hopefully, I maintain at least 90% of the sites!

Liftover should be pretty easy to manage at this point. Here's the

> Blade14:  /mnt/iscsi/vnx_gliu_7/ruminant_project

```bash
/home/dbickhart/bin/liftOver cattle_bc_tfbs_btau4.bed bosTau4ToBosTau6.over.chain cattle_bc_tfbs_umd3_liftover.bed cattle_bc_tfbs_umd3_liftover.unmapped 

wc -l cattle_bc_tfbs_btau4.bed cattle_bc_tfbs_umd3_liftover.bed cattle_bc_tfbs_umd3_liftover.unmapped
  379333 cattle_bc_tfbs_btau4.bed
  377281 cattle_bc_tfbs_umd3_liftover.bed
    4104 cattle_bc_tfbs_umd3_liftover.unmapped

# Well! That's better than I hoped!
# Just taking a brief look at the unmapped coordinates quickly
grep '#Deleted' cattle_bc_tfbs_umd3_liftover.unmapped | wc -l
	1676
grep '#Partially' cattle_bc_tfbs_umd3_liftover.unmapped | wc -l
	349

# OK, I didn't want those sites in any case, it looks like!
```