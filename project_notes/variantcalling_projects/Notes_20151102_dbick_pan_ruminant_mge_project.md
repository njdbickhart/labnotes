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