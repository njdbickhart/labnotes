# Bird CNV calling
---
*9/4/2015*

These are my commands and notes on the CNV calling on a wild bird dataset.

## Table of Contents
* [BAM file summary statistics](#stats)
* [Annotation entries](#annotation)
* [Running RAPTR-SV](#raptr)
* [Running CNVNator on the bird data](#cnvnator)
* [Running JaRMS on the bird data](#jarms)
* [Summary of current results](#summary1)
* [Consolidating Tandem Dup information for Morgan](#consol)
	* [Strong Tandem Candidates](#strongtandem)
	* [Other Tandem Candidates](#othertandem)

<a name="stats"></a>
## BAM file summary statistics

Let's get a sense for what I'm working with here.

> Blade14: /mnt/nfs/nfs2/dbickhart/bird_data

```bash
samtools view -H sj.cor_trimmed_paired_mapping.bam | grep '@SQ' | wc -l
	34,314  <- chromosome count

# Unsorted bam
# Unindexed

# Let's make life easier by sorting and indexing
samtools sort -@ 5 -o sj.cor_trimmed_paired_mapping.sorted.bam -T sj.cor_temp samtools sort -@ 5 -o /mnt/iscsi/vnx_gliu_7/bird_data/sj.cor_trimmed_paired_mapping.sorted.bam -T /mnt/iscsi/vnx_gliu_7/bird_data/sj.cor_temp sj.cor_trimmed_paired_mapping.bam

# Damn, there was an error!
samtools view -H sj.cor_temp.1114.bam > sam.header
samtools merge -h sam.header sj.cor_trimmed_paired_mapping.sorted.bam $bams

# Damn! It's the file handle limit that's screwing things up here!
bams=`for i in $(seq 1 1020); do num=$(printf '%04d' $i); echo -n "sj.cor_temp.${num}.bam "; done; echo`
samtools merge -h sam.header sj.cor_trimmed_paired_mapping.sorted.1020.bam $bams

for i in $(seq 1 1020); do num=$(printf '%04d' $i); rm sj.cor_temp.${num}.bam; done
samtools merge -h sam.header sj.cor_trimmed_paired_mapping.sorted.rest.bam sj.cor_temp.*.bam

samtools merge -h sam.header sj.cor_trimmed_paired_mapping.sorted.final.bam sj.cor_trimmed_paired_mapping.sorted.1020.bam sj.cor_trimmed_paired_mapping.sorted.rest.bam

samtools index sj.cor_trimmed_paired_mapping.sorted.final.bam

```

Let's look at the gff file. I need to look at the gene specific info first.

```bash
# I need to grep out the GFF data as this is a combo GFF + fasta
grep -n FASTA MorganFinal1KbReplicate_v3.gff
	5205233:##FASTA
head -n 5205232 MorganFinal1KbReplicate_v3.gff > MorganFinal1KbReplicate_v3.gff.real.gff

perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -e '#' -f MorganFinal1KbReplicate_v3.gff.real.gff -c 1 -m
```

|Entry           |   Count|
|:---------------|-------:|
|.               |   34314|
|augustus_masked |  191578|
|blastx          | 1071841|
|cdna2genome     |  336284|
|maker           |  326989|
|protein2genome  |  875969|
|repeatmasker    | 1271694|
|repeatrunner    |     606|
|tblastx         | 1025406|

Maker and augustus are the annotation programs. I'm guessing that Chris also did a blastx and used two transcript assembly comparison programs (cdna2genome & protein2genome) to identify gene regions. Finally, he used a repeatmasker profile to identify repetitive elements. 

First, let's get the fasta sequence because I need that for RAPTR-SV

```bash
# extracting the fasta sequence
perl -e 'my $start = 0; while(<>){if($_ =~ /\##FASTA/){$start = 1;}elsif($start){print $_;}}' < MorganFinal1KbReplicate_v3.gff > MorganFinal1KbReplicate_v3.fa

# Gzipping it for transfer
gzip MorganFinal1KbReplicate_v3.fa
```
Now, let's pull out the repeatmasker entries in the gff and use that to mask the fasta.

> blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
perl -lane 'if($F[0] =~ /^#/ || $F[0] eq ""){next;} if($F[1] eq "repeatmasker"){print "$F[0]\t$F[3]\t$F[4]";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | sortBedFileSTDIN.pl | mergeBed -i stdin > fasta/sj.corr.repeatmasker.bed

wc -l fasta/sj.corr.repeatmasker.bed
	629631 fasta/sj.corr.repeatmasker.bed

# masking the fasta using the repeatmasker information
~/bedtools-2.17.0/bin/maskFastaFromBed -fi fasta/MorganFinal1KbReplicate_v3.fa -bed fasta/sj.corr.repeatmasker.bed -fo fasta/MorganFinal1KbReplicate_v3.repeatmasked.fa

# Indexing with mrsfast
mrsfast --index fasta/MorganFinal1KbReplicate_v3.repeatmasked.fa

```
<a name="annotation"></a>
## Annotation entries

I need to pull the annotation from the gff, but it looks like Chris used several different prediction methods. Thankfully, my annotation program will let me use each one individually. Let's separate them from the GFF file into individual BED files.

> blade14: /mnt/iscsi/vnx_gliu_7/bird_data/

```bash
mkdir annotation

# Now to pull individual sets from the second column identifier
# Maker is going to be a super-set of the other methods, most likely.
# I've noticed that augustus has "match" and "match_part" annotations
perl -lane 'if($F[1] eq "augustus_masked" && $F[2] eq "match"){print $_;}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | wc -l
	21,301 <- match annotations
perl -lane 'if($F[1] eq "augustus_masked" && $F[2] eq "match_part"){print $_;}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | wc -l
	170,277 <- match part annotations
perl -lane 'if($F[1] eq "augustus_masked" && $F[2] eq "expressed_sequence_match"){print $_;}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | wc -l
	0 <- expressed sequence match annotations

# Let's pull both match and matchpart annotations and put them in the same bed
# Ok, this is why I hate gffs. 
# The match_part entries have the ID within the "TARGET=" tag, and the match entries have the ID within the "NAME=" tag
perl -lane 'if($F[1] eq "augustus_masked" ){ $name; if($F[2] eq "match"){ ($name) = $F[8] =~ /Name=(.+)/;}elsif($F[2] eq "match_part"){($name) = $F[8] =~ /Target=(\S+)/;} print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | sortBedFileSTDIN.pl > annotation/augustus_morgan_annotation.bed

# now protein2genome
grep 'protein2genome' /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 2
	Entry   Count
	match_part      674697
	protein_match   201272

perl -lane 'if($F[1] eq "protein2genome" ){ $name; if($F[2] eq "protein_match"){ ($name) = $F[8] =~ /Name=(.+)/;}elsif($F[2] eq "match_part"){($name) = $F[8] =~ /Target=(\S+)/;} print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | sortBedFileSTDIN.pl > annotation/protein2genome_morgan_annotation.bed

# now cdna2genome
grep 'cdna2genome' /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 2                                                                                	Entry   Count
	expressed_sequence_match        51390
	match_part      284894

perl -lane 'if($F[1] eq "cdna2genome" ){ $name; if($F[2] eq "expressed_sequence_match"){ ($name) = $F[8] =~ /Name=(.+)/;}elsif($F[2] eq "match_part"){($name) = $F[8] =~ /Target=(\S+)/;} print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | sortBedFileSTDIN.pl > annotation/cdna2genome_morgan_annotation.bed

# Finally, maker
grep 'maker' /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 2
	Entry   Count
	CDS     144940
	exon    146180
	five_prime_UTR  777
	gene    17132
	mRNA    17132
	three_prime_UTR 828

# OK, I want to be selective here. I think that the gene ID is going to be best here
perl -lane 'if($F[1] eq "maker" && $F[2] eq "gene"){ ($name) = $F[8] =~ /Name=(.+)/; print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | sortBedFileSTDIN.pl > annotation/maker_gene_morgan_annotation.bed
# NOTE: the maker gene IDs are unique. This is key for my annotation program

# Let's see if we can get any info from the blastx searches
grep 'blastx' /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | grep -v 'tblastx' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 2
	Entry   Count
	match_part      903631
	protein_match   168210
grep 'tblastx' /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 2
	Entry   Count
	match_part      975400
	translated_nucleotide_match     50006

# Let's only pull out the complete entries from these two files
perl -lane 'if($F[1] eq "blastx" && $F[2] eq "protein_match"){ ($name) = $F[8] =~ /Name=(.+)/; print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff > annotation/blastx_morgan_annotation.bed
perl -lane 'if($F[1] eq "tblastx" && $F[2] eq "translated_nucleotide_match"){ ($name) = $F[8] =~ /Name=(.+)/; print "$F[0]\t$F[3]\t$F[4]\t$name";}' < /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.gff.real.gff > annotation/tblastx_morgan_annotation.bed


```

<a name="raptr"></a>
## Running RAPTR-SV

I need to prepare everything first, run the preprocess step and then identify if I need to scale this back because of the shear number of scaffolds.


#### Setting up necessary files
> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
# OK, setting up the working directory
mkdir raptr

# Now to set up the reference genome fasta directory
mkdir fasta
cp /mnt/nfs/nfs2/dbickhart/bird_data/MorganFinal1KbReplicate_v3.fa.gz ./fasta/
gunzip fasta/MorganFinal1KbReplicate_v3.fa.gz

# moving into the fasta directory briefly to index the fasta with mrsfast
mrsfast --index MorganFinal1KbReplicate_v3.fa

# there were some gaps in the scaffolds -- identifying them for the clustering step
~/jdk1.8.0_05/bin/java -jar ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -f MorganFinal1KbReplicate_v3.fa -o MorganFinal1KbReplicate_v3.gaps.bed -s MorganFinal1KbReplicate_v3.gaps.stats

perl -e '$c = 0; $s = 0; <>;  while(<>){chomp; @s = split(/\t/); $c++; $s += $s[7];} print ($s / $c); print "\n";' < MorganFinal1KbReplicate_v3.gaps.stats
	11.1033588494158	<- 11.1% of the assembly, on average, is comprised of gaps
```

One note from the gap analysis: if needed, I can remove most scaffolds that are mostly comprised of gaps if the number of scaffolds is too high.

<a name="raptr"></a>
#### Running the RAPTR-SV pipeline
> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
~/jdk1.8.0_05/bin/java -Xmx40g -jar ~/RAPTR-SV/store/RAPTR-SV.jar preprocess -i sj.cor_trimmed_paired_mapping.sorted.final.bam -o raptr/sj.cor_trimmed.raptr.preprocess -r fasta/MorganFinal1KbReplicate_v3.fa -g -t 10 -p ../tmp

# Damn! Samtools HTSLIB merger generated thousands of read groups. I'm going to treat both libraries as one for now
~/jdk1.8.0_05/bin/java -Xmx40g -jar ~/RAPTR-SV/store/RAPTR-SV.jar preprocess -i sj.cor_trimmed_paired_mapping.sord.final.bam -o raptr/sj.cor_trimmed.raptr.preprocess -r fasta/MorganFinal1KbReplicate_v3.fa -t 10 -p ../tmp
	[MAIN] Setting temporary file directory to: ../tmp
	[MAIN] Setting ForkJoin thread ceiling to: 10
	[METADATA] Identified more than one read length in BAM! Using most frequent sampled read length to filter split reads: 124
	[PREPROCESS] Read input file and calculated sample thresholds.
	Sample: D Avg Ins size: 1648.8695 Stdev Ins size: 1918.163622644829

# Hmm... read length sampling suggests way too much of a distribution stdev. I might have to clean up the header afterall
# I'm going to pull the nightly build of samtools from github and try using that for the merger instead.
# NOTE: I needed to pull the nightly build of htslib as well!
samtools sort -m 2G -o sj.cor_trimmed_paired_mapping.resorted.bam -T sj.cor.temp -@ 10 /mnt/nfs/nfs2/dbickhart/bird_data/sj.cor_trimmed_paired_mapping.bam

# rerunning on the resorted, repeatmasked data
~/jdk1.8.0_05/bin/java -jar ~/RAPTR-SV/store/RAPTR-SV.jar preprocess -i sj.cor_trimmed_paired_mapping.resorted.bam -o raptr/sj.cor_resorted.raptr.preprocess -r fasta/MorganFinal1KbReplicate_v3.repeatmasked.fa -g -t 10 -p ../tmp
	Sample: 3436b63b-fe72-471e-8018-58e591fe2e60 Avg Ins size: 3401.7501 Stdev Ins size: 1781.3147379026404
	Sample: e83dd82e-8ff0-4ea6-9b38-25c34eef94e8 Avg Ins size: 384.0294 Stdev Ins size: 179.4941774866667

```

There was a problem with the split read alignment portion -- I think that there were NO alignments selected for reprocessing for some reason. I'll create a dummy bam file and proceed regardless.

```bash
samtools view -b -H raptr/sj.cor_resorted.raptr.preprocess.e83dd82e-8ff0-4ea6-9b38-25c34eef94e8.anchor.bam > raptr/sj.cor_resorted.raptr.preprocess.e83dd82e-8ff0-4ea6-9b38-25c34eef94e8.dummysplit.bam
# Then I added the split read file to the preprocess list, and removed the other read group from consideration

# OK, going to attempt this without the split reads
~/jdk1.8.0_05/bin/java -jar ~/RAPTR-SV/store/RAPTR-SV.jar cluster -s raptr/sj.cor_resorted.raptr.preprocess.flat -g fasta/MorganFinal1KbReplicate_v3.gaps.bed -o raptr/sj.cor_resorted.raptr.cluster -t 10 -p ../tmp/

```
<a name="cnvnator"></a>
#### Running CNVNator on the bird data

I was able to get CNVNator to work! I had to compile within my cnvnator directory on blade14, and here are the environmental variables that I needed to set:

> $LD_LIBRARY_PATH = /usr/lib64/:/usr/lib64/root/

> $ROOTSYS = /usr/lib64/root/

Because the system is "headless" with respect to X11 libraries, I cannot install a default ROOT package. Since CNVNator doesn't really use the X11 libraries from ROOT, these are optional to include, but the install fails early if it doesn't find them. Now to use my newly compiled CNVNator executable to call cnvs on this bam.

> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
# Creating histogram
~/CNVnator_v0.2.7/src/cnvnator -root cnvnator/sj.corr.combined.root -genome fasta/MorganFinal1KbReplicate_v3.fa -tree sj.cor_trimmed_paired_mapping.resorted.bam

# Damn, I'm getting a segmentation fault and a stack trace
The lines below might hint at the cause of the crash.
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
#5  0x00000031cb48636a in strlen () from /usr/lib64/libc.so.6
#6  0x0000003ad8195a85 in TString::TString(char const*) () from /usr/lib64/root/libCore.so.5.34
#7  0x00007f76020a6305 in TH1::FFT(TH1*, char const*) () from /usr/lib64/root/libHist.so.5.34
#8  0x0000000000423fd4 in HisMaker::produceTrees(std::string*, int, std::string*, int, bool) ()
#9  0x0000000000407f04 in main ()
===========================================================

```

<a name="jarms"></a>
## JaRMS data collection

OK, I got frustrated with CNVNator and I wrote my own program. Let's see how this works.

> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
# Threads aren't fully implemented yet in this version, but I don't want to limit the innate fork-join pool
~/jdk1.8.0_05/bin/java -Xmx56g -jar ~/JaRMS/store/JaRMS.jar call -i sj.cor_trimmed_paired_mapping.resorted.bam -f fasta/MorganFinal1KbReplicate_v3.repeatmasked.fa -o sj.cor.jarms.calls.bed -t 10

# NOTE: I need better temp file directory handling. I use the output directory path literal instead of a proper directory
```

Damn, I think that the program hangs on contigs that have zero reads (ie. the nearly completely repeatmasked regions). I was able to parse my binary file with the GC corrected read depth, so let's try to call CNVs off of the data we've got.

> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data

```bash
#	Let's get the mean and stdev of the GC corrected data
# From the log file:
grep 'mean' JaRMS.call.2015_10_19_12_43_09.0.0.log
	Oct 19,2015 12:43 -- INFO -- jarms.modes.CallMode -> SingleThreadRun -- [CALLMODE] Unadjusted RD values: mean-> 1259.5969244443722 sd-> 5427.133162033951

# OK, CNVnator takes 1/2 the mean and 2 * the mean to get the values needed for the stdev estimate
# Because of several huge outliers, I'm going to do that as well
perl -lane 'if($F[3] < (1259 / 2) || $F[3] > (1259 * 2)){next;}else{ print $F[3];}' < sj.cor.jarms.calls.bed.gccorr.corr.bed | statStd.pl
total   1655880
	Minimum 629.5015845566844
	Maximum 2517.9903990167227
	Average 1261.007442
	Median  1186.76349951635
	Standard Deviation      391.614905
	Mode(Highest Distributed Value) 998.0583212449338

# OK, let's now try to co-opt Can Alkan's scripts to get the CNV windows
# Duplication cutoff = 2435.837442
perl ~/wssd-package/wssd_picker.pl -f sj.cor.jarms.calls.bed.gccorr.corr.bed -w 7 -s 6 -c 2435 -b 3 -n 1 -o sj.cor.jarms.calls.gccorr.dups.bed
wc -l sj.cor.jarms.calls.gccorr.dups.bed
	273 sj.cor.jarms.calls.gccorr.dups.bed

# OK! Let's see how many bases this comprises
perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[2] - $s[1];} print "$c\n";' < sj.cor.jarms.calls.gccorr.dups.bed
1,146,647	<- about 1% of the data

# Now for the deletions
# Deletion cutoff = 479
perl ~/wssd-package/wssd_picker.pl -f sj.cor.jarms.calls.bed.gccorr.corr.bed -w 7 -s 6 -c 479 -m -b 3 -n 1 -o sj.cor.jarms.calls.gccorr.dels.bed
wc -l sj.cor.jarms.calls.gccorr.dels.bed
	11136 sj.cor.jarms.calls.gccorr.dels.bed

# I knew there would be quite a few. Let's remove the repeat-masked regions first
subtractBed -a sj.cor.jarms.calls.gccorr.dels.bed -b fasta/sj.corr.repeatmasker.bed | perl -lane 'if($F[2] - $F[1] < 200){next;}else{print $_;}' > sj.cor.jarms.calls.gccorr.dels.norpt.bed
wc -l sj.cor.jarms.calls.gccorr.dels.norpt.bed
	22982 sj.cor.jarms.calls.gccorr.dels.norpt.bed

cat sj.cor.jarms.calls.gccorr.dels.norpt.bed | perl ~/bin/bed_length_sum.pl
        Interval Numbers:       22982
        Total Length:           36816588
        Length Average:         1601.97493690714
        Length Median:          1084
        Min Length:             2926
        Max Length:             3999
        Length Stdev:           1556.22709580478

# Let's see how these stack up compared to the sex chromosomes
cat sj.cor.jarms.calls.gccorr.dels.norpt.bed | cut -f1 | sort | uniq > deletion_chrs.txt
wc -l deletion_chrs.txt
	7070 deletion_chrs.txt

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl deletion_chrs.txt /mnt/nfs/nfs2/dbickhart/bird_data/Morgan_ZW_scaffolds.csv
File Number 1: deletion_chrs.txt
File Number 2: /mnt/nfs/nfs2/dbickhart/bird_data/Morgan_ZW_scaffolds.csv
Set     Count
1       6608
1;2     462
2       1088

# Not very convincing
# Let's compare this to the dups
cat sj.cor.jarms.calls.gccorr.dups.bed | cut -f1 | sort | uniq > dup_chrs.txt
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl dup_chrs.txt /mnt/nfs/nfs2/dbickhart/bird_data/Morgan_ZW_scaffolds.csv
File Number 1: dup_chrs.txt
File Number 2: /mnt/nfs/nfs2/dbickhart/bird_data/Morgan_ZW_scaffolds.csv
Set     Count
1       215
1;2     15
2       1535

# Fewer, at least!
# I've thought of a way to estimate the likelihood of RD windows being derived from the Z/W chromosome scaffolds
# First, let's get the statistics on the ZW scaffolds
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]}) && $s[3] > 0){print join("\t", @s); print "\n";}} close IN;' /mnt/nfs/nfs2/dbickhart/bird_data/Morgan_ZW_scaffolds.csv sj.cor.jarms.calls.bed.gccorr.corr.bed > putative_zw_scaffolds.gccorr.rd.bed

cat putative_zw_scaffolds.gccorr.rd.bed | cut -f4 | statStd.pl
	total   140710
	Minimum 0.8756060903062457
	Maximum 3621139.9085821067
	Average 909.008232
	Median  552.614753483193
	Standard Deviation      17523.841203
	Mode(Highest Distributed Value) 15.109298926334505

# The average value is lower, but there are some contigs that have significantly heightened estimates of read depth.
# These are likely repetitive fragments common to the sex chromosomes -- if birds are anything like mammals here.
```

Just to give Chris something to use, I'm going to intersect the duplication regions (ie. the WSSDs) with gene locations.

```bash
# Annotation of WSSD regions
~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d annotation/genedb.list -i sj.cor.jarms.calls.gccorr.dups.bed -o sj.cor.jarms.wssd
~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d annotation/genedb.list -i sj.cor.jarms.calls.gccorr.dups.bed -o sj.cor.jarms.wssd -t

# Now, annotation of the deletions...
~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d annotation/genedb.list -i sj.cor.jarms.calls.gccorr.dels.norpt.bed -o sj.cor.jarms.dels
~/jdk1.8.0_05/bin/java -jar ~/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d annotation/genedb.list -i sj.cor.jarms.calls.gccorr.dels.norpt.bed -o sj.cor.jarms.dels -t
```
<a name="summary1"></a>
## Summary of current results

#### Methods
All data was processed in the following stepwise manner:

1. Created a set of non-overlapping, 500 bp windows for each scaffold
2. Calculated global mean and stdev from a fitted Guassian distribution
3. Estimated G+C % effect on read depth values
4. Normalized G+C % effects on read depth values
5. Reestimated mean and stdev -- fitted values from simplistic mean cutoff
6. Established cutoffs for duplications and deletions
	1. Duplication cutoff read depth: > 2435
	2. Deletion cutoff read depth: < 479
7. Called CNVs if, out of 7 consecutive windows, 6 windows were higher or lower than the cutoffs
8. Removed Deletion windows that intersected with repeat-masked regions
8. Annotated CNVs with gene lists that were provided in the original gff file

#### List of files and descriptions

| File | Description|
| :--- | :--- |
sj.cor.jarms.dels_regions.tab | tab delimited list of deletions or potential sites of misassembly/misalignment
sj.cor.jarms.wssd_regions.tab | tab delimited list of predicted Segmental Duplications
sj.cor.jarms.wssd_anno.xls | excel spreadsheet with gene intersections of the Segmental Duplications
sj.cor.jarms.dels_anno.xls | excel spreadsheet with gene intersections of the misassemblies/deletions

<a name="consol"></a>
## Consolidating Tandem Dup information for Morgan

I'm now going to take the RAPTR-SV divet file and attempt to find evidence of tandem duplications near previously detected WSSD in Morgan that match up with Joe (the other sample). My strategy is pretty straight-forward:

* Convert RAPTR-SV divet file into a bedpe and intersect bed coordinates of remapped WSSD segments with everted bedpe
* Count the support and define tandem dup breakpoints

Let's subsection the divet so we can rapidly interrogate it.

> Blade14: /mnt/iscsi/vnx_gliu_7/bird_data/raptr

```bash
# Selecting only the scaffolds that we need
grep 'eversion' sj.cor_resorted.raptr.preprocess.3436b63b-fe72-471e-8018-58e591fe2e60.divet | perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %keep; while(<IN>){chomp; $keep{$_} = 1;} close IN; while(<STDIN>){chomp; @s = split(/\t/); if(exists($keep{$s[1]})){print "$s[1]\t$s[2]\t$s[3]\t$s[5]\t$s[6]\t$s[7]\t$s[0]\n";}}' scaffs.list > sj.cor.wssd.firstselection.eversion.bedpe

wc -l *.bedpe
	95815 sj.cor.wssd.firstselection.eversion.bedpe 

# Well, that's a reasonable number! Let's see how close it compares to our bed coordinates from the WSSD blasting
# First, let's generate some end sequence
perl -lane '$s1 = $F[1] - 100; $e1 = $F[1] + 100; if($s1 < 1){$s1 = 1;} $s2 = $F[2] - 100; $e2 = $F[2] + 100; print "$F[0]\t$s1\t$e1\n$F[0]\t$s2\t$e2";' < sj.cor_refined_wssd.bed > sj.cor_refined_wssd.ends.bed
perl -lane '$s1 = $F[1] - 100; $e1 = $F[1] + 100; if($s1 < 1){$s1 = 1;} $s2 = $F[2] - 100; $e2 = $F[2] + 100; print "$F[0]\t$s1\t$e1\n$F[0]\t$s2\t$e2";' < sj.cor_original_wssd.bed > sj.cor_original_wssd.ends.bed


# After a few successive tests, I'm worried about repeats. Some of the "ends" are in repetitive regions
# Let's still gather the evidence and deal with it later
intersectBed -a sj.cor_original_wssd.ends.bed -b sj.cor.wssd.firstselection.eversion.bedpe -c | perl -lane 'if($F[3]){print $_;}'
	jcf7180005232807        11401   11601   2
>	jcf7180005232807        14400   14600   221
	jcf7180005232807        65401   65601   382
	jcf7180005232807        65401   65601   382
	jcf7180005232807        65401   65601   382
	jcf7180005234354        215401  215601  2
	jcf7180005234354        215401  215601  2
	jcf7180005234811        196901  197101  188
	jcf7180005234811        200400  200600  4
	jcf7180005219577        9400    9600    1
	jcf7180005234681        16401   16601   19
	jcf7180005234681        19900   20100   7
	jcf7180005234416        470901  471101  13
	jcf7180005232807        6900    7100    2
	jcf7180005198788        1       101     194
	jcf7180005198788        1       101     194

intersectBed -a sj.cor_refined_wssd.ends.bed -b sj.cor.wssd.firstselection.eversion.bedpe -c | perl -lane 'if($F[3]){print $_;}'
	jcf7180005232807        11401   11601   2
	jcf7180005232807        13500   13700   5
	jcf7180005232807        65841   66041   10
>	jcf7180005232807        66577   66777   338
	jcf7180005232807        72282   72482   174
	jcf7180005232807        72342   72542   174
	jcf7180005234675        10410   10610   4
	jcf7180005234811        198682  198882  4
	jcf7180005198788        9153    9353    793
	jcf7180005198788        9620    9820    671
	jcf7180005232807        5688    5888    37
	jcf7180005198788        1       101     194
	jcf7180005198788        5436    5636    201
	jcf7180005198788        5875    6075    411

# Now to eliminate putative repetitive regions
intersectBed -a sj.cor_original_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005232807'
	jcf7180005232807        14568   14600   jcf7180005232807        14568   14761
intersectBed -a sj.cor_refined_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005232807'
	jcf7180005232807        65903   65923   jcf7180005232807        65903   65923

intersectBed -a sj.cor_refined_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005234811'
intersectBed -a sj.cor_original_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005234811'
# jcf7180005234811 is clean

intersectBed -a sj.cor_original_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005198788'
	jcf7180005198788        8455    8601    jcf7180005198788        8455    8653
intersectBed -a sj.cor_refined_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005198788'
# Interesting fact: 388 reads support this tandem dup:
# jcf7180005198788        73      194     jcf7180005198788        9176    9209

intersectBed -a sj.cor_original_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005234681'
intersectBed -a sj.cor_refined_wssd.ends.bed -b ../fasta/sj.corr.repeatmasker.bed -wb | grep 'jcf7180005234681'
# jcf7180005234681 is clean

# Just to test this, going to generate some BAM slices and view them with SVS Genome browse
```

That makes three potential Tandem Dups with good supporting evidence:
<a name="strongtandem"></a>
#### Strong Tandem candidates
*	jcf7180005234811        197001  200500
*	jcf7180005198788	73	9209                    <- IGV signature suggests tandem dup + dispersed dup near 5,608
*	jcf7180005234681        16501   20000


Now I'm going to try a more aggressive approach that includes the majority of the window interval; however, I am going to actively filter repetitive regions.

```bash
intersectBed -a sj.cor.wssd.firstselection.eversion.bedpe -b ../fasta/sj.corr.repeatmasker.bed -v > sj.cor.wssd.firstselection.eversion.norepeats.bedpe
wc -l *.bedpe
   95815 sj.cor.wssd.firstselection.eversion.bedpe
   92099 sj.cor.wssd.firstselection.eversion.norepeats.bedpe

# That removed far fewer than I expected!
# Let's see how the data intersects now
# NOTE: the refined coordinates were a subset of the following selection
intersectBed -a sj.cor_original_wssd.bed -b sj.cor.wssd.firstselection.eversion.norepeats.bedpe -c | perl -lane 'if($F[3]){print $_;}' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -m
intersectBed -a sj.cor_original_wssd.ends.bed -b sj.cor.wssd.firstselection.eversion.norepeats.bedpe -c | perl -lane 'if($F[3]){print $_;}' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -m
```
<a name="othertandem"></a>
#### Other Tandem Candidates
#### Original 
|Entry            | Count|
|:----------------|-----:|
|jcf7180005198788 |     3|
|jcf7180005219577 |     1|
|jcf7180005228792 |     2|
|jcf7180005232807 |     5|
|jcf7180005234354 |     2|
|jcf7180005234416 |     1|
|jcf7180005234601 |     1|
|jcf7180005234675 |     2|
|jcf7180005234681 |     1|
|jcf7180005234717 |     1|
|jcf7180005234811 |     1|

#### Original ends
|Entry            | Count|
|:----------------|-----:|
|jcf7180005198788 |     2|
|jcf7180005219577 |     1|
|jcf7180005232807 |     6|
|jcf7180005234354 |     2|
|jcf7180005234416 |     1|
|jcf7180005234681 |     2|
|jcf7180005234811 |     2|

So, there is strong evidence for the prior three tandem dups I detected, and some very complex evidence for the tandem dups in the table above (minus the three which are in both lists).