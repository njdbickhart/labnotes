# Bird CNV calling
---
*9/4/2015*

These are my commands and notes on the CNV calling on a wild bird dataset.

## Table of Contents
* [BAM file summary statistics](#stats)
* [Running RAPTR-SV](#raptr)

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
That should be enough to get started

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