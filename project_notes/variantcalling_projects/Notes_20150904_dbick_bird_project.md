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
```
