# Bird CNV calling
---
*9/4/2015*

These are my commands and notes on the CNV calling on a wild bird dataset.

## Table of Contents
* [BAM file summary statistics](#stats)

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
