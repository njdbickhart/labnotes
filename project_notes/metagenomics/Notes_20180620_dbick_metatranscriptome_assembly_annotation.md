# Metagenome assembly and classification
---
**6/20/1982**

These are my notes for generating an RNA-seq de novo assembly from Wenli's data and running classification on it.

## Table of Contents


## Guiding principles

Here is my plan:

* Assemble de novo transcripts from ALL data using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity).
* Annotate using HMMR3 with Pfam HMMs [ideas here:](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5594617/ )
* Annotate using blastp searches against nr database


## Generating the de novo assemblies

These are my notes for using Trinity and the filtered, processed fastq reads from Wenli's dataset. 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/Rumen_epimural_RNAseq/rumen_microbial_RNAseq

```bash
# Let's prepare the read files for input into Trinity
perl -e '@f = `ls *.fastq`; %h; chomp(@f); foreach my $c (@f){@csegs = split(/[_\.]/, $c); $h{$csegs[0]}->{$csegs[-2]} = $c;} $left = ""; $right = ""; foreach my $c (keys(%h)){ $left .= $h{$c}->{"R1"} . ","; $right .= $h{$c}->{"R2"} . ",";} chop($left); chop($right); print "--left $left --right $right\n";'

module load samtools bowtie2/2.3.0 jellyfish/2.2.3

# And now to run Trinity
# Running it with the jaccard-clip feature first (I am expecting ALLOT of genes in this data!)
sbatch --nodes=1 --ntasks-per-node=20 --mem=91000 -p assemble2 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/Trinityrnaseq-v2.6.6/Trinity --seqType fq --max_memory 90G --left rumen6775Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6773Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6771Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6768Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6793_2Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6766Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6765Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq,rumen6792_2Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R1.fastq --right rumen6775Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6773Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6771Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6768Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6793_2Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6766Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6765Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq,rumen6792_2Unmapped_against_Cattle_Unmapped_against_rRNA_database.out_R2.fastq --CPU 20 --jaccard_clip --output trinity_all_jaccard"
```