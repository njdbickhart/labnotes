# Cattle assembly manuscript work
---
*1/29/2019*

These are my notes on work done to support the cattle assembly manuscript.

## Table of Contents


## Summary comparison statistics

I want to draw a few brief comparisons between assemblies in order to highlight the differences. My first priority is establishing the rank order correlation between the optical maps and the cattle assemblies.

### Rank order comparisons

The cattle optical maps were aligned using the Bionano proprietary software, I think. The alignment files are located here on Ceres:

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/OpticalMap

First order of business is to identify simple counts of alignments and blocks.

```bash
wc -l *.aln
   18674 OM-ARS-UCD1.2.full.aln
   18629 OM-UMD3.1.1.full.aln

# I think that Ben included the unplaced contigs as well, so this makes things very difficult.

cat OM-ARS-UCD1.2.full.aln | cut -f1 | sort | uniq | wc -l
18674	<- queries = the optical map contigs
cat OM-ARS-UCD1.2.full.aln | cut -f2 | sort | uniq | wc -l
39		<- reference = assembled chromosomes

# OK, this is a bad situation. The OM was split into windows prior to mapping
# So, I need to take into account the OM contig as well as the window numbers

# I wrote a script to generate rank order from these alignment files
# It uses the reference chromosomes as the basis for the preliminary sort, so it assumes that the assembly scaffolds are larger than the OM scaffolds
sbatch sortAndRankOMData.pl OM-ARS-UCD1.2.full.aln OM-ARS-UCD1.2.full.ranks OM-ARS-UCD1.2.full.om.consensus
sbatch sortAndRankOMData.pl OM-UMD3.1.1.full.aln OM-UMD3.1.1.full.ranks OM-UMD3.1.1.full.om.consensus
```

Now to check the rank order correlations.

```R
# ARS UCD rank order test
data <- read.delim("OM-ARS-UCD1.2.full.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 6.7575e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.3772741

# UMD rank order test
data <- read.delim("OM-UMD3.1.1.full.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 9.9289e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
0.07837418
```

These are atrocious for both assemblies! I could determine ref chromosome consensus order, but I'm not sure if that would even help at this point! I need to see how disjointed the alignments are between optical maps first. 

Here's the consensus approach:

```bash
sbatch sortAndRankOMDataCons.pl OM-ARS-UCD1.2.full.aln OM-ARS-UCD1.2.cons.ranks OM-ARS-UCD1.2.cons.om.consensus
sbatch sortAndRankOMDataCons.pl OM-UMD3.1.1.full.aln OM-UMD3.1.1.cons.ranks OM-UMD3.1.1.cons.om.consensus
```

```R
data <- read.delim("OM-ARS-UCD1.2.cons.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 2.1526e+10, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.9801629

data <- read.delim("OM-UMD3.1.1.cons.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 4.1647e+10, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.9613426
```

OK, so that was a huge improvement. The magnitude of misalignments was higher than expected, so that threw off my ranking of the OM windows. 

## RepeatMasker analysis

I am going to try to run RepeatMasker on the cattle reference using the default CERES module (no special libraries). I know that this will upset allot of people, but this is exploratory for now.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/repeatmasker

```bash
sbatch --nodes=1 --mem=30000 --ntasks-per-node=10 -p medium --wrap="RepeatMasker -pa 10 -q -species cow -no_is -gff ../ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna"
# Building species libraries in: /home/derek.bickharhth/.RepeatMaskerCache/dc20170127/bos_taurus
#   - 175 ancestral and ubiquitous sequence(s) for bos taurus   <- ruh-roh! That ain't good!
#   - 0 lineage specific sequence(s) for bos taurus
```

## Unique sequence in the cattle assembly

OK, now I'm going to start scraping the reads from the cattle assembly. I'm going to be using [this manuscript](https://www.nature.com/articles/s41588-018-0273-y#MOESM1) as a guide for the workflow.

Here's their generic command listing so that I don't have to go into the supplemental every time I want to check on how to queue this:

```bash
---------------------------------------
Bowtie 2 alignment, per sample
---------------------------------------
bowtie2-build [GRCh38_no_alt] [GRCh38_no_alt_idx]
bowtie2 -x [GRCh38_no_alt_idx] [reads1] [reads2] > [alignments.bam]
---------------------------------------------------------------------------------
Extraction of unaligned reads (and mates) via samtools, per sample
---------------------------------------------------------------------------------
samtools fastq –f 12 [alignments.bam] -1 [mateUnmapped_R1.fq] -2 [mateUnmapped_R2.fq]
samtools fastq –f 68 –F 8 [alignments.bam] > [R1_mateMapped.fq]
samtools fastq –f 132 –F 8 [alignments.bam] > [R2_mateMapped.fq]
samtools view –f 8 –F 4 [alignments.bam] > [GRCh38Links.bam]
-----------------------------------------
MaSuRCA assembly, per sample
-----------------------------------------
masurca_config.txt: **********************************************************
DATA
PE= pe 300 50 [mateUnmapped_R1.fq] [mateUnmapped_R2.fq]
PE= s1 300 50 [R1_mateMapped.fq]
PE= s2 300 50 [R2_mateMapped.fq]
END
PARAMETERS
GRAPH_KMER_SIZE=auto
USE_LINKING_MATES=1
KMER_COUNT_THRESHOLD = 1
NUM_THREADS=24
JF_SIZE=200000000
DO_HOMOPOLYMER_TRIM=0
END
*********************************************************
masurca masurca_config.txt && ./assemble.sh
----------------------------
Centrifuge, per sample
----------------------------
centrifuge --report-file [centrifuge.report] -x [centrifugedb] -k 1 --host-taxids 9606
-f [masurca_contigs_over1kb.fa] > [centrifuge.output]
centrifuge-kreport -x [centrifugedb] [centrifuge.output] --min-score 0 --min-length 0
> [centrifuge.krakenOut]
*** centrifuge.krakenOut was used to filter any non-chordate identified reads. ***
---------------------------------------------------------
RepeatMasker on assembly contigs, per sample
--------------------------------------------------------- 
RepeatMasker -nolow -species human [filteredContigs.fa]
-------------------------------------------------------------
Bowtie 2 alignment of reads to contigs, per sample
-------------------------------------------------------------
bowtie2-build [filteredContigs.fa.masked] [contigIdx]
bowtie2 -x [contigIdx] -U [mateUnmapped_R2.fq],[mateUnmapped_R2.fq] -S
[readContigAlignment.sam]
-------------------------------------------------------------------------------------
Linking mates to implicated region, and aligning to region, per sample
-------------------------------------------------------------------------------------
samtools view -h -F 256 [readContigAlignment.sam] | samtools sort - -n –O bam |
bedtools bamtobed -i stdin | awk '{OFS=”\t”}{print $4,$1,$6,$2,$3}' | sort >
[readContigAlignment.txt]
samtools view -H [GRCh38Links.bam] | cat - <(awk 'FNR==NR{main[$1]=$0;next} $1 in main
{print main[$1]}' <(samtools view [GRCh38Links.bam]) [readContigAlignment.txt]) |
samtools sort -n –O bam | bedtools bamtobed -i stdin | awk '{OFS=”\t”}{print
$4,$1,$6,$2,$3}' | sed -e 's/\/[1-2]//g' | sort > [matchedMates.txt]
join -j 1 [readContigAlignment.txt] [matchedMates.txt] > [mateLinks.txt]
*** Filtering was performed here using python scripts to examine links to contig ends only, and filter based on
described unambiguity criteria (see methods). Contig ends and GRCh38 regions meeting criteria were extracted with
samtools faidx ***
nucmer --maxmatch -l 15 -b 1 -c 15 -p [deltaFile] [GRCh38Regions.fa]
[filteredContigEnds.fa]
----------------------------------
Clustering of placed contigs
----------------------------------
bedtools merge -d 100 -c 4 -o distinct [placedCtgLocations.bed] > [mergedClusters.bed]
nucmer -p [deltaFile] [repCtg.fa] [restOfClusterCtgs.fa]
nucmer -p [deltaFile] [verifiedClusterCtgs.fa] [unplacedCtgs.fa]
----------------------------------------------------------------------------
Left/Right one end placement merging into two end placement
----------------------------------------------------------------------------
nucmer --maxmatch --nosimplify -p [deltaFile] [leftEndedPlaced.fa] [rightEndPlaced.fa]
show-coords -H -T -l -c -o [deltaFile] > [coordsFile]
------------------------------------------
Removal of redundant placements
------------------------------------------
nucmer --maxmatch --nosimplify -p [deltaFile] [allPlaced.fa] [allPlaced.fa]
-------------------------------------
Clustering of unplaced contigs
-------------------------------------
nucmer --maxmatch --nosimplify -l 31 -c 100 -p [deltaFile] [unPlaced.fa] [unPlaced.fa]
show-coords -H -T -l -c -o [deltaFile] > [coordsFile]
*** Additional analysis was performed on the alignments to find and remove contigs contained within two contigs
with the ends overlapping (see methods) ***
----------------------
Further screening
----------------------
kraken --db [database] [APG_Sequences.fa]
blastn -db [nt] -query [kraken_nonMamalHits.fa] -outfmt "6 qseqid sseqid pident length
mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qcovs qcovhsp
staxids sscinames" -max_hsps 1 -max_target_seqs 1 –out [blastOutput]
bwa index [GRCh38.p10_primaryChrs]
bwa index [GRCh38.p10]
bwa mem [GRCh38.p10_primaryChrs] [APG_Sequences_noContamiants.fa]
bwa mem [GRCh38.p10] [APG_Sequences_noContamiants.fa]
```

First, let's prepare a script for the scraping. We'll then assemble the reads per individual as indicated in the manuscript supplementary information.

> Ceres: /home/derek.bickharhth/bostauruscnv/novel_seq

```bash
# Getting the list of bam IDs
ls ../bam/*.bam | cut -d'/' -f3 | cut -d'.' -f1 > bam_id_nums.list

# Scraping the unmapped reads
for i in `cat bam_id_nums.list`; do echo $i; sbatch scrap_unmapped_reads.sh ../bam $i scrap; done
```