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
module load repeatmasker/4.0.7
sbatch --nodes=1 --mem=30000 --ntasks-per-node=10 -p medium --wrap="RepeatMasker -pa 10 -q -species cow -no_is -gff ../ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna"
# Building species libraries in: /home/derek.bickharhth/.RepeatMaskerCache/dc20170127/bos_taurus
#   - 175 ancestral and ubiquitous sequence(s) for bos taurus   <- ruh-roh! That ain't good!
#   - 0 lineage specific sequence(s) for bos taurus

# OK, the process terminated without finding repeats! I think there are two reasons why:
# The Scinet team have an outdated repeat library installed
# Also, the fasta I used was soft-masked by NCBI

# First, let's get the repeat library updated
wget http://www.repeatmasker.org/libraries/RepeatMaskerMetaData-20181026.tar.gz
# Hmm, it's not so simple to install this as a local library. I could try hacking @INC to point to this directory, but it's going to screw up a bunch of things
# Let's try the run again but with the proper assembly fasta this time.

sbatch --nodes=1 --mem=35000 --ntasks-per-node=20 -p medium --wrap="RepeatMasker -pa 20 -q -species cow -no_is -gff ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa"

# The repeatmasker module was improperly installed! No wonder why I was getting nothing!
sbatch --nodes=1 --mem=45000 --ntasks-per-node=50 -p msn --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 50 -q -species cow -no_is -gff ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa"

# OK, now I need the UMD3 reference genome
for i in `seq 1 29` X Y U; do echo $i; wget /pub/data/assembly/Bos_taurus/Bos_taurus_UMD_3.1/Chr${i}.fa.gz; done
for i in `seq 1 29` X Y U; do echo $i; gunzip Chr${i}.fa.gz; cat Chr${i}.fa >> umd3_reference_genome.fasta; done

module load bwa samtools
module load java/64/1.8.0_121

sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index umd3_reference_genome.fasta"

# Identifying gaps
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p msn --wrap="perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o umd3_reference_genome.fasta -s ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa -g ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -j /software/7/apps/java/1.8.0_121/bin/java -d umd3_gaps_on_arsucd.tab"

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f umd3_gaps_on_arsucd.tab -c 0 -m -d '\t'
|Entry    | Value|
|:--------|-----:|
|Closed   | 58110|
|Trans    |  8508|
|Unmapped |  5570|
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

# Checking the output from the slurm STDOUT reports:
echo -e "Animal\tUnmapReads\tR1UnmapReads\tR2UnmapReads"; for i in slurm-*.out; do perl -e '@d; $n; while(<>){chomp; if($_ =~ /^mkdir/){next;}elsif($_ =~ /^\[/){@b = split(/\s+/); push(@d, $b[2]);}elsif($_ =~ /^Done/){@b = split(/\s+/); $n = $b[3];}} print "$n\t$d[1]\t$d[3]\t$d[5]\n";' < $i; done

# Crap! there are no completely unmapped read pairs! I think it's part of Bob's processing pipeline

# Let's calculate the number of links per chromosome to see how often the one-end anchors could be expected.
# Setting up 10 mb windows across the genome.
perl -lane 'for($x = 1; $x < $F[1]; $x += 10000000){$e = $x + 10000000; if($x + 10000000 > $F[1]){$e = $F[1];} print "$F[0]\t$x\t$e";}' < ../assembly/ARS-UCD1.2_Btau5.0.1Y.fa.fai > ARSUCD1.2.links.bed

module load bedtools; for i in scrap/*; do name=`echo $i | cut -d'/' -f2`; echo $name; perl -lane '$e = $F[3] + 1; print "$F[2]\t$F[3]\t$e";' < $i/$name.links.bam | bedtools intersect -a ARSUCD1.2.links.bed -b stdin -c > $name.wins.intersect.bed; done

# I'm going to go out on a limb and try to intersect these
python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -c 0 -c 1 -c 2 -m 3 -o left_join_total.wins.intersect.bed -f 1289.wins.intersect.bed -f 1315.wins.intersect.bed -f 1510.wins.intersect.bed -f 18394.wins.intersect.bed -f 18395.wins.intersect.bed -f 18474.wins.intersect.bed -f 18476.wins.intersect.bed -f 18557.wins.intersect.bed -f 185615.wins.intersect.bed -f 185616.wins.intersect.bed -f 185617.wins.intersect.bed -f 185618.wins.intersect.bed -f 185619.wins.intersect.bed -f 185620.wins.intersect.bed -f 185621.wins.intersect.bed -f 185622.wins.intersect.bed -f 185623.wins.intersect.bed -f 185624.wins.intersect.bed -f 185625.wins.intersect.bed -f 185626.wins.intersect.bed -f 185627.wins.intersect.bed -f 185628.wins.intersect.bed -f 185629.wins.intersect.bed -f 185630.wins.intersect.bed -f 185631.wins.intersect.bed -f 185632.wins.intersect.bed -f 185633.wins.intersect.bed -f 185634.wins.intersect.bed -f 185635.wins.intersect.bed -f 185636.wins.intersect.bed -f 185637.wins.intersect.bed -f 185639.wins.intersect.bed -f 185641.wins.intersect.bed -f 185643.wins.intersect.bed -f 185644.wins.intersect.bed -f 185646.wins.intersect.bed -f 185647.wins.intersect.bed -f 185648.wins.intersect.bed -f 185649.wins.intersect.bed -f 185650.wins.intersect.bed -f 185651.wins.intersect.bed -f 185652.wins.intersect.bed -f 185653.wins.intersect.bed -f 185654.wins.intersect.bed -f 185655.wins.intersect.bed -f 185656.wins.intersect.bed -f 185657.wins.intersect.bed -f 185658.wins.intersect.bed -f 185659.wins.intersect.bed -f 185660.wins.intersect.bed -f 185661.wins.intersect.bed -f 185662.wins.intersect.bed -f 185663.wins.intersect.bed -f 185664.wins.intersect.bed -f 185665.wins.intersect.bed -f 185666.wins.intersect.bed -f 185667.wins.intersect.bed -f 185668.wins.intersect.bed -f 185669.wins.intersect.bed -f 185670.wins.intersect.bed -f 185671.wins.intersect.bed -f 185672.wins.intersect.bed -f 185673.wins.intersect.bed -f 185674.wins.intersect.bed -f 185675.wins.intersect.bed -f 185676.wins.intersect.bed -f 185677.wins.intersect.bed -f 185678.wins.intersect.bed -f 185679.wins.intersect.bed -f 185680.wins.intersect.bed -f 185681.wins.intersect.bed -f 185682.wins.intersect.bed -f 185683.wins.intersect.bed -f 185684.wins.intersect.bed -f 185685.wins.intersect.bed -f 185686.wins.intersect.bed -f 185687.wins.intersect.bed -f 185688.wins.intersect.bed -f 185689.wins.intersect.bed -f 185690.wins.intersect.bed -f 185691.wins.intersect.bed -f 185692.wins.intersect.bed -f 185693.wins.intersect.bed -f 185694.wins.intersect.bed -f 185695.wins.intersect.bed -f 185696.wins.intersect.bed -f 185697.wins.intersect.bed -f 185698.wins.intersect.bed -f 185699.wins.intersect.bed -f 185700.wins.intersect.bed -f 185701.wins.intersect.bed -f 185702.wins.intersect.bed -f 185703.wins.intersect.bed -f 185704.wins.intersect.bed -f 18682.wins.intersect.bed -f 19020.wins.intersect.bed -f 19225.wins.intersect.bed -f 193865.wins.intersect.bed -f 193866.wins.intersect.bed -f 193867.wins.intersect.bed -f 193868.wins.intersect.bed -f 193869.wins.intersect.bed -f 193870.wins.intersect.bed -f 193871.wins.intersect.bed -f 193872.wins.intersect.bed -f 193873.wins.intersect.bed -f 193874.wins.intersect.bed -f 193875.wins.intersect.bed -f 193876.wins.intersect.bed -f 193877.wins.intersect.bed -f 193878.wins.intersect.bed -f 193879.wins.intersect.bed -f 193880.wins.intersect.bed -f 193881.wins.intersect.bed -f 193882.wins.intersect.bed -f 193883.wins.intersect.bed -f 193884.wins.intersect.bed -f 193885.wins.intersect.bed -f 193886.wins.intersect.bed -f 193887.wins.intersect.bed -f 193888.wins.intersect.bed -f 193889.wins.intersect.bed -f 193890.wins.intersect.bed -f 193891.wins.intersect.bed -f 193892.wins.intersect.bed -f 193893.wins.intersect.bed -f 193894.wins.intersect.bed -f 193895.wins.intersect.bed -f 193896.wins.intersect.bed -f 193897.wins.intersect.bed -f 193898.wins.intersect.bed -f 193899.wins.intersect.bed -f 193900.wins.intersect.bed -f 193901.wins.intersect.bed -f 193902.wins.intersect.bed -f 193903.wins.intersect.bed -f 193904.wins.intersect.bed -f 193905.wins.intersect.bed -f 193906.wins.intersect.bed -f 193907.wins.intersect.bed -f 193908.wins.intersect.bed -f 193909.wins.intersect.bed -f 193910.wins.intersect.bed -f 193911.wins.intersect.bed -f 193912.wins.intersect.bed -f 193913.wins.intersect.bed -f 193914.wins.intersect.bed -f 193915.wins.intersect.bed -f 193916.wins.intersect.bed -f 193917.wins.intersect.bed -f 193918.wins.intersect.bed -f 193919.wins.intersect.bed -f 193920.wins.intersect.bed -f 193921.wins.intersect.bed -f 193922.wins.intersect.bed -f 193923.wins.intersect.bed -f 193924.wins.intersect.bed -f 193925.wins.intersect.bed -f 193926.wins.intersect.bed -f 193927.wins.intersect.bed -f 193928.wins.intersect.bed -f 193929.wins.intersect.bed -f 193930.wins.intersect.bed -f 193931.wins.intersect.bed -f 193932.wins.intersect.bed -f 193933.wins.intersect.bed -f 193934.wins.intersect.bed -f 193935.wins.intersect.bed -f 193936.wins.intersect.bed -f 193937.wins.intersect.bed -f 193938.wins.intersect.bed -f 193939.wins.intersect.bed -f 193940.wins.intersect.bed -f 193941.wins.intersect.bed -f 193942.wins.intersect.bed -f 193943.wins.intersect.bed -f 193944.wins.intersect.bed -f 193945.wins.intersect.bed -f 193946.wins.intersect.bed -f 193947.wins.intersect.bed -f 193948.wins.intersect.bed -f 193949.wins.intersect.bed -f 193950.wins.intersect.bed -f 193951.wins.intersect.bed -f 193952.wins.intersect.bed -f 193953.wins.intersect.bed -f 193954.wins.intersect.bed -f 193955.wins.intersect.bed -f 193956.wins.intersect.bed -f 193957.wins.intersect.bed -f 193958.wins.intersect.bed -f 193959.wins.intersect.bed -f 193960.wins.intersect.bed -f 193961.wins.intersect.bed -f 193962.wins.intersect.bed -f 193963.wins.intersect.bed -f 193964.wins.intersect.bed -f 193965.wins.intersect.bed -f 193966.wins.intersect.bed -f 193967.wins.intersect.bed -f 193968.wins.intersect.bed -f 193969.wins.intersect.bed -f 193970.wins.intersect.bed -f 193971.wins.intersect.bed -f 193972.wins.intersect.bed -f 193973.wins.intersect.bed -f 193974.wins.intersect.bed -f 193975.wins.intersect.bed -f 193976.wins.intersect.bed -f 193977.wins.intersect.bed -f 193978.wins.intersect.bed -f 193979.wins.intersect.bed -f 193980.wins.intersect.bed -f 193981.wins.intersect.bed -f 193982.wins.intersect.bed -f 193983.wins.intersect.bed -f 193984.wins.intersect.bed -f 193985.wins.intersect.bed -f 193986.wins.intersect.bed -f 193987.wins.intersect.bed -f 193988.wins.intersect.bed -f 193989.wins.intersect.bed -f 193990.wins.intersect.bed -f 193991.wins.intersect.bed -f 193992.wins.intersect.bed -f 193993.wins.intersect.bed -f 193994.wins.intersect.bed -f 193995.wins.intersect.bed -f 193996.wins.intersect.bed -f 193997.wins.intersect.bed -f 193998.wins.intersect.bed -f 193999.wins.intersect.bed -f 194000.wins.intersect.bed -f 194001.wins.intersect.bed -f 194002.wins.intersect.bed -f 19423.wins.intersect.bed -f 19599.wins.intersect.bed -f 19628.wins.intersect.bed -f 2089.wins.intersect.bed -f 2133.wins.intersect.bed -f 22009.wins.intersect.bed -f 22049.wins.intersect.bed -f 22111.wins.intersect.bed -f 22168.wins.intersect.bed -f 22224.wins.intersect.bed -f 22448.wins.intersect.bed -f 22640.wins.intersect.bed -f 22641.wins.intersect.bed -f 22726.wins.intersect.bed -f 22856.wins.intersect.bed -f 23085.wins.intersect.bed -f 23150.wins.intersect.bed -f 23307.wins.intersect.bed -f 23574.wins.intersect.bed -f 23737.wins.intersect.bed -f 24004.wins.intersect.bed -f 24031.wins.intersect.bed -f 24045.wins.intersect.bed -f 24612.wins.intersect.bed -f 27164.wins.intersect.bed -f 27249.wins.intersect.bed -f 35808.wins.intersect.bed -f 47325.wins.intersect.bed -f 47332.wins.intersect.bed -f 47381.wins.intersect.bed -f 6606.wins.intersect.bed -f 6980.wins.intersect.bed -f 88651.wins.intersect.bed -f 88653.wins.intersect.bed -f 88655.wins.intersect.bed -f 986.wins.intersect.bed

# Wow! It worked the first try!!!! I just need to make a header
echo -ne "chr\tstart\tend" > left_join_header; for i in *.intersect.bed; do name=`echo $i | cut -d'.' -f1`; echo -ne "\t$name"; done >> left_join_header; echo >> left_join_header
cat left_join_header left_join_total.wins.intersect.bed > left_join_header.wins.intersect.bed

mkdir oea_link_counts
mv *.intersect.bed ./oea_link_counts/
```

Amateur attempts to try to find interesting outliers by scaling:

```R
library(dplyr)
links <- read.delim("oea_link_counts/left_join_header.wins.intersect.bed", header=TRUE)
links.trim <- head(links, n=267)

# Outliers by column
links.trim %>% filter_all(any_vars(. >= quantile(., 0.99, na.rm = TRUE)))

# The scale function has several attributes that I can grab. I'm going to take the scale factor here from the outliers
scale_factors <- attributes(scale(links.trim %>% filter_all(any_vars(. >= quantile(., 0.99, na.rm = TRUE)))))$`scaled:scale`

# Visual inspection of the top 5%
scale_factors[scale_factors >= quantile(scale_factors, 0.95)]
```