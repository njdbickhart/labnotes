# Generating the red clover reference assembly
---
*12/14/2018*

These are my notes on generating the data, performing the assembly and providing some analysis for the red clover genome assembly project.

## Table of Contents


## Diagnostic assembly

I want to try to assemble a portion of the data (~ 15 Gbp) of the first batch of red clover nanopore reads for diagnosis of how much more data we need. I am going to queue up a full job and hopefully it will run to conclusion over the weekend.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
module load canu/1.8 java/64/1.8.0_121
# Run with all nanopore reads, with parameters designed to reduce influence of systematic error overestimation
sbatch --nodes=1 --ntasks-per-node=20 --mem=64G -p medium --wrap="canu -p clover_hen_test -d clover_hen_test genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14'  -nanopore-raw ./*/*.fastq"

# Note, the meryl routine had substantial problems with the number of filtered short reads in our dataset! 
# Some datasets had >90% filtered reads because of reads < 1000 bp. 
# It appears to be running still, since we're still above 30 X coverage for the genome size.

# OK, so Canu submitted everything on the "short" queue, and one job went over the two day limit and was lost
```

I just uploaded the rest of the clover data and now I'm going to try the whole assembly.

```bash
module load canu/1.8 java/64/1.8.0_121
sbatch --nodes=1 --ntasks-per-node=20 --mem=64G -p medium --wrap="canu -p clover_hen_full -d clover_hen_full genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p medium' -nanopore-raw ./*/*.fastq"
```

OK, the output of the unitiging program appears to show a much lower coverage than expected. Only 36% of reads as "unique" and 13 X coverage at that! I'm going to let the current run finish, but Serge told me that there's another method to try to fix the input data: manual correction prior to another round of manual correction.

His suggestion is to run with "-correction corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100" and then use the output reads in the full pipeline run afterwards (same settings).

He didn't see any signs of heterozygosity or other issues. 

Yeah, the output just produced an assembly with an NG50 contig size of 51 kbp and a full assembly length of 392,742,549 and 10,975 contigs. I'm going to try Serge's approach first.

```bash
# First, the correction
sbatch --nodes=1 --ntasks-per-node=20 --mem=64G -p medium --wrap="canu -correct -p clover_hen_retry -d clover_hen_retry genomeSize=420m corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 'gridOptions=-p medium' -nanopore-raw ./*/*.fastq"

sbatch --nodes=1 --ntasks-per-node=20 --mem=64G -p medium --wrap="canu -p clover_hen_postcor -d clover_hen_postcor genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p medium' -nanopore-raw clover_hen_retry/clover_hen_retry.correctedReads.fasta.gz"
```


## KAT analysis of Diagnostic assembly

I am using the "clover_hen_postcor" assembly version in this analysis, as it had the most assembled bases and the largest average contig lengths. 

I will be downloading several SRA datasets for future analysis, but I will be using the one with the highest mapping rate for kmer-based analysis of assembly completion. If the assembly is fairly complete, then we just need more coverage to increase the contiguity. If we're missing big chunks, then there's a heterozygosity problem or a bias in the sequencing. 

Here are the SRA accessions I will be downloading (all Aberystwyth University, UK):
* ERR3063534
* ERR3063535
* ERR3063536
* ERR3063537
* ERR3063538
* ERR3063539
* ERR3063540
* ERR3063541

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/red_clover_public

```bash
sbatch ~/bin/download_sra_dataset.sh clover_uk_sra_list.txt
```

Preparing assembly fasta for alignment.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_hen_postcor

```bash
module load bwa
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p short --wrap="bwa index clover_hen_postcor.contigs.fasta"
```

Now to align the fastqs to the assembly, separately:

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/red_clover_public

```bash
module load samtools
ls *.fastq.gz | xargs -I {} sbatch --nodes=1 --ntasks-per-node=10 --mem=25000 -p short --wrap="bwa mem -t 8 -M ../clover_hen_postcor/clover_hen_postcor.contigs.fasta {} | samtools sort -m 2G -o {}.sorted.bam -T {}.temp -"

ls *.bam | xargs -I {} sbatch --nodes=1 --ntasks-per-node=1 --mem=10000 -p short --wrap="samtools index {}"

# Checking the mapping percentages:
for i in *.sorted.bam; do echo -ne "$i\t"; samtools idxstats $i | perl -e '$un = 0; $m = 0; while(<STDIN>){chomp; @s = split(/\t/); $m += $s[2]; $un += $s[3];} $perc = $m / ($m + $un); $perc *= 100; print "$m\t$un\t$perc\n";' ; done
ERR3063534_1.fastq.gz.sorted.bam        159503126       25574908        86.1815541005801
ERR3063535_1.fastq.gz.sorted.bam        199324740       26029001        88.4497142650053
ERR3063536_1.fastq.gz.sorted.bam        169998247       31918785        84.1921284777997
ERR3063537_1.fastq.gz.sorted.bam        190936056       32079594        85.6155413308438
ERR3063538_1.fastq.gz.sorted.bam        198148738       32623653        85.863277292993
ERR3063539_1.fastq.gz.sorted.bam        189012804       28185700        87.0230689986705
ERR3063540_1.fastq.gz.sorted.bam        201118225       29134382        87.3467743190417
ERR3063541_1.fastq.gz.sorted.bam        197161358       37262635        84.1045984572066

# Now checking for contig coverage using bedtools
for i in red_clover_public/*.sorted.bam; do echo $i; sbatch --nodes=1 --ntasks-per-node=2 --mem=16000 -p short --wrap="bedtools genomecov -ibam $i -g clover_hen_postcor/clover_hen_postcor.contigs.fasta > $i.cov"; done

# I'm noticing a high degree of zero coverage bases on contigs. Like so:
for i in red_clover_public/*.cov; do echo -ne "$i\t"; avg=`perl -lane 'if($F[1] == 0){print $F[4];}' < $i | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl | grep 'Average'`; echo $avg; done
red_clover_public/ERR3063534_1.fastq.gz.sorted.bam.cov  Average 0.854841
red_clover_public/ERR3063535_1.fastq.gz.sorted.bam.cov  Average 0.835845
red_clover_public/ERR3063536_1.fastq.gz.sorted.bam.cov  Average 0.840202
red_clover_public/ERR3063537_1.fastq.gz.sorted.bam.cov  Average 0.849392
red_clover_public/ERR3063538_1.fastq.gz.sorted.bam.cov  Average 0.892309
red_clover_public/ERR3063539_1.fastq.gz.sorted.bam.cov  Average 0.853269
red_clover_public/ERR3063540_1.fastq.gz.sorted.bam.cov  Average 0.882853
red_clover_public/ERR3063541_1.fastq.gz.sorted.bam.cov  Average 0.853424

# I think that the sequence data may be really bonkers...
for i in red_clover_public/*.cov; do echo -ne "$i\t"; max=`cat $i | cut -f2 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl | grep 'Maximum'`; echo $max; done
red_clover_public/ERR3063534_1.fastq.gz.sorted.bam.cov  Maximum 1253838
red_clover_public/ERR3063535_1.fastq.gz.sorted.bam.cov  Maximum 1498571
red_clover_public/ERR3063536_1.fastq.gz.sorted.bam.cov  Maximum 1294225
red_clover_public/ERR3063537_1.fastq.gz.sorted.bam.cov  Maximum 2098104
red_clover_public/ERR3063538_1.fastq.gz.sorted.bam.cov  Maximum 1704903
red_clover_public/ERR3063539_1.fastq.gz.sorted.bam.cov  Maximum 1787357
red_clover_public/ERR3063540_1.fastq.gz.sorted.bam.cov  Maximum 1525081
red_clover_public/ERR3063541_1.fastq.gz.sorted.bam.cov  Maximum 1365049

# I'm going to map the short-read ENA red clover assembly (2015) onto the new assembly to see just how different things are here.
sbatch --nodes=1 --ntasks-per-node=3 --mem=20000 -p short --wrap="minimap2 -x asm10 clover_hen_postcor/clover_hen_postcor.contigs.fasta short_read_red_clover_ena.fasta > clover_hen_postcor.vs.short_read.paf"
sbatch --nodes=1 --ntasks-per-node=3 --mem=20000 -p short --wrap="minimap2 -x asm10 short_read_red_clover_ena.fasta clover_hen_postcor/clover_hen_postcor.contigs.fasta > short_read.vs.clover_hen_postcor.paf"

perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); if($s[11] > 0){ $c += $s[9];}} print "$c\n";' < clover_hen_postcor.vs.short_read.paf
88,982,406
perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); if($s[11] > 0){ $c += $s[9];}} print "$c\n";' < short_read.vs.clover_hen_postcor.paf
156,057,180
```

#### Aligning SSR markers to the Red clover assembly

I just want to see how many of Heathcliffe/John's SSR markers align to the current iteration of red clover. Let's see...

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
dos2unix ssr_marker_est.fasta

# Don't try this at home kids! I queued up an interactive shell prior to running this short command
module load bwa
bwa mem clover_hen_postcor/clover_hen_postcor.contigs.fasta ssr_marker_est.fasta > ssr_marker_est.clover_hen_postcor.sam

# Now to turn this into a human-readable format using one of my scripts
echo -e "Marker\tSamFlag\tContig\tStart\tAlignEnd\tAlignLen\tMappingQual" > ssr_marker_est.clover_hen_postcor.tab; perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_scripts/BriefSamOutFormat.pl -s ssr_marker_est.clover_hen_postcor.sam >> ssr_marker_est.clover_hen_postcor.tab
```

## Second-round draft assembly

OK, we've generated about 30 Gigabases of more data and I think we're ready to try to polish the assembly again. I'll try two assembly approaches: 1. the normal Canu run and 2. the post-correct Canu run.

#### Unpacking the data and modifying it

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
tar -xvf Clover.tar.gz
sbatch --nodes=1 --ntasks-per-node=1 --mem=2000 -p msn --wrap="tar -xvf Clover2.tar.gz"
sbatch --nodes=1 --mem=2000 --ntasks-per-node=1 -p msn --wrap="tar -xvf CloverB.tar.gz"

mkdir clover12_fastqs
mv Clover/20190508_2255_GA10000_FAK10494_8301799c/fastq_pass/*.fastq ./clover12_fastqs/
sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="mkdir clover13_fastqs; mv Clover2/*/fastq_pass/*.fastq ./clover13_fastqs/"
sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="mkdir clover14_fastqs; mv CloverB/*/fastq_pass/*.fastq ./clover14_fastqs/"

# Mash sketching
for i in `seq 1 14`; do name=clover${i}_fastqs; echo $name; cat clover${i}_fastqs/*.fastq > clover${i}_fastqs/clover${i}.combined.fastq; done
for i in `seq 1 14`; do name=clover${i}_fastqs; echo $name; sbatch --nodes=1 --mem=10000 --ntasks-per-node=4 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -o clover${i}.msh -p 4 -s 100000 -r -m 4 -g 420M $name/clover${i}.combined.fastq"; done
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste clover_combined clover10.msh clover11.msh clover12.msh clover13.msh clover14.msh clover1.msh clover2.msh clover3.msh clover4.msh clover5.msh clover6.msh clover7.msh clover8.msh clover9.msh
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -t clover_combined.msh clover10.msh clover11.msh clover12.msh clover13.msh clover14.msh clover1.msh clover2.msh clover3.msh clover4.msh clover5.msh clover6.msh clover7.msh clover8.msh clover9.msh > clover_combined.dist

perl -ne '$_ =~ s/(clover\d+)_fastqs\/clover\d+\.combined\.fastq/$1/g; print $_;' < clover_combined.dist > clover_combined.rfmt.dist

for i in 12 13 14; do perl -e 'chomp(@ARGV); open(IN, $ARGV[0]); while(<IN>){$s = <IN>; chomp($s); print "$ARGV[1]\t" . length($s) . "\n"; <IN>; <IN>;} close IN;' clover${i}_fastqs/clover${i}.combined.fastq $i; done > clover_read_lengths.new.tab

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f clover_read_lengths.new.tab -c 0 -d '\t' -m

perl -e '$h = <>; @s = split(/\s+/, $h); print join("\t", @s) . "\n"; %d; while($t = <>){chomp($t); @s = split(/\s+/, $t); shift(@s); $d{$s[0]} = [$s[1], $s[2]];} foreach my $f (sort {$a <=> $b} keys(%d)){print "$f\t" . join("\t", @{$d{$f}}) . "\n";}' < read_count_summary.tab > read_count_summary.sorted.tab
```

|Entry |   Value|
|:-----|-------:|
|14    | 1007247|
|13    |  781812|
|12    |  746246|

Now to plot the summary statistics to share with collaborators.

```R
library(dplyr)
library(ggplot2)
library(reshape)
dist <- read.delim("clover_combined.rfmt.dist", header=TRUE)
rownames(dist) <- dist$X.query
dist <- dist[, c(2:15)]
dist.m <- melt(as.matrix(dist))

pdf(file="clover_dataset_distance_heatmap.pdf", useDingbats=FALSE)
ggplot(dist.m, aes(X1, X2)) + geom_tile(aes(fill=value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") + labs(title = "Clover Dataset Mash Distances")
dev.off()

reads <- read.delim("clover_read_lengths.tab", header=FALSE)
colnames(reads) <- c("Flowcell", "ReadLen")

library(gridExtra)
library(ggridges)
library(viridis)
pdf("clover_read_distribution.pdf")
ggplot(reads, aes(x=ReadLen, y=Flowcell, fill=..x..)) + geom_density_ridges_gradient(scale=3, rel_min_height = 0.01) + scale_fill_viridis(name = "Read Length (bp)", option= "C") + labs(title= "Read Length Density per Flowcell") + xlim(0, 50000)
dev.off()

cbar <- ggplot(read.counts, aes(x=Flowcell, y=count)) + geom_bar(stat="identity", color="blue", fill="white") + labs(title = "A. Read counts per flowcell", ylab="Number of Reads")
bbar <- ggplot(read.counts, aes(x=Flowcell, y=bases)) + geom_bar(stat="identity", color="red", fill="white") + labs(title = "B. Base pair counts per flowcell", ylab="Sum of DNA bases (bp)", xlab = "Flowcells (in numerical order)")
pdf("read_length_bases_barchart.pdf")
grid.arrange(cbar, bbar, nrow=2)
dev.off()

new <- read.delim("clover_read_lengths.new.tab", header=FALSE)
colnames(new) <- c("Flowcell", "ReadLen")

reads <- reads %>% mutate(Tag = "Old")
new <- new %>% mutate(Tag = "New")
new$Flowcell <- as.factor(new$Flowcell)


rcombine <- bind_rows(reads, new)
rcombine$Flowcell <- as.factor(rcombine$Flowcell)
rcombine$Tag <- as.factor(rcombine$Tag)

pdf("clover_read_dist_new.pdf", useDingbats=FALSE)
ggplot(rcombine, aes(x=ReadLen, y=Flowcell, fill=..x..)) + geom_density_ridges_gradient(scale=3, rel_min_height = 0.01) + scale_fill_viridis(name = "Read Length (bp)", option= "C") + scale_x_log10() + labs(title= "Read Length Density per Flowcell")
dev.off()

read.counts <- group_by(reads, Flowcell) %>% summarize(count = n(), bases = sum(as.numeric(ReadLen)))
rcombine.counts <- group_by(reads, Flowcell) %>% summarize(count = n(), bases = sum(as.numeric(ReadLen)))

sorted <- read.delim("read_count_summary.sorted.tab", header=TRUE)
sorted$Tag <- c(rep("Old", 11), rep("New", 3))
sorted$Flowcell <- as.factor(sorted$Flowcell)
sorted$Tag <- as.factor(sorted$Tag)

cbar <- ggplot(sorted, aes(x=Flowcell, y=count, fill=Tag)) + geom_bar(stat="identity", color="blue") + labs(title = "A. Read counts per flowcell", ylab="Number of Reads") + scale_fill_manual(values=c("blue", "white"))
bbar <- ggplot(sorted, aes(x=Flowcell, y=bases, fill=Tag)) + geom_bar(stat="identity", color="red") + labs(title = "B. Base pair counts per flowcell", ylab="Sum of DNA bases (bp)", xlab = "Flowcells (in numerical order)") + scale_fill_manual(values=c("red", "white"))
pdf("clover_new_read_distribution.pdf", useDingbats=FALSE)

```

#### The assembly

```bash
module load canu/1.8 java/1.8.0_121

### De novo, no-pre correct
sbatch --nodes=1 --ntasks-per-node=30 --mem=75G -p msn --wrap="canu -p clover_no_pre -d clover_no_pre genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p msn' -nanopore-raw ./*/*.fastq"

### De novo, pre-correct
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short --wrap="canu -correct -p clover_pre -d clover_pre genomeSize=420m corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 saveReadCorrections=true 'gridOptions=-p msn' -nanopore-raw ./*/*combined.fastq"

sbatch --nodes=1 --ntasks-per-node=30 --mem=10G -p msn --wrap="canu -p clover_correct -d clover_correct genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p msn,short' -nanopore-raw clover_pre/clover_pre.correctedReads.fasta.gz"
```


## Nanopolish

I was able to install nanopolish in my home directory on Ceres. It looks like it requires allot of information that wasn't saved from our previous runs (or was a pain to save!), so I will limit my polishing to the three new cells (~30 gigabases) we just produced for Clover. 

```bash
module load hdf5/1.8.19 eigen/3.2.6

```

## Checking the validity of larger assembly sizes. 

So, our two clover assemblies finished, but they're both predicting an assembly size of 700 megabases! That isn't right... Let's test a couple of hypotheses here. 

#### Potential sources for larger assemblies
1. Contamination
	* I can test this by using centrifuge on each contig
	* If any contig is closer to human, then that's a potential source of contamination. 
	* KAT is also a natural comparison here.
2. Bias in the libraries
	* It could be that the earlier data is causing issues due to the different basecalling
	* What if we left it out of the assembly?
	* To test, let's see how many contigs are only supported by the earlier data through read alignments
	* Also, low coverage contigs would be interesting to see
3. Clover hybrid
	* This is more difficult, but could be a possibility? 
	* If all else fails to show significant differences, this could be a cryptic hybrid
	* How to test?


Let's queue up tests for the first two possibilities in the meantime. It should be easy enough to find anomalies for both.

> Ceres:  /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
module load bwa samtools minimap2 jellyfish2/2.2.9

### The read alignment test ###
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index clover_correct/clover_correct.contigs.fasta"
for i in red_clover_public/*.fastq.gz; do echo $i; sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p msn --dependency=afterok:662826 --wrap="bwa mem clover_correct/clover_correct.contigs.fasta $i | samtools sort -T ${i}.tmp -o ${i}.correct.sort.bam -"; done

# Now to align the "new" and "old" reads to the assembly separately
# I think that I can get away with just aligning the "new" and removing contigs that don't have any coverage
sbatch --nodes=1 --mem=300 --ntasks-per-node=1 -p msn --wrap="cat ./clover12_fastqs/clover12.combined.fastq ./clover13_fastqs/clover13.combined.fastq ./clover14_fastqs/clover14.combined.fastq > clover_new_combined_fastqs.fastq"

sbatch --nodes=1 --mem=15000 --ntasks-per-node=4 -p msn --wrap="minimap2 -a -x map-ont clover_correct/clover_correct.contigs.fasta clover_new_combined_fastqs.fastq | samtools sort -T clov_combined.temp -o clover_new_combined_fastqs.sort.bam -"


### The kmer composition test ###
# Let's start by making Jellyfish hashes from everything in existence
for i in clover_new_combined_fastqs.fastq clover_correct/clover_correct.contigs.fasta short_read_red_clover_ena.fasta; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=5 --wrap="jellyfish count -m 21 -s 90G -t 5 -o ${i}.jf -C $i"; done
```