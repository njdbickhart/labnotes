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

# Gathering the information from canu on each contig
grep '>' clover_correct/clover_correct.contigs.fasta | perl -e 'print "Contig\tLen\tReads\tCovStat\tgappedBases\tclass\tsuggestRepeat\tsuggestCircular\n"; while(<STDIN>){chomp; @data = (); @segs = split(/\s+/); $ctg = shift(@segs); $ctg =~ s/\>//; push(@data, $ctg); foreach my $r (@segs){($t) = $r =~ m/.+\=(.+)/; push(@data, $t . "");} print join("\t", @data); print "\n";}' > clover_correct/clover_correct.contig_header.info.tab

# Checking to see if any are circular or repetitive
perl -lane 'if($F[-1] eq "yes" || $F[-2] eq "yes"){print $_;}' < clover_correct/clover_correct.contig_header.info.tab | wc -l
483

# Now, just circular
perl -lane 'if($F[-1] eq "yes"){print $_;}' < clover_correct/clover_correct.contig_header.info.tab | wc -l
3

# Interesting... let's see what they are
perl -lane 'if($F[-1] eq "yes"){print $_;}' < clover_correct/clover_correct.contig_header.info.tab
tig00002968     97363   84      22.22   no      contig  no      yes	<- hits to 18S rRNA and to medicago chloroplast
tig00004516     108709  35      60.66   no      contig  yes     yes <- medicago chr3, receptor-like protein genes
tig00072474     98236   22      69.68   no      contig  no      yes <- no significant similarity results

# I wonder what the nr database thinks of these? Ran each one on nucleotide blast with the megablast option
samtools faidx clover_correct/clover_correct.contigs.fasta
```

So, we've got a partial chloroplast (should be 200 kb?) and one potential virus? Maybe?


```bash
### The read alignment test ###
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index clover_correct/clover_correct.contigs.fasta"
for i in red_clover_public/*.fastq.gz; do echo $i; sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p msn --dependency=afterok:662826 --wrap="bwa mem clover_correct/clover_correct.contigs.fasta $i | samtools sort -T ${i}.tmp -o ${i}.correct.sort.bam -"; done

# Now to align the "new" and "old" reads to the assembly separately
# I think that I can get away with just aligning the "new" and removing contigs that don't have any coverage
sbatch --nodes=1 --mem=300 --ntasks-per-node=1 -p msn --wrap="cat ./clover12_fastqs/clover12.combined.fastq ./clover13_fastqs/clover13.combined.fastq ./clover14_fastqs/clover14.combined.fastq > clover_new_combined_fastqs.fastq"

sbatch --nodes=1 --mem=15000 --ntasks-per-node=4 -p msn --wrap="minimap2 -a -x map-ont clover_correct/clover_correct.contigs.fasta clover_new_combined_fastqs.fastq | samtools sort -T clov_combined.temp -o clover_new_combined_fastqs.sort.bam -"

sbatch --nodes=1 --mem=2000 --ntasks-per-node=1 -p msn --wrap="samtools index clover_new_combined_fastqs.sort.bam"

# Identify contigs with zero mappings
samtools idxstats clover_new_combined_fastqs.sort.bam | perl -lane 'if($F[2] == 0){print $_;}' | head
*       0       0       36138  <- Only the unaligned reads

# Now, get a distribution of the count of reads to contig length (not very indicative, but an approximation)
samtools idxstats clover_new_combined_fastqs.sort.bam | perl -lane 'if($F[1] != 0){print ($F[2] / $F[1]);}' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3132
Sum:    142.036951768705
Minimum 0.000349545590732048
Maximum 27.8661461238148
Average 0.045350
Median  0.00455151813437383
Standard Deviation      0.726313
Mode(Highest Distributed Value) 0.00455793719598611

# Number of contigs with read/contig len ratios greater than the average
samtools idxstats clover_new_combined_fastqs.sort.bam | perl -lane 'if($F[1] != 0){print ($F[2] / $F[1]);}' | perl -lane 'if($F[0] > 0.45){print $_;}' | wc -l
28

# Let's align the old reads and do the comparison next
sbatch --nodes=1 --mem=800 --ntasks-per-node=2 -t 1-0 -p msn --wrap="for i in `seq 1 11`; do ls clover${i}_fastqs/clover${i}.combined.fastq; cat clover${i}_fastqs/clover${i}.combined.fastq >> clover_old_combined_fastqs.fastq; done"

sbatch --nodes=1 --mem=15000 --ntasks-per-node=4 -t 1-0 -p msn --wrap="minimap2 -a -x map-ont clover_correct/clover_correct.contigs.fasta clover_old_combined_fastqs.fastq | samtools sort -T clov_combined.old.temp -o clover_old_combined_fastqs.sort.bam -"

### The kmer composition test ###
# Let's start by making Jellyfish hashes from everything in existence
module load miniconda

for i in clover_new_combined_fastqs.fastq clover_correct/clover_correct.contigs.fasta short_read_red_clover_ena.fasta; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=5 --wrap="jellyfish count -m 21 -s 90G -t 5 -o ${i}.jf -C $i"; done

# Now to load a virtual environment in which I've installed Kat
source activate /KEEP/rumen_longread_metagenome_assembly/kat

sbatch --nodes=1 --mem=50000 --ntasks-per-node=8 -p msn --wrap="source activate /KEEP/rumen_longread_metagenome_assembly/kat; kat comp -t 8 -o assembly_kmer_comps clover_new_combined_fastqs.fastq.jf short_read_red_clover_ena.fasta.jf clover_correct/clover_correct.contigs.fasta.jf"

# Now to print out a histogram for each dataset respectively
sbatch --nodes=1 --mem=50000 --ntasks-per-node=15 -t 1-0 -p msn --wrap="source activate /KEEP/rumen_longread_metagenome_assembly/kat; kat hist -t 15 -o clover_new_kat_hist clover_new_combined_fastqs.fastq.jf"
sbatch --nodes=1 --mem=50000 --ntasks-per-node=15 -t 1-0 -p msn --wrap="source activate /KEEP/rumen_longread_metagenome_assembly/kat; kat hist -t 15 -o clover_srasm_kat_hist short_read_red_clover_ena.fasta.jf"
sbatch --nodes=1 --mem=50000 --ntasks-per-node=15 -t 1-0 -p msn --wrap="source activate /KEEP/rumen_longread_metagenome_assembly/kat; kat hist -t 15 -o clover_lrasm_kat_hist clover_correct/clover_correct.contigs.fasta.jf"


### The contamination test
# Let's run this across the read sets and the assemblies and see which is likely to be the source of contamination.
# I'm guessing that human contamination is the most likely reason, but there could be other plant species.
for i in clover_new_combined_fastqs.fastq clover_old_combined_fastqs.fastq clover_correct/clover_correct.contigs.fasta short_read_red_clover_ena.fasta; do echo $i; sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -t 1-0 -p msn --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/centrifuge/centrifuge --report-file ${i}.ctfg.reports -x /beegfs/project/bostauruscnv/novel_seq/cddr_pipeline/refseq_abv_hum_mus_bos -k 1 -f $i > ${i}.ctfg.out"; done
```


```bash
module load canu/latest
### New assembly New reads only
sbatch --nodes=1 --ntasks-per-node=30 --mem=10G -p msn --wrap="canu -p clover_no_old -d clover_no_old genomeSize=420m correctedErrorRate=0.105 overlapper=mhap utgReAlign=true 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p msn' -nanopore-raw clover_new_combined_fastqs.fastq"
```

## New test with Flye and more reads

OK, the previous assemblies weren't fantastic and I think it's primarily due to the fact that we had less data  in the new read only version and we had to tweak tons of settings in canu. Let's shake things up with more data (one more flowcell: clover15) and with the flye assembly.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
module load miniconda
module load samtools
conda create --prefix /KEEP/rumen_longread_metagenome_assembly/flye
source activate /KEEP/rumen_longread_metagenome_assembly/flye
conda install flye

tar -xvf clover15_fastqs.tar.gz
cat Clover15/20190712_1708_GA10000_FAK17200_f03b9ef2/pass/*.fastq >> clover_new_combined_fastqs.fastq

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn flye --nano-raw clover_new_combined_fastqs.fastq -g 420m -t 70 -i 2 -o clover_new_flye

# It completed and it looks like there are 3000 scaffolds!
samtools faidx scaffolds.fasta
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i scaffolds.fasta.fai
N50 length:     367289526
N50 value:      555657
L50 value:      361

perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < scaffolds.fasta.fai
733752078

# So about the same (less contiguous though) than Canu
```


## Red clover WGS alignment and polishing

> Ceres: /home/derek.bickharhth/forage/analysis/clover_assemblies

```bash
module load bwa samtools
# Gathering our sequence data
ls /project/forage_assemblies/sequence_data/CloverGenome-137246120/*/*/*.gz > clover_hen_fastqs.tab

# Ran it on /project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_correct/clover_correct.contigs.fasta -> cont_check

# Clover_no_old
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="bwa index /project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_no_old/clover_no_old.contigs.fasta"
sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p msn --dependency=afterany:733572 --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b clover_no_old -t clover_hen_fastqs.tab -f /project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_no_old/clover_no_old.contigs.fasta -p msn -m"
```


## New idea: look at kmers and/or read coverage to resolve highly heterozygous regions

I previously tried to look at KAT plots with little success. Let's try a different kmer based approach with [Smudgeplot](https://github.com/KamilSJaron/smudgeplot)

> Ceres: /home/derek.bickharhth/forage/analysis/clover_kmers

```bash
export LD_LIBRARY_PATH=/usr/lib64:/lib64:$LD_LIBRARY_PATH
module load jellyfish2/2.2.9;
# Creating the Jellyfish DB
ls /project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/*/*.gz | perl -lane 'print "gunzip -c $F[0]";' > red_clover_reads.generators

sbatch --nodes=1 --mem=300000 --ntasks-per-node=60 -p msn --wrap="jellyfish count -C -m 21 -s 100M -t 60 --bf-size 50G -g red_clover_reads.generators -o red_clover_21mer_illumina"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p msn --wrap='jellyfish histo red_clover_21mer_illumina > red_clover_21mer_illumina.hist'

# Determining cutoff
~/forage/binaries/smudgeplot/exec/smudgeplot.py cutoff red_clover_21mer_illumina.hist L
Running smudgeplot v0.2.1
Task: cutoff
11
Done!

~/forage/binaries/smudgeplot/exec/smudgeplot.py cutoff red_clover_21mer_illumina.hist U
Running smudgeplot v0.2.1
Task: cutoff
1500
Done!

# Running the script
jellyfish dump -c -L 11 -U 1500 red_clover_21mer_illumina | ~/forage/binaries/smudgeplot/exec/smudgeplot.py hetkmers -o red_clover_21mer_illumina.kmer_pairs

## It did not complete quickly. I'm going to try KMC instead
ls /project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/*/*.gz > red_clover_reads.list

mkdir tmp; sbatch --nodes=1 --mem=200000 --ntasks-per-node=60 -p msn --wrap='/project/forage_assemblies/binaries/KMC/bin/kmc -k21 -t60 -m200 -cs10000 @red_clover_reads.list red_clover_kmc tmp'
sbatch --nodes=1 --mem=50000 --ntasks-per-node=2 -p msn --wrap='/project/forage_assemblies/binaries/KMC/bin/kmc_tools transform red_clover_kmc histogram red_clover_kmc_k21.hist -cx10000'

~/forage/binaries/smudgeplot/exec/smudgeplot.py cutoff red_clover_kmc_k21.hist L
~/forage/binaries/smudgeplot/exec/smudgeplot.py cutoff red_clover_kmc_k21.hist U

# same as above, so pretty consistent!
sbatch --nodes=1 --mem=50000 --ntasks-per-node=2 -p msn --wrap='/project/forage_assemblies/binaries/KMC/bin/kmc_tools transform red_clover_kmc -ci11 -cx1500 dump -s red_clover_kmc_k21.dump'
sbatch --nodes=1 --mem=300000 --ntasks-per-node=2 -p msn --wrap='~/forage/binaries/smudgeplot/exec/smudgeplot.py hetkmers -o red_clover_kmc.kmer_pairs < red_clover_kmc_k21.dump'

~/forage_assemblies/binaries/smudgeplot/exec/smudgeplot.py plot red_clover_kmc.kmer_pairs_coverages.tsv
```

That generated two plots that are very, very confused! The heterozygosity of the genome confuses these models quite a bit. I ran the kmer profile in [Genomescope](http://qb.cshl.edu/genomescope/analysis.php?code=Yws2xXiXHsLcbIK4TQbm) and it looks like there is a clear peak at the 25X coverage region that could be the unique sequence. 

## Pseudohaploid correction

Michael Schatz's group just published the Pineapple genome with a [new algorithm](https://github.com/schatzlab/pseudohaploid) that resolves highly heterozygous assemblies. Let's test it out. I will run it with the default parameters for the clover correct (best Canu assembly) and the clover flye (flye assembly).

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
module load mummer/3.23

sbatch --nodes=1 --mem=60000 --ntasks-per-node=2 -p msn /project/rumen_longread_metagenome_assembly/binaries/pseudohaploid/create_pseudohaploid.sh clover_correct/clover_correct.contigs.fasta clover_correct.pseudohaplo

```

Damn! It didn't work because the cluster guys didn't compile mummer with 64-bit capacity! I'll have to readapt the program to minimap2 or something else. 

OK, I've rewritten things as a snakemake pipeline using minimap2. Let's see if I can run this and make it work!

```bash
module load minimap2/2.6
cp ~/rumen_longread_metagenome_assembly/binaries/pseudohaploid/default.json ./

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p msn snakemake -s /project/rumen_longread_metagenome_assembly/binaries/pseudohaploid/snakeFile --cluster "sbatch --nodes=1 --ntasks-per-node=3 --mem=20000 --partition=short" --jobs 999

# Minimap was not giving the granularity that the program required, so I compiled a 64-bit version of nucmer and I am running that now
sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/MUMmer3.23/nucmer --maxmatch -c 100 -l 500 clover_correct/clover_correct.contigs.fasta clover_correct/clover_correct.contigs.fasta -p clover_selfnucmer"

sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p msn --wrap='/project/rumen_longread_metagenome_assembly/binaries/MUMmer3.23/delta-filter -l 1000 -i 90 clover_selfnucmer.delta > clover_selfnucmer.filter.delta'

sbatch --nodes=1 --mem=25000 --ntasks-per-node=3 -p msn --wrap='/project/rumen_longread_metagenome_assembly/binaries/MUMmer3.23/show-coords -rclH clover_selfnucmer.filter.delta > clover_selfnucmer.filter.coords'

perl /project/rumen_longread_metagenome_assembly/binaries/pseudohaploid/pseudohaploid.chains.pl clover_selfnucmer.filter.coords 90 93 20000 > clover_selfnucmer.filter.chains

grep '^#' clover_selfnucmer.filter.chains | awk -v cut=93 '{if ($4 >= cut){print ">"$2}}' > clover_selfnucmer.contained.ids

wc -l clover_selfnucmer.contained.ids
376 clover_selfnucmer.contained.ids

/project/rumen_longread_metagenome_assembly/binaries/pseudohaploid/filter_seq -v clover_selfnucmer.contained.ids clover_correct/clover_correct.contigs.fasta > clover_correct.pseudohap.fa
samtools faidx clover_correct.pseudohap.fa

wc -l clover_correct.pseudohap.fa.fai
2759 clover_correct.pseudohap.fa.fai
wc -l clover_correct/*.fai
3135 clover_correct/clover_correct.contigs.fasta.fai

perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < clover_correct.pseudohap.fa.fai
721484706
perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < clover_correct/clover_correct.contigs.fasta.fai
748498712

# Only a reduction in 20 megabases...
```

OK, last attempt with the pipelines: let's try the [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) alignment suite.

```bash
module load minimap2/2.6 r/3.5.2 bedtools/2.25.0 samtools
ln -s /project/rumen_longread_metagenome_assembly/binaries/purge_haplotigs/bin/purge_haplotigs ~/bin/purge_haplotigs

mkdir purge_haplotigs
sbatch --nodes=1 --mem=14000 --ntasks-per-node=4 -p short --wrap='minimap2 -t 4 -ax map-pb clover_correct/clover_correct.contigs.fasta clover_correct/clover_correct.correctedReads.fasta.gz --secondary=no | samtools sort -m 1G -o purge_haplotigs/aligned.bam -T tmp.align'

# Step 1
sbatch --nodes=1 --ntasks-per-node=10 --mem=32000 -p msn --wrap='purge_haplotigs hist -b purge_haplotigs/aligned.bam -g clover_correct/clover_correct.contigs.fasta -t 10'

# I had to modify the script to print to pdf, but then I saw that the distribution of read depth was unimodal
# Odd... there should be bimodal peaks if this is true. perhaps I smashed everything too much with the clover_correct assembly? Let's run two tests.
```

Quickly making a flye assembly of all Clover reads.

> Ceres: /project/forage_assemblies/analysis/clover_assemblies

```bash
module load miniconda/3.6
cat /project/forage_assemblies/sequence_data/clover_old_combined_fastqs.fastq /project/forage_assemblies/sequence_data/clover_new_combined_fastqs.fastq > /project/forage_assemblies/sequence_data/clover_total_combined_fastqs.fastq

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn --wrap='source activate /KEEP/rumen_longread_metagenome_assembly/flye; flye --nano-raw /project/forage_assemblies/sequence_data/clover_total_combined_fastqs.fastq -g 420m -t 70 -i 2 -o clover_total_flye'


```

> Ceres: /project/forage_assemblies/assemblies/red_clover

```bash
module load miniconda
source activate /KEEP/rumen_longread_metagenome_assembly/flye
# Trying with the metrics that Serge used for the mammalian genomes:
sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn -q msn flye --nano-raw /project/forage_assemblies/sequence_data/clover_total_combined_fastqs.fastq -g 420m -t 70 -i 2 -m 10000 --asm-coverage 40 -o clover_limit_flye
```

Now, I'm going to try the purge haplotigs pipeline again before trying a modified version of pseudohaploid again. I think that the limit flye assembly actually comes out better than the canu assembly with fewer small unplaced contigs. We will see! The Pseudohaploid run will use the less strict -k19 -w19 settings found in [Ragoo](https://github.com/malonge/RaGOO/blob/master/ragoo.py). 

```bash
module load minimap2/2.6 r/3.5.2 bedtools/2.25.0 samtools

mkdir flye_purgehaplotigs
sbatch --nodes=1 --mem=14000 --ntasks-per-node=9 -p short -q memlimit --wrap='minimap2 -t 6 -ax map-pb clover_limit_flye/assembly.fasta /project/forage_assemblies/sequence_data/clover_total_combined_fastqs.fastq --secondary=no | samtools sort -m 1G -@ 2 -o flye_purgehaplotigs/flye_aligned.bam -T tmp.align'

# Step 1
sbatch --nodes=1 --ntasks-per-node=10 --mem=48000 -p short -q memlimit --wrap='samtools index flye_purgehaplotigs/flye_aligned.bam; purge_haplotigs hist -b flye_purgehaplotigs/flye_aligned.bam -g clover_limit_flye/assembly.fasta -t 10'

# OK! This got me somewhere! I now have a low and midpoint cutoff rate from the pdf file produced by the program
# Low point: 35
# Mid point: 95
# high point: 160

# Step 2
sbatch --nodes=1 --ntasks-per-node=1 --mem=8000 -p short -q memlimit --wrap='purge_haplotigs cov -i flye_aligned.bam.gencov -l 35 -h 160 -m 95 -o flye_clover_stats.csv'
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f flye_clover_stats.csv -d ',' -c 1 -m -i '#'
|Entry | Value|
|:-----|-----:|
|s     |  2718|
|      |   380|	<- this is where the "junk" ratio is not met and the "haplotig" ratio is not met. 90% diploid
|j     |   124|


# Step 3
sbatch --nodes=1 --ntasks-per-node=10 --mem=48000 -p short -q memlimit --wrap="purge_haplotigs purge -g clover_limit_flye/assembly.fasta -c flye_clover_stats.csv -t 10 -o flye_clover_curated -d -b flye_purgehaplotigs/flye_aligned.bam"

# note that there was a pipe failure at the end of some of the multicov bedtools operations for some reason
samtools faidx flye_clover_curated.artefacts.fasta
samtools faidx flye_clover_curated.fasta
samtools faidx flye_clover_curated.haplotigs.fasta

wc -l *.fai
   124 flye_clover_curated.artefacts.fasta.fai
  1597 flye_clover_curated.fasta.fai
  1501 flye_clover_curated.haplotigs.fasta.fai
  3222 total

for i in *.fai; do echo $i;  perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i $i; done
flye_clover_curated.artefacts.fasta.fai
	N50 length:     3119077
	N50 value:      73472
	L50 value:      30
flye_clover_curated.fasta.fai
	N50 length:     216376573
	N50 value:      514218
	L50 value:      239
flye_clover_curated.haplotigs.fasta.fai
	N50 length:     133405088
	N50 value:      357870
	L50 value:      228

for i in *.fai; do echo -ne "$i\t"; perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < $i; done
	flye_clover_curated.artefacts.fasta.fai 	6122171
	flye_clover_curated.fasta.fai   			431899544
	flye_clover_curated.haplotigs.fasta.fai 	266559604

for i in *.fai; do perl -lane 'print $F[0];' < $i > $i.ctg.list; done
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl flye_clover_curated.artefacts.fasta.fai.ctg.list flye_clover_curated.fasta.fai.ctg.list flye_clover_curated.haplotigs.fasta.fai.ctg.list
	File Number 1: flye_clover_curated.artefacts.fasta.fai.ctg.list
	File Number 2: flye_clover_curated.fasta.fai.ctg.list
	File Number 3: flye_clover_curated.haplotigs.fasta.fai.ctg.list
	Set     Count
	1       124
	2       1597
	3       1501

# So they're all unique. Let's compare that to the Schatz pipeline later
```

## Attempting to hack the trio canu pipeline

I think that I can pull apart reads from the heterozygous site linkage groups using some informatics analysis of read kmers. My goal is to try to identify kmers that are present in a "lower peak" (ie. haplotype specific) and then calculate their co-occurrence on error-prone long reads to try to estimate linkage-group-specific reads. Then I assemble both pools of reads separately.

From what I've found in the canu source code, [splithaplotype.c](https://github.com/marbl/canu/blob/master/src/haplotyping/splitHaplotype.C) determines the minima and maxima kmer occurrence peaks from a histogram and then pulls kmers from each parental haplotype from the pool in between. This avoids dealing with "noise" kmers and kmers that are shared by both parents. 

Let's start by plotting the exact jellyfish kmer profile to try to select the best frequency cutoffs.

> Ceres: /home/derek.bickharhth/forage/analysis/clover_kmers

```R
library(ggplot2)
data <- read.delim("red_clover_21mer_illumina.hist", sep=" ", header=FALSE)

pdf(file="red_clover_21mer_illumina.hist.pdf", useDingbats=FALSE)
ggplot(data, aes(x = V1, y= V2)) + geom_col(color="darkblue", fill="lightblue") + xlim(0,100) + ylim(0,8000000)
dev.off()
```

Basically produced the same plot as Genomescope (see above) but a little clearer on the kmer frequencies. It looks like a minimum cutoff of 15 and a maximum of 30 are warranted by visual inspection. Let's dump the kmers to see how much we're working with here.

```bash
/project/forage_assemblies/binaries/KMC/bin/kmc_dump -ci15 -cx30 red_clover_kmc red_clover_kmc_k21.uniq_het.peak.kmers

wc -l red_clover_kmc_k21.uniq_het.peak.kmers
99725171 red_clover_kmc_k21.uniq_het.peak.kmers

# That's quite alot!
```

I need to use [directed acyclic graphs](https://algs4.cs.princeton.edu/41graph/) to organize the data, and calculate the [minimizer values](http://yorke.umd.edu/papers/genome_papers/minimizerpaper.pdf) for kmers and use a concept of [granularity in hierarchically organized kmers](https://medium.com/@infoecho/constructing-a-graph-for-genome-comparison-swiftly-d47dcd7eae5d) to try to compare and sort reads. I think that most of what I'm thinking is far too complex and that I should start with the generation of a frequency 'graph' to print edge weights (edges are counts of observances).

Let's start out by pulling the "one bp difference" kmers from smudgeplot out by using the het kmer counts we just found above.

```bash
sbatch --nodes=1 --mem=30000 --ntasks-per-node=2 -p priority -q msn intersect_hetpeak_with_kmerpairs.pl red_clover_kmc_k21.uniq_het.peak.kmers red_clover_kmc.kmer_pairs_sequences.tsv red_clover_uniq_het.kmer_pairs.sequences
```

#### Advice from Arang:
> meryl-lookup -memory 12 -include -mers $meryl -sequence $read1 -sequence2 $read2 -r2 $out.R2.fastq.gz
> From Meryl version 1

## Testing purge_dups

OK, I'm going to try to use [purge_dups](https://github.com/dfguan/purge_dups) instead. It looks like the pipeline manager tool that the author created won't work on Ceres (no QOS terms), so I need to create manual scripts to queue everything up.

> Ceres: /lustre/project/forage_assemblies/analysis/clover_kmers/

```bash

```

## Alfalfa side project

> CERES: /project/forage_assemblies/analysis/alfalfa_project

```bash
module load minimap2/2.6

# mapping the reads to the target sequence
minimap2 -x map-ont t_dna_pmls485.fa /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq > gusalfb1_nanopore_mappings.paf

# Now I'm going to steal my viral overlap script and modify it so that a bed file can be generated
cp ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAnalysisScripts/selectLikelyViralOverhangs.pl ./
cp ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAnalysisScripts/filterViralOverhangsAndGenerateSeq.pl ./

perl selectLikelyViralOverhangs.pl gusalfb1_nanopore_mappings.paf gusalfb1_overhangs

# OK, there aren't many reads! I will start out with the whole list (22 reads) and then move onto the filtered list if needed
# Full set of reads that had ANY alignments whatsoever
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %reads; while(<IN>){chomp; @s = split(/\t/); $reads{$s[0]} = 1;} close IN; open(IN, "< $ARGV[1]"); while($h = <IN>){$h =~ s/^\@//; if(exists($reads{$h})){print "\@$h\n"; $s = <>; $p = <>; $q = <>; print "$s$p$q";}else{<>; <>; <>;}} close IN;' gusalfb1_nanopore_mappings.paf /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq > full_align_filtered_reads.fastq

# I had to queue this on the cluster because it was taking a long time
sbatch --nodes=1 --mem=6000 --ntasks-per-node=2 -p msn -q msn --wrap="perl filter_fastq_file.pl gusalfb1_nanopore_mappings.paf /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq full_align_filtered_reads.fastq"

# OK, let's try assembling these with canu
module load canu/1.9
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p short -q memlimit --wrap="canu -p gus_alfalfa -d canu_alfalfa genomeSize=10k corMhapSensitivity=high corMinCoverage=0 contigFilter='2 0 1.0 0.5 0' 'gridOptions=-p short -q memlimit' -nanopore-raw full_align_filtered_reads.fastq"

grep '>' canu_alfalfa/gus_alfalfa.contigs.fasta
>tig00000001 len=52946 reads=10 class=contig suggestRepeat=no suggestCircular=no

# Jackpot? Let's align the old segment to this one to check
sbatch --nodes=1 --mem=12000 --ntasks-per-node=3 -p msn -q msn --wrap='minimap2 -x asm10 canu_alfalfa/gus_alfalfa.contigs.fasta t_dna_pmls485.fa > gus_alfalfa_to_tdna.paf'
# The output
T-DNA_from_pMLS485      5497    80      5441    -       tig00000001     52946   30491   35791   4348    5366    60      tp:A:P  cm:i:403        s1:i:4343       s2:i:91 dv:f:0.0166
# Right smack in the middle! Lots of errors though.

# Let's try 3 billion iterations (2) of racon to polish this thing
# First, align all reads to this larger contig
sbatch --nodes=1 --mem=12000 --ntasks-per-node=3 -p msn -q msn --wrap='minimap2 -x map-ont canu_alfalfa/gus_alfalfa.contigs.fasta /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq > gus_alfalfa_contig_map_it1.paf'

# There were allot more reads aligning than I expected! I think that maybe Canu joined a chimeric contig here?
module load bedtools/2.25.0
perl -lane 'print "$F[5]\t$F[7]\t$F[8]";' < gus_alfalfa_contig_map_it1.paf | bedtools sort -i stdin | bedtools intersect -a gus_alfalfa.contigs.len.bed -b stdin -wao | less
...
tig00000001     1       52946   tig00000001     16374   17096   722
tig00000001     1       52946   tig00000001     16374   17132   758
tig00000001     1       52946   tig00000001     20410   52937   32527
tig00000001     1       52946   tig00000001     22789   38127   15338
tig00000001     1       52946   tig00000001     23465   38519   15054
...

# So, let's crop the contig at 20410 bp onwards
module load samtools
samtools faidx canu_alfalfa/gus_alfalfa.contigs.fasta tig00000001:20410-52946 > gus_alfalfa.contigs.cropped.fasta
# I renamed the head of the fasta file to "target_region"
# Let's realign then!
sbatch --nodes=1 --mem=12000 --ntasks-per-node=3 -p msn -q msn --wrap='minimap2 -x asm10 gus_alfalfa.contigs.cropped.fasta t_dna_pmls485.fa > gus_alfalfa_cropped_to_tdna.paf'
# Now, 10kb from the start
T-DNA_from_pMLS485      5497    80      5441    -       target_region   32537   10082   15382   4348    5366    60      tp:A:P  cm:i:403        s1:i:4343       s2:i:91 dv:f:0.0166

sbatch --nodes=1 --mem=12000 --ntasks-per-node=3 -p msn -q msn --wrap='minimap2 -x map-ont gus_alfalfa.contigs.cropped.fasta /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq > gus_alfalfa_contig_map_it2.paf'

# There's still a ton of repetitive chaff aligning to the contig. Going to filter it
perl -lane 'if($F[9] > 1000){print $_;}' < gus_alfalfa_contig_map_it2.paf > gus_alfalfa_contig_map_it2.filt.paf
# It's racon time!
module load racon/1.3.0
sbatch --nodes=1 --mem=30000 --ntasks-per-node=10 -p msn -q msn --wrap='racon -t 10 /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq gus_alfalfa_contig_map_it2.filt.paf gus_alfalfa.contigs.cropped.fasta > gus_alfalfa.contigs.racon1.fasta'

# Now for racon iteration 2!
sbatch --nodes=1 --mem=12000 --ntasks-per-node=3 -p msn -q msn --wrap='minimap2 -x map-ont gus_alfalfa.contigs.racon1.fasta /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq > gus_alfalfa_contig_map_it3.paf'

# note: did not filter the paf this time. fingers crossed!
sbatch --nodes=1 --mem=30000 --ntasks-per-node=10 -p msn -q msn --wrap='racon -t 10 /project/forage_assemblies/sequence_data/gusalfb1_nanopore_raw.fastq gus_alfalfa_contig_map_it3.filt.paf gus_alfalfa.contigs.racon1.fasta > gus_alfalfa.contigs.racon2.fasta'

# That was horrible and collapsed the fasta! we're sticking with the first round results.
```

We increased the X coverage by double in the meantime. I am going to rerun the pipeline and pick out only the premium reads.

```bash
sbatch --nodes=1 --mem=16000 --ntasks-per-node=3 -p msn -q msn --wrap="minimap2 -x map-ont t_dna_pmls485.fa /project/forage_assemblies/sequence_data/gusalfb1_total_nanopore_raw.fastq > gusalfb1_nanopore_mappings_total.paf"

perl -lane 'if($F[9] > 2000 && $F[10] > 5000 && $F[11] >= 60){print $_;}' < gusalfb1_nanopore_mappings_total.paf | wc -l
17  # These are the reads with good alignment to the insert, good mapping quality and reasonable mapping rates

perl -lane 'if($F[9] > 2000 && $F[10] > 5000 && $F[11] >= 60){print $_;}' < gusalfb1_nanopore_mappings_total.paf | perl -lane 'print $F[0];' > gusalfb1_nanopore_mappings_total.hqreadnames.list

perl filter_fastq_file.pl gusalfb1_nanopore_mappings_total.hqreadnames.list /project/forage_assemblies/sequence_data/gusalfb1_total_nanopore_raw.fastq gusalfb1_nanopore_mappings_total.hqreads.fastq

module load canu/1.9
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="canu -p gus_alfalfa -d canu_totalfalfa genomeSize=10k corMhapSensitivity=high corMinCoverage=0 contigFilter='2 0 1.0 0.5 0' 'gridOptions=-p msn -q msn' -nanopore-raw gusalfb1_nanopore_mappings_total.hqreads.fastq"

# Wow! It assembled two contigs. I'll check to see which one contains the insert
minimap2 -x asm10 canu_totalfalfa/gus_alfalfa.contigs.fasta t_dna_pmls485.fa > canu_totalalfalfa.tdna.align.paf

grep '>' canu_totalfalfa/gus_alfalfa.contigs.fasta
>tig00000007 len=49247 reads=8 class=contig suggestRepeat=no suggestCircular=no
>tig00000008 len=26084 reads=1 class=contig suggestRepeat=no suggestCircular=no

# I think that the first contig is probably the best. Maybe the second is a chimera?
# Ah, a clustalomega alignment shows that the second contig is a subset of the first that was separated for some reason. The only sequence divergence is at the end of the read.


## Let's see how the pileups look
samtools faidx canu_totalfalfa/gus_alfalfa.contigs.fasta tig00000007 > gus_total_fasta_tig7.fasta

minimap2 -a -x map-ont gus_total_fasta_tig7.fasta canu_totalfalfa/gus_alfalfa.correctedReads.fasta.gz | samtools sort -T gus -o gus_alfalfa_corrected_tig7.bam -
samtools index gus_alfalfa_corrected_tig7.bam

samtools tview gus_alfalfa_corrected_tig7.bam gus_total_fasta_tig7.fasta

## Doing this on the uncorrected reads
minimap2 -a -x map-ont gus_total_fasta_tig7.fasta gusalfb1_nanopore_mappings_total.hqreads.fastq | samtools sort -T gus -o gus_alfalfa_uncorrected_tig7.bam
samtools index gus_alfalfa_uncorrected_tig7.bam

```


## Polishing the assembly for public release

I want to run a quick pilon polish to make sure that the assembly is ready for use by other groups. I will use the clover hen corrected canu assembly for this.

> Ceres: /project/forage_assemblies/assemblies/red_clover

```bash
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p priority -q msn --wrap="bwa index clover_limit_flye/assembly.fasta"

ls /project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/*/*fastq.gz > clover_hen_illumina_sequence.tab

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b pilon_first -t clover_hen_illumina_sequence.tab -f /project/forage_assemblies/assemblies/red_clover/clover_limit_flye/assembly.fasta -m -q msn -p priority

python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f /project/forage_assemblies/assemblies/red_clover/pilon_first/clover/clover.sorted.merged.bam -g /project/forage_assemblies/assemblies/red_clover/clover_limit_flye/assembly.fasta -q msn -p priority -o pilon_first_correction -m

# Now to realign to the new assembly!
cat pilon_first_correction/*.fasta > pilon_first_correction/clover_pilonone_corrected.fasta
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn --wrap="module load bwa; bwa index pilon_first_correction/clover_pilonone_corrected.fasta"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b pilon_second -t clover_hen_illumina_sequence.tab -f /project/forage_assemblies/assemblies/red_clover/pilon_first_correction/clover_pilonone_corrected.fasta -m -q msn -p priority
python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f /project/forage_assemblies/assemblies/red_clover/pilon_second/clover/clover.sorted.merged.bam -g /project/forage_assemblies/assemblies/red_clover/pilon_first_correction/clover_pilonone_corrected.fasta -q msn -p priority -o pilon_second_correction -m

# Creating the last version of the assembly
perl -e '@f = `ls pilon_second_correction/*.fasta`; open(OUT, "> clover_second_pilon_assembly.fasta"); foreach $j (@f){open(IN, "< $j"); while(<IN>){if($_ =~ /^>/){$_ =~ s/_pilon//g;} print OUT $_;} close IN;} close OUT;'

# I'm just curious on how the polishing performed. Let's generate some stats
sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p priority -q msn --wrap="nucmer -maxmatch -l 100 -c 500 clover_limit_flye/assembly.fasta clover_second_pilon_assembly.fasta -prefix pre_vs_post_corr_clover"
```

## Polished assembly QC

Queueing up my assembly QC pipeline to provide stats on Clover

> Ceres: /lustre/project/forage_assemblies/assemblies/red_clover

```bash
module load miniconda/3.6

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```

## Clover recalling

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
for i in CloverB Clover2 Clover; do echo $i; sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/bin/guppy_basecaller -i $i/*/fast5_pass -s $i/g3_6_f -c /lustre/project/rumen_longread_metagenome_assembly/binaries/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --cpu_threads_per_caller 14 --num_callers 5"; done
```

## Clover CCS

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/flye/

sbatch -N 1 -n 80 --mem=1100000 -p priority-mem -q msn-mem -t 5-0 flye --pacbio-hifi /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_200929_165102.Q20.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_201002_155306.Q20.fastq -g 800m -t 78 -o clover_hiflye

vim default.json

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 4-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 4-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda

### IPA assembly
conda activate /KEEP/rumen_longread_metagenome_assembly/ipa

sbatch -N 1 -n 2 --mem=30000 -p priority -q msn -t 3-0 ipa dist --nthreads 24 --njobs 100 -i /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_200929_165102.Q20.fastq -i /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_201002_155306.Q20.fastq --run-dir clover_ipa --cluster-args "sbatch -N 1 -n {params.num_threads} -o logs/ --mem=40000 -p priority -q msn -t 3-0"

# It created a config file, but it didn't point to the correct snakemake file. I had to dig through and requeue manually using the suggested params in the slurm-out stdout file
sbatch -N 1 -n 2 --mem=15000 -p priority -q msn -t 4-0 --wrap="/KEEP/rumen_longread_metagenome_assembly/ipa/bin/python3 -m snakemake -j 100 -d clover_ipa -p -s /KEEP/rumen_longread_metagenome_assembly/ipa/etc/ipa.snakefile --configfile clover_ipa/config.yaml --reason --cluster 'sbatch -N 1 -n {params.num_threads} -o stdout.{rule}.out --mem=85000 -p priority -q msn -t 4-0' --latency-wait 60 --rerun-incomplete"
```

#### Generating meryl venn diagram for comparison


> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/meryl

sbatch -N 1 -n 4 --mem=30000 -p priority -q msn -t 2-0 --wrap="python ~/python_toolchain/sequenceData/merylVennUpset.py -m /lustre/project/rumen_longread_metagenome_assembly/binaries/meryl/build/bin/meryl -o clover_hifi_comp --d merqury/altIPA/altIPA.meryl -d merqury/primaryIPA/primaryIPA.meryl -d merqury/hifiASM/hifiASM.meryl -d merqury/hifiFlye/hifiFlye.meryl -d merqury/ontFlye/ontFlye.meryl"
```

#### Testing out website generations script

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs

```bash
conda env create -f ~/python_toolchain/snakeMake/assemblyValidation/envs/base.yaml -p ninja
conda activate /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/ninja

python ~/python_toolchain/snakeMake/assemblyValidation/scripts/createWebpage.py -f final -o clover_summary -c hifiASM_altIPA -c hifiASM_hifiFlye -c hifiASM_primaryIPA -c hifiFlye_altIPA -c hifiFlye_primaryIPA -c ontFlye_altIPA -c ontFlye_hifiASM -c ontFlye_hifiFlye -c ontFlye_primaryIPA -c primaryIPA_altIPA -s clover_limit_flye.fasta -s red_clover_2cells_hifiasm.p_ctg.fasta -s clover_hiflye.fasta -s final.p_ctg.fasta -s final.a_ctg.fasta -a ontFlye -a hifiASM -a hifiFlye -a primaryIPA -a altIPA

mkdir package

wget https://zenodo.org/record/17232/files/redclover_v2.1.fasta
```


#### Clover hic qc

```bash
module load samtools bwa miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/hic_qc
sbatch -N 1 -n 32 --mem=55000 -p priority -q msn --wrap='bwa mem -t 28 /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/fastas/primaryIPA.fa LIB109270_S1_L001_R1_001.fastq.gz LIB109270_S1_L001_R2_001.fastq.gz | samtools sort -T clover_hic -n -@ 4 -o clover_hic.sorted.bam -'

sbatch -N 1 -n 2 --mem=12000 -p priority -q msn --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/hic_qc/hic_qc.py -b clover_hic.sorted.bam --outfile_prefix  clover_hic "

## Running on the additional data
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/fastas/primaryIPA.fa -b omnic_ipa_paligns -t omnic_data.tab -q msn -p priority -n -m

sbatch -N 1 -n 2 --mem=12000 -p priority -q msn --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/hic_qc/hic_qc.py -b omnic_ipa_paligns/clover/clover.sorted.merged.bam --outfile_prefix  clover_omnic "



sbatch -N 1 -n 1 --mem=45000 -p priority -q msn --wrap="bwa index clover_primary_ipa.purged.fa"
sbatch -N 1 -n 1 --mem=9000 --dependency=afterok:5851618 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_primary_ipa.purged.fa -b omnic_purge_paligns -t omnic_data.tab -q msn -p priority -n -m"
```

#### Clover salsa scaffolding

```bash
module load bedtools/2.25.0 miniconda/3.6 minimap2
conda activate /KEEP/rumen_longread_metagenome_assembly/salsa

bamToBed -i omnic_ipa_paligns/clover/clover.sorted.merged.bam > omnic_ipa_paligns/clover.sorted.bed
sbatch -N 1 -n 2 --mem=280000 -p priority -q msn --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/SALSA/run_pipeline.py -a fastas/primaryIPA.fa -l fastas/primaryIPA.fa.fai -e DNASE -b omnic_ipa_paligns/clover.sorted.bed -o omnic_ipa_paligns/salsa_scaffolds"

# Now for ASM-ASM comparisons
sbatch -N 1 -n 3 --mem=50000 -p priority -q msn --wrap="minimap2 -x asm5 fastas/primaryIPA.fa omnic_ipa_paligns/salsa_scaffolds/scaffolds_FINAL.fasta > omnic_ipa_paligns/rc_pipa_to_salsa.paf"

conda activate /KEEP/rumen_longread_metagenome_assembly/r
Rscript /lustre/project/rumen_longread_metagenome_assembly/binaries/Themis-ASM/scripts/pafDotPlotly.R -i omnic_ipa_paligns/rc_pipa_to_salsa.paf -v -l -s -o rc_pipa_to_salsa

# Running scaffolding on purged primary IPA assembly
sbatch -N 1 -n 1 --mem=10000 -p priority -q msn --wrap="bamToBed -i omnic_purge_paligns/clover/clover.sorted.merged.bam > omnic_purge_paligns/clover.sorted.bed"

sbatch -N 1 -n 2 --mem=280000 -p priority -q msn --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/SALSA/run_pipeline.py -a clover_primary_ipa.purged.fa -l clover_primary_ipa.purged.fa.fai -e DNASE -b omnic_purge_paligns/clover.sorted.bed -o omnic_purge_paligns/salsa_scaffolds"

# Now for an ASM-ASM comparison

```

The scaffolding wasn't the best. Maybe I should purge dups first?

## Purge dups

```bash
module load miniconda/3.6 "java/11.0.2" bedtools minimap2

conda activate /KEEP/rumen_longread_metagenome_assembly/purge_dups/

for i in /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/*Q20.fastq; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 15 -p priority -q msn --mem=55000 --wrap="minimap2 -t 15 -x asm20 clover_ipa/19-final/final.p_ctg.fasta $i > $name.asm20.paf"; done

for i in *.paf; do sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gzip $i"; done

sbatch -N 1 -n 2 --mem=65000 -p priority -q msn --wrap="pbcstat *asm20.paf.gz; calcuts PB.stat > cutoffs 2> calcults.log"

sbatch -N 1 -n 4 --mem=55000 -p priority -q msn --wrap="split_fa clover_ipa/19-final/final.p_ctg.fasta > ipa_primary.fasta.split; minimap2 -xasm5 -DP ipa_primary.fasta.split ipa_primary.fasta.split | gzip -c - > ipa_primary.split.self.paf.gz"

sbatch -N 1 -n 4 --mem=55000 -p priority -q msn --wrap="get_seqs -l 5000 -e dups.bed clover_ipa/19-final/final.p_ctg.fasta"

sbatch -N 1 -n 1 --mem=20000 -p primary -q msn --wrap="python3 ~/python_toolchain/sequenceData/reformatFasta.py -f clover_primary_ipa.purged.fa -o clover_hifi_temp.fasta -v"

mv clover_hifi_temp.fasta clover_primary_ipa.purged.fa
```

## Themis run

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_ccs_themis

```bash
module load miniconda/3.6 "java/11.0.2" bedtools

PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/redclover_v2.1.fasta -n shortRead -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_primary_ipa.purged.fa -n hifiContigs -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/omnic_purge_paligns/salsa_scaffolds/scaffolds_FINAL.fasta -n hifiScaffolds -b eudicots_odb10 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L001-ds.ac229a7a209e4f809c4ac2bda32e46fe/Hen17Cv_S1_L001_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L001-ds.ac229a7a209e4f809c4ac2bda32e46fe/Hen17Cv_S1_L001_R2_001.fastq.gz -s lane1 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L002-ds.e4d313d6db7a45ec88dba60164d0d717/Hen17Cv_S1_L002_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L002-ds.e4d313d6db7a45ec88dba60164d0d717/Hen17Cv_S1_L002_R2_001.fastq.gz -s lane2 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L003-ds.20e6182d55ac40e2b889f55011234126/Hen17Cv_S1_L003_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L003-ds.20e6182d55ac40e2b889f55011234126/Hen17Cv_S1_L003_R2_001.fastq.gz -s lane3 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L004-ds.8b23d704e97542d982a7416ef5ce3c8a/Hen17Cv_S1_L004_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L004-ds.8b23d704e97542d982a7416ef5ce3c8a/Hen17Cv_S1_L004_R2_001.fastq.gz -s lane4 -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'
```

## Download bac ends

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/bac_ends

```bash
conda activate /KEEP/rumen_longread_metagenome_assembly/eutils
PERL5LIB=/KEEP/rumen_longread_metagenome_assembly/eutils/lib/perl5:/KEEP/rumen_longread_metagenome_assembly/eutils/lib/5.26.2/x86_64-linux-thread-multi/:/KEEP/rumen_longread_metagenome_assembly/eutils/lib/site_perl/:$PERL5LIB

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='for i in `seq 235466 298279`; do echo $i; esearch -db nucleotide -query HR${i} | efetch -format fasta >> bac_ends_unsorted.fasta; done'

perl generate_pairs.pl bac_ends_unsorted.fasta sorted_bac_ends

grep '>' sorted_bac_ends.f.fasta > for.list
grep '>' sorted_bac_ends.r.fasta > rev.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl for.list rev.list
File Number 1: for.list
File Number 2: rev.list
Set     Count
1       2497
1;2     29069
2       2175

# looks like they screwed up their submission! I'm going to have to filter the file to keep only the paired bac ends
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_2 for.list rev.list | perl -lane 'if($F[0] =~ /File/){next;}else{$F[0] =~ s/>//; print "$F[0]";}' > shared.list

perl generate_pairs.pl bac_ends_unsorted.fasta sorted_bac_ends shared.list

module load minimap2
minimap2 -ax sr /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/omnic_purge_paligns/salsa_scaffolds/scaffolds_FINAL.fasta sorted_bac_ends.f.fasta sorted_bac_ends.r.fasta > hifi_bacends_paired.sam
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] & 2048){}else{print $_;}' < hifi_bacends_paired.sam | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 6 -i '@' -d '\t' | perl -e '%h; <>; while(<>){chomp; @s = split(/\t/); if($s[0] eq "="){$h{"same"} += $s[1];}elsif($s[0] eq "*"){$h{"unmap"} += $s[1];}else{$h{"different"} += $s[1];}} foreach $i ("same", "different", "unmap"){print "$i\t$h{$i}\n";}'
same    27984
different       29188
unmap   966

```

|Entry        | Value|
|:------------|-----:|
|same		  |27,984|
|different    |29,188|
|unmap        |   966|

```bash
minimap2 -ax sr /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/redclover_v2.1.fasta sorted_bac_ends.f.fasta sorted_bac_ends.r.fasta > short_bacends_paired.sam

perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] & 2048){}else{print $_;}' < short_bacends_paired.sam | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 6 -i '@' -d '\t' | perl -e '%h; <>; while(<>){chomp; @s = split(/\t/); if($s[0] eq "="){$h{"same"} += $s[1];}elsif($s[0] eq "*"){$h{"unmap"} += $s[1];}else{$h{"different"} += $s[1];}} foreach $i ("same", "different", "unmap"){print "$i\t$h{$i}\n";}'
same    14714
different       42456
unmap   968
```

|Entry                      | Value|Pairs |
|:--------------------------|-----:|-----:|
|same                       |14,714| 7,357|
|different					|42,456|21,228|
|unmap                      |   968|   484|

So, the Salsa scaffold HIFI assembly better resolves BAC end placement on the same scaffolds than the short read assembly! Big time! Almost twice as likely to resolve it!

#### Super scaffold run

```bash
minimap2 -ax sr ../manuscript_themis/fastas/ARSDFRC1.fa sorted_bac_ends.f.fasta sorted_bac_ends.r.fasta > ARSRC1_bacends_paired.sam

perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] & 2048){}else{print $_;}' < ARSRC1_bacends_paired.sam | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 6 -i '@' -d '\t' | perl -e '%h; <>; while(<>){chomp; @s = split(/\t/); if($s[0] eq "="){$h{"same"} += $s[1];}elsif($s[0] eq "*"){$h{"unmap"} += $s[1];}else{$h{"different"} += $s[1];}} foreach $i ("same", "different", "unmap"){print "$i\t$h{$i}\n";}'
same    31590
different       25582
unmap   966
```

|Entry                      | Value|Pairs |
|:--------------------------|-----:|-----:|
|same                       |31,590|15,795|
|different					|25,582|12,791|
|unmap                      |   966|   483|

Let's create a graph so we can compare our scaffold-scaffold aligments with additional evidence.

```bash
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] & 2048){}else{if($F[6] eq "="){$F[6] = $F[2];} print "$F[2]-$F[6]";}' < hifi_bacends_paired.sam | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -d '\t' | perl -e '<>; while(<>){chomp; @s = split(/[-\t]/); if($s[0] eq $s[1]){next;} print join("\t", @s) . "\n";}' > hifi_bacends_discord.bigraph

perl -lane 'print $F[-1];' < hifi_bacends_discord.bigraph | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3962
Sum:    29188
Minimum 1
Maximum 286
Average 7.366986
Median  2
Standard Deviation      16.919970

perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] & 2048){}else{if($F[6] ne "=" && $F[3] < 250000){print "$F[2]-$F[6]";}' < hifi_bacends_paired.sam | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -d '\t'

# I think that the BAC end alignments don't work as well due to my lack of repeatmasking. Let's try again.
sbatch -N 1 -n 50 --mem=55000 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 50 -q -species arabidopsis -no_is -gff scaffolds_FINAL.fasta"
```

## Super scaffolding run

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_ccs_themis

```bash
grep 'LG' mapped/mapshortRead_hifiScaffolds.paf | perl -e '%h; %l; while(<>){chomp; @s = split(/\t/); $h{$s[0]}->{$s[5]} += $s[9]; $l{$s[0]} = $s[1];} foreach $q (sort{$l{$b} <=> $l{$a}} keys(%l)){$num = scalar(keys(%{$h{$q}})); print "$q\t$l{$q}\t$num"; foreach $r (sort{$h{$q}->{$b} <=> $h{$q}->{$a}} keys(%{$h{$q}})){if($h{$q}->{$r} < 1000){next;}print "\t$r:" . $h{$q}->{$r};} print "\n";}' > scaffold_lg_linkage.tab

# I think that one of my old perl scripts can help here
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -g mapped/mapshortRead_hifiScaffolds.paf -m 5000 -o shortread_superscaffold_order


# Manual inspection SOP
perl -lane 'print "$F[0]\t$F[1]\t$F[2]";' < shortread_superscaffold_order.alnstats | grep 'LG6' | head

grep 'LG' shortread_superscaffold_order.segs | perl -lane 'if($F[-1] == 0){next;}else{print $_;}' | less
```

|Group | Scaffold list |
|:--- | :--- |
| LG6 | scaffold_19:-; scaffold_5:-; scaffold_17:-|
|LG3 | scaffold_3:+; scaffold_14:+; scaffold_10:-|
|LG5 | scaffold_9:-; scaffold_18:-; scaffold_25:+|
|LG4 | scaffold_26:-; scaffold_8:+; scaffold_1:+|
|LG7 | scaffold_11:-; scaffold_12:-; scaffold_7:+|
|LG1 | scaffold_4:-; scaffold_6:-|
|LG2 | scaffold_2:+|

OK, The alignments were all over the place so this needs confirmation. Let's use these as a base to see how the bigraph alignment of BAC ends looks. Nevermind that! It is crap! I don't know which to trust, so I'm going to align the ONT reads as an orthogonal dataset to compare against.

```bash
module load minimap2

sbatch -N 1 -n 30 --mem=75000 -p priority -q msn --wrap='minimap2 -t 30 -x map-ont /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/omnic_purge_paligns/salsa_scaffolds/scaffolds_FINAL.fasta /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/RC_ULR.total.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_new_combined_fastqs.fastq.gz > clover_scaffolds_ont_aligns.paf'

sbatch -N 1 -n 2 --mem=55000 -p priority -q msn --wrap="perl ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAnalysisScripts/selectLikelyViralOverhangs.pl clover_scaffolds_ont_aligns.paf clover_scaffolds_ont_overhangs"

sbatch -N 1 -n 2 --mem=10000 -p priority -q msn --wrap='cat /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_new_combined_fastqs.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/RC_ULR.total.fastq | perl convert_fastq_to_fasta.pl; perl ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAnalysisScripts/filterViralOverhangsAndGenerateSeq.pl clover_scaffolds_ont_overhangs.subread.bed temp.fa 200 clover_scaffolds_ont_overhangs.fa; rm temp.fa'

sbatch -N 1 -n 3 --mem=36000 -p priority -q msn --wrap='minimap2 -x map-ont /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/omnic_purge_paligns/salsa_scaffolds/scaffolds_FINAL.fasta clover_scaffolds_ont_overhangs.fa > clover_scaffolds_overhangs.paf'

sbatch -N 1 -n 2 --mem=55000 -p priority -q msn --wrap="perl ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAnalysisScripts/filterOverhangAlignments.pl clover_scaffolds_overhangs.paf clover_scaffolds_ont_aligns.paf clover_scaffolds_filtered_association_table.tab"

perl -e '%h; while(<>){chomp; @s = split(/\t/); shift(@s); @s = sort {$a cmp $b} @s; if($s[0] eq $s[1]){next;} $h{"$s[0]\t$s[1]"} += 1;} foreach $k (sort{$h{$b} <=> $h{$a}} keys(%h)){print "$k\t$h{$k}\n";}' < clover_scaffolds_filtered_association_table.tab
scaffold_10     scaffold_34     107
scaffold_11     scaffold_7      77
scaffold_14     scaffold_34     59
scaffold_2      scaffold_54     57
scaffold_4      scaffold_44     55
scaffold_2      scaffold_32     53
scaffold_3      scaffold_64     45
scaffold_36     scaffold_5      42
scaffold_1      scaffold_141    41
scaffold_1      scaffold_89     35
scaffold_20     scaffold_24     33
scaffold_140    scaffold_18     32
scaffold_41     scaffold_46     28
scaffold_15     scaffold_82     27
scaffold_19     scaffold_42     27
scaffold_13     scaffold_29     22
scaffold_16     scaffold_97     20
scaffold_11     scaffold_20     20
scaffold_48     scaffold_9      20
scaffold_14     scaffold_3      19
scaffold_111    scaffold_9      19
scaffold_20     scaffold_49     17
scaffold_138    scaffold_9      17
scaffold_26     scaffold_35     16
scaffold_1      scaffold_50     16
scaffold_46     scaffold_68     16
scaffold_53     scaffold_82     15


# Trying this with the contigs
sbatch -N 1 -n 30 --mem=75000 -p priority -q msn --wrap='minimap2 -t 30 -x map-ont clover_ccs_themis/fastas/hifiContigs.fa /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/RC_ULR.total.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_new_combined_fastqs.fastq.gz > clover_contigs_ont_aligns.paf'
```

## Ragtag scaffolding

OK, the orthogonal data failed to help me resolve scaffolds. Let's see how RagTag performs. I will try to create superscaffolds using the LG groups from red clover 2, and then using ragtag's Hi-C error correction to try to fix them. Let's start with the contigs from the purged IPA assembly first and then test them against the Salsa scaffolds.

```bash
module load miniconda/3.6 minimap2

conda activate /KEEP/rumen_longread_metagenome_assembly/ragtag

sbatch -N 1 -n 20 -p priority -q msn --mem=105000 --wrap="ragtag.py scaffold -o clover_contigs_rcv2 -t 20 /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/redclover_v2.1.fasta /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_primary_ipa.purged.fa"

sbatch -N 1 -n 20 -p priority -q msn --mem=105000 --wrap="ragtag.py scaffold -o clover_scaffolds_rcv2 -t 20 /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/redclover_v2.1.fasta /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/omnic_purge_paligns/salsa_scaffolds/scaffolds_FINAL.fasta"

# Creating reads for assembly correction
sbatch -N 1 -n 2 --mem=10000 -p priority -q msn --wrap='cat /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano/clover_new_combined_fastqs.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/protist_and_clover/RC_ULR.total.fastq | perl convert_fastq_to_fasta.pl'

mv temp.fa redclover_ont_combined.fasta
gzip redclover_ont_combined.fasta

# Making a comparative alignment plot
conda activate /KEEP/rumen_longread_metagenome_assembly/r

sbatch -N 1 -n 3 --mem=30000 -p priority -q msn --wrap='minimap2 -x asm10 clover_contigs_rcv2/ragtag.scaffold.fasta clover_scaffolds_rcv2/ragtag.scaffold.fasta > ragtag_scaffoldvscontig.paf; Rscript /lustre/project/rumen_longread_metagenome_assembly/binaries/Themis-ASM/scripts/pafDotPlotly.R -i ragtag_scaffoldvscontig.paf -o ragtag_svsc -v -s -l'
```

## Contingency: HifiASM

I am going to try to run HifiASM on the HiFi reads + Hi-C reads to see if we can best resolve the assembly that way. It will make Heng Li happy if it does!

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs

```bash
sbatch -N 1 -n 70 -p priority-mem -q msn-mem --mem=1000000 --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/hifiasm/hifiasm -o redclover.hifi.asm -t 70 --h1 /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L001_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L002_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L003_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L004_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S1_L001_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L001_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L002_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L003_R1_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L004_R1_001.fastq.gz --h2 /lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L001_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L002_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L003_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S10_L004_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S1_L001_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L001_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L002_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L003_R2_001.fastq.gz,/lustre/project/rumen_longread_metagenome_assembly/sequence_data/red_clover_omniC/LIB109270_S9_L004_R2_001.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_200929_165102.Q20.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/clover_ccs/m54337U_201002_155306.Q20.fastq"


perl -ne 'chomp; @F = split(/\t/); if($F[0] eq "S"){$F[2] =~ s/(.{60})/$1\n/g; print ">$F[1]\n$F[2]\n";}' < redclover.hifi.asm.hic.hap1.p_ctg.gfa > redclover.hifi.asm.hic.hap1.p_ctg.fasta

samtools faidx redclover.hifi.asm.hic.hap1.p_ctg.fasta
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i redclover.hifi.asm.hic.hap1.p_ctg.fasta.fai
Total length:   457806091
Num contigs:    1502
N50 length:     233830670
N50 value:      7234713
L50 value:      18
Max:    24557875
Min:    7799

conda activate /KEEP/rumen_longread_metagenome_assembly/r
sbatch -N 1 -n 3 --mem=30000 -p priority -q msn --wrap='minimap2 -x asm10 /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_primary_ipa.purged.fa redclover.hifi.asm.hic.hap1.p_ctg.fasta > hifiasm_ipactg_align.paf; Rscript /lustre/project/rumen_longread_metagenome_assembly/binaries/Themis-ASM/scripts/pafDotPlotly.R -i hifiasm_ipactg_align.paf -o hifiasm_ipactg_align -v -s -l'
```

## Validation EST map

```bash
module load bwa samtools
# Downloading the EST sequences
perl -lane 'print $F[0]; open(IN, "wget -O - http://marker.kazusa.or.jp/Red_clover/sequence/$F[0] |"); open(OUT, "> $F[0].fasta"); $pre = 0; while(<IN>){chomp; if($_ =~ /\<pre\>/){$pre = 1; next;}elsif($_ =~ /\<\/pre\>/){$pre = 0;} if($pre){if($_ =~ /^>/){$_ = ">$F[0]";} print {OUT} $_;}} close IN; close OUT;' < linkage_map.tab

# Now to filter the linkage map to exclude reads that weren't identified
perl -lane 'if( -s "est_linkage/$F[0].fasta"){print $_;}' < est_linkage/linkage_map.tab > est_linkage/linkage_map_filtered.tab

# Moving individual fasta sequence to a new directory to keep things clean
mkdir est_linkage/indiv_fastas
mv est_linkage/*.fasta est_linkage/indiv_fastas/

cat est_linkage/indiv_fastas/*.fasta > est_linkage/linkage_map_filtered.fasta

# converting to a format that gives linkage information
perl -e 'chomp(@ARGV); %data; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $s[2] =~ s/\.\d+$//; $data{$s[0]} = [$s[1], $s[2]]; } close IN; open(IN, " < $ARGV[1]"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; @d = @{$data{$_}}; unshift(@d, $_); print ">" . join(".", @d) . "\n";}elsif($_ =~ /^\s+$/){next;}else{print $_;}} close IN;' est_linkage/linkage_map_filtered.tab est_linkage/linkage_map_filtered.fasta > est_linkage/linkage_map_filtered.fmt.fasta

# Running the order script
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a clover_ccs_themis/fastas/hifiScaffolds.fa -p est_linkage/linkage_map_filtered.fmt.fasta -o est_linkage/hifi_scaffolds

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a clover_ccs_themis/clover_scaffolds_rcv2/ragtag.scaffold.fasta -p est_linkage/linkage_map_filtered.fmt.fasta -o est_linkage/ragtag_scaffolds

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a clover_ccs_themis/clover_contigs_rcv2/ragtag.scaffold.fasta -p est_linkage/linkage_map_filtered.fmt.fasta -o est_linkage/ragtag_contigs

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a clover_ccs_themis/fastas/hifiContigs.fa -p est_linkage/linkage_map_filtered.fmt.fasta -o est_linkage/hifi_contigs

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a clover_ccs_themis/fastas/shortRead.fa -p est_linkage/linkage_map_filtered.fmt.fasta -o est_linkage/short_read

# Identifying best assembly based on consensus statistics
for i in est_linkage/*.stats; do echo $i; perl -lane 'if($F[0] eq 'Mapping'){next;}else{print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < $i; done

## If we base things off of consensus incorporation of sequence that has at least 10 mapping markers in succession, it looks like the ragtag contigs are best. 
## LG2: LG2 scaffold + ctg.000054F_1_1  <- actually didn't pan out might be heterozygous alignment?
## LG5: LG5 scaffold + ctg.000021F_1_1 (<- also interior to the scaffold) + ctg.000066F_1_1 + ctg.000085F_1_1
```

So the only modification I should make appears to be:

LG5: ctg.000066F_1_1:-, ctg.000085F_1_1:+, LG5

Let's make that modification in a new agp file and generate what I feel is a suitable penultimate assembly. Then, maybe I can run ragtag with hi-c? Let's just make the fasta file and see what the collaborators think first.

```bash
perl -lane '$t = $F[0]; $t =~ s/Tp57577_TGAC_v2_//; $t =~ s/_RagTag//; print "$F[0]\t1\t$F[1]\t$t\t1\t+";' < clover_ccs_themis/clover_contigs_rcv2/ragtag.scaffold.fasta.fai > est_linkage/ragtag_contigs.plan.bed

vim est_linkage/ragtag_contigs.plan.bed

module load java
java -Xmx15g -jar ~/rumen_longread_metagenome_assembly/binaries/CombineFasta/store/CombineFasta.jar agp2fasta -f clover_ccs_themis/clover_contigs_rcv2/ragtag.scaffold.fasta -b est_linkage/ragtag_contigs.plan.bed -o est_linkage/ARS_DFRC.1.fasta
samtools faidx est_linkage/ARS_DFRC.1.fasta

# I'm not too happy because only 260 Megs are in superscaffolds!
```


### Network analysis

OK, I want to see if I can consolidate all of the linkage information that I have into networks and analyze the topology. I will use the EST data to assess contig/scaffold assignment in general (to linkage groups) but then create edges based on each method. I will then try to consolidate all networks and weight the edges based on discovery. 

Let's first organize all of the data I have:

#### Data types
* BAC_ends: bac_ends/hifi_bacends_discord.bigraph
* EST_linkages: est_linkage/hifi_scaffolds.tab
* ONT longreads: clover_ccs_themis/clover_scaffolds_filtered_association_table.tab


I just checked assembly sizes, and the Hifi scaffolds (from Salsa) have more bases assigned to linkage groups than the base contigs. I will try to use the scaffolds in the first round of network analysis and omit the Hi-C linkage data (as it will be highly biased!). 

All data files will need to be converted to a bi-graph format.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs

```bash
mkdir graph_comps

# Removing duplicate entries from the bacend analysis
perl -e '$c = 0; while(<>){if($c % 2 == 0){print $_;} $c++;}' < bac_ends/hifi_bacends_discord.bigraph > graph_comps/hifiScaffolds_bacends.bigraph

# Taking raw associations between ESTS into the bigraph
perl -e '%data; $lastc = "None"; $lasts = "None"; while(<>){chomp; @s = split(/\t/); if($s[4] ne $lastc){$lastc = $s[4];}elsif($s[1] eq "*"){next;}else{@cats = sort { $a cmp $b} ($lasts, $s[1]); $data{$cats[0]}->{$cats[1]} += 1;} $lasts = $s[1];} foreach $first (keys(%data)){foreach $second (keys(%{$data{$first}})){print "$first\t$second\t" . $data{$first}->{$second} . "\n";}}' < est_linkage/hifi_scaffolds.tab > graph_comps/hifiScaffolds_estlinks.bigraph

# converting the ont read overlaps into a bigraph
perl -e '%data; while(<>){chomp; @s = split(/\t/); @y = sort{$a cmp $b} ($s[1], $s[2]); $data{$y[0]}->{$y[1]} += 1;} foreach $first (keys(%data)){foreach $second (keys(%{$data{$first}})){print "$first\t$second\t" . $data{$first}->{$second} . "\n";}}' < clover_ccs_themis/clover_scaffolds_filtered_association_table.tab > graph_comps/hifiScaffolds_ont.bigraph

# I need attribute data for each scaffold including LG assignment based on EST maps. I'll keep it simple and just consolidate scaffold length and LG assignment for now
perl -e 'chomp(@ARGV); @list = ("LG", "Length"); %data; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[0]}->{"Length"} = $s[1]; $data{$s[0]}->{"LG"} = "Unassigned";} close IN; open(IN, "< $ARGV[0]"); <IN>; while(<IN>){@s = split(/\s+/); $lg = $s[1] . ""; for($x =0; $x < 3; $x++){shift(@s);} foreach my $d (@s){@bsegs = split(/:/, $d); if($bsegs[0] eq "*"){next;} $data{$bsegs[0]}->{"LG"} = $lg;} } close IN; print "Scaffold\t" . join("\t", @list) . "\n"; foreach my $k (keys(%data)){print "$k"; foreach my $l (@list){print "\t" . $data{$k}->{$l}} print "\n";}' est_linkage/hifi_scaffolds.stats clover_ccs_themis/fastas/hifiScaffolds.fa.fai > graph_comps/hifiScaffolds_attributes.tab


```

## Finalizing scaffolds

OK, my network analysis resulted in a strategy for finishing this up. I need to figure out how several of the unassigned scaffolds fit into the mix, but once I do, I can finish everything up.

I am going to use the Salsa HiFi scaffolds to start as my base. I will then take all EST assigned scaffolds and place them in the best order possible. Then I will add 8 otherwise unassigned scaffolds to the linkage groups based on the network analysis. 

```bash
module load bwa samtools
# First, generate a consensus for each scaffold and a range for values
perl -e '%data; while(<>){chomp; @s = split(/\t/); if($s[1] eq "*"){next;}  $data{$s[1]}->{$s[4]} += 1;} foreach $s (sort {$a cmp $b} keys(%data)){foreach $k (sort {$data{$s}->{$b} <=> $data{$s}->{$a}} keys(%{$data{$s}})){print "$s\t$k\t" . $data{$s}->{$k} . "\n"; last;}}' < est_linkage/hifi_scaffolds.tab > est_linkage/hifi_scaffolds_consensus.tab

perl -e 'chomp(@ARGV); %cons; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $cons{$s[0]} = $s[1];} close IN; %range; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[4] ne $cons{$s[1]}){next;} if(!exists($range{$s[1]})){$range{$s[1]} = [$cons{$s[1]}, $s[5], $s[5], $s[2], $s[2]];}else{$range{$s[1]}->[2] = $s[5]; $range{$s[1]}->[4] = $s[2];}} close IN; foreach my $k (sort {$range{$a}->[0] cmp $range{$b}->[0] || $range{$a}->[1] <=> $range{$b}->[1]} keys(%range)){print "$k\t" . join("\t", @{$range{$k}}) . "\n";}' est_linkage/hifi_scaffolds_consensus.tab est_linkage/hifi_scaffolds.tab > est_linkage/hifi_scaffolds.prelim.plan.tab

# I'm going to use the minimum placing for each scaffold in the meantime. Then I will pad entries by 1000 bp in the final scaffolds
perl -e 'chomp(@ARGV); %lens; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); $c = 1; $last = "None"; while(<IN>){chomp; @s = split(/\t/); if($last ne $s[1]){$last = $s[1]; $c = 1;} $end = $lens{$s[0]}; $orient = ($s[4] - $s[5] <= 0)? "+" : "-"; print "$s[0]\t1\t$end\t$s[1]\t$c\t$orient\n"; $c++;} close IN;' clover_ccs_themis/fastas/hifiScaffolds.fa.fai est_linkage/hifi_scaffolds.prelim.plan.tab > est_linkage/hifi_scaffolds.prelim.plan.bed

mkdir placed_scaffs
cp est_linkage/hifi_scaffolds.prelim.plan.bed placed_scaffs/hifi_scaffolds.plan.bed

# I only added in four scaffolds that were 100kb+ for a total addition of around 1 Mbp to the final assembly
perl -e '%data; while(<>){chomp; @s = split(/\t/); $data{$s[3]} += $s[2];} foreach $k (keys(%data)){print "$k\t$data{$k}\n";}' < placed_scaffs/hifi_scaffolds.plan.bed
LG6     46133763
LG2     74472114
LG7     59238277
LG5     57895441
LG1     53631259
LG4     60492793
LG3     54406022

# I need to now add in the remainder at the end of the assembly
perl -lane 'print $F[0];' < placed_scaffs/hifi_scaffolds.plan.bed > placed_scaffs/hifi_scaffolds.plan.scaff.list
perl -lane 'print $F[0];' < clover_ccs_themis/fastas/hifiScaffolds.fa.fai > hifiScaffolds.scaff.list
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl hifiScaffolds.scaff.list placed_scaffs/hifi_scaffolds.plan.scaff.list
File Number 1: hifiScaffolds.scaff.list
File Number 2: placed_scaffs/hifi_scaffolds.plan.scaff.list
Set     Count
1       81
1;2     62

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1 hifiScaffolds.scaff.list placed_scaffs/hifi_scaffolds.plan.scaff.list > hifiScaffolds.missing.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f clover_ccs_themis/fastas/hifiScaffolds.fa.fai -l hifiScaffolds.missing.list -d '\t' -c 0 | perl -lane 'print "$F[0]\t1\t$F[1]\t$F[0]\t1\t+";' > hifiScaffolds.missing.bed
cat placed_scaffs/hifi_scaffolds.plan.bed hifiScaffolds.missing.bed > placed_scaffs/hifi_scaffolds.plan.final.bed

# Now let's make a fasta!
module load java
java -Xmx15g -jar ~/rumen_longread_metagenome_assembly/binaries/CombineFasta/store/CombineFasta.jar agp2fasta -f clover_ccs_themis/fastas/hifiScaffolds.fa -b placed_scaffs/hifi_scaffolds.plan.final.bed -o placed_scaffs/ARS_DFRC_1.1.fasta -i 1000

samtools faidx placed_scaffs/ARS_DFRC_1.1.fasta
```

#### New themis run

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis

```bash
module load python_3/3.6.6 miniconda/3.6 "java/11.0.2" bedtools

PATH=$PATH:/software/7/apps/merqury/1.1/:/software/7/apps/meryl/1.0/Linux-amd64/bin/
MERQURY=/software/7/apps/merqury/1.1

ln -s /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_ccs_themis/busco_downloads busco_downloads
mkdir logs

sbatch -N 1 -n 2 --mem=6000 -p priority -q msn --wrap='python3 ~/rumen_longread_metagenome_assembly/binaries/Themis-ASM/themisASM.py -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/redclover_v2.1.fasta -n shortRead -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/clover_primary_ipa.purged.fa -n hifiContigs -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/placed_scaffs/ARS_DFRC_1.1.fasta -n ARSDFRC1 -a /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/ncbi-genomes-2021-08-04/GCA_000219495.2_MedtrA17_4.0_genomic.fna -n MedTr4 -b eudicots_odb10 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L001-ds.ac229a7a209e4f809c4ac2bda32e46fe/Hen17Cv_S1_L001_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L001-ds.ac229a7a209e4f809c4ac2bda32e46fe/Hen17Cv_S1_L001_R2_001.fastq.gz -s lane1 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L002-ds.e4d313d6db7a45ec88dba60164d0d717/Hen17Cv_S1_L002_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L002-ds.e4d313d6db7a45ec88dba60164d0d717/Hen17Cv_S1_L002_R2_001.fastq.gz -s lane2 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L003-ds.20e6182d55ac40e2b889f55011234126/Hen17Cv_S1_L003_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L003-ds.20e6182d55ac40e2b889f55011234126/Hen17Cv_S1_L003_R2_001.fastq.gz -s lane3 -f /lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L004-ds.8b23d704e97542d982a7416ef5ce3c8a/Hen17Cv_S1_L004_R1_001.fastq.gz,/lustre/project/forage_assemblies/sequence_data/CloverGenome-137246120/FASTQ_Generation_2019-06-16_03_19_42Z-188429241/Hen17Cv_L004-ds.8b23d704e97542d982a7416ef5ce3c8a/Hen17Cv_S1_L004_R2_001.fastq.gz -s lane4 -c "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -j 100 --resume'


module load busco4

busco_configurator.py /software/7/apps/busco4/4.0.2/config/config.ini /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_config.ini
export BUSCO_CONFIG_FILE=/lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_config.ini

cp -Rp /software/7/apps/augustus/3.3.2/config/ /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/AUGUSTUS_CONFIG
export AUGUSTUS_CONFIG_PATH=/lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/AUGUSTUS_CONFIG

conda activate .snakemake/conda/a6cbb216
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="busco -m genome -i fastas/ARSDFRC1.fa -f --offline -l /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_downloads/lineages/eudicots_odb10 -c 70 -o busco_arsdfrc_rerun"
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="busco -m genome -i fastas/hifiContigs.fa -f --offline -l /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_downloads/lineages/eudicots_odb10 -c 70 -o busco_hificontigs_rerun"
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="busco -m genome -i fastas/MedTr4.fa -f --offline -l /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_downloads/lineages/eudicots_odb10 -c 70 -o busco_medtr4_rerun"
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="busco -m genome -i fastas/shortRead.fa -f --offline -l /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/manuscript_themis/busco_downloads/lineages/eudicots_odb10 -c 70 -o busco_shortread_rerun"
```

#### Star annotation

```bash
module load star/2.7.9a

# creating the index
sbatch -n 20 -N 1 --mem=36000 -p priority -q msn --wrap="mkdir fastas/star_ars; STAR --runThreadN 20 --runMode genomeGenerate --genomeDir fastas/star_ars --genomeFastaFiles fastas/ARSDFRC1.fa"

sbatch -n 20 -N 1 --mem=36000 -p priority -q msn --wrap="mkdir fastas/star_ars; STAR --runThreadN 20 --genomeDir fastas/star_ars --readFilesIn /lustre/project/rumen_longread_metagenome_assembly/plants/red_clover_IsoSeq/collapse_isoforms.fastq --outFileNamePrefix star_ars_dfrc"

# The star aligner complains about quality score string lengths...
module load minimap2

sbatch -n 20 -N 1 --mem=36000 -p priority -q msn --wrap='minimap2 -x splice fastas/ARSDFRC1.fa /lustre/project/rumen_longread_metagenome_assembly/plants/red_clover_IsoSeq/hq_transcripts.fastq > ars_dfrc_hq_transcripts.paf'

perl -lane 'print "$F[5]\t$F[7]\t$F[8]\t$F[0]";' < ars_dfrc_hq_transcripts.paf | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o distinct -delim ';' | head
```


#### Creating circos diagram from PAF

> /mnt/c/SharedFolders/sequencing_projects/red_clover/hifi_assemblies 

```bash
python3 /mnt/c/SharedFolders/test_software/CircosAlignmentPlotter/p2c.py mapshortRead_ARSDFRC1.paf shortRead.fa.fai ARSDFRC1.fa.fai targets.bed circos_TGAC_ARSDFRC

cd circos_TGAC_ARSDFRC

conda activate circos 

circos
```

Note: I had to remove a ton of link data because the script failed to filter links that were not present in the "targets.bed" file. I also had to modify the configuration file to remove extraneous stuff, and I also had to manipulate the config files to make the picture look a bit better. Still, not a bad start!


#### NCBI region check

```bash
module load samtools minimap2 bedtools/2.25.0 java

perl -e 'while(<>){chomp; @s = split(/\t/); $type = pop(@s); $chr = shift(@s); shift(@s); @bsegs = split(/,/, $s[0]); foreach $k (@bsegs){@c = split(/\.\./, $k); $len = $c[1] - $c[0]; print "$chr\t$c[0]\t$c[1]\t$len\t$type\n";}}' < ncbi_regions.txt > ncbi_regions.bed

sbatch -N 1 -n 30 -p priority -q msn -t 1-0 --mem=50000 --wrap='minimap2 -ax asm20 -t 20 ARS_RC_1.1.fasta ../m54337U_200929_165102.Q20.fastq | samtools sort -T m54 -o m54337U_200929_165102.Q20.bam -@ 10 -'

perl -lane 'open(IN, "samtools view m54337U_200929_165102.Q20.bam $F[0]:$F[1]-$F[2] |"); $h = 0; $c = 0; while(<IN>){chomp; @s = split(/\t/); if($s[5] =~ /H/){$h += 1;} $c++;} close IN; $r = $h/$c; print join("\t", @F) . "\t$h\t$c\t$r";' < ncbi_regions.bed > ncbi_regions.hardclipped.tab

# Now to mask all sequence that was mitochondrial and clearly misassembled
perl -lane 'if($F[5] > 0){print $F[3];}' < ncbi_regions.hardclipped.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   47
Sum:    23044
Minimum 119
Maximum 8907
Average 490.297872
Median  196
Standard Deviation      1328.009444
Mode(Highest Distributed Value) 119

perl -lane 'if($F[5] > 0){print "$F[0]\t$F[1]\t$F[2]";}' < ncbi_regions.hardclipped.tab > ncbi_regions.tomask.bed
bedtools maskfasta -fi ARS_RC_1.1.fasta -fo ARS_RC_1.1.v2.fasta -bed ncbi_regions.tomask.bed

perl -lane 'if($F[5] < 1){print "$F[0]:$F[1]-$F[2]";}' < ncbi_regions.hardclipped.tab > ncbi_regions_numt.tab

perl -lane 'print "$F[0]\t1\t$F[1]";' < ARS_RC_1.1.fasta.fai > ARS_RC_1.1.fasta.lens.bed
java -Xmx15g -jar ~/rumen_longread_metagenome_assembly/binaries/CombineFasta/store/CombineFasta.jar agp2fasta -f  ARS_RC_1.1.fasta -b ncbi_regions.tomask.newplan.bed -o ARS_RC_1.1.v2.fasta -i 50

perl -lane 'system("samtools faidx ARS_RC_1.1.fasta $F[0]");' < ncbi_regions_numt.tab > ncbi_regions_numt.fa
```

```bash
minimap2 -x asm5 -N 50 --secondary=no /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/placed_scaffs/ARS_RC_1.1.v2.fasta final_IsoSeq_analysis/hq_transcripts.fasta > ARS_RC_1.1.v2.transcript.paf

minimap2 -ax splice --cs /lustre/project/rumen_longread_metagenome_assembly/assemblies/clover_ccs/placed_scaffs/ARS_RC_1.1.v2.fasta final_IsoSeq_analysis/hq_transcripts.fasta | samtools sort -T temp -o ARS_RC_1.1.v2.transcript.bam -

bedtools bamtobed -bed12 -i ARS_RC_1.1.v2.transcript.bam > ARS_RC_1.1.v2.transcript.bed
/lustre/project/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/bedToGenePred ARS_RC_1.1.v2.transcript.bed ARS_RC_1.1.v2.transcript.genepred
/lustre/project/rumen_longread_metagenome_assembly/binaries/kentUtils/bin/linux.x86_64/genePredToGtf 'file' ARS_RC_1.1.v2.transcript.genepred ARS_RC_1.1.transcript.gtf
```