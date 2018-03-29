# Rumen microbiome technology bias check
---
*8/9/2017*

These are my notes and commands for testing the bias of different sequencing platforms in the collection and determining if there are major limitations to each.

## Table of Contents
* Mash Profile comparisons
* Mash Sketch Generation of OutGroups
* [Illumina downsampling test](#downsampling)
* [Mash NMDS for Dataset Comparison](#mash_nmdsv1)

## Mash profile comparisons

I am going to test for kmer cardinality differences in the methods using MASH. My concern is that certain kmers may be underrepresented in the different methods and this is the best method for finding out!

Some oddities of MASH: it doesn't read gzipped files and each separate file is considered an "island" in each sketch comparison. I need to merge files to keep them in the same profile for comparison. The great news is that the sketches I create can be transferred and used on other systems, making them great portable compressed profiles.

#### Fastq merger and preparation

```bash
cat *.fastq > nanopore_yu_and_morrison_3.combined.fq
for i in illuminaR*/*.fastq.gz; do echo $i; gunzip $i; done

# I wrote a python script to filter the reads
python3 ../../../programs_source/python_toolchain/sequenceData/filterNextseqFastqFiles.py -f YMPrepCannula_S1_L001_R1_001.fastq -r YMPrepCannula_S1_L001_R2_001.fastq -o YMRewriteFilterGbase
Processed 29316700 reads and identified 884938 single G-base and 3990662 both G-base artifacts

# Automating the rest
for i in illuminaR1/YMPrepCannula_S1_L002_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L002_R2_001.fastq illuminaR1/YMPrepCannula_S1_L003_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L003_R2_001.fastq illuminaR1/YMPrepCannula_S1_L004_R1_001.fastq,illuminaR1/YMPrepCannula_S1_L004_R2_001.fastq illuminaR2/YMPrepCannula_S1_L001_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L001_R2_001.fastq illuminaR2/YMPrepCannula_S1_L002_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L002_R2_001.fastq illuminaR2/YMPrepCannula_S1_L003_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L003_R2_001.fastq illuminaR2/YMPrepCannula_S1_L004_R1_001.fastq,illuminaR2/YMPrepCannula_S1_L004_R2_001.fastq; do file1=`echo $i | cut -d',' -f1`; file2=`echo $i | cut -d',' -f2`; folder=`echo $i | cut -d'/' -f1`; echo $folder; lane=`basename $file1 | cut -d'_' -f3`; python3 ../../programs_source/python_toolchain/sequenceData/filterNextseqFastqFiles.py -f $file1 -r $file2 -o ${folder}"/YMRewriteFilterGbase."${lane} ; done
illuminaR1
Processed 52887277 reads and identified 1863215 single G-base and 26864810 both G-base artifacts
illuminaR1
Processed 26327191 reads and identified 572981 single G-base and 1135473 both G-base artifacts
illuminaR1
Processed 28095512 reads and identified 793048 single G-base and 3164439 both G-base artifacts
illuminaR2
Processed 63812154 reads and identified 808176 single G-base and 34846 both G-base artifacts
illuminaR2
Processed 63221220 reads and identified 892873 single G-base and 85939 both G-base artifacts
illuminaR2
Processed 63346742 reads and identified 777543 single G-base and 17110 both G-base artifacts
illuminaR2
Processed 62625319 reads and identified 818123 single G-base and 34208 both G-base artifacts
```

#### Mash sketch generation

I want to create Mash sketches that conform with Serge's metagenomics estimates and are comparable to the sketches I've generated [previously](https://github.com/njdbickhart/labnotes/blob/master/project_notes/metagenomics/Notes_20161219_dbick_metagenomics_software_test_notes.md#mash). So I will use the same settings as before.

```bash
for i in *.fastq; do name=`echo $i | cut -d'.' -f1`; echo $name; mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name $i; done

for i in cheryl/*/*.subreads.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; samtools fastq $i > $name.fq; mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name $name.fq; rm $name.fq; done

# Sketching all Illumina reads from stdin
for i in illuminaR*/YMRewriteFilterGbase*fq; do cat $i; done | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o illumina_filtered_non-interleaved -

# Sketching pacbio reads separately from bam
for i in *.bam; do name=`echo $i | cut -d'.' -f1`; echo $name; samtools fastq $i | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name - ; done

for i in ./*.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; samtools fastq $i | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o $name - ; done

# It's also only fair that I create a combined sketch from each facility's instrument for comparison
# Tim's first (RSII and Sequel)
cat LIB*.fastq | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o pacbio_rsii_summary -
for i in ./*/m54033*.bam; do samtools fastq $i; done | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o pacbio_sequel_tim_summary -

# Cheryl's second (Sequel and CCS'ed reads)
find ./ -name '*.bam' | grep -v 'scraps' | grep -v 'm54033' | xargs -I{} samtools fastq {} | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o pacbio_sequel_cheryl_summary -
cat ./*/*.fastq | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o ../pacbio_ccs_cheryl_summary -

# Nanopore read sketching is relatively straightforward; all of my generated data will be sketched together
for i in `find ./nanopore -name "*.fastq"`; do cat $i ; done | mash sketch -p 3 -k 21 -s 10000 -r -m 2 -o nanopore_yu_morrison_fastq -

```


## Mash sketch generation of outgroups

I want to include the cattle manure sample studies and the [Hess et al](https://www.ncbi.nlm.nih.gov/pubmed/21273488) dataset as "outgroups" to compare with our new dataset MASH profiles.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/

```bash
# The Hess et al. survey
sbatch --nodes=1 --ntasks-per-node=6 --mem=10000 --partition=assemble3 --wrap='for i in SRR094166_1.fastq.gz SRR094166_2.fastq.gz SRR094403_1.fastq.gz SRR094403_2.fastq.gz SRR094405_1.fastq.gz SRR094405_2.fastq.gz SRR094415_1.fastq.gz SRR094415_2.fastq.gz SRR094416_1.fastq.gz SRR094416_2.fastq.gz SRR094417_1.fastq.gz SRR094417_2.fastq.gz SRR094418_1.fastq.gz SRR094418_2.fastq.gz SRR094419_1.fastq.gz SRR094419_2.fastq.gz SRR094424_1.fastq.gz SRR094424_2.fastq.gz SRR094427_1.fastq.gz SRR094427_2.fastq.gz SRR094428_1.fastq.gz SRR094428_2.fastq.gz SRR094429_1.fastq.gz SRR094437_1.fastq.gz SRR094437_2.fastq.gz SRR094926_1.fastq.gz SRR094926_2.fastq.gz; do gunzip -c datasources/$i ; done | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -s 10000 -r -m 2 -o illumina_hess_combined - '

# Manure sample survey
sbatch --nodes=1 --ntasks-per-node=6 --mem=10000 --partition=assemble3 --wrap='for i in SRR2329878_1.fastq.gz SRR2329878_2.fastq.gz SRR2329910_1.fastq.gz SRR2329910_2.fastq.gz SRR2329939_1.fastq.gz SRR2329939_2.fastq.gz SRR2329962_1.fastq.gz SRR2329962_2.fastq.gz; do gunzip -c datasources/$i; done | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -s 10000 -r -m 2 -o illumina_manure_combined - '

# I transferred all of the files to my linux VM for distance comparison
```

## Mash distance estimation and concatenation

OK, now to generate the distance matrix and to create ordination plots.

> pwd: /home/dbickhart/share/metagenomics/pilot_project/mash_profiles

```bash
ls illumina_filtered_non-interleaved.msh illumina_hess_combined.msh illumina_manure_combined.msh FNFAE24738.msh nanopore_yu_morrison_fastq.msh pacbio_ccs_cheryl_summary.msh pacbio_rsii_summary.msh pacbio_sequel_cheryl_summary.msh pacbio_sequel_tim_summary.msh > combined_summary_sketches.list

mash paste -l combined_summary_sketches combined_summary_sketches.list

for i in `cat combined_summary_sketches.list`; do echo $i; mash dist -p 3 -t combined_summary_sketches.msh $i > $i.combined.dist; done
mkdir combined_dist
mv *.dist ./combined_dist/

# Because most of my sketches were done based on standard input, I have to be creative in how I paste the distance matricies together
# You can see the file IDs in the mash "info" field. 
perl -lane 'open(IN, "< combined_dist/$F[0].combined.dist"); $h = <IN>; $d = <IN>; chomp $d; @dsegs = split(/\t/, $d); $dsegs[0] = $F[0]; print join("\t", @dsegs);' < combined_summary_sketches.list > combined_dist/combined_distance.matrix

# Now to read it in and do some plotting in R
```

My goal here is a rapid NMDS with the following simplistic categories:

* Illumina
* Nanopore
* PacBio
* Outgroup (combined manure profile)

```R
library(MASS)
library(vegan)

# just in case I need them
library(dplyr)
library(tidyr)

mash_data <- read.delim("combined_dist/combined_distance.matrix", header=FALSE)
entry_names <- mash_data[,1]
mash_data.format <- select(mash_data, -V1)
rownames(mash_data.format) <- entry_names
colnames(mash_data.format) <- entry_names

ord <- metaMDS(as.dist(mash_data.format), trymax=100)
# Only took 20 iterations

# Base Vegan plotting polymorphic functions did not work on the data. Had to hack out the NMDS coords
plot(scores(ord), col=c("red", "blue", "blue", "green", "red", "purple", "purple", "purple", "purple"), pch= c(15, 15, 16, 15, 16, 16, 17, 15, 15))
legend("bottomleft", c("UKNano", "USIllum", "Hess", "Manure", "USNano", "PBCCS", "PBRSII", "PBCheryl", "PBTim"), col=c("red", "blue", "blue", "green", "red", "purple", "purple", "purple", "purple"), pch= c(15, 15, 16, 15, 16, 16, 17, 15, 15))
dev.copy2pdf(file="combined_profile_nmds_vegan.pdf", useDingbats=FALSE)
```
<a name="downsampling"></a>
## Illumina downsampling test

I want to see how reproducible the mash profile is for different proportions of the Illumina data. Hopefully, I can get a good reproducibility curve from a progressive downsampling that shows a cutoff for identity.

> pwd: /home/dbickhart/share/metagenomics/pilot_project/downsample

```bash
# a simple test on one fastq file
perl ~/share/programs_source/Perl/perl_toolchain/metagenomics_scripts/downsampleIlluminaReads.pl ../illuminaR1/YMPrepCannula_S1_L001_R1_001.fastq.gz

# Now for the whole shebang
perl ~/share/programs_source/Perl/perl_toolchain/metagenomics_scripts/downsampleIlluminaReads.pl ../illuminaR1/YMPrepCannula_S1_L001_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L001_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L002_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L002_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L003_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L003_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L004_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L004_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L001_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L001_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L002_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L002_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L003_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L003_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L004_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L004_R2_001.fastq.gz

# Now let's generate the distance matrix so that we can estimate the percent similarity to the original dataset
for i in downsample*.msh; do name=`echo $i | cut -d'.' -f1,2`; echo $name; mash dist -t ../mash_profiles/illumina_filtered_non-interleaved.msh $i > $name.dist; done

# The 20% only got down to 3% dissimilarity. I'm worried that the majority of mash sketches are picking up the same species and that abundance is killing things
# Testing out with more gradients (5% intervals)
perl ~/share/programs_source/Perl/perl_toolchain/metagenomics_scripts/downsampleIlluminaReads.pl ../illuminaR1/YMPrepCannula_S1_L001_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L001_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L002_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L002_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L003_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L003_R2_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L004_R1_001.fastq.gz ../illuminaR1/YMPrepCannula_S1_L004_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L001_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L001_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L002_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L002_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L003_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L003_R2_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L004_R1_001.fastq.gz ../illuminaR2/YMPrepCannula_S1_L004_R2_001.fastq.gz

for i in downsample*.msh; do name=`echo $i | cut -d'.' -f1,2`; echo $name; mash dist -t downsampleSketch_1.0.msh $i > $name.dist; done
# This time, the 5% dataset only reached 9% dissimilarity

# Creating a table for loading to R
perl -e '@f = `ls *.dist`; chomp(@f); for($x = scalar(@f) - 1; $x >= 0; $x--){open($IN, "< $f[$x]"); <$IN>; $v = <$IN>; chomp $v; $v =~ s/^\s+//g; ($b) = $f[$x] =~ m/downsampleSketch_(.+)\.dist/; $b *= 100; $v = 1 - $v; print "$b\t$v\n"; close $IN;}' > downsampleSimilarityMatrix.tab
```

Now I'm going to load it into R and generate a line plot.

```R
simMatrix <- read.delim(file="downsampleSimilarityMatrix.tab", header=FALSE)
colnames(simMatrix) <- c("sampling", "similarity")
library(ggplot2)

ggplot(data=simMatrix, aes(x = sampling, y = similarity)) + geom_line() + geom_point() + expand_limits(y=0.875) + xlab(label="Sampling percentage") + ylab(label="Similarity propotion (1 - Mash distance)")

# I'd like to calculate the derivative as well to see if things are slowing down as we reach higher sampling
```

I just ran Tim's PCR free library and will be testing it out against the Nextera prep.

```bash
for i in illuminaR3PCRFree/*.fastq.gz; do echo -n "$i "; done; echo
perl ~/share/programs_source/Perl/perl_toolchain/metagenomics_scripts/downsampleIlluminaReads.pl illuminaR3PCRFree/YMPrepCannula_S1_L001_R1_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L001_R2_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L002_R1_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L002_R2_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L003_R1_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L003_R2_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L004_R1_001.fastq.gz illuminaR3PCRFree/YMPrepCannula_S1_L004_R2_001.fastq.gz

for i in downsample*.msh; do name=`echo $i | cut -d'.' -f1,2`; echo $name; mash dist -t downsampleSketch_1.0.msh $i > $name.dist; done
perl -e '@f = `ls *.dist`; chomp(@f); for($x = scalar(@f) - 1; $x >= 0; $x--){open($IN, "< $f[$x]"); <$IN>; $v = <$IN>; chomp $v; $v =~ s/^\s+//g; ($b) = $f[$x] =~ m/downsampleSketch_(.+)\.dist/; $b *= 100; $v = 1 - $v; print "$b\t$v\n"; close $IN;}' > downsampleSimilarityMatrixR3.tab
```

It's the same profile -- low distances between the 5% and the 100% samples. I wonder if we have a substantial bias here? I will test it out on some of the assemblies other groups have prepared.

```bash
for i in ../../hungate/assembly_fasta/*.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; perl -e 'chomp(@ARGV); open(IN, "gunzip -c $ARGV[0] |"); while(<IN>){chomp; if($_ =~ /\>(.+)/){print ">$ARGV[1]\_$1\n";}else{print "$_\n";}}' $i $name >>  hungate_concatenated.fa; done

samtools faidx hungate_concatenated.fa
bwa index hungate_concatenated.fa

```

## GenomeScope test

Using a python script to download the tarballs of the reads from my Gdrive, I was able to get the data I needed to start the characterization. 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina/illumina

```bash
ls *.fastq.gz | xargs -n 1 echo gunzip -c > generators
sbatch --mem=200000 --ntasks-per-node=20 --nodes=1 -p assemble3 --wrap="jellyfish count -C -m 21 -s 1000000000 -t 10 -o YMPrepCannula_r1_21mer.jf -g generators -G 2"

# Increased "high bound" of histos to 1,000,000.
# If you don't increase this value, then all kmers with counts > 10,000 get slotted into the same bin
sbatch --mem=30000 --ntasks-per-node=10 --nodes=1 -p assemble3 --wrap="jellyfish histo -t 10 --o YMPrepCannula_r1_21mer.histo -h 1000000 YMPrepCannula_r1_21mer.jf"

# Downloading the other files from run3: 
perl -e '%h = ("0BxbRXPCzWa5-Mkd6WmtLTDF2Sk0" => "YMPrepCannula_run3_L1_R1.fastq.gz", "0BxbRXPCzWa5-SE91UkVqU1hVeE0" => "YMPrepCannula_run3_L2_R1.fastq.gz", "0BxbRXPCzWa5-QVRCVFVRZXM1SkU" => "YMPrepCannula_run3_L2_R2.fastq.gz", "0BxbRXPCzWa5-WWJPcUJtcVY0eXc" => "YMPrepCannula_run3_L3_R1.fastq.gz", "0BxbRXPCzWa5-TlVoajNBN3l2TWs" => "YMPrepCannula_run3_L3_R2.fastq.gz", "0BxbRXPCzWa5-UHNocU8zR094d2s" => "YMPrepCannula_run3_L4_R1.fastq.gz", "0BxbRXPCzWa5-UWNPSmI3VWlIV1E" => "YMPrepCannula_run3_L4_R2.fastq.gz"); foreach my $k (keys(%h)){print "$h{$k}\n"; system("python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py $k /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina/$h{$k}")}'
```

I think that the GenomeScope model will break down at the high frequency of kmers in the datasets. Let's see if I can align error corrected Illumina reads to the data to identify variant sites within the error corrected reads.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/error_corrected_reads

```bash
sbatch --mem=15000 --nodes=1 --ntasks-per-node=1 --wrap="bwa index rumen_nanopore_corrected.fasta"
```

OK, that won't work because the read names are all overlapping due to the shear number of reads. Let's try with the original assemblies.

```bash
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 0BxbRXPCzWa5-MGNXaWhlRXgwR00 rumen_pacbio_pilot_asm.fasta.gz
gunzip rumen_pacbio_pilot_asm.fasta.gz
sbatch --mem=15000 --nodes=1 --ntasks-per-node=1 --wrap="bwa index rumen_pacbio_pilot_asm.fasta"

python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 0BxbRXPCzWa5-Q1VMa1RLOTNhWm8 rumen_nanopore_pilot_asm.fasta.gz
gunzip rumen_nanopore_pilot_asm.fasta.gz
sbatch --mem=15000 --nodes=1 --ntasks-per-node=1 --wrap="bwa index rumen_nanopore_pilot_asm.fasta"

perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b pacbio -t illumina_reads.tab -f rumen_pacbio_pilot_asm.fasta -m
perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b nanopore -t illumina_reads.tab -f rumen_nanopore_pilot_asm.fasta -m


# Now to accumulate stats
perl ~/sperl/sequence_data_scripts/getBamStats.pl -b nanopore/rumen/rumen.sorted.merged.bam
Determining raw x coverage from 1 bams...
BamName TotalReads      MappedReads     UnmappedReads   RawXCov MapXCov AvgRawChrcov    AvgMapChrcov
nanopore/rumen/rumen.sorted.merged.bam  1,068,428,447      336,440,793       731,987,654       970.349316651647        305.556347266827        317.55390201709   284.210371661399

perl ~/sperl/sequence_data_scripts/getBamStats.pl -b pacbio/rumen/rumen.sorted.merged.bam
Determining raw x coverage from 1 bams...
BamName TotalReads      MappedReads     UnmappedReads   RawXCov MapXCov AvgRawChrcov    AvgMapChrcov
pacbio/rumen/rumen.sorted.merged.bam    1,028,994,158      437,357,636       591,636,522       328.668224704145        139.695212715602        152.607518056304  137.119932416167

# Generating GC percentage data on each contig
python3 calcGCcontentFasta.py -f rumen_pacbio_pilot_asm.fasta -o rumen_pacbio_pilot_asm.gc.tab -t 10

# Gathering free data from each contig's name
perl -ne 'chomp; @s = split(/\s+/); $j = $s[0]; $k = $s[1]; $k =~ s/len=//; $l = $s[2]; $l =~ s/reads=//; $m = $s[3]; $m =~ s/covStat=//; $o = $s[8]; print "$j\t$k\t$l\t$m\t$o\n";' < rumen_pacbio_pilot_asm.gc.tab > rumen_pacbio_pilot_asm.combined.stats.tab
# columns are: 
# tigname	length	numreads	covStat	gcperc

samtools idxstats pacbio/rumen/rumen.sorted.merged.bam | grep -v '*' > rumen_pacbio_pilot_asm.idxstats.tab
```

I'm going to take this to R to start drawing plots and identifying correlations.

```R
data <- read.delim("rumen_pacbio_pilot_asm.combined.stats.tab", header=FALSE)
idxdata <- read.delim("rumen_pacbio_pilot_asm.idxstats.tab", header=FALSE)

colnames(data) <- c("contig", "length", "pbreads", "pbcov", "gcperc")
colnames(idxdata) <- c("contig", "length", "ilmap", "ilunmap")

library(dplyr)
fulldata <- left_join(data, idxdata, by="contig")
fulldata.filter <- select(fulldata, -length.y)

# First, let's see how the data behaves in a correlation matrix
cor(fulldata.filter[,2:7])
            length    pbreads      pbcov     gcperc     ilmap    ilunmap
length  1.00000000 0.50511970 0.97630081 0.10015184 0.1203805 0.01315528
pbreads 0.50511970 1.00000000 0.51619721 0.04209997 0.1579580 0.01937891
pbcov   0.97630081 0.51619721 1.00000000 0.08682802 0.1226578 0.01346729
gcperc  0.10015184 0.04209997 0.08682802 1.00000000 0.0265572 0.01609518
ilmap   0.12038054 0.15795800 0.12265784 0.02655720 1.0000000 0.98150715
ilunmap 0.01315528 0.01937891 0.01346729 0.01609518 0.9815071 1.00000000

# So, pacbio coverage is a stronger determinant of contig length and is only mildly biased by GC percentage
# There's a very low correlation between the contigs size, number of reads and pbcoverage with illumina reads though
cov(fulldata.filter[,2:7])
              length      pbreads        pbcov       gcperc        ilmap
length  1.160104e+08 3.478463e+05 4.541154e+07 9.581425e+01 1.863762e+08
pbreads 3.478463e+05 4.087796e+03 1.425261e+05 2.390832e-01 1.451683e+06
pbcov   4.541154e+07 1.425261e+05 1.864955e+07 3.330556e+01 7.614040e+07
gcperc  9.581425e+01 2.390832e-01 3.330556e+01 7.889427e-03 3.390710e+02
ilmap   1.863762e+08 1.451683e+06 7.614040e+07 3.390710e+02 2.066195e+10
ilunmap 3.719469e+06 3.252418e+04 1.526675e+06 3.752761e+01 3.703493e+09
             ilunmap
length  3.719469e+06
pbreads 3.252418e+04
pbcov   1.526675e+06
gcperc  3.752761e+01
ilmap   3.703493e+09
ilunmap 6.890721e+08

# gcpercentage has a minor covariance with the illumina datasets, but is nothing to write home about
# Printing out the correlations to show later
library(corrplot)
pdf(file="pacbio_asm_correlation_plot.pdf", useDingbats=FALSE)
corrplot(cor(fulldata.filter[,2:7]), method="shade")
dev.off()

# Printing out a matrix of scatterplots
pdf(file="pacbio_asm_scatterplots.pdf", useDingbats=FALSE)
pairs(~length+pbreads+pbcov+gcperc+ilmap+ilunmap, data=fulldata.filter, main="PacBio summary stats plot")
dev.off()

# I should really check the correlation of illumina read coverage instead of counts
# Also adding a proportion of unmapped illumina reads (based on the number of mapped reads)
fulldata.filter <- mutate(fulldata.filter, ilcov = (ilmap * 150) / length, ilumapprop = ilunmap / ilmap)

# Some final summary statistics
summary(fulldata.filter[,2:9])
     length          pbreads            pbcov              gcperc
 Min.   :  1014   Min.   :   2.00   Min.   : -5113.3   Min.   :0.003254
 1st Qu.:  6304   1st Qu.:   3.00   1st Qu.:   494.4   1st Qu.:0.433027
 Median :  9145   Median :   6.00   Median :  1487.8   Median :0.491760
 Mean   : 11805   Mean   :  14.13   Mean   :  2767.4   Mean   :0.475510
 3rd Qu.: 14041   3rd Qu.:  12.00   3rd Qu.:  3381.8   3rd Qu.:0.534069
 Max.   :247015   Max.   :7905.00   Max.   :101547.3   Max.   :0.991419
     ilmap             ilunmap            ilcov            ilumapprop
 Min.   :      31   Min.   :      0   Min.   :     0.9   Min.   :0.00000
 1st Qu.:    2039   1st Qu.:    180   1st Qu.:    41.3   1st Qu.:0.05700
 Median :    4571   Median :    410   Median :    70.5   Median :0.09077
 Mean   :   10994   Mean   :    884   Mean   :   137.1   Mean   :0.10406
 3rd Qu.:   10055   3rd Qu.:    843   3rd Qu.:   122.6   3rd Qu.:0.13656
 Max.   :28080943   Max.   :5221720   Max.   :365098.5   Max.   :0.75493

# Holy crap! There are some contigs with high mapped reads!
fulldata.filter[fulldata.filter$ilcov > 10000,]
           contig length pbreads   pbcov    gcperc    ilmap ilunmap     ilcov
13108 tig00034561  11537     157 2552.65 0.9914189 28080943 5221720 365098.50
25384 tig00058960   2477       3  304.61 0.2006459   201633   76230  12210.31
38144 tig00092998   2571       6  308.03 0.1979774   261847   30359  15276.95
      ilumapprop
13108  0.1859524
25384  0.3780631
38144  0.1159418

# OK, I checked in samtools and tig00034561 is almost all C's
# tig00058960 has closest hit in a Musa acuminata kinase
# tig00092998 has closest hit in a Macaca mulatta BAC

fulldata.filter[fulldata.filter$ilumapprop > 0.5,]

# tig00050777 had the highest proportion of unmapped reads at 0.75 
# the last 200bp of the contig hit Gossypium hirsutum, which suggests plant contamination?

# Let's try to cluster the dataset
distances <- dist(fulldata.filter[,2:9])
summary(distances)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
       9     4942     9476    19500    18780 28560000

nrow(fulldata.filter)
[1] 39781

# We don't want to cluster every contig into a unique cluster. Let's stop with a k of 40 to start
twcss <- sapply(1:kmax, function(k){kmeans(distances, k)$tot.withinss})
# It took too long.
```

<a name="binning"></a>
## Metagenomic binning test

In order to help Serge, I promised to bin the error corrected Pacbio and Nanopore reads to try to separate them into separate clusters. First, to gather all of the data I need and to prepare the files.

#### MetaBat binning

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/error_corrected_reads

```bash
# creating indicies for the reads -- may take a long time!
for i in rumen_nanopore_corrected.fasta rumen_pacbio_corrected.fasta; do echo $i; sbatch --mem=10000 --ntasks-per-node=1 --nodes=1 --wrap="module load bwa; bwa index $i"; sbatch --mem=4000 --ntasks-per-node=1 --nodes=1 --wrap="module load samtools; samtools faidx $i"; done

# Damn, I think that the script ran out of memory
for i in rumen_nanopore_corrected.fasta rumen_pacbio_corrected.fasta; do echo $i; sbatch --mem=30000 --ntasks-per-node=1 --nodes=1 --wrap="module load bwa; bwa index $i"; done

# The error corrected pacbio reads indexing is taking a long time for BWA!
# Going to start the nanopore correction first
perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b nanopore_reads -t illumina_reads.tab -f rumen_nanopore_corrected.fasta -m

# It's taking far too long, going to try minimap2 instead
# Indexing reference
/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -d rumen_pacbio_corrected.mmi rumen_pacbio_corrected.fasta

# And running the illumina alignments
perl -lane '@fsegs = split(/\//, $F[0]); @nsegs = split(/\./, $fsegs[-1]); $cmd = "/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -x sr rumen_pacbio_corrected.mmi $F[0] $F[1] > $nsegs[0].out.paf"; system(qq{sbatch --nodes=1 --ntasks-per-node=1 --mem=25000 --wrap="$cmd"});' < illumina_reads.tab

# Finally, I'm going to count the alignments and generate the matrix of read depth
sbatch -p assemble3 pafChrCounter.pl pacbio_corrected_run3_reads.count YMPrepCannula_run3_L2_R1.out.paf YMPrepCannula_run3_L3_R1.out.paf YMPrepCannula_run3_L4_R1.out.paf
sbatch -p assemble3 pafChrCounter.pl pacbio_corrected_run1_reads.count YMPrepCannula_S1_L001_R1_001.out.paf YMPrepCannula_S1_L003_R1_001.out.paf YMPrepCannula_S1_L004_R1_001.out.paf

# The PAF files do not have best alignments, so they are not completely worthwhile. 
# Serge recommended that I generate BAMs instead
module load samtools; perl -lane '@fsegs = split(/\//, $F[0]); @nsegs = split(/\./, $fsegs[-1]); $cmd = "/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -ax sr rumen_pacbio_corrected.mmi $F[0] $F[1] | samtools sort -m 2G -o $nsegs[0].out.bam -T $nsegs[0].temp -"; system(qq{sbatch --nodes=1 --ntasks-per-node=2 --mem=27000 --wrap="$cmd"});' < illumina_reads.tab

# Damn, it doesn't output sam headers, so samtools sort crashed!
module load samtools; perl -lane '@fsegs = split(/\//, $F[0]); @nsegs = split(/\./, $fsegs[-1]); $cmd = "/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -ax sr rumen_pacbio_corrected.mmi $F[0] $F[1] > $nsegs[0].out.sam"; system(qq{sbatch --nodes=1 --ntasks-per-node=2 --mem=27000 --wrap="$cmd"});' < illumina_reads.tab

# The larger run3 file has taken over 5 days and has generated 2.5 billion mappings. I'm going to try to condense one of the files into a sorted bam to test out the compression
# I don't think this will hurt the final results, as I may have mapped a large proportion of the reads already
sbatch -p assemble1 --mem=28000 --ntasks-per-node=2 --nodes=1 --wrap="module load samtools; samtools view -bht rumen_pacbio_corrected.fasta.fai YMPrepCannula_run3_L2_R1.out.sam | samtools sort -T YMPrepCannula_run3_L2_R1.temp -o YMPrepCannula_run3_L2_R1.sorted.bam -m 15G -"
sbatch -p assemble1 --mem=28000 --ntasks-per-node=2 --nodes=1 --wrap="module load samtools; samtools view -bht rumen_pacbio_corrected.fasta.fai YMPrepCannula_run3_L3_R1.out.sam | samtools sort -T YMPrepCannula_run3_L3_R1.temp -o YMPrepCannula_run3_L3_R1.sorted.bam -m 15G -"
sbatch -p assemble1 --mem=28000 --ntasks-per-node=2 --nodes=1 --wrap="module load samtools; samtools view -bht rumen_pacbio_corrected.fasta.fai YMPrepCannula_run3_L4_R1.out.sam | samtools sort -T YMPrepCannula_run3_L4_R1.temp -o YMPrepCannula_run3_L4_R1.sorted.bam -m 15G -"

for i in YMPrepCannula_run3_L2_R1.sorted.bam YMPrepCannula_run3_L3_R1.sorted.bam YMPrepCannula_run3_L4_R1.sorted.bam YMPrepCannula_S1_L003_R1_001.sorted.bam YMPrepCannula_S1_L004_R1_001.sorted.bam; do echo $i; sbatch --mem=6000 --ntasks-per-node=1 --nodes=1 --wrap="module load samtools; samtools index $i;"; done

/mnt/nfs/nfs2/bickhart-users/binaries/metabat/jgi_summarize_bam_contig_depths --outputDepth YMPrepCannula_jgi_depth_file.txt YMPrepCannula_run3_L2_R1.sorted.bam YMPrepCannula_run3_L3_R1.sorted.bam YMPrepCannula_run3_L4_R1.sorted.bam YMPrepCannula_S1_L001_R1_001.sorted.bam YMPrepCannula_S1_L003_R1_001.sorted.bam YMPrepCannula_S1_L004_R1_001.sorted.bam

sbatch --nodes=1 -p assemble3 --ntasks-per-node=3 --mem=100000 --wrap="metabat2 -i rumen_pacbio_corrected.fasta.gz -a YMPrepCannula_jgi_depth_file.txt -o YMP/bin -v"
MetaBAT 2 (v2.12.1) using minContig 2500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, and maxEdges 200.
[00:16:56] Finished reading 5931682 contigs and 6 coverages from YMPrepCannula_jgi_depth_file.txt
[00:17:10] Number of target contigs: 427181 of large (>= 2500) and 9357 of small ones (>=1000 & <2500).
[00:22:51] Finished TNF calculation.
[00:45:55] Finished Preparing TNF Graph Building [pTNF = 97.06]
[01:22:23] Finished Building TNF Graph (405472 vertices and 17317295 edges) [46.6Gb / 503.9Gb]                            
[01:23:08] Building SCR Graph and Binning (40583 vertices and 37093 edges) [P = 9.50%; 46.5Gb / 503.9Gb]                  [01:23:13] Building SCR Graph and Binning (81166 vertices and 106764 edges) [P = 19.00%; 46.5Gb / 503.9Gb]                [01:23:17] Building SCR Graph and Binning (121747 vertices and 209892 edges) [P = 28.50%; 46.5Gb / 503.9Gb]               [01:23:24] Building SCR Graph and Binning (162329 vertices and 341147 edges) [P = 38.00%; 46.5Gb / 503.9Gb]               [01:23:47] Building SCR Graph and Binning (202911 vertices and 545071 edges) [P = 47.50%; 46.5Gb / 503.9Gb]               [01:24:39] Building SCR Graph and Binning (243494 vertices and 875160 edges) [P = 57.00%; 46.5Gb / 503.9Gb]               [01:25:20] Building SCR Graph and Binning (284076 vertices and 1230131 edges) [P = 66.50%; 46.6Gb / 503.9Gb]              [01:26:01] Building SCR Graph and Binning (324658 vertices and 1595284 edges) [P = 76.00%; 46.6Gb / 503.9Gb]              [01:26:38] Building SCR Graph and Binning (365240 vertices and 2156485 edges) [P = 85.50%; 46.6Gb / 503.9Gb]              [01:27:30] Building SCR Graph and Binning (376812 vertices and 2478338 edges) [P = 95.00%; 46.6Gb / 503.9Gb]              
[01:27:37] 0.39% (8077936 bases) of large (>=2500) contigs were re-binned out of small bins (<200000).
[01:28:29] 47.60% (2071784572 bases) of large (>=2500) and 0.81% (132587 bases) of small (<2500) contigs were binned.
1563 bins (2071917159 bases in total) formed.

export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin:$PATH
module load hmmer/3.1b1
checkm data setRoot /mnt/nfs/nfs2/bickhart-users/binaries/CheckM/databases

checkm lineage_wf -f YMP/CheckM.txt -t 8 -x fa YMP YMP/SCG
 Determining marker sets for each genome bin.
    Finished processing 1563 of 1563 (100.00%) bins (current: bin.796).

  Marker set written to: YMP/SCG/lineage.ms

  { Current stage: 0:04:01.640 || Total: 1:05:23.517 }
  Calculating AAI between multi-copy marker genes.

  Reading HMM info from file.
  Parsing HMM hits to marker genes:
    Finished parsing hits for 1563 of 1563 (100.00%) bins.

  QA information written to: YMP/CheckM.txt

  { Current stage: 0:05:38.476 || Total: 2:36:20.053 }

```

Now to check the benchmarking in R.

```R
# First run
source('http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/benchmark.R')
printPerf(list(calcPerfBySCG("./CheckM.txt", removeStrain=F)), rec=c(seq(.1,.9,.1),.95), prec=c(seq(.6,.9,.1),.95,.99)) [[1]]
         Recall
Precision 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
     0.6   10  10   5   0   0   0   0   0   0    0
     0.7    6   6   2   0   0   0   0   0   0    0
     0.8    5   5   1   0   0   0   0   0   0    0
     0.9    0   0   0   0   0   0   0   0   0    0
     0.95   0   0   0   0   0   0   0   0   0    0
     0.99   0   0   0   0   0   0   0   0   0    0

printPerf(list(calcPerfBySCG("./CheckM.txt", removeStrain=T)), rec=c(seq(.1,.9,.1),.95), prec=c(seq(.6,.9,.1),.95,.99))
[[1]]
         Recall
Precision 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
     0.6   18  18  13   5   1   0   0   0   0    0
     0.7   17  17  12   4   1   0   0   0   0    0
     0.8   12  12   7   1   0   0   0   0   0    0
     0.9    6   6   2   0   0   0   0   0   0    0
     0.95   3   3   1   0   0   0   0   0   0    0
     0.99   0   0   0   0   0   0   0   0   0    0


# Not very good! Our higher precision reads are all from contaminants

```

#### MetaProb testing

The attractive features of MetaProb include the ability of the program to bin without requiring read depth. Unfortunately, mapping short reads to longer reads takes allot of time!

> pwd: /home/dbickhart/share/metagenomics/read_binning

```bash
../../test_software/MetaProb/Release/MetaProb -si rumen_pacbio_corrected.fasta -numSp 5000 -feature 2 -m 45 -dirOutput pacbio_metaprob_binned

Directory output: pacbio_metaprob_binned
Files: rumen_pacbio_corrected.fasta
N. Cluster: 5000
Parameter q: 111111111111111111111111111111
Parameter m: 45
Parameter SeedSize: 9000
Parameter lmerfreq: 4
Parameter Kmeans iteration max: 100
Norm: D2star_Group_Prob_Bernulli_Lmer_Euclidian
Graph Type: Paired

```


## Alignment percentage checks

In this project, I am going to check the alignment percentage of reads to clusters in our assembly and Mick's assemblies. It will be a "cross-alignment" validation of read composition. I am going to queue up jobs on each cluster individually, and align the 

> /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project

```bash
ls /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina/YMPrepCannula_run3_L* > raw_illumina_fastas.tab
ls /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/micks_reads/* >> raw_illumina_fastas.tab

# Testing the run to see how it progresses first
sbatch process_individual_clusters_forRD.pl raw_illumina_fastas.tab best_genome_clusters/cluster.100.fasta testrun
# took 7 hours

# Now testing with minimap2 to see if we can speed up alignment
sbatch process_individual_clusters_forRD.pl raw_illumina_fastas.tab best_genome_clusters/cluster.100.fasta testrunM

# Looks good. Running on all of the USDA samples
ls best_genome_clusters/*.fasta | xargs -I {} sbatch process_individual_clusters_forRD.pl raw_illumina_fastas.tab {} {}

# And mick's clusters
ls micks_clusters/*.fasta | xargs -I {} sbatch -p assemble1 process_individual_clusters_forRD.pl raw_illumina_fastas.tab {} {}

# It looks like the analysis completed successfully. I need to grab some more stats from the files to improve my comparisons
# Let's get the overal GC content and fasta lengths first
module load samtools

for i in best_genome_clusters/*.final.bam; do name=`echo $i | cut -d'.' -f1,2`; echo $name; samtools faidx $name.fasta; python3 error_corrected_reads/calcGCcontentFasta.py -f $name.fasta -o $name.gc -t 10; done
for i in micks_clusters/*.final.bam; do name=`echo $i | cut -d'.' -f1,2`; echo $name; samtools faidx $name.fasta; python3 error_corrected_reads/calcGCcontentFasta.py -f $name.fasta -o $name.gc -t 10; done

perl -e 'my $fold = "best_genome_clusters"; <>; %f; while(<>){chomp; @s = split(/\t/); $f{$s[0]} = 1;} print "Cluster\tContig\tGC\tLen\n"; foreach my $t (sort{$a cmp $b} keys(%f)){my %l; my %g; open(IN, "< $fold/$t.fai"); while(<IN>){chomp; @s = split(/\t/); $l{$s[0]} = $s[1];} close IN; $h = $t; $h =~ s/\.fasta//; open(IN, "< $fold/$h.gc"); while(<IN>){chomp; @s = split(/\t/); $g{$s[0]} = $s[1];} close IN; foreach my $k (sort{$a cmp $b} keys(%l)){print "$t\t$k\t$g{$k}\t$l{$k}\n";}}' < usda_clusters_rg_counts.tab > usda_clusters_gclen.tab
perl -e 'my $fold = "micks_clusters"; <>; %f; while(<>){chomp; @s = split(/\t/); $f{$s[0]} = 1;} print "Cluster\tContig\tGC\tLen\n"; foreach my $t (sort{$a cmp $b} keys(%f)){my %l; my %g; open(IN, "< $fold/$t.fai"); while(<IN>){chomp; @s = split(/\t/); $l{$s[0]} = $s[1];} close IN; $h = $t; $h =~ s/\.fasta//; open(IN, "< $fold/$h.gc"); while(<IN>){chomp; @s = split(/\t/); $g{$s[0]} = $s[1];} close IN; foreach my $k (sort{$a cmp $b} keys(%l)){print "$t\t$k\t$g{$k}\t$l{$k}\n";}}' < micks_clusters_rg_counts.tab > micks_clusters_gclen.tab

# I also have a hunch that cross-mapping validation may be a good sign of repetitive or contaminant quality
# I updated the script to count map quality and zeromapQ reads
sbatch -p assemble3 calculate_cluster_read_depth.pl -f best_genome_clusters -o usda_clusters_rg_counts.ext.tab
sbatch -p assemble3 calculate_cluster_read_depth.pl -f micks_clusters -o micks_clusters_rg_counts.ext.tab
```

I am plotting these on my local VM to allow me to interrogate the plots in real time.

```R
library(dplyr)
library(ggplot2)

mickdata <- read.delim("micks_clusters_rg_counts.tab")
# Mutating data for text plotting
mickdata <- mutate(mickdata, outlier=ifelse(Count > 1000000, Cluster, as.numeric(NA)))
# Damn R classes!
mickdata <- mutate(mickdata, Label=ifelse(outlier, as.character(Cluster), ""))

ggplot(mickdata, aes(x=ReadGroup, y=Count)) + geom_boxplot() + scale_y_log10() + labs(y="Read Count (log10)") + ggtitle("Alignments to Mick's Clusters") + geom_text(aes(label=Label), na.rm=TRUE, hjust= -0.3)
dev.copy2pdf(file="mick_contig_aligns.boxplot.pdf", useDingbats=FALSE)


# Let's try the same with my dataset
usdadata <- read.delim("usda_clusters_rg_counts.tab")
usdadata <- mutate(usdadata, Label=ifelse((ReadGroup == "MICK" & Count > 1000000) | (ReadGroup == "USDA" & Count > 5000000), Cluster, as.numeric(NA)))
usdadata <- mutate(usdadata, Label=ifelse(Label, as.character(Cluster), ""))

ggplot(usdadata, aes(x=ReadGroup, y=Count)) + geom_boxplot() + scale_y_log10() + labs(y="Read Count (log10)") + ggtitle("Alignments to USDA Clusters") + geom_text(aes(label=Label), na.rm=TRUE, hjust= -0.3)
dev.copy2pdf(file="usda_contig_aligns.boxplot.pdf", useDingbats=FALSE)

# Now I'm going to estimate how the read counts vary by GC percentage and contig length
usdagc <- read.delim("usda_clusters_gclen.tab")
library(tidyr)

# Condensing the data and associating with GC and len
usdadata <- usdadata %>% spread(ReadGroup, Count) %>% mutate(Diff = USDA - MICK)
usdadata.join <- inner_join(usdadata, usdagc, by=c("Cluster", "Contig"))

# Damn joining program didn't associate all values!
usdadata.filter <- usdadata.join %>% filter(!is.na(Diff))
cor(usdadata.filter[,c(6,7,8)])
# So, a correlation between contig length and the proportion of USDA reads mapping, but Mick's reads are not mapping
           Diff          GC         Len
Diff  1.0000000 -0.01417650  0.33446263
GC   -0.0141765  1.00000000 -0.04456602
Len   0.3344626 -0.04456602  1.00000000

library("PerformanceAnalytics")
chart.Correlation(usdadata.filter[,c(6,7,8)], histogram=TRUE, pch=19)
dev.copy2pdf(file="usda_correlation_plots_alignsvscontigs.pdf", useDingbats=FALSE)

# The longer contigs with high difference values may be contaminants. Let's plot them separately
# This plot showed similar results to the correlation plots, so I used it to orient my tag printing
ggplot(usdadata.filter, aes(x = Diff, y = Len, col = GC)) + geom_jitter() + geom_smooth(lwd = 2, se = FALSE, method="loess") + theme_bw()

usdadata.filter <- mutate(usdadata.filter, Label = ifelse(Diff > 1500000 | Diff < -500000, Cluster, as.numeric(NA)))
usdadata.filter <- mutate(usdadata.filter, Label=ifelse(Label, as.character(Cluster), ""))

ggplot(usdadata.filter, aes(y = Diff, x = Len, col = GC)) + geom_jitter() + geom_smooth(lwd = 2, se = FALSE, method="loess", col = "red") + geom_text(aes(label=Label), na.rm=TRUE, hjust=-0.3) + theme_bw() + labs(x="Contig Length (bp)", y = "Read Count Difference (USDA - MICK)") + ggtitle("USDA Contig length vs ReadCount Difference")
dev.copy2pdf(file="usda_len_vs_diff_scatter.pdf", useDingbats=FALSE)

# The right-most value was very interesting, so I wanted to take a look at that cluster to see if it was a contaminant
usdadata.filter[usdadata.filter$Cluster == "cluster.621.fasta",]
                Cluster      Contig             Label    MICK    USDA    Diff
27860 cluster.621.fasta tig00122470              <NA>     780   19254   18474
27861 cluster.621.fasta tig00125079              <NA>     449   15582   15133
27862 cluster.621.fasta tig00178536              <NA>    1610   13046   11436
27863 cluster.621.fasta tig00497612 cluster.621.fasta 2241956 6051057 3809101
27864 cluster.621.fasta tig00497842              <NA>   22123   67783   45660
             GC    Len
27860 0.5412320   4026
27861 0.5303644   3211
27862 0.5101818   2750
27863 0.5707110 438390
27864 0.5222579  44838

# Looks like it may be legit! Or a contaminant that's shared between datasets.

# Now to look at Mick's data
mickgc <- read.delim("micks_clusters_gclen.tab")
mickdata.temp <- mickdata %>% spread(ReadGroup, Count) %>% mutate(Diff = MICK - USDA)
mickdata.join <- inner_join(mickdata.temp, mickgc, by=c("Cluster", "Contig"))

mickdata.filter <- mickdata.join %>% filter(!is.na(Diff))
cor(mickdata.filter[,c(7,8,9)])
          Diff         GC        Len
Diff 1.0000000 0.10804427 0.27144416
GC   0.1080443 1.00000000 0.05880973
Len  0.2714442 0.05880973 1.00000000

# Similar to the USDA contigs, but Mick has a slight change with the GC bias
chart.Correlation(mickdata.filter[,c(7,8,9)], histogram=TRUE, pch=19)
dev.copy2pdf(file="mick_correlation_plots_alignsvscontigs.pdf", useDingbats=FALSE)

# The difference values in Mick's data are more evenly distributed. I'll plot a comparison histogram of them against our dataset last
mickdata.filter <- mutate(mickdata.filter, Label = ifelse(Diff > 300000 | Diff < -300000, Cluster, as.numeric(NA)))
mickdata.filter <- mutate(mickdata.filter, Label = ifelse(Label, as.character(Cluster), ""))

# Given the wider range of "Difference values" in our datasets, it's becoming far clearer that the USDA sample has far more diversity
# Mick's sample is relatively even, by comparison!
ggplot(mickdata.filter, aes(y= Diff, x = Len, col = GC)) + geom_jitter() + geom_smooth(lwd = 2, se = FALSE, method="loess", col = "red") + geom_text(aes(label=Label), na.rm=TRUE, hjust=-0.3) + theme_bw() + labs(x="Contig Length (bp)", y = "Read Count Difference (MICK - USDA)") + ggtitle("MICK Contig length vs ReadCount Difference")
dev.copy2pdf(file="mick_len_vs_diff_scatter.pdf", useDingbats=FALSE)

# Now to try to plot a histogram of the difference values from the two respective datasets
join.diff <- mickdata.filter %>% select(Diff) %>% mutate(Data = c("MICK"))
join.diff <- bind_rows(join.diff, usdadata.filter %>% select(Diff) %>% mutate(Data = c("USDA")))

ggplot(join.diff, aes(x = Diff, fill = Data)) + geom_density(alpha = 0.30) + theme_bw() + xlim(c(-50000, 50000)) + labs(x = "Differences in Read Depth per Contig", y = "Density") + ggtitle("Read Count Differences in Mick's Clusters vs USDA's Clusters")

```

Now to combine all of the data so that I can regenerate the plots with additional mapping quality data.

```bash
perl -e 'my $fold = "best_genome_clusters"; <>; %f; while(<>){chomp; @s = split(/\t/); $f{$s[0]}->{$s[1]}->{$s[2]} = [$s[3], $s[4], $s[5]];} print "Cluster\tContig\tReadGroup\tCount\tAvgQual\tNonZeroQ\tGC\tLen\n"; foreach my $t (sort{$a cmp $b} keys(%f)){my %l; my %g; open(IN, "< $fold/$t.fai"); while(<IN>){chomp; @s = split(/\t/); $l{$s[0]} = $s[1];} close IN; $h = $t; $h =~ s/\.fasta//; open(IN, "< $fold/$h.gc"); while(<IN>){chomp; @s = split(/\t/); $g{$s[0]} = $s[1];} close IN; foreach my $k (sort{$a cmp $b} keys(%l)){foreach my $rg (keys(%{$f{$t}->{$k}})){my $orig = join("\t", @{$f{$t}->{$k}->{$rg}}); print "$t\t$k\t$rg\t$orig\t$g{$k}\t$l{$k}\n";}}}' < usda_clusters_rg_counts.ext.tab > usda_clusters_rg_counts.ext.full.tab

perl -e 'my $fold = "micks_clusters"; <>; %f; while(<>){chomp; @s = split(/\t/); $f{$s[0]}->{$s[1]}->{$s[2]} = [$s[3], $s[4], $s[5]];} print "Cluster\tContig\tReadGroup\tCount\tAvgQual\tNonZeroQ\tGC\tLen\n"; foreach my $t (sort{$a cmp $b} keys(%f)){my %l; my %g; open(IN, "< $fold/$t.fai"); while(<IN>){chomp; @s = split(/\t/); $l{$s[0]} = $s[1];} close IN; $h = $t; $h =~ s/\.fasta//; open(IN, "< $fold/$h.gc"); while(<IN>){chomp; @s = split(/\t/); $g{$s[0]} = $s[1];} close IN; foreach my $k (sort{$a cmp $b} keys(%l)){foreach my $rg (keys(%{$f{$t}->{$k}})){my $orig = join("\t", @{$f{$t}->{$k}->{$rg}}); print "$t\t$k\t$rg\t$orig\t$g{$k}\t$l{$k}\n";}}}' < micks_clusters_rg_counts.ext.tab > micks_clusters_rg_counts.ext.full.tab
```


Finally, to complete my analysis with the addition of mapping data.

```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(PerformanceAnalytics)

mickfull <- read.delim("micks_clusters_rg_counts.ext.full.tab")
# Let's get a more meaningful proportion of non-zero quality reads
mickfull <- mutate(mickfull, ZeroProp = ifelse(Count > 0, 1 - (NonZeroQ / Count), 0))

# Drawing out only the metrics I will plot
mickfull.metrics <- select(mickfull, -NonZeroQ)
mickfull.metrics.mick <- filter(mickfull, ReadGroup == "MICK") %>% select(-ReadGroup)
mickfull.metrics.usda <- filter(mickfull, ReadGroup == "USDA") %>% select(-ReadGroup)

cor(mickfull.metrics.mick[,c(3,4,6,7,8)])
               Count     AvgQual          GC         Len     ZeroProp
Count     1.00000000 -0.01697127  0.15250212 0.430468302  0.021225964
AvgQual  -0.01697127  1.00000000 -0.18075584 0.057938474 -0.227404773
GC        0.15250212 -0.18075584  1.00000000 0.060272939 -0.010140248
Len       0.43046830  0.05793847  0.06027294 1.000000000  0.003730061
ZeroProp  0.02122596 -0.22740477 -0.01014025 0.003730061  1.000000000

cor(mickfull.metrics.usda[,c(3,4,6,7,8)])
               Count     AvgQual           GC         Len    ZeroProp
Count    1.000000000  0.02204298  0.009192094  0.02702579  0.03856566
AvgQual  0.022042977  1.00000000 -0.147042915 -0.14880451 -0.13634277
GC       0.009192094 -0.14704292  1.000000000  0.05881794 -0.01373861
Len      0.027025790 -0.14880451  0.058817943  1.00000000  0.02369836
ZeroProp 0.038565663 -0.13634277 -0.013738608  0.02369836  1.00000000

# So our read alignment qualities decrease on Mick's larger contigs
chart.Correlation(mickfull.metrics.mick[,c(3,4,6,7,8)], histogram=TRUE, pch=19, title="MICK alignment to MICK clusters")
dev.copy2pdf(file="mick_extend_mick_correlation.pdf", useDingbats=FALSE)

chart.Correlation(mickfull.metrics.usda[,c(3,4,6,7,8)], histogram=TRUE, pch=19)
dev.copy2pdf(file="mick_extend_usda_correlation.pdf", useDingbats=FALSE)

png("mick_extend_mick_correlation.png", width=1200, height=800)
chart.Correlation(mickfull.metrics.mick[,c(3,4,6,7,8)], histogram=TRUE, pch=19)
There were 50 or more warnings (use warnings() to see the first 50)
dev.off()

png("mick_extend_usda_correlation.png", width=1200, height=800)
chart.Correlation(mickfull.metrics.usda[,c(3,4,6,7,8)], histogram=TRUE, pch=19)
There were 50 or more warnings (use warnings() to see the first 50)
dev.off()

# Now to do the same for the USDA clusters
usdafull <- read.delim("usda_clusters_rg_counts.ext.full.tab")
usdafull <- mutate(usdafull, ZeroProp = ifelse(Count > 0, 1 - (NonZeroQ / Count), 0))

summary(usdafull) # Interesting... the ZeroProp read average is far higher than in Mick's data
usdafull.metrics <- select(usdafull, -NonZeroQ)

usdafull.metrics.mick <- filter(usdafull.metrics, ReadGroup == "MICK") %>% select(-ReadGroup)
usdafull.metrics.usda <- filter(usdafull.metrics, ReadGroup == "USDA") %>% select(-ReadGroup)

cor(usdafull.metrics.mick[,c(3,4,5,6,7)])
              Count     AvgQual           GC         Len     ZeroProp
Count    1.00000000  0.20159698  0.053087937  0.24155250  0.046192754
AvgQual  0.20159698  1.00000000 -0.097552829  0.01416813 -0.183985431
GC       0.05308794 -0.09755283  1.000000000 -0.05433986 -0.008441778
Len      0.24155250  0.01416813 -0.054339858  1.00000000  0.054342640
ZeroProp 0.04619275 -0.18398543 -0.008441778  0.05434264  1.000000000

cor(usdafull.metrics.usda[,c(3,4,5,6,7)])
              Count     AvgQual           GC         Len     ZeroProp
Count     1.0000000 -0.04782600  0.010521697  0.23900193  0.124856845
AvgQual  -0.0478260  1.00000000 -0.258912729  0.05120062 -0.401512812
GC        0.0105217 -0.25891273  1.000000000 -0.03894586 -0.003804244
Len       0.2390019  0.05120062 -0.038945857  1.00000000  0.060269000
ZeroProp  0.1248568 -0.40151281 -0.003804244  0.06026900  1.000000000

# It looks like some of our contigs have issues with repetition and the higher GC stuff tends to have more issues
png("usda_extend_mick_correlation.png", width=1200, height=800)
chart.Correlation(usdafull.metrics.mick[,c(3,4,5,6,7)], histogram=TRUE, pch=19)
dev.off()

png("usda_extend_usda_correlation.png", width=1200, height=800)
chart.Correlation(usdafull.metrics.usda[,c(3,4,5,6,7)], histogram=TRUE, pch=19)
dev.off()

# Interesting, I found that several of our clusters have no reads mapping from Micks' dataset, but there are no clusters in our dataset that have no reads from our dataset
nrow(setdiff(usdafull.metrics.usda[,c(1,2)], usdafull.metrics.mick[,c(1,2)]))
[1] 6033
usdamickzeroclusters <- setdiff(usdafull.metrics.usda[,c(1,2)], usdafull.metrics.mick[,c(1,2)]) %>% select(Cluster) %>% distinct()

# It turns out that there are contigs in each cluster that aren't represented in Mick's data though
usdamickzeroclusters <- setdiff(usdafull.metrics.usda[,c(1,2)], usdafull.metrics.mick[,c(1,2)])
summary(filter(usdafull.metrics.usda, Cluster %in% usdamickzeroclusters$Cluster & Contig %in% usdamickzeroclusters$Contig)[,c(3,4,5,6,7)])
     Count          AvgQual             GC              Len        
 Min.   :    6   Min.   : 3.618   Min.   :0.1551   Min.   :  2000  
 1st Qu.:  556   1st Qu.:45.404   1st Qu.:0.3995   1st Qu.:  5509  
 Median : 1138   Median :51.656   Median :0.4593   Median :  8270  
 Mean   : 2124   Mean   :49.391   Mean   :0.4603   Mean   : 10149  
 3rd Qu.: 2320   3rd Qu.:55.652   3rd Qu.:0.5191   3rd Qu.: 12343  
 Max.   :93871   Max.   :59.689   Max.   :0.7456   Max.   :116968  
    ZeroProp      
 Min.   :0.00000  
 1st Qu.:0.00000  
 Median :0.00000  
 Mean   :0.01478  
 3rd Qu.:0.00000  
 Max.   :0.87259
a
# It looks like most of the stuff is slightly lower length than what we had, but a correlation plot shows it is the best "behaved" of our data with very little repetitive content
png("usda_extend_notinmick_correlation.png", width=1200, height=800)
chart.Correlation(filter(usdafull.metrics.usda, Cluster %in% usdamickzeroclusters$Cluster & Contig %in% usdamickzeroclusters$Contig)[,c(3,4,5,6,7)], histogram=TRUE, pch=19)
dev.off()

# Let's plot the length and quality distributions to compare
usdafull.metrics.subset <- filter(usdafull.metrics.usda, Cluster %in% usdamickzeroclusters$Cluster & Contig %in% usdamickzeroclusters$Contig)
usdafull.subcomp <- bind_rows(usdafull.subcomp, select(usdafull.metrics.usda, Count, AvgQual, GC, Len, ZeroProp) %>% mutate(Data = c("Full")))
library(gridExtra)

p.density <- ggplot(usdafull.subcomp, aes(x = Count, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Read Depth") + xlim(c(0, 50000))
p.qual <- ggplot(usdafull.subcomp, aes(x = AvgQual, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Average Map Qual Score")
p.gc <- ggplot(usdafull.subcomp, aes(x = GC, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig GC Percentage")
p.len <- ggplot(usdafull.subcomp, aes(x = log10(Len), fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Length Log10(bp)")
p.zero <- ggplot(usdafull.subcomp, aes(x = ZeroProp, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Proportion of Zero MapQ Reads per Contig") + xlim(c(0.0, 0.10))

png("usda_contigs_distributions_not_inmicks.png", height = 1200, width = 800)
grid.arrange(p.density, p.qual, p.gc, p.len, p.zero, ncol=2, nrow=3)
dev.off()

# Now let's find the contigs that did the USDA dataset did not map to in Mick's clusters
mickusdazeroclusters <- setdiff(mickfull.metrics.mick[,c(1,2)], mickfull.metrics.usda[,c(1,2)])
nrow(mickusdazeroclusters)
[1] 4513

mickfull.metrics.subset <- filter(mickfull.metrics.mick, Cluster %in% mickusdazeroclusters$Cluster & Contig %in% mickusdazeroclusters$Contig)
mickfull.subcomp <- select(mickfull.metrics.subset, Count, AvgQual, GC, Len, ZeroProp) %>% mutate(Data = c("Subset"))
mickfull.subcomp <- bind_rows(mickfull.subcomp, select(mickfull.metrics.mick, Count, AvgQual, GC, Len, ZeroProp) %>% mutate(Data = c("Full")))

p.density <- ggplot(mickfull.subcomp, aes(x = Count, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Read Depth") + xlim(c(0,50000))
p.qual <- ggplot(mickfull.subcomp, aes(x = AvgQual, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Average Map Qual Score")
p.gc <- ggplot(mickfull.subcomp, aes(x = GC, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig GC Percentage")
p.len <- ggplot(mickfull.subcomp, aes(x = log10(Len), fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Length Log10(bp)")
p.zero <- ggplot(mickfull.subcomp, aes(x = ZeroProp, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Proportion of Zero MapQ Reads per Contig") + xlim(c(0.0, 0.10))

png("mick_contigs_distributions_not_inusda.png", height = 1200, width = 800)
grid.arrange(p.density, p.qual, p.gc, p.len, p.zero, ncol=2, nrow=3)
dev.off()
```

Some notes before I start digging into specific contigs/clusters for analysis:
* Our reads have a GC problem (NextSeq specific?)
* There appears to be far more repetitive content within our clusters (higher proportion of zero-mapped reads)
* Some of our longer contigs have higher proportions of zeromapped reads than Mick's, suggesting repetitive content in our base reads

Here are some of the "winners" of the top percentiles of each dataset for the USDA clusters

```R
# USDA clusters and USDA reads
write.table(usdafull.metrics.usda[usdafull.metrics.usda$Count > 5000000,], "usdacontigs.usdareads.highcount.tab", sep="\t", row.names=FALSE, quote=FALSE)
write.table(usdafull.metrics.usda[usdafull.metrics.usda$ZeroProp > 0.9,], "usdacontigs.usdareads.lowqual.tab", sep="\t", row.names=FALSE, quote=FALSE)

# Now to grep out the contigs that don't map with Mick's reads
write.table(filter(usdafull.metrics.usda, Cluster %in% usdamickzeroclusters$Cluster & Contig %in% usdamickzeroclusters$Contig) %>% select(Contig) %>% distinct(), "usdacontigs.usdareads.nomickreads.contigs.list", row.names=FALSE, col.names=FALSE, quote=FALSE)

```

I am just going to grep out the cluster tax IDs from the Phase genomics data quickly and then send the data to the group for digestion.

> pwd: /home/dbickhart/share/metagenomics/pilot_project/cluster_mapping

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %c; while(<IN>){chomp; @s = split(/\t/); $s[0] =~ s/\.fasta//; push(@{$c{$s[1]}}, $s[0]);} close IN; %l; open(IN, "< $ARGV[1]"); while(<IN>){chomp; foreach my $j (@{$c{$_}}){$l{$j} += 1;} } close IN; open(IN, "< $ARGV[2]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if(exists($l{$s[0]})){print "$s[0]\t$l{$s[0]}\t$s[3]\t$s[5]\t$s[6]\n";}}; close IN;' usda_clusters_rg_counts.ext.full.tab usdacontigs.usdareads.nomickreads.contigs.list best_genomes_report.tsv > usdacontigs.usdareads.nomickreads.contigs.clusterphylolist.tab

# There were 323 clusters out of 643 that had contigs with no mappings, so this was not an amazing filter
```

<a name="mash_nmdsv1"></a>
## Mash NMDS for dataset comparison

I am going to start with the existing datasets that I had for my previous agnostic SRA dataset mash sketch comparison and expand them with Micks' illumina data, our illumina data and our pacbio reads.

Unfortunately, the dataset fasta files are in various states of gzip compression, so I need to progress carefully here.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects

```bash
for i in `ls *.gz | cut -d'_' -f1 | sort | uniq`; do echo $i; sbatch --nodes=1 -p assemble3 --ntasks-per-node=2 --mem=5000 --wrap="zcat $i*.fastq.gz | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 2 -o $i.21k.1Ms -k 21 -s 1000000 -r -m 2 -"; done

for i in `ls *.fastq | cut -d'_' -f1 | sort | uniq`; do echo $i; sbatch --nodes=1 -p assemble3 --ntasks-per-node=2 --mem=5000 --wrap="cat $i*.fastq | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 2 -o $i.21k.1Ms -k 21 -s 1000000 -r -m 2 -"; done

# Mick's reads
sbatch --nodes=1 -p assemble2 --ntasks-per-node=10 --mem=5000 --wrap="zcat ../pilot_project/micks_reads/*.fastq.gz | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 10 -o Mick.21k.1Ms -k 21 -s 1000000 -r -m 2 -"

# USDA's reads
sbatch --nodes=1 -p assemble2 --ntasks-per-node=10 --mem=5000 --wrap="zcat ../pilot_project/illumina/YMPrepCannula_run3*.fastq.gz | /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 10 -o USDA.21k.1Ms -k 21 -s 1000000 -r -m 2 -"

# Moving data to a new location to keep it organized
mkdir pilot_project/first_msh_sketch
mv datasources/*.msh pilot_project/first_msh_sketch/
```

Now to go to the other folder and to start the comparison.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/first_msh_sketch

```bash
ls *.msh > combined_summary_sketches.list

# Pasting the sketches together and generating a distance matrix
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash paste -l combined_sra_sketches combined_summary_sketches.list

for i in `cat combined_summary_sketches.list`; do echo $i; /mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash dist -p 3 -t combined_sra_sketches.msh $i > $i.combined.dist; done

perl -lane 'open(IN, "< $F[0].combined.dist"); $h = <IN>; $d = <IN>; chomp $d; @dsegs = split(/\t/, $d); $dsegs[0] = $F[0]; print join("\t", @dsegs);' < combined_summary_sketches.list > first_sketches_combineddist.matrix

perl -ne '$_ =~ s/\.21k\.1Ms\.msh//g; print $_;' < first_sketches_combineddist.matrix > first_sketches_combineddist.reformat.matrix
```

Now to make the plot locally so that I can tweak it.

> pwd: /home/dbickhart/share/metagenomics/pilot_project/sra_mash_sketches

```R
library(vegan)
library(dplyr)
library(tidyr)

# Generating the formatted data matrix
data <- read.delim("first_sketches_combineddist.reformat.matrix", header=FALSE)
entries <- data[,1]
data.format <- select(data, -V1)
row.names(data.format) <- entries
colnames(data.format) <- entries

# Now grabbing sample information and listings
samples <- read.delim("sra_file_accession_list.csv", sep=",")
samples.filtered <- filter(samples, Run %in% entries)
samples.condensed <- select(samples.filtered, Run, Tissue, Sampling.Method, Breed)
samples.final <- rbind(samples.condensed, data.frame(Run = c("USDA", "Mick"), Tissue = c("Rumen", "Rumen"), Sampling.Method = c("Fiber", "Stomach Tube"), Breed = c("Dairy", "Dairy")))

# Setting sample colors
samples.final$col <- c("#1b9e77")
samples.final[samples.final$Tissue == "Rumen", 5] <- c("#d95f02")
samples.final[samples.final$Tissue == "Mammary", 5] <- c("#7570b3")
samples.final[samples.final$Tissue == "Manure", 5] <- c("#e7298a")
samples.final[samples.final$Tissue == "Oral", 5] <- c("#66a61e")
row.names(samples.final) <- samples.final[,1]
samples.final.o <- samples.final[match(rownames(data.format), samples.final$Run),]
samples.final.o[samples.final.o$Tissue == "Nasal", 5] <- c("#000000")

samples.final.o$pch <- c(1)
samples.final.o[samples.final.o$Run == "Mick",6] <- 17
samples.final.o[samples.final.o$Run == "USDA",6] <- 18


# Now for the plotting
ord <- metaMDS(as.dist(data.format), trymax=100)
# ...
# Run 67 stress 0.1439099 
# ... Procrustes: rmse 6.207925e-05  max resid 0.0003494643

plot(scores(ord), col=samples.final.o$col, pch=samples.final.o$pch)
legend("bottomleft", legend=c("Foot", "Nasal", "Rumen", "Mammary", "Manure", "Oral", "Mick", "USDA"), col=c("#1b9e77", "#000000", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#d95f02", "#d95f02"), pch=c(1,1,1,1,1,1,17,18), pt.cex=2.5)
dev.copy2pdf(file="vegan_nmds_mick_and_usda.pdf", useDingbats=FALSE)
```

## Generating read depth matricies for downstream data plots

I need to generate a data file format that can be used in binning and/or dataset partitioning. I have my  own format (with additional GC and other information), but I know that I am still missing nucleotide composition data that might be useful down the road.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project

```bash
# This will take care of the rest of the fastas that did not finish before
for i in `seq 411 588; seq 41 57`; do name="best_genome_clusters/cluster."${i}".fasta"; echo $name; sbatch -p assemble3 process_individual_clusters_forRD.pl raw_illumina_fastas.tab $name $name; done

# Now to start processing the clusters from the remade USDA Illumina assemblies
for i in `seq 1 1181`; do name="illumina_usda_clusters/cluster."${i}".fasta"; echo $name; sbatch -p assemble1 process_individual_clusters_forRD.pl raw_illumina_fastas.tab $name $name; done
for i in `seq 1182 2362`; do name="illumina_usda_clusters/cluster."${i}".fasta"; echo $name; sbatch -p assemble3 process_individual_clusters_forRD.pl raw_illumina_fastas.tab $name $name; done

# Finishing up the pacbio USDA clusters with illumina data
module load samtools; for i in `seq 411 588; seq 41 57`; do name="best_genome_clusters/cluster."${i}; echo $name; samtools faidx $name.fasta; python3 /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/sequenceData/calcGCcontentFasta.py -f $name.fasta -o $name.gc -t 25; done


sbatch -p assemble2 calculate_cluster_read_depth.pl -f best_genome_clusters -o new_usda_pacbio_clusters.rd.ext.tab
sbatch -p assemble2 calculate_cluster_read_depth.pl -f micks_clusters -o new_mick_ilm_clusters.rd.ext.tab

# Condensing all of the information together
perl condense_cluster_information.pl micks_clusters new_mick_ilm_clusters.rd.ext.tab new_mick_ilm_clusters.rd.ext.full.tab
```

##### Correcting runtime errors

Some of the alignment jobs failed to complete. I'm going to ID the ones that still need to be run and requeue them.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project

```bash
for i in illumina_usda_clusters/*.final.bam; do name=`basename $i | cut -d'.' -f1,2`; echo $name; done > usda_ilmn_clusters.done.list
for i in illumina_usda_clusters/*.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; done > usda_ilmn_clusters.total.list
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl usda_ilmn_clusters.done.list usda_ilmn_clusters.total.list
File Number 1: usda_ilmn_clusters.done.list
File Number 2: usda_ilmn_clusters.total.list
Set     Count
1       21
1;2     2281
2       81

# I'm actually surprised to see several bam files that are unique to the "done" list!
# Ah, I know why: it's because those cluster numbers actually aren't included in the folder!
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o usda_ilmn_clusters.done.list usda_ilmn_clusters.total.list

for i in `cat group_1.txt`; do echo $i; rm illumina_usda_clusters/${i}.fasta.final.bam; done
# Removing temporary bam files
for i in `cat group_2.txt`; do echo $i; rm illumina_usda_clusters/${i}.fasta.*.bam; done

for i in `cat group_2.txt`; do name="illumina_usda_clusters/"${i}".fasta"; echo $name; sbatch -p assemble2 process_individual_clusters_forRD.pl raw_illumina_fastas.tab $name $name; done
for i in illumina_usda_clusters/*.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/sequenceData/calcGCcontentFasta.py -f illumina_usda_clusters/${name}.fasta -o illumina_usda_clusters/${name}.gc -t 10; done

sbatch -p assemble2 calculate_cluster_read_depth.pl -f illumina_usda_clusters -o new_usda_ilm_clusters.rd.ext.tab
for i in illumina_usda_clusters/*.fasta; do echo $i; sbatch -p assemble1 --nodes=1 --ntasks-per-node=1 --mem=2000 --wrap="module load samtools; samtools faidx $i;"; done

perl condense_cluster_information.pl illumina_usda_clusters new_usda_ilm_clusters.rd.ext.tab new_usda_ilm_clusters.rd.ext.full.tab
```

#### Running R statistics on the Illumina data

```R
library(dplyr)
library(tidyr)
library(ggplot2)

data <- read.delim("new_usda_ilm_clusters.rd.ext.full.tab")
data <- mutate(data, ZeroProp = ifelse(Count > 0, 1 - (NonZeroQ / Count), 0))

data.metrics <- select(data, -NonZeroQ)
data.metrics.mick <-filter(data, ReadGroup == "MICK") %>% select(-ReadGroup)
data.metrics.usda <-filter(data, ReadGroup == "USDA") %>% select(-ReadGroup)

# We have ALLOT more clusters that have no reads from MICKs dataset
nrow(setdiff(data.metrics.usda[,c(1,2)], data.metrics.mick[,c(1,2)]))
[1] 129888

usdamickzeroclusters <- setdiff(data.metrics.usda[,c(1,2)], data.metrics.mick[,c(1,2)])
usdadatasubset <- filter(data.metrics.usda, Cluster %in% usdamickzeroclusters$Cluster & Contig %in% usdamickzeroclusters$Contig)

data.subcomp <- select(usdadatasubset, Count, AvgQual, GC, Len, ZeroProp) %>% mutate(Data = c("Unique"))
data.subcomp <- bind_rows(data.subcomp, select(data.metrics.usda, Count, AvgQual, GC, Len, ZeroProp) %>% mutate(Data = c("Full")))

p.density <- ggplot(data.subcomp, aes(x = Count, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Read Depth") + xlim(c(0, 10000))
p.qual <- ggplot(data.subcomp, aes(x = AvgQual, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Average Map Qual Score")
p.gc <- ggplot(data.subcomp, aes(x = GC, fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig GC Percentage")
p.len <- ggplot(data.subcomp, aes(x = log10(Len), fill=Data)) + geom_density(alpha=0.30) + theme_bw() + labs(x = "Contig Length Log10(bp)")

png("usda_contigs_not_found_in_mick.png", height=1200, width = 800)
library(gridExtra)
grid.arrange(p.density, p.qual, p.gc, p.len, ncol=2, nrow=2)
dev.off()

# The PNG image looked ugly, so I rewrote it as a pdf
pdf("usda_contigs_not_found_in_mick.pdf", useDingbats=FALSE)
grid.arrange(p.density, p.qual, p.gc, p.len, ncol=2, nrow=2)
dev.off()
```

## Running Pilon on PacBio clusters


```bash
mkdir pacbio_usda_pilon

# Testing it out on one sample
sbatch -p assemble1 pilon_process_reads.pl cluster.100 pacbio_usda_clusters/cluster.100.fasta pacbio_usda_clusters/cluster.100.fasta.final.bam pacbio_usda_pilon

# I think that it worked fairly well. Now to queue the rest up
for i in pacbio_usda_clusters/*.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -p assemble1 pilon_process_reads.pl $name $i $i.final.bam pacbio_usda_pilon ; done

mkdir pacbio_usda_corrected
cp pacbio_usda_pilon/*.fasta pacbio_usda_corrected/
tar -czvf pacbio_usda_corrected.tar.gz pacbio_usda_corrected

# Let's estimate how many indels were corrected
module load samtools; for i in pacbio_usda_pilon/*.fasta; do echo $i; samtools faidx $i; done

# So, Pilon considers indels that are PASS filter that but considers deletions with a Del filter event
perl calculate_pilon_stats.pl pacbio_usda_pilon_correction_stats.tab

# Now to condense all of the pacbio clusters and to run a long-term alignment of illumina data to them
mkdir pacbio_pilon_accumulated

cat pacbio_usda_pilon/cluster.*.fasta | perl -ne '$_ =~ s/_pilon//; print $_;' > pacbio_pilon_accumulated/total_contigs.fasta
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 --wrap="bwa index pacbio_pilon_accumulated/total_contigs.fasta"
```

## Generating blobtools plots and some summary statistics

I am going to run the blobtools pipeline to generate some stats on the pacbio dataset and the Illumina dataset before sending it all to Mick and company.

####Note that the Diamond database is located here: 
* **/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/diamond**


```bash
# Generating a combined alignment to the pilon-corrected reads
mkdir pacbio_pilon_accumulated
cat pacbio_usda_pilon/cluster.*.fasta | perl -ne '$_ =~ s/_pilon//; print $_;' > pacbio_pilon_accumulated/total_contigs.fasta
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 --wrap="bwa index pacbio_pilon_accumulated/total_contigs.fasta"

perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b aligns -t usda_illumina_fastas.tab -f total_contigs.fasta -m

perl ~/sperl/sequence_data_scripts/getBamStats.pl -b pacbio_pilon_accumulated/aligns/USDA/USDA.sorted.merged.bam
Determining raw x coverage from 1 bams...
BamName TotalReads      MappedReads     UnmappedReads   RawXCov MapXCov AvgRawChrcov    AvgMapChrcov
USDA.sorted.merged.bam     860211818       373701639       486510179       152.532528997474        66.2646744608792  68.8191602730411        61.5598801926716 <- 43.44% mapping rate for illumina reads. 

cd pacbio_pilon_accumulated
sbatch -p assemble2 --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/bin/diamond blastx --query total_contigs.fasta --db /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/diamond/uniprot_ref_proteosomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25"

# AH, the program prints everything to STDOUT for some reason. Good thing I kept the slurm log!
mv slurm-785603.out pacbio_accum_diamond_uniprot.tsv

/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools taxify -f pacbio_accum_diamond_uniprot.tsv -m /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/diamond/uniprot_ref_proteomes.taxids -s 0 -t 2 -o pacbio_accum_diamond_uniprot

# Now to generate the plots
/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools create -i total_contigs.fasta -b aligns/USDA/USDA.sorted.merged.bam -t pacbio_accum_diamond_uniprot.pacbio_accum_diamond_uniprot.tsv.taxified.out -o usda_pacbio_pilon_blobplot
[+] Parsing FASTA - total_contigs.fasta
[ERROR:5]       : Sequence header tig00126445 is not unique.

# Doh! Looks like I need to repartition the fasta to make sure there are no duplicates
# The entries are duplicates of the same contigs, so let's try this using samtools
perl -lane 'system("samtools faidx total_contigs.fasta $F[0] >> total_unique_contigs.fasta");' < total_contigs.fasta.fai

/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools create -i total_unique_contigs.fasta -b aligns/USDA/USDA.sorted.merged.bam -t pacbio_accum_diamond_uniprot.pacbio_accum_diamond_uniprot.tsv.taxified.out -o usda_pacbio_pilon_blobplot
[+] Parsing FASTA - total_unique_contigs.fasta
[+] names.dmp/nodes.dmp not specified. Retrieving nodesDB from /mnt/nfs/nfs2/bickhart-users/binaries/blobtools/data/nodesDB.txt
[%]     100%
[+] Parsing tax0 - /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated/pacbio_accum_diamond_uniprot.pacbio_accum_diamond_uniprot.tsv.taxified.out
[+] Computing taxonomy using taxrule(s) bestsum
[%]     100%
[+] Parsing bam0 - /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated/aligns/USDA/USDA.sorted.merged.bam
[+]     Checking with 'samtools flagstat'
[+]     Mapping reads = 373,701,639, total reads = 860,211,818 (mapping rate = 43.4%)
[%]     98%
[-] Based on samtools flagstat: expected 373701639 reads, 367974149 reads were parsed
[+]     Writing usda_pacbio_pilon_blobplot.USDA.sorted.merged.bam.cov
[+] Generating BlobDB and writing to file usda_pacbio_pilon_blobplot.blobDB.json
[+]     Writing usda_pacbio_pilon_blobplot.USDA.sorted.merged.bam.cov
[+] Generating BlobDB and writing to file usda_pacbio_pilon_blobplot.blobDB.json

# Generating a phylum-level (I think!) plot
/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools plot -i usda_pacbio_pilon_blobplot.blobDB.json --notitle -o usda_pacbio_phylum

/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools plot -i usda_pacbio_pilon_blobplot.blobDB.json --notitle -r genus --sort_first "no-hit,other,undef" -p 14 -o usda_pacbio_genus

# Just for kicks, but it'll be a huge mess:
/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools plot -i usda_pacbio_pilon_blobplot.blobDB.json --notitle -r species --sort_first "no-hit,other,undef" -p 14 -o usda_pacbio_species

# I saw in the tardigrade reference paper that there was a way to remove contaminants this way:
/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools plot -i usda_pacbio_pilon_blobplot.blobDB.json --notitle -r superkingdom -o usda_pacbio_supkingdom
```

So, some notes on the blobplots:
* The superkingdom plot appears to have highlighted some potential protist contigs in the upper left hand corner! (low GC, high coverage)
* We also have 1.23 megabases of putative viral origin contigs!
* Other Eukaryotes appear to make up a smaller proportion (but larger contig sizes) of the dataset.

Now to generate some tables to try to grep out other stats.

```bash
# Species level view to see contamination
/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools view -i usda_pacbio_pilon_blobplot.blobDB.json -r species -b -o usda_pacbio_species

# This gave me a good rundown of genus names from that list in the meantime
perl -lane 'if($F[0] =~ /^#/){next;}else{print $F[5];}' < usda_pacbio_species.usda_pacbio_pilon_blobplot.blobDB.table.txt | perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -c 0 -f stdin

/mnt/nfs/nfs2/bickhart-users/binaries/blobtools/blobtools view -i usda_pacbio_pilon_blobplot.blobDB.json -r superkingdom -b -o usda_pacbio_supkingdom

grep 'tig00498169' usda_pacbio_species.usda_pacbio_pilon_blobplot.blobDB.table.txt               tig00498169     105788  0.4073  0       436.551 Pseudomonas virus phiKZ 443.0   0       tax0=Pseudomonas virus phiKZ:443.0; # A really cool large virus! Typically 200kb genome size
```


#### Blobplots on Illumina contigs

Let's see how the Illumina assemblies stand up compared to the PacBio reads.

```bash
mkdir illumina_usda_accumulated
sbatch --mem=15000 --nodes=1 --ntasks-per-node=1 -p assemble1 --wrap="cat illumina_usda_clusters/*.fasta > illumina_usda_accumulated/prefilter_usda_illumina_contigs.fa; samtools faidx illumina_usda_accumulated/prefilter_usda_illumina_contigs.fa; perl unique_fasta.pl illumina_usda_accumulated/prefilter_usda_illumina_contigs.fa.fai illumina_usda_accumulated/prefilter_usda_illumina_contigs.fa illumina_usda_accumulated/unique_usda_illumina_contigs.fa; samtools faidx illumina_usda_accumulated/unique_usda_illumina_contigs.fa; module load bwa; bwa index illumina_usda_accumulated/unique_usda_illumina_contigs.fa;"

```

#### Redoing pilon correction and Illumina assembly blobs

The Pacbio pilon correction was faulty because not all contigs were present in the cluster dataset. I'm going to fix that now.

```bash
# Downloading Illumina contigs
python ../../binaries/download_from_gdrive.py 1BEVP_ft3Tk24faPprYtnfW19EPgaHVP1 ./illumina_usda_accumulated/mick_megahit_final_full.fasta

sbatch --nodes=1 --mem=20000 -p assemble1 --ntasks-per-node=1 --wrap="bwa index mick_megahit_final_full.fasta"

sleep 4h; perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b aligns -t ../usda_illumina_fastas.tab -f mick_megahit_final_full.fasta -m -p assemble1
```

And the redo of pilon correction.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_retry_pilon

```bash
bwa index rumen_pacbio_multiround.asm.fasta

perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b aligns -t ../usda_illumina_fastas.tab -f rumen_pacbio_multiround.asm.fasta -m -p assemble1
```