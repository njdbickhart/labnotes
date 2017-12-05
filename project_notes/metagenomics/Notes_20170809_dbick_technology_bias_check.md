# Rumen microbiome technology bias check
---
*8/9/2017*

These are my notes and commands for testing the bias of different sequencing platforms in the collection and determining if there are major limitations to each.

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

# Base Vegan plotting polymorhpic functions did not work on the data. Had to hack out the NMDS coords
plot(scores(ord), col=c("red", "blue", "blue", "green", "red", "purple", "purple", "purple", "purple"), pch= c(15, 15, 16, 15, 16, 16, 17, 15, 15))
legend("bottomleft", c("UKNano", "USIllum", "Hess", "Manure", "USNano", "PBCCS", "PBRSII", "PBCheryl", "PBTim"), col=c("red", "blue", "blue", "green", "red", "purple", "purple", "purple", "purple"), pch= c(15, 15, 16, 15, 16, 16, 17, 15, 15))
dev.copy2pdf(file="combined_profile_nmds_vegan.pdf", useDingbats=FALSE)
```

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

