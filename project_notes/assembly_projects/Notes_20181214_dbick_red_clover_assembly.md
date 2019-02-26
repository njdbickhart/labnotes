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