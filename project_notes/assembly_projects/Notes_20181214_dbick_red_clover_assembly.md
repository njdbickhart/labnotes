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


# KAT analysis of Diagnostic assembly

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