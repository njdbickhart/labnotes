# Generating the red clover reference assembly
---
*12/14/2018*

These are my notes on generating the data, performing the assembly and providing some analysis for the red clover genome assembly project.

## Table of Contents


## Diagnostic assembly

I want to try to assemble a portion of the data (~ 15 Gbp) of the first batch of red clover nanopore reads for diagnosis of how much more data we need. I am going to queue up a full job and hopefully it will run to conclusion over the weekend.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
# Run with all nanopore reads, with parameters designed to reduce influence of systematic error overestimation
sbatch --nodes=1 --ntasks-per-node=20 --mem=64G -p medium --wrap="canu -p clover_hen_test -d clover_hen_test genomeSize=420m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' -nanopore-raw ./*/*.fastq"

# Note, the meryl routine had substantial problems with the number of filtered short reads in our dataset! 
# Some datasets had >90% filtered reads because of reads < 1000 bp. 
# It appears to be running still, since we're still above 30 X coverage for the genome size.
```