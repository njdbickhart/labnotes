# Dorper White Romanov cross assembly
---
*10/28/2019*

These are my notes on assisting with the Dorper x White Romanov assembly project.

## Table of contents


## Generating trio canu corrected reads

I first want to bin all of the reads from Tim's pacbio run.

> Ceres:

```bash
module load canu/1.8 samtools
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -q msn -p msn --wrap="canu -p dorp_wroman -d dorp_wroman genomeSize=2800m corMhapSensitivity=normal corOutCoverage=200 saveReadCorrections=true 'gridOptions=-q msn -p msn' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"
```