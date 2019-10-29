# Dorper White Romanov cross assembly
---
*10/28/2019*

These are my notes on assisting with the Dorper x White Romanov assembly project.

## Table of contents


## Generating trio canu corrected reads

I first want to bin all of the reads from Tim's pacbio run.

> Ceres: /project/forage_assemblies/assemblies/romanov_whitedorp

```bash
module load canu/1.8 samtools
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -q msn -p msn --wrap="canu -p dorp_wroman -d dorp_wroman genomeSize=2800m -haplotypeSire /project/gaur_genome_assembly/WhiteDorper_x_Romanov/Dorper_sire.fastq  -haplotypeDam /project/gaur_genome_assembly/WhiteDorper_x_Romanov/Romanov_dam.fastq corMhapSensitivity=normal corOutCoverage=200 saveReadCorrections=true 'gridOptions=-q msn -p msn' -pacbio-raw /project/gaur_genome_assembly/WhiteDorper_x_Romanov/pacbio_data/all_reads.fq"
```