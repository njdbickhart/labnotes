# Sheep fecal project
---
*4/19/2019*

These are my notes from working on Tim's Sheep feces project.

## Table of Contents

## Starting out

Serge has prepared an assembly and was interested in some statistics on it. I will try to run CheckM on the contigs (even though they are not already binned). I will use the existing sheep CCS reads as a coverage dataset for downstream analysis.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sheep_poop

```bash
module load bwa checkm/v1.0.11 metabat/2.12.1 samtools

sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa index asm.contigs.fasta"

sbatch --nodes=1 --ntasks-per-node=4 --mem=15000 -p short --wrap="bwa mem -t 3 asm.contigs.fasta sheep_poop_CCS.fastq.gz | samtools sort -T sheep.temp -o sheep_poop_ctg_ccs.sort.bam -"
```