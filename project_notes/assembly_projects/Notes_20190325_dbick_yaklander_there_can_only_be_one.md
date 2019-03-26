# Yaklander genome assembly analysis
---
*3/25/2019*

These are my notes on running some analysis on the Yaklander (Yak x Highland cross) triobinned assembly.

## Table of Contents


## Preparing the assembly fastas

I just need to gather the materials I need for the analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander

```bash
module load bwa
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=14000 --ntasks-per-node=1 -p short --wrap="bwa index yaklander.sire.gapfilled.arrow2.fasta"
```

 