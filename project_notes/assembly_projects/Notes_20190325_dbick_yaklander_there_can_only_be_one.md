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


## Running the repeat analysis

OK, now that everything is ready, let's queue up RepeatMasker and generate the files I need. I've already done repeatmasking on the ARS-UCDv1.2 assembly and the UMD3.1 assembly, so those comparative repeat lengths can be added to the plots if needed for comparison.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/yaklander
 
```bash
module unload perl/5.24.1
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.dam.gapfilled.arrow2.fasta"
sbatch --nodes=1 --mem=40000 --ntasks-per-node=25 -p medium --wrap="/beegfs/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 25 -q -species cow -no_is -gff yaklander.sire.gapfilled.arrow2.fasta"
```