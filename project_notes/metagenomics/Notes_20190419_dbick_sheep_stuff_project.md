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

sbatch --nodes=1 --ntasks-per-node=4 --mem=15000 -p msn --wrap="bwa mem -t 3 asm.contigs.fasta sheep_poop_CCS.fastq.gz | samtools sort -T sheep.temp -o sheep_poop_ctg_ccs.sort.bam -"
```

## HiFi

I want to call variants on the data using the [PacBio HiFi](https://www.biorxiv.org/content/biorxiv/early/2019/01/23/519025.full.pdf) pipeline. I will be using minimap2 with the defaults recommended in this [github repo](https://github.com/PacificBiosciences/pbmm2).

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sheep_poop

```bash
module load minimap2

sbatch --nodes=1 --mem=15000 --ntasks-per-node=3 -p msn --wrap="minimap2 -a -k 19 -w 10 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 ccs_all/asm.contigs.fasta sheep_poop_CCS.fastq.gz > sheep_canu_normal.ccsalgn.sam"
sbatch --nodes=1 --mem=15000 --ntasks-per-node=3 -p msn --wrap="minimap2 -a -k 19 -w 10 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 asm.contigs.fasta sheep_poop_CCS.fastq.gz > sheep_canu_compressed.ccsalgn.sam"
sbatch --nodes=1 --mem=15000 --ntasks-per-node=3 -p msn --wrap="minimap2 -a -k 19 -w 10 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 ccs_flye/scaffolds.fasta sheep_poop_CCS.fastq.gz > sheep_flye.ccsalgn.sam"

# Now to sort all the files into bams
module load samtools; for i in sheep_canu_compressed.ccsalgn.sam sheep_canu_normal.ccsalgn.sam sheep_flye.ccsalgn.sam; do name=`echo $i | cut -d'.' -f1,2`; echo $name; sbatch --nodes=1 --mem=3000 --ntasks-per-node=1 -p short --wrap="samtools sort -T $name.tmp -o $name.sorted.bam $i"; done
```