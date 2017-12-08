# Haplotype variant discovery and filtering
---
*12/6/2017*

This is my notes file for generating haplotype regions for full scale variant calling and analysis.

## Table of contents


## Generating initial data files and testing code

I need to convert the files from UMD3.1 coordinates to the new ARS-UCDv1.0.14 assembly. First, let's try a quick work-around with the recombination map. I plan to take recmap segment mappings to locate the "ends" of existing haplotypes and then check to see if they are contiguous or on the same chromosome. I will split those regions into separate categories, then visually inspect them for fidelity.

> Assembler2: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/haplotype_segs

```bash
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a ../ARS-UCD1.0.14.clean.wIGCHaps.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/rcmap_manifest_correct.sorted.fa -o ../ARS-UCD1.0.14.clean.wIGCHaps.recmap

# I cleaned up some segments that were small misassemblies present -- but eventually fixed -- in this version
vim ARS-UCD1.0.14.clean.wIGCHaps.recmap.cor.segs

```