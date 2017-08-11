# SNP probe remapping on cattle assemblies
---
*8/10/2017*

These are my notes for remapping UMD3-based SNP probes onto the new cattle assembly. My hope is to do as thorough of a job as possible and to resolve >99% of marker locations, if possible. 

## Table of contents


## Preparing reference fastas and generating metadata

First, let's talk about what I need. I want to approach this from a hierarchical standpoint, where I:

* Identify rearranged segments of the genome using liftover
* Confirm simple rearrangements using a fast aligner (BWA MEM)
* Perform specialized local alignment using a more precise alignment tool to resolve recalcitrant probes

I suspect that 98% of markers will be resolved by the fast aligner, whereas the remaining 2% will require far more care. If needed, I may even need to go to the LD information on the marker to try to find regions of sequence to align to in a very targeted approach!

OK, first things first, I need to generate the requisite data for liftover. In order to assist with this, I am going to use RepeatMasker to reduce the complexity of the genome (I need the repeat information for our subsequent analysis in any case, so this will be killing two birds with one stone). Reducing the complexity of the genome will make the blat alignment much faster by eliminating ambiguous alignments. 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/RepeatMasker/RepeatMasker -pa 40 -q -species cow -no_is -dir ars_ucd_14_igc_rmask ARS-UCD1.0.14.clean.wIGCHaps.fasta