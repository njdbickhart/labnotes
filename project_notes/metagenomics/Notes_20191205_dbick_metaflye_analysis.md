# MetaFlye analysis of sheep dataset
---
*12/5/2019*

These are my notes and plans for validating the MetaFlye data for use in the reviewer rebuttals.

## Table of Contents


## Analysis plan

I will need to divy up tasks and run them in order to generate the results needed.

* Run CheckM on both assemblies
* Run BlobTools on both assemblies
	* Taxify contigs and derive contig taxonomic profiles
	* Derive a consensus taxonomic affiliation per bin and plot CheckM results
* Run Minimap2 using CCS reads on both assemblies
* Align Hi-C reads to both assemblies
	* Run the viral association script
	* Count number of split-read or hard-clipped CCS reads
	* Adapt DESMAN to pick strains from the alignments