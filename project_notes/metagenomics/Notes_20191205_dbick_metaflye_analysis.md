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


## Running the analysis

#### CCS read alignments

> Ceres: /project/forage_assemblies/sheep_project

```bash
module load minimap2/2.6
for i in *.tar.gz; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="tar -xvf $i"; done

# OK, so the contigs are not properly binned, but given their sizes, this is probably OK.
# Let me start the long-term alignment first so that I can analyze this next week
# I will use this for blobtools "coverage" estimates as well
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch --nodes=1 --mem=30000 --ntasks-per-node=3 -p msn --wrap="minimap2 -x map-pb $i /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq.gz > $i.ccs.paf"; done

# The blobtools coverage estimates didn't work very well with the CCS reads in sam format. Maybe I can resolve this with scripting?
```