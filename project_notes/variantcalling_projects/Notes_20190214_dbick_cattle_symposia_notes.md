# ADSA Cattle Symposia Talk Notes
---
*2/14/2018*

These are my notes and commands for generating datasets for use in my talk at the 2019 ASDA Cattle Genetics symposia.


## Table of Contents

## Lit review and structure

### Genetic and physical mapping

[James Womack](https://genome.cshlp.org/content/15/12/1699.full.html) has a great review on the progress up until the early 00's. The first gene map was by [Heuertz et al](https://www.karger.com/Article/Abstract/131601) but [Womack](https://www.ncbi.nlm.nih.gov/pubmed/3082971) expanded this map greatly by using cattle-hamster hybrid somatic cell lines in 1986. In-situ hybridization was added to cattle-hamster hybrid somatic cell lines to make a comprehensive map of [many cattle genes](https://link.springer.com/article/10.1007/BF00296815) in the early 90's. A [multi-national group](https://www.nature.com/articles/ng0394-227) and [USDA-MARC](https://www.genetics.org/content/136/2/619.abstract?ijkey=c8173747b99141827e87e12dd6281489f54c1196&keytype2=tf_ipsecsha) published a linkage map almost simultaneously in 1994 spanned 2464 centimorgans of the cattle genome. **NOTE: the USDA-MARC map could be fun to plot on ARS-UCDv1.2.** The USMARC map used 280 out of 468 microsatellites to develop linkage maps for cattle. Next, Womack led the development of a [radiation-hybrid map](https://link.springer.com/article/10.1007%2Fs003359900593?LI=true) for cattle. 

James Womack predicted that the bovine genome [would never be sequenced](https://www.annualreviews.org/doi/pdf/10.1146/annurev-animal-020518-114902) in the 90's, in order to justify grant proposals for comparative mapping. 

### Molecular basis of testing

The earliest markers used by cattle geneticists were probed southern blots and PCR-based markers. Term "satellite" comes from [Saul Kit](https://www.sciencedirect.com/science/article/pii/S0022283661800752?via%3Dihub) who found bands from CsCl gradient centrifugation. Microsatellites are 1 to 10 bp repetitive, low complexity sequence. Useful because they are multi-allelic and are frequently found outside of coding regions (no positive selection). 

### Reduced representation sequencing for SNP50 (Van Tassel Nat. Methods)

Curt, Tad and company first used DraI (TTT^AAA) in an *in silico* digest of Btau3.1 but found that it did not perform as well *in vitro*. This may be due to their use of the incomplete assembly (Btau3.1 = 2.4 Gbp assembly size; 

### Genome assemblies

All [Baylor-led](https://www.hgsc.bcm.edu/other-mammals/bovine-genome-project) genome assembly projects are hosted online and have links to the NCBI repositories to download the data. 

## Assessing relative genome resolution of different technologies

I want to draw comparative maps of the different technologies used to track genetic variants in the cattle genome. I intend this to visualize the progress that each method has had on the resulting assembly map and to show the differences in resolution that each method provides.

## Comparing improvements of genome assemblies throughout the years

I see the following figures as being integral to the presentation and subsequent manuscript:

1. A Circos plot showing marker locations on the assembly from physical maps
2. A histogram showing the assembly size, and continuity of the assemblies, from Btau1.0 to ARS-UCDv1.2
3. A Chord diagram showing the differences between Btau1.0 to ARS-UCDv1.2