# Simulation test of Mash profile vs normal G matrix
---
*10/6/2017*

## Table of Contents


## Setting up the test environment

Here are the components of the selection that I am hoping to emulate:

* A pedigree system that generates random matings from a set of founders
* Haplotype tracking via mating simulations
* Two autosomes and one sex chromosome (to simplify things and allow for far more comparisons over time)
* Output of low density and HD markers (both fixed sites; low density = 1 per haplotype and HD = several per haplotype)
* Output of sequence reads directly to MASH for sketch creation

In order to simulate the reads, I need a "hard-point" low level program to sample my genome. Here are some ideas to fill the gap while I consider what I want:

Pedigree simulation:
* AIPL genosim
* Plink (doesn't actually output markers and locations)

WGS variant generation:
* seqan
