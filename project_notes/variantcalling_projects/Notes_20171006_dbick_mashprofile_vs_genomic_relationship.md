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


I have generated data from genosim and will use the following chromosomes from ARS-UCDv14 as a "seed" for the program:

|Chr|size|Newname|
|:---|---:|:----|
|1|158915630|1|
|12|87180481|2|
|28|45965725|3|

While I wait for Paul to test my version of his markersim program, I will try to alter the chromosome.data file for use in this program:

> assembler2: /mnt/nfs/nfs2/bickhart-users/binaries/aip_genosim

```bash
perl -lane 'if($F[0] == 1 || $F[0] == 12 || $F[0] == 28){print $_;}' < /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta.fai > test_run_chrs.fai

python3 ../python_toolchain/simulations/markersimChrEvents.py -c chromosome.data -f test_run_chrs.fai -o expandedchr.data

# OK, so I understand that the markersim program generates a subset of true variants vs false
# The full set of markers (5000 in my test run) is downsampled to provide the lower density chip estimates
# Let's increase the size of the simulated marker pool to a variant every 5 kb or so: 
./pedsim
./markersim
python3 ../python_toolchain/simulations/markersimChrEvents.py -c chromosome.data -f test_run_chrs.fai -o expandedchr.data
./genosim
```

Now I have the pedigrees, the genotypes and the chromosome.data table to proceed. I just need to write a program to generate the sequence data from this input file. 