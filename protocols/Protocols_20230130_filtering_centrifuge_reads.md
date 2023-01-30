# Filtering metagenome reads from Centrifuge
---
*01/30/2023*

These are my notes on creating an environment and testing my python script to filter reads from fastq files for use in metagenome assembly. 

## Setting up the environment

I tried to reduce the number of requirements for the script, but the python3 networkX module had too many useful features already built in. Let's first create the environment that will be used for running the script.

**NOTE: you only need to create the environment once! Afterwards you can activate it at any time to get the required dependencies! I will highlight this in the comment below**

> Ceres: /project/rumen_longread_metagenome_assembly/analysis/centrifuge

```bash
# Create interactive session to run shell commands
srun -N 1 -n 3 --mem=16000 -p priority -q msn -t 1-0 --pty bash

# Load the conda package manager
module load miniconda

# Creating the environment (ONLY NEEDS TO BE DONE THIS ONCE!)
conda create -p /project/rumen_longread_metagenome_assembly/analysis/centrifuge/fastq_filter networkx gzip

# Activating the environment (the only part of the proceedure that you need to do each time)
conda activate /project/rumen_longread_metagenome_assembly/analysis/centrifuge/fastq_filter
```

You will also need to download my [python toolchain](https://github.com/njdbickhart/python_toolchain) repository from github. Here is how to do so:

> Ceres: ~

```bash
# Go to your home folder if you aren't there already
cd ~

# Download the repo
git clone https://github.com/njdbickhart/python_toolchain.git
```

## Creating test datasets

I want to run a limited test on the data to debug my script. So I will subsection 1000 reads from the fastq and the ".cout" file.

```bash
# retrieving 1000 reads (4000 lines) from the fastq file
head -n 4000 /90daydata/rumen_longread_metagenome_assembly/rumen_pools_11cells.HiFi.fastq > test_data.fastq

# Actually, I just checked the cout file and the centrifuge output isn't ordered in the same way that the fastq file is ordered! This means that I will have to load the entire cout file each time or filter the reads. 
# OK, it takes more work, but let's filter the reads from the cout file
grep -P '^@m' test_data.fastq | perl -ne '$_ =~ s/^@//; print $_;' > reads.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f rumen_pools_11cells.cout -c 0 -l reads.list -d '\t' > temp.cout
head -n 1 rumen_pools_11cells.cout > header.cout

cat header.cout temp.cout > test.cout
rm header.cout temp.cout
```

There are several other essential files that were moved or copied to that directory. Here is a full listing of the required files:

* test_data.fastq (-f)
* test.cout (-c)
* nodes.dmp (-n)

## Assessing taxids for recursive filtering

The goal is to take a higher level taxid (ie. Phylum level) and use that to filter all taxonomic entries below it (ie. to Species level). Our filtration won't be perfect, but we want to take all entries below a certain taxid out of our dataset. I have written the script to allow you to filter multiple taxIDs, but it will also recursively eliminate everything below those selections. This way, you can be more precise by entering multiple taxIds for Genus- or Species-level entries.

The Centrifuge program produces "Kraken-style" reports (.kfile) that allow you to skim the phylogeny of assignments in your dataset quite easily. Using that file, I found a substantial proportion of the data was contained within the Chordata:

```
# %     Nreads  R@level Tax	 	Taxid		Taxonomic name
 36.25  5804190 209665  D       2759        Eukaryota
 31.58  5057049 0       -       33154         Opisthokonta
 30.24  4843163 91988   K       33208           Metazoa
 29.66  4749192 0       -       6072              Eumetazoa
 29.61  4741642 0       -       33213               Bilateria
 26.41  4229551 0       -       33511                 Deuterostomia
 26.39  4226442 24294   P       7711                    Chordata
 26.22  4199361 0       -       89593                     Craniata
 26.22  4199361 0       -       7742                        Vertebrata
 26.22  4199214 0       -       7776                          Gnathostomata
 26.21  4197170 0       -       117570                          Teleostomi
 26.21  4197170 0       -       117571                            Euteleostomi
 25.31  4052958 0       -       8287                                Sarcopterygii
 25.30  4051849 0       -       1338369                               Dipnotetrapodomorpha
 25.30  4051841 0       -       32523                                   Tetrapoda
 25.28  4047779 0       -       32524                                     Amniota
 24.94  3994460 10595   C       40674                                       Mammalia
 24.87  3982660 0       -       32525                                         Theria
 24.84  3977888 0       -       9347                                            Eutheria
 24.82  3973824 0       -       1437010                                           Boreoeutheria
 23.71  3797065 0       -       314145                                              Laurasiatheria
 23.60  3779961 0       -       91561                                                 Cetartiodactyla
 23.55  3770487 0       -       9845                                                    Ruminantia
 23.55  3770487 0       -       35500                                                     Pecora
 17.27  2766238 1178    F       9895                                                        Bovidae
 16.15  2586218 0       -       9963                                                          Caprinae
 16.13  2582340 146     G       9935                                                            Ovis
 16.05  2569531 979     S       37174                                                             Ovis canadensis
 16.04  2568552 2568552 -       112262                                                              Ovis canadensis canadensis
  0.08  12663   7602    S       9940                                                              Ovis aries
  0.03  5061    5061    -       9938                                                                Ovis aries musimon
```

Chordata represents a decent entry point for our filtration, and it should get rid of a large proportion of reads that are likely from host DNA here. 

## Running the script

Using the test files and our environment, we can now test out the script and filter the fastq file.

> Ceres: /project/rumen_longread_metagenome_assembly/analysis/centrifuge

```bash
# If it is not already active:
conda activate /project/rumen_longread_metagenome_assembly/analysis/centrifuge/fastq_filter

# running the script
python3 ~/python_toolchain/metagenomics/filterRecursiveFastqsTaxid.py -f test_data.fastq -t 7711 -c test.cout -n nodes.dmp -o test_output.fastq.gz
	Identified 118542 descendants from taxid 7711
	Total set entries comprise 118542 elements.
	Retained 944 from 1000 reads (94.40)%)
	Filtered 56 reads out of a total of 1000 (0.06)
```

## Running on the full dataset

I will run this in an interactive session, but it should be able to be queued using SBATCH as well.

```bash
srun -N 1 -n 3 --mem=48000 -p priority -q msn -t 1-0 --pty bash

conda activate /project/rumen_longread_metagenome_assembly/analysis/centrifuge/fastq_filter

python3 ~/python_toolchain/metagenomics/filterRecursiveFastqsTaxid.py -f /90daydata/rumen_longread_metagenome_assembly/rumen_pools_11cells.HiFi.fastq -t 7711 -c rumen_pools_11cells.cout -n nodes.dmp -o rumen_pools_11cells.HiFi.filt.fastq.gz
```