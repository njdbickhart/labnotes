# Rumen Metagenome CCS assembly
---
*4/27/2020*

This is a new take on generating a rumen metagenome assembly with CCS reads.

## Table of Contents

## Assembly

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/rumen_ccs

```bash
module load miniconda/3.6

conda activate /KEEP/rumen_longread_metagenome_assembly/flye/

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p priority -q msn -t 4-0 flye -g 1.0g --pacbio-hifi /lustre/project/rumen_longread_metagenome_assembly/sequence_data/rumen_CCS/m54337U_200422_195119.Q20.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/rumen_CCS/m54337U_200525_191807.Q20.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/rumen_CCS/m54337U_200527_075128.Q20.fastq /lustre/project/rumen_longread_metagenome_assembly/sequence_data/rumen_CCS/m54337U_200627_185709.Q20.fastq -t 70 --meta -o flye_rumen_css

```