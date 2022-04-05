# Sheep sample number 2
---
*2/15/2022*

These are my notes on the processing of the second HiFi sheep sample

## Table of contents


## Setup

Concatenating the Hi-C data for use in binning.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/sequence_data/sheep_hic

```bash
sbatch -N 1 -n 2 --mem=10000 -p priority -q msn -t 1-0 --wrap='gunzip -c *_R1_001.fastq.gz > combined_sheep2_hic.R1.fastq'
sbatch -N 1 -n 2 --mem=10000 -p priority -q msn -t 1-0 --wrap='gunzip -c *_R2_001.fastq.gz > combined_sheep2_hic.R2.fastq'
```

## Running the project

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/analysis/sheep2

```bash
module load python_3/3.6.6 miniconda/3.6 samtools

cp ~/python_toolchain/snakeMake/hifiMAGManuscript/Snakefile ./
cp -r ~/python_toolchain/snakeMake/hifiMAGManuscript/envs ./
cp -r ~/python_toolchain/snakeMake/hifiMAGManuscript/scripts ./
cp -r ~/python_toolchain/snakeMake/hifiMAGManuscript/cluster.json ./

# Note, I had to edit the Snakefile to add metabat binning rules and to add other necessary parameters

sbatch -N 1 -n 2 --mem=10000 -p priority -q msn -t 2-0 snakemake -s Snakefile --cluster-config cluster.json --cluster "sbatch -N 1 --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} -p priority -q msn -o {cluster.stdout} -t 2-0" -p --use-conda --jobs 250 --verbose --latency-wait 40
```