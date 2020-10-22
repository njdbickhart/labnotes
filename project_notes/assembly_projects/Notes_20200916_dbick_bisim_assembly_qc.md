# Bison Simmental analysis 
---
*9/16/2020*

These are my notes on quality control and publication of the Bison x Simmental manuscript

## Table of Contents


## Setting up analysis pipeline

I am going to queue up these assemblies for analysis in separate sets. The bison and simmental in separate folders.

#### Bison

The Bison's name is Woody.

> Ceres: /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_x_simmental_qc/bison_qc

```bash
module load miniconda/3.6
cp ~/python_toolchain/snakeMake/assemblyValidation/default.json ./

ls /lustre/project/cattle_genome_assemblies/bison_x_simmental/illumina_data/Woody/*.fastq.gz | perl -lane 'print "    \"$F[0]\",";'

vim default.json

sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 8-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 8-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```

#### Simmental

And the Simmental's name is Hollary (cow)

> Ceres: /lustre/project/cattle_genome_assemblies/bison_x_simmental/bison_x_simmental_qc/simmental_qc

```bash
module load miniconda/3.6
cp ~/python_toolchain/snakeMake/assemblyValidation/default.json ./

ls /lustre/project/cattle_genome_assemblies/bison_x_simmental/illumina_data/Hollary/*.fastq.gz | perl -lane 'print "    \"$F[0]\",";'

vim default.json
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn -t 5-0 snakemake --cluster-config ~/python_toolchain/snakeMake/assemblyValidation/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout} -t 5-0" -p --jobs 250 -s ~/python_toolchain/snakeMake/assemblyValidation/assemblyValidation --use-conda
```