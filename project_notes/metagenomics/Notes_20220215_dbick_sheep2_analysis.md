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

## Generating information on the Phase MAGs

> Ceres: /project/rumen_longread_metagenome_assembly/analysis/sheep2

```bash
module load miniconda snakemake/6.12.1

git clone https://github.com/WatsonLab/MAGpy.git

# Had problems with ete3 install
conda activate /project/rumen_longread_metagenome_assembly/environments/ete3
conda env export > ete3_replacement.yml
conda activate /project/rumen_longread_metagenome_assembly/environments/sourmash
conda env export > ../sourmash_replacement.yaml

cd MAGpy
cp ../ete3_replacement.yml envs/ete3.yaml
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="snakemake -rp -s MAGpy --cores 1 --use-conda test"
sbatch -N 1 -n 48 -p priority -q msn --mem=300000 --dependency=afterok:7763667 --wrap="rm -rf magpy_dbs; snakemake -rp -s MAGpy --cores 48 --use-conda setup"

# copying mags
for i in /project/rumen_longread_metagenome_assembly/assembly/sheep2/microbial_genomes/smith_sheep_63_new_rerun_clusters/*.fasta; do name=`basename $i | cut -d'.' -f1`; echo $name; cp $i mags/$name.fa; done

sbatch -N 1 -n 2 --mem=9000 -p priority -q msn -t 8-0 snakemake -rp -s MAGpy --use-conda all --cluster "sbatch -N 1 --ntasks-per-node=16 --mem=100000 -p priority -q msn -t 8-0" --jobs 100
```

## Running blobtools on the full original sheep datast

> Ceres: /project/rumen_longread_metagenome_assembly/assembly/sheep1

```bash
conda activate /project/rumen_longread_metagenome_assembly/environments/blobtools

sbatch -N 1 -n 72 --mem=300000 -p priority -q msn --wrap="diamond blastx --threads 72 --query edges.fasta --max-target-seqs 1 --db /project/rumen_longread_metagenome_assembly/binaries/blobtools/data/uniprot_referenceproteomes_202202.dmnd --evalue 1e-25 --outfmt 6 --tmpdir /project/rumen_longread_metagenome_assembly/assembly/sheep1/temp  --out edges.vs.uniprot.out"
```