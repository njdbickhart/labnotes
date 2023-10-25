# Creating a genetic background for causal variant screening
---

## Table of Contents

## Generating alignment files for variant calling

I am testing quite a few things here, so I want to progress cautiously while I work out the kinks of the variant calling workflow. I will be using [octopus](https://luntergroup.github.io/octopus/docs/tutorials/germline) as that seems to provide the best quality to resources value from [a recent publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08365-3). 

> Anunna: /scratch/HG/bickha01/turkey_background/ 


```bash
module load SHARED/sra-toolkit
ln -s /lustre/shared/Liftover_Pipeline/Mgal_WUR_HG_1.0.fa Mgal_WUR_HG_1.0.fa

conda activate /lustre/backup/HG/bickha01/conda_envs/sraGeneticBackground
module load gatk/4.2.6.1

# Dry run test:
snakemake -s ~/python_toolchain/snakeMake/sraGeneticBackground/sraGeneticBackground.snk -n -r

# DAG printing
snakemake -s ~/python_toolchain/snakeMake/sraGeneticBackground/sraGeneticBackground.snk --forceall --rulegraph | dot -Tpdf > dag.pdf

sbatch -N 1 -n 1 --mem=5000 -p main -q std -t 5-0 --wrap='snakemake -s ~/python_toolchain/snakeMake/sraGeneticBackground/sraGeneticBackground.snk --cluster-config ~/python_toolchain/snakeMake/sraGeneticBackground/cluster.json --cluster "sbatch -N 1 -n {cluster.ntasks-per-node} --mem={cluster.mem} -p main -q std -J {cluster.jobname} -o logs/{cluster.output} -t {cluster.time}" -p --jobs 25 --verbose --latency-wait 40 -T 1'
```

### NOTE: All of the Wageningen samples were incorrectly submitted as "paired-end" without read pairing information. 

> Anunna: /scratch/HG/bickha01/chicken_background

```bash
perl -ne 'chomp; if($_ =~ /^>/){if($_ =~ /unlocalized/ || $_ =~ /unplaced/){@s = split(/\s+/, $_); print "$s[0]\n";}elsif($_ =~ /mitochondrion/){print ">mitochondrion\n";}else{($v) = $_ =~ /chromosome (.{1,2}),/; print ">chr$v\n";}}else{print "$_\n";}' <  GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna > GalGal1.mat.broiler.renamed.fasta

module load samtools
samtools faidx GalGal1.mat.broiler.renamed.fasta

samtools faidx GalGal1.mat.broiler.renamed.fasta chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29 chr30 chr31 chr32 chr33 chr34 chr35 chr36 chr37 chr38 chr39 chrW chrZ mitochondrion > GalGal1.mat.broiler.renamed.kary.fasta
samtools faidx GalGal1.mat.broiler.renamed.kary.fasta

perl -lane 'print "  \"$F[0]\" : [\n   \"P\",\n   \"$F[1]\"\n  ],";' < samples_raw.tab > samples_fmt.json
cp samples_fmt.json default.json

conda activate /lustre/backup/HG/bickha01/conda_envs/sraGeneticBackground
module load gatk/4.2.6.1

# Dry run test:
snakemake -s ~/python_toolchain/snakeMake/sraGeneticBackground/sraGeneticBackground.snk -n -r

sbatch -N 1 -n 2 --mem=5000 -p main -q std --wrap="picard CreateSequenceDictionary R=/scratch/HG/bickha01/chicken_background/GalGal1.mat.broiler.renamed.kary.fasta O=/scratch/HG/bickha01/chicken_background/GalGal1.mat.broiler.renamed.kary.dict"

sbatch -N 1 -n 1 --mem=5000 -p main -q std -t 5-0 --wrap='snakemake -s ~/python_toolchain/snakeMake/sraGeneticBackground/sraGeneticBackground.snk --cluster-config ~/python_toolchain/snakeMake/sraGeneticBackground/cluster.json --cluster "sbatch -N 1 -n {cluster.ntasks-per-node} --mem={cluster.mem} -p main -q std -J {cluster.jobname} -o logs/{cluster.output} -t {cluster.time}" -p --jobs 25 --verbose --latency-wait 40 -T 1 freebayes_only'
```

## Use of a new snpcalling workflow

> Anunna: /scratch/HG/bickha01/side_projects/turkey_eye

```bash
# to format all of the samples for inclusion in the default.json file
ls /scratch/HG/bickha01/side_projects/turkey_eye/raw_files/*.gz | perl -e 'use File::Basename; %samples; while($f = <STDIN>){chomp $f; $b = basename($f); @bsegs = split(/_/, $b); push(@{$samples{$bsegs[0]}}, $f);} foreach $k (keys(%samples)){print "    \"$k\" : [\n"; print "      [\"" . join("\",\n      \"", @{$samples{$k}}) . "\"],\n    ],\n";}'

conda activate /lustre/backup/HG/bickha01/conda_envs/sraGeneticBackground/

sbatch -N 1 -n 1 -p main -q std --mem=6000 -t 5-0 snakemake -s ~/python_toolchain/snakeMake/snpCalling/snpCalling --cluster-config ~/python_toolchain/snakeMake/snpCalling/cluster.json --cluster "sbatch -N 1 -n {cluster.ntasks-per-node} --mem={cluster.mem} -p main -q std -J {cluster.jobname} -o logs/{cluster.output} -t {cluster.time}" -p --jobs 20 --verbose --latency-wait 40
```