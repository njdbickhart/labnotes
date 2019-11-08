# Vetch Genome assembly
---
*5/20/2019*

## Table of Contents

## Testing the variance of the different Vetch libraries

OK, so Lisa tried out several different methods for diluting the DNA preps prior to library creation, and we need to see if there are any obvious biases in these methods that could cause us to choose to avoid them in the future!

We're just going to do the Mash screen to compare them (and to compare them against a red clover sample as an outgroup!).

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/red_clover_nano

```bash
# Creating a merged fastq file
for i in VetchBooger VetchMN VetchVortex VetchZymo; do echo $i; cat $i/*/fastq_pass/*.fastq > $i/$i.combined.pass.fastq; done

for i in VetchBooger VetchMN VetchVortex VetchZymo; do echo $i; sbatch --nodes=1 --mem=10000 --ntasks-per-node=4 -p msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -o $i/$i.combined.msh -p 4 -s 100000 -r -m 4 -g 420M $i/$i.combined.pass.fastq"; done
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste vetch_combined VetchBooger/VetchBooger.combined.msh VetchMN/VetchMN.combined.msh VetchZymo/VetchZymo.combined.msh VetchVortex/VetchVortex.combined.msh clover14.msh
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -t vetch_combined.msh VetchBooger/VetchBooger.combined.msh VetchMN/VetchMN.combined.msh VetchZymo/VetchZymo.combined.msh VetchVortex/VetchVortex.combined.msh clover14.msh > vetch_combined.dist

perl -ne '$_ =~ s/(Vetch.+)\/.+fastq/$1/g; print $_;' < vetch_combined.dist > vetch_combined.rfmt.dist

```

```R
library(dplyr)
library(ggplot2)
library(reshape)


vdist <- read.delim("vetch_combined.rfmt.dist", header=TRUE)
rownames(vdist) <- vdist$query
vdist <- vdist[,c(2:6)]
vdist.m <- melt(as.matrix(vdist))

pdf(file="vetch_distance_heatmap.pdf", useDingbats=FALSE)
ggplot(vdist.m, aes(X1, X2)) + geom_tile(aes(fill=value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") + labs(title = "Vetch Prep Dataset Mash Distances")
dev.off()
```

## Starting the assembly process

I'm going to unpack what I have and run an assembly over the long weekend.

> Ceres: /home/derek.bickharhth/forage_assemblies/sequence_data

```bash
for i in Vetch*.tar.gz; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="tar -xvf $i"; done

# Now to concatenate the fastqs
sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p msn --wrap='for i in Vetch*/*/fastq_pass/*.fastq; do cat $i; done > vetch_combined_reads.fastq'
```

Now to run Flye on Vetch to test it out.

> Ceres: /project/forage_assemblies/assemblies/vetch

```bash
module load miniconda
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=300000 --ntasks-per-node=70 -p msn flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye

## Ran out of memory! Going to try to increase a slight bit, but if it fails I'll have to try an alternate strategy
sbatch --nodes=1 --mem=370000 --ntasks-per-node=70 -p msn flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye --resume

## That also ran out of memory. Going to queue up on the MEM nodes for the long haul
sbatch --nodes=1 --mem=900000 --ntasks-per-node=70 -p mem flye --nano-raw /home/derek.bickharhth/forage_assemblies/sequence_data/vetch_combined_reads.fastq -g 2000m -t 70 -i 2 -o vetch_flye --resume

## I am going to try a co-assembly approach with the Noble vetch reads as well
sbatch --nodes=1 --mem=900000 --ntasks-per-node=70 -p mem -q memlimit flye --nano-raw /project/forage_assemblies/sequence_data/combined_both_vetch.fastq -g 2000m -t 70 -i 2 -m 10000 --asm-coverage 40 -o vetch_limit_flye
```

Now trying Canu

```bash
module load canu/1.8 java/1.8.0_121

# De novo correction
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short --wrap="canu -correct -p vetch_pre -d vetch_pre genomeSize=2000m corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 saveReadCorrections=true 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"

### NOT RUN YET replace fastq###
sbatch --nodes=1 --ntasks-per-node=30 --mem=10G -p short --wrap="canu -p vetch_canu -d vetch_canu genomeSize=2000m correctedErrorRate=0.120 'corMhapOptions=--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14' 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"

# No correction; assume "polyploid"
sbatch --nodes=1 --ntasks-per-node=10 --mem=10G -p short --wrap="canu -p vetch_nocorrect -d vetch_nocorrect genomeSize=2000m corMhapSensitivity=normal corOutCoverage=200 'batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50' saveReadCorrections=true 'gridOptions=-p short' -nanopore-raw /project/forage_assemblies/sequence_data/vetch_combined_reads.fastq"
```