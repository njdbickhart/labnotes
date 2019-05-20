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