# Cattle assembly manuscript work
---
*1/29/2019*

These are my notes on work done to support the cattle assembly manuscript.

## Table of Contents


## Summary comparison statistics

I want to draw a few brief comparisons between assemblies in order to highlight the differences. My first priority is establishing the rank order correlation between the optical maps and the cattle assemblies.

### Rank order comparisons

The cattle optical maps were aligned using the Bionano proprietary software, I think. The alignment files are located here on Ceres:

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/OpticalMap

First order of business is to identify simple counts of alignments and blocks.

```bash
wc -l *.aln
   18674 OM-ARS-UCD1.2.full.aln
   18629 OM-UMD3.1.1.full.aln

# I think that Ben included the unplaced contigs as well, so this makes things very difficult.

cat OM-ARS-UCD1.2.full.aln | cut -f1 | sort | uniq | wc -l
18674	<- queries = the optical map contigs
cat OM-ARS-UCD1.2.full.aln | cut -f2 | sort | uniq | wc -l
39		<- reference = assembled chromosomes

# OK, this is a bad situation. The OM was split into windows prior to mapping
# So, I need to take into account the OM contig as well as the window numbers

# I wrote a script to generate rank order from these alignment files
# It uses the reference chromosomes as the basis for the preliminary sort, so it assumes that the assembly scaffolds are larger than the OM scaffolds
sbatch sortAndRankOMData.pl OM-ARS-UCD1.2.full.aln OM-ARS-UCD1.2.full.ranks OM-ARS-UCD1.2.full.om.consensus
sbatch sortAndRankOMData.pl OM-UMD3.1.1.full.aln OM-UMD3.1.1.full.ranks OM-UMD3.1.1.full.om.consensus
```

Now to check the rank order correlations.

```R
# ARS UCD rank order test
data <- read.delim("OM-ARS-UCD1.2.full.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 6.7575e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.3772741

# UMD rank order test
data <- read.delim("OM-UMD3.1.1.full.ranks", header=TRUE)
cor.test(data$RefRank, data$OMRank, alternative="two.sided", method="spearman", exact = TRUE)

        Spearman's rank correlation rho

data:  data$RefRank and data$OMRank
S = 9.9289e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho
0.07837418
```

These are atrocious for both assemblies! I could determine ref chromosome consensus order, but I'm not sure if that would even help at this point! I need to see how disjointed the alignments are between optical maps first. 