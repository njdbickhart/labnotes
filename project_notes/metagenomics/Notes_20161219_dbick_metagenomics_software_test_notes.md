# Testing Metagenomics datasets with new methods
---
*12/19/2016**

These are my notes on tinkering around on previously generated metagenomics and ribotyping experiments. I will try to slowly advance towards the use of new techniques on the data so that I can garner some publications in the field.

## Table of Contents
* [The Austrian Acidosis dataset summary](#austriansummary)
* [Meta-analysis of metagenomics WGS data](#metanalysis)
	* [MASH sketching](#mash)
	* [MetaSpades assembly](#metaspades)


<a name="austriansummary"></a>
## The Austrian Acidosis dataset summary

I found this dataset on [SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB9353) and thought it would be a good test case for the following reasons:

* The Austrians did not label the data well
* They tested the epimural layer
* They tested different animal timepoint stages
* They considered animals in acidosis

Here is their project summary:

> This study aimed to investigate the impact of diet induced subacute rumen acidosis (SARA) in a transient feeding model on the ruminal bacterial epithelial microbiome in dairy cattle. Eight dry dairy cows were fed high-grain diet, mainly consisting of barley grain, wheat, corn, and rapeseed meal. Baseline diet was a forage-only diet consisting of hay and grass silage Feeding model was conducted as follow: one week adaptation (10% daily from 0 to 60%, on dry matter basis), one week SARA challenge, one week only forage-feeding and two weeks SARA challenge again (this time using a two-days adaptation phase). DNA was isolated from the rumen mucosa and used for Illumina MiSeq 16S rRNA gene amplicon sequencing.

My goal is to: 

1. Perform a cursory ribotyping experiment using their data
2. Use principal component analysis to isolate putative timepoint layers

It should be a straightforward exercise, but the differentiation of the different samples will be difficult as they have two SARA challenges, and I don't know how well they performed their sample collection.

Let's collect the data from a list of accessions that I've gathered from NCBI.

> pwd: /home/dbickhart/share/metagenomics/austrian_metagenomics_test_resources

```bash
# Test download
~/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --split-files ERR920940

# There was only one fastq file... the Austrians  must have combined their reads from the MiSeq into a long, overlapping read
# Getting all the rest
for i in `cat SraAccList.txt`; do echo $i; ~/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --split-files $i; done

# Moving them all to a separate folder
mkdir fastqs
mv ./*.fastq ./fastqs/
```

<a name="metanalysis"></a>
## Meta-analysis of metagenomics WGS data

I think that there are a number of studies on SRA that consist of WGS experiments on cattle rumen. I am hoping to download just the cattle data, process it through [MASH](https://github.com/marbl/Mash) and then identify inter-individual variability/WGS coverage bias using the dataset.

I have selected a large group of WGS experiments from SRA, and I have removed mischaracterized datasets that point to chicken, pig or other ruminants. The list of files is in SharedFolders/metagenomics/metanalysis/sra_file_accession_list.csv


<a name="mash"></a>
### MASH sketching 

I am now going to download the list of fastqs from SRA to perform the downstream analysis.

> Fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/datasources

```bash
# sending off the list of accession entries to fastq-dump
cat accession_list.txt | xargs -I {} sbatch --nodes 1 --tasks-per-node 1 --mem 500 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/sratoolkit.2.8.1-centos_linux64/bin/fastq-dump.2 -I --split-files {}"
```

Now I need to try to generate Mash sketches. First, I'd like to see how kmer size impacts MASH sketch attributes. Let's test this out on a range of sizes.

```bash
mkdir sketch_test
#noise filter = 2, 15mer
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 15 -r -m 2 -o sketch_test/DRR017219_15 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"

# std output follows
Sketching datasources/DRR017219_1.fastq...
Sketching datasources/DRR017219_2.fastq...
Estimated genome size: 2.1573e+08
Estimated coverage:    7.202
Estimated genome size: 2.29383e+08
Estimated coverage:    7.105
Writing to sketch_test/DRR017219_15.msh...

mv slurm-653995.out ./sketch_test/DRR017219_15.out

sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 18 -r -m 2 -o sketch_test/DRR017219_18 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -r -m 2 -o sketch_test/DRR017219_21 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 24 -r -m 2 -o sketch_test/DRR017219_24 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 32 -r -m 2 -o sketch_test/DRR017219_32 datasources/DRR017219_1.fastq datasources/DRR017219_2.fastq"


# In order to estimate the distance between sketches, they must be the same kmer size
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 15 -r -m 2 -o sketch_test/DRR019503_15 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 18 -r -m 2 -o sketch_test/DRR019503_18 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -r -m 2 -o sketch_test/DRR019503_21 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 24 -r -m 2 -o sketch_test/DRR019503_24 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"
sbatch --mem=20000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 32 -r -m 2 -o sketch_test/DRR019503_32 datasources/DRR019503_1.fastq datasources/DRR019503_2.fastq"

```

OK, I have an email from Serge that gives me good advice for clustering the data I've gathered.

> 1. You can use the sketch here:
RefSeqSketchesDefaults.msh.gz
However, if you are comparing mash sketches of a metagenome to a reference it isn’t going to do the right thing. It is using a similarity score (so full sequence represented by metagenome compared to full refseq genome) whereas you want a contains version instead. This sketch should work fine for clonal samples. This has been on our todo list but we haven’t finished it.
>
>2. The option is -r -m 2 (or a larger number) which is an exact filter vs the bloom filter and does a better job removing noise.
>
> 3. We’ve generally used sketch = 10000, kmer=21 for clustering Illumina metagenomic datasets.


I need to sketch all of the data files, then get a good outgroup and then do a dist comparison.

> fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects

```bash
mkdir meta_sketches
for i in `cat datasources/accession_list.txt`; do echo $i; sbatch --mem=10000 --nodes=1 --ntasks-per-node=5 --wrap="../binaries/mash-Linux64-v1.1.1/mash sketch -p 5 -k 21 -s 10000 -r -m 2 -o meta_sketches/${i}_21mer_10k datasources/${i}_1.fastq datasources/${i}_2.fastq"; done

# Combine meta_sketches
ls meta_sketches/*.msh > combined_meta_sketchs.list
../binaries/mash-Linux64-v1.1.1/mash paste -l meta_sketches/combined_profile combined_meta_sketchs.list

# Now to automate Mash distance estimation for each sketch
for i in `ls *_10k.msh`; do echo $i; sbatch --mem=8000 --nodes=1 --ntasks-per-node=3 --wrap="../../binaries/mash-Linux64-v1.1.1/mash dist -p 3 -t combined_profile.msh $i > $i.dist"; done

# now to generate a phylip formatted distance matrix
perl ~/sperl/metagenomics_scripts/convertMashToPhylip.pl -f meta_sketches -o meta_sketches/phylip_dist.matrix

# Finally, I need to learn how to generate a phylip tree!
module load phylip/3.696
neighbor

# The "outtree" file is now useable for plotting trees
```

```R
library(ape)
library(Biostrings)
library(ggtree)
tree <- read.tree(file="outtree")
ls <- data.frame(x1=0, x2=0.5, y1=0, y2=0)
ggtree(tree) + geom_segment(data=ls, aes(x=x1, xend=x2, y=y1, yend=y2))

# the tree looks good, but I want to also generate a heatmap
library(gplots)
matrix_data <- read.delim("phylip_dist.matrix", header=FALSE, sep=" ")
library(dplyr)
library(tidyr)

# Now to reformat the dataframe to remove some of the artifacts needed by phylip
samplenames <- filter(matrix_data, V1 != 170) %>% select(V1)
samplenames <- samplenames[,1] # To make it into a vector

matrix_data.format <- filter(matrix_data, V1 != 170) %>% select(-V1)
row.names(matrix_data.format) <- samplenames
colnames(matrix_data.format) <- samplenames

library(colorspace)
library(dendextend)
library(stringr)

# Now to assign colors based on sample categories
metadata <- read.delim("sra_file_accession_list.csv", header=TRUE, sep=",")
metadata.categories <- select(metadata, Run, Platform, Tissue, Sampling.Method, Breed)
summary(metadata.categories)
metadata.categories$Abbv <- str_sub(metadata.categories$Run, start=-8)
# Getting the replication count
abbr_run_names <- data.frame(name=str_sub(samplenames, start=-10, end=8))
abbr_run_names <- group_by(abbr_run_names, name)
summarize(abbr_run_names, n())
colnames(abbr_run_names) <- c("Abbv", "Count")
metadata.categories.final <- na.omit(left_join(filter(metadata.categories, Platform == "ILLUMINA"), abbr_run_names, by="Abbv"))

# I had  not included the Ion torrent and 454 data into the dataset, so these were removed
# Note that dplyr throws stuff into a new "table" format object that does not respond to R's base summarize function

# Picking the tissue type first
# Duplicates the colors because each sample is duplicated in paired end fashion!
cat_col <- rep(rev(rainbow_hcl(6))[as.numeric(metadata.categories.final$Tissue)], metadata.categories.final$Count)

# Testing clustering methods for correlation
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
iris_dendlist <- dendlist()
for(i in seq_along(hclust_methods)){ hc_meta <- hclust(as.dist(matrix_data.format), method = hclust_methods[i]); iris_dendlist <- dendlist(iris_dendlist, as.dendrogram(hc_meta))}
names(iris_dendlist) <- hclust_methods

# Now to check the correlation of clustering methods
meta_dendlist_corr <- cor.dendlist(iris_dendlist)
library(corrplot)
corrplot(meta_dendlist_corr, "pie", "lower")
dev.copy2pdf(file="hclust_meta_corrplot.pdf", useDingbats=FALSE)

# Average and Ward.D look quite different, let's plot a tanglegram to see the differences
iris_dendlist %>% dendlist(which = c(1,4)) %>% ladderize %>% set("rank_branches") %>% tanglegram(common_subtrees_color_branches = TRUE) # 1 is Ward.D and 4 is Average
dev.copy2pdf(file="ward_vs_avg_meta_clusttanglegram.pdf", useDingbats=FALSE)

# Plotting all methods to compare
par(mfrow = c(2,4))
for(i in 1:8){iris_dendlist[[i]] %>% set("branches_k_color", k=2) %>% plot(axes= FALSE, horiz = TRUE); title(names(iris_dendlist)[i])}

# Finally, let's try to compare these methods against the number of clusters we hoped to get
# Note: the clusters assume that the order of the data is the cluster order -- not the case here!
get_ordered_5_clusters <- function(dend){cutree(dend, k=5)[order.dendrogram(dend)]}
dend_5_clusters <- lapply(iris_dendlist, get_ordered_5_clusters)
compare_clusters_to_meta <- function(clus){FM_index(clus, rep(1:5, c(36, 26, 2, 34, 72)), assume_sorted_vectors = TRUE)}
clusters_performance <- sapply(dend_5_clusters, compare_clusters_to_meta)
dotchart(sort(clusters_performance), xlim = c(0.4,1), xlab = "Fowlkes-Mallows Index", main = "Performance of clustering algorithms on metagenomics hashes", pch = 19)

# So now that I have the routine, let's try a better reassignment
library(hashmap)
catLookup <- hashmap(c("Rumen", "Foot", "Manure", "Nasal", "Oral"), c(1,2,3,4,5))
ordered_cats <- rep(metadata.categories.final$Tissue, each=2)
get_ordered_5_clusters <- function(dend){catLookup[[ordered_cats[order.dendrogram(dend)]]]}
dend_5_clusters <- lapply(iris_dendlist, get_ordered_5_clusters)
clusters_performance <- sapply(dend_5_clusters, compare_clusters_to_meta)
dotchart(sort(clusters_performance), xlim = c(0.2,1), xlab = "Fowlkes-Mallows Index", main = "Performance of clustering algorithms on metagenomics hashes", pch = 19)

#So the "single" method is the best
dev.copy2pdf(file="fm_index_comparison_methods.pdf", useDingbats=FALSE)

## Generating the heatmap

singleClust <- hclust(as.dist(matrix_data.format), method="single")
singleDend <- as.dendrogram(singleClust)
labels_colors(singleDend) <- rainbow_hcl(5)[get_ordered_5_clusters(singleDend)]
singleDend <- color_branches(singleDend, k=5)

library(gplots)
some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
gplots::heatmap.2(as.matrix(matrix_data.format), main = "Metagenome Tissue MASH Clustering", dendrogram = "row", Rowv = singleDend, Colv = "NA", trace = "none", denscol = "grey", density.info = "density", RowSideColors = rev(labels_colors(singleDend)), col = some_col_func)
dev.copy2pdf(file="MASH_heatmap.pdf", useDingbats=FALSE)
gplots::heatmap.2(as.matrix(matrix_data.format), main = "Metagenome Tissue MASH Clustering", dendrogram = "none", Rowv = "NA", Colv = "NA", trace = "none", denscol = "grey", density.info = "density", RowSideColors = rev(labels_colors(singleDend)), col = some_col_func)
dev.copy2pdf(file="MASH_heatmap_noclust.pdf", useDingbats=FALSE)
# Colors for final legend
rainbow_hcl(5)
	[1] "#E495A5" "#BDAB66" "#65BC8C" "#55B8D0" "#C29DDE"

# Rumen = pink = #E495A5
# Foot = tan = #BDAB66
# Manure = light green = #65BC8C
# Nasal = light blue = #55B8D0
# Oral = light purple = #C29DDE

# One last test, let's make a plot colored by sampling method
sampLookup <- hashmap(c("Biopsy", "Digester", "Fecal", "Fiber", "Oral", "Stomach Tube", "Swab"), c(1,2,3,4,5,6,7))
ordered_samp_cats <- rep(metadata.categories.final$Sampling.Method, 2)
get_ordered_7_clusters <- function(dend){sampLookup[[ordered_samp_cats[order.dendrogram(dend)]]]}
labels_colors(singleDend) <- rainbow_hcl(7)[get_ordered_7_clusters(singleDend)]
singleDend <- color_branches(singleDend, k=7)

gplots::heatmap.2(as.matrix(matrix_data.format), main = "Metagenome Tissue MASH Clustering", dendrogram = "row", Rowv = singleDend, Colv = "NA", trace = "none", denscol = "grey", density.info = "density", RowSideColors = rev(labels_colors(singleDend)), col = some_col_func)
gplots::heatmap.2(as.matrix(matrix_data.format), main = "Metagenome Tissue MASH Clustering", dendrogram = "none", Rowv = "NA", Colv = "NA", trace = "none", denscol = "grey", density.info = "density", RowSideColors = rev(labels_colors(singleDend)), col = some_col_func)

rainbow_hcl(7)
[1] "#E495A5" "#CEA472" "#9CB469" "#56BD96" "#46BAC8" "#99A9E2" "#D497D3"
# Biopsy	= pink = #E495A5
# Digester = tan = #CEA472
# Fecal  = light green = #9CB469
# Fiber = light blue-green = #56BD96
# Oral = light cyan = #46BAC8
# Stomach Tube = light purple = #99A9E2
# Swab
```

Finally, I'm going to attempt an Nonmetric MDS to see how the data fall out.

```R
library(MASS)
fit <- isoMDS(as.dist(matrix_data.format), k=2)
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", main = "MASH distance Nonmetric MDS", pch = c(15,16,17,18,25)[catLookup[[ordered_cats]]], col = c("blue", "dark green", "yellow", "orange", "red")[catLookup[[ordered_cats]]])
legend("topleft", c("Rumen", "Foot", "Manure", "Nasal", "Oral"), pch = c(15,16,17,18,25), col = c("blue", "dark green", "yellow", "orange", "red"), cex = 1.25)

# It looked good, but I had a stress value above 19. Anything above 0.3 is considered suspect!
library(vegan)
ord <- metaMDS(as.dist(matrix_data.format), trymax=100)
# It took 94 cycles to reach convergence

plot(ord, type="n")
orditorp(ord, display="sites", col=c("blue", "dark green", "purple", "orange", "red")[catLookup[[ordered_cats]]], pch = c(15,16,17,18,25)[catLookup[[ordered_cats]]])
legend("bottomleft", c("Rumen", "Foot", "Manure", "Nasal", "Oral"), pch = c(15,16,17,18,25), col = c("blue", "dark green", "purple", "orange", "red"), cex = 1.25, text.col = c("blue", "dark green", "purple", "orange", "red"))

# The 
```

I am also generating distance profiles with SourMash. Should be comparable and generates easy comparison heatmaps.

```bash
# I was only able to install sourmash thus far.
```

<a name="metaspades"></a>
### MetaSpades assembly

First, a test run of spades on the cluster with the meta tag runtime option:

> fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms

```bash
sbatch --nodes=1 --mem=200000 --ntasks-per-node=16 -p assemble2 ../../binaries/SPAdes-3.10.1-Linux/bin/metaspades.py -o DRR017219 -pe1-1 ../datasources/DRR017219_1.fastq -pe1-2 ../datasources/DRR017219_2.fastq
	ImportError: No module named spades_init
	# I need to add a Python path environmental variable here

sbatch --nodes=1 --mem=200000 --ntasks-per-node=16 -p assemble2 --wrap="PYTHONPATH=$PYTHONPATH:/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms; ../../binaries/SPAdes-3.10.1-Linux/bin/metaspades.py -o DRR017219 -1 ../datasources/DRR017219_1.fastq -2 ../datasources/DRR017219_2.fastq "

# Just testing a few things while it runs:
grep '>' /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms/DRR017219//K21/final_contigs.fasta | wc -l
12174565  <- that's 12 million contigs 
sbatch --nodes=1 --ntasks-per-node=1 --mem=1000 --wrap="module load samtools; samtools faidx /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms/DRR017219//K21/final_contigs.fasta"
cat /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms/DRR017219//K21/final_contigs.fasta.fai | cut -f2 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
	total   6357025
	Minimum 22
	Maximum 12160
	Average 60.678315
	Median  24
	Standard Deviation      98.195571
	Mode(Highest Distributed Value) 22

# OK, so this is to be expected. The assembly only takes me halfway the binning is what's needed next.
# Spades finished after 3 days of running. Here's the final tally:
grep '>' /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms/DRR017219/contigs.fasta | wc -l
	1,677,688

grep '>' /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms/DRR017219/scaffolds.fasta | wc -l
	1,677,530	<- only a few fewer than the contigs! Not much scaffolding information

# Lets gear up for the whole project
head -n 47 ../datasources/accession_list.txt > assembler2.list
tail -n 46 ../datasources/accession_list.txt > assembler3.list

# Assemble2
for i in `cat assembler2.list`; do echo $i; sbatch --nodes=1 --mem=200000 --ntasks-per-node=30 -p assemble2 -o ${i}.out -e ${i}.err --wrap="PYTHONPATH=$PYTHONPATH:/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms; ../../binaries/SPAdes-3.10.1-Linux/bin/metaspades.py -o $i -t 30 -m 200 -1 ../datasources/${i}_1.fastq -2 ../datasources/${i}_2.fastq "; done

# Assemble3
for i in `cat assembler3.list`; do echo $i; sbatch --nodes=1 --mem=200000 --ntasks-per-node=30 -p assemble3 -o ${i}.out -e ${i}.err --wrap="PYTHONPATH=$PYTHONPATH:/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/spades_asms; ../../binaries/SPAdes-3.10.1-Linux/bin/metaspades.py -o $i -t 30 -m 200 -1 ../datasources/${i}_1.fastq -2 ../datasources/${i}_2.fastq "; done
```

OK, several of these assemblies finished, but Serge advised me against my current path: consensus overlap of separate assemblies. I am now going to use MegaHit to assemble all of the sequence together in one go.

<a name="megahit"></a>
## MegaHit Assembly

Going to queue up the MegaHit run to go over the weekend.

> fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/metahit

```bash
sbatch --nodes=1 --ntasks-per-node=60 --mem=510000 -p assemble2 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/megahit_v1.1.1_LINUX_CPUONLY_x86_64-bin/megahit --k-min 27 --kmin-1pass -t 60 -1 ../datasources/DRR017219_1.fastq,../datasources/DRR019503_1.fastq,../datasources/DRR019504_1.fastq,../datasources/DRR019505_1.fastq,../datasources/DRR019506_1.fastq,../datasources/DRR019507_1.fastq,../datasources/DRR019508_1.fastq,../datasources/DRR019509_1.fastq,../datasources/DRR019510_1.fastq,../datasources/DRR019511_1.fastq,../datasources/ERR571345_1.fastq,../datasources/ERR571346_1.fastq,../datasources/ERR571347_1.fastq,../datasources/ERR571348_1.fastq,../datasources/ERR571349_1.fastq,../datasources/ERR571350_1.fastq,../datasources/ERR571351_1.fastq,../datasources/ERR571352_1.fastq,../datasources/ERR833046_1.fastq,../datasources/ERR833047_1.fastq,../datasources/ERR833048_1.fastq,../datasources/ERR833049_1.fastq,../datasources/ERR833050_1.fastq,../datasources/ERR833051_1.fastq,../datasources/ERR833052_1.fastq,../datasources/ERR833211_1.fastq,../datasources/ERR833212_1.fastq,../datasources/ERR833213_1.fastq,../datasources/ERR833214_1.fastq,../datasources/ERR833215_1.fastq,../datasources/ERR833216_1.fastq,../datasources/ERR833217_1.fastq,../datasources/ERR983262_1.fastq,../datasources/ERR983263_1.fastq,../datasources/ERR983264_1.fastq,../datasources/ERR983265_1.fastq,../datasources/ERR983266_1.fastq,../datasources/ERR983267_1.fastq,../datasources/ERR983268_1.fastq,../datasources/ERR983269_1.fastq,../datasources/SRR094166_1.fastq,../datasources/SRR094403_1.fastq,../datasources/SRR094405_1.fastq,../datasources/SRR094415_1.fastq,../datasources/SRR094416_1.fastq,../datasources/SRR094417_1.fastq,../datasources/SRR094418_1.fastq,../datasources/SRR094419_1.fastq,../datasources/SRR094424_1.fastq,../datasources/SRR094427_1.fastq,../datasources/SRR094428_1.fastq,../datasources/SRR094429_1.fastq,../datasources/SRR094437_1.fastq,../datasources/SRR094926_1.fastq,../datasources/SRR1725612_1.fastq,../datasources/SRR1776513_1.fastq,../datasources/SRR1805214_1.fastq,../datasources/SRR1805215_1.fastq,../datasources/SRR1810823_1.fastq,../datasources/SRR1810824_1.fastq,../datasources/SRR2016834_1.fastq,../datasources/SRR2016841_1.fastq,../datasources/SRR2016847_1.fastq,../datasources/SRR2016852_1.fastq,../datasources/SRR2145287_1.fastq,../datasources/SRR2329878_1.fastq,../datasources/SRR2329910_1.fastq,../datasources/SRR2329939_1.fastq,../datasources/SRR2329962_1.fastq,../datasources/SRR2329984_1.fastq,../datasources/SRR3166066_1.fastq,../datasources/SRR3166067_1.fastq,../datasources/SRR3166092_1.fastq,../datasources/SRR3166102_1.fastq,../datasources/SRR4435556_1.fastq,../datasources/SRR4435557_1.fastq,../datasources/SRR4435558_1.fastq,../datasources/SRR4435559_1.fastq,../datasources/SRR4435560_1.fastq,../datasources/SRR4435561_1.fastq,../datasources/SRR4435562_1.fastq,../datasources/SRR4435563_1.fastq,../datasources/SRR4435564_1.fastq,../datasources/SRR4435565_1.fastq,../datasources/SRR4435566_1.fastq,../datasources/SRR4435567_1.fastq,../datasources/SRR4435578_1.fastq,../datasources/SRR4435579_1.fastq,../datassqueueources/SRR4435580_1.fastq,../datasources/SRR4435581_1.fastq,../datasources/SRR493932_1.fastq,../datasources/SRR493933_1.fastq,../datasources/SRR493934_1.fastq -2 ../datasources/DRR017219_2.fastq,../datasources/DRR019503_2.fastq,../datasources/DRR019504_2.fastq,../datasources/DRR019505_2.fastq,../datasources/DRR019506_2.fastq,../datasources/DRR019507_2.fastq,../datasources/DRR019508_2.fastq,../datasources/DRR019509_2.fastq,../datasources/DRR019510_2.fastq,../datasources/DRR019511_2.fastq,../datasources/ERR571345_2.fastq,../datasources/ERR571346_2.fastq,../datasources/ERR571347_2.fastq,../datasources/ERR571348_2.fastq,../datasources/ERR571349_2.fastq,../datasources/ERR571350_2.fastq,../datasources/ERR571351_2.fastq,../datasources/ERR571352_2.fastq,../datasources/ERR833046_2.fastq,../datasources/ERR833047_2.fastq,../datasources/ERR833048_2.fastq,../datasources/ERR833049_2.fastq,../datasources/ERR833050_2.fastq,../datasources/ERR833051_2.fastq,../datasources/ERR833052_2.fastq,../datasources/ERR833211_2.fastq,../datasources/ERR833212_2.fastq,../datasources/ERR833213_2.fastq,../datasources/ERR833214_2.fastq,../datasources/ERR833215_2.fastq,../datasources/ERR833216_2.fastq,../datasources/ERR833217_2.fastq,../datasources/ERR983262_2.fastq,../datasources/ERR983263_2.fastq,../datasources/ERR983264_2.fastq,../datasources/ERR983265_2.fastq,../datasources/ERR983266_2.fastq,../datasources/ERR983267_2.fastq,../datasources/ERR983268_2.fastq,../datasources/ERR983269_2.fastq,../datasources/SRR094166_2.fastq,../datasources/SRR094403_2.fastq,../datasources/SRR094405_2.fastq,../datasources/SRR094415_2.fastq,../datasources/SRR094416_2.fastq,../datasources/SRR094417_2.fastq,../datasources/SRR094418_2.fastq,../datasources/SRR094419_2.fastq,../datasources/SRR094424_2.fastq,../datasources/SRR094427_2.fastq,../datasources/SRR094428_2.fastq,../datasources/SRR094437_2.fastq,../datasources/SRR094926_2.fastq,../datasources/SRR1725612_2.fastq,../datasources/SRR1776513_2.fastq,../datasources/SRR1805214_2.fastq,../datasources/SRR1805215_2.fastq,../datasources/SRR1810823_2.fastq,../datasources/SRR1810824_2.fastq,../datasources/SRR2016834_2.fastq,../datasources/SRR2016841_2.fastq,../datasources/SRR2016847_2.fastq,../datasources/SRR2016852_2.fastq,../datasources/SRR2145287_2.fastq,../datasources/SRR2329878_2.fastq,../datasources/SRR2329910_2.fastq,../datasources/SRR2329939_2.fastq,../datasources/SRR2329962_2.fastq,../datasources/SRR2329984_2.fastq,../datasources/SRR3166066_2.fastq,../datasources/SRR3166067_2.fastq,../datasources/SRR3166092_2.fastq,../datasources/SRR3166102_2.fastq,../datasources/SRR4435556_2.fastq,../datasources/SRR4435557_2.fastq,../datasources/SRR4435558_2.fastq,../datasources/SRR4435559_2.fastq,../datasources/SRR4435560_2.fastq,../datasources/SRR4435561_2.fastq,../datasources/SRR4435562_2.fastq,../datasources/SRR4435563_2.fastq,../datasources/SRR4435564_2.fastq,../datasources/SRR4435565_2.fastq,../datasources/SRR4435566_2.fastq,../datasources/SRR4435567_2.fastq,../datasources/SRR4435578_2.fastq,../datasources/SRR4435579_2.fastq,../datasources/SRR4435580_2.fastq,../datasources/SRR4435581_2.fastq,../datasources/SRR493932_2.fastq,../datasources/SRR493933_2.fastq,../datasources/SRR493934_2.fastq -o rumen_megahit_asm"

  [ERROR] [sequence_manager.cpp     : 123]: File(s) ../datasources/SRR1776513_1.fastq,../datasources/SRR1776513_2.fastq: Number of sequences not the same in paired files. Abort.
Error occurs when running "megahit_asm_core buildlib"; please refer to rumen_megahit_asm/log for detail

wc -l ../datasources/SRR1776513_1.fastq ../datasources/SRR1776513_2.fastq
  3020000 ../datasources/SRR1776513_1.fastq
  3019996 ../datasources/SRR1776513_2.fastq
```

Damn! Gotta remove that dataset because of the number of reads. 

An important consideration: I did not remove adaptor/faulty bases! I need to run sickle on this data first and then selectively add files to the assembly.


> fry: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/datasources

```bash
mkdir sickle

for i in `cat accession_list.txt`; do echo $i; sbatch --nodes=1 --ntasks-per-node=2 --mem=5000 -o sickle/${i}_slurm.out --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/sickle/sickle pe -f ${i}_1.fastq -r ${i}_2.fastq -o sickle/${i}_sickle.1.fastq.gz -p sickle/${i}_sickle.2.fastq.gz -s sickle/${i}_singles.fastq.gz -t sanger -q 15 -l 36 -g"; done
```
