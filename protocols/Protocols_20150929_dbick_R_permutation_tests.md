# Running Permutation tests with R
---
*9/29/2015*

These are my notes on how to run a permutation test of interval intersections with R package ["regioneR"](http://www.bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.pdf).

Lingyang's manuscript needs a permutation test to see if identified CNVRs overlap with genes more frequently than would be predicted at random. Let's work with this data in order to find out how to run the test.

> pwd: /home/dbickhart/share/side_projects/lingyang_projects

```bash
# Preparing the files for loading into R
# I need bed files for the GRanges object creation process
cp ../../umd3_data/genedb/ucsc_umd3_otherspec_refgene.bed ./
cp ../../umd3_data/genedb/ensGene_umd3_coords.bed ./
cp ../../umd3_data/genedb/refGeneName_umd3_coords.bed ./

# Let's prepare an "all gene" bed file for the "gene universe" option of regioneR as well
cat ensGene_umd3_coords.bed refGeneName_umd3_coords.bed ucsc_umd3_otherspec_refgene.bed | sortBedFileSTDIN.pl > gene_universe.bed

# Finally, I need the fasta index file in order to generate the "genome" object
cp ../../umd3_data/umd3_kary_extend_hgap.fa.fai ./
```

OK, it's time to start tinkering around in R

```R
# Loading the CNVs
cnvs.dels <- toGRanges("deletion_cnvs.bed")
cnvs.all <- toGRanges("all_cnvs.bed")

# Loading gene lists
genes.all <- toGRanges("gene_universe.bed")
genes.ensgene <- toGRanges("ensGene_umd3_coords.bed")
genes.xeno <- toGRanges("ucsc_umd3_otherspec_refgene.bed")
genes.refgene <- toGRanges("refGeneName_umd3_coords.bed")

# Making the genome length object
umd3.faidx <- read.table("umd3_kary_extend_hgap.fa.fai")
umd3.genome <- umd3.faidx[,c("V1", "V2")]
colnames(umd3.genome) <- c("chr", "length")

# OK, since Lingyang used the HD chip genotypes, and did not filter away any regions in particular, I'm not
# going to use a mask this time.

# First, refgene association
pt.all.refgene <- overlapPermTest(A=cnvs.all, B=genes.refgene, ntimes=1000, genome=umd3.genome, non.overlapping=FALSE)
pt.del.refgene <- overlapPermTest(A=cnvs.dels, B=genes.refgene, ntimes=1000, genome=umd3.genome, non.overlapping=FALSE)

# Data from the tests
pt.all.refgene
	P-value: 0.000999000999000999
	Z-score: -3.145
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 67
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
pt.del.refgene
	P-value: 0.000999000999000999
	Z-score: -2.36
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 49
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions

# Testing if "non-overlapping" changes the statistics significantly
pt.all.refgene <- overlapPermTest(A=cnvs.all, B=genes.refgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.all.refgene
	P-value: 0.000999000999000999
	Z-score: -3.1296
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 67
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
# Slight reduction of Z-score

# OK, so at least the refGene works! Let's try ensgene Ids now
pt.all.ensgene <- overlapPermTest(A=cnvs.all, B=genes.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.del.ensgene <- overlapPermTest(A=cnvs.dels, B=genes.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)

pt.all.ensgene
	P-value: 0.215784215784216
	Z-score: 0.7393
	Number of iterations: 1000
	Alternative: greater
	Evaluation of the original region set: 230
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
pt.del.ensgene
	P-value: 0.210789210789211
	Z-score: 0.6446
	Number of iterations: 1000
	Alternative: greater
	Evaluation of the original region set: 158
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
# Not significant
# Let's double-check my ensgene listing to see how many gene entries there are.
```

I'm going to see if I need to make sure that my gene bed files are non-overlapping!

```bash
../../programs_source/Perl/backup/table_bed_length_sum.pl ensGene_umd3_coords.bed
FName   IntNum  TotLen  LenAvg  LenStdev        LenMedian       SmallestL       LargestL
ensGene_umd3_coords.bed 26740   873,599,724       32670.1467464473        71526.8187146969        9293    26      1851199

# That's 873 megabases! Let's see what merging the bed does to the length
cat ensGene_umd3_coords.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -nms > ensGene_umd3_merged.bed

../../programs_source/Perl/backup/table_bed_length_sum.pl ensGene_umd3_coords.bed ensGene_umd3_merged.bed
FName   IntNum  TotLen  LenAvg  LenStdev        LenMedian       SmallestL       LargestL
ensGene_umd3_coords.bed 26740   873,599,724       32670.1467464473        71526.8187146969        9293    26      1851199
ensGene_umd3_merged.bed 22313   776,950,769       34820.542688119 73954.7187827353        10463   26      1851199

# A reduction of 100 Mb -- pretty significant!
# Let's merge the other files and attempt it again with R
cat ucsc_umd3_otherspec_refgene.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -nms > ucsc_umd3_otherspec_merged.bed
cat gene_universe.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -nms > gene_universe_merged.bed
```

OK, going to try the permutations again with R.

```R
genes.all <- toGRanges("gene_universe_merged.bed")
genes.ensgene <- toGRanges("ensGene_umd3_merged.bed")
genes.xeno <- toGRanges("ucsc_umd3_otherspec_merged.bed")

# Now the ensgene test again
pt.all.ensgene <- overlapPermTest(A=cnvs.all, B=genes.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.all.ensgene
	P-value: 0.24975024975025
	Z-score: 0.5444
	Number of iterations: 1000
	Alternative: greater
	Evaluation of the original region set: 192
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions

# Even less than before! Let's try xenorefgene instead
pt.all.xeno <- overlapPermTest(A=cnvs.all, B=genes.xeno, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.all.xeno
	P-value: 0.001998001998002
	Z-score: -2.7148
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 155
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions

# That was significant, but let's see the all vs all comparison
pt.all.all <- overlapPermTest(A=cnvs.all, B=genes.all, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.all.all
	P-value: 0.363636363636364
	Z-score: -0.3521
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 200
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions

# Not significant
# I'm wondering if I downloaded a weird copy of ensemble genes now...
```

Taking a look at my files...

```bash
perl ../../programs_source/Perl/backup/table_bed_length_sum.pl ensGene_umd3_merged.bed gene_universe_merged.bed ucsc_umd3_otherspec_merged.bed refGeneName_umd3_coords.bed
FName   IntNum  TotLen  LenAvg  LenStdev        LenMedian       SmallestL       LargestL
ensGene_umd3_merged.bed 22313   776,950,769       34820.542688119 73954.7187827353        10463   26      1851199
gene_universe_merged.bed        15948   1,363,532,605      85498.6584524705        206189.500051045        16663   265987469
ucsc_umd3_otherspec_merged.bed  14392   1,345,143,423      93464.6625208449        212594.458884406        22038   305987469
refGeneName_umd3_coords.bed     14176   620,495,211       43770.8247037246        92854.9904886637        15685   342345944

# So weird! There's apparently 20 megabases between the xeno and merged dataset that are 
# contained within the ensgene file that account for most of the CNV regions!

# downloaded the 1.81 Umd3.1 gff3 file from Ensembl
perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f Bos_taurus.UMD3.1.81.gff3 -c 2 -e '#' -m
```
#### Ensembl 1.81 annotation table

|Entry                |  Count|
|:--------------------|------:|
|CDS                  | 214704|
|RNA                  |    175|
|chromosome           |     31|
|exon                 | 227731|
|five_prime_UTR       |  17713|
|gene                 |  20167|
|miRNA                |   1153|
|miRNA_gene           |   1153|
|mt_gene              |     22|
|processed_pseudogene |    171|
|pseudogene           |   1252|
|rRNA                 |    405|
|rRNA_gene            |    405|
|snRNA                |   1222|
|snRNA_gene           |   1222|
|snoRNA               |    846|
|snoRNA_gene          |    846|
|supercontig          |   3286|
|three_prime_UTR      |  15442|
|transcript           |  22317|

OK, let's pull out only the genes.

```bash
perl -lane 'if($F[0] =~ /^#/){next;}elsif($F[2] eq "gene"){$F[8] =~ /ID=gene:(.{16,19});/; print "chr$F[0]\t$F[3]\t$F[4]\t$1";}' < Bos_taurus.UMD3.1.81.gff3 > ensgene.1.81.named.bed

wc -l ensgene.1.81.named.bed
	20167 ensgene.1.81.named.bed
cat ensgene.1.81.named.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -nms | wc -l
	19340

# Let's merge it to make the permutation testing fair
cat ensgene.1.81.named.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -nms > ensgene.1.81.named.merged.bed

perl ../../programs_source/Perl/backup/table_bed_length_sum.pl ensgene.1.81.named.merged.bed ensGene_umd3_merged.bed
	FName   IntNum  TotLen  LenAvg  LenStdev        LenMedian       SmallestL       LargestL
	ensgene.1.81.named.merged.bed   19340   776,178,447       40133.321975181 78088.2242353799        14594   98      1851198
	ensGene_umd3_merged.bed 22313   776,950,769       34820.542688119 73954.7187827353        10463   26      1851199

# Similar in terms of length... Oh well, let's test it in R
```

Working in the same directory in R:

```R
genes.ensgene <- toGRanges("ensgene.1.81.named.merged.bed")

pt.all.ensgene <- overlapPermTest(A=cnvs.all, B=genes.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)

pt.all.ensgene
	P-value: 0.218781218781219
	Z-score: 0.7188
	Number of iterations: 1000
	Alternative: greater
	Evaluation of the original region set: 182
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
# Nope! Still not significant!

pt.dels.ensgene <- overlapPermTest(A=cnvs.dels, B=genes.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)

pt.dels.ensgene
	P-value: 0.230769230769231
	Z-score: 0.5866
	Number of iterations: 1000
	Alternative: greater
	Evaluation of the original region set: 125
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions
# Not significant either
```

Let's summarize the results so far:

| CNV set | Gene Set | p value | z score |
| :--- | :--- | ---: | ---: |
all | refgene | 0.0009 | -3.145
dels | refgene | 0.0009 | -2.36
all | ensgene | 0.218 | 0.7188
dels | ensgene | 0.230 | 0.5866
all | xeno | 0.0019 | -2.7148
dels | xeno | 0.0029 | -2.3743
all | all | 0.3636 | -0.3521


Testing one last thing! Let's see just how much those 20 megabases not found in the xeno regions influence the ensgene statistics.

```bash
intersectBed -b ucsc_umd3_otherspec_merged.bed -a ensgene.1.81.named.merged.bed  > ensgene.1.81.trimmed.bed
```

```R
test.ensgene <- toGRanges("ensgene.1.81.trimmed.bed")
pt.test.all.ensgene <- overlapPermTest(A=cnvs.all, B=test.ensgene, ntimes=1000, genome=umd3.genome, non.overlapping=TRUE)
pt.test.all.ensgene
	P-value: 0.141858141858142
	Z-score: -1.0248
	Number of iterations: 1000
	Alternative: less
	Evaluation of the original region set: 140
	Evaluation function: numOverlaps
	Randomization function: randomizeRegions

# Hmm... about 40 CNVRs intersect with ensgenes within a 17 megabase region? What genes are there?
```

```bash
subtractBed -a ensgene.1.81.named.merged.bed -b ucsc_umd3_otherspec_merged.bed > ensegene.1.81.unique.bed

intersectBed -a ensegene.1.81.unique.bed -b all_cnvs.bed | less
# I found out the reason:
# There are several large, nearly overlapping CNVs on chromosome 10 that overlap dozens of ensembl gene IDs 
# that are otherwise not annotated.
```

The problematic region is [here](http://genome.ucsc.edu/cgi-bin/hgTracks?db=bosTau6&position=chr10%3A23180741-25173443&hgsid=442136149_dzdZINHHQr3Nv998N8qGGAraOeo4) and there are tons of gaps along with lots of predicted (but not fully annotated) genes. They all have hits to "T-cell receptor" genes, but most of the data for this comes from ESTs, not sequence homology.