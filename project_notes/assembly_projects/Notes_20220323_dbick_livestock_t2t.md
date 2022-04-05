# T2T livestock genome assemblies
---
*3/23/2022*

These are my notes on the analysis and validation of the T2T livestock genome assemblies.

## Table of Contents

## rDNA array analysis

I want to write a script/program that identifies rDNA arrays automatically from either repeatmasker or some other program's output.

> Ceres: /project/rumen_longread_metagenome_assembly/analysis/rdna

```bash
module load python_3/3.6.6 trnascanse/2.0.5

sbatch -N 1 -n 30 --mem=150000 -p priority -q msn --wrap='tRNAscan-SE -o verkko_sheep.trnascan.out --thread 30 -m verkko_sheep.trnascan.stats -b verkko_sheep.trnascan.bed /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta'
```

Nope! That just produced tRNA assignments contrary to what the BioStars people suggested. Instead, let's work with Ben's RM file and plot the rRNA gene locations first to make sure that they are in the right places.

```bash
module load python_3/3.6.6 miniconda/3.6

conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/

# first to extract the rRNA locations to a bed file
perl -e '<>; <>;<>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[10] eq "rRNA"){print("$s[4]\t$s[5]\t$s[6]\t$s[9]\n");}}' < assembly.trio_resolved.fasta.out > sheep_rrna_locs.bed

python3 ~/python_toolchain/sequenceData/plotBedIdeogram.py -f renamed_gapped.fasta.fai -b sheep_rrna_locs.bed -o sheep_rrna_locs.pdf -t 50

# the tool works, but the plot looks really noisy. Let's try a filter on the repeatmasker output

```

OK, having read the T2T manuscript, I think that I know how to generate more accurate representations of the rDNA clusters using read recruitment and sparse graph construction. Here is the order of operations:

#### The strategy

1. Filter RDNA annotations to identify likely consensus regions for read recruitment.
2. Select HiFi reads that align to the rDNA regions and generate a sparse graph with MBG
3. Resolve the sparse graph with UL ONT reads that are graph aligned
4. Pull reads that map to the combined graph that appear to be in chromosome-specific morphs and separate into pools
5. Assemble read pools to reconstruct the morphs for each array

#### Data sources

* Hifi reads: /90daydata/sheep_genome_assemblies/sergek/hifi/
* ONT reads: /90daydata/sheep_genome_assemblies/sergek/nano/PAH86096/
* Assembly: /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta

> Ceres: /90daydata/rumen_longread_metagenome_assembly/analysis/rdna

```bash
module load python_3/3.6.6 miniconda/3.6 r/4.1.2
# Some stats and distributions
# divergence
perl -e '<>; <>;<>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[10] eq "rRNA"){print "$s[1]\n";}}' < /project/rumen_longread_metagenome_assembly/analysis/rdna/assembly.trio_resolved.fasta.out | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3889
Sum:    39205.9
Minimum 0.0
Maximum 38.2
Average 10.081229
Median  10.0
Standard Deviation      8.435337
Mode(Highest Distributed Value) 0.8

# Length
perl -e '<>; <>;<>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[10] eq "rRNA"){$l = $s[6] - $s[5]; print "$l\n";}}' < /project/rumen_longread_metagenome_assembly/analysis/rdna/assembly.trio_resolved.fasta.out | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3889
Sum:    2745168
Minimum 28
Maximum 4928
Average 705.880175
Median  120
Standard Deviation      1338.690503
Mode(Highest Distributed Value) 120

# Generating data for distribution plotting
perl -e 'print "Div\tLen\tChrom\tType\n"; <>; <>;<>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[10] eq "rRNA"){$l = $s[6] - $s[5]; print "$s[1]\t$l\t$s[4]\t$s[9]\n";}}' < /project/rumen_longread_metagenome_assembly/analysis/rdna/assembly.trio_resolved.fasta.out > sheep_rrna_data.tab

```

```R
library(ggplot2)
library(dplyr)

data <- read.delim("sheep_rrna_data.tab", sep="\t", head=TRUE)

pdf(file="sheep_rrna_data.bytype.pdf", useDingbats=FALSE)
ggplot(data=data, aes(x=Div, y=Len, color=Type)) + geom_point() + theme_bw() + scale_color_brewer(palette="Dark2") + facet_grid(Type ~ ., scales="free_y")
dev.off()

summtable <- data %>% group_by(Chrom,Type) %>% summarize(count = n(), lenAvg = mean(Len), lenMax = max(Len), divAvg = mean(Div), divMin = min(Div))
write.table(x=summtable, file="sheep_rrna_group_summary.tab", sep="\t", quote=FALSE)
```

Running through distributions of length, divergence and count, I think that I have narrowed down unitigs that are suitable for read recruitment. Let's start with some of the unused unitigs.

#### LSU

|Chrom	|Type	|count	|lenAvg	|lenMax	| divAvg	| divMin |
|:--- | :--- | :--- | :--- | :--- | :--- | :--- |
|unused_utig4-1265_len_158503	| LSU-rRNA_Hsa	| 5	|4906.2	|4908	|5.14	|5.1|
|unused_utig4-1757_len_59982	|LSU-rRNA_Hsa	|3	|4908.333333	|4910	|5.2	|5.2|
|unused_utig4-1947_len_47983	| LSU-rRNA_Hsa	| 3	| 4910.333333	| 4916	| 5.2	| 5.2|
|unused_utig4-729_len_81941	|LSU-rRNA_Hsa	|3	|4909	|4911	|5.266666667	|5.2|
|unused_utig4-731_len_119185	|LSU-rRNA_Hsa	|4	|4908.75	|4911	|5.275	|5.2|
|unused_utig4-732_len_89511	|LSU-rRNA_Hsa	|3	|4909	|4911	|5.266666667	|5.2|

#### SSU

| Chrom                        | Type         | count | lenAvg      | lenMax | divAvg      | divMin |
|------------------------------|--------------|-------|-------------|--------|-------------|--------|
| unused_utig4-1265_len_158503 | SSU-rRNA_Hsa | 5     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-1267_len_74267  | SSU-rRNA_Hsa | 3     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-1446_len_72491  | SSU-rRNA_Hsa | 3     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-1757_len_59982  | SSU-rRNA_Hsa | 3     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-1947_len_47983  | SSU-rRNA_Hsa | 3     | 1800        | 1868   | 0.566666667 | 0.5    |
| unused_utig4-719_len_71352   | SSU-rRNA_Hsa | 4     | 1866.75     | 1868   | 0.5         | 0.5    |
| unused_utig4-721_len_61239   | SSU-rRNA_Hsa | 3     | 1867.666667 | 1868   | 0.5         | 0.5    |
| unused_utig4-729_len_81941   | SSU-rRNA_Hsa | 3     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-731_len_119185  | SSU-rRNA_Hsa | 5     | 1868        | 1868   | 0.5         | 0.5    |
| unused_utig4-732_len_89511   | SSU-rRNA_Hsa | 3     | 1868        | 1868   | 0.5         | 0.5    |


#### 5S

| Chrom                       | Type | count | lenAvg      | lenMax | divAvg      | divMin |
|-----------------------------|------|-------|-------------|--------|-------------|--------|
| path_from_utig4-1326        | 5S   | 66    | 119.1515152 | 120    | 1.215151515 | 0.8    |
| path_from_utig4-1936        | 5S   | 286   | 117.9090909 | 121    | 1.76048951  | 0.8    |
| path_from_utig4-2188        | 5S   | 191   | 117.2722513 | 121    | 2.042408377 | 0.8    |
| path_from_utig4-815         | 5S   | 62    | 119.1451613 | 120    | 1.461290323 | 0.8    |
| unused_utig4-1013_len_84705 | 5S   | 52    | 120.25      | 121    | 1           | 0.8    |
| unused_utig4-668_len_154406 | 5S   | 112   | 119.9910714 | 122    | 0.88125     | 0.8    |

So, I think that the units should be all of the LSU unused unitigs (perhaps not the first with length 158,503, but we'll see), and subsections of the 5S path unitigs if they are chains of 5S sequences. The following path coordinates appear to be legitimate in terms of arrays of 5S spaced in arrays.

path_from_utig4-1936    42714874	43277189
path_from_utig4-2188	1099	606302
path_from_utig4-1326	9180112	9306592
path_from_utig4-815	971	119344

#### Running the alignment to generate the read list

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2 r/4.1.2

samtools faidx /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta unused_utig4-1265_len_158503 unused_utig4-1757_len_59982 unused_utig4-1947_len_47983 unused_utig4-729_len_81941 unused_utig4-731_len_119185 unused_utig4-732_len_89511 unused_utig4-1013_len_84705 unused_utig4-668_len_154406 path_from_utig4-1936:42714874-43277189 path_from_utig4-2188:1099-606302 path_from_utig4-1326:9180112-9306592 path_from_utig4-815:971-119344 > sheep_verkko_rdna_clusters.fasta
samtools faidx sheep_verkko_rdna_clusters.fasta

for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do echo -n "$i "; done; echo
sbatch -N 1 -n 30 -p priority -q msn --mem=180000 --wrap='minimap2 -x map-hifi -t 30 -o sheep_verkko_rdna_alignments.paf /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210601_165737.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210610_143100.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210611_233100.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210621_195743.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210623_154734.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210716_180058.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210718_025906.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210726_175325.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_210728_025148.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211026_163028.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211028_032534.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211104_191436.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211106_060943.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211118_170215.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m54337U_211120_035803.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m64015e_211119_011043.hifi_reads.fastq.gz /90daydata/sheep_genome_assemblies/sergek/hifi/m64015e_211120_120456.hifi_reads.fastq.gz'

# That was taking forever. Let's queue each read job separately so we can scatter-gather later
mkdir hifi_rdna_aligns
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="minimap2 -x map-hifi -t 5 -o hifi_rdna_aligns/$name.alignments.paf /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/gapped/renamed_gapped.fasta $i"; done;

# I want to see the distribution of read alignments and make sure that I need to filter the list in the first place!
perl -e 'print "LenRatio\tMapRatio\n"; %data; while(<>){chomp; @s = split(/\t/); $lr = $s[10] / $s[1]; $mr = $s[9] / $s[10]; if(exists($data{$s[0]})){if($data{$s[0]}->[0] < $lr){$data{$s[0]}->[0] = $lr; $data{$s[0]}->[1] = $mr;}}else{$data{$s[0]} = [$lr, $mr];}} foreach my $k (keys(%data)){print "$data{$k}->[0]\t$data{$k}->[1]\n";}' < hifi_rdna_aligns/m54337U_210726_175325.alignments.paf > test_alignment_qualities.tab
(base) [derek.bickharhth@ceres-login rdna]$ head test_alignment_qualities.tab
```

```R
library(dplyr)
library(ggplot2)
library(gridExtra)
data <- read.delim("test_alignment_qualities.tab", sep="\t", header=TRUE)

summary(data)
    LenRatio           MapRatio
 Min.   :0.003551   Min.   :0.005176
 1st Qu.:0.992278   1st Qu.:0.766556
 Median :0.997706   Median :0.857084
 Mean   :0.989266   Mean   :0.806359
 3rd Qu.:0.999537   3rd Qu.:0.928715
 Max.   :2.224475   Max.   :1.000000

data <- mutate(data, sumRatio = LenRatio + MapRatio)
p1 <- ggplot(data, aes(x=LenRatio, y=MapRatio)) + geom_point() + scale_fill_brewer(palette="Dark2") + theme_bw()
p2 <- ggplot(data, aes(x="sumRatio", y=sumRatio)) + geom_violin(trim=FALSE, fill="red") + geom_boxplot(width=0.1) + theme_bw()

pdf(file="test_alignment_qualities.pdf", useDingbats=FALSE)
grid.arrange(p1, p2, ncol=1, nrow=2)
dev.off()

pdf(file="test_alignment_qualities_density.pdf", useDingbats=FALSE)
ggplot(data, aes(x=LenRatio, y=MapRatio)) + geom_bin2d() + scale_fill_continuous(type = "viridis") + theme_bw()
dev.off()
```

It looks like the sweet spot that contains most of the alignments is between a Length ratio of 90-100% and a map ratio of 80-100%. I will filter the reads and then I will try to set up three different MBG runs with different kmer values. First will be 2001, then 3501 (the paper), and then 5001. Window sizes will be (k-1)/10 and max kmer sizes will be 15000 uniformly.

```bash
module load python_3/3.6.6 miniconda/3.6 samtools minimap2 r/4.1.2

python3 filterHifiRDNAPafs.py hifi_rdna_aligns/m54337U_210726_175325.alignments.paf test_reads.list
wc -l test_reads.list
	72771 test_reads.list
wc -l test_alignment_qualities.tab
	99977 test_alignment_qualities.tab

# SO about a 30% reduction based on those criteria!
# Let's queue it up
mkdir hifi_rdna_lists
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="python3 filterHifiRDNAPafs.py hifi_rdna_aligns/$name.alignments.paf hifi_rdna_lists/$name.reads.list"; done

conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/
mkdir hifi_rdna_fastqs
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=25000 --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q $i -l hifi_rdna_lists/$name.reads.list -o hifi_rdna_fastqs/$name.rdna.fastq; gzip hifi_rdna_fastqs/$name.rdna.fastq"; done


#### Queuing up MBG runs
conda activate /project/rumen_longread_metagenome_assembly/environments/verkko
for i in hifi_rdna_fastqs/*; do echo -n "-i $i "; done; echo

mkdir hifi_rdna_mbg
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_2001_200.slurm.out --wrap="MBG -t 35 -k 2001 -w 200 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_2001_200_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_2001_200_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"

sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_3501_350.slurm.out --wrap="MBG -t 35 -k 3501 -w 350 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_3501_350_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_3501_350_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"

sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_5001_500.slurm.out --wrap="MBG -t 35 -k 5001 -w 500 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_5001_500_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_5001_500_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"

# Damn, they're too lengthy. The repeat structure is too big. Let's just pull the reads that map exclusively to rDNA sequences.
mkdir hifi_paragon_fastqs

```