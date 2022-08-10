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

|Chrom							|Type			|count	|lenAvg	|lenMax	| divAvg| divMin |
|:--- 							| :--- 			| :--- 	| :--- | :--- | :--- | :--- |
|unused_utig4-1265_len_158503	|LSU-rRNA_Hsa	|5		|4906.2	|4908	|5.14	|5.1|
|unused_utig4-1757_len_59982	|LSU-rRNA_Hsa	|3		|4908.3	|4910	|5.2	|5.2|
|unused_utig4-1947_len_47983	|LSU-rRNA_Hsa	|3		|4910.3	| 4916	|5.2	|5.2|
|unused_utig4-729_len_81941		|LSU-rRNA_Hsa	|3		|4909	|4911	|5.26	|5.2|
|unused_utig4-731_len_119185	|LSU-rRNA_Hsa	|4		|4908.7	|4911	|5.275	|5.2|
|unused_utig4-732_len_89511		|LSU-rRNA_Hsa	|3		|4909	|4911	|5.26	|5.2|
	
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

> Ceres: /90daydata/rumen_longread_metagenome_assembly/analysis/rdna

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
# The consensus 45S sequence is here: 45s_array_consensus.fasta
mkdir hifi_paragon_aligns

for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="minimap2 -x map-hifi -t 5 -o hifi_paragon_aligns/$name.alignments.paf 45s_array_consensus.fasta $i"; done

# Separate sam run for read depth estimation
mkdir hifi_sam_aligns
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="minimap2 -ax map-hifi -t 5 -o hifi_sam_aligns/$name.alignments.sam 45s_array_consensus.fasta $i"; done

for i in hifi_sam_aligns/*.sam; do echo $i; name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 -t 1-0 --wrap="samtools sort --reference 45s_array_consensus.fasta -o hifi_sam_aligns/$name.bam --threads 5 -T $name $i"; done 

perl -e 'print "LenRatio\tMapRatio\tMaxMQ\n"; %data; while(<>){chomp; @s = split(/\t/); $lr = $s[10] / $s[1]; $mr = $s[9] / $s[10]; if(exists($data{$s[0]})){if($data{$s[0]}->[0] < $lr){$data{$s[0]}->[0] = $lr; $data{$s[0]}->[1] = $mr; $data{$s[0]}->[2] = $s[11];}}else{$data{$s[0]} = [$lr, $mr, $s[11]];}} foreach my $k (keys(%data)){print "$data{$k}->[0]\t$data{$k}->[1]\t$data{$k}->[2]\n";}' < hifi_paragon_aligns/m54337U_211106_060943.alignments.paf > paragon_test_aligns.tab

mkdir hifi_paragon_lists
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="python3 filterHifiRDNAByQPafs.py hifi_paragon_aligns/$name.alignments.paf hifi_paragon_lists/$name.reads.list"; done

mkdir hifi_paragon_flists
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 5 -p priority -q msn --mem=22000 --wrap="python3 filterHifiRDNAByQPafsF.py hifi_paragon_aligns/$name.alignments.paf hifi_paragon_flists/$name.reads.list"; done

conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/
mkdir hifi_paragon_fastqs
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=25000 --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q $i -l hifi_paragon_lists/$name.reads.list -o hifi_paragon_fastqs/$name.rdna.fastq; gzip hifi_paragon_fastqs/$name.rdna.fastq"; done

conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/
mkdir hifi_paragon_ffastqs
for i in /90daydata/sheep_genome_assemblies/sergek/hifi/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=25000 --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q $i -l hifi_paragon_flists/$name.reads.list -o hifi_paragon_ffastqs/$name.rdna.fastq; gzip hifi_paragon_ffastqs/$name.rdna.fastq"; done

mkdir hifi_paragon_mbg
conda activate /project/rumen_longread_metagenome_assembly/environments/verkko
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_paragon_3501_350.slurm.out --wrap="MBG -t 35 -k 3501 -w 350 -r 15000 --output-sequence-paths hifi_paragon_mbg/sheep_rdna_3501_350_paths.gaf  --out hifi_paragon_mbg/sheep_rdna_3501_350_graph.gfa -i hifi_paragon_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_paragon_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_paragon_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_paragon_fastqs/m64015e_211120_120456.rdna.fastq.gz"

mkdir hifi_fparagon_mbg
conda activate /project/rumen_longread_metagenome_assembly/environments/verkko
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_2501_250.slurm.out --wrap="MBG -t 35 -k 2501 -w 250 -r 15000 --output-sequence-paths hifi_fparagon_mbg/sheep_rdna_2501_250_paths.gaf  --out hifi_fparagon_mbg/sheep_rdna_2501_250_graph.gfa -i hifi_paragon_ffastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211120_120456.rdna.fastq.gz"
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_5000_500.slurm.out --wrap="MBG -t 35 -k 5001 -w 500 -r 15000 --output-sequence-paths hifi_fparagon_mbg/sheep_rdna_5001_500_paths.gaf  --out hifi_fparagon_mbg/sheep_rdna_5001_500_graph.gfa -i hifi_paragon_ffastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211120_120456.rdna.fastq.gz"
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_3501_350.slurm.out --wrap="MBG -t 35 -k 3501 -w 350 -r 15000 --output-sequence-paths hifi_fparagon_mbg/sheep_rdna_3501_350_paths.gaf  --out hifi_fparagon_mbg/sheep_rdna_3501_350_graph.gfa -i hifi_paragon_ffastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_paragon_ffastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_paragon_ffastqs/m64015e_211120_120456.rdna.fastq.gz"
```

OK, finally found a filtration method that worked and graphs that appear normal. There are supposedly 5 45S arrays in sheep, but only two clear ones appear in the lower k value graphs. Let's pull reads and see if there are similarities between the graphs in the two lowest k values

```bash
sbatch -N 1 -n 2 --mem=8000 -p priority -q msn -o hifi_fparagon_mbg/2501_250_cluster_read_assigns.tab --wrap="python3 pull_reads_from_gaf.py hifi_fparagon_mbg/2501_250_cluster1_nodes.list hifi_fparagon_mbg/2501_250_cluster2_nodes.list hifi_fparagon_mbg/2501_250_cluster3_nodes.list hifi_fparagon_mbg/2501_250_cluster4_nodes.list hifi_fparagon_mbg/sheep_rdna_2501_250_paths.gaf"

sbatch -N 1 -n 2 --mem=8000 -p priority -q msn -o hifi_fparagon_mbg/3501_350_cluster_read_assigns.tab --wrap="python3 pull_reads_from_gaf.py hifi_fparagon_mbg/3501_350_cluster1_nodes.list hifi_fparagon_mbg/3501_350_cluster2_nodes.list hifi_fparagon_mbg/3501_350_cluster3_nodes.list hifi_fparagon_mbg/3501_350_cluster4_nodes.list hifi_fparagon_mbg/sheep_rdna_3501_350_paths.gaf"

mkdir read_lists
for i in 0 1 2 3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] == $ARGV[1]){print "$s[0]\n";}} close IN;' hifi_fparagon_mbg/3501_350_cluster_read_assigns.tab $i > read_lists/3501_350_cluster_${i}.list; done
for i in 0 1 2 3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] == $ARGV[1]){print "$s[0]\n";}} close IN;' hifi_fparagon_mbg/2501_250_cluster_read_assigns.tab $i > read_lists/2501_250_cluster_${i}.list; done

for i in 0 1 2 3; do for j in 0 1 2 3; do echo "ITERATION: $i $j"; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl read_lists/2501_250_cluster_${i}.list read_lists/3501_350_cluster_${j}.list; done; done
```

| cluster | 1 | 2 | 3 |  4|
|---------|---|---|---|---|
| 1		  | 50| 25|  0|  1|
| 2		  |  1| 10| 25| 25|
| 3		  |  0| 25|  0|  0|
| 4 	  |  0|  0|  0|  1|

Based on graph visualization, there are a couple anomalous nodes that could serve as breakpoints for separating clusters into distinct fragments. Cluster 4 in the k=3501 graph appears completely unique -- maybe add a different smaller cluster to compare as well? 

```bash
mkdir cluster_lists_two
mkdir cluster_lists_two/read_lists

sbatch -N 1 -n 2 --mem=8000 -p priority -q msn -o cluster_lists_two/read_lists/2501_250_cluster_read_assigns.tab --wrap="python3 pull_reads_from_gaf.py cluster_lists_two/2501_250_cluster1_nodes.list cluster_lists_two/2501_250_cluster2_nodes.list cluster_lists_two/2501_250_cluster3_nodes.list cluster_lists_two/2501_250_cluster4_nodes.list cluster_lists_two/2501_250_cluster5_nodes.list cluster_lists_two/2501_250_cluster6_nodes.list hifi_fparagon_mbg/sheep_rdna_2501_250_paths.gaf"

sbatch -N 1 -n 2 --mem=8000 -p priority -q msn -o cluster_lists_two/read_lists/3501_350_cluster_read_assigns.tab --wrap="python3 pull_reads_from_gaf.py cluster_lists_two/3501_350_cluster1_nodes.list cluster_lists_two/3501_350_cluster2_nodes.list cluster_lists_two/3501_350_cluster3_nodes.list cluster_lists_two/3501_350_cluster4_nodes.list cluster_lists_two/3501_350_cluster5_nodes.list cluster_lists_two/3501_350_cluster6_nodes.list hifi_fparagon_mbg/sheep_rdna_3501_350_paths.gaf"

for i in 0 1 2 3 4 5; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] == $ARGV[1]){print "$s[0]\n";}} close IN;' cluster_lists_two/read_lists/2501_250_cluster_read_assigns.tab $i > cluster_lists_two/read_lists/2501_250_cluster_${i}.list; done
for i in 0 1 2 3 4 5; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[1] == $ARGV[1]){print "$s[0]\n";}} close IN;' cluster_lists_two/read_lists/3501_350_cluster_read_assigns.tab $i > cluster_lists_two/read_lists/3501_350_cluster_${i}.list; done

for i in 0 1 2 3 4 5; do for j in 0 1 2 3 4 5; do echo "ITERATION: $i $j"; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl cluster_lists_two/read_lists/2501_250_cluster_${i}.list cluster_lists_two/read_lists/3501_350_cluster_${j}.list; done; done 
```

| cluster | 1 | 2 | 3 |  4|  5|  6|
|---------|---|---|---|---|---|---|
| 1		  | 30| 50|  5|  0|  0|  1|
| 2		  |  0|  0|  0| 10|  0|  0|
| 3		  |  0|  0|  1| 25| 60|  0|
| 4 	  |  0|  1|  1|  1|  0| 75|
| 5 	  |  0|  0|  1| 60|  0|  0|
| 6 	  |  0|  0|  0|  0|  0|  1|


## Retry with read mappings to the graph

OK, my previous attempt probably filtered away too many reads. I need to find the rDNA tigs in the graph and pull reads that were used in their construction. Then I can create the sub graph separately and try to tease it apart as done in the T2T manuscript.

The TIG to ID files are likely faulty because of naming conversion changes throughout the verkko cycle. Let's list where these files are located:

* /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/packages.tigName_to_ID.map
* /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/7-consensus/packages.tigName_to_ID.map


And the graph that I am using is here:

* /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/8-trio/unitig-popped-unitig-normal-connected-tip.noseq.gfa


First, let's see where my current crop of reads align:

> Ceres: /90daydata/rumen_longread_metagenome_assembly/analysis/rdna/graph_analysis

```bash
# preparing files
perl -e '<>; <>; while(<>){$_ =~ s/^\s+//; chomp; @s = split(/\s+/); print join("\t", @s) . "\n";}' < /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/7-consensus/packages.readName_to_ID.map > cons.readname_to_id.tab

perl -e '<>; <>; while(<>){$_ =~ s/^\s+//; chomp; @s = split(/\s+/); print join("\t", @s) . "\n";}' < /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/7-consensus/packages.tigName_to_ID.map > cons.tigname_to_id.tab

# grepping the read names
cat ../hifi_paragon_flists/*.list > combined_align_read_lists.list
wc -l combined_align_read_lists.list
21377 combined_align_read_lists.list

# These are the reads that are guaranteed to be in the rDNA regions. Let's use them as seeds for tig identification
module load miniconda samtools

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f cons.readname_to_id.tab -c 2 -l combined_align_read_lists.list -d '\t' | perl -lane 'print $F[1];' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -o contig_counts.tab

perl -e '<>; while(<>){print $_;}' < contig_counts.tab > contig_counts.fmt.tab
perl -lane 'print $F[0];' < contig_counts.fmt.tab > contig_counts.fmt.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f cons.tigname_to_id.tab -l contig_counts.fmt.list -c 0 -d '\t' -o cons.tigname_to_id.grepped.tab

python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f cons.tigname_to_id.grepped.tab -f contig_counts.fmt.tab -c 0 -m 1 -o rdna_contigs_grepped.tab
perl -lane 'print "$F[1]\t$F[0]\t$F[2]";' < rdna_contigs_grepped.tab > rdna_contigs_grepped.fmt.tab
python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f /90daydata/sheep_genome_assemblies/sergek/verkko_beta2/7-consensus/unitig-popped.fasta.fai -f rdna_contigs_grepped.fmt.tab -c 0 -m 1 -m 2 -o rdna_contigs_grepped.fmt.plslen.tab
perl -lane 'if($F[3] eq "-"){next;}else{print "$F[0]\t$F[1]\t$F[3]\t$F[4]";}' < rdna_contigs_grepped.fmt.plslen.tab | sort -k4nr > rdna_contigs_grepped.fmt.plslen.fixed.tab

# This may miss quite a few potential connections, but let's run with a criteria of a unitig less than 100kb and having at least 2 candidate conservative rRNA reads mapping
perl -lane 'if($F[1] < 160000 && $F[3] > 1){print $F[0];}' < rdna_contigs_grepped.fmt.plslen.fixed.tab > candidate_rdna_unitigs.list  #<- ~ 1/3rd of the total

# Now let's pull all of the reads that were used to construct these nodes
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f cons.tigname_to_id.tab -c 1 -l candidate_rdna_unitigs.list | perl -lane 'print $F[0];' > candidate_rdna_unitigs.ids.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f cons.readname_to_id.tab -l candidate_rdna_unitigs.ids.list -c 1 | perl -lane 'print $F[2];' | sort | uniq > candidate_rdna_reads.list

wc -l combined_align_read_lists.list candidate_rdna_reads.list
  21377 combined_align_read_lists.list
   9105 candidate_rdna_reads.list
  30482 total
 
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl combined_align_read_lists.list candidate_rdna_reads.list
File Number 1: combined_align_read_lists.list
File Number 2: candidate_rdna_reads.list
Set     Count
1       16699
1;2     4678
2       4427

# Creating read association list with the main assembly:
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f cons.readname_to_id.tab -l candidate_rdna_unitigs.ids.list -c 1 | perl -lane 'print "$F[1]\t$F[2]";' > candidate_rdna_reads.plusctgids.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $n = $h{$s[0]}; print "$s[1]\t$n\n";} close IN;' cons.tigname_to_id.tab candidate_rdna_reads.plusctgids.tab > candidate_rdna_reads.plusctgnames.tab

for i in hifi_rdna_mbg/*.gaf; do name=`basename $i | cut -d'_' -f1,2,3,4`; echo $name; perl -lane '@s = split(/[<>]/, $F[5]); print "$F[0]\t$s[1]";' < $i > ${name}.rnameassoc.tab; done

for i in hifi_rdna_mbg/*.gaf; do name=`basename $i | cut -d'_' -f1,2,3,4`; echo $name; python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f ${name}.rnameassoc.tab -f candidate_rdna_reads.plusctgnames.tab -c 0 -m 1 -o ${name}.allassoc.tab; done
```

Hmm, that's a very small amount of constitutive reads and I'm very surprised to see more unique reads in the candidate dataset. Lets try to assemble only with the candidate read set, and then with the superset to see if that resolves the graph.


```bash
module load miniconda samtools minimap2
conda activate /project/rumen_longread_metagenome_assembly/environments/seaborn/

mkdir hifi_rdna_fastqs
for i in `ls /90daydata/sheep_genome_assemblies/Churro_x_Friesian/hifi/*.gz`; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=25000 --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q $i -l candidate_rdna_reads.list -o hifi_rdna_fastqs/$name.rdna.fastq; gzip hifi_rdna_fastqs/$name.rdna.fastq"; done

# OK, now time to try to make the graphs again
mkdir hifi_rdna_mbg
conda activate /project/rumen_longread_metagenome_assembly/environments/verkko
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_2501_250.slurm.out --wrap="MBG -t 35 -k 2501 -w 250 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_2501_250_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_2501_250_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_5000_500.slurm.out --wrap="MBG -t 35 -k 5001 -w 500 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_5001_500_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_5001_500_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"
sbatch -N 1 -n 35 --mem=300000 -p priority -q msn -o mbg_fparagon_3501_350.slurm.out --wrap="MBG -t 35 -k 3501 -w 350 -r 15000 --output-sequence-paths hifi_rdna_mbg/sheep_rdna_3501_350_paths.gaf  --out hifi_rdna_mbg/sheep_rdna_3501_350_graph.gfa -i hifi_rdna_fastqs/m54337U_210601_165737.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210610_143100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210611_233100.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210621_195743.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210623_154734.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210716_180058.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210718_025906.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210726_175325.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_210728_025148.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211026_163028.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211028_032534.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211104_191436.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211106_060943.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211118_170215.rdna.fastq.gz -i hifi_rdna_fastqs/m54337U_211120_035803.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211119_011043.rdna.fastq.gz -i hifi_rdna_fastqs/m64015e_211120_120456.rdna.fastq.gz"
```


### Gaur Pied analysis

> Ceres: /90daydata/gaur_genome_assembly/Gaur_x_Pied/verkko_DC_test3

```bash
module load repeatmasker/4.1.0 java/8 minimap2
for i in haplotype1 haplotype2; do echo $i; sbatch -N 1 -n 50 --mem=55000 -p priority -q msn --wrap="RepeatMasker -pa 50 -q -species cattle -no_is -gff assembly.$i.fasta"; done

for i in haplotype1 haplotype2; do echo $i; sbatch -N 1 -n 3 -p priority -q msn --mem=45000 --wrap="java -Xmx45g -jar /project/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f assembly.$i.fasta -o $i.gaps.bed -s $i.gaps.stats"; done

for i in assembly.haplotype?.fasta.out; do echo $i; perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < $i > $i.bed; done

python3 ~/python_toolchain/sequenceData/checkTelomereFromRMOutput.py -b assembly.haplotype1.fasta.out.bed -f assembly.haplotype1.fasta.fai -o assembly.haplotype1.fasta.out.telomeres
python3 ~/python_toolchain/sequenceData/checkTelomereFromRMOutput.py -b assembly.haplotype2.fasta.out.bed -f assembly.haplotype2.fasta.fai -o assembly.haplotype2.fasta.out.telomeres

for i in haplotype1 haplotype2; do echo $i; perl -lane 'if($F[1] > 20000000){print "$F[2];$F[3]";}' < assembly.$i.fasta.out.telomeres | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m; done
haplotype1
|Entry       | Value|
|:-----------|-----:|
|True;False  |    19|
|False;False |     4|
|True;True   |     4|
|False;True  |     3|
haplotype2
|Entry       | Value|
|:-----------|-----:|
|False;False |    14|
|True;False  |     8|
|False;True  |     7|
|True;True   |     1|

for i in haplotype1 haplotype2; do echo $i; sbatch -N 1 -n 3 --mem=9000 -p priority -q msn --wrap="minimap2 -x asm20 cattle_centromere.fa assembly.$i.fasta > centromere.$i.paf"; done
```

### Churro x Friesian analysis

> Ceres: /90daydata/sheep_genome_assemblies/Churro_x_Friesian/ 

```bash
module load repeatmasker/4.1.0 java/8
sbatch -N 1 -n 50 --mem=55000 -p priority -q msn --wrap="RepeatMasker -pa 50 -q -species sheep -no_is -gff assembly.fasta"

sbatch -N 1 -n 3 -p priority -q msn --mem=45000 --wrap="java -Xmx45g -jar /project/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f assembly.fasta -o assembly.gaps.bed -s assembly.gaps.stats"

perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < assembly.fasta.out > assembly.fasta.out.bed
python3 ~/python_toolchain/sequenceData/checkTelomereFromRMOutput.py -b assembly.fasta.out.bed -f assembly.fasta.fai -o assembly.fasta.out.telomeres

perl -lane 'if($F[1] > 20000000){print "$F[2];$F[3]";}' < assembly.fasta.out.telomeres | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m
|Entry       | Value|
|:-----------|-----:|
|False;False |    91|
|True;False  |     7|
|False;True  |     2|
|True;True   |     1|
```

### Kiko x Saanen analysis

> Ceres: /90daydata/sheep_genome_assemblies/Kiko_x_Saanen/verkko_DCken4

```bash
module load repeatmasker/4.1.0 java/8
for i in haplotype1 haplotype2; do echo $i; sbatch -N 1 -n 50 --mem=55000 -p priority -q msn --wrap="RepeatMasker -pa 50 -q -species goat -no_is -gff assembly.$i.fasta"; done

for i in haplotype1 haplotype2; do echo $i; sbatch -N 1 -n 3 -p priority -q msn --mem=45000 --wrap="java -Xmx45g -jar /project/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f assembly.$i.fasta -o $i.gaps.bed -s $i.gaps.stats"; done

for i in assembly.haplotype?.fasta.out; do echo $i; perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < $i > $i.bed; done

for i in haplotype1 haplotype2; do echo $i; samtools faidx assembly.$i.fasta; python3 ~/python_toolchain/sequenceData/checkTelomereFromRMOutput.py -b assembly.$i.fasta.out.bed -f assembly.$i.fasta.fai -o assembly.$i.fasta.out.telomeres; perl -lane 'if($F[1] > 20000000){print "$F[2];$F[3]";}' < assembly.$i.fasta.out.telomeres | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -m; done
haplotype1
|Entry       | Value|
|:-----------|-----:|
|True;False  |    19|
|False;False |     7|
|True;True   |     4|
|False;True  |     1|
haplotype2
|Entry       | Value|
|:-----------|-----:|
|True;False  |    24|
|False;False |     6|
|True;True   |     1|
```