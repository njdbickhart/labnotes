# MetaFlye analysis of sheep dataset
---
*12/5/2019*

These are my notes and plans for validating the MetaFlye data for use in the reviewer rebuttals.

## Table of Contents
* [Analysis plan](#analysis)
* [Running the analysis](#running)
	* [CCS read alignments](#ccs)
	* [WGS read alignments](#wgs)
	* [ORF annotation](#orf)
	* [Blobtools](#blobtools)
	* [Network plot and read alignment](#network)
	* [Hi-C and metabat binning](#binning)
	* [Checkm plotting](#checkm)
	* [Bin3c off-diagonal analysis and Hi-C inter-contig links](#offdiag)
	* [BUSCO analysis of Eukaryotic contigs](#busco)
	* [DESMAN strain inference](#desman)
* [Generating figures and tables](#generating)
	* [Kingdom-level and GC bin contig size plot](#kingdomlevel)
	* [ORF prediction classification](#orfprediction)

<a name="analysis"></a>
## Analysis plan

I will need to divy up tasks and run them in order to generate the results needed.

* Run CheckM on both assemblies
* Run BlobTools on both assemblies
	* Taxify contigs and derive contig taxonomic profiles
	* Derive a consensus taxonomic affiliation per bin and plot CheckM results
* Run Minimap2 using CCS reads on both assemblies
* Align Hi-C reads to both assemblies
	* Run the viral association script
	* Count number of split-read or hard-clipped CCS reads
	* Adapt DESMAN to pick strains from the alignments
* Run prodigal on the contigs

<a name="running"></a>
## Running the analysis

<a name="ccs"></a>
#### CCS read alignments

> Ceres: /project/forage_assemblies/sheep_project

```bash
module load minimap2/2.6
for i in *.tar.gz; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=1 -p msn --wrap="tar -xvf $i"; done

# OK, so the contigs are not properly binned, but given their sizes, this is probably OK.
# Let me start the long-term alignment first so that I can analyze this next week
# I will use this for blobtools "coverage" estimates as well
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch --nodes=1 --mem=30000 --ntasks-per-node=3 -p msn --wrap="minimap2 -x map-pb $i /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq.gz > $i.ccs.paf"; done

# The blobtools coverage estimates didn't work very well with the CCS reads in sam format. Maybe I can resolve this with scripting?

# Total size of the CCS reads
perl -e '$c = 0; while(<>){$s = <>; chomp($s); $c += length($s); <>; <>;} print "$c\n";' < /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq
51,928,785,035		<- 51X coverage of a assembly size of 1 Gbp

# Rerunning on v2 assembly
sbatch --nodes=1 --mem=30000 --ntasks-per-node=3 -p msn -q msn --wrap="minimap2 -x map-pb flye2.contigs.fasta /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq > flye2.contigs.fasta.ccs.paf"

# Some exploratory data analysis on CCS read alignments
# Just getting all reads that map multiple times and counts of their mapping start and end locations (to try to track repeats)
perl -e '%data; %keeps; while(<>){chomp; @s = split(/\t/); $data{$s[0]}->{$s[2]}->{$s[3]} += 1; $keeps{$s[0]} += 1;} @rs; foreach $j (sort {$keeps{$b} <=> $keeps{$a}} keys(%keeps)){if($keeps{$j} > 1){push(@rs, $j);}} foreach $chr (@rs){foreach $start (keys(%{$data{$chr}})){foreach $end (keys(%{$data{$chr}->{$start}})){$val = $data{$chr}->{$start}->{$end}; print "$chr\t$start\t$end\t$val\n";}}}' < canu.contigs.fasta.ccs.paf > canu.contigs.fasta.ccs.mult.algn.bed

cat canu.contigs.fasta.ccs.mult.algn.bed | cut -f4 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   2672241
Sum:    4019981
Minimum 1
Maximum 12
Average 1.504348
Median  1
Standard Deviation      1.009274
Mode(Highest Distributed Value) 1

# OK, time to translate this to a new script.
perl identify_multimapped_ccs_reads.pl canu.contigs.fasta.ccs.paf canu.contigs.fasta.ccs.mult.algn.bed canu.contigs.fasta.ccs.mult.algn.graph

cat canu.contigs.fasta.ccs.mult.algn.graph | cut -f3 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   34001
Sum:    325368
Minimum 1
Maximum 4555
Average 9.569366
Median  2
Standard Deviation      66.898076
Mode(Highest Distributed Value) 1

# Let's do this on Flye now
perl identify_multimapped_ccs_reads.pl flye2.contigs.fasta.ccs.paf flye2.contigs.fasta.ccs.mult.algn.bed flye2.contigs.fasta.ccs.mult.algn.graph

cat flye2.contigs.fasta.ccs.mult.algn.graph | cut -f3 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   26069
Sum:    148533
Minimum 1
Maximum 1302
Average 5.697687
Median  2
Standard Deviation      21.533184
Mode(Highest Distributed Value) 1

# Graph connections greater than 51 X coverage
perl -lane 'if($F[2] > 51){print $_;}' < flye2.contigs.fasta.ccs.mult.algn.graph > flye2.contigs.fasta.ccs.mult.gt51.graph
perl -lane 'if($F[2] > 51){print $_;}' < canu.contigs.fasta.ccs.mult.algn.graph > canu.contigs.fasta.ccs.mult.gt51.graph

## The CCS read counts were highly inflated by repetitive region alignments. Let's see if I can remove those multi-mapped alignments
sort -k1 flye2.contigs.fasta.ccs.mult.algn.bed | perl -lane 'print "$F[0]-$F[1]-$F[2]";' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 -o flye2.contigs.fasta.ccs.mult.ctg.counts

# Actually, I found errors in the logic in the script above. I rewrote it and ran it. I think it's working better now.
# I added a third output file format that links the read mapping coords together in a link for later filtering
perl identify_multimapped_ccs_reads.pl canu.contigs.fasta.ccs.paf canu.contigs.fasta.ccs.mult.algn.bed canu.contigs.fasta.ccs.mult.algn.graph canu.contigs.fasta.ccs.mult.algn.links

perl -lane 'if($F[-1] > 5){print $_;}' < canu.contigs.fasta.ccs.mult.algn.links | wc -l
107

perl identify_multimapped_ccs_reads.pl flye2.contigs.fasta.ccs.paf flye2.contigs.fasta.ccs.mult.algn.bed flye2.contigs.fasta.ccs.mult.algn.graph flye2.contigs.fasta.ccs.mult.algn.links

# Removing links between legitimate bin3c bins
for i in flye2 canu; do echo $i; perl -e 'chomp @ARGV; %bbins; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); push(@{$bbins{$s[1]}}, $s[0]);} close IN; my %graph; foreach $c (keys(%bbins)){ @working = @{$bbins{$c}}; for($i = 0; $i < scalar(@working); $i++){for($j = 0; $j < scalar(@working); $j++){$graph{$working[$i]}->{$working[$j]} = 1;}}} open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($graph{$s[0]}->{$s[3]}) || exists($graph{$s[1]}->{$s[0]})){next;}else{print join("\t", @s) . "\n";}} close IN;' $i.bin3c.bins.tab $i.contigs.fasta.ccs.mult.algn.links > $i.contigs.fasta.ccs.mult.algn.filt.links; done

# A script to link the reads with hi-c
for i in flye2 canu; do sbatch --nodes=1 --mem=1000 --ntasks-per-node=2 -p msn -q msn --wrap="perl add_hic_evidence.pl $i.contigs.fasta.ccs.mult.algn.filt.links ${i}_hic/63/63.sorted.merged.bam $i.contigs.fasta.ccs.mult.algn.links.evd"; done

# sorting and condensing
perl sort_evidence.pl flye2.contigs.fasta.ccs.mult.algn.links.evd flye2.contigs.fasta.ccs.mult.algn.links.cond.evd
perl sort_evidence.pl canu.contigs.fasta.ccs.mult.algn.links.evd canu.contigs.fasta.ccs.mult.algn.links.cond.evd

# HIC read X coverage
python3 ~/python_toolchain/sequenceData/bamStats.py -f flye2_hic/63/63.sorted.merged.bam
#Genome size: {'1,469,466,257'} and avg read len: {'150'}. Cropping similar columns
File    TotalReads      MappedReads     MapPerc RawXCov MapXCov AvgMapChrXCov   AvgChrMapPerc
63.sorted.merged.bam    221,409,724     191,177,768     0.863   22.601  19.515  15.632  0.977
```

Now to try to find relevant statistics on the data I've collated.

```R
library(dplyr)
fevd <- read.delim("flye2.contigs.fasta.ccs.mult.algn.links.cond.evd", header=FALSE)
fevd <- mutate(fevd, LEN = V3 - V2) %>% mutate(COV = (V8 * 150) / LEN)
nrow(fevd[fevd$V7 > 5 & fevd$COV > 1,])
222

# Using Hi-C X coverage estimates above for filtering
nrow(fevd[fevd$V7 > 5 & fevd$COV > 20,])
[1] 9
write.table(fevd[fevd$V7 > 5 & fevd$COV > 20,], file="flye2.ccs.fevd.likelies", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

cevd <- read.delim("canu.contigs.fasta.ccs.mult.algn.links.cond.evd", header=FALSE)
cevd <- mutate(cevd, LEN = V3 - V2) %>% mutate(COV = (V8 * 150) / LEN)
nrow(cevd[cevd$V7 > 5 & cevd$COV > 20,])

write.table(cevd[cevd$V7 > 5 & cevd$COV > 20,], file="canu.ccs.fevd.likelies", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
```

Finally, let's validate each one. I want to see how many CCS reads support the contig as-is and how many are divergent. This will show if there's a problem with chimerism or repeats.

```bash
perl count_read_map_ratio.pl canu.ccs.fevd.likelies canu.contigs.fasta.ccs.paf canu.ccs.fevd.rcounts
for i in *.rcounts*; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 5 -d '\t' -m; done
canu.ccs.fevd.rcounts_1.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |  5352|
|PARTIAL |   987|
|NORM    |     1|
canu.ccs.fevd.rcounts_2.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |  4794|
|PARTIAL |   982|
|NORM    |     1|
canu.ccs.fevd.rcounts_3.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |   145|
|PARTIAL |    56|
|NORM    |    20|
canu.ccs.fevd.rcounts_4.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |   172|
|PARTIAL |    46|
|NORM    |    20|
canu.ccs.fevd.rcounts_5.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |   105|
|PARTIAL |    44|
|NORM    |    20|

perl count_read_map_ratio.pl flye2.ccs.fevd.likelies flye2.contigs.fasta.ccs.paf flye2.ccs.fevd.rcounts
for i in flye2*.rcounts*; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 5 -d '\t' -m; done
flye2.ccs.fevd.rcounts_1.out
|Entry   | Value|
|:-------|-----:|
|PARTIAL |  2159|
|MAPQ0   |    62|
flye2.ccs.fevd.rcounts_2.out
|Entry   | Value|
|:-------|-----:|
|PARTIAL | 29326|
|MAPQ0   |   182|
flye2.ccs.fevd.rcounts_3.out
|Entry   | Value|
|:-------|-----:|
|PARTIAL | 43376|
|MAPQ0   |   208|
flye2.ccs.fevd.rcounts_4.out
|Entry   | Value|
|:-------|-----:|
|NORM    |  1826|
|PARTIAL |  1380|
|MAPQ0   |   568|
flye2.ccs.fevd.rcounts_5.out
|Entry   | Value|
|:-------|-----:|
|PARTIAL |  6913|
|NORM    |  5904|
|MAPQ0   |     1|
flye2.ccs.fevd.rcounts_6.out
|Entry   | Value|
|:-------|-----:|
|NORM    | 29018|
|PARTIAL | 16273|
|MAPQ0   |  2148|
flye2.ccs.fevd.rcounts_7.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |   171|
|PARTIAL |    61|
|NORM    |    35|
flye2.ccs.fevd.rcounts_8.out
|Entry   | Value|
|:-------|-----:|
|MAPQ0   |  1787|
|NORM    |  1521|
|PARTIAL |   952|
flye2.ccs.fevd.rcounts_9.out
|Entry   | Value|
|:-------|-----:|
|PARTIAL |  2304|
|NORM    |  1477|
|MAPQ0   |   821|

```

Let's condense this all down into something that's useful.

#### Canu likely chimeric

| Chr | Start | End | Chr | Start | End | CCs reads | Hi-C link cov | MapQ0 | Partial | Norm |
| :---| :---| :----| :---| :--- | :---- | :--------| :-------------| :---- | :----| :----|
|tig00022560|	26340|	39630|	tig00064035|	2508|	18524|	24 |20.1 | 5352 | 987 | 1|
| tig00053493|	84365|	87068|	tig00080874 |	4808|	7510|	9 | 90 | 172 | 56 | 20 | 

#### Flye2 likely chimeric

| Chr | Start | End | Chr | Start | End | CCs reads | Hi-C link cov | MapQ0 | Partial | Norm |
| :---| :---| :----| :---| :--- | :---- | :--------| :-------------| :---- | :----| :----|
| contig_1050 |	1996938 |	1997074	| contig_11674	|323	|21004	|96 | 29.7 | 62 | 2159 | 0  |
| contig_15841|	164213|	165317|	contig_6761|	30609|	40086|	443 | 167 | 182 | 29700 | 0 |
| contig_15841|	209118|	210297|	contig_6761|	30608|	40086|	4200 | 22.5 | 208 | 43000| 0 |
| contig_14057|	17|	5617|	contig_9755|	0|	8524|	121 | 97.0 | 568 | 1380 | 1826 |
| contig_12232|	0|	531|	contig_468|	8|	18407|	1091 | 171.0 | 1 | 6900 | 5900|
| contig_14333|	18|	2456|	contig_9516|	6|	7074|	1869 | 20.7 | 2148 | 16000 | 29000|
| contig_12644|	5942|	10660|	contig_3716|	10|	6954|	17 | 75.2 | 171 | 61 | 35 |
| contig_15873|	12|	562|	contig_9755|	177|	8524|	15 | 162.0 | 1787 | 952 | 1521 |
| contig_14854	|15	|1850	|contig_9755	|44	|8524	|42 | 114.4 | 821 | 2304 | 1477 |


<a name="wgs"></a>
#### WGS read alignments

```bash
# BWA indexing
for i in canu.contigs.fasta flye_contigs.fasta ; do echo $i; sbatch --nodes=1 --mem=15000 --ntasks-per-node=1 -p msn -q msn --wrap="module load bwa; bwa index $i"; done

# queueing alignment jobs
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye_contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/hic_links.tab -b flye_hic -p short -q memlimit -m
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye_contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/wgs_reads.tab -b flye_wgs -p short -q memlimit -m

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/canu.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/hic_links.tab -b canu_hic -p short -q memlimit -m
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/canu.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/wgs_reads.tab -b canu_wgs -p short -q memlimit -m

# Rerunning on v2 assembly
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn -q msn --wrap="bwa index flye2.contigs.fasta"

sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p msn -q msn --wrap="bwa mem -5SP -t 3 flye2.contigs.fasta /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R1_001.fastq.gz /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R2_001.fastq.gz | samtools view -F 0x904 -bS - | samtools sort -n -T flye2hic.tmp -o flye2.hiclinks.bam -@ 3 -"
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye2.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/wgs_reads.tab -b flye2_wgs -p short -q memlimit -m
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -f /project/forage_assemblies/sheep_project/flye2.contigs.fasta -t /project/rumen_longread_metagenome_assembly/sheep_poop/hic_links.tab -b flye2_hic -p short -q memlimit -m
```
<a name="orf"></a>
#### ORF annotation

I am going to run both assemblies through Prodigal. The ORF data will be extremely useful downstream.

```bash
module load prodigalorffinder/2.6.3

for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p msn --wrap="prodigal -a $i.prod.prottrans -c -d $i.prod.genenuc -f gff -i $i -o $i.prod.out -p meta"; done

for i in canu flye; do grep '>' $i.contigs.fasta.prod.prottrans |  perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > $i.contigs.prod.shortform.tab; done

# Rerunning on v2 assembly
sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p msn --wrap="prodigal -a flye2.contigs.fasta.prod.prottrans -c -d flye2.contigs.fasta.prod.genenuc -f gff -i flye2.contigs.fasta -o flye2.contigs.fasta.prod.out -p meta"
```
<a name="blobtools"></a>
#### Blobtools

```bash
module load miniconda
conda activate /KEEP/rumen_longread_metagenome_assembly/blobtools

# Diamond run before blobtools
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query $i --db /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o $i.diamondout.tsv"; done

# NOTE: I had to clone the blobtools repo because Conda screwed up the libraries for the program
# Blobtools taxify
for i in canu.contigs.fasta flye_contigs.fasta; do echo $i; sbatch -t 2-0 --nodes=1 --mem=20000 --ntasks-per-node=3 -p msn -q msn --wrap="./blobtools/blobtools taxify -f $i.diamondout.tsv -m /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.taxids -s 0 -t 2 -o ${i}_unip"; done

# Copying the database file
cp /project/rumen_longread_metagenome_assembly/assemblies/protists/blob_ncbi.db ./

# Blobtools create
for i in canu flye; do sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools create -i $i.contigs.fasta -b ${i}_wgs/63/63.sorted.merged.bam -t $i.contigs.fasta_unip.$i.contigs.fasta.diamondout.tsv.taxified.out -o $i.blobtools --db blob_ncbi.db"; done

# Blobtools view and plot
for i in canu flye; do sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i $i.blobtools.blobDB.json --notitle -r superkingdom --format pdf -o ${i}_supkingdom"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i $i.blobtools.blobDB.json --notitle --format pdf -r phylum -o ${i}_phylum"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools view -i $i.blobtools.blobDB.json -o ${i}_table -r all"; done

mkdir blob_plots
mv *.stats.txt ./blob_plots/
mv *.png ./blob_plots/

# Taking stock of superkingdoms
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f canu_table.canu.blobtools.blobDB.table.txt -c 5 -i '#' -d '\t' -m
|Entry     | Value|
|:---------|-----:|
|Bacteria  | 18821|
|Archaea   |   881|
|Eukaryota |   296|
|Viruses   |   208|
|no-hit    |    71|

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f flye_table.flye.blobtools.blobDB.table.txt -c 5 -i '#' -d '\t' -m
|Entry     | Value|
|:---------|-----:|
|Bacteria  | 15751|
|Archaea   |   584|
|no-hit    |   249|
|Eukaryota |   205|
|Viruses   |   123|


# Rerun on v2 assembly
sbatch -t 2-0 -p msn -q msn --nodes=1 --ntasks-per-node=30 --mem=100000 --wrap="diamond blastx --query flye2.contigs.fasta --db /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.diamond.dmnd --threads 29 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o flye2.contigs.fasta.diamondout.tsv";

sbatch -t 2-0 --nodes=1 --mem=20000 --ntasks-per-node=3 -p msn -q msn --wrap="./blobtools/blobtools taxify -f flye2.contigs.fasta.diamondout.tsv -m /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.taxids -s 0 -t 2 -o flye2.contigs.fasta_unip"

sbatch -t 2-0 --dependency=afterok:1427747 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools create -i flye2.contigs.fasta -b flye2_wgs/63/63.sorted.merged.bam -t flye2.contigs.fasta_unip.flye2.contigs.fasta.diamondout.tsv.taxified.out -o flye2.blobtools --db blob_ncbi.db"

sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i flye2.blobtools.blobDB.json --notitle --format pdf -r superkingdom -o flye2_supkingdom"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools plot -i flye2.blobtools.blobDB.json --notitle -r phylum --format pdf -o flye2_phylum"; sbatch -t 2-0 --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn -q msn --wrap="./blobtools/blobtools view -i flye2.blobtools.blobDB.json -o flye2_table -r all"
```
<a name="network"></a>
#### Network plot and read alignment

I added a rudimentary network plot to my pipeline script. Let's see if it works!

```bash
source activate /KEEP/rumen_longread_metagenome_assembly/seaborn/
module load samtools/1.9 minimap2/2.6

perl -lane 'if($F[0] =~ /^#/){next;} if($F[5] eq "Viruses"){print "$F[0]\t$F[1]";}' < canu_table.canu.blobtools.blobDB.table.txt > canu.viral.contigs.list
perl -lane 'if($F[0] =~ /^#/){next;} if($F[5] eq "Viruses"){print "$F[0]\t$F[1]";}' < flye_table.flye.blobtools.blobDB.table.txt > flye.viral.contigs.list

# Preparing list of viral contig fastas
perl -lane 'print $F[0];' < canu.viral.contigs.list | xargs -I {} samtools faidx canu.contigs.fasta {} >> canu.viral.contigs.fa
perl -lane 'print $F[0];' < flye.viral.contigs.list | xargs -I {} samtools faidx flye.contigs.fasta {} >> flye.viral.contigs.fa

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p msn -q msn -J canuflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a canu.contigs.fasta -g canu.viral.contigs.fa -b canu_table.canu.blobtools.blobDB.table.txt -i canu_hic/63/63.sorted.merged.bam -v canu.viral.contigs.list -l /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq -m /software/7/apps/minimap2/2.6/minimap2 -o canu.contigs.vassoc"

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p msn -q msn -J flyeflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a flye.contigs.fasta -g flye.viral.contigs.fa -b flye_table.flye.blobtools.blobDB.table.txt -i flye_hic/63/63.sorted.merged.bam -v flye.viral.contigs.list -l /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq -m /software/7/apps/minimap2/2.6/minimap2 -o flye.contigs.vassoc"

# For flye version 2
perl -lane 'if($F[0] =~ /^#/){next;} if($F[5] eq "Viruses"){print "$F[0]\t$F[1]";}' < flye2_table.flye2.blobtools.blobDB.table.txt > flye2.viral.contigs.list

perl -lane 'print $F[0];' < flye2.viral.contigs.list | xargs -I {} samtools faidx flye2.contigs.fasta {} >> flye2.viral.contigs.fa

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p msn -q msn -J flyeflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a flye2.contigs.fasta -g flye2.viral.contigs.fa -b flye2_table.flye2.blobtools.blobDB.table.txt -i flye2_hic/63/63.sorted.merged.bam -v flye2.viral.contigs.list -l /project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_CCS.fastq -m /software/7/apps/minimap2/2.6/minimap2 -o flye2.contigs.vassoc"
```

<a name="binning"></a>
#### Hi-C and metabat binning

OK, I think that I will run Aaron Darling's [Bin3c tool](https://github.com/cerebis/bin3C) to save some time. The instructions listed a very important point: mate-pairing Hi-C reads is likely to fail because of the distances between pairs! I should use their recommended alignment strategy to avoid faulty alignments.

> Ceres:

```bash
module load bwa/0.7.17 samtools/1.9 miniconda

for i in canu flye; do sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p msn -q msn --wrap="bwa mem -5SP -t 3 $i.contigs.fasta /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R1_001.fastq.gz /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R2_001.fastq.gz | samtools view -F 0x904 -bS - | samtools sort -n -T $i.tmp -o $i.hiclinks.bam -@ 3 -"; done

# Getting HIC read statistics
for i in flye2 canu; do echo $i; sbatch --nodes=1 --mem=1000 --ntasks-per-node=2 -p msn -q msn --wrap="perl calc_proportion_intercontig.pl ${i}_hic/63/63.sorted.merged.bam $i.hic_readcounts.tab"; done

head *readcounts.tab
==> canu.hic_readcounts.tab <==
Total:  215424778
Mapped: 185004106
Unmapped:       30420672
MapRatio:       0.8588
INTRACTG:       170551154
INTERCTG:       11523430
INTRARatio:     0.7917
INTERRatio:     0.0535

==> flye2.hic_readcounts.tab <==
Total:  215424778
Mapped: 185192822
Unmapped:       30231956
MapRatio:       0.8597
INTRACTG:       170609094
INTERCTG:       11654976
INTRARatio:     0.7920
INTERRatio:     0.0541

# Setup bin3C
git clone --recursive https://github.com/cerebis/bin3C

conda create --prefix=/KEEP/rumen_longread_metagenome_assembly/bin3C python=2.7
conda activate /KEEP/rumen_longread_metagenome_assembly/bin3C


# OK, I don't know which RE was used, but I'm guessing it's Sau3AI for educational purposes here
python bin3C/bin3C.py mkmap -e Sau3AI -v canu.contigs.fasta canu.hiclinks.bam canu.hiclinks.bin3c_out
	IOError: BAM file must be sorted by read name

# URRRRRRGH! 
# Sorting by read name
for i in canu flye; do echo $i; sbatch --nodes=1 --mem=10000 --ntasks-per-node=4 -p msn -q msn --wrap="samtools sort -T $i.tmp -n -o $i.hiclinks.rname.bam -@ 3 $i.hiclinks.bam"; done

python bin3C/bin3C.py mkmap -e Sau3AI -v canu.contigs.fasta canu.hiclinks.rname.bam canu_bin3c

sbatch --nodes=1 --mem=16000 --ntasks-per-node=2 -p msn -q msn --wrap="python bin3C/bin3C.py cluster -v --no-spades canu_bin3c/contact_map.p.gz canu_bin3c_cluster"

sbatch --nodes=1 --mem=16000 --ntasks-per-node=2 -p msn -q msn --wrap="python bin3C/bin3C.py mkmap -e Sau3AI -v flye.contigs.fasta flye.hiclinks.rname.bam flye_bin3c"

sbatch --nodes=1 --dependency=afterany:1397741 --mem=16000 --ntasks-per-node=2 -p msn -q msn --wrap="python bin3C/bin3C.py cluster -v --no-spades flye_bin3c/contact_map.p.gz flye_bin3c_cluster"

#### Metabat2
module load metabat/2.12.1 miniconda/3.6 usearch/11.0.667 hmmer3/3.2.1
# Getting coverage files
for i in flye canu; do sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p brief-low -q memlimit --wrap="jgi_summarize_bam_contig_depths --outputDepth $i.wgs.cov ${i}_wgs/63/63.sorted.merged.bam"; done

for i in flye canu; do sbatch --nodes=1 -p brief-low -q memlimit --ntasks-per-node=15 --mem=20000 --wrap="metabat2 -i $i.contigs.fasta -a $i.wgs.cov -o ${i}_meta/$i.bin -t 15 -v"; done


#### DAS_TOOL
# First associate bins with contigs
for i in flye canu; do echo $i; for j in ${i}_meta/*.fa; do bin=`basename $j | cut -d'.' -f2,3`; echo $bin; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; if($_ =~ /^>/){$_ =~ s/>//; print "$_\t$ARGV[1]\n";}} close IN;' $j $bin >> $i.metabat.bins.tab; done; done
for i in canu flye; do echo $i; perl -lane 'print "$F[2]\t$F[0]";' < ${i}_bin3c_bins.ctgassoc > $i.bin3c.bins.tab; done

module load diamond/0.9.28 usearch/11.0.667
conda activate /KEEP/rumen_longread_metagenome_assembly/das_tool/
for i in canu flye; do echo $i; mkdir ${i}_dastool; sbatch --nodes=1 --mem=25000 --ntasks-per-node=4 -p msn -q msn --wrap="DAS_Tool --search_engine 'diamond' -i ${i}.bin3c.bins.tab,${i}.metabat.bins.tab -l bin3c,metabat -c $i.contigs.fasta -o ${i}_dastool/$i.das -t 4 --write_bins 1"; done


# For the second version
sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p msn -q msn --wrap="bwa mem -5SP -t 3 flye2.contigs.fasta /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R1_001.fastq.gz /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R2_001.fastq.gz | samtools view -F 0x904 -bS - | samtools sort -n -T flye2.tmp -o flye2.hiclinks.bam -@ 3 -"

sbatch --nodes=1 --mem=16000 --ntasks-per-node=2 -p msn -q msn --wrap="python bin3C/bin3C.py mkmap -e Sau3AI -v flye2.contigs.fasta flye2.hiclinks.bam flye2_bin3c"

sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p brief-low -q memlimit --wrap="jgi_summarize_bam_contig_depths --outputDepth flye2.wgs.cov flye2_wgs/63/63.sorted.merged.bam"

sbatch --dependency=afterok:1444378 --nodes=1 -p brief-low -q memlimit --ntasks-per-node=15 --mem=20000 --wrap="metabat2 -i flye2.contigs.fasta -a flye2.wgs.cov -o flye2_meta/flye2.bin -t 15 -v"

# I'm going to try to use DasTool only on the bin3c bins to save time
python bin3C/bin3C.py cluster -v --no-spades flye2_bin3c/contact_map.p.gz flye2_bin3c_cluster
perl pull_bin3c_ctg_name_tax.pl flye2_bin3c_cluster/fasta/ flye2_table.flye2.blobtools.blobDB.table.txt flye2_bin3c_bins
perl -lane 'print "$F[2]\t$F[0]";' < flye2_bin3c_bins.ctgassoc > flye2.bin3c.bins.tab

mkdir flye2_dastool; sbatch --nodes=1 --mem=25000 --ntasks-per-node=4 -p msn -q msn --wrap="DAS_Tool --search_engine 'diamond' -i flye2.bin3c.bins.tab -l bin3c -c flye2.contigs.fasta -o flye2_dastool/flye2.das -t 4 --write_bins 1"

# CheckM runs
for i in canu flye2; do echo $i; sbatch --nodes=1 --mem=50000 --ntasks-per-node=10 -p msn -q msn checkm lineage_wf -f $i.checkm.tab --tab_table -t 10 ${i}_bin3c_cluster/fasta/ ${i}_bin3c_checkm; done


# Summary stats on bin3c binning:
# Contigs per bin
for i in flye2_bin3c_cluster/cluster_report.csv canu_bin3c_cluster/cluster_report.csv; do echo $i; perl -e '$c = 0; $sum = 0; <>; while(<>){chomp; @s = split(/,/); $c++; $sum += $s[2];} print "$c\t$sum\t" . ($sum / $c) . "\n\n";' < $i; done
flye2_bin3c_cluster/cluster_report.csv
589     6150    10.4414261460102

canu_bin3c_cluster/cluster_report.csv
446     4765    10.6838565022422

# bases per bin
for i in flye2_bin3c_cluster/cluster_report.csv canu_bin3c_cluster/cluster_report.csv; do echo $i; perl -e '$c = 0; $sum = 0; <>; while(<>){chomp; @s = split(/,/); $c++; $sum += $s[3];} print "$c\t$sum\t" . ($sum / $c) . "\n\n";' < $i; done
flye2_bin3c_cluster/cluster_report.csv
589     1084505102      1841265.02886248

canu_bin3c_cluster/cluster_report.csv
446     946219639       2121568.69730942
``` 

<a name="checkm"></a>
#### Checkm plotting

```R
library(doMC)
library(data.table)
library(ggplot2)
library(dplyr)

flye <- read.delim("flye2.checkm.tab", header=TRUE)
canu <- read.delim("canu.checkm.tab", header=TRUE)

asms <- list(Canu = canu, metaFlye = flye)
asmnames <- names(asms)
for(i in 1:length(asmnames)){
asms[[i]]$ASM <- asmnames[i]
}

result_table <- do.call(rbind.data.frame, asms)
tmp_wide <- result_table %>% group_by(ASM) %>% filter(Contamination < 5) %>% summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90), `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70))

melt(tmp_wide,id.vars = 'ASM', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table <- melt(tmp_wide,id.vars = 'ASM', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table$title <- "CheckM Bin Statistics (< 5% Contamination)"
packageVersion("ggplot2")
plot_table$Completeness <- factor(plot_table$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%')))
colors <- rev(c("#08306B","#1664AB","#4A97C9","#93C4DE"))
pdf("flye2_canu_bin3c_checkm_stats.pdf", useDingbats=FALSE)
ggplot(plot_table, aes(ASM, Bins, fill=Completeness)) + geom_bar(stat="identity", position="stack")  + facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_manual(values = colors) +  scale_x_discrete(limits = asmnames)
dev.off()
```
<a name="offdiag"></a>
#### Bin3c off-diagonal analysis and Hi-C inter-contig links

I want to quantify the number of off-diagonal hits on the bin3c plots and to count the overall number of inter- and intra- contig links to compare against assemblies. The contact map is a pickled instance of the ContactMap class in the "mzd" folder of the bin3c repo. I think that I can pull a matrix from that and then run through the off-diagonal upper-triangle to identify bad hits.

```bash
module load samtools/1.9 miniconda
conda activate /KEEP/rumen_longread_metagenome_assembly/bin3C
```

```python
import pandas as pd
import numpy as np
import sys
sys.path.insert(1, '/project/forage_assemblies/sheep_project/bin3C')
# hmm.. no cPickle on conda. Will have to edit those documents
import pickle
import gzip

from mzd.io_utils import load_object
contact = load_object('canu_bin3c/contact_map.p.gz')

# didn't work from here -- I need to run the whole clustering pipeline I think!
AttributeError: ContactMap instance has no attribute 'seqinfo'
```

OK, try #2 -- let's add a print-out of the matrix or a threshold contact-contact table in the clustering script. I editted the plotting function in bin3C/mzd/contact_map.py

I was able to print out the matrix and labels for each cluster. The clusters take up a large proportion of cells of the matrix each, so I will have to devise labels for them. I've written a script to try to parse out all of the off-diagonal associations.

```bash
perl check_mat_offdiagonal.pl flye2_bin3c_cluster/cluster_plot.png.mat flye2_bin3c_cluster/cluster_plot.png.ticks flye2_bin3c_cluster/flye2_offdiagonal_hits.tab > flye2_bin3c_cluster/flye2_offdiagonal_stats.txt

perl check_mat_offdiagonal.pl canu_bin3c_cluster/cluster_plot.png.mat canu_bin3c_cluster/cluster_plot.png.ticks canu_bin3c_cluster/canu_offdiagonal_hits.tab > canu_bin3c_cluster/canu_offdiagonal_stats.txt

wc -l *_bin3c_cluster/*_hits.tab
  840 canu_bin3c_cluster/canu_offdiagonal_hits.tab
 1549 flye2_bin3c_cluster/flye2_offdiagonal_hits.tab
 2389 total

grep 'CL' flye2_bin3c_cluster/flye2_offdiagonal_stats.txt > flye2_bin3c_cluster/flye2_offdiag_digraph.tab
grep 'CL' canu_bin3c_cluster/canu_offdiagonal_stats.txt > canu_bin3c_cluster/canu_offdiag_digraph.tab
```

<a name="busco"></a>
#### BUSCO analysis of Eukaryotic contigs

OK, the BUSCO scores are likely to be low on the Eukaryotes, but we'll give it a try!

```bash
module load busco3 samtools/1.9

# First, I need to classify the bins and then pull the contigs from each. I've written a script to do this
perl pull_bin3c_ctg_name_tax.pl canu_bin3c_cluster/fasta/ canu_table.canu.blobtools.blobDB.table.txt canu_bin3c_bins
perl pull_bin3c_ctg_name_tax.pl flye_bin3c_cluster/fasta/ flye_table.flye.blobtools.blobDB.table.txt flye_bin3c_bins

# No unique Eukaryotic bins... interesting. I guess that I'll have to do this after DAS_tool dereplication?
# In the meantime, I'll pull each fasta individually and run BUSCO on them
perl -ne 'if($_ =~ /^#/){next;} chomp; @F = split(/\t/); if($F[5] eq "Eukaryota"){print "$F[0]\n";}' < canu_table.canu.blobtools.blobDB.table.txt > canu_table.canu.blobtools.euk_ctg.list
perl -ne 'if($_ =~ /^#/){next;} chomp; @F = split(/\t/); if($F[5] eq "Eukaryota"){print "$F[0]\n";}' < flye_table.flye.blobtools.blobDB.table.txt > flye_table.flye.blobtools.euk_ctg.list

for i in canu flye; do mkdir ${i}_euk; for j in `cat ${i}_table.$i.blobtools.euk_ctg.list`; do echo $j; samtools faidx $i.contigs.fasta $j > ${i}_euk/$j.ctg.fa; done; done
cat flye_euk/*.fa > flye_euk/full_flye_euk.ctg.fa
cat canu_euk/*.fa > canu_euk/full_canu_euk.ctg.fa

# OK, there are complex rules for working with BUSCO that are uncovered by using the "module show busco3" command
cp /software/7/apps/busco3/3.1.0/config/config.ini.default ./busco_config.ini

export BUSCO_CONFIG_FILE=/project/forage_assemblies/sheep_project/busco_config.ini
cp -Rp /software/apps/augustus/gcc/64/3.2.3/config /project/rumen_longread_metagenome_assembly/AUGUSTUS_CONFIG
export AUGUSTUS_CONFIG_PATH=/project/rumen_longread_metagenome_assembly/AUGUSTUS_CONFIG

# Testing
run_BUSCO.py -i flye_euk/contig_1003.ctg.fa -l /reference/data/BUSCO/latest/nematoda_odb9 -o busco.contig_1003 -m geno

# Queuing up the rest of the scripts
mkdir flye_busco
cd flye_busco/
for i in /project/forage_assemblies/sheep_project/flye_euk/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap="run_BUSCO.py -i $i -l /reference/data/BUSCO/latest/nematoda_odb9 -o busco.${name} -m geno"; done
cd ..

mkdir canu_busco
cd canu_busco/
for i in /project/forage_assemblies/sheep_project/canu_euk/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap="run_BUSCO.py -i $i -l /reference/data/BUSCO/latest/nematoda_odb9 -o busco.${name} -m geno"; done
cd ..

# The nematode database had few hits. Let's try the eukaryota database?
cd flye_busco/
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap="run_BUSCO.py -i /project/forage_assemblies/sheep_project/flye_euk/full_flye_euk.ctg.fa -l /reference/data/BUSCO/latest/eukaryota_odb9 -o busco.eukaryota.full -m geno"
cd ..

cd canu_busco/
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p msn -q msn --wrap="run_BUSCO.py -i /project/forage_assemblies/sheep_project/canu_euk/full_canu_euk.ctg.fa -l /reference/data/BUSCO/latest/eukaryota_odb9 -o busco.eukaryota.full -m geno"
cd ..
```

<a name="desman"></a>
#### DESMAN strain inference

NOTE: I've installed DESMAN in my Blobtools virtual environment because they require the same version of python and they're used in pretty similar circumstances.

> Ceres:

```bash
module load r/3.6.1 prodigalorffinder/2.6.3 samtools/1.9 miniconda
conda activate /KEEP/rumen_longread_metagenome_assembly/blobtools/

mkdir desman
for i in flye canu; do mkdir desman/${i}_bins; mkdir desman/${i}_beds; done

# Creating the lists from the bin3c bins
for i in flye canu; do echo $i; perl -e 'chomp(@ARGV); $outfold = "$ARGV[1]_bins"; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); open(OUT, ">> desman/$outfold/$F[1].list"); print OUT "$F[0]\n"; close OUT;} close IN;' ${i}_dastool/${i}.das_DASTool_scaffolds2bin.txt $i; done

# And the bed lists
for i in flye canu; do echo $i; cat ${i}_dastool/*.scg > ${i}_dastool/$i.combined.scg; done

for i in canu flye; do grep '>' ${i}_dastool/${i}.das_proteins.faa |  perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > ${i}_dastool/${i}.das_proteins.shortform.tab; done

for i in canu flye; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ${i}_dastool/${i}.das_proteins.shortform.tab -c 0 -l ${i}_dastool/$i.combined.scg | perl -lane '$r = $F[0]; $r =~ s/_\d{1,3}$//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > ${i}_dastool/$i.prod.proteins.scg.loc.bed; done

for i in flye canu; do echo $i; for j in desman/${i}_bins/*.list; do name=`basename $j | cut -d'.' -f1,2`; echo $name; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ${i}_dastool/$i.prod.proteins.scg.loc.bed -c 0 -l $j > desman/${i}_beds/${name}.scg.bed; done; done

### Testing it out
WORKING=/project/forage_assemblies/sheep_project/
sbatch -p msn -q msn ~/python_toolchain/metagenomics/desmanStrainInference.py -a $WORKING/flye.contigs.fasta -c $WORKING/desman/flye_bins/bin.611.list -g $WORKING/desman/flye_beds/bin.611.scg.bed -d /project/rumen_longread_metagenome_assembly/binaries/DESMAN -b $WORKING/flye_wgs/63/63.sorted.merged.bam -o $WORKING/desman/test

### Now queuing up the real thing
for j in canu flye; do echo $j; for i in desman/${j}_bins/*.list; do echo $WORKING/$i; name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -p msn -q msn ~/python_toolchain/metagenomics/desmanStrainInference.py -a $WORKING/$j.contigs.fasta -c $WORKING/$i -g $WORKING/desman/${j}_beds/$name.scg.bed -d /project/rumen_longread_metagenome_assembly/binaries/DESMAN -b $WORKING/${j}_wgs/63/63.sorted.merged.bam -o $WORKING/desman/${j}_output/$name; done; done

### Testing desman strain prediction
sbatch -p msn -q msn ~/python_toolchain/metagenomics/desmanStrainPrediction.py -a /project/forage_assemblies/sheep_project/canu.contigs.fasta -c /project/forage_assemblies/sheep_project/canu_dastool/canu.das_DASTool_bins/bin.1159.contigs.fa -o /project/forage_assemblies/sheep_project/desman/canu_output/bin.1159 -d /project/rumen_longread_metagenome_assembly/binaries/DESMAN -g /project/forage_assemblies/sheep_project/desman/canu_beds/bin.1159.scg.bed

### Running on version2
mkdir desman/flye2_beds; mkdir desman/flye2_bins

perl -e 'chomp(@ARGV); $outfold = "$ARGV[1]_bins"; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); open(OUT, ">> desman/$outfold/$F[1].list"); print OUT "$F[0]\n"; close OUT;} close IN;' flye2_dastool/flye2.das_DASTool_scaffolds2bin.txt flye2

cat flye2_dastool/*.scg > flye2_dastool/flye2.combined.scg

i=flye2; grep '>' ${i}_dastool/${i}.das_proteins.faa |  perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > ${i}_dastool/${i}.das_proteins.shortform.tab

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ${i}_dastool/${i}.das_proteins.shortform.tab -c 0 -l ${i}_dastool/$i.combined.scg | perl -lane '$r = $F[0]; $r =~ s/_\d{1,4}$//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > ${i}_dastool/$i.prod.proteins.scg.loc.bed

for i in flye2; do echo $i; for j in desman/${i}_bins/*.list; do name=`basename $j | cut -d'.' -f1,2`; echo $name; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ${i}_dastool/$i.prod.proteins.scg.loc.bed -c 0 -l $j > desman/${i}_beds/${name}.scg.bed; done; done

WORKING=/project/forage_assemblies/sheep_project/
for j in flye2; do echo $j; for i in desman/${j}_bins/*.list; do echo $WORKING/$i; name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -p msn -q msn ~/python_toolchain/metagenomics/desmanStrainInference.py -a $WORKING/$j.contigs.fasta -c $WORKING/$i -g $WORKING/desman/${j}_beds/$name.scg.bed -d /project/rumen_longread_metagenome_assembly/binaries/DESMAN -b $WORKING/${j}_wgs/63/63.sorted.merged.bam -o $WORKING/desman/${j}_output/$name; done; done
```


<a name="generating"></a>
## Generating figures and tables

It's time to start summarizing the stats from this run.

<a name="kingdomlevel"></a>
#### Kingdom-level and GC bin contig size plot

First I need to prepare the data for the plots.

```bash
# Preparing ctg length, GC and superkingdom tables
for i in canu flye; do echo $i; perl -ne 'chomp; @F = split(/\t/); if($F[0] =~ /^#/){next;} print "$F[1]\t$F[2]\t$F[5]\n";' < ${i}_table.$i.blobtools.blobDB.table.txt > $i.blobtools.lenbygcbyking.tab; done

perl -ne 'chomp; @F = split(/\t/); if($F[0] =~ /^#/){next;} print "$F[1]\t$F[2]\t$F[5]\n";' < flye2_table.flye2.blobtools.blobDB.table.txt > flye2.blobtools.lenbygcbyking.tab
```

```R
library(ggridges)
library(ggplot2)
library(dplyr)

canu <- read.delim("canu.blobtools.lenbygcbyking.tab", header=FALSE)
flye <- read.delim("flye.blobtools.lenbygcbyking.tab", header=FALSE)

canu <- mutate(canu, ASM=c("canu"))
flye <- mutate(flye, ASM=c("flye"))
colnames(canu) <- c("LEN", "GC", "KING", "ASM")
colnames(flye) <- c("LEN", "GC", "KING", "ASM")

combined <- bind_rows(canu, flye)
combined$ASM <- as.factor(combined$ASM)
combined$KING <- factor(combined$KING, levels = c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))

pdf(file="superkingdom_length.pdf", useDingbats=FALSE)
ggplot(combined, aes(y=KING, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(1000, 10000, 100000, 1000000, 6000000), limits=c(100, 6000000), labels=c("1000", "10,000", "100,000", "1,000,000", "6,000,000")) + xlab(label = "Log10 Contig Lengths (bp)") + ylab(label= "Contig Superkingdom Taxonomic Assignment")
dev.off()

# now for the GC%
pdf(file="gc_assembly_counts.pdf", useDingbats=FALSE)
ggplot(combined, aes(x=ASM, y=GC, fill=ASM)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white") + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + labs(x="Assembly", y = "Avg GC ratio per contig")
dev.off()

### Flye v2 plots
canu <- read.delim("canu.blobtools.lenbygcbyking.tab", header=FALSE)
flye <- read.delim("flye2.blobtools.lenbygcbyking.tab", header=FALSE)

canu <- mutate(canu, ASM=c("canu"))
flye <- mutate(flye, ASM=c("flye"))
colnames(canu) <- c("LEN", "GC", "KING", "ASM")
colnames(flye) <- c("LEN", "GC", "KING", "ASM")

combined <- bind_rows(canu, flye)
combined$ASM <- as.factor(combined$ASM)
combined$KING <- factor(combined$KING, levels = c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))

pdf(file="superkingdom_flye2_length.pdf", useDingbats=FALSE)
ggplot(combined, aes(y=KING, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(1000, 10000, 100000, 1000000, 6000000), limits=c(100, 6000000), labels=c("1000", "10,000", "100,000", "1,000,000", "6,000,000")) + xlab(label = "Log10 Contig Lengths (bp)") + ylab(label= "Contig Superkingdom Taxonomic Assignment")
dev.off()

pdf(file="gc_assembly_counts_flye2.pdf", useDingbats=FALSE)
ggplot(combined, aes(x=ASM, y=GC, fill=ASM)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white") + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + labs(x="Assembly", y = "Avg GC ratio per contig")
dev.off()
```

<a name="orfprediction"></a>
#### ORF prediction classification

I think that I have everything that I need for this in my *.shortform.tab files. Let's check.

```R
library(ggplot2)
library(dplyr)
library(ggridges)

canu <- read.delim("canu.contigs.prod.shortform.tab")
flye <- read.delim("flye.contigs.prod.shortform.tab")

cdata <- mutate(canu, ASM=c("canu"), LEN=End - Start) %>% select(ASM, LEN)
fdata <- mutate(flye, ASM=c("flye"), LEN=End - Start) %>% select(ASM, LEN)

combined <- bind_rows(cdata, fdata)
pdf(file="orf_lengths.pdf", useDingbats=FALSE)
ggplot(combined, aes(y=ASM, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(100, 1000, 10000, 50000), limits=c(50, 55000), labels=c("100", "1000", "10,000", "50,000")) + xlab(label = "Log10 ORF lengths (bp)") + ylab(label="Assembly")
```

# TODO:

## Redo plot of contig sizes per superkingdom
## Prepare table of blobplots summary stats (contig N50, total size, etc)
## Prepare bin3c table of bins
## Run CheckM on canu and Bin3C bins
## Finish writeup of methods
## Quantify chimeric long-reads
