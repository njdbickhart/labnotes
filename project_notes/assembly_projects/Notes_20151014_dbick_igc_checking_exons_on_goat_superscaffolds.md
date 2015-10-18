# Identifying IGC regions that are resolved from optical mapping
---
*10/14/2015*

These are my notes on doing a very rapid identification of immune gene cluster (IGC) exons that span different PacBio contigs, but are linked by optical mapping technologies -- namely BioNano Genomics Irys!


## Preparing materials for the analysis

I'm going to have to get exons of the genes that I'm interested in interrogating. The first step is getting Human and Cattle exon sequence. Human is likely better annotated, but Cattle will be the closest/highest quality ruminant reference genome sequence to Goat. I will then map both sets to the Goat reference assembly and try to identify exons that span different contigs. Then I will see if any of those exons lie within Super-scaffolds.

First, let's get the fasta files for the list of contigs from Ensembl.

OK, ran into a problem: to restrict bandwidth, Ensembl restricts gene accession to EnsGene and EnsTranscript IDs. Here is the XML of how I retrieved the Human genes with GO type "Immune response:"

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "transcript_source" value = "ensembl"/>
		<Filter name = "go_parent_name" value = "immune response"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>
```

There were 400+ genes with that GO annotation. Now for Cattle.

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "btaurus_gene_ensembl" interface = "default" >
		<Filter name = "transcript_source" value = "ensembl"/>
		<Filter name = "go_parent_name" value = "immune response"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>
```

Strange, cattle had 600 genes this time. Oh well, I'm going to try pulling gene lists from the UCSC.

The UCSC genome browser does have a nice list of exon coordinates for the assembly, but it does not contain the exon fastas themselves. I'll have to generate that myself with a custom script. I wrote a script [here](https://github.com/njdbickhart/perl_toolchain/blob/master/bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl) that should take care of the problem.

> pwd: /home/dbickhart/share/goat_assembly_paper/igc_annotation

```bash
perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\n";}' < cattle_ensid_immune_response.txt > cattle_enstranscript_immune_response.txt

perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl -b umd3EnsGene.txt -f ../../umd3_data/umd3_kary_unmask_ngap.fa -e cattle_enstranscript_immune_response.txt -o cattle_enstranscript_immune_response.fa
```

I'm really curious to see how this turns out. Let's start the alignment and see how the data piles up. I'm going to use the frozen BNG version 5 assembly in order to test this out -- I should have the super scaffolds nicely broken up and tagged already.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/annotation

```bash
bwa mem ../papadum-v5bng-pilon-split2kbgap.fa cattle_enstranscript_immune_response.fa | samtools view -bS - | samtools sort -T cattle.temp -o cattle_igc_goat_align.bam -

samtools index cattle_igc_goat_align.bam
samtools flagstat cattle_igc_goat_align.bam
	5275 + 0 in total (QC-passed reads + QC-failed reads)
	0 + 0 secondary
	35 + 0 supplementary
	0 + 0 duplicates
	5210 + 0 mapped (98.77%:-nan%)
	0 + 0 paired in sequencing
	0 + 0 read1
	0 + 0 read2
	0 + 0 properly paired (-nan%:-nan%)
	0 + 0 with itself and mate mapped
	0 + 0 singletons (-nan%:-nan%)
	0 + 0 with mate mapped to a different chr
	0 + 0 with mate mapped to a different chr (mapQ>=5)

# So there are 35 that have multiple mappings
# Lots have aligned to the super-reads -- let's see if we can find alternative mappings

samtools view cattle_igc_goat_align.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); $s[0] =~ s/\_\d{1,3}$//; $h{$s[0]}->{$s[2]} = 1;} foreach my $contig (keys %h){foreach my $chr (keys %{$h{$contig}}){ print "$contig\t$chr\n";}}' | less
ENSBTAT00000007438      Scaffold_154.1
ENSBTAT00000007438      Scaffold_94.1

ENSBTAT00000007743      Scaffold_94.1
ENSBTAT00000007743      Scaffold_134.1

# Very good example Interleukin 1 alpha
ENSBTAT00000013665      Scaffold_22.3
ENSBTAT00000013665      Scaffold_22.1
```

Let's put this into a formalized script and generate a summary file. I've written a script that goes along with this pipeline [here](https://github.com/njdbickhart/perl_toolchain/blob/master/assembly_scripts/identifyCrossContigExonAlignments.pl).

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/annotation

```bash
perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b cattle_igc_goat_align.bam -o cattle_igc_goat_align.summary

ENSBTAT00000013665      Scaffold_22.1   7       Scaffold_22.3   2,4,5,6,7,7
```

Let's do this with the human annotations now. Unfortunately, the UCSC genome browser doesn't handle all annotations in a uniform format. I'm going to have to convert their file by changing the first column to the actual ensembl gene IDs.

> pwd: /home/dbickhart/share/goat_assembly_paper/igc_annotation

```bash
# Converting with the UCSC conversion database text file
perl -e '%c; chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\.\d{1,3}//; $c{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($c{$s[0]})){$s[0] = $c{$s[0]};} print "1\t" . join("\t", @s); print "\n";} close IN;' humanknownToEnsembl.txt humanknownGene.txt > human_ensgeneconverted.txt

# Now I can run my script!
perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\n";}' < hg_ensid_immune_response.txt > hg_enstrans_immune_response.txt

perl ../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl -b human_ensgeneconverted.txt -f hg38.fa -e hg_enstrans_immune_response.txt -o human_ensgene_immune.fa
```

Now to transfer everything to the server so that I can do the alignment and subsequent summary.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/annotation

```bash
# Aligning the exons
bwa mem ../papadum-v5bng-pilon-split2kbgap.fa human_ensgene_immune.fa | samtools view -bS - | samtools sort -T human.temp -o human_igc_goat_align.bam -

# Now to summarize them.
perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b human_igc_goat_align.bam -o human_igc_goat_align.summary


# Quite allot of exons were completely unmapped. It would be interesting to see how many were unmapped in our assembly
# Vs the BGI assembly.
bwa index chi_ref_CHIR_1.0_all.fa
bwa mem chi_ref_CHIR_1.0_all.fa human_ensgene_immune.fa | samtools view -bS - | samtools sort -T human.temp -o human_igc_goat_chi_align.bam -
perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b human_igc_goat_chi_align.bam -o human_igc_goat_chi_align.summary

bwa mem chi_ref_CHIR_1.0_all.fa cattle_enstranscript_immune_response.fa | samtools view -bS - | samtools sort -T cattle.temp -o cattle_igc_goat_chi_align.bam -
perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b cattle_igc_goat_chi_align.bam -o cattle_igc_goat_chi_align.summary

```

From a brief scan of the chinese and pacbio assemblies, I can see that there are discrepencies in the number of scaffolds that contain the IGCs. Let's count them up.

```bash
# First, let's get a listing of the number of unique scaffold mappings (by exons) and unmapped regions
for i in *.summary; do echo $i; perl -lane '%h; %ccount; for($x = 1; $x < scalar(@F); $x += 2){$contig = $F[$x]; $exons = $F[$x + 1]; @s = split(/,/, $exons); if($contig ne "unmapped"){$h{"mapped"} += scalar(@s); $ccount{$contig} = 1;}else{$h{"unmapped"} += scalar(@s);}} if(!exists($h{"mapped"})){ $h{"mapped"} = 0;} if(!exists($h{"unmapped"})){$h{"unmapped"} = 0;} print "$F[0]\t" . (scalar(keys(%ccount))) . "\t$h{mapped}\t$h{unmapped}"; %h = (); %ccount = ();' < $i > $i.stats; done
```

Here are the columns of the stats files:

* Ensembl ID
* Number of mappable contigs that contain the ID
* Number of exons on mappable contigs
* Number of exons that were unmapped

```bash
# Summarizing the number of scaffolds that contain the exons
cat cattle_igc_goat_align.summary.stats | cut -f2 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
	total   655
	Minimum 0
	Maximum 7
	Average 1.129771
	Median  1
	Standard Deviation      0.521656
	Mode(Highest Distributed Value) 1

cat cattle_igc_goat_chi_align.summary.stats | cut -f2 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
	total   655
	Minimum 1
	Maximum 15
	Average 1.325191
	Median  1
	Standard Deviation      1.070693
	Mode(Highest Distributed Value) 1

# Interesting, so the BGI assembly has more exon distribution, but the PacBio assembly actually has completely unmapped genes.
# Let's see what gene that was for fairness sake
cat cattle_igc_goat_align.summary.stats | perl -lane 'if($F[1] == 0){print $_;}'
	ENSBTAT00000066072      0       0       1
# Checking the BGI assembly
grep 'ENSBTAT00000066072' cattle_igc_goat_chi_align.summary.stats
	ENSBTAT00000066072      1       1       0

# So one exon out of everything. Let's see the gene alignment
samtools view cattle_igc_goat_chi_align.bam | grep 'ENSBTAT00000066072'
ENSBTAT00000066072_1    0       gi|541128976|ref|NC_022303.1|   103748862       60      710S176M ...

# Pretty heavily soft-clipped, but let's count this as a victory for BGI. 
# Let's make a count of transcripts
for i in *.stats; do echo -n -e "$i\t"; perl -e '$c = 0; $t = 0; while(<>){chomp; @s = split(/\t/); if($s[1] == 1){$c++;} $t++;} print "$c\t$t\n";' < $i; done

# OOPs! I need to align to the "pseudo chromosome" assembly as well!
bwa index chi_ref_CHIR_1.0_chrall.fa
bwa mem chi_ref_CHIR_1.0_chrall.fa human_ensgene_immune.fa | samtools view -bS - | samtools sort -T human.temp -o human_igc_goat_chichr_align.bam -
bwa mem chi_ref_CHIR_1.0_chrall.fa cattle_enstranscript_immune_response.fa | samtools view -bS - | samtools sort -T cattle.temp -o cattle_igc_goat_chichr_align.bam -

perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b human_igc_goat_chichr_align.bam -o human_igc_goat_chichr_align.summary
perl ~/perl_toolchain/assembly_scripts/identifyCrossContigExonAlignments.pl -b cattle_igc_goat_chichr_align.bam -o cattle_igc_goat_chichr_align.summary

# Let's make a count of transcripts
for i in *.stats; do echo -n -e "$i\t"; perl -e '$c = 0; $t = 0; while(<>){chomp; @s = split(/\t/); if($s[1] == 1){$c++;} $t++;} print "$c\t$t\n";' < $i; done
```
#### Ensembl transcripts that are placed on only one scaffold per assembly

| Stat file | Assembly | transcripts on one scaffold | total transcripts| ratio placed |
| :--- | :--- | ---: | ---: | ---: |
cattle_igc_goat_align.summary.stats  | PacBio |   601  |   655 | 0.917
cattle_igc_goat_chi_align.summary.stats | CHI_1.0 |  544  |   655 |  0.831
cattle_igc_goat_chichr_align.summary.stats | CHI_1.0  |   606  |   655 | 0.925
human_igc_goat_align.summary.stats  | PacBio |    475  |   581 |  0.818
human_igc_goat_chi_align.summary.stats | CHI_1.0 |  394  |   581 |  0.678
human_igc_goat_chichr_align.summary.stats | CHI_1.0 |   469  |   581 |  0.807

So, on a per-scaffold basis, PacBio is much better than CHI_1.0. Only the Pseudochromosomes save the data here. Now let's tabulate the number of super-scaffolds that have transcript mappings

```bash
perl -e '%h = ("scaffold" => 0, "contig" => 0); $t = 0; while(<>){chomp; @s = split(/\t/); $scaffold = 0; $contig = 0; $t++; for($i = 1; $i < scalar(@s); $i += 2){if($s[$i] =~ /Scaffold/){$scaffold = 1;}elsif($s[$i] ne "unmapped"){$contig = 1;}} if($scaffold){$h{"scaffold"} +=1;}elsif($contig){$h{"contig"} += 1;}} print "$h{scaffold}\t$h{contig}\t$t\n";' < cattle_igc_goat_align.summary
perl -e '%h = ("scaffold" => 0, "contig" => 0); $t = 0; while(<>){chomp; @s = split(/\t/); $scaffold = 0; $contig = 0; $t++; for($i = 1; $i < scalar(@s); $i += 2){if($s[$i] =~ /Scaffold/){$scaffold = 1;}elsif($s[$i] ne "unmapped"){$contig = 1;}} if($scaffold){$h{"scaffold"} +=1;}elsif($contig){$h{"contig"} += 1;}} print "$h{scaffold}\t$h{contig}\t$t\n";' < human_igc_goat_align.summary
```

| dataset | scaffold mappings | contig mappings | total | ratio scaffold |
| :--- | ---: | ---: | ---: | ---: |
cattle ensid | 624  |   30  |    655 | 0.952
human ensid | 523  |   8    |   581  | 0.900

Finally, I'm going to create bed files of the alignments so that I can create R plots using [Sushi](https://www.bioconductor.org/packages/release/bioc/html/Sushi.html). 

```bash
# Creating the bed files with a one-liner
samtools view cattle_igc_goat_align.bam | perl -lane 'if($F[2] =~ /\*/){next;} $orient = ($F[1] & 0x10)? -1 : 1; $F[0] =~ s/\.\d{1,3}//; $end = length($F[9]); print "$F[2]\t$F[3]\t$end\t$F[0]\t\.\t$orient\texon";' > cattle_igc_goat_align.bed
samtools view cattle_igc_goat_chi_align.bam | perl -lane 'if($F[2] =~ /\*/){next;} $orient = ($F[1] & 0x10)? -1 : 1; $F[0] =~ s/\.\d{1,3}//; $end = length($F[9]); print "$F[2]\t$F[3]\t$end\t$F[0]\t\.\t$orient\texon";' > cattle_igc_goat_chi_align.bed

samtools view human_igc_goat_align.bam | perl -lane 'if($F[2] =~ /\*/){next;} $orient = ($F[1] & 0x10)? -1 : 1; $F[0] =~ s/\.\d{1,3}//; $end = length($F[9]); print "$F[2]\t$F[3]\t$end\t$F[0]\t\.\t$orient\texon";' > human_igc_goat_align.bed
samtools view human_igc_goat_chi_align.bam | perl -lane 'if($F[2] =~ /\*/){next;} $orient = ($F[1] & 0x10)? -1 : 1; $F[0] =~ s/\.\d{1,3}//; $end = length($F[9]); print "$F[2]\t$F[3]\t$end\t$F[0]\t\.\t$orient\texon";' > human_igc_goat_chi_align.bed
```

#### Target genes
| Ens transcript ID | Gene name | Description | Comments |
| :--- | :--- | :--- | :--- |
ENSBTAT00000007129 | NLRC5 | Involved in Cytokine response | Full contig in all assemblies; PacBio has more mapped exons 
ENSBTAT00000022979 | C3-201 | Complement component 3 | Nearly full contig in PacBio; 14 scaffolds contain exons in BGI
ENSBTAT00000011795 | MHCI JSP.1 | MHC class I receptor | Spread across 4 scaffolds BGI

OK, time to draw the plots in R using Sushi.

```R
library(Sushi)
base.cattle.pacbio <- read.table("cattle_igc_goat_align.bed", header=FALSE)
colnames(base.cattle.pacbio) <- c("chrom","start","stop","gene","score","strand","type")

base.cattle.bgi <- read.table("cattle_igc_goat_chi_align.bed", header=FALSE)
colnames(base.cattle.bgi) <- c("chrom","start","stop","gene","score","strand","type")

# Starting with the best example first: the complement component 3 gene
ENSBTAT00000022979.cattle.pacbio <- base.cattle.pacbio[base.cattle.pacbio$gene == "ENSBTAT00000022979", ]
ENSBTAT00000022979.cattle.bgi <- base.cattle.bgi[base.cattle.bgi$gene == "ENSBTAT00000022979", ]

# Important note: chrom, gene and score all need to be character types!
ENSBTAT00000022979.cattle.pacbio$chrom <- as.character(ENSBTAT00000022979.cattle.pacbio$chrom)
ENSBTAT00000022979.cattle.pacbio$gene <- as.character(ENSBTAT00000022979.cattle.pacbio$gene)
ENSBTAT00000022979.cattle.pacbio$score <- as.character(ENSBTAT00000022979.cattle.pacbio$score)

pg <- plotGenes(ENSBTAT00000022979.cattle.pacbio, "Scaffold_330.3", 1170000, 1209000, plotgenetype="arrow", bentline=FALSE,fontsize=1.2, arrowlength = 0.025, labeltext=TRUE)
labelgenome("Scaffold_330.3", 1170000, 1209000, n=3, scale="Kb")
dev.copy2pdf(file="ENSBTAT00000022979.cattle.pacbio.pdf", useDingbats=FALSE)

# The BGI copy is harder to plot -- I'm going to try to use par to try to plot multiple chromosomes
ENSBTAT00000022979.cattle.bgi$chrom <- as.character(ENSBTAT00000022979.cattle.bgi$chrom)
ENSBTAT00000022979.cattle.bgi$gene <- as.character(ENSBTAT00000022979.cattle.bgi$gene)

# Well, the figure margins are too large, so I'll have to print them out separately and concatenate later
pg <- plotGenes(ENSBTAT00000022979.cattle.bgi, "gi|541128980|ref|NC_022299.1|", 14272000, 14306000, plotgenetype="arrow", bentline=FALSE,fontsize=1.2, arrowlength = 0.025, labeltext=TRUE)
labelgenome("gi|541128980|ref|NC_022299.1|", 14272000, 14306000, n=3, scale="Kb")
dev.copy2pdf(file="ENSBTAT00000022979.cattle.bgi.1.pdf", useDingbats=FALSE)
```

I think that the NCBI pipe ("|") characters are screwing things up here. I'll go with the nuclear option and draw a [parasight](http://eichlerlab.gs.washington.edu/jeff/parasight/) image instead.

```bash
# Preparing the extra and showseq file for the complement component gene
echo -e "seqname\tlength\tbegin\tend\nScaffold_330.3\t39000\t1170000\t1209000\nutg3912\t5000\t1700\t2200" > ENSBTAT00000022979.pacbio.seq
grep 'ENSBTAT00000022979' cattle_igc_goat_align.bed | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\tblue\t-10";' > ENSBTAT00000022979.pacbio.extra
