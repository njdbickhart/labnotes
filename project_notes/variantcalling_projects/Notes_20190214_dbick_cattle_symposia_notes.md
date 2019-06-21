# ADSA Cattle Symposia Talk Notes
---
*2/14/2018*

These are my notes and commands for generating datasets for use in my talk at the 2019 ASDA Cattle Genetics symposia.


## Table of Contents

## Lit review and structure

### Genetic and physical mapping

[James Womack](https://genome.cshlp.org/content/15/12/1699.full.html) has a great review on the progress up until the early 00's. The first gene map was by [Heuertz et al](https://www.karger.com/Article/Abstract/131601) but [Womack](https://www.ncbi.nlm.nih.gov/pubmed/3082971) expanded this map greatly by using cattle-hamster hybrid somatic cell lines in 1986. In-situ hybridization was added to cattle-hamster hybrid somatic cell lines to make a comprehensive map of [many cattle genes](https://link.springer.com/article/10.1007/BF00296815) in the early 90's. A [multi-national group](https://www.nature.com/articles/ng0394-227) and [USDA-MARC](https://www.genetics.org/content/136/2/619.abstract?ijkey=c8173747b99141827e87e12dd6281489f54c1196&keytype2=tf_ipsecsha) published a linkage map almost simultaneously in 1994 spanned 2464 centimorgans of the cattle genome. **NOTE: the USDA-MARC map could be fun to plot on ARS-UCDv1.2.** The USMARC map used 280 out of 468 microsatellites to develop linkage maps for cattle. Next, Womack led the development of a [radiation-hybrid map](https://link.springer.com/article/10.1007%2Fs003359900593?LI=true) for cattle. 

James Womack predicted that the bovine genome [would never be sequenced](https://www.annualreviews.org/doi/pdf/10.1146/annurev-animal-020518-114902) in the 90's, in order to justify grant proposals for comparative mapping. 

### Molecular basis of testing

The earliest markers used by cattle geneticists were probed southern blots and PCR-based markers. Term "satellite" comes from [Saul Kit](https://www.sciencedirect.com/science/article/pii/S0022283661800752?via%3Dihub) who found bands from CsCl gradient centrifugation. Microsatellites are 1 to 10 bp repetitive, low complexity sequence. Useful because they are multi-allelic and are frequently found outside of coding regions (no positive selection). 

### Reduced representation sequencing for SNP50 (Van Tassel Nat. Methods)

Curt, Tad and company first used DraI (TTT^AAA) in an *in silico* digest of Btau3.1 but found that it did not perform as well *in vitro*. This may be due to their use of the incomplete assembly (Btau3.1 = 2.4 Gbp assembly size; 

### SNP chips

The [Affy 10K chip](https://www.affymetrix.com/support/technical/byproduct.affx?product=bo-10ksnp) was released in 2005 during the cattle genome sequencing project. 92% of markers were identified from the sequence data of four cattle breeds, Holstein, Angus, Limousin and Hereford) and 8% were contributed by CSIRO were identified as being within genes (though not necessarily coding. The [Bovine SNP50](https://www.illumina.com/Documents/products/datasheets/datasheet_bovine_snp5O.pdf) was released in ~2009 and improved on the 10k chip. It includes some of the 10k markers, but also has 23k markers from Curt's RRS libraries and ~10k markers from the [Bovine Hapmap project](https://science.sciencemag.org/content/324/5926/528.long). The [Bovine HD](https://www.illumina.com/Documents/products/datasheets/datasheet_bovineHD.pdf) has 93% novel content from a separate 180X coverage WGS alignment of 20 breeds of cattle against the UMD3.1 reference genome

### Genome assemblies

All [Baylor-led](https://www.hgsc.bcm.edu/other-mammals/bovine-genome-project) genome assembly projects are hosted online and have links to the NCBI repositories to download the data. 

Let's download them all for further analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/symposium_comparison

```bash
wget ftp://ftp.hgsc.bcm.edu/Btaurus/fasta/Btau20040927-freeze/contigs/Btau20040927-freeze-contigs.fa.gz
gunzip Btau20040927-freeze-contigs.fa.gz

# There are 4 iterations of Btau4.0 that were released on a yearly basis, and incorporated new polishing information. I'm going to pick the "best" one that was hosted on NCBI for this and future comparisons
# I need to remove the sill NCBI annotations in the fasta header first
mkdir btau4_6
cd btau4_6/
for i in `seq 1 29` X Y; do file="chr${i}.fa.gz"; name="chr${i}"; perl -e 'chomp(@ARGV); open(IN, "gunzip -c $ARGV[0] |"); <IN>; print ">$ARGV[1]\n"; while(<IN>){print "$_";} close IN;' $file $name; done > btau_4.6.karytype.fa
cd ..
mv btau4_6/btau_4.6.karytype.fa ./

# UMD3
cp ../repeatmasker/umd3_reference_genome.fasta ./

# And ARS-UCDv1.2
cp ../ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa ./
```

L1 Dominette is the reference cow, and her sire is L1 Domino. Both are line 1 Hereford cattle, which is an inbred line of Beef cattle developed in 1934. The goal was to replicate the successes of inbred plant line development to fix positive traits and to enable heterosis by inbred line crossing. 

## Assessing relative genome resolution of different technologies

I want to draw comparative maps of the different technologies used to track genetic variants in the cattle genome. I intend this to visualize the progress that each method has had on the resulting assembly map and to show the differences in resolution that each method provides.

## Comparing improvements of genome assemblies throughout the years

I see the following figures as being integral to the presentation and subsequent manuscript:

1. A Circos plot showing marker locations on the assembly from physical maps
2. A histogram showing the assembly size, and continuity of the assemblies, from Btau1.0 to ARS-UCDv1.2
3. A Chord diagram showing the differences between Btau1.0 to ARS-UCDv1.2

Let's start with minimap alignments of each assembly so that we can generate chord diagrams.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/symposium_comparison

```bash
module load minimap2/2.6 samtools bedtools/2.25.0

# Getting the chromosome/contig sizes
for i in ARS-UCD1.2_Btau5.0.1Y.fa Btau20040927-freeze-contigs.fa btau_4.6.karytype.fa; do echo $i; samtools faidx $i; done

# Getting the UMD3 gap sizes
java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f umd3_reference_genome.fasta -o umd3_reference_genome.gaps.bed -s umd3_reference_genome.gaps.stats

# ARS-UCD vs Btau1.0
sbatch --nodes=1 --mem=10000 --ntasks-per-node=3 -p short --wrap="minimap2 -x asm5 ARS-UCD1.2_Btau5.0.1Y.fa Btau20040927-freeze-contigs.fa > btau_1_on_ars-ucd.paf"
# ARS-UCD vs Btau4.6.1
sbatch --nodes=1 --mem=10000 --ntasks-per-node=3 -p short --wrap="minimap2 -x asm5 ARS-UCD1.2_Btau5.0.1Y.fa btau_4.6.karytype.fa > btau_4.6_on_ars-ucd.paf"
# ARS-UCD vs UMD3.1
sbatch --nodes=1 --mem=10000 --ntasks-per-node=3 -p short --wrap="minimap2 -x asm5 ARS-UCD1.2_Btau5.0.1Y.fa umd3_reference_genome.fasta > umd_3_on_ars-ucd.paf"

# Converting PAFs into sorted, merged beds so that I can calculate the alignment stats
for i in *.paf; do name=`echo $i | cut -d'_' -f1,2`; echo $name; perl -lane 'print "$F[0]\t$F[2]\t$F[3]";' < $i | bedtools sort -i stdin | bedtools merge -i stdin > ${name}.${name}.mapped.bed; perl -lane 'print "$F[5]\t$F[7]\t$F[8]";' < $i | bedtools sort -i stdin | bedtools merge -i stdin > ${name}.ARSUCD.mapped.bed; done


wc -l *.bed
  704318 btau_1.ARSUCD.mapped.bed
  811541 btau_1.btau_1.mapped.bed
    9725 btau_4.6.ARSUCD.mapped.bed
   16798 btau_4.6.btau_4.6.mapped.bed
   23692 umd_3.ARSUCD.mapped.bed
   26636 umd_3.umd_3.mapped.bed
 1592710 total

# Now to count the tallies of bases per mapped file
for i in *.bed; do echo -ne "$i\t"; perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[2] - $s[1];} print "$c\n";' < $i; done
btau_1.ARSUCD.mapped.bed        2222427688
btau_1.btau_1.mapped.bed        2221700003
btau_4.6.ARSUCD.mapped.bed      2470182979
btau_4.6.btau_4.6.mapped.bed    2525002031
umd_3.ARSUCD.mapped.bed 2665236846
umd_3.umd_3.mapped.bed  2659719035

# Now to get the total number of bases in the original assembly
for i in *.fai; do echo -ne "$i\t"; perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < $i; done
ARS-UCD1.2_Btau5.0.1Y.fa.fai    2759153975
Btau20040927-freeze-contigs.fa.fai      2262579180
btau_4.6.karytype.fa.fai        2673141463
umd3_reference_genome.fasta.fai 2670123310
```

| ASM  | Btau1     |  ARS-UCD  |
| :--- | ---:      |  ----:    |
| Btau1| 40879177  | 2221700003|
|ARSUCD| 2222427688| 536726287 |


And plotting the chord diagrams with the above information.

```bash
# OK, new strategy, let's try an upset plot using a script I wrote to make the comparisons easy
python3 ~/python_toolchain/sequenceData/alignmentCoordUpsetLists.py -f ARS-UCD1.2_Btau5.0.1Y.fa.fai -p btau_1_on_ars-ucd.paf -p btau_4.6_on_ars-ucd.paf -p umd_3_on_ars-ucd.paf -o ars-ucd.list -o btau_1_on_ars-ucd.list -o btau_4.6_on_ars-ucd.list -o umd_3_on_ars-ucd.list
```

> F:\SharedFolders\side_projects\cattle_symposium

```R
library(chorddiag)
library(RColorBrewer)

m <- matrix(c(40879177, 2221700003, 2222427688, 536726287), byrow=TRUE, nrow=2, ncol=2)
assemblies <- c("Btau 1.0", "ARS-UCDv1.2")
dimnames(m) <- list(origin = assemblies, destination = assemblies)
colors <- brewer.pal(3, "Dark2")
chorddiag(m, groupColors = colors, groupnamePadding = 20)

# That did not look good! I'll try a different method like upset plots
```

## Generating the assembly comparison histogram plots

OK, I'm going to more simplistic, but hopefully, still illuminating, plots of contiguity and assembly size. Let's start out by gathering the statistics I need for the Histograms.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette/symposium_comparison

```bash
module load bedtools/2.25.0 samtools
# Getting gap regions in the scaffolded assemblies
java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f ARS-UCD1.2_Btau5.0.1Y.fa -o ARS-UCD1.2_Btau5.0.1Y.gaps.bed -s ARS-UCD1.2_Btau5.0.1Y.gaps.stats
java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f btau_4.6.karytype.fa -o btau_4.6.karytype.gaps.bed -s btau_4.6.karytype.gaps.stats

wc -l *.gaps.bed
    404 ARS-UCD1.2_Btau5.0.1Y.gaps.bed	
  66293 btau_4.6.karytype.gaps.bed
  74464 umd3_reference_genome.gaps.bed

# Getting the contig regions
for i in *.fai; do echo $i; perl -lane 'print "$F[0]\t1\t$F[1]";' < $i > $i.bed; done
for i in ARS-UCD1.2_Btau5.0.1Y btau_4.6.karytype umd3_reference_genome; do echo $i; bedtools subtract -a $i.fa*.fai.bed -b $i.gaps.bed > $i.ctgs.bed; done

wc -l *.ctgs.bed
   2614 ARS-UCD1.2_Btau5.0.1Y.ctgs.bed
  65423 btau_4.6.karytype.ctgs.bed
  77554 umd3_reference_genome.ctgs.bed

# And calculating the number of contigs per other assemblies (this is easier)
for i in Btau20040927-freeze-contigs.fa.fai Btau20050310.contigs.fa.fai Btau20060815.contigs.fa.fai; do wc -l $i; done
795212 Btau20040927-freeze-contigs.fa.fai
321107 Btau20050310.contigs.fa.fai
131728 Btau20060815.contigs.fa.fai

# Now to calculate assembly length for everyone (contig/scaffold n50 is the second number)
for i in *.fai; do echo -ne "$i\t"; perl -e '@j = (); $c = 0; while(<>){chomp; @s = split(/\t/); push(@j, $s[1]); $c += $s[1];} @j = sort{$b <=> $a} @j; $n = $c / 2; $n50 = 0; $csum = 0; foreach $k (@j){$csum += $k; if($csum >= $n){$n50 = $k; last;}} print "$c\t$n50\n";' < $i; done
ARS-UCD1.2_Btau5.0.1Y.fa.fai    2759153975      103308737
Btau20040927-freeze-contigs.fa.fai      2262579180      4185
Btau20050310.contigs.fa.fai     2642506507      18884
Btau20060815.contigs.fa.fai     3247500072      48721
btau_4.6.karytype.fa.fai        2673141463      105982576
umd3_reference_genome.fasta.fai 2670123310      105708250
```

Now let's plot this in R.

> F:/SharedFolders/side_projects/cattle_symposium/

```R
library(ggplot2)
library(scales)
data <- data.frame(Assembly = c("Btau_1.0", "Btau_2.0", "Btau_4.0", "Btau_4.6.1", "UMD3.1", "ARS-UCDv1.2"), Length = c(2262579180, 2642506507, 3247500072, 2673141463, 2670123310, 2759153975), Contigs = c(795212, 321107, 131728, 65423, 77554, 2614))

data$Assembly <- as.factor(data$Assembly)
data$Assembly <- factor(data$Assembly, levels(data$Assembly)[c(2:6,1)])

# Assembly lengths
pdf(file="assembly_lengths_plot.pdf", width = 7, height = 4, useDingbats = FALSE)
ggplot(data, aes(x=Assembly, y=Length, fill=Assembly)) + geom_bar(stat="identity") + scale_color_brewer(palette="Dark2") + theme_classic() + geom_text(aes(label=Length, fontface="bold"), vjust=1.6, color="white") + theme(legend.position="none", axis.text.x = element_text(angle=40, hjust=1, size=12), axis.text=element_text(size=12), axis.title=element_text(size=14, face="bold")) + scale_y_continuous(label = unit_format(unit = "Gbp", scale = 1e-9)) + xlab("Assemblies") + ylab("Assembly Lengths (Gbp)")
dev.off()

# Contig counts
pdf(file="assembly_contig_counts.pdf", width = 7, height = 4, useDingbats=FALSE)
ggplot(data, aes(x=Assembly, y=Contigs, color=Assembly)) + geom_bar(stat="identity", fill="white") + scale_color_brewer(palette="Dark2") + theme_classic() + geom_text(aes(label=Contigs, fontface="bold"), vjust=1.6, color="black") + theme(legend.position="none", axis.text.x = element_text(angle=40, hjust=1, size=12), axis.text=element_text(size=12), axis.title=element_text(size=14, face="bold")) + scale_y_continuous(label = comma) + xlab("Assemblies") + ylab("Number of Contigs")
dev.off()
```

## Bipartite plot generation

OK, here's the idea: the [bipartite](https://cran.r-project.org/web/packages/bipartite/bipartite.pdf) package in R has a "plotweb" function that can be used to plot the order and orientation of SNP probes. Maybe I can take the "segs" output from my alignment perl script, create a matrix that shows the interaction of segments against chromosomes, and then 


First, let's remap all of the SNP coordinates and then find the best way to present the data.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/symposium

```bash
module load samtools bwa
# First, let's grab the information we need
mv 5-0108\ Bovine\ Batch\ 1\ 10K\ Panel\ R1_5.txt bovine_10k_panel.txt
dos2unix bovine_10k_panel.txt

wget http://hgdownload.soe.ucsc.edu/goldenPath/bosTau2/bigZips/bosTau2.softmask.fa.gz
gunzip bosTau2.softmask.fa.gz
samtools faidx bosTau2.softmask.fa

perl extract_10k_snps.pl bosTau2.softmask.fa bovine_10k_panel.txt bovine_10k_panel.fa
grep '>' bovine_10k_panel.fa | wc -l
5040  <- not quite 10k placed!

# sorting
perl -e '%d; while($h = <>){chomp $h; $s = <>; @hsegs = split(/\./, $h); $d{$hsegs[1]}->{$hsegs[2]} = [$h, $s];} foreach my $chr (sort {$a <=> $b} keys(%d)){ foreach my $pos (sort {$a <=> $b} keys(%{$d{$chr}})){ $ref = $d{$chr}->{$pos}; print "$ref->[0]\n$ref->[1]";}}' < bovine_10k_panel.fa > bovine_10k_panel.sorted.fa

# And processing
perl ~/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a ../ncbi/ARSUCD1.2.current_ref.fa -p bovine_10k_panel.sorted.fa -o affy10k_to_arsucd1.2

# baseline counts
head -n 1 affy10k_to_arsucd1.2.stats
Mapping probes: 4954    Unmapped probes: 1 <- not too shabby!

# Preparing the BovineSNP50 fasta for remapping
dos2unix snp50_manifest.tab
perl -e '$h = <>; while(<>){chomp; @s = split(/\t/); print ">$s[0].$s[2].$s[3]\n$s[1]\n";}' < snp50_manifest.tab > BovineSNP50.v2.probseq.coords.fa

# Let's queue up the rest, then I have an idea on how to plot the interactions
perl ~/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a ../ncbi/ARSUCD1.2.current_ref.fa -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.rev1coords.fa -o bovinehd_to_arsucd1.2

perl ~/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a ../ncbi/ARSUCD1.2.current_ref.fa -p BovineSNP50.v2.probseq.coords.fa -o bovinesnp50_to_arsucd1.2


# Mapped vs unmapped probe lists
head -n 1 *.stats
==> affy10k_to_arsucd1.2.stats <==
Mapping probes: 5036    Unmapped probes: 2

==> bovinehd_to_arsucd1.2.stats <==
Mapping probes: 734630  Unmapped probes: 1006

==> bovinesnp50_to_arsucd1.2.stats <==
Mapping probes: 54022   Unmapped probes: 35

# Still pretty good. This suggests that most of the SNP probe sequence was designed from the unique regions of the genome, and the transferability of that sequence is easy
# Now to start crafting matrices for bipartite network plotting.
# The Bipartite library in R allows for "tripartite" plots from two different matricies. We could do a comparison of affy vs UMD3 locations vs ARS-ucd
perl createChrMappingMatrix.pl -i affy10k_to_arsucd1.2.tab -o affy10k_to_arsucd1.2.matrix

# OK, it works well enough! Let's run it on all three separately and then concatenate
perl createChrMappingMatrix.pl -i bovinehd_to_arsucd1.2.tab -o bovinehd_to_arsucd1.2.matrix
perl createChrMappingMatrix.pl -i bovinesnp50_to_arsucd1.2.tab -o bovinesnp50_to_arsucd1.2.matrix

# Now to concatenate for the larger plot
perl createChrMappingMatrix.pl -i affy10k_to_arsucd1.2.tab,bovinesnp50_to_arsucd1.2.tab,bovinehd_to_arsucd1.2.tab -o combined_to_arsucd1.2.matrix
```

Let's try to plot the networks as bipartite plots first, and then transition to heatmaps if things get too messy.

```R
library(bipartite)

data <- read.delim("combined_to_arsucd1.2.matrix", header=TRUE, row.names= 1, check.names = FALSE)
matrix <- data.matrix(data, rownames.force = TRUE)

pdf(file="combined_plotweb.pdf", useDingbats = FALSE)
plotweb(matrix, method = "normal")
dev.off()
# Unfortunately, it was a huge mess. Onto the heatmap!
pdf(file="combined_visweb.pdf", useDingbats=FALSE)
visweb(matrix)
dev.off()

# OK, let's try with the smaller chips. There are likely to be more recognizable differences
data <- read.delim("affy10k_to_arsucd1.2.matrix", header=TRUE, row.names= 1, check.names = FALSE)
matrix <- data.matrix(data, rownames.force = TRUE)

pdf(file="affy10k_plotweb.pdf", useDingbats = FALSE)
plotweb(matrix, method = "normal")
dev.off()

pdf(file="affy10k_visweb.pdf", useDingbats = FALSE)
visweb(matrix, type="none")
dev.off()

# And the SNP50
data <- read.delim("bovinesnp50_to_arsucd1.2.matrix", header=TRUE, row.names= 1, check.names = FALSE)
matrix <- data.matrix(data, rownames.force = TRUE)
plotweb(matrix, method = "normal")
pdf("bovinesnp50_plotweb.pdf", useDingbats=FALSE)
plotweb(matrix, method = "normal")
dev.off()

pdf("bovinesnp50_visweb.pdf", useDingbats=FALSE)
visweb(matrix, type="none")
dev.off()
```