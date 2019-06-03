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