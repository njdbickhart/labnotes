# Polled TF investigation
---
*3/27/2019*

These are my notes on trying to identify transcription factor binding sites in the regions of the polled gene for future gene editing experiments.

## Table of contents


## Recap of conversation

We want to try to identify any functional element within the polled locus that may be causal for the phenotype. One possibility is to identify TF factors in the region. I will try to do a liftover from my previous analysis. Worst case -- I may need to rerun TFLOC or find human sequence conservation sites to do this.

I just reviewed the Nguyen et al. manuscript and it appears that their analysis overlapped fairly well with mine. I think that their results are going to be a good approximation of what I did albeit without all of the tuning that I would need to do on my analysis pipeline! I will use their data as a stepping off point. 

| Locus | UMD3 coords | Btau4 coords |
|:--- | :--- | :--- |
|Celtic | chr1:1705834-1706045 | chr1:1,517,223-1,517,434 |


One more bit of information: Musk Deer has a frameshift in RXFP2 which appears to be the gene target for the lncRNA that is downstream of the Celtic mutation. The lncRNA expression is inversely proportional to the expression of RXFP2.

From Tad:

```
RXFP2 is the only the headgear-specific expressed gene (Fig. S11a) and was pseudogenization in both species, with one bp deletion and insertion in the twelfth and twenty-first codon, respectively, result in frame shift (Fig. S10).
```

There may also be a ZEB2 site that is deleted in a new polled mutation. Here's a reference that may contain the [consensus sequence for ZEB2](http://jem.rupress.org/content/212/12/2041).

## Preparing files

I want to gather all of the materials I need. Unfortunately, there's a paucity of chromatin accessibility data in cattle. I want to use this data to validate associations among TF's found in the Australian publication, but the only tissue that the FAANG group uploaded was for [CD-4 cells!](https://www.ebi.ac.uk/ena/data/view/ERX2628476) Oh well... maybe something sticks out? 

I'm going to download a bunch of data for analysis.

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/unmapped_scraping/side_test

```bash
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/000/ERR2611830/ERR2611830_1.fastq.gz"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=1 -p msn --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/000/ERR2611830/ERR2611830_2.fastq.gz"
```

I ran some ContraTF analysis on the entire region and on the Celtic locus. I'm looking for any TF that could be in the region. The 212 bp duplicated in the Celtic locus has some multi-species conservation sites.

The mongolian locus has an enhancer site located in this region: chr1:1974983-1978220 .

More information: This sequence GCCACGCTAGAATTAGATGTCATGAGTGTTTTTATTTAATGACAGTG|CAACTCAATTTAAGAGCAGACTGTGCAGGGCCAAGGGG is conserved in cattle and buffalo but is split by a 1.8 kb insertion. I ran it in [MirBase](http://www.mirbase.org/cgi-bin/blast.pl) and found the following miRNA binding site right on the split site:

```
Query: 35-57 bta-miR-12020: 1-23 score: 79 evalue: 0.43
UserSeq             35  uuuaaugacagugcaacucaauu  57  
                        ||||||  |||||| | ||||||
bta-miR-12020        1  uuuaauaccagugccaaucaauu  23


Query: 38-56 csa-miR-281: 5-23 score: 68 evalue: 3.6
UserSeq             38  aaugacagugcaacucaau  56  
                        ||||| || ||||||| ||
csa-miR-281         23  aaugagagagcaacuccau  5 
```


Check lncRNA Cs1-orf62-as1 binding to RXFP2

OK, let's convert the coordinates to ARS-UCDv1.2 so that I have a good baseline to start with.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver polled_loci.bed liftovers/UMD3.11_to_ARS-UCD1.2/umd3_kary_unmask_ngap_to_ARS-UCD1.2.mmap.liftover.chain polled_loci.ars-ucd.bed polled_loci.ars-ucd.unmapped

```

#### Coordinates

| Locus | UMD3 coords | ARS-UCDv1.2 coords|
| :--- | :--- | :--- |
| Celtic| chr1:1705834-1706045 | NKLS02000001.1:2429109-2429320|
| Mongolian | chr1:1974983-1978220 | NKLS02000001.1:2694744-2697981 |
| RXFP2 | chr12:29234118-29294783 | NKLS02000012.1:29212314-29274747 |

Let's look for lncRNA within the celtic and mongolian regions using the [LGC](http://bigd.big.ac.cn/lgc/calculator) database.

```bash
module load samtools
samtools faidx ncbi/ARS-UCD1.2.PlusY.fa NKLS02000001.1:2429109-2697981 > celtic_mongolian_loci.fa
```

The LGC database does not output coordinates or sequence. I had to blast through an R package to try to find the sequence.

> pwd: F:/SharedFolders/side_projects/tad_polledTF/

```R
library("LncFinder")
library(seqinr)

celtic <- read.fasta("celtic_mongolian_loci.fa")
result_1 <- LncFinder::lnc_finder(celtic)

result_1
                                 Pred Coding.Potential ORF.Max.Len ORF.Max.Cov Seq.lnc.Dist Seq.pct.Dist
NKLS02000001.1:2429109-2697981 Coding        0.8650677        1122 0.004172974    -8.170863    -7.130609
                               Seq.Dist.Ratio Signal.Peak        SNR Signal.Min Signal.Q1 Signal.Q2 Signal.Max
NKLS02000001.1:2429109-2697981       1.145886    160.4746 0.04058561   285.6901  322.7192  377.0941   12555.11

# Psh, what is up with these LncRNA predicting programs? Sure, there's something here, but we're not going to tell you where!

# This found the sequence:
find_orfs(celtic)
```

OK, now that we have the LncRNA sequence, let's check to see how it aligns to ARS-UCDv1.2

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/dominette

```bash
module load minimap2
minimap2 -x asm5 ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa just_testing.fa > just_testing.paf


celticLncRNA    1122    7       992     +       1       158534110       2638171 2639156 243     985     36      tp:A:P  cm:i:17 s1:i:243        s2:i:41   dv:f:0.0030     rl:i:968
flanking        85      36      81      +       12      87216183        29212939        29212984        45      45      5       tp:A:P  cm:i:3  s1:i:45   s2:i:0  dv:f:0.0269     rl:i:0  <- this is the flanking sequence I aligned above
```

This is unexpected, as I thought the insertion sequence site should have been on chr1 because of the mongolian locus. Let's try to get the exact coordinates again from the UMD3 reference using minimap.

```bash
minimap2 -x asm5 repeatmasker/umd3_reference_genome.fasta just_testing.fa > just_testing.umd.paf
celticLncRNA    1122    7       992     +       Chr1    158337067       1918410 1919395 243     985     36      tp:A:P  cm:i:17 s1:i:243        s2:i:41   dv:f:0.0030     rl:i:968
flanking        85      36      81      +       Chr12   91163125        29234743        29234788        45      45      5       tp:A:P  cm:i:3  s1:i:45   s2:i:0  dv:f:0.0269     rl:i:0
```

It's the same. Ah, I see now! The flanking sequence that Tad gave me was from the 5' UTR of cattle RXFP2 and that sheep has a 1.8 kb insertion in that region. Let's see if we can "sheepize" the cattle UTR and see if there are any lncRNA binding sites or other anomalies we can find. 

First, align to the **sheep** to see if its there.

```bash
minimap2 -x asm5 GCA_000298735.2_Oar_v4.0_genomic.fna ../dominette/just_testing.fa > align_to_sheep.paf

cat align_to_sheep.paf
flanking        85      36      81      +       CM001591.2      86377204        29434912        29434957        45      45      5       tp:A:P  cm:i:3    s1:i:45 s2:i:0  dv:f:0.0269     rl:i:0  <- sheep chr10

samtools faidx GCA_000298735.2_Oar_v4.0_genomic.fna CM001591.2:29433000-29434957
# This is pretty much the region! 

# finding the location of the mir site:
perl -e '<>; $s = ""; while(<>){chomp; $s .= uc($_);} $s =~ /AAAACCTTCAGAAGGAAAGGA/; print $-[0] . "\n";' < sheep_region.fa
1729  <- that means that it is at chr10:29434729-29434749
```

Running the region through mirbase (piece by piece because of the 1000 bp limit on queries), I found the following miRNA species with very high identities:

```
Query: 830-850 eca-miR-703: 1-21 score: 105 evalue: 0.013
UserSeq            830  aaaaccuucagaaggaaagga  850 
                        |||||||||||||||||||||
eca-miR-703          1  aaaaccuucagaaggaaagga  21  
Query: 830-850 mmu-miR-703: 1-21 score: 96 evalue: 0.075
UserSeq            830  aaaaccuucagaaggaaagga  850 
                        ||||||||||||||||||| |
mmu-miR-703          1  aaaaccuucagaaggaaagaa  21  

...
Query: 984-1001 bta-miR-2325c: 3-20 score: 72 evalue: 7.4
UserSeq            984  gaaaaaaaaaaaaaaaaa  1001
                        |||||| |||||||| ||
bta-miR-2325c       20  gaaaaagaaaaaaaacaa  3   

# A lower hit bos taurus mir gene
```
