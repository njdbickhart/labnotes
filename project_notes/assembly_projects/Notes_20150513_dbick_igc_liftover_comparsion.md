# Goat assembly Liftover
---
*5/13/2015*

These are my notes on processing Ben's liftover coordinates for the Goat assembly. I will try to resolve some of the issues that he has with the files by scripting custom commands to prevent loss of data from GFF to BED and vice versa.

## Table of Contents
* [Locations](#Locations)
* [Automating split region resolution](#automate)


<a name="Locations"></a>
## Locations
Most of the files are located in this directory:

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover

From Ben's email, he's trying to generate comparative coordinates for the gene models created by the BGI assembly team. Basically -- using the LastZ alignment "chain file," he's trying to map the gene exons back to the BGI assembly for comparison.

Here are the relevant files in this directory:

| File | Comments |
| :--- | :--- |
|NCBI_USDA.over.clean.chain | The chain file to convert BGI coords to PacBio coords |
|ref_CHIR “top_level”.gff | The top BGI gene model calls |
|ref_CHIR “gnomon_top_level”.gff| The less confident BGI gene model calls |
|NCBI_…bed | Ben's GFF to BED conversion |
|USDA_...bed | Ben's liftover to PacBio coords |

My job is to try to get the exon coordinates lifted over, since Ben's bed files are getting trimmed of exon coordinates before and after conversion from GFF.

## Liftover of gene models
OK, let's start digging into the guts of the files here.

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover

```bash
# First, the NCBI GFF
head ref_CHIR_1.0_top_level.gff3
	##gff-version 3
	#!gff-spec-version 1.20
	#!processor NCBI annotwriter
	##sequence-region NC_022293.1 1 155011307
	##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9925
	NC_022293.1     RefSeq  region  1       155011307       .       +       .       ID=id0;Name=1;Dbxref=taxon:9925;breed=Yunnan black goat;chromosome=1;collection-date=Sep-2009;country=China: Kunming;gbkey=Src;genome=chromosome;mol_type=genomic DNA
	NC_022293.1     Gnomon  gene    52833   111810  .       -       .       ID=gene0;Name=LOC102190458;Dbxref=GeneID:102190458;gbkey=Gene;gene=LOC102190458;part=1%2F1
	NC_022293.1     Gnomon  mRNA    52833   111810  .       -       .       ID=rna0;Name=XM_005674738.1;Parent=gene0;Dbxref=GeneID:102190458,Genbank:XM_005674738.1;gbkey=mRNA;gene=LOC102190458;product=chloride intracellular channel protein 6-like;transcript_id=XM_005674738.1
	NC_022293.1     Gnomon  exon    111682  111810  .       -       .       ID=id1;Parent=rna0;Dbxref=GeneID:102190458,Genbank:XM_005674738.1;gbkey=mRNA;gene=LOC102190458;product=chloride intracellular channel protein 6-like;transcript_id=XM_005674738.1
	NC_022293.1     Gnomon  exon    110132  110335  .       -       .       ID=id2;Parent=rna0;Dbxref=GeneID:102190458,Genbank:XM_005674738.1;gbkey=mRNA;gene=LOC102190458;product=chloride intracellular channel protein 6-like;transcript_id=XM_005674738.1

# OK, Let's check out the other GFF files.
# I think I found Ben's problem:
head *.gff3
==> USDA_merged_exon.gff3 <==
utg0    liftOver        exon    12890   13051   0       -       .       ID=id379254
utg0    liftOver        exon    12890   13051   0       -       .       ID=id460240
utg0    liftOver        exon    13867   13896   0       -       .       ID=id454967
utg0    liftOver        exon    17153   17264   0       +       .       ID=id454964
utg0    liftOver        exon    18958   19133   0       +       .       ID=id454965
utg0    liftOver        exon    19501   19653   0       +       .       ID=id454966
utg0    liftOver        exon    45176   45205   0       +       .       ID=id320860
utg0    liftOver        exon    45176   45205   0       +       .       ID=id454972
utg0    liftOver        exon    61213   61374   0       +       .       ID=id320861
utg0    liftOver        exon    61213   61374   0       +       .       ID=id454973

==> USDA_merged_gene-exon.gff3 <==
utg0    liftOver        gene    12890   13051   0       -       .       ID=gene25642
utg0    liftOver        gene    12890   13051   0       -       .       ID=gene57148
utg0    liftOver        gene    90888   107609  0       +       .       ID=gene40688
utg0    liftOver        gene    109288  123269  0       -       .       ID=gene18729
utg0    liftOver        gene    109288  123269  0       -       .       ID=gene40685
utg0    liftOver        gene    123442  134377  0       +       .       ID=gene40684
utg0    liftOver        gene    181741  182236  0       -       .       ID=gene40681
utg0    liftOver        gene    240546  241511  0       -       .       ID=gene18726
utg0    liftOver        gene    240546  241511  0       -       .       ID=gene40678
utg0    liftOver        gene    243295  244484  0       +       .       ID=gene40677
``` 

OK, so Ben is missing the Parent tags on the GFF files. Basically, he's missing ALL the tags! NCBI GFF3 hierarchy looks like it follows this format for gene models:

1. Region
2. gene
3. mRNA
4. exon

##### Looking at existing liftover results

Ben has already tried lifting over the coordinates from a "minimalist" bed file that he generated from the GFF files. Here is an example from the directory.

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover

```bash
wc -l USDA_top_level_exon.bed USDA_top_level_exon_unMapped USDA_top_level_exon_multiple.bed USDA_top_level_exon_unMapped_multiple
  308851 USDA_top_level_exon.bed
   27034 USDA_top_level_exon_unMapped			27034 / 2 = 13,517 unmapped exons
  310757 USDA_top_level_exon_multiple.bed
   24512 USDA_top_level_exon_unMapped_multiple	24512 / 2 = 12,256 unmapped exons, even after multiple mapping
  671154 total
```

So, very few are being lost with the conversion -- that's good! Let's see what bed file Ben is using to perform the liftover.

```bash
head NCBI_merged_exon.bed
	NC_005044.2     0       68      id399503        0       +
	NC_005044.2     68      637     id399504        0       +
	NC_005044.2     1024    1091    id399505        0       +
	NC_005044.2     1091    2664    id399506        0       +
```

OK, so the Bed file only has the ID in the name field! We can do better. Normally I would sock the entire information string into the "name" section of the bed file and parse it later, but I imagine that the memory overhead is going to be huge if I do that (or there are hardcoded string length sizes). Let's deconvolute the gff3 by creating two files (indexable by "ID" field). 

Actually, I think that Ben has already done most of the deconvolution. Let's instead use the base gff3 information field to generate the hierarchy for the subsequent gff3 conversion. Some stats on his files:

```bash
wc -l NCBI_top_level_CDS.bed
	291737 NCBI_top_level_CDS.bed

grep 'CDS' ref_CHIR_1.0_top_level.gff3 | wc -l
	310754	# I thought this should be the same as the line count for the bed file above

grep 'CDS' ref_CHIR_1.0_gnomon_top_level.gff3 | wc -l
	456966

grep 'CDS' ref_CHIR_1.0_merge.gff3 | wc -l
	86039	# How is this smaller??
```

OK, that's really concerning. **There are 19,017 CDS entries missing from Ben's bed compared to a straight pull from the BGI GFF3 file.** Time to do a brute force check.

```bash
# Creating my own bed from the CDS entries of the GFF3 file
grep 'CDS' ref_CHIR_1.0_top_level.gff3 | perl -lane '($id) = $F[8] =~ /ID=(.{4,10});.+/; print "$F[0]\t$F[3]\t$F[4]\t$id\t$F[7]\t$F[6]";' | ~/bin/sortBedFileSTDIN.pl > BGI.top_level.CDS.bed

# Resorting Ben's bed file for bedtools comparison
cat NCBI_top_level_CDS.bed | sortBedFileSTDIN.pl > NCBI_top_level_CDS.lexsort.bed

# Checking if Ben's bed has any unique entries that overlap my entries
intersectBed -a NCBI_top_level_CDS.lexsort.bed -b BGI.top_level.CDS.bed -v | wc -l
	0	# None

# Checking if there are differences in area between Ben's bed and mine
subtractBed -a NCBI_top_level_CDS.lexsort.bed -b BGI.top_level.CDS.bed | bed_length_sum.pl
        Interval Numbers:       271623
        Total Length:           271623
        Length Average:         1
        Length Median:          1
        Min Length:             1
        Max Length:             1
        Length Stdev:           0

# That's just a one base difference per interval -- perhaps his conversion program was assuming one-base status per entry?
# OK, now the reverse intersect -v
intersectBed -a BGI.top_level.CDS.bed -b NCBI_top_level_CDS.lexsort.bed -v | wc -l
	365 # OK, this is far less than I expected! Let's look at the output

intersectBed -a BGI.top_level.CDS.bed -b NCBI_top_level_CDS.lexsort.bed -v | head
	NC_022305.1     3161658 3162149 id162557        .       +
	NC_022305.1     3183380 3183471 id162558        .       +
	NC_022305.1     15082374        15086202        id163131        .       +
	NC_022305.1     41743720        41744057        id166409        .       -
	NC_022305.1     50310337        50310508        id167359        .       +

# Hmm... and if I return to the original gff file, here are some of the results I find:
grep 'id162557' ref_CHIR_1.0_top_level.gff3
	NC_022305.1     .... ID=id162557;Parent=rna14728;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred complete CDS: inserted 5 bases in 3 codons%3B deleted 1 base in 1 codon;Dbxref=GeneID:102189493,Genbank:XM_005687833.1;exception=unclassified transcription discrepancy;gbkey=mRNA;gene=SLX4IP;product=SLX4 interacting protein;transcript_id=XM_005687833.1

grep 'id162558' ref_CHIR_1.0_top_level.gff3
	NC_022305.1     .... ID=id162558;Parent=rna14728;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred complete CDS: inserted 5 bases in 3 codons%3B deleted 1 base in 1 codon;Dbxref=GeneID:102189493,Genbank:XM_005687833.1;exception=unclassified transcription discrepancy;gbkey=mRNA;gene=SLX4IP;product=SLX4 interacting protein;transcript_id=XM_005687833.1

grep 'id167359' ref_CHIR_1.0_top_level.gff3
	NC_022305.1     .... ID=id167359;Parent=rna15173;Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred complete CDS: inserted 2 bases in 1 codon;Dbxref=GeneID:102183079,Genbank:XM_005688207.1;exception=unclassified transcription discrepancy;gbkey=mRNA;gene=LOC102183079;product=leucine zipper putative tumor suppressor 3-like;transcript_id=XM_005688207.1

# OK, so these all have "exception" tags within the info field of the gff
# Let's count how many CDS entries have "exceptions"
grep 'CDS' ref_CHIR_1.0_top_level.gff3 | grep 'exception' | wc -l
	36,698 # This is still more entries than are missing in Ben's bed file... there's something screwy going on here
```

##### Starting from scratch

OK, I think that the best way forward is to naively trust the NCBI files and process them with one-off scripts in order to ensure that there are no data fidelity issues with the output. I think that I can create a one-shot script that can process the liftover bed files using the GFF3 file info fields. 

First order of business is to generate the bed files that I need-- I've already created the CDS file from the commands above.

```bash
grep 'exon' ref_CHIR_1.0_top_level.gff3 | perl -lane '($id) = $F[8] =~ /ID=(.{4,10});.+/; print "$F[0]\t$F[3]\t$F[4]\t$id\t$F[7]\t$F[6]";' | ~/bin/sortBedFileSTDIN.pl > BGI.top_level.exon.bed

grep 'gene' ref_CHIR_1.0_top_level.gff3 | perl -lane '($id) = $F[8] =~ /ID=(.{4,10});.+/; print "$F[0]\t$F[3]\t$F[4]\t$id\t$F[7]\t$F[6]";' | ~/bin/sortBedFileSTDIN.pl > BGI.top_level.gene.bed

grep 'mRNA' ref_CHIR_1.0_top_level.gff3 | perl -lane '($id) = $F[8] =~ /ID=(.{4,10});.+/; print "$F[0]\t$F[3]\t$F[4]\t$id\t$F[7]\t$F[6]";' | ~/bin/sortBedFileSTDIN.pl > BGI.top_level.mRNA.bed

grep 'tRNA' ref_CHIR_1.0_top_level.gff3 | perl -lane '($id) = $F[8] =~ /ID=(.{4,10});.+/; print "$F[0]\t$F[3]\t$F[4]\t$id\t$F[7]\t$F[6]";' | ~/bin/sortBedFileSTDIN.pl > BGI.top_level.tRNA.bed

wc -l BGI.top_level.*.bed
  310754 BGI.top_level.CDS.bed
  322615 BGI.top_level.exon.bed
  671363 BGI.top_level.gene.bed
  333609 BGI.top_level.mRNA.bed
    8689 BGI.top_level.tRNA.bed
```

I'm going to keep this run tidy by creating sub directory structures. Top level directory will contain these parsed bed files and the lower level ones will contain the liftover results. 

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover/BGI_top_level

```bash
mkdir CDS_liftover
liftOver BGI.top_level.CDS.bed ../NCBI_USDA.over.clean.chain CDS_liftover/PacBio.top_level.CDS.over.bed CDS_liftover/PacBio.top_level.CDS.over.bed.unmapped

liftOver -multiple BGI.top_level.CDS.bed ../NCBI_USDA.over.clean.chain CDS_liftover/PacBio.top_level.CDS.over.mult.bed CDS_liftover/PacBio.top_level.CDS.over.mult.bed.unmapped

wc -l CDS_liftover/*
  298152 CDS_liftover/PacBio.top_level.CDS.over.bed
   25204 CDS_liftover/PacBio.top_level.CDS.over.bed.unmapped
  302069 CDS_liftover/PacBio.top_level.CDS.over.mult.bed
   22700 CDS_liftover/PacBio.top_level.CDS.over.mult.bed.unmapped
  648125 total
```

**Wait one second, I noticed in the liftover help menu that it takes gff files.** Let's see if we can do the conversion directly.

```bash
liftOver -gff ../ref_CHIR_1.0_top_level.gff3 ../NCBI_USDA.over.clean.chain CDS_liftover/PacBio.top_level.gff.over.bed CDS_liftover/PacBio.top_level.gff.over.bed.unmapped
	# This was the STDERR I received
	Reading liftover chains
	Mapping coordinates
	WARNING: -gff is not recommended.
	Use 'ldHgGene -out=<file.gp>' and then 'liftOver -genePred <file.gp>'

# OK!!! It worked, but I guessed wrong: the file format is gff, not bed
mv CDS_liftover/PacBio.top_level.gff.over.bed CDS_liftover/PacBio.top_level.over.gff
mv CDS_liftover/PacBio.top_level.gff.over.bed.unmapped CDS_liftover/PacBio.top_level.over.gff.unmapped

wc -l CDS_liftover/PacBio.top_level.over.gff*
   842681 CDS_liftover/PacBio.top_level.over.gff
   155722 CDS_liftover/PacBio.top_level.over.gff.unmapped  155722 / 2 = 77861
	# 77861 + 842681 = 920,542

wc -l ../ref_CHIR_1.0_top_level.gff3
	920,542 ../ref_CHIR_1.0_top_level.gff3
```

OK, the liftover worked well. Let's do the same command, but with multiple mapping

```bash
liftOver -multiple -gff ../ref_CHIR_1.0_top_level.gff3 ../NCBI_USDA.over.clean.chain CDS_liftover/PacBio.top_level.over.mult.gff CDS_liftover/PacBio.top_level.over.mult.gff.unmapped

wc -l CDS_liftover/PacBio.top_level.over.*
   842681 CDS_liftover/PacBio.top_level.over.gff
   155722 CDS_liftover/PacBio.top_level.over.gff.unmapped
   842681 CDS_liftover/PacBio.top_level.over.mult.gff
   155722 CDS_liftover/PacBio.top_level.over.mult.gff.unmapped

# They're the same
```

I suspect that liftover ignores the **-multiple** command when non-bed input is found (overriding options?). Well, even despite the lack of multiple map entries, we have a good start. 

I'm still seeing issues with numbers between Ben's count of gene models and information from the gff files. Here's an example:

```bash
# Ben said that the top_level file should have 25,789 gene models
for i in ../ref_CHIR*.gff3; do echo -n "$i "; grep 'gene' $i | wc -l; done
../ref_CHIR_1.0_gnomon_top_level.gff3 	1016475
../ref_CHIR_1.0_merge.gff3 				557002
../ref_CHIR_1.0_top_level.gff3 			671363 # vs 25,789

for i in ../ref_CHIR*.gff3; do echo -n "$i "; grep 'region' $i | wc -l; done
../ref_CHIR_1.0_gnomon_top_level.gff3 	0
../ref_CHIR_1.0_merge.gff3 				72883
../ref_CHIR_1.0_top_level.gff3 			155301 # vs 25,789

perl -e '%h; while(<>){chomp; @s = split(/\t/); $h{$s[2]} +=1;} foreach my $k (keys(%h)){print "$k\t$h{$k}\n";}' < ../ref_CHIR_1.0_top_level.gff3
```
|gff tag | count |
| :--- | ---: |
region|  77022
CDS    | 291737
**Null field**       | 154046
D_loop | 1
V_gene_segment | 265
cDNA_match     | 17731
transcript     | 1549
match   |333
gene    |25789
mRNA    |27947
exon    |322368
tRNA    |1664
rRNA    |2
C_gene_segment  |88

Ah ok, so it's just that the gff file is not "grep safe." So my "grep" of **gene** and **CDS** tags were pretty sloppy with respect to the field fidelity. 

*5/14/2015*

--

OK, now to quickly run the liftover on the gnomon results as well. 

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover

```bash
# Just moving to a new working directory to keep things organized.
mkdir BGI_gnomon
cd BGI_gnomon/

# The liftover commands
liftOver -gff ../ref_CHIR_1.0_gnomon_top_level.gff3 ../NCBI_USDA.over.clean.chain PacBio.gnomon.over.gff PacBio.gnomon.over.gff.unmapped
liftOver -multiple -gff ../ref_CHIR_1.0_gnomon_top_level.gff3 ../NCBI_USDA.over.clean.chain PacBio.gnomon.over.mult.gff PacBio.gnomon.over.mult.gff.unmapped

``` 
I just wrote a new script to process the gff file to generate automatic counts of gff tab column entries. Should be much more maintainable than my one-liners. 

```bash
# Here's how to invoke the script
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f PacBio.gnomon.over.gff -c 2 -e '#' -m
```
|Entry          |  Count|
|:--------------|------:|
|CDS            | 281020|
|C_gene_segment |     76|
|V_gene_segment |    187|
|cDNA_match     |  16380|
|exon           | 308851|
|gene           |  16753|
|mRNA           |  17924|
|match          |    290|
|rRNA           |      1|
|region         |  44477|
|tRNA           |   1592|
|transcript     |   1084|

**Even better: I added multiple file comparison options:**
```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f PacBio.gnomon.over.gff,PacBio.gnomon.over.mult.gff -c 2 -e '#' -m
```
File1:  PacBio.gnomon.over.gff
File2:  PacBio.gnomon.over.mult.gff

|File1          |  Count|File2          |  Count|
|:--------------|------:|:--------------|------:|
|CDS            | 281020|CDS            | 281020|
|C_gene_segment |     76|C_gene_segment |     76|
|V_gene_segment |    187|V_gene_segment |    187|
|cDNA_match     |  16380|cDNA_match     |  16380|
|exon           | 308851|exon           | 308851|
|gene           |  16753|gene           |  16753|
|mRNA           |  17924|mRNA           |  17924|
|match          |    290|match          |    290|
|rRNA           |      1|rRNA           |      1|
|region         |  44477|region         |  44477|
|tRNA           |   1592|tRNA           |   1592|
|transcript     |   1084|transcript     |   1084|


##### Comparison of liftover gff results

```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f BGI_gnomon/PacBio.gnomon.over.gff,BGI_top_level/CDS_liftover/PacBio.top_level.over.gff -c 2 -e '#' -m
```

File1:  BGI_gnomon/PacBio.gnomon.over.gff
File2:  BGI_top_level/CDS_liftover/PacBio.top_level.over.gff

|File1      |  Count|File2          |  Count|
|:----------|------:|:--------------|------:|
|CDS        | 406755|CDS            | 281020|
|exon       | 436058|C_gene_segment |     76|
|gene       |  38423|V_gene_segment |    187|
|mRNA       |  43119|cDNA_match     |  16380|
|transcript |   4188|exon           | 308851|
|           |       |gene           |  16753|
|           |       |mRNA           |  17924|
|           |       |match          |    290|
|           |       |rRNA           |      1|
|           |       |region         |  44477|
|           |       |tRNA           |   1592|
|           |       |transcript     |   1084|

## Resolution of Unmapped liftover results

The goal here is to take a look at the unmapped results to get basic statistics on what carried over from the BGI annotation, and what didn't. In order to get an eagle-eye view of the situation first, let's count the number of split vs deleted entries.

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/liftover/BGI_top_level/CDS_liftover

```bash
perl -e '%c; while($h = <>){chomp $h; $h =~ s/\# //g; $g = <>; chomp $g; @s = split(/\t/, $g); $c{$h}->{$s[2]} += 1;} foreach $t (sort{$a cmp $b} keys %c){foreach $n (sort{$a cmp $b} keys %{$c{$t}}){print "$t\t\|$n\t\|". $c{$t}->{$n}; print "\n";}}' < PacBio.top_level.over.gff.unmapped
```

|liftover | type | count|
|:--- | :--- | ---:|
Deleted in new  |CDS    |2295
Deleted in new  |C_gene_segment |7
Deleted in new  |V_gene_segment |41
Deleted in new  |cDNA_match     |186
**Deleted in new**  |**exon**   |**2731**
**Deleted in new**  |**gene**   |**103**
Deleted in new  |mRNA   |33
Deleted in new  |match  |8
Deleted in new  |region |6373
Deleted in new  |tRNA   |19
Deleted in new  |transcript     |4
Partially deleted in new        |CDS    |7802
Partially deleted in new        |C_gene_segment |4
Partially deleted in new        |D_loop |1
Partially deleted in new        |V_gene_segment |29
Partially deleted in new        |cDNA_match     |1097
**Partially deleted in new**        |**exon**   |**9838**
**Partially deleted in new**        |**gene**   |**3310**
Partially deleted in new        |mRNA   |3228
Partially deleted in new        |match  |25
Partially deleted in new        |region |23125
Partially deleted in new        |tRNA   |43
Partially deleted in new        |transcript     |183
Split in new    |CDS    |620
Split in new    |C_gene_segment |1
Split in new    |V_gene_segment |8
Split in new    |cDNA_match     |68
**Split in new**    |**exon**   |**948**
**Split in new**    |**gene**   |**5623**
Split in new    |mRNA   |6762
Split in new    |match  |10
Split in new    |rRNA   |1
Split in new    |region |3047
Split in new    |tRNA   |10
Split in new    |transcript     |278

OK, lots of split genes and exons. Not as many deleted genes and exons -- **I wonder if Ben also ran LastZ on the "deg" contig file?** Let's start with the split genes/exons. If they lie on contig borders (less likely as there are nearly as many genes as there are contig borders) then they might be evidence for scaffolding. Let's look at exons first, as I can do a BWA MEM on the shorter sequence of them to see finer-tuned alignments.

```bash
perl -e '%c; while($h = <>){chomp $h; $h =~ s/\# //g; $g = <>; chomp $g; @s = split(/\t/, $g); if($h =~ /Split in new/ && $s[2] eq "exon"){print "$g\n";} $c{$h}->{$s[2]} += 1;}' < PacBio.top_level.over.gff.unmapped | head
```

> NC_022293.1     Gnomon  exon    110132  110335  .       -       .       ID=id2;Parent=rna0;Dbxref=GeneID:102190458,Genbank:XM_005674738.1;gbkey=mRNA;gene=LOC102190458;product=chloride intracellular channel protein 6-like;transcript_id=XM_005674738.1
NC_022293.1     Gnomon  exon    4172533 4173007 .       -       .       ID=id671;Parent=rna65;Dbxref=GeneID:102179518,Genbank:XM_005674697.1;gbkey=mRNA;gene=LOC102179518;product=keratin-associated protein 6-2-like;transcript_id=XM_005674697.1
NC_022293.1     Gnomon  exon    4187183 4187943 .       +       .       ID=id674;Parent=rna66;Dbxref=GeneID:102179810,Genbank:XM_005674698.1;gbkey=mRNA;gene=LOC102179810;product=keratin-associated protein 6-2-like;transcript_id=XM_005674698.1
NC_022293.1     Gnomon  exon    4753646 4754511 .       -       .       ID=id682;Parent=rna72;Dbxref=GeneID:102181431,Genbank:XM_005674703.1;gbkey=mRNA;gene=LOC102181431;product=keratin-associated protein 13-1-like;transcript_id=XM_005674703.1
NC_022293.1     Gnomon  exon    24615276        24616184        .       -       .       ID=id1406;Parent=gene113;Dbxref=GeneID:102183454;exon_number=1;gbkey=exon;gene=LOC102183454;part=1%2F1

**NC_022293.1 is chr1 from an NCBI search.** So let's faidx chr1 and start going to town.

> Blade14:/mnt/nfs/nfs2/GoatData/goat_annotation/dbickhart/bwa_comp

```bash
cp ../../genomes/NCBI_BGI/unmasked/chi_ref_CHIR_1.0_chr1.fa ./

# Getting the BGI sequence for the coordinates
samtools faidx chi_ref_CHIR_1.0_chr1.fa
head chi_ref_CHIR_1.0_chr1.fa.fai
	# gi|541128986|ref|NC_022293.1|   155011307       122     70      71

samtools faidx chi_ref_CHIR_1.0_chr1.fa 'gi|541128986|ref|NC_022293.1|:110132-110335'
	>gi|541128986|ref|NC_022293.1|:110132-110335
	CTTCACGAAGAGGGTGATGTCGTGGTCCTGCCCCGGGGCCCCCTCCTCGGACGCCTCGCC
	GCCCTCCTGGCCGCCGTTCTCGGGCACTGTCCCGCCGCCGCGCTCCGCGCTGCCCTCGCC
	CGCCCGGTGCCTCCCCCCCCTCCCCCCCCCAGCCACATCCATTAGGCCCCTTGTTCCCAC
	AACCTGCTGGCTGGGCACCTTGTC

samtools faidx chi_ref_CHIR_1.0_chr1.fa 'gi|541128986|ref|NC_022293.1|:110132-110335' > top_nolift_exon_id2.fa

# Now to pull the PacBio assembly and generate a BWA index for it
# Actually, it looks like Ben has already done this:
ls ../../genomes/USDA/USDA_V3.fasta*
	#../../genomes/USDA/USDA_V3.fasta      ../../genomes/USDA/USDA_V3.fasta.ann  ../../genomes/USDA/USDA_V3.fasta.pac
	#../../genomes/USDA/USDA_V3.fasta.amb  ../../genomes/USDA/USDA_V3.fasta.bwt  ../../genomes/USDA/USDA_V3.fasta.sa

# OK, BWA MEM time
# Ben is hitting the mount a bit hard, so disk IO is way too slow
# copying over to a new share
cp ../../genomes/USDA/USDA_V3.fasta /mnt/iscsi/vnx_gliu_7/goat_assembly/
```
> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/

```bash
bwa index USDA_V3.fasta
samtools faidx USDA_V3.fasta

bwa mem USDA_V3.fasta /mnt/nfs/nfs2/GoatData/goat_annotation/dbickhart/bwa_comp/top_nolift_exon_id2.fa
	...
	gi|541128986|ref|NC_022293.1|:110132-110335     16      utg2462        294089  60 80S86M1I37M     *       0       0       GACAAGGTGCCCAGCCAGCAGGTTGTGGGAACAAGGGGCCTAATGGATGTGGCTGGGGGGGGGAGGGGGGGGAGGCACCGGGCGGGCGAGGGCAGCGCGGAGCGCGGCGGCGGGACAGTGCCCGAGAACGGCGGCCAGGAGGGCGGCGAGGCGTCCGAGGAGGGGGCCCCGGGGCAGGACCACGACATCACCCTCTTCGTGAAG      *       NM:i:5  MD:Z:5A4C7T2A101        AS:i:96   XS:i:20 SA:Z:utg3802_len=1658306_reads=6343_status=X_microHet=1.00_covStat=1.00_quiver,1039863,+,129S75M,60,2;
	gi|541128986|ref|NC_022293.1|:110132-110335     2048    utg3802       1039863 60 129H75M *       0       0       CCTCCCCCCCCTCCCCCCCCCAGCCACATCCATTAGGCCCCTTGTTCCCACAACCTGCTGGCTGGGCACCTTGTC     *       NM:i:2  MD:Z:5T8A60     AS:i:65   XS:i:22 SA:Z:utg2462_len=532208_reads=2354_status=X_microHet=1.00_covStat=1.00_quiver,294089,-,80S86M1I37M,60,5;
	
# INTERESTING NOTE: 
# BWA Mem's SA:Z tags are circular references to MEM split reads!
```

It's pretty obvious to me that this is a contig break. Here's why: 

**CCTCCCCCCCCTCCCCCCCCC** 

Pretty easy to spot simple repeat and they probably just patched a bandaid over the region without regard for the consequences. When I put the 128 bases that map to utg2462 in a NCBI nr blast, I get these results:

> PREDICTED: Capra hircus chloride intracellular channel protein 6-like (LOC102190458), mRNA
Sequence ID: ref|XM_005674738.1|Length: 1192Number of Matches: 1
Related Information
Gene-associated gene details
Map Viewer-aligned genomic context
Range 1: 206 to 333GenBankGraphicsNext MatchPrevious Match
Alignment statistics for match #1
Score	Expect	Identities	Gaps	Strand
237 bits(128)	3e-59	128/128(100%)	0/128(0%)	Plus/Plus

>PREDICTED: Ovis aries musimon chloride intracellular channel protein 6 (LOC101112021), transcript variant X2, mRNA
Sequence ID: ref|XM_012170865.1|Length: 4828Number of Matches: 1
Related Information
Gene-associated gene details
Range 1: 2046 to 2173GenBankGraphicsNext MatchPrevious Match
Alignment statistics for match #1
Score	Expect	Identities	Gaps	Strand
198 bits(107)	1e-47	121/128(95%)	0/128(0%)	Plus/Plus

Hahaha! A blast with the full sequence against the Non-redundant database shows that the BGI assembly is split among closely related ruminants too! So that explains the discrepancy. **I need to find a way to automate this search by pulling sequence from all of the split exons and by mapping them back to the assembly with BWA mem.**

--
##### Quick gff conversion

Ben needs a quick modification of a gff file in order to give "dummy" exon coordinates. I need to create a duplicate "exon" coordinate that is the child of the "CDS" coordinate in each case.

> Blade14: /mnt/nfs/nfs2/GoatData/goat_annotation/EVM/working

```bash
head gene_predictions.gff3
	utg0    AUGUSTUS_masked gene    12203   13931   0.15    -       .       ID=g1;
	utg0    AUGUSTUS_masked mRNA    12203   13931   0.15    -       .       ID=g1.t1;Parent=g1
	utg0    AUGUSTUS_masked stop_codon      12203   12205   .       -       0       Parent=g1.t1;
	utg0    AUGUSTUS_masked CDS     12203   12653   0.8     -       1       ID=g1.t1.CDS1;Parent=g1.t1

# NOTE: the gff doesn't have '#' header tags like the NCBI version
# I can easily write a one-liner for this
perl -lane 'if($F[2] eq "CDS"){$F[8] =~ /ID=(.+);Parent=(.+$)/; $id = $1; $parent = $2; $parent = $id; $id =~ s/CDS/exon/; print join("\t", @F); print "$F[0]\t$F[1]\texon\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\tID=$id;Parent=$parent";}else{print $_;}' < gene_predictions.gff3 | head -n 20
	utg0    AUGUSTUS_masked gene    12203   13931   0.15    -       .       ID=g1;
	utg0    AUGUSTUS_masked mRNA    12203   13931   0.15    -       .       ID=g1.t1;Parent=g1
	utg0    AUGUSTUS_masked stop_codon      12203   12205   .       -       0       Parent=g1.t1;
	utg0    AUGUSTUS_masked CDS     12203   12653   0.8     -       1       ID=g1.t1.CDS1;Parent=g1.t1
	utg0    AUGUSTUS_masked exon    12203   12653   0.8     -       1       ID=g1.t1.exon1;Parent=g1.t1.CDS1
	utg0    AUGUSTUS_masked intron  12654   12730   0.82    -       .       Parent=g1.t1;

# Bingo. Now to make the file.
```

<a name="automate"></a>
## Automating split region resolution

*5/20/2015*

--
I need to generate a script that reads the unmapped liftover coordinates and remaps the split regions to try to find assembly discrepancies. 

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/liftover_comp

```bash
# First, let's get the files that I need in this directory
cp /mnt/nfs/nfs2/GoatData/goat_annotation/liftover/BGI_top_level/CDS_liftover/PacBio.top_level.over.gff.unmapped ./

# I need to create an indexable reference genome with chromosome names that match the liftover chromosomes
# NCBI makes this as painful as possible, sadly
for i in /mnt/nfs/nfs2/GoatData/goat_annotation/genomes/NCBI_BGI/unmasked/chi_ref_CHIR_1.0_chr*.fa; do perl -e 'while(<>){chomp; if($_ =~ />/){($c) = $_ =~ />gi\|\d+\|ref\|(NC_.+)\| .+/; print ">$c\n";}else{print "$_\n";}}' < $i; done > chir_bgi_pseudo_chr.fa

# The above script should remove the pipes and other silly characters that they use in their scaffold naming scheme
# Now to get the stats and index the chromosome
samtools faidx chir_bgi_pseudo_chr.fa
# NC_005044.2 is chrMt
bwa index chir_bgi_pseudo_chr.fa

# I'm going to pull all of the unmapped, split exons, generate fastas for them, and then remap them back.
perl -e '%c; while($h = <>){chomp $h; $h =~ s/\# //g; $g = <>; chomp $g; @s = split(/\t/, $g); ($id) = $s[8] =~ /ID=(.+);Parent/; if($h =~ /Split in new/ && $s[2] eq "exon"){open(IN, "samtools faidx chir_bgi_pseudo_chr.fa $s[0]:$s[3]-$s[4] |"); $fh = <IN>; chomp $fh; if($fh eq ""){next;} print "$fh\_$id\n"; while($line = <IN>){print $line;} close IN;}}' < PacBio.top_level.over.gff.unmapped > split_in_new_exons.fa

# The read names should contain the exon IDs from the gff file. 
# This will allow me to refer back to the annotation if needed.
bwa mem -M -t 5 /mnt/nfs/nfs2/GoatData/goat_annotation/genomes/USDA/USDA_V3.fasta split_in_new_exons.fa > split_in_new_exons.sam

samtools view -bS split_in_new_exons.sam >  split_in_new_exons.bam

```
*5/21/2015*

--

Some really stupid mistakes that BGI made with the annotation:

**NC_022293.1:140352356-140353005_id11935** is a "vegetative cell wall protein, gp1-like" from their annotation, and it was clearly lifted from a mine-drainage metagenomics study. When you look at their assembled sequence, it's clearly a simple repetitive element.

I'm going to filter the alignments to try to get tabular data on each exon region.

```bash
perl ~/perl_toolchain/sequence_data_scripts/resolveLiftoverDifferences.pl -i split_in_new_exons.bam -o split_in_new_exons_format.tab

# This file has the remapped locations for these regions.
# Now, I'm going to get some basic stats about the bam, like number of split reads, etc.

samtools view -f '0x100' split_in_new_exons.bam | wc -l
	125	# This is the number of secondary alignments from the realignment
samtools view -f '0x100' split_in_new_exons.bam | perl -lane '($h) = $F[5] =~ /(\d+)H/; if($h){print $h;}' | statStd.pl
	total   125
	Minimum 13
	Maximum 5004
	Average 753.776000   # Average number of Hardclipped bases from alignment
	Median  377
	Standard Deviation      1045.663456
	Mode(Highest Distributed Value) 130

```

## Comparative genomics alignment with Sheep

Since I had such a great success with that one test case, I'm interested in taking ALL the Sheep exon nucleotide sequence, blasting it against BGI, and blasting it against the PacBio assembly. THEN, I plan to see which exons fall out in BGI vs PacBio and vice versa.

>Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/liftover_comp

```bash
# Downloading the data from Ensembl
wget ftp://ftp.ensembl.org/pub/release-79/gtf/ovis_aries/Ovis_aries.Oar_v3.1.79.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-79/fasta/ovis_aries/dna/Ovis_aries.Oar_v3.1.dna.toplevel.fa.gz

# Here is a test script to help me pull the exon fastas from sheep
head Ovis_aries.Oar_v3.1.79.gtf | perl -ne '@F = split(/\t/); if($F[0] =~ /^#/){next;}else{($i) = $F[8] =~ /exon_id \"(.{15,20})\"\;/; print $i; print "\n";}'

# Now to pull the ovis exon fastas
samtools faidx Ovis_aries.Oar_v3.1.dna.toplevel.fa
cat Ovis_aries.Oar_v3.1.79.gtf | perl -ne '@F = split(/\t/); if($F[0] =~ /^#/){next;}elsif($F[2] eq "exon"){($i) = $F[8] =~ /exon_id \"(.{12,20})\"\;/; open(IN, "samtools faidx Ovis_aries.Oar_v3.1.dna.toplevel.fa $F[0]:$F[3]-$F[4] |"); my $h = <IN>; chomp $h; $h =~ s/>//g; $h .= "_" . $i; print ">$h\n"; while(<IN>){print $_;} close IN;;}' > ovis.aries.exon.nuc.fa

# Aligning first to BGI
bwa mem -M chir_bgi_pseudo_chr.fa ovis.aries.exon.nuc.fa > oar_bgi_exon_align.sam
samtools view -bS oar_bgi_exon_align.sam > oar_bgi_exon_align.bam

bwa mem -M USDA_V3.fasta ovis.aries.exon.nuc.fa > oar_pacbio_exon_align.sam
samtools view -bS oar_pacbio_exon_align.sam > oar_pacbio_exon_align.bam
```

OK, now to calculate the stats of each alignment.

```bash
# PacBio dataset
perl ~/perl_toolchain/sequence_data_scripts/getSamFlagStats.pl oar_pacbio_exon_align.bam
```
|Flag                     |  Count|
|:------------------------|------:|
|TotalReads               | 237632|
|SegUnmapped              |   9930|
|ThisSeqReverse           | 113665|
|SecondaryAlignment       |   3127|

```bash
# BGI dataset
perl ~/perl_toolchain/sequence_data_scripts/getSamFlagStats.pl oar_bgi_exon_align.bam
```

|Flag                     |  Count|
|:------------------------|------:|
|TotalReads               | 239552|
|SegUnmapped              |  19172|
|ThisSeqReverse           |  38866|
|SecondaryAlignment       |   5047|

OK, this is good! The PacBio dataset mapped more (~10,000 more mapped) and split fewer exons (~ 2,000 fewer split) from Sheep. Let's create venns to display these differences. I'm going to produce simple lists of the Ensembl exon IDs for each category, and then run my perl script to create the Venn comparison.

```bash
# Just keeping things tidy
mkdir venn_comparisons
samtools view oar_bgi_exon_align.bam | perl -lane 'if(($F[1] & 4) == 4){@b = split(/_/, $F[0]); print $b[1];}' > venn_comparisons/oar_bgi_unmapped_exons.txt
samtools view oar_bgi_exon_align.bam | perl -lane 'if(($F[1] & 256) == 256){@b = split(/_/, $F[0]); print $b[1];}' > venn_comparisons/oar_bgi_split_exons.txt

samtools view oar_pacbio_exon_align.bam | perl -lane 'if(($F[1] & 4) == 4){@b = split(/_/, $F[0]); print $b[1];}' > venn_comparisons/oar_pacbio_unmapped_exons.txt
samtools view oar_pacbio_exon_align.bam | perl -lane 'if(($F[1] & 256) == 256){@b = split(/_/, $F[0]); print $b[1];}' > venn_comparisons/oar_pacbio_split_exons.txt

wc -l venn_comparisons/*
  5047 venn_comparisons/oar_bgi_split_exons.txt
 19172 venn_comparisons/oar_bgi_unmapped_exons.txt
  3127 venn_comparisons/oar_pacbio_split_exons.txt
  9930 venn_comparisons/oar_pacbio_unmapped_exons.txt

# OK, everything matches up so far! Let's go to town with the Venn!

# First, the split reads
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl venn_comparisons/oar_bgi_split_exons.txt venn_comparisons/oar_pacbio_split_exons.txt
	File Number 1: venn_comparisons/oar_bgi_split_exons.txt
	File Number 2: venn_comparisons/oar_pacbio_split_exons.txt
	Set     Count
	1       1930
	1;2     1528
	2       592

# Important note!!! I forgot that the split reads had mulitple mappings
# My Venn script reduces down the names to unique entries, so the numbers sum up to be less here

# Now, the unmapped reads
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl venn_comparisons/oar_bgi_unmapped_exons.txt venn_comparisons/oar_pacbio_unmapped_exons.txt
	File Number 1: venn_comparisons/oar_bgi_unmapped_exons.txt
	File Number 2: venn_comparisons/oar_pacbio_unmapped_exons.txt
	Set     Count
	1       9534
	1;2     9140
	2       564

# I added a "printout" function to the venn script that prints the names in each set to text files
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o venn_comparisons/oar_bgi_unmapped_exons.txt venn_comparisons/oar_pacbio_unmapped_exons.txt

# I'm going to work on the BGI exclusive unmapped exons first
# Let's see how they map in the PacBio file
samtools view oar_pacbio_exon_align.bam | grep 'ENSOARE00000135090'
	11:26767193-26767365_ENSOARE00000135090 0       utg3197_len=152401_reads=517_status=X_microHet=1.00_covStat=1.00_quiver 38336   60      8M1I30M1I9M1I83M1I39M     *       0       0       TGCAGGTTCCGGGGGAGGAGGTGGGTACTCACGTTTTTTGGGGGGCCGTGGGGGAGCAGTGTAAGTGATGACTGAGGCGGCTTGGGCCCCCCCTGTGCTGGCCGGCCCACTCCCTGCTCCACCACTGCTGCTGCCCCCATGTACAATGACTGAGGCAGCCGGGGGGGCAAAAG     *       NM:i:7  MD:Z:87T63G7T9  AS:i:126        XS:i:22

# Hmmm, a few insertions in the exon, which is odd, but hopefully not too indicative of the dataset
	14:55312643-55312930_ENSOARE00000135137 16      utg5374_len=444093_reads=1434_status=X_microHet=1.00_covStat=1.00_quiver        347153  60      93M1I87M1I106M    *       0       0       AGGACACATTGAGCGTGACGGTCCTCTCCGCGCTCACAGAAGCCCTGCGGAAGGTCACCCGGCAGGTGAGGCTTGAGCCGTGATCCTGAGGGCGGGGGGTGAGCAGGATCTCCGAGGAGTTATAGGCCTCCTTGGAGTCCAGTCCCGGGGGTCTGAGGGCAGCCCCGGTCCAGGAGATGGTGGGGGGCGTGGCCCAGTCACAGGCCTCAGGCAGGGAGCATTTGAGGTGGCTGCTGGAGCCGGACTCCAGGGGCTCCTCGATGTGGACATCAGGGGTCCATGTTAGCG    *       NM:i:7  MD:Z:19A53C5A8C56A140   AS:i:247        XS:i:23

# Similar pattern. Maps but a few insertions
# Let's be fair and test the opposite situation (maps in BGI but not ours)
samtools view oar_bgi_exon_align.bam | grep 'ENSOARE00000187526'
	3:858498-858734_ENSOARE00000187526      0       NC_022303.1     103905526       60      195S42M *       0       0       GCGGGGCGCGGCGGGCTCCGCTGGGGGCGCCTCCGGGCGGGCCTCCATGGCCACCCCCTCCAGCGCCTCCAGGAAATCGAAGCTTCGACAGTGGCGCTTATTGAAAGGCGCGGGCACTGCAGGCCGGGCCATGGGGCCCTCACATGGCGGGGGCCTCGACACGGGCCCCGGGGCCGCCACCTTGATGTCCTGGTACACGGTCGACACCAGCAGATCCGGCGGGTCGGTGCGGGTCAT       *       NM:i:0  MD:Z:42 AS:i:42 XS:i:0

# With a cigar score like that, you can hardly call it mapping!
samtools view oar_bgi_exon_align.bam | grep 'ENSOARE00000053917'
	15:27555669-27555698_ENSOARE00000053917 0       NC_022307.1     26339931        60      30M     *       0       0       AACCAGACTAAAGAGGTATCCTTCTCTCTC    *       NM:i:0  MD:Z:30 AS:i:30 XS:i:0
samtools view oar_bgi_exon_align.bam | grep 'ENSOARE00000130790'
	X:14481160-14481307_ENSOARE00000130790  16      NC_022322.1     114826454       60      148M    *       0       0       CACATTTAATACTCTCTATCCGCACAGGGAGGCCAGACTGTGCCGCAGCAATTAATTTCAGGGCAACGTAAAACTGTGTCAGACCAAAATAACCAACCCGCTTGGCACCACACAATTCTGTGATCTGGAAACAGACAAGAGAAGGGAT      *       NM:i:1  MD:Z:127A20       AS:i:143        XS:i:20

# Those aligned fair and square -- I wonder if the regions are on the deg contig list of Sergey's assembly?
# Saving the lists
mv group_1.txt venn_comparisons/bgi_exclusive_unmapped_exons.txt
mv group_2.txt venn_comparisons/pacbio_exclusive_unmapped_exons.txt
mv group_1_2.txt venn_comparisons/shared_exclusive_unmapped_exons.txt

# Now for the split reads
samtools view oar_pacbio_exon_align.bam | grep 'ENSOARE00000105192'
	14:53127459-53129584_ENSOARE00000105192 0       utg3570_len=2632531_reads=9134_status=X_microHet=1.00_covStat=1.00_quiver       2300967 60      209M1I95M1I454M1I298M2D70M2D75M1I103M1I311M65I62M1I99M86D92M1D187M        *       0       0       (long sequence)        *       NM:i:191        MD:Z:78G406G13C19T54T4G53A194T144G13C68^GG33G31G4^GC36C129G42C4C28A30T77G18G9A5G96C17C1G145^GGGCATGCCTCCGCAGTGGAGAAGCCACTAAGGTACCCTGGGCGTGGGGCGTGGGGTGGAGGCCCGCCTGCCGCACAGTCCTCCGG92^T5T28A112T4G34       AS:i:1682       XS:i:0
# That shows an 86 bp insertion in the middle of the exon!
# Let's compare to the BGI alignments
samtools view oar_bgi_exon_align.bam | grep 'ENSOARE00000105192'
	14:53127459-53129584_ENSOARE00000105192 0       NC_022310.1     52501884        60      582S477M2D70M2D462M535S *       0       0       (long sequence)        *       NM:i:23 MD:Z:53A177A17T144G13C68^GG33G31G4^GC36C107C22G43C4C28A30T77G18G9A5G56C15 AS:i:898        XS:i:0  SA:Z:NC_022310.1,52500222,+,455M1671S,60,2;NC_022310.1,52503159,+,1847S92M1D187M,60,5;
	14:53127459-53129584_ENSOARE00000105192 256     NC_022310.1     52500222        60      455M1671H       *       0       0       AGTCAGCCACAACAACAAGGTGAACCTCATGACCAGCGAGAACCTGTCCATCTGCTTCTGGCCCACCCTGATGAGGCCAGACTTCAGCACCATGGACGCCCTCACGGCCACGCGCACCTACCAGACAATCATCGAGCTCTTCATCCAGCAGTGCCCCTTCTTCTTCCACAACCGGCCCATCAGCGAGCCCCCTGGTGCCACGCCCAGCTCCCCCTCTGCCCTGGCCTCCACCGTGCCCTTCCTCACCTCCACGCCCGTCACAAGTCAGCCGTCGCCGCCCCAGTCGCCTCCGCCCACCCCCCAGTCCCCAATGCAGCCGCTGCTCCCCTCCCAACTTCAAGCCGAACACACGCTGTGAGCCCCCGAGGCCTCGGGCTGAAGGAGACAAGGTCTTCTCTCTTGACGGGGTGGCACTGGGCCTTGAACAGAACCAGGCCTGCTGGGGTGGGGCGGGG*NM:i:2  MD:Z:78G32A343  AS:i:445        XS:i:119        SA:Z:NC_022310.1,52501884,+,582S477M2D70M2D462M535S,60,23;NC_022310.1,52503159,+,1847S92M1D187M,60,5;
	14:53127459-53129584_ENSOARE00000105192 256     NC_022310.1     52503159        60      1847H92M1D187M  *       0       0       GGGTCTCTTTGCACTTTAAGAAGCAGCACCTGTGGTATGTTCGAGGCCGCCTGACAAAGTGAGAGACGCCTCAGAACAGTCTCGGCCACCCCGGCCCCAGACGCTTCCCACTGGCTGTTGGAGGGCGGGCCAGGCCCTGTGCTCACCGACTCTTCCTCTGGATAAGCTGGATGTCGCACTGTCTCTAGTTCTTGTACAGCAAAATAAGGTGCTTGGAAAAACTACATGTCCCCTGCACGCTTCACCGCCTCCTTAAGGACACAAGTCAGAAAACCTAAG     *       NM:i:5  MD:Z:92^T5T28A112T4G34    AS:i:252        XS:i:0  SA:Z:NC_022310.1,52501884,+,582S477M2D70M2D462M535S,60,23;NC_022310.1,52500222,+,455M1671S,60,2;

# OK, in BGI, they have a 455 alignment at 52,500,222, a 1009 alignment at 52,501,884 and a 279 alignment at 52,503,159. 
# So the alignment is split three ways 
# NCBI blast assigns this to a Rho GTPase

# One more
'ENSOARE00000015203' # is split on the same BGI chr, about 1kb apart from its other half
# In the PacBio assembly, it is all contiguous

# Now, the opposite circumstance again (Pacbio split -- BGI whole)
'ENSOARE00000139505' # is split on the same PacBio chr, about 200 bp apart from its other half
'ENSOARE00000139505' # is unmapped on BGI

'ENSOARE00000068956' # is split twice, 200 bp apart again on the same PacBio contig
# It maps once to BGI, but it is split by hard clipping right in the middle

# Saving the files
mv group_1.txt venn_comparisons/bgi_exclusive_split_exons.txt
mv group_2.txt venn_comparisons/pacbio_exclusive_split_exons.txt
mv group_1_2.txt venn_comparisons/shared_exclusive_split_exons.txt
```

In general, it looks like the PacBio assembly recovers far more of the sheep exons than the BGI assembly could. Let's test the hypothesis that the PacBio assembly is better here by taking one of the clear competitive test cases : **ENSOARE00000105192**

The signature that I'm looking for here, is the presence of sufficient coverage of reads that have NO significant soft-clipped bases and/or discordant pairing profiles from the test animals that I downloaded from the SRA. I'll pull the fasta sequence from the BGI and PacBio genomes in the area, then align the reads to both and pull out the non-split reads.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/test_animals

```bash
# Downloading the SRA data I need
wget wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR405/ERR405776/ERR405776.sra
```

Now, to generate the test regions so that I can align the reads from the SRA and check fidelity.

| Genome | Coords |
| :--- | ---: |
| BGI | NC_022310.1:52,500,222-52,503,438 |
| PacBio | utg3570_len=2632531_reads=9134_status=X_microHet=1.00_covStat=1.00_quiver:2300967-2303093 | 

```bash
mkdir test_region
samtools faidx chir_bgi_pseudo_chr.fa NC_022310.1:52,500,222-52,503,438 > test_region/bgi_ENSOARE00000105192_region.fa
# NOTE: The BGI region has a huge gap in it!
samtools faidx USDA_V3.fasta utg3570_len=2632531_reads=9134_status=X_microHet=1.00_covStat=1.00_quiver:2300967-2303093 > pacbio_ENSOARE00000105192_region.fa

cd test_region/
bwa index bgi_ENSOARE00000105192_region.fa
bwa index pacbio_ENSOARE00000105192_region.fa
```

I'll have to wait for the SRA data to finish downloading. In the meantime, here's a summary table of comparative gene mappings (total exons: 234505):

| O. aries exons | Exclusive PacBio | Exclusive BGI | Shared |
| :--- | :--- | :--- | :--- |
| Total Mappings| 9,534 | 564 | 225,365 |
| Unmapped | 564 | 9,534 | 9,140 |
| Split Exons | 592 | 1,930 | 1,528 |

*5/22/2015*

--

OK, the files downloaded, its time to decompress them and then start the alignment. 

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/test_animals

```bash
~/sratoolkit.2.4.5-2-ubuntu64/bin/fastq-dump ERR405776.sra

# This file is pretty damn big, but I'm not sure if I want to split it
rm ERR405776.sra
```

OK, this will take a long time, but hopefully I can run the unmapped read filtering simultaneously here.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/liftover_comp/test_region

```bash
# Running with 5 threads and piping to a samtools process that will filter unmapped reads by default
bwa mem -M -t 5 bgi_ENSOARE00000105192_region.fa ../../test_animals/ERR405776.fastq | samtools view -bS -F 4 - > bgi_ERR405776.bam

bwa mem -M -t 5 pacbio_ENSOARE00000105192_region.fa ../../test_animals/ERR405776.fastq | samtools view -bS -F 4 - > pacbio_ERR405776.bam

samtools sort bgi_ERR405776.bam bgi_ERR405776.sorted
samtools sort pacbio_ERR405776.bam pacbio_ERR405776.sorted

samtools index pacbio_ERR405776.sorted.bam
samtools index bgi_ERR405776.sorted.bam
```
