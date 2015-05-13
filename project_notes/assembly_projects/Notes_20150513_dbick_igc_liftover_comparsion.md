# Goat assembly Liftover
---
*5/13/2015*

These are my notes on processing Ben's liftover coordinates for the Goat assembly. I will try to resolve some of the issues that he has with the files by scripting custom commands to prevent loss of data from GFF to BED and vice versa.

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