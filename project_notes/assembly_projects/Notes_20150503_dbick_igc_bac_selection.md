# Selecting BAC clones for IGC grant
---
*05/03/2015*

### IGC BAC selection coordinate refinement

George Liu remembers that many of the BAC plates for the RPCI42 and CHORI-240 cell lines were dropped in the move. I will be selecting BAC clones that cover three regions that contain the IGCs of interest and reordering them from CHORI bacpac:

| IGC type | UMD3 region | Size |
| :---: | :---: | :---: |
| MHC I | chr23:28,500,000-28,750,000 | 250 kb |
| LRC | chr18:63,100,000-63,400,000 | 300 kb |
| NKC | chr5:99,500,000-99,850,000 | 350 kb |
| MHC refined| chr23:28,300,000-28,750,000 | 450 kb |

These coordinates were approximated from the maps that John Hammond provided, and are likely to be inexact. I'm going to check the UCSC genome browser for each of these coordinate windows just to make sure that the genes I'm interested in are actually there.

  * [MHC I UCSC Link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=bosTau6&position=chr23%3A28300000-28750000&hgsid=424706521_EdAAYfF9G1BVcOsr6Jlvyha9s30M)
  * [LRC UCSC Link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=bosTau6&position=chr18%3A63100000-63400000&hgsid=424706521_EdAAYfF9G1BVcOsr6Jlvyha9s30M)
  * [NKC UCSC Link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=bosTau6&position=chr5%3A99500000-99850000&hgsid=424706521_EdAAYfF9G1BVcOsr6Jlvyha9s30M)

I enlarged the search window for the MHC class I receptor because the BOLA refseq annotation extended out past the initial region that was highlighted in John's figures. 

### BAC clone selection 

I have downloaded the RPCI BAC end sequences from NCBI, but I need to get the approximate coordinates for each pairing. I'm going to try using BWA MEM as a first-pass alignment -- I know it isn't going to be the best alignment, but it's fast and automated.

Actually, it looks like I've already done most of this work in one of my previous labnotes files (./textpad_archive/Lab_notes_20141009_dbick_john_chr18_bac_selection.note). I'm going to transfer most of the commands over to this document for better notekeeping.

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing/RP42

```bash
# This script turned each Sanger read into a pseudo fastq read
perl process_fasta_to_fastq_bac_libraries.pl ncbi_genbank_full_bac_end.fasta
bwa mem -t 2 /home/dbickhart/share/umd3_data/umd3_kary_unmask_ngap.fa ncbi_genbank_full_bac_end.fasta.fq > ncbi_genbank_full_bac_end.sam
samtools view -bS ncbi_genbank_full_bac_end.sam | samtools sort - ncbi_genbank_full_bac_end_sorted
samtools index ncbi_genbank_full_bac_end_sorted.bam

# I searched NCBI using the following query: GCC (in the dropbox) (RPCI-42) AND Bos taurus[Organism] 
perl process_fasta_to_fastq_bac_libraries.pl ncbi_everything_gss_bac.fasta
bwa mem -t 2 /home/dbickhart/share/umd3_data/umd3_kary_unmask_ngap.fa ncbi_everything_gss_bac.fasta.fq > ncbi_everything_genbank_bac_end.sam
samtools view -bS ncbi_everything_genbank_bac_end.sam | samtools sort - ncbi_everything_genbank_bac_end_sorted
samtools index ncbi_everything_genbank_bac_end_sorted.bam

# This final bam file (ncbi_everything_genbank_bac_end_sorted.bam) should contain all of the RPCI-42 BAC end sequence traces in NCBI currently
```

I could improve this by making it a "pseudo-paired-end" fastq, but this is fine for now. I will manually check overlapping BAC ends to ensure that they're mapping to the same region.

##### RPCI42 MHC class I region

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing/RP42

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr23:28,300,000-28,750,000 | wc -l
# 20 BAC ends in the region
# Let's reformat this so that it is more easily parseable
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr23:28,300,000-28,750,000 | perl -lane 'print "$F[2]\t$F[3]\t$F[0]";'
```
> chr23	28321536	RPCI42_133J13.TV <- If there wasn't an issue with the other end, this would be perfect
chr23	28327765	RP42-47B01
chr23	28328466	RPCI42_113K1.TV <- good candidate
chr23	28343013	RPCI42_148M10.TJ
chr23	28343014	RPCI42_148O8.TJ		<- good candidate
chr23	28352169	RPCI42_101L23.TV
chr23	28423747	RPCI42_160H4.TV
chr23	28444302	RPCI42_154E17.TV
chr23	28447844	RPCI42_133J13.TJ	<- Not sure why this is the same as the next one
chr23	28448607	RPCI42_133J13.TJ
chr23	28452390	RPCI42_113K1.TJ
chr23	28465336	RPCI42_122O1.TJ <- good candidate
chr23	28548725	RPCI42_148O8.TV
chr23	28548735	RPCI42_148M10.TV
chr23	28639463	RP42-164F10
chr23	28657595	RPCI42_122O1.TV
chr23	28725589	RPCI42_113O16.TJ
chr23	28726368	RPCI42_154E4.TV
chr23	28740559	RPCI42_135I11.TJ
chr23	28749847	RP42-14B12

Let's check RPCI42_133J13.TJ to see what's going on real quick.

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr23:28,300,000-28,750,000 | grep 'RPCI42_133J13.TJ'
```
>RPCI42_133J13.TJ	16	chr23	28447844	60	71M1I14M1I520M**107S**	*	0	0 {...}	NM:i:41	AS:i:396	XS:i:57	XP:Z:chr23,-28448607,609S105M,60,6;
RPCI42_133J13.TJ	16	chr23	28448607	60	609S**105M**	*	0	0	{...}	NM:i:6	AS:i:79	XS:i:33	XP:Z:chr23,-28447844,71M1I14M1I520M107S,60,41;

OK, the CIGAR scores indicate that this is a split read on the assembly! Really interesting! So, this is a definite keeper. 

I'm worried that I'm not capturing the last half of the region sufficiently with the current paired BAC clones. Let's expand the region and see if some of the clones near the 28.6 MB region fit in with my criteria.

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr23:28,600,000-28,950,000 | perl -lane 'print "$F[2]\t$F[3]\t$F[0]";'
```
> chr23	28639463	RP42-164F10
chr23	28657595	RPCI42_122O1.TV
chr23	28725589	RPCI42_113O16.TJ <- This one pairs up now
chr23	28726368	RPCI42_154E4.TV
chr23	28740559	RPCI42_135I11.TJ
chr23	28749847	RP42-14B12
chr23	28874325	RP42-164F10	<- This one pairs up with the previous results now!
chr23	28886699	RPCI42_113O16.TV
chr23	28916246	RP42-23N10
chr23	28931479	RPCI42_135I11.TV

OK, let' summarize for MHC:

| Clone Name | End Coords | Notes |
| :--- | :---: | :--- |
| RPCI42_133J13 | chr23:28,321,536-28,448,607 | The 3' end is a split alignment on this assembly |
| RPCI42_113K1 | chr23:28,328,466-28,452,390 | |
| RPCI42_148O8 | chr23:28,343,014-28,548,725 | |
| RPCI42_122O1 | chr23:28,465,336-28,657,595 | |
| RP42-164F10 | chr23:28,639,463-28,874,325 | This is probably the smallest overlap region |
| RPCI42_113O16 | chr23:28,725,589-28,886,699 | This extends 100kb beyond the region |

##### RPCI42 LRC Region

OK, this should be the same drill as with the MHC class I region.

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing/RP42

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr18:63,100,000-63,400,000 | wc -l
# 1 clone in the region

# Oh boy, I was fearing something like this. Let's extend the search by 150kb in both directions
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr18:62,950,000-63,550,000 | wc -l
# 6 clones now
# OK, let's get the info to see what we have
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr18:62,950,000-63,550,000 | perl -lane 'print "$F[2]\t$F[3]\t$F[0]";'
```
> chr18	63005919	RPCI42_145P10.TV
chr18	63288030	RP42-168O11	<- this is a natural
chr18	63408180	RP42-168O11
chr18	63478958	RP42-165K18
chr18	63538092	RP42-41H04
chr18	63539979	RPCI42_148B8.TV

There are some clones that *might* be in the same region, but it's looking pretty dire. Let's grep out some of the possible contenders. When I pulled **145P10**, the ".TS" end sequence mapped to chr26 originally. When I blasted it against the NCBI nr database; however, the top hits were for LILR genes and pseudogenes. So, this might be an anchor! **165K18** mapped further downstream (chr18:63594540). 

OK, so only two clones for this one? That's a shame, but we'll press on.

| Clone Name | End Coords | Notes |
| :--- | :---: | :--- |
| RPCI42_145P10 | chr18:63,005,919-? | The end sequence might not currently be on this assembly |
| RP42-168O11 | chr18:63288030-63408180 | |

##### RPCI42 NKC Region

OK, third time's a charm.

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing/RP42

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr5:99,500,000-99,850,000 | wc -l
# 3 clones in the region

# Sometimes I hate my job
# Let's expand the region to see if that gets us more to work with
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr5:99,450,000-100,000,000 | wc -l
# 7 clones in the region

samtools view ncbi_everything_genbank_bac_end_sorted.bam chr5:99,450,000-100,000,000 | perl -lane 'print "$F[2]\t$F[3]\t$F[0]";'
```
> chr5	99665696	RP42-161F13 <- good one
chr5	99735207	RPCI42_154D6.TV <- good one
chr5	99780339	RP42-161F13
chr5	99866305	RP42-162P15 <- good one
chr5	99922628	RPCI42_157F9.TV
chr5	99968671	RPCI42_154D6.TJ
chr5	99991606	RP42-162P15

We've got 3 good ones, but nothing that covers the 99.5 Mb region. I'm expanding the search in the 5' direction to try to find any possible leads.

```bash
samtools view ncbi_everything_genbank_bac_end_sorted.bam chr5:99,250,000-99,600,000 | perl -lane 'print "$F[2]\t$F[3]\t$F[0]";'
```
>chr5	99319620	RP42-15O08
chr5	99336873	RPCI42_156K1.TJ

Damn, I'm going to have to grep to see if these two are potential anchors. **RP42-15O08** mapped downstream (chr5:99,121,246). **RPCI42_156K1** also mapped downstream (chr5:99,197,978). So, we're missing coverage from chr5:99,450,000-99,665,969 for NKC right now. Maybe we can create probes for BAC selection down the line instead of ramming our heads against putative assembly issues.

So, to summarize the NKC region:

| Clone Name | End Coords | Notes |
| :--- | :---: | :--- |
| RP42-161F13 | chr5:99,665,696-99,780,339 | Shorter than the average BAC... Also earliest one for this region |
| RPCI42_154D6 | chr5:99,735,207-99,968,671 | |
| RP42-162P15 | chr5:99,866,305-99,991,606 | This one is 10kb upstream of the region, but it could anchor |

##### John's chr18 region

Just for my notes:

> tag SNP is ARS-BFGL-NGS-109285 at 57,589,121 bp on BTA 18

I had already selected BACs for John's project by using similar methods. Here is a table of the clones to use:

| Clone Name | End Coords | Notes |
| :--- | :---: | :--- |
| RPCI42_154M1 | chr18:57,371,161-57,531,510 | Flanking 5' region|
| RPCI42_118F24 | chr18:57,548,063-57,704,285 | Overlaps SNP |
| RPCI42_118B22 | chr18:57,548,064-57,704,269 | Perhaps a clone of the previous region? |
| RPCI42_3D15 | chr18:57,657,380-57,806,510 | Flanking 3' region|

So we've got two clones that overlap the SNP, and two clones that flank it on either side. Pretty good to start with.

--

### RPCI-42 BAC selection summary table and notes

Here are the full list of clones that we can pull for RPCI42 for this region and project

| Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :---: | :--- |
| MHC | RPCI42_133J13 | chr23:28,321,536-28,448,607 | The 3' end is a split alignment on this assembly |
| MHC | RPCI42_113K1 | chr23:28,328,466-28,452,390 | |
| MHC | RPCI42_148O8 | chr23:28,343,014-28,548,725 | |
| MHC | RPCI42_122O1 | chr23:28,465,336-28,657,595 | |
| MHC | RP42-164F10 | chr23:28,639,463-28,874,325 | This is probably the smallest overlap region |
| MHC | RPCI42_113O16 | chr23:28,725,589-28,886,699 | This extends 100kb beyond the region |
| LRC | RPCI42_145P10 | chr18:63,005,919-? | The end sequence might not currently be on this assembly |
| LRC | RP42-168O11 | chr18:63288030-63408180 | |
| NKC | RP42-161F13 | chr5:99,665,696-99,780,339 | Shorter than the average BAC... Also earliest one for this region |
| NKC | RPCI42_154D6 | chr5:99,735,207-99,968,671 | |
| NKC | RP42-162P15 | chr5:99,866,305-99,991,606 | This one is 10kb upstream of the region, but it could anchor |
| chr18 | RPCI42_154M1 | chr18:57,371,161-57,531,510 | Flanking 5' region|
| chr18 | RPCI42_118F24 | chr18:57,548,063-57,704,285 | Overlaps SNP |
| chr18 | RPCI42_118B22 | chr18:57,548,064-57,704,269 | Perhaps a clone of the previous region? |
| chr18 | RPCI42_3D15 | chr18:57,657,380-57,806,510 | Flanking 3' region|

I'm missing a BAC that covers chr5:99,450,000-99,665,969 for NKC, and the LRC region is indeterminately covered (one BAC has only a partial mapping to the region!). Not the best start, but it's a start and we have 3 years at least. John will get 4 BACs that overlap or flank the tag SNP -- we can start from there.

--

*05/04/2015*

### CHORI BAC selection for IGC

I have already selected John's chr18 BACs, so now I need to identify the BACs that I need for John Hammond's project. Thankfully, these BACs from the Domino library are incredibly well characterized. So, let's work with the bed file containing their locations on UMD3 that I prepared previously.

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing

Actually, windows remote desktop is far too slow for me to work on my linux virtualbox. Instead, I will use a modification of my java BedTools package to do the intersection.

| Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :---: | :--- |
| MHC | CH240-503O23 | chr23:28,268,254-28,468,033 | |
| MHC | CH240-264A9 | chr23:28,631,696-28,807,313 | |
| MHC | CH240-203K5 | chr23:28,672,885-28,900,932 | |
| LRC | CH240-392A13 | chr18:63,415,417-63,594,513 | The only BAC nearby! |
| NKC | CH240-239G9 | chr5:99,854,599-100,001,447 | Another flanking BAC! |


Damn, not that many. I'm going to have to dig to find some of these. I'll use liftover with far more permissive settings in order to find these locations.

##### Working with difficult BACs for IGC regions

> pwd: /home/dbickhart/share/side_projects/john_sequencing/bac_sequencing

```bash
# First, 80% minimum match
../../../umd3_data/liftover/liftOver -minMatch=0.80 bac_placed_coord_sorted.bed ../../../umd3_data/liftover/bosTau4ToBosTau6.over.chain bac_permissive_umd3.bed bac_permissive_umd3.unmapped

wc -l *.bed
  19189 bac_full_list_ambiguous_sorted.bed
  11123 bac_full_list_ambiguous_sorted_umd3.bed	<- old file
    132 bac_oneandunknown_bacsorted.bed
  14727 bac_permissive_umd3.bed	<- about 3600 more clones mapping
  18626 bac_placed_coord_sorted.bed
    431 bac_unknown_bacsorted.bed
```

Still nothing. OK, the problem is that these regions are not getting proper BAC coverage because they're highly repetitive. Instead, let's try to back-liftover the regions to btau4 and select mapped clones from the list directly. 

UCSC only lifts over from UMD3 -> baylor 7 -> baylor 4. I had to set the list to multiple hits, and give it a "named bed" of the original coordinates. I then took only the coordinates that mapped on placed chromosomes in the region. Minimum query hit of 1000bp worked for LRC. 

> chr18	62990861	62991955	LRC	1

| Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :---: | :--- |
| LRC | CH240-370M3 | chr18:62,794,084-63,033,719 | Only ID'ed from permissive liftover. Btau4 coords |

		
> Here are the UMD3 coordinates for the 1kb segment
> chr18	63,129,043	63,130,130

OK, we're going to have to get even more permissive with NKC. I tried a 500bp minimum match, but this is the best I could get:

> chr5	107222161	107223401	NKC	1
chr5	107211670	107212619	NKC	1

This translates to the 100 Mb region on chromosome 5. Good enough for a flanking BAC, but not good enough to cover the whole region. Maybe I can work with this though, if I get another flanking BAC clone on the other side. I went to the browser and selected the leftmost region before a large gap in sequence: **chr5:99,493,589-99,513,114**. I'll try to lift over this region instead.

OK! That worked! The new coordinates are: 
**Btau4: chr5:106356067-106375125** 

Note, I shifted the coordinates over by 50kb on both sides to get the next bac clone results:

| Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :---: | :--- |
| NKC | CH240-49I22	| chr5:106,206,907-106,445,457 | Btau4 coords |
| NKC | CH240-60G5 | chr5:106,309,704-106,537,019 | Btau4 coords |
| NKC | CH240-274P11 | chr5:106,412,271-106,652,564 | Btau4 coords |

Now, I think that I have a good set of BACs to sequence. Let's list John's remaining BACs and then I'll create a final listing of coordinates and positions.

##### John's BACs (chosen by liftOver)

I had to liftOver John's SNP coordinate to find overlapping BACs (the region is far too shuffled for straight liftOver). The SNP coordinates on Btau4:

> chr18	57082908	57082909

| Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :---: | :--- |
| chr18 |  CH240-389P14 | chr18:56,954,654-57,129,335 | Btau4 coords |
| chr18 | CH240-234E12 | chr18:57,058,248-57,236,865 | Btau4 coords |
| chr18 | CH240-280L6 | chr18:57,092,237-57,268,067 | Btau4 coords |
| chr18 | CH240-34N7 | chr18:57,129,383-57,288,223 | Btau4 coords; Flanking 5' end |

OK, now I've got everything. Time to get the master table together.

---

# Final table of selected BACs for sequencing

|Library | Region Name | Clone Name | End Coords | Notes |
| :--- | :--- | :--- | :---: | :--- |
| RPCI-42 | MHC | RPCI42_133J13 | chr23:28,321,536-28,448,607 | The 3' end is a split alignment on this assembly |
| RPCI-42 | MHC | RPCI42_113K1 | chr23:28,328,466-28,452,390 | |
| RPCI-42 | MHC | RPCI42_148O8 | chr23:28,343,014-28,548,725 | |
| RPCI-42 | MHC | RPCI42_122O1 | chr23:28,465,336-28,657,595 | |
| RPCI-42 | MHC | RP42-164F10 | chr23:28,639,463-28,874,325 | This is probably the smallest overlap region |
| RPCI-42 | MHC | RPCI42_113O16 | chr23:28,725,589-28,886,699 | This extends 100kb beyond the region |
| RPCI-42 | LRC | RPCI42_145P10 | chr18:63,005,919-? | The end sequence might not currently be on this assembly |
| RPCI-42 | LRC | RP42-168O11 | chr18:63288030-63408180 | |
| RPCI-42 | NKC | RP42-161F13 | chr5:99,665,696-99,780,339 | Shorter than the average BAC... Also earliest one for this region |
| RPCI-42 | NKC | RPCI42_154D6 | chr5:99,735,207-99,968,671 | |
| RPCI-42 | NKC | RP42-162P15 | chr5:99,866,305-99,991,606 | This one is 10kb upstream of the region, but it could anchor |
| RPCI-42 | chr18 | RPCI42_154M1 | chr18:57,371,161-57,531,510 | Flanking 5' region|
| RPCI-42 | chr18 | RPCI42_118F24 | chr18:57,548,063-57,704,285 | Overlaps SNP |
| RPCI-42 | chr18 | RPCI42_118B22 | chr18:57,548,064-57,704,269 | Perhaps a clone of the previous region? |
| RPCI-42 | chr18 | RPCI42_3D15 | chr18:57,657,380-57,806,510 | Flanking 3' region|
| CHORI-240| MHC | CH240-503O23 | chr23:28,268,254-28,468,033 | |
| CHORI-240| MHC | CH240-264A9 | chr23:28,631,696-28,807,313 | |
| CHORI-240| MHC | CH240-203K5 | chr23:28,672,885-28,900,932 | |
| CHORI-240| LRC | CH240-392A13 | chr18:63,415,417-63,594,513 | The only BAC nearby! |
| CHORI-240| LRC | CH240-370M3 | chr18:62,794,084-63,033,719 | Only ID'ed from permissive liftover. Btau4 coords |
| CHORI-240| NKC | CH240-239G9 | chr5:99,854,599-100,001,447 | Another flanking BAC! |
| CHORI-240| NKC | CH240-49I22	| chr5:106,206,907-106,445,457 | Btau4 coords |
| CHORI-240| NKC | CH240-60G5 | chr5:106,309,704-106,537,019 | Btau4 coords |
| CHORI-240| NKC | CH240-274P11 | chr5:106,412,271-106,652,564 | Btau4 coords |
| CHORI-240| chr18 |  CH240-389P14 | chr18:56,954,654-57,129,335 | Btau4 coords |
| CHORI-240| chr18 | CH240-234E12 | chr18:57,058,248-57,236,865 | Btau4 coords |
| CHORI-240| chr18 | CH240-280L6 | chr18:57,092,237-57,268,067 | Btau4 coords |
| CHORI-240| chr18 | CH240-34N7 | chr18:57,129,383-57,288,223 | Btau4 coords; Flanking 5' end |

28 Total BACs in this list. The LRC region is woefully undercovered. Let me ask John Hammond for suggestions/advice. Perhaps he's found an unplaced contig that has LRC genes? I will also ask him for probe sequence to see if I can map the probe sequence to new coordinates on several assemblies to try to resolve this region.

-- 

*05/05/2015*

# Increasing LRC BAC coverage

John Hammond sent me a list of LRC BACs that the pirbright team selected for sequencing, previously. MHC Class I is appropriately covered in this dataset, but we could use more LRC BACs. From [the excel spreadsheet](E:\SharedFolders\grants\immune_gene_cluster_grant\bac_selection\list_of_clones.xlsx) I am going to select clones that have end sequencing and that are on different chromosomes. 

| chromosome | clone |
| :---: | :---: |
|4	|CH240-312G3|
|4|	CH240-342K15|
4	|CH240-372D22
4	|CH240-371O5
4	|CH240-467K1
4	|CH240-386N21
10	|CH240-273N20
18	|CH240-373G10
18	|CH240-436K15
18	|CH240-399M16
18	|CH240-358F12
18	|CH240-269F7
18	|CH240-364C4
18	|CH240-379C3
X	|CH240-391K10
X	|CH240-315A1
X	|CH240-455F3
1	|CH240-387G13
1	|CH240-194J1
1	|CH240-330I21
1	|CH240-378O5
18	|CH240-448F6
18	|CH240-388F19
18	|CH240-237P10

I'm going to pull the coordinates for each clone and select a few that overlap for each chromosome segment. 


