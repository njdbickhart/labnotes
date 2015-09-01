# Identifying Immune gene region differences between assemblies
---
*8/31/2015*

I need to see if overly complex regions of the genome are resolved better through the use of PacBio + lachesis + BNG information.

I just learned from Ben that the INDELs that we corrected via PILON completely borked his gene predictions. He'll have to recall the genes using his pipeline from the beginning. 

## Table of Contents
* [Summary statistics](#stats)

<a name="stats"></a>
## Summary statistics

*9/1/2015*

I'm going to compile some statistics from the pre-corrected and post-corrected assemblies. Here is the location of Steve's PILON corrected data:

> Blade14: /mnt/nfs/nfs2/GoatData/Pilon

Here is a decription of the contents of the directory according to Steve's README:

> **Goat333HybScaffolds1242contigs0723** - the BNG scaffolded assembly (based on Papadum-v4 : split contigs from v3 - no degenerates)
> 
>**papadum-v3-pilon** - Pilon corrected contigs (utg's only) from Papadum-v3 (no degenerates)
>
>**papadum-v4s-pilon** - Pilon corrected everything (utg's & dtg's) based on papadum-v3 with split contigs as predicted by Irys scaffolds
>
>**papadum-BNG-pilon** - Pilon corrected BNG scaffolding - no degenerates


I'm going to generate the following statistics for each file:
* Number of N gaps
* Total count of N's
* Total number of scaffolds
* Total count of non-N bases

> Blade14: /mnt/nfs/nfs2/GoatData/Pilon

```bash
# NOTE: papadum-BNG-pilon was not ready by the time I started this analysis
# Getting the number of contigs
wc -l *.fai
   1575 Goat333HybScaffolds1242contigs0723.fa.fai
   3074 papadum-v3-pilon.fa.fai
  33803 papadum-v4s-pilon.fa.fai

# Now processing gap regions
~/jdk1.8.0_05/bin/java -jar ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -f Goat333HybScaffolds1242contigs0723.fa -o Goat333HybScaffolds1242contigs0723.gaps.bed -s Goat333HybScaffolds1242contigs0723.gaps.stats

~/jdk1.8.0_05/bin/java -jar ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -f papadum-v3-pilon.fa -o papadum-v3-pilon.gaps.bed -s papadum-v3-pilon.gaps.stats

~/jdk1.8.0_05/bin/java -jar ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -f papadum-v4s-pilon.fa -o papadum-v4s-pilon.gaps.bed -s papadum-v4s-pilon.gaps.stats

# OK, I noticed a discrepancy between utg0 in the 4s and 3 file. I'm going to check it out using Mummer.
```

#### Side project: check reference fasta fidelity

It's always bugged me that they found discrepancies between my split fast and the v3 assembly fasta. Let's see what's going on here from an alignment standpoint.

```bash
# extracting the reference and query fasta
samtools faidx papadum-v3-pilon.fa utg0_pilon:1-2942271 > /mnt/iscsi/vnx_gliu_7/goat_assembly/v3_pilon_test.fa
samtools faidx papadum-v4s-pilon.fa utg0:1-2942312 > /mnt/iscsi/vnx_gliu_7/goat_assembly/v4_pilon_test.fa
```

Now I'm going to migrate to my ISCSI share to work on the nucmer aligns

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/

```bash
# Running nucmer
~/MUMmer3.23/nucmer -mumref -l 100 -c 1000 -p v3_vs_v4_test v3_pilon_test.fa v4_pilon_test.fa
~/MUMmer3.23/delta-filter -i 95 -o 95 v3_vs_v4_test.delta > v3_vs_v4_test.best.delta
~/MUMmer3.23/dnadiff -p v3_vs_v4_test. -d v3_vs_v4_test.best.delta

# Making the image
# Note: the fasta has to be the query
~/MUMmer3.23/mummerplot -p v3_vs_v4_test --large --fat --png -Q v4_pilon_test.fa v3_vs_v4_test..1delta
cat v3_vs_v4_test.gp |awk '{if (match($0, "len=")) { print substr($1, 1, index($1, "_")-1)"\" "$2" "$3} else print $0}'  > v3_vs_v4_test.fixed.gp

# I had to change the terminal to "emf" and comment ('#') out the mouse settings
gnuplot -geometry 900x900+0+0 -title mummerplot v3_vs_v4_test.fixed.gp
```

There were no major changes detected... the ends are the same. I'm not sure why the middle was different.

Time to brute force the issue with BWA:

```bash
bwa mem v3_pilon_test.fa v4_pilon_test.fa > v3_vs_v4_test.sam

# Now to process the cigar string
perl -lane 'if($F[0] =~ /^@/){next;}else{%h; while($F[5] =~ /(\d+)([MID])/g){$h{$2} += $1;} foreach my $k (keys(%h)){print "$k:\t$h{$k}";}}' < v3_vs_v4_test.sam
M:      2942230
I:      82
D:      41
```

Hmm... there are twice as many Insertions as deletions, could Pilon be randomly correcting alignments here? 

MD:Z tag:
> MD:Z:33993^G21G6740T211G226T11G8454C162T11G25A12T67G2T0G7^A16C41A29G27G13T30T16G28C0T12C27T39G11^A10117G87G22945^G288^G971091^G147^C15643A14701^C92690A191T38G31653A105G59A241^C1808^C189^G62017A7G61T4956A9G96G35C259G91205^G261^G40447^C401^C58324^A508^G38343G126^C55592G51T189^G128047^G388^C89585A14C26A75A15G2C64717C0G248823C19385G245A35064A178098^G246835A4C232^C26T6C1^CACAT8T13G27566T139G83^T68989^C30G224^G26245G78^C14863^TT3T0T281^T22172A9C154^T20796G34A34028T379A1621^T4058G267^G8313G212C20221^T56C9032^G414^G83321T307A16683G14A99A4735

Cigar string:
> 33993M1D7245M1I8624M2I93M1D241M1I8M1I52M1D10076M1I465M1I490M1I382M1I19825M1I338M1I1575M1D288M1D231993M1I577M1I424250M1I334M1I31393
7M1D147M1D30345M1D10M5I92747M1I32226M1D1808M1D189M1D66931M1I91722M1D261M1D40447M1D401M1D58324M1D5
08M1D38470M1D55834M1D128047M1D388M1D27063M1I337M1I362730M1I539M1I12840M1I54085M1I478M1I177822M1I1
67M1D66769M1I431M1I25595M1I567M1I152080M1I1631M1D3M6I4M1I4M1I7M1I4M1I4M1I5M13I4M5D10M1I8M5I27795M
1D68989M1D255M1D26324M1D14863M2D286M1D22337M1D56862M1D4326M1D8219M1I146M14I338M1I20045M1D9089M1D4
14M1D62102M1I43062M 

There are allot of SNPs too. Let's call variants here so that we have an easily indexable version:

```bash
samtools view -bS v3_vs_v4_test.sam > v3_vs_v4_test.bam
samtools index v3_vs_v4_test.bam

samtools mpileup -ugf v3_pilon_test.fa v3_vs_v4_test.bam | bcftools call -vmO v -o v3_vs_v4_test.vcf

grep 'INDEL' v3_vs_v4_test.vcf | wc -l
22

grep -v '#' v3_vs_v4_test.vcf | wc -l
105  <- 105 - 22 = 83 SNPs
```

So there are 83 SNPs and 22 INDELs compared between the two contigs after PILON correction. This I have to check with my original fasta!

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/contig_check

```bash
grep 'utg0' /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa.fai
utg0    2936728 6       60      61

grep 'utg0' /mnt/nfs/nfs2/GoatData/BioNano_conflict_res/goat_split_36ctg_assembly.fa.fai
utg0    2936728 250176000       60      61

# they look the same... let's try aligning them quickly
samtools faidx /mnt/nfs/nfs2/GoatData/BioNano_conflict_res/goat_split_36ctg_assembly.fa utg0:1-2936728 > split_36_ctg_utg0test.fa
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v3/papadum-v3s.ctg.fa utg0:1-2936728 > original_ctg_utg0test.fa

bwa index original_ctg_utg0test.fa
bwa mem original_ctg_utg0test.fa split_36_ctg_utg0test.fa > split_original_test.sam
```

Cigar string:
> 2936728M

Identical! So that means that the inclusion of the degenerate contigs in Steve's papadumv4 most likely impacted the correction of the contig. So, it's not a deterministic correction unless you have the same input and output data.

#### End side project

To recap:
* Both the split contig fasta and the original v3 assembly fasta had the same exact sequence for utg0
* The v3 pilon correction without the degenerate contigs had 82 fewer deleted bases and 41 more inserted bases than the v4 pilon corrected assembly
* The only notable differences in these two assemblies for utg0 is the fact that v4 had degenerate contigs present.

