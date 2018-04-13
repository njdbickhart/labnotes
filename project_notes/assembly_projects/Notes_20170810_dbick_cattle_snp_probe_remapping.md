# SNP probe remapping on cattle assemblies
---
*8/10/2017*

These are my notes for remapping UMD3-based SNP probes onto the new cattle assembly. My hope is to do as thorough of a job as possible and to resolve >99% of marker locations, if possible. 

## Table of contents
* [Preparing reference fastas and generating metadata](#prepwork)
	* [Coordinate changes](#coordchanges)
	* [HD probe remapping](#hdremapping)
* [Testing regions of the assembly for completion](#testcompletion)
* [Non-repeatmasked chain of assembly file](#nonrepeatmasked)
	* [Chaining v25 to UMD3](#chain25umd)
* [Mapping Bob's markers to v25 from NCBI](#mapbobsmarkers)

<a name="prepwork"></a>
## Preparing reference fastas and generating metadata

First, let's talk about what I need. I want to approach this from a hierarchical standpoint, where I:

* Identify rearranged segments of the genome using liftover
* Confirm simple rearrangements using a fast aligner (BWA MEM)
* Perform specialized local alignment using a more precise alignment tool to resolve recalcitrant probes

I suspect that 98% of markers will be resolved by the fast aligner, whereas the remaining 2% will require far more care. If needed, I may even need to go to the LD information on the marker to try to find regions of sequence to align to in a very targeted approach!

OK, first things first, I need to generate the requisite data for liftover. I am following the general guide, listed [here](http://genomewiki.ucsc.edu/index.php/LiftOver_Howto). In order to assist with this, I am going to use RepeatMasker to reduce the complexity of the genome (I need the repeat information for our subsequent analysis in any case, so this will be killing two birds with one stone). Reducing the complexity of the genome will make the blat alignment much faster by eliminating ambiguous alignments. 

I had to download a mirror of the compiled [kent utilities](https://github.com/ENCODE-DCC/kentUtils) in order to process lastz alignments. I found the sequence conservation parameters recommended for cow alignments [here](http://genomewiki.ucsc.edu/index.php/Hg19_conservation_lastz_parameters). I would like to note that the parameters are hg19 to cow... but given the improvement of the quality of the reference, perhaps these are appropriate?

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/RepeatMasker/RepeatMasker -pa 40 -q -species cow -no_is -dir ars_ucd_14_igc_rmask ARS-UCD1.0.14.clean.wIGCHaps.fasta

# Now to segregate all fasta entries into separate files for parallel alignment
cd ars_ucd_14_igc_rmask/
module load samtools
samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta.masked

mkdir masked_chr_fastas
for i in `cat ARS-UCD1.0.14.clean.wIGCHaps.fasta.masked.fai | cut -f1`; do echo $i; samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta.masked $i > masked_chr_fastas/${i}.fa ; done

# Now to index the masked UMD3.1 reference
module load blat/0.35
faToTwoBit /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa umd3_kary_unmask_ngap.2bit
blat umd3_kary_unmask_ngap.2bit /dev/null /dev/null -makeOoc=11.ooc -repMatch=1024

# Generating alignment files now
# Note: larger chromosomes are too large to use the -fastMap option, so I need to speed the process up by using the ooc file
mkdir psl
for i in masked_chr_fastas/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --partition=assemble1 --ntasks-per-node=1 --mem=5000 --wrap="blat umd3_kary_unmask_ngap.2bit $i -tileSize=11  -minIdentity=98 psl/${name}.psl -ooc=11.ooc -noHead -minScore=100"; done

# Well, it looks like the larger chromosomes and contigs all failed due to memory constraints
grep 'CANCELLED' *.out | cut -d':' -f1 | xargs -I{} cat {} | grep 'Query' | perl -ane 'print "$F[2].fa ";'; echo
for i in masked_chr_fastas/10.fa masked_chr_fastas/11.fa masked_chr_fastas/12.fa masked_chr_fastas/13.fa masked_chr_fastas/14.fa masked_chr_fastas/15.fa masked_chr_fastas/16.fa masked_chr_fastas/17.fa masked_chr_fastas/18.fa masked_chr_fastas/19.fa masked_chr_fastas/1.fa masked_chr_fastas/20.fa masked_chr_fastas/21.fa masked_chr_fastas/22.fa masked_chr_fastas/23.fa masked_chr_fastas/24.fa masked_chr_fastas/25.fa masked_chr_fastas/26.fa masked_chr_fastas/27.fa masked_chr_fastas/28.fa masked_chr_fastas/29.fa masked_chr_fastas/2.fa masked_chr_fastas/3.fa masked_chr_fastas/4.fa masked_chr_fastas/5.fa masked_chr_fastas/6.fa masked_chr_fastas/7.fa masked_chr_fastas/8.fa masked_chr_fastas/9.fa masked_chr_fastas/Super-Scaffold_1723_ScbfJmS_2085.fa masked_chr_fastas/X.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --partition=assemble1 --ntasks-per-node=1 --mem=50000 --wrap="blat umd3_kary_unmask_ngap.2bit $i -tileSize=11  -minIdentity=98 psl/${name}.psl -ooc=11.ooc -noHead -minScore=100"; done

# Man, it keeps on running out of alloted memory! I think that the query size is too large for blat
# Moving on to lastz; Nevermind! Got an error message telling me that I need to use individual fastas for the reference!
#mkdir lav; for i in masked_chr_fastas/10.fa masked_chr_fastas/11.fa masked_chr_fastas/12.fa masked_chr_fastas/13.fa masked_chr_fastas/14.fa masked_chr_fastas/15.fa masked_chr_fastas/16.fa masked_chr_fastas/17.fa masked_chr_fastas/18.fa masked_chr_fastas/19.fa masked_chr_fastas/1.fa masked_chr_fastas/20.fa masked_chr_fastas/21.fa masked_chr_fastas/22.fa masked_chr_fastas/23.fa masked_chr_fastas/24.fa masked_chr_fastas/25.fa masked_chr_fastas/26.fa masked_chr_fastas/27.fa masked_chr_fastas/28.fa masked_chr_fastas/29.fa masked_chr_fastas/2.fa masked_chr_fastas/3.fa masked_chr_fastas/4.fa masked_chr_fastas/5.fa masked_chr_fastas/6.fa masked_chr_fastas/7.fa masked_chr_fastas/8.fa masked_chr_fastas/9.fa masked_chr_fastas/X.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --partition=assemble1 --ntasks-per-node=1 --mem=25000 --wrap="lastz umd3_kary_unmask_ngap.2bit $i Y=9400 L=3000 K=3000 > lav/${name}.lav "; done

mkdir masked_chr_fastas/chr_sub_chunks
sbatch --partition=assemble1 faSplit_blat_align.sh umd3_kary_unmask_ngap.2bit masked_chr_fastas/1.fa
for i in masked_chr_fastas/10.fa masked_chr_fastas/11.fa masked_chr_fastas/12.fa masked_chr_fastas/13.fa masked_chr_fastas/14.fa masked_chr_fastas/15.fa masked_chr_fastas/16.fa masked_chr_fastas/17.fa masked_chr_fastas/18.fa masked_chr_fastas/19.fa masked_chr_fastas/20.fa masked_chr_fastas/21.fa masked_chr_fastas/22.fa masked_chr_fastas/23.fa masked_chr_fastas/24.fa masked_chr_fastas/25.fa masked_chr_fastas/26.fa masked_chr_fastas/27.fa masked_chr_fastas/28.fa masked_chr_fastas/29.fa masked_chr_fastas/2.fa masked_chr_fastas/3.fa masked_chr_fastas/4.fa masked_chr_fastas/5.fa masked_chr_fastas/6.fa masked_chr_fastas/7.fa masked_chr_fastas/8.fa masked_chr_fastas/1.fa masked_chr_fastas/9.fa masked_chr_fastas/X.fa; do sbatch --partition=assemble1 faSplit_blat_align.sh umd3_kary_unmask_ngap.2bit $i; done

# Because the OLD files are so numerous, I'm going to get rid of them
ls psl/chr_sub_chunks/OLD*.psl | xargs -I{} rm {}

# Now to turn the psl files into chains. I'll condense the chr_sub_chunks first
faToTwoBit ARS-UCD1.0.14.clean.wIGCHaps.fasta.masked ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit
sbatch --partition=assemble3 psl_merge.sh umd3_kary_unmask_ngap.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit

# I screwed up the chain sort algorithm and had to run this step in a separate loop
for i in `seq 1 29` X; do /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/psl/chr_sub_chunks/${i}_*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/psl/${i}_chain stdin; done

# It produced separate chromosome chain files for each major chromosome in that directory.
# Trying to merge them into a superset of chains
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/psl/*_chain/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/combined_chain stdin

# Now for the unscaffolded elements
for i in psl/*.psl; do base=`basename $i | cut -d'.' -f1`; echo $base; /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl $i umd3_kary_unmask_ngap.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit psl/${base}.chain; done
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/psl/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/unscaffolded_chain stdin

# Final merge
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/combined_chain/*.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/unscaffolded_chain/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/final_merged_chain stdin

# Preparing for the net step
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/twoBitInfo ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/twoBitInfo umd3_kary_unmask_ngap.2bit umd3_kary_unmask_ngap.2bit.info

mkdir net
cat final_merged_chain/*.chain > final_merged_chain/combined.umd3_to_ars-ucd.chain
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort final_merged_chain/combined.umd3_to_ars-ucd.chain final_merged_chain/combined.umd3_to_ars-ucd.sorted.chain

# Generating the net
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet final_merged_chain/combined.umd3_to_ars-ucd.sorted.chain umd3_kary_unmask_ngap.2bit.info ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info net/combined.umd3_to_ars-ucd.net /dev/null
# And finally, the liftover chain file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset net/combined.umd3_to_ars-ucd.net final_merged_chain/combined.umd3_to_ars-ucd.sorted.chain combined.umd3_to_ars-ucd.liftover.chain

```


Generating the data to compare against the reference will require at least two different data files (bed and fasta). I will try to test out my new liftover files against the GGP chip data that Geneseek has sent me.

```bash
perl -e 'for($x = 0; $x < 8; $x++){<>;} while(<>){chomp; $_ =~ s/\r//g; @s = split(/,/); print "chr$s[9]\t$s[10]\t$s[10]\t$s[1]\n";}' < GGP-F250_20000778_A1.csv > GGP-F250_20000778.umd3.snplocs.bed

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver GGP-F250_20000778.umd3.snplocs.bed combined.umd3_to_ars-ucd.liftover.chain GGP-F250_20000778.ars-ucd14.snplocs.bed GGP-F250_20000778.ars-ucd14.snplocs.unmapped

# Everything was unmapped. I wonder if the fact that I didn't have an end coordinate of +1 mattered here...
perl -lane '$F[2] += 1; print join("\t", @F);' < GGP-F250_20000778.umd3.snplocs.bed > GGP-F250_20000778.umd3.snplocs.onebase.bed
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver GGP-F250_20000778.umd3.snplocs.onebase.bed combined.umd3_to_ars-ucd.liftover.chain GGP-F250_20000778.ars-ucd14.snplocs.bed GGP-F250_20000778.ars-ucd14.snplocs.unmapped

# That was it
wc -l GGP-F250_20000778.ars-ucd14.snplocs.unmapped GGP-F250_20000778.ars-ucd14.snplocs.bed
   35090 GGP-F250_20000778.ars-ucd14.snplocs.unmapped	<- is really 17,545 unmapped
  203570 GGP-F250_20000778.ars-ucd14.snplocs.bed
```

Just as predicted, ~8% of probes were not mapped (due to "deletions" in new assembly). This may be due to repeatmasked segments. I will try my mapping approach and compare the two predictive methods to see how the coordinates pan out.

```bash
perl -e 'for($x = 0; $x < 8; $x++){<>;} while(<>){chomp; $_ =~ s/\r//g; @s = split(/,/); print "\>$s[1]\n$s[5]\n";}' < GGP-F250_20000778_A1.csv > GGP-F250_20000778_A1.probeseq.fa

# The main difference here is that I am going to map to an unrepeatmasked reference
bwa mem ../ARS-UCD1.0.14.clean.wIGCHaps.fasta GGP-F250_20000778_A1.probeseq.fa > GGP-F250_20000778_A1.probeseq.sam
perl -lane 'if($F[0] =~ /^@/){next;}elsif($F[1] & 2048 == 2048){next;}else{$s = 0; $e = 0; if($F[1] & 16 == 16){ $s = $F[3] - 1; $e = $s + 1;}else{$s = $F[3] + 50; $e = $s + 1;} print "$F[2]\t$s\t$e\t$F[0]";}' < GGP-F250_20000778_A1.probeseq.sam > GGP-F250_20000778_A1.probeseq.bwa.ars-ucd14.bed

# Easy enough, but is it sensitive enough?
```

Now to test the overlap of coordinates

```bash
module load bedtools/2.26.0
intersectBed -a GGP-F250_20000778_A1.probeseq.bwa.ars-ucd14.bed -b GGP-F250_20000778.ars-ucd14.snplocs.bed | wc -l
87027 	<- not that amazing
```

I'm going to do this quickly with R, as there are built-in functions in the tidyverse that will do the comparison quickly.

```R
umd3_coords <- read.delim("GGP-F250_20000778.umd3.snplocs.onebase.bed", header=FALSE)
colnames(umd3_coords) <- c("umdChr", "umdStart", "umdEnd", "snpname")

bwa_coords <- read.delim("GGP-F250_20000778_A1.probeseq.bwa.ars-ucd14.bed", header=FALSE)
colnames(bwa_coords) <- c("bwaChr", "bwaStart", "bwaEnd", "snpname")

liftover_coords <- read.delim("GGP-F250_20000778.ars-ucd14.snplocs.bed", header=FALSE)
colnames(liftover_coords) <- c("liftChr", "liftStart", "liftEnd", "snpname")

library(dplyr)
full_data <- left_join(umd3_coords, liftover_coords, by ="snpname")
full_data <- left_join(full_data, bwa_coords, by = "snpname")
write.csv(full_data, file="F250_joined_test_datafield.tab", quote=FALSE)
```

Now some quick summary stats.

```bash
perl -e '<>; while(<>){chomp; @s = split(/,/); $s[1] =~ s/chr//g; if($s[1] ne $s[8] && length($s[8]) <= 2){print join("\t", @s); print "\n";}}' < F250_joined_test_datafield.tab | wc -l
	2141 <- changed chromosome from UMD3 to ARS-UCD14
```

<a name="coordchanges"></a>
#### Coordinate changes

Bob just informed me of a big flaw in the infinium assays: Two types of transversions will cause problems with my assumptions of variant base position. In variants: [A/T], [T/A], [C/G] and [G/C], the variant base is exactly at the end of the probe, whereas in all other probes the variant is at the 3' + 1 position.

I will write a python script to process the sam files with this in mind.

First, let's extract the probe sequence fasta, but with the variant allele notation in the snp name:

```bash
perl -e 'for($x = 0; $x < 8; $x++){<>;} while(<>){chomp; $_ =~ s/\r//g; @s = split(/,/); print "\>$s[1]_$s[3]\n$s[5]\n";}' < GGP-F250_20000778_A1.csv > GGP-F250_20000778_A1.probeseq.fa

bwa mem ../ARS-UCD1.0.14.clean.wIGCHaps.fasta GGP-F250_20000778_A1.probeseq.fa > GGP-F250_20000778_A1.probeseq.sam
python3 ../../../binaries/python_toolchain/snpMapping/getProbeSeqCoords.py -s GGP-F250_20000778_A1.probeseq.sam -o GGP-F250_20000778_A1.probeseq.python.tab

intersectBed -a GGP-F250_20000778_A1.probeseq.python.tab -b GGP-F250_20000778.ars-ucd14.snplocs.bed | wc -l
95477 	<- that improved it by quite a bit
intersectBed -a GGP-F250_20000778_A1.probeseq.python.tab -b GGP-F250_20000778_A1.probeseq.bwa.ars-ucd14.bed | wc -l
11279	<- wow! That's only half of the initial markers!

grep '*' GGP-F250_20000778_A1.probeseq.python.sorted.tab | wc -l
721 <- number of unaligned markers from BWA mem

# Now to check to see if liftover helped to place unaligned markers
grep '*' GGP-F250_20000778_A1.probeseq.python.sorted.tab | perl -lane '$F[3] =~ s/_\[.\/.\]//g; system("grep $F[3] GGP-F250_20000778.ars-ucd14.snplocs.sorted.bed");' | wc -l
62	<- only 62 placed by liftover out of 721

grep '*' GGP-F250_20000778_A1.probeseq.python.sorted.tab | perl -lane '$F[3] =~ s/_\[.\/.\]//g; system("grep $F[3] GGP-F250_20000778.ars-ucd14.snplocs.unmapped");' | wc -l
661	<- about 3.7% of the unmapped in the liftover file. 

# Checking liftover coords against bwa coords for ~ 1bp proximity
closestBed -a GGP-F250_20000778_A1.probeseq.python.sorted.tab -b GGP-F250_20000778.ars-ucd14.snplocs.sorted.bed -d | perl -lane 'if($F[0] eq "*"){next;}elsif($F[8] <= 2){print $_;}' | wc -l
198617	<- 89.8% of the BWA mem alignments, and 98% of the liftover coords

# Checking the number of unmapped liftover coordinates that account for this discrepency
closestBed -a GGP-F250_20000778_A1.probeseq.python.sorted.tab -b GGP-F250_20000778.ars-ucd14.snplocs.sorted.bed -d | perl -lane 'if($F[0] eq "*"){next;}elsif($F[8] > 2){$F[3] =~ s/_\[.\/.\]//g; system("grep $F[3] GGP-F250_20000778.ars-ucd14.snplocs.unmapped");}' | wc -l
15370	<- 15370 / 17545 = 87.6% of the unmapped probes

# Just to make sure that all of my comparisons are unbiased, let's prepare a flat venn cross-section to make sure that everything is on point
closestBed -a GGP-F250_20000778_A1.probeseq.python.sorted.tab -b GGP-F250_20000778.ars-ucd14.snplocs.sorted.bed -d | perl -lane 'if($F[8] > 2){$F[3] =~ s/_\[.\/.\]//g; print $F[3];}' > GGP_bwa_liftover_not_same.list
perl -lane 'if($F[0] eq "*"){$F[3] =~ s/_\[.\/.\]//g; print $F[3];}' < GGP-F250_20000778_A1.probeseq.python.sorted.tab > GGP_bwa_unmapped.list
perl -lane 'if($F[0] ne "*"){$F[3] =~ s/_\[.\/.\]//g; print $F[3];}' < GGP-F250_20000778_A1.probeseq.python.sorted.tab > GGP_bwa_mapped.list
perl -lane 'print $F[3]' < GGP-F250_20000778.ars-ucd14.snplocs.sorted.bed > GGP_liftover_mapped.list
grep -v '#' GGP-F250_20000778.ars-ucd14.snplocs.unmapped | perl -lane 'print $F[3];' > GGP_liftover_unmapped.list

perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl GGP_bwa_mapped.list GGP_bwa_unmapped.list GGP_liftover_mapped.list GGP_liftover_unmapped.list GGP_bwa_liftover_not_same.list
File Number 1: GGP_bwa_mapped.list
File Number 2: GGP_bwa_unmapped.list
File Number 3: GGP_liftover_mapped.list
File Number 4: GGP_liftover_unmapped.list
File Number 5: GGP_bwa_liftover_not_same.list
Set     Count
1;3     196381
1;3;5   7129
1;4     1574
1;4;5   15310
2;3     60
2;4     661

# This mirrors the numbers I've gathered above. 661 markers were unrecoverable and 7129 markers represent a problem
# Grepping out the problematic marker names
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_3_5 GGP_bwa_mapped.list GGP_bwa_unmapped.list GGP_liftover_mapped.list GGP_liftover_unmapped.list GGP_bwa_liftover_not_same.list &> GGP_bwa_liftover_major_coorddiffs.list
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 2_4 GGP_bwa_mapped.list GGP_bwa_unmapped.list GGP_liftover_mapped.list GGP_liftover_unmapped.list GGP_bwa_liftover_not_same.list &> GGP_bwa_liftover_unrecovered.list

# I used vim to remove the file number header from each list file
# I am also going to standardize the SNP names so that I can load them into R for quick comparisons
perl -ne '$_ =~ s/_\[.\/.\]//g; print $_;' < GGP-F250_20000778_A1.probeseq.python.sorted.tab > GGP-F250_20000778_A1.probeseq.python.sorted.reformat.tab
perl -lane 'if($F[0] =~ /^@/){next;}elsif($F[1] & 2048 == 2048){next;}else{$F[0] =~ s/_\[.\/.\]//g; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[5]";}' < GGP-F250_20000778_A1.probeseq.sam > GGP-F250_20000778_A1.probeseq.brief.sam

#TODO: compare the data in R to see if any common trends pop up
```

Now to add them all together in R to try to identify trends and to merge everything together.

```R
library(dplyr)
bwa_coords <- read.delim("GGP-F250_20000778_A1.probeseq.python.sorted.reformat.tab", header=FALSE)
colnames(bwa_coords) <- c("bwaChr", "bwaStart", "bwaEnd", "snpname")

# This dataframe was loaded automatically from a saved workspace image
liftover_coords <- read.delim("GGP-F250_20000778.ars-ucd14.snplocs.bed", header=FALSE)
colnames(liftover_coords) <- c("liftChr", "liftStart", "liftEnd", "snpname")

head GGP-F250_20000778_A1.probeseq.brief.sam
colnames(sam_aligns) <- c("snpname", "samFlag", "samChr", "samPos", "samCigar")

full_data <- left_join(umd3_coords, liftover_coords, by ="snpname")
full_data <- left_join(full_data, bwa_coords, by = "snpname")
full_data <- left_join(full_data, sam_aligns, by = "snpname")

# Now to slice out the profiles of variants that differ in coordinates
differ_onebase <- filter(full_data, abs(bwaStart - liftStart) == 1)
count(differ_onebase) # 96170
differ_twobase <- filter(full_data, abs(bwaStart - liftStart) == 2)
count(differ_twobase) # 4766
differ_threebase <- filter(full_data, abs(bwaStart - liftStart) == 3)
count(differ_threebase) # 4132
exact_match <- filter(full_data, bwaStart == liftStart)
count(exact_match) # 94500

# Nearly all exact_matches are reverse oriented reads
# Nearly all differing coordinates are forward oriented reads
pdf(file="bwa_vs_liftover_sam_flag_comp.pdf", useDingbats=FALSE)
full_data <- mutate(full_data, cdiff = abs(bwaStart - liftStart))
boxplot(cdiff ~ samFlag, data=full_data, xlab = "SAM read orientation", ylab = "Coordinate bp difference", main = "Base bias in alignment coordinate algorithm conversion")
dev.off()

# The plot was too compressed to really tell me much, unfortunately
differ_onebase <- mutate(differ_onebase, cdiff = liftStart - bwaStart)
summary(differ_onebase)
# The majority (99.73%) of events were negative, so that means that the bwaStart is typically 1 bp more than the liftover start
summary(as.factor(differ_onebase$samFlag))
    0    16
95896   274
summary(as.factor(differ_onebase$cdiff))
   -1     1
96039   131

# They're pretty similar, let's try a chi square test
chisq.test(as.factor(differ_onebase$samFlag), as.factor(differ_onebase$cdiff), correct=FALSE)
# X-squared = 34738, df = 1, p-value < 2.2e-16

# Supposedly they're dependent, and a quick test of 50 columns does show dependency
# Solution: always subtract one base from the forward alignment position in the script
```

Testing my conclusions from above:

```bash
python3 ../../../binaries/python_toolchain/snpMapping/getProbeSeqCoords.py -s GGP-F250_20000778_A1.probeseq.sam -o GGP-F250_20000778_A1.probeseq.python.tab -c
sortBed -i GGP-F250_20000778_A1.probeseq.python.tab > GGP-F250_20000778_A1.probeseq.python.sorted.tab

intersectBed -a GGP-F250_20000778_A1.probeseq.python.sorted.tab -b GGP-F250_20000778.ars-ucd14.snplocs.bed | wc -l
190775 	<- Thats 93.7% of all markers
```

Now to test this against the GH2 chip that George sent me. I will run the liftover and bwa mem approaches to see how they line up.

```bash
# Converting the GH2 chip to fasta and UMD3 SNP coords
perl -e '<>; while(<>){chomp; @s = split(/\t/); @bsegs = split(/\[/, $s[6]); print ">$s[1]_\[$s[7]\/$s[8]\]\n$bsegs[0]\n";}' < SNP_GH2.txt > SNP_GH2.fa
perl -e '<>; while(<>){chomp; @s = split(/\t/); if($s[4] == 0){next;} $e = $s[5] + 1; print "chr$s[4]\t$s[5]\t$e\t$s[1]_\[$s[7]\/$s[8]\]\n";}' < SNP_GH2.txt > SNP_GH2.umd3_coords.bed

# Generating liftover coords
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver SNP_GH2.umd3_coords.bed combined.umd3_to_ars-ucd.liftover.chain SNP_GH2.ars-ucd14_coords.bed SNP_GH2.ars-ucd14_coords.unmapped
wc -l SNP_GH2.ars-ucd14_coords.bed SNP_GH2.ars-ucd14_coords.unmapped
  159 SNP_GH2.ars-ucd14_coords.bed		<- should be 322 total markers, but is less because of a lack of chr assignments
    8 SNP_GH2.ars-ucd14_coords.unmapped	<- only 4 unmapped

# Generating BWA alignments
bwa mem ../ARS-UCD1.0.14.clean.wIGCHaps.fasta SNP_GH2.fa > SNP_GH2.ars-ucd14.sam
python3 ../../../binaries/python_toolchain/snpMapping/getProbeSeqCoords.py -s SNP_GH2.ars-ucd14.sam -o SNP_GH2.ars-ucd14.tab

# Note that I added a new "-c" flag to toggle the 1bp correction algorithm in the above python script between now and then.
```

<a name="hdremapping"></a>
#### HD probe remapping

I will try to remap the BovineHD probes using both liftover and BWA to compare positions and profiles.

```bash
perl -e 'for($x = 0; $x < 8; $x++){<>;} while(<>){chomp; $_ =~ s/\r//g; @s = split(/,/); print "\>$s[1]_$s[3]\n$s[5]\n";}' < bovinehd-manifest-b.csv > bovinehd_illumina.probeseq.fa
perl -e 'for($x = 0; $x < 8; $x++){<>;} while(<>){chomp; $_ =~ s/\r//g; @s = split(/,/); $e = $s[10] + 1; print "chr$s[9]\t$s[10]\t$e\t$s[1]\n";}' < bovinehd-manifest-b.csv > bovinehd_illumina.snplocs.bed

# I removed the terminal controls from the end of the file manually with Vim

# Now to run the comparison!
bwa mem ../ARS-UCD1.0.14.clean.wIGCHaps.fasta bovinehd_illumina.probeseq.fa > bovinehd_illumina.probeseq.sam
python3 /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/snpMapping/getProbeSeqCoords.py -s bovinehd_illumina.probeseq.sam -o bovinehd_illumina.probeseq.bwa.bed -c

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver bovinehd_illumina.snplocs.bed combined.umd3_to_ars-ucd.liftover.chain bovinehd_illumina.snplocs.ars-ucd14.bed bovinehd_illumina.snplocs.ars-ucd14.unmapped
wc -l bovinehd_illumina.snplocs.ars-ucd14.bed bovinehd_illumina.snplocs.ars-ucd14.unmapped
  646503 bovinehd_illumina.snplocs.ars-ucd14.bed
  262918 bovinehd_illumina.snplocs.ars-ucd14.unmapped	<- to be expected, given the repeatmasking!

intersectBed -a bovinehd_illumina.probeseq.bwa.bed -b bovinehd_illumina.snplocs.ars-ucd14.bed | wc -l
608711	<- OK, so we have about 39,000 that don't match up

# Testing to see how many BWA unmapped were lifted over
perl -lane 'if($F[0] eq "*"){print $F[3];}' < bovinehd_illumina.probeseq.bwa.bed > bovinehd_illumina.probeseq.bwa.unmappedprobes.list
perl -lane 'print $F[3];' < bovinehd_illumina.snplocs.ars-ucd14.bed > bovinehd_illumina.snplocs.ars-ucd14.liftovermappedprobes.list

perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl bovinehd_illumina.probeseq.bwa.unmappedprobes.list bovinehd_illumina.snplocs.ars-ucd14.liftovermappedprobes.list
File Number 1: bovinehd_illumina.probeseq.bwa.unmappedprobes.list
File Number 2: bovinehd_illumina.snplocs.ars-ucd14.liftovermappedprobes.list
Set     Count
1       3282
1;2     136		<- 136 were recovered by liftover
2       646367
```

Loading them into R to check coordinates

```bash
library(dplyr)
hd_bwa <- read.delim("bovinehd_illumina.probeseq.bwa.bed", header=FALSE)
colnames(hd_bwa) <- c("bwaChr", "bwaStart", "bwaEnd", "snpname")

hd_lift <- read.delim("bovinehd_illumina.snplocs.ars-ucd14.bed", header=FALSE)
colnames(hd_lift) <- c("liftChr", "liftStart", "liftEnd", "snpname")

hd_full <- left_join(hd_bwa, hd_lift, by="snpname")
hd_full <- mutate(hd_full, cdiff = bwaStart - liftStart)
hd_full <- mutate(hd_full, abscdiff = abs(bwaStart - liftStart))
hd_full <- mutate(hd_full, bwaunmap = bwaChr == "*")

summary(hd_full)
     bwaChr          bwaStart             bwaEnd            snpname
 1      : 46637   Min.   :        0   Min.   :        1   Length:777962
 2      : 40004   1st Qu.: 21406256   1st Qu.: 21406256   Class :character
 X      : 37486   Median : 43790343   Median : 43790344   Mode  :character
 3      : 35631   Mean   : 49762810   Mean   : 49762811
 4      : 34898   3rd Qu.: 72059510   3rd Qu.: 72059511
 5      : 34799   Max.   :167313065   Max.   :167313066
 (Other):548507
    liftChr         liftStart            liftEnd              cdiff
 1      : 37957   Min.   :      282   Min.   :      283   Min.   :-151069026
 2      : 32935   1st Qu.: 21791194   1st Qu.: 21791195   1st Qu.:         0
 8      : 32051   Median : 43860955   Median : 43860956   Median :         0
 4      : 28931   Mean   : 49709112   Mean   : 49709113   Mean   :    206191
 3      : 28852   3rd Qu.: 71866844   3rd Qu.: 71866844   3rd Qu.:         0
 (Other):485777   Max.   :167302265   Max.   :167302266   Max.   : 157140187
 NA's   :131459   NA's   :131459      NA's   :131459      NA's   :131459
    abscdiff          bwaunmap
 Min.   :        0   Mode :logical
 1st Qu.:        0   FALSE:774544
 Median :        0   TRUE :3418
 Mean   :   304871   NA's :0
 3rd Qu.:        0
 Max.   :157140187
 NA's   :131459

count(filter(hd_full, abscdiff >= 1)) # 37934 
```
<a name="testcompletion"></a>
## Testing regions of the assembly for completion

I am going to check several canu haplotypes for alignment to the v14 reference.

> assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
# Here is one tig section that Bob thinks is missing from the assembly
samtools faidx /mnt/nfs/nfs2/dbickhart/dominette_asm/canu.mhap.all.fasta tig00000805:1691178-1744155 > tig00000805_1691178-1744155.fa

bwa mem ../ARS-UCD1.0.14.clean.wIGCHaps.fasta tig00000805_1691178-1744155.fa > tig00000805_1691178-1744155.sam
perl ~/sperl/sequence_data_scripts/BriefSamOutFormat.pl -s tig00000805_1691178-1744155.sam
tig00000805:1691178-1744155     16      15      49017753        49022770        52978   60

# Now to generate nucmer plots
# liftover for the region from chr4 umd3
perl check_coordinates.pl final_merged_chain/chr4.chain 96770424 96823498 tig00000805_chr4_coords.tab
# That has it on chr4 in the new assembly as well. Going to make plots for both regions

samtools faidx ../ARS-UCD1.0.14.clean.wIGCHaps.fasta 4:95900318-96000000 > ars_ucd_sub_4_95900318_96000000.fa
samtools faidx ../ARS-UCD1.0.14.clean.wIGCHaps.fasta 15:49017753-49022770 > ars_ucd_sub_15_49017753_49022770.fa

sh /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh tig00000805_1691178-1744155.fa ars_ucd_sub_4_95900318_96000000.fa
sh /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh tig00000805_1691178-1744155.fa ars_ucd_sub_15_49017753_49022770.fa

# Both gave errors for no alignment data to plot!
rm ars_ucd_sub_*
# We're trying whole chromosomes now
samtools faidx ../ARS-UCD1.0.14.clean.wIGCHaps.fasta 4 > ars_ucd_sub_4.fa
sh /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ars_ucd_sub_4.fa tig00000805_1691178-1744155.fa
# No alignment data here as well!

samtools faidx /mnt/nfs/nfs2/dbickhart/dominette_asm/canu.mhap.all.fasta tig00000805 > tig00000805.fa
sh /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ars_ucd_sub_4.fa tig00000805.fa
# I resized the gnuplot image, but only 3kb aligned!

samtools faidx ../ARS-UCD1.0.14.clean.wIGCHaps.fasta 15 > ars_ucd_sub_15.fa
sh /mnt/nfs/nfs2/bickhart-users/binaries/run_nucmer_plot_automation_script.sh ars_ucd_sub_15.fa tig00000805.fa
# 48027402 51365201

ARS-BFGL-NGS-95780      4       95965226
ARS-BFGL-NGS-95780      4       95965226        +       4       96754893
Hapmap28499-BTA-142459  4       96007967        +       4       96861591


samtools faidx ../ARS-UCD1.0.14.clean.wIGCHaps.fasta 4:95971185-96013926 > ars_ucd_probregion_4.fa
samtools faidx /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa chr4 > umd3_chr4.fa
sh ../../../binaries/run_nucmer_plot_automation_script.sh umd3_chr4.fa ars_ucd_probregion_4.fa

````

<a name="nonrepeatmasked"></a>
## Non-repeatmasked chain of assembly file

I  need to check to see if I can get this working on my non-repeatmasked fastas for full liftover of v14 to umd3.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
mkdir chrchunks
perl -lane 'system("samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0] > chrchunks/$F[0].fa");' < ARS-UCD1.0.14.clean.wIGCHaps.fasta.fai

# Converting to 2bit assemblies
for i in chrchunks/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/faToTwoBit $i chrchunks/${name}.2bit; done

module load blat
for i in chrchunks/*.fa; do sbatch -p assemble3 faSplit_blat_align.sh umd3_kary_unmask_ngap.2bit $i; done

# testing minimap2 to see if it is faster than blat
ln -s /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa umd3_kary_unmask_ngap.fa
sbatch faSplit_minimap2_align.sh umd3_kary_unmask_ngap.fa chrchunks/10.fa

# I tried using the "genome to genome alignment" metric for 95% similarity, but the output of minimap2 was too disjointed
# testing this as if the read is just a very long pacbio read
sbatch --nodes=1 --mem=18000 --ntasks-per-node=4 --exclude=vm-agil-[233,235] --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -ax map-pb umd3_kary_unmask_ngap.fa.mmi /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/chrchunks/test_sub_chunks/10_100000.fa > OLD.10_100000.fa.sam; python /mnt/nfs/nfs2/bickhart-users/binaries/fusioncatcher/bin/sam2psl.py -i OLD.10_100000.fa.sam -o OLD.10_100000.fa.psl"

# That was even worse. Returning to the 95% asm comparison run and trying the rest of the pipeline
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/faToTwoBit ARS-UCD1.0.14.clean.wIGCHaps.fasta ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit
sbatch -p assemble2 psl_merge.sh umd3_kary_unmask_ngap.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit
```

<a name="chain25umd"></a>
#### Chaining v25 to UMD3

> Lewis: /home/bickhartd/bickhartd/ars

```bash
# Because I don't want to piss off the Lewis cluster admins, I'm going to write scripts for most of this processing
sbatch splitChrs.pl ARS-UCDv1.0.25.fasta chrchunks

sbatch --nodes=1 --mem=8000 --ntasks-per-node=1 --wrap="../kentUtils/bin/linux.x86_64/faToTwoBit umd3_kary_unmask_ngap.fa umd3_kary_unmask_ngap.2bit"

# Making the ooc file
sbatch --nodes=1 --mem=25000 --ntasks-per-node=1 --wrap="blat ../umd3/umd3_kary_unmask_ngap.2bit /dev/null /dev/null -makeOoc=11.ooc -repMatch=1024"

for i in chrchunks/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --account=biocommunity -p BioCompute ../bin/faSplit_blat_align.sh ../umd3/umd3_kary_unmask_ngap.2bit $i; done

```

<a name="mapbobsmarkers"></a>
## Mapping Bob's markers to v25 from NCBI

Using Bob's script to download and standardize the assembly.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi

##### download_process_reference.sh

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=25000
#SBATCH --ntasks-per-node=8

module load bwa
module load samtools
wget 'ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/NK/LS/NKLS02/*'
unpigz *nt.gz

wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/Assembled_chromosomes/seq/bt_alt_Btau_5.0.1_chrY.fa.gz'
unpigz bt_alt_Btau_5.0.1_chrY.fa.gz

cat *.fsa_nt >ARS-UCD1.2.ORIG.fa
# replace the genbank chr name with just the accession
sed s/'gi|355477162|ref|'// <bt_alt_Btau_5.0.1_chrY.fa >bt_alt_Btau_5.0.1_chrY.fa2
sed s/'|'// <bt_alt_Btau_5.0.1_chrY.fa2 >bt_alt_Btau_5.0.1_chrY.fa3

cat ARS-UCD1.2.ORIG.fa bt_alt_Btau_5.0.1_chrY.fa3 >ARS-UCD1.2.PlusY.fa

grep '>' ARS-UCD1.2.PlusY.fa >ARS-UCD1.2.PlusY.fa.NAMES

samtools faidx ARS-UCD1.2.PlusY.fa
bwa index ARS-UCD1.2.PlusY.fa
```

Downloading markers and queuing analysis.

```bash
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 11ouBkPXCNr0blgkFsAcmMGIGt9daH3iM 9913_CHIP_DEREK_A.csv.1.fasta.gz
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1-n9iFE1qUMy4OqZrySVclEBahYfo5qYG 9913_CHIP_DEREK_A.csv.2.fasta.gz
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1dC6zHFv31RC4podsUZaFdZhz7sRT5HAW 9913_CHIP_DEREK_A.csv.3.fasta.gz
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1uXAuaQe-QTuHulu-53owIOqq1qhWRvfu 9913_CHIP_DEREK_B.csv.1.fasta.gz
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1pvRxVQMQZZyEWo25h-rb70H7uLSjuwnW 9913_CHIP_DEREK_B.csv.2.fasta.gz
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1gIeU5X5Hw4OqaQ8IEYSqhlwmEAPY084N 9913_CHIP_DEREK_B.csv.3.fasta.gz

for i in 9913_CHIP_DEREK_*.gz; do echo $i; gunzip $i; done

for i in 9913_CHIP_DEREK_*.fasta; do echo $i; sbatch --mem=20000 -p assemble1 --ntasks-per-node=1 --nodes=1 --wrap="perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a ARS-UCD1.2.PlusY.fa -p $i -o $i; bwa mem ARS-UCD1.2.PlusY.fa $i > $i.sam"; done

# And for UMD3
module load bwa; for i in snp_remappings/*.fasta; do name=`basename $i`; echo $name; sbatch --mem=20000 --ntasks-per-node=1 --nodes=1 --wrap="perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a bt_ref_Bos_taurus_UMD_3.1.1.fasta -p $i -o $name; bwa mem bt_ref_Bos_taurus_UMD_3.1.1.fasta $i > $name.sam"; done

mv 9913_CHIP* ./umd3_mappings/
```