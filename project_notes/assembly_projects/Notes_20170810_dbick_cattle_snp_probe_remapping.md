# SNP probe remapping on cattle assemblies
---
*8/10/2017*

These are my notes for remapping UMD3-based SNP probes onto the new cattle assembly. My hope is to do as thorough of a job as possible and to resolve >99% of marker locations, if possible. 

## Table of contents


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
