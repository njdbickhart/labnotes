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