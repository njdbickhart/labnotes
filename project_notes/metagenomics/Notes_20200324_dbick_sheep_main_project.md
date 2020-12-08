# Sheep dataset main analysis
---
*3/24/2020*

These are my notes on analyzing the sheep fecal dataset for eventual publication.

## Table of Contents


## Setting up analysis pipeline

I've been asked to analyze many metagenome assemblies lately and I need an automated way to process them. My previous work on the rumen dataset generated many useful tables and summary stats, but I should automate this so I can crank the gears of the HPC while working on other projects.

I'll be using snakemake again. *shudder*

> Ceres: /project/forage_assemblies/sheep_project/pipeline_test

```bash
module load miniconda/3.6
/KEEP/rumen_longread_metagenome_assembly/python3

cp ~/python_toolchain/snakeMake/metaGenomeAnalysis/config.json ./

# Making sure the environments are all created before plowing ahead!
snakemake -n -s ~/python_toolchain/snakeMake/metaGenomeAnalysis/Snakefile --use-conda --create-envs-only

# OK, a billion bugs later, I think that I'm ready to test out the run
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/metaGenomeAnalysis/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/metaGenomeAnalysis/Snakefile --use-conda

# I tested out the taxify pipeline, but the same mucked up conda environment is screwing me up!
# Blobtools2 is apparently designed from the ground-up for conda, so I will need to refactor to use that in my pipeline
# Actually, it isn't! Damn...

```

## Eukaryote analysis

Just for some quick comparative analysis, I want to see how many contigs map to Eukaryotic parasite genomes. I will test against the following genomes (based on their identification using 18S sequence):

* [Blastocystis sp.](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Blastocystis_hominis/all_assembly_versions/GCA_000151665.1_ASM15166v1)
* [Eimeria Tenella](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Eimeria_tenella/latest_assembly_versions/GCA_000499545.1_ETH001)
* Eimeria Zuemii
* Cyniclomyces guttulatus
* [Haemonchus contortus](ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/haemonchus_contortus/PRJNA205202/haemonchus_contortus.PRJNA205202.WBPS14.genomic.fa.gz)

> Ceres: /project/forage_assemblies/sheep_project/eukcheck

```bash
module load minimap2 samtools
# Preparing the reference fasta files
for i in *.tar; do echo $i; tar -xvf $i; done
for i in ncbi-genomes-2020-04-16/*.gz; do echo $i; gunzip $i; done

mv ncbi-genomes-2020-04-16/GCA_000151665.1_ASM15166v1_genomic.fna ncbi-genomes-2020-04-16/bhominis.1_ASM15166v1_genomic.fna
mv ncbi-genomes-2020-04-16/GCA_000499545.1_ETH001_genomic.fna ncbi-genomes-2020-04-16/etenella.1_ETH001_genomic.fna
mv ncbi-genomes-2020-04-16/GCA_007637855.2_ASM763785v2_genomic.fna ncbi-genomes-2020-04-16/hcontortus.2_ASM763785v2_genomic.fna

# I will separate out the "eukaryotic" origin contigs, align them separately, and then align the whole shebang for full alignment. I'll then compare the alignment matches
perl -ne '@s = split(/\t/); print "$s[0] ";' < ../flye2_table.flye2.blobtools.blobDB.table.euk.txt; echo

samtools faidx ../flye2.contigs.fasta contig_1004 contig_10176 contig_1050 contig_10526 contig_10980 contig_1110 contig_11135 contig_1145 contig_11454 contig_1152 contig_11739 contig_11829 contig_11951 contig_12004 contig_121 contig_12108 contig_12145 contig_12223 contig_12289 contig_12460 contig_12647 contig_12653 contig_12692 contig_12705 contig_12896 contig_12925 contig_12946 contig_12955 contig_12994 contig_13209 contig_13318 contig_13408 contig_13453 contig_13497 contig_136 contig_13665 contig_13716 contig_14043 contig_14097 contig_14610 contig_14620 contig_14699 contig_14718 contig_14932 contig_1522 contig_15258 contig_15295 contig_15363 contig_15371 contig_15442 contig_15462 contig_1551 contig_1559 contig_1561 contig_15694 contig_15754 contig_15786 contig_1579 contig_15845 contig_1587 contig_16 contig_16084 contig_1609 contig_16136 contig_16142 contig_16148 contig_16437 contig_16462 contig_1653 contig_16569 contig_16592 contig_16601 contig_16750 contig_16848 contig_16977 contig_171 contig_17110 contig_17341 contig_17512 contig_17658 contig_17688 contig_17791 contig_17798 contig_17890 contig_18150 contig_18516 contig_1864 contig_18642 contig_18709 contig_18834 contig_18846 contig_18866 contig_18940 contig_18956 contig_19012 contig_2040 contig_2193 contig_2217 contig_2270 contig_2287 contig_2300 contig_2360 contig_242 contig_2454 contig_2582 contig_2588 contig_2631 contig_2634 contig_2667 contig_2726 contig_2768 contig_2794 contig_2810 contig_2815 contig_2880 contig_2906 contig_2926 contig_3023 contig_3124 contig_3132 contig_3356 contig_368 contig_3694 contig_3695 contig_3708 contig_3843 contig_3889 contig_3894 contig_4117 contig_4162 contig_4491 contig_4516 contig_4702 contig_4967 contig_5034 contig_506 contig_5084 contig_509 contig_5096 contig_5146 contig_515 contig_5318 contig_5321 contig_5417 contig_542 contig_5474 contig_5793 contig_6222 contig_6249 contig_634 contig_6341 contig_6435 contig_6443 contig_6452 contig_6555 contig_659 contig_6667 contig_6668 contig_6671 contig_679 contig_6810 contig_6870 contig_7032 contig_7034 contig_7046 contig_7070 contig_7074 contig_7148 contig_7175 contig_7372 contig_7428 contig_7579 contig_7595 contig_7604 contig_7613 contig_77 contig_770 contig_7747 contig_7753 contig_786 contig_7962 contig_7993 contig_8041 contig_837 contig_8537 contig_8597 contig_8724 contig_8786 contig_8791 contig_8871 contig_8950 contig_90 contig_9006 contig_9093 contig_9216 contig_9256 contig_932 contig_9471 contig_9679 contig_9886 contig_9887 contig_9934 > flye2.euk_contigs.fa

for i in flye2.euk_contigs.fa ../flye2.contigs.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; for j in ncbi-genomes-2020-04-16/*.fna; do query=`basename $j | cut -d'.' -f1`; echo $query; outfn=$name.$query.paf; echo $outfn; sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p priority -q msn --wrap="minimap2 -x asm10 -t 10 -o $outfn $j $i"; done; done

for i in *.paf; do echo -ne "$i\t"; perl -lane 'if($F[11] == 60){print $_;}' < $i | wc -l; done
flye2.contigs.bhominis.paf      3
flye2.contigs.etenella.paf      69
flye2.contigs.hcontortus.paf    188
flye2.euk_contigs.bhominis.paf  2
flye2.euk_contigs.etenella.paf  2
flye2.euk_contigs.hcontortus.paf        10

# Very strange! Only small regions mapping... Let's check this on the raw reads
for j in ncbi-genomes-2020-04-16/*.fna; do query=`basename $j | cut -d'.' -f1`; echo $query; sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p priority -q msn --wrap="minimap2 -x map-ont -t 10 -o $query.ccsmaps.paf $j sheep_fecal_combined_ccs.fastq"; done

# Preparing for R plots
perl -lane '$d = "ccsreads"; if($F[11] != 0){$r = $F[9] / $F[1]; print "$r\t$F[9]\t$F[10]\t$d";}' < hcontortus.ccsmaps.paf > hcontortus.ccsmaps.mapdata
perl -lane '$d = "asm"; if($F[11] != 0){$r = $F[9] / $F[1]; print "$r\t$F[9]\t$F[10]\t$d";}' < flye2.contigs.hcontortus.paf > hcontortus.asm.mapdata
```

```R
library(dplyr)
library(ggplot2)
library(gridExtra)

ccs <- read.delim("hcontortus.ccsmaps.mapdata", header=FALSE)
asm <- read.delim("hcontortus.asm.mapdata", header=FALSE)

colnames(ccs) <- c("qratio", "qlen", "rlen", "Data")
colnames(asm) <- c("qratio", "qlen", "rlen", "Data")

total <- bind_rows(ccs, asm)
total$Data <- as.factor(total$Data)

p1 <- ggplot(total, aes(y = qlen, x=Data, fill=Data)) + geom_boxplot(width=0.1) + geom_violin(trim=FALSE) + scale_fill_brewer(palette="Dark2") + labs(title="Query Length",x="Data", y = "Length (bp)") + scale_y_log10()
p2 <- ggplot(total, aes(y = rlen, x=Data, fill=Data)) + geom_boxplot(width=0.1) + geom_violin(trim=FALSE) + scale_fill_brewer(palette="Dark2") + labs(title="Reference Length",x="Data", y = "Length (bp)") + scale_y_log10()
p3 <- ggplot(total, aes(y = qratio, x=Data, fill=Data)) + geom_boxplot(width=0.1) + geom_violin(trim=FALSE) + scale_fill_brewer(palette="Dark2") + labs(title="Ratio query alignment",x="Data", y = "Proportion of reads")

pdf(file="alignment_to_hcontortus.pdf", useDingbats=FALSE)
grid.arrange(p1, p2, p3, nrow=2)
dev.off()
```

```bash
# Mash Screen analysis
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w RefSeq88n.msh flye2.euk_contigs.fa > flye2.euk_contigs.mashscreen.tab"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w RefSeq88n.msh ../flye2.contigs.fasta > flye2.contigs.mashscreen.tab"

# It mostly caught bacterial genomes... but I noticed that the other Eukaryotic genomes were missing
for i in ncbi-genomes-2020-04-16/*.fna; do echo $i; sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -S 0 -o $i $i"; done
mv ncbi-genomes-2020-04-16/*.msh ./
/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste RefSeq88n_andfecaleuk RefSeq88n.msh bhominis.1_ASM15166v1_genomic.fna.msh etenella.1_ETH001_genomic.fna.msh hcontortus.2_ASM763785v2_genomic.fna.msh

sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w RefSeq88n_andfecaleuk.msh flye2.euk_contigs.fa > flye2.euk_contigs.new.mashscreen.tab"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w RefSeq88n_andfecaleuk.msh ../flye2.contigs.fasta > flye2.contigs.new.mashscreen.tab"

### Interesting, most of the Nematodes map to Trichuris trichiura. Let's download that sequence to see how much we've assembled
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/trichuris_suis/PRJNA208415/trichuris_suis.PRJNA208415.WBPS14.genomic.fa.gz
ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS14.genomic.fa.gz

for i in trichuris_suis.PRJNA208415.WBPS14.genomic.fa.gz trichuris_trichiura.PRJEB535.WBPS14.genomic.fa.gz; do gunzip $i; done
for i in trichuris_suis.PRJNA208415.WBPS14.genomic.fa trichuris_trichiura.PRJEB535.WBPS14.genomic.fa; do ref=`echo $i | cut -d'.' -f1`; echo $ref; for j in flye2.euk_contigs.fa ../flye2.contigs.fasta; do query=`basename $j | cut -d'.' -f1,2`; echo $query; outfn=$ref.$query.paf; echo $outfn; sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p priority -q msn --wrap="minimap2 -x asm10 -t 10 -o $outfn $i $j"; done; done

# Now let's compare contig names to see how many are shared
for i in t*.paf; do name=`echo $i | cut -d'.' -f1,2,3`; echo $name; perl -lane 'if($F[11] > 10 && $F[10] > 5000 && $F[9] > 1000){print $F[0];}' < $i > $name.list; done
grep 'Nematoda' ../flye2_table.flye2.blobtools.blobDB.table.txt | perl -lane 'print $F[0];' > flye2.blobtools.nematoda.list
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o flye2.blobtools.nematoda.list trichuris_suis.flye2.euk_contigs.list trichuris_trichiura.flye2.euk_contigs.list trichuris_suis.flye2.contigs.list trichuris_trichiura.flye2.contigs.list
File Number 1: flye2.blobtools.nematoda.list
File Number 2: trichuris_suis.flye2.euk_contigs.list
File Number 3: trichuris_trichiura.flye2.euk_contigs.list
File Number 4: trichuris_suis.flye2.contigs.list
File Number 5: trichuris_trichiura.flye2.contigs.list
Set     Count
1       85
1;3;5   5
3;5     10
5       18

# Based on this, it's closest to Trichuris trichiura
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1 flye2.blobtools.nematoda.list trichuris_suis.flye2.euk_contigs.list trichuris_trichiura.flye2.euk_contigs.list trichuris_suis.flye2.contigs.list trichuris_trichiura.flye2.contigs.list > flye2.blobtools.uniq.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../flye2_table.flye2.blobtools.blobDB.table.txt -d '\t' -c 0 -l flye2.blobtools.uniq.list | perl -ne '@s = split(/\t/); print "$s[20]\n";' | sort | uniq -c
      1 Caenorhabditis remanei
      1 Diploscapter pachys
      2 Haemonchus placei
     55 Trichuris trichiura

grep 'Nematoda' ../flye2_table.flye2.blobtools.blobDB.table.txt | perl -ne '@s = split(/\t/); print "$s[20]\n";' | sort | uniq -c
      1 Angiostrongylus costaricensis
      1 Caenorhabditis remanei
      1 Diploscapter pachys
      1 Elaeophora elaphi
      2 Haemonchus placei
     84 Trichuris trichiura

# So, it looks like the H placei, C remanei and D pachys have no alignments to the other genomes and are likely not part of this species. Let's pull them and do comparative alignments with them.
samtools faidx ../flye2.contigs.fasta contig_1004 contig_1110 contig_1145 contig_11454 contig_1152 contig_121 contig_12108 contig_12460 contig_12692 contig_12925 contig_12946 contig_12994 contig_13209 contig_13497 contig_136 contig_13716 contig_14043 contig_14097 contig_14610 contig_1522 contig_15462 contig_1559 contig_15754 contig_15845 contig_16084 contig_16136 contig_16142 contig_1653 contig_16601 contig_16848 contig_16977 contig_17341 contig_17512 contig_18150 contig_18516 contig_18642 contig_18834 contig_18846 contig_18956 contig_2040 contig_2193 contig_2270 contig_2454 contig_2588 contig_2634 contig_2667 contig_2726 contig_2768 contig_2794 contig_2906 contig_3132 contig_368 contig_3694 contig_3695 contig_4117 contig_4516 contig_506 contig_5096 contig_5146 contig_5417 contig_542 contig_5474 contig_6222 contig_6443 contig_6452 contig_6555 contig_6668 contig_6671 contig_7034 contig_7070 contig_7074 contig_7148 contig_7372 contig_7613 contig_7747 contig_7753 contig_8041 contig_837 contig_8597 contig_8871 contig_90 contig_9216 contig_9471 contig_9679 > t.trichiura.contigs.fasta

samtools faidx t.trichiura.contigs.fasta

samtools faidx ../flye2.contigs.fasta scaffold_4292 contig_9854 contig_6834 scaffold_7509 contig_253 contig_11646 contig_4710 contig_2508 contig_12083 contig_511 contig_14949 contig_10992 contig_11974 contig_851 contig_13275 contig_2661 contig_5603 contig_10975 contig_8724 contig_7993 contig_16 contig_1579 contig_2287 contig_8791 contig_10980 contig_13408 contig_5793 contig_3708 contig_15754 contig_7747 contig_7070 contig_1561 contig_13497 contig_136 contig_16977 contig_12004 contig_13716 contig_13209 contig_2768 contig_18846 contig_15845 contig_2454 contig_837 contig_3694 contig_17341 contig_506 contig_7034 contig_6222 contig_3695 contig_9934 contig_14610 contig_2906 contig_2794 contig_90 contig_6443 contig_18516 contig_16136 contig_9216 contig_14043 contig_9471 contig_15694 contig_12108 contig_5096 contig_5034 contig_542 contig_7372 contig_7613 contig_7753 contig_1559 contig_16848 contig_2270 contig_2193 contig_6668 contig_1653 contig_368 contig_5417 contig_1145 contig_16601 contig_7148 contig_1004 contig_18956 contig_8041 contig_11454 contig_12946 contig_18834 contig_18150 contig_7074 contig_5474 contig_14097 contig_2726 contig_12692 contig_17512 contig_4516 contig_1522 contig_4117 contig_16084 contig_9679 contig_1152 contig_7595 contig_121 contig_8871 contig_12994 contig_6452 contig_18642 contig_2588 contig_16142 contig_2040 contig_2634 contig_6671 contig_12925 contig_5146 contig_6555 contig_3132 contig_15462 contig_2667 contig_1110 contig_8597 contig_12460 > t.trichiura.extended.fasta

samtools faidx t.trichiura.extended.fasta

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i t.trichiura.contigs.fasta.fai
Total length:   22389789
Num contigs:    84
N50 length:     11890218
N50 value:      724621
L50 value:      10
Max:    1952461
Min:    6484

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i t.trichiura.extended.fasta.fai
Total length:   27587529
Num contigs:    118
N50 length:     14285003
N50 value:      580020
L50 value:      13
Max:    1952461
Min:    5234

module load busco4/4.0.2
for i in t_contigs t_extended; do busco_configurator.py /software/7/apps/busco4/4.0.2/config/config.ini $i/config.ini; done
cp -Rp /software/7/apps/augustus/3.3.2/config/ ./AUGUSTUS_CONFIG

for i in contigs extended; do echo $i; fold="t_"${i}; asm="t.trichiura."${i}".fasta"; export BUSCO_CONFIG_FILE=/project/forage_assemblies/sheep_project/eukcheck/$fold/config.ini; sbatch --nodes=1 --mem=45000 --ntasks-per-node=20 -p priority -q msn --wrap="busco -i $asm -c 20 -o "t_busco_"${i} -m geno --offline -l /project/reference/data/BUSCO/v4/lineages/nematoda_odb10 --config $BUSCO_CONFIG_FILE"; done

# Aligning CCS reads to the refseq genome 
sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p priority -q msn --wrap="minimap2 -x map-ont -t 10 -o t_trichiura.ccsmaps.paf trichuris_trichiura.PRJEB535.WBPS14.genomic.fa sheep_fecal_combined_ccs.fastq"

perl -e '%h; while(<>){chomp; @s = split(/\t/); if($s[11] > 0 && ($s[9] / $s[1] > 0.5 || $s[10] / $s[6] > 0.5)){$h{$s[0]} = 1;}} foreach my $k (keys(%h)){print "$k\n";}' < t_trichiura.ccsmaps.paf > t_trichiura.reads.list

# Attempting a rescue assembly separate from the metagenome
python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -q sheep_fecal_combined_ccs.fastq -l t_trichiura.reads.list -o sheep_fecal_combined_ccs.ttrichuria.fastq

sbatch --nodes=1 --mem=200000 --ntasks-per-node=30 -p priority -q msn flye --pacbio-corr sheep_fecal_combined_ccs.ttrichuria.dedup.fastq -g 70m -o ttrichuria_subasm -t 30
```

## Final sheep assembly

#### Eukaryote mash screen test

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w /project/forage_assemblies/sheep_project/eukcheck/RefSeq88n_andfecaleuk.msh flye3.contigs.fasta > flye3.contigs.mashscreen.tab"
```

I ran into issues with the classification prior to this, but let's continue with the analysis using the raw reads and the new flye4 reference.

```bash
# First the Q20 CCS reads
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash screen -p 15 -w /lustre/project/forage_assemblies/sheep_project/eukcheck/RefSeq88n_andfecaleuk.msh /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/*Q20.fasta.gz > sheep_ccs_read_screen.tab"

# Next let's try a distance strategy for the contigs...
sbatch --nodes=1 --mem=20000 --ntasks-per-node=15 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -p 15 -i -v 0.05 -d 0.10 /lustre/project/forage_assemblies/sheep_project/eukcheck/RefSeq88n_andfecaleuk.msh flye4.contigs.fasta > flye4_mash_dist_euk.tab"
```

#### Running a megahit assembly for comparison

> Ceres: /luster/project/forage_assemblies/sheep_project/complete_flye

```bash
sbatch --nodes=1 --mem=500000 --ntasks-per-node=60 -p priority-mem -q msn-mem --wrap="megahit -m 500000000000 -t 60 --min-contig-len 1000 --out-prefix sheep_megahit -1 /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz -2 /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz"
```

#### Running the pipeline on both assemblies

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/metaGenomeAnalysis/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/metaGenomeAnalysis/Snakefile --use-conda

# URG! I forgot that megahit doesn't use proper fasta format! I need to reformat to 60 bp per line
python3 ~/python_toolchain/snakeMake/metaGenomeAnalysis/scripts/convertMegahitContigNames.py megahit.contigs.fasta megahit.reformat.fasta

# Rerunning for flye4 assembly
module load miniconda/3.6; sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/metaGenomeAnalysis/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/metaGenomeAnalysis/Snakefile --use-conda

# Copying clr datasets for the run
for i in clr1 clr2 clr3; do name=`echo "sheep_"$i`; echo $name; cp /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/$name/output/${name}_pilon_2.fasta.gz /lustre/project/forage_assemblies/sheep_project/complete_flye/$i.contigs.fasta.gz; done

module load miniconda/3.6; sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p priority -q msn snakemake --cluster-config ~/python_toolchain/snakeMake/metaGenomeAnalysis/cluster.json --cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}" -p --jobs 250 -s ~/python_toolchain/snakeMake/metaGenomeAnalysis/Snakefile --use-conda
```

## Generating CLR reads for CCS comparison

> Ceres: /lustre/project/rumen_longread_metagenome_assembly

```bash
module load pbsuite_new/15.8.24
for i in m54337U_200203_184822.subreads.bam m54337U_200211_145859.subreads.bam; do echo $i; sbatch --nodes=1 --mem=35000 --ntasks-per-node=3 -p priority -q msn -t 1-0 --wrap="bamToFastq $i > $i.fastq"; done

# OK, let's see if we can pull CCS reads from the fastqs. It looks like the run name is the first part of the read, with the ZMW being between '/' brackets. I want to calculate the number of reads present and then the average lengths
# Getting CCS readnames first
gunzip -c m54337U_200213_024013.Q20.fasta.gz | grep '>' | perl -lane '$_ =~ s/>//; print $_;' > m54337U_200213_024013.Q20.fasta.rnames

# Calculating the subread length statistics
perl calc_ccs_stats.pl m54337U_200213_024013.Q20.fasta.rnames m54337U_200213_024013.fastq.gz > m54337U_200213_024013.stats
3204385
```

```R
library(dplyr)
library(ggplot2)
library(gridExtra)
data <- read.delim("m54337U_200213_024013.stats", header=FALSE)
colnames(data) <- c("Read", "Length")

# Unique count of reads
length(unique(data$Read))
[1] 2699414

# Creating a dataframe to summarize
rsummaries <- data %>% group_by(Read) %>% mutate(Ct = n(), Avg = mean(Length)) %>% summarize(RCt = n(), AvgLen = mean(Avg))
rsummaries %>% summarize(Avg = mean(RCt), Min = min(RCt), Max = max(RCt))
# A tibble: 1 x 3
    Avg   Min   Max
  <dbl> <dbl> <dbl>
1  15.5     4  2179

# Let's see how many reads fit within the minimum
tally(rsummaries[rsummaries$RCt == 4,])
# A tibble: 1 x 1
      n
  <int>
1     3

rsummaries[rsummaries$RCt == 4,]
# A tibble: 3 x 3
       Read   RCt AvgLen
      <int> <int>  <dbl>
1  37685763     4 16090.
2  87163605     4 10753.
3 177343308     4 11350.

data[data$Read %in% c(37685763, 87163605, 177343308),]
              Read Length
8533757   37685763  15689
8533758   37685763  16273
8533759   37685763  16119
8533760   37685763  16277
20080347  87163605  10827
20080348  87163605  10695
20080349  87163605  10928
20080350  87163605  10561
41037503 177343308  11444
41037504 177343308  11319
41037505 177343308  11390
41037506 177343308  11245

# Just plotting the summaries for viewing later
pdf(file="m54337U_200213_024013.stats.pdf", useDingbats=FALSE)
p1 <- ggplot(rsummaries, aes(x="", y=RCt)) + geom_boxplot() + labs(title="SubRead counts")
p2 <- ggplot(rsummaries, aes(x="", y=AvgLen)) + geom_boxplot() + labs(title="Average SubRead Length")
grid.arrange(p1, p2, nrow=2)
dev.off()
```

#### Sequel I reads

First, to extract the reads from the bams.

```bash
module load pbsuite_new/15.8.24 samtools

for i in subreads.bam/*.bam; do echo $i; sbatch --nodes=1 --mem=35000 --ntasks-per-node=3 -p priority -q msn -t 1-0 --wrap="bamToFastq $i > $i.fastq"; done

# Now to get the stats
samtools faidx sheep_poop_sequel1_CCS.fastq
perl -lane 'print $F[1];' < sheep_poop_sequel1_CCS.fastq.fai > sheep_poop_sequel1_CCS.fastq.rnames

perl calc_ccs_stats.pl sheep_poop_sequel1_CCS.fastq.rnames subreads.bam/m54033_180827_171729.subreads.bam.fastq.gz > m54033_180827_171729.clr.stats

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f m54033_180919_161442.clr.stats -c 0 -d '\t' | perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\n";}' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   206553
Sum:    2286804
Minimum 3
Maximum 158
Average 11.071270
Median  10
Standard Deviation      5.647489
Mode(Highest Distributed Value) 6

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f m54033_180919_161442.clr.stats -c 0 -d '\t' | perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[1]\n";}' | perl -lane 'if($F[0] == 3){print $F[0];}' | wc -l
7

# Converting this to fasta
perl -e '%read; while($r = <>){chomp($r); $s = <>; chomp($s); <>; <>; $r =~ s/^@//; if(exists($read{$r})){next;} $s =~ s/(.{1,60})/$1\n/g; $read{$r} = 1; print ">$r\n$s";}' < sheep_poop_sequel1_CCS.fastq > sheep_poop_sequel1_CCS.fasta
```

So, I will keep the "skip the first read" philosophy to avoid screwing up the pipeline and sampling shorter CCS segments. Looks like it won't matter in most cases anyways!

### Generating the CLR subreads

So, I think that I should ignore the first read in all cases and then just extract up to the fourth subread in all cases. I wrote a python script to do this on a file-by-file basis.

```bash
# testing this out on the first read set
sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 -p priority -q msn --wrap="python3 ~/python_toolchain/assembly/extractPacbioCLRFromCCSData.py -f m54337U_200213_024013.Q20.fasta.gz -q m54337U_200213_024013.fastq.gz -o m54337U_200213_024013_clr"

# Had some weird issues with downloading files from gdrive. Half didn't make it. Working on the remainder for now
for i in m54337U_200203_184822 m54337U_200211_145859 m54337U_200220_195233; do echo $i; sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 -p priority -q msn --wrap="python3 ~/python_toolchain/assembly/extractPacbioCLRFromCCSData.py -f $i.Q20.fasta.gz -q $i.fastq.gz -o ${i}_clr"; done

# And the rest
for i in m54337U_200214_212611 m54337U_200222_075656 m54337U_200223_200916 m54337U_200227_201603; do echo $i; sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 -p priority -q msn --wrap="python3 ~/python_toolchain/assembly/extractPacbioCLRFromCCSData.py -f $i.Q20.fasta.gz -q $i.fastq.gz -o ${i}_clr"; done
```

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/

```bash
# Now, let's test a metaflye run on this
module load miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200203_184822_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200211_145859_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200213_024013_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200214_212611_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200220_195233_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200222_075656_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200223_200916_clr.1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/m54337U_200227_201603_clr.1.fastq.gz -t 70 -i 2 -m 4000 --asm-coverage 40 --meta -o sheep_clr1
```

#### Sequel I reads

Generating the pseudo CLR fastqs needed

```bash
for i in subreads.bam/*.fastq.gz; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=35000 --ntasks-per-node=2 -p priority -q msn --wrap="python3 ~/python_toolchain/assembly/extractPacbioCLRFromCCSData.py -f sheep_poop_sequel1_CCS.fasta.gz -q $i -o ${name}_clr"; done

# Some libraries must not have made it to CCS. There were several that were empty after processing
# Most had losses of 1-10 reads in the third file, though some were perfectly equal
for i in *_clr*.fastq; do echo $i; sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p priority -q msn --wrap="gzip $i"; done

```

And now, time for the assembly of everything! Let's test one assembly out first

```bash
module load miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.1.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.1.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.1.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.1.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.1.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.1.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.1.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.1.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.1.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.1.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.1.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.1.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.1.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.1.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.1.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.1.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.1.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.1.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.1.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.1.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.1.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.1.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.1.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.1.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.1.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.1.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.1.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.1.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.1.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.1.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.1.fastq.gz -o sheep_clr1 --meta -t 70

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.2.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.2.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.2.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.2.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.2.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.2.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.2.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.2.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.2.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.2.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.2.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.2.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.2.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.2.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.2.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.2.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.2.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.2.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.2.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.2.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.2.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.2.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.2.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.2.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.2.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.2.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.2.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.2.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.2.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.2.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.2.fastq.gz -o sheep_clr2 --meta -t 70

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.3.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.3.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.3.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.3.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.3.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.3.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.3.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.3.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.3.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.3.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.3.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.3.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.3.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.3.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.3.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.3.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.3.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.3.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.3.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.3.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.3.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.3.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.3.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.3.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.3.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.3.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.3.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.3.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.3.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.3.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.3.fastq.gz -o sheep_clr3 --meta -t 70
```

Empty libraries for CCS:
* m54033_180827_171729
* m54033_181218_133601

#### Quick stats

```bash
perl -e '<>; <>; <>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[-3] > 5 && $s[-2] < 10){print "$s[0]\t$s[-3]\n";}}' < ../checkm/qa_table.txt > flye3_checkm_contigs.tab
perl -e '<>; <>; <>; while(<>){$_ =~ s/^\s+//; @s = split(/\s+/); if($s[-3] > 5 && $s[-2] < 10){print "$s[0]\t$s[-3]\n";}}' < flye_sheep_2/checkm_table.txt > flye4_checkm_contigs.tab

perl -e '$asm = "flye4"; chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = [$s[1]];} close IN; open(IN, "< $ARGV[1]"); <IN>; print "Contig\tCompleteness\tSize\tAvgCov\tASM\n"; while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){push(@{$h{$s[0]}}, ($s[1], $s[2]));}} close IN; foreach my $c (keys(%h)){print "$c\t" . join("\t", @{$h{$c}}) . "\t$asm\n";}' flye4_checkm_contigs.tab binning/metabat2/flye4.contigs/jgi_abund.txt > flye4_checkm_contigs.extended.tab
```

```R
library(dplyr)
library(ggplot2)
flye3 <- read.delim("flye3_long_contigs.tab", header=FALSE)
flye4 <- read.delim("flye4_long_contigs.tab", header=FALSE)
combined <- bind_rows(flye3, flye4)
head(combined)
colnames(combined) <- c("Len", "Cov", "Phylum", "ASM")

pdf("Xcoverages.pdf", useDingbats=FALSE)
ggplot(combined, aes(x=Phylum, y=Cov, fill=Phylum)) + geom_jitter(aes(colour=Phylum)) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() + theme(axis.title=element_blank(), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + theme(legend.position = "none") + scale_fill_brewer(palette="Dark2") + labs(title="XCoverage of 1Mbp+ contigs") + facet_wrap(~ASM)
dev.off()

pdf("len_byX.pdf", useDingbats=FALSE)
ggplot(combined, aes(x=Len, y=Cov, colour=Phylum)) + geom_point(shape=18) + geom_smooth(method=lm, linetype="dashed", color="darkred", fill="coral") + labs(title="X coverage by contig length") + facet_wrap(~ASM)
dev.off()

flye4 <- read.delim("flye4_checkm_contigs.extended.tab")
pdf("flye4_cov_by_completeness.pdf", useDingbats=FALSE)
ggplot(flye4, aes(x=Completeness, y=AvgCov, color="blue")) + geom_point(aes(size=Size), alpha=0.10) + labs(title="Coverage by CheckM Completeness", ylab = "Average Coverage (log10)") + scale_y_log10()
dev.off()
```

And quick stats for the CLR assemblies.

```bash
tail -n 8 sheep_clr*/*.log | perl -e '%h; while(<>){chomp; if($_ =~ /[=\[]/){next;} $_ =~ s/^\s+//; @s = split(/\t/); if(scalar(@s) < 2){next;} push(@{$h{$s[0]}}, $s[1]);} foreach $k (keys(%h)){print "$k\t" . join("\t", @{$h{$k}}) . "\n";}'
```

#### Pilon correction of CLR assemblies

I coopted a new snakemake workflow to automatically pilon correct assemblies. Let's test it out.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep

```bash
# Preparing for the pipeline
for i in sheep_clr1 sheep_clr2 sheep_clr3; do echo $i; mkdir $i/asm1; cp $i/assembly.fasta $i/asm1/$i.fa; done

cp ~/python_toolchain/snakeMake/pilonCorrection/default.json ./
for i in sheep_clr1 sheep_clr2 sheep_clr3; do echo $i; cp default.json $i/; done

module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/python3/
# Because the Ceres staff never update their modules, I had to update my virtual environment to contain all the dependencies I need. 

# Now I need to run the snakemake pipeline in each folder like so
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/pilonCorrection/Snakefile --cluster "sbatch --nodes=1 --ntasks-per-node={threads} --mem=33000 -p priority -q msn -o logs/{rule}.stdout" --latency-wait 40 --jobs 50


#NOTE: the script must have the same file extension for all sequence files. Since pilon automatically affixes a "fasta", that means everything else must have one!
```

## Statistical calculations

Here are some notes on generating statistics on the polished assembly files.

#### Read depth ratios (long and short reads)

Pavel raised a good point about the balance of the short and long reads to each other on the assembled contigs. I will align both files and then calculate the per-base depth on each assembly. I will also use proportion of total bases aligned instead of an absolute read value (but I will output this for comparisons anyways!).

First I need to align the CCS reads to each assembly.

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
for i in sheep_poop_sequel1_CCS sheep_poop_sequelII_CCS; do for j in flye4.contigs.fasta clr1.contigs.fasta clr2.contigs.fasta clr3.contigs.fasta; do echo $i; echo $j; sbatch --nodes=1 --mem=45000 --ntasks-per-node=10 -p priority -q msn --wrap="minimap2 -x asm10 -t 10 -N 50 --secondary=no $j /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/$i.fasta > $i.$j.paf"; done; done

for i in clr1.contigs clr2.contigs clr3.contigs flye4.contigs; do echo $i; for j in sheep_poop_sequel1_CCS sheep_poop_sequelII_CCS; do echo $j; sbatch --nodes=1 --mem=50000 --ntasks-per-node=3 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/calcReadBalanceASM.py -o $j.$i.fasta.rb -b mapping/$i/Lib101_L1.bam -b mapping/$i/Lib101_L2.bam -b mapping/$i/Lib101_L3.bam -b mapping/$i/Lib101_L4.bam -p $j.$i.fasta.paf -l /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/$j.fasta"; done; done
```

## Tim methylation analysis


```bash
perl -e '$max = 0; $ctg = ""; while(<>){chomp; @s = split(/\t/); if($s[1] > $max){$max = $s[1]; $ctg = $s[0];}} print "$max\t$ctg\n";' < flye4.contigs.fasta.fai
5546585 contig_65409

cat sheep_poop_sequel1_CCS.flye4.contigs.ids sheep_poop_sequelII_CCS.flye4.contigs.ids | perl -lane 'if($F[1] eq "contig_65409"){print $F[0];}' > contig_65409.ccs.read.ids

# I want to make separate bins for read IDs for each contig. I'll also bin them into sizes that are > 5 mb, 1mb 100kb and 10kb, and separate them by read sets
mkdir sequel1contigreads
mkdir sequel2contigreads

# This should sort them all into the highest bin that they fit
perl -e '$odir = "sequel1contigreads"; @bins = (5000000, 1000000, 100000, 10000); foreach $i (@bins){mkdir("$odir/$i");} chomp(@ARGV); open(IN, "< $ARGV[0]"); %ctgs; while(<IN>){chomp; @s = split(/\t/); $ctgs{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $z = $ctgs{$s[1]}; $fb = -1; foreach $i (@bins){if($z >= $i){$fb = $i; last;}} if($fb == -1){next;} open(OUT, ">> $odir/$fb/$s[1].ids"); print OUT "$s[0]\n"; close OUT;} close IN;' flye4.contigs.fasta.fai sheep_poop_sequel1_CCS.flye4.contigs.ids

perl -e '$odir = "sequel2contigreads"; @bins = (5000000, 1000000, 100000, 10000); foreach $i (@bins){mkdir("$odir/$i");} chomp(@ARGV); open(IN, "< $ARGV[0]"); %ctgs; while(<IN>){chomp; @s = split(/\t/); $ctgs{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $z = $ctgs{$s[1]}; $fb = -1; foreach $i (@bins){if($z >= $i){$fb = $i; last;}} if($fb == -1){next;} open(OUT, ">> $odir/$fb/$s[1].ids"); print OUT "$s[0]\n"; close OUT;} close IN;' flye4.contigs.fasta.fai sheep_poop_sequelII_CCS.flye4.contigs.ids
```

#### TODO: Try aligning long reads and use SNP callers to identify variants. Stratify into potential strains based on graph of assembly.

```bash
module load samtools minimap2 

for i in flye3.contigs.fasta flye4.contigs.fasta clr1.contigs.fasta clr2.contigs.fasta clr3.contigs.fasta; do echo $i; sbatch -N 1 -n 3 -p priority -q msn --mem=55000 -t 3-0 --wrap="minimap2 -ax asm20 -R '@RG\tID:CCS\tSM:CCS' $i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel1_CCS.fasta /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta | samtools sort -T $i.tmp -o $i.ccs.bam -"; done

for i in clr1.contigs.fasta.ccs.bam clr2.contigs.fasta.ccs.bam clr3.contigs.fasta.ccs.bam flye4.contigs.fasta.ccs.bam; do echo $i; sbatch -N 1 -n 1 --mem=10000 -p priority -q msn --wrap="samtools index $i"; done

# Freebayes
module load miniconda/3.6
conda activate /lustre/project/forage_assemblies/sheep_project/complete_flye/freebayes

for i in flye3.contigs.fasta flye4.contigs.fasta clr1.contigs.fasta clr2.contigs.fasta clr3.contigs.fasta; do name=`echo $i | cut -d'.' -f1`; echo $name; sbatch run_freeparallel.sh $i $i.fai 100000 $i.ccs.bam $name.ccs.vcf; done

###### Desman analysis first
cp /lustre/project/rumen_longread_metagenome_assembly/analysis/desman/illumina_rest/run_desman_pipeline_illumina_rest.sh ./

mkdir desman/bin_lists
mkdir desman/bed_lists
for i in flye4 clr1 clr2 clr3; do mkdir desman/bin_lists/$i; mkdir desman/bed_lists/$i; done

for i in flye4 clr1 clr2 clr3; do echo $i; cat b3c_${i}_dastool/$i.das_proteins.faa.archaea.scg b3c_${i}_dastool/$i.das_proteins.faa.bacteria.scg > desman/$i.das_proteins.faa.faa.combined.scg; cat b3c_${i}_dastool/$i.das_proteins.faa | grep '>' | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' >  desman/$i.das_prodigal.shortform.tab ; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$r = $s[0]; $r =~ s/_\d{1,3}//; print "$h{$s[0]},$r,$s[1],$s[2],$s[0],$s[3]\n";}} close IN;' desman/$i.das_proteins.faa.faa.combined.scg desman/$i.das_prodigal.shortform.tab > desman/$i.das_prodigal.master_cogs.csv; done

# Now I've generated the master files for partitioning into separate bins. Let's generate the individual bins for parallelization

for i in flye4 clr1 clr2 clr3; do echo $i; cat b3c_${i}_dastool/$i.das_proteins.faa.archaea.scg b3c_${i}_dastool/$i.das_proteins.faa.bacteria.scg | perl -lane 'print $F[0];' > desman/$i.scg.list; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f desman/$i.das_prodigal.shortform.tab -c 0 -l desman/$i.scg.list | perl -lane '$r = $F[0]; $r  =~ s/_\d{1,3}$//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > desman/$i.scg.loc.bed; perl -lane 'open(OUT, ">> desman/bin_lists/$F[1].hqdas.bin.list"); print OUT "$F[0]"; close OUT;' < b3c_${i}_dastool/${i}.das_DASTool_scaffolds2bin.txt; mv desman/bin_lists/*.list desman/bin_lists/$i/; for j in desman/bin_lists/$i/*.list; do name=`basename $j | cut -d'.' -f1,2`; echo " $name"; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f desman/$i.scg.loc.bed -c 0 -l $j > desman/bed_lists/$i/${name}.scg.bed; done; done


# Now to cross my fingers and queue it up!

conda activate /KEEP/rumen_longread_metagenome_assembly/desman

# Note, I have to do this individually for each assembly as there are more than 1000 jobs in total
for i in flye4 clr1 clr2 clr3; do echo $i; mkdir desman/results_${i}; cd desman/results_${i}/; for j in ../bin_lists/$i/*.list; do name=`echo $j | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "bin_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $j | cut -d'.' -f1,2`; echo $bed; sbatch /lustre/project/forage_assemblies/sheep_project/complete_flye/run_desman_pipeline_illumina_rest.sh $name $bed $i; done; cd ../../; done

# Ugh, the virtual environments are letting me down!
conda activate /KEEP/rumen_longread_metagenome_assembly/desman
for i in flye4; do echo $i; mkdir desman/results_${i}; cd desman/results_${i}/; for j in ../bin_lists/$i/*.list; do name=`echo $j | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "bin_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $j | cut -d'.' -f1,2`; echo $bed; RUMEN=/lustre/project/rumen_longread_metagenome_assembly; FORAGE=/lustre/project/forage_assemblies/sheep_project/complete_flye; sbatch -N 1 -n 5 --mem=28000 -p priority -q msn --wrap=" python $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $FORAGE/${i}.contigs.fasta -c $FORAGE/desman/bin_lists/${i}/${bed}.hqdas.bin.list -g $FORAGE/desman/bed_lists/${i}/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -b $FORAGE/${i}.contigs.fasta.ccs.bam -o $name; python $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a $FORAGE/${i}.contigs.fasta -c $FORAGE/desman/bin_lists/${i}/${bed}.hqdas.bin.list -g $FORAGE/desman/bed_lists/${i}/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -o $name"; done; cd ../../; done

for i in clr1 clr2 clr3; do echo $i; mkdir desman/results_${i}; cd desman/results_${i}/; for j in ../bin_lists/$i/*.list; do name=`echo $j | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "bin_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $j | cut -d'.' -f1,2`; echo $bed; RUMEN=/lustre/project/rumen_longread_metagenome_assembly; FORAGE=/lustre/project/forage_assemblies/sheep_project/complete_flye; sbatch -N 1 -n 5 --mem=28000 -p priority -q msn --wrap=" python $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $FORAGE/${i}.contigs.fasta -c $FORAGE/desman/bin_lists/${i}/${bed}.hqdas.bin.list -g $FORAGE/desman/bed_lists/${i}/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -b $FORAGE/${i}.contigs.fasta.ccs.bam -o $name; python $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a $FORAGE/${i}.contigs.fasta -c $FORAGE/desman/bin_lists/${i}/${bed}.hqdas.bin.list -g $FORAGE/desman/bed_lists/${i}/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -o $name"; done; cd ../../; done

cat desman/results_flye4/*/desman_dic.fits.strain.count > desman/results_flye4/total_strain_counts.tab
for i in clr1 clr2 clr3; do cat desman/results_${i}/*/desman_dic.fits.strain.count > desman/results_${i}/total_strain_counts.tab; done

# Whatshap
conda activate /KEEP/rumen_longread_metagenome_assembly/whatshap


#### Testing
freebayes -f flye4.contigs.fasta -r contig_11399:1-228080 flye4.contigs.fasta.ccs.bam > flye4.contig_11399.example.vcf
perl -lane 'if($F[0] =~ /#/){print $_;}elsif($F[5] > 1){print $_;}' < flye4.contig_11399.example.vcf > flye4.contig_11399.filtered.vcf

freebayes -f flye4.contigs.fasta -r contig_46389:1-2593 flye4.contigs.fasta.ccs.bam > flye4.contig_4638.example.vcf
perl -lane 'if($F[0] =~ /#/){print $_;}elsif($F[5] > 1){print $_;}' < flye4.contig_4638.example.vcf > flye4.contig_4638.filtered.vcf
```

Creating a quick Desman strain count plot to show collaborators.

```R
# Created in perl
# perl -lane 'if($F[1] eq "None"){$F[1] = 0;} print "$F[0]\t$F[1]\tHiFi";' < desman/results_flye4/total_strain_counts.tab > desman/initial_strain_counts.tab
# perl -lane 'if($F[1] eq "None"){$F[1] = 0;} print "$F[0]\t$F[1]\tClr1";' < desman/results_clr1/total_strain_counts.tab >> desman/initial_strain_counts.tab

library(dplyr)
library(ggplot2)

sdata <- read.delim("desman/initial_strain_counts.tab", header=FALSE, sep="\t")
colnames(sdata) <- c("BName", "Strains", "ASM")

pdf(file="desman/initial_strain_counts.pdf", useDingbats=FALSE)
ggplot(sdata, aes(x=Strains, fill=ASM)) + geom_histogram(alpha=0.6, position="identity") + theme_bw() + ggtitle("DESMAN strain count CLR vs HIFI")
dev.off()

```


#### Viral association analysis

```bash
module load samtools/1.9 minimap2/2.6 miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/seaborn/

sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p priority -q msn -J flyeflye --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a flye4.contigs.fasta -b blobtools/table.flye4.contigs.blobDB.table.txt -i /lustre/project/forage_assemblies/sheep_project/complete_flye/mapping/flye4.contigs/hic_Sau3AI.bam -v flye4.viral.contigs.list -l /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta -m /software/7/apps/minimap2/2.6/minimap2 -o flye4.contigs.vassoc"

for i in clr1 clr2 clr3; do echo $i; perl -ne '@s = split(/\t/); if($s[29] eq "Viruses"){print "$s[0]\t$s[1]\n";}' < blobtools/table.${i}.contigs.blobDB.table.txt > ${i}.viral.contigs.list; done

for i in clr1 clr2 clr3; do echo $i; sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p priority -q msn -J ${i}_virus --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a ${i}.contigs.fasta -b blobtools/table.${i}.contigs.blobDB.table.txt -i /lustre/project/forage_assemblies/sheep_project/complete_flye/mapping/${i}.contigs/hic_Sau3AI.bam -v ${i}.viral.contigs.list -l /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta -m /software/7/apps/minimap2/2.6/minimap2 -o ${i}.contigs.vassoc"; done
```

That worked, but I tried to print a bipartite graph with no success. Much of the issue is that viral-viral associations are present. I want to test to see if I exclude these and sort the graph by viral class if I can solve this issue.


#### Rerunning das_tool

```bash
module load diamond/0.9.28 usearch/11.0.667 miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/das_tool/

perl -ne 'chomp; @s = split(/,/); print "$s[0]\tgb2.$s[1]\n";' < flye4.gb2.csvgraphbin2_wo_inconsistent_2.csv > flye4.gb2.csvgraphbin2_wo_inconsistent_2.tab

mkdir gb2_dastool; sbatch --nodes=1 --mem=25000 --ntasks-per-node=4 -p priority -q msn --wrap="DAS_Tool --search_engine 'diamond' -i flye4.gb2.csvgraphbin2_wo_inconsistent_2.tab,binning/bin3c/flye4.contigs/bin3c.full.clusters.tab,binning/metabat2/flye4.contigs/metabat2.full.clusters.tab -l graph2bin,bin3c,metabat -c flye4.contigs.fasta -o gb2_dastool/flye4.das -t 4 --write_bins 1"

perl -ne 'chomp; @s = split(/,/); print "$s[0]\tgb2.$s[1]\n";' < flye4.gb2.csvgraphbin2_wo_unsupported_1.csv > flye4.gb2.csvgraphbin2_wo_unsupported_1.tab

mkdir gb1_dastool; sbatch --nodes=1 --mem=25000 --ntasks-per-node=4 -p priority -q msn --wrap="DAS_Tool --search_engine 'diamond' -i flye4.gb2.csvgraphbin2_wo_unsupported_1.tab,binning/bin3c/flye4.contigs/bin3c.full.clusters.tab,binning/metabat2/flye4.contigs/metabat2.full.clusters.tab -l graph2bin1,bin3c,metabat -c flye4.contigs.fasta -o gb1_dastool/flye4.gb1.das -t 4 --write_bins 1"


### I need to rerun dastool on just the bin3c bins to get this ready for DEsman
for i in flye4 clr1 clr2 clr3; do echo $i; mkdir b3c_${i}_dastool; sbatch --nodes=1 --mem=25000 --ntasks-per-node=4 -p priority -q msn --wrap="DAS_Tool --search_engine 'diamond' -i binning/bin3c/${i}.contigs/bin3c.full.clusters.tab -l bin3c -c $i.contigs.fasta -o b3c_${i}_dastool/${i}.das -t 4 --write_bins 1"; done
```