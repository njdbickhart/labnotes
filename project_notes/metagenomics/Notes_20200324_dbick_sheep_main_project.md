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
```