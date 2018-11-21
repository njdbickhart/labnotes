# Rumen taxonomic analysis and major hypothesis testing
---
*7/20/2018*

These are my notes on testing the major hypotheses in the datasets. For starters, we want to know what is unique about the pacbio reads and assembly, and we want to know how much data is needed to finish the assembly.

## Table of contents


## Salivary microbiome presence in the rumen contents

Kevin's 2-way genera comparison venn containing the PacBio 16S long read data showed five genera that were not present in the amplicon data from the same sample. What are they, are they present in the assembly and what is unique about them that prevents them from detection in the amplicon data?

Here are the genera and brief descriptions of them:

#### PacBio 16S unique genera
| Genera | Description | PacbioCtgs? | IlluminaCtgs? |
| :--- | :--- | ---: | ---: |
|Mucinivorans | an anaerobic, mucin-degrading bacterium isolated from the digestive tract of the medicinal leech | 26 | 427 |
|Anaerorhabdus | Nonsaccharolytic; Acetic and lactic acids producer; Isolated from infected appendix, lung abscesses, and abdominal abscesses. Infrequently isolated from human and pig feces. | 0 | 0 |
|Aestuariispira | a Gram-negative, aerobic and non-motile bacteria from the genus of Aestuariispira which has been isolated from tidal flat; only known species | 0 | 0 |
|Coprobacter | Gram-negative obligate anaerobic bacterium belonging to the phylum Bacteroidetes. C. fastidiosus strain NSB1 isolated from human infant feces. | 0 | 0 |
|Microbacter | non-motile bacterium from the genus of Microbacter which has been isolated from sediments from the Tinto River |  5 | 182 |

OK, so allot of false positives (likely from 16S profiling from a faulty database?). Let's see what's present and what's not in the blob dbs from a simple venn comparison.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools

```bash
# Generating the lists
# NOTE: I'm removing any "undef" entries by default
perl -e '$prefix = "pacbio"; %taxa = ("phylum" => 27, "order" => 31, "family" => 35, "genus" => 39, "species" => 43); %files; foreach my $k (keys(%taxa)){open(my $TEMP, "> $prefix\_$k\.list"); $files{$k} = \*$TEMP;} while(<>){chomp; if($_ =~ /^\#/ || $_ =~ /\-undef/){next;}else{my @s = split(/\t/); foreach my $k (keys(%taxa)){my $idx = $taxa{$k}; my $fh = $files{$k}; print {$fh} "$s[$idx]\n";}}} foreach my $k (keys(%files)){close $files{$k};}' < pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt

perl -e '$prefix = "illumina"; %taxa = ("phylum" => 27, "order" => 31, "family" => 35, "genus" => 39, "species" => 43); %files; foreach my $k (keys(%taxa)){open(my $TEMP, "> $prefix\_$k\.list"); $files{$k} = \*$TEMP;} while(<>){chomp; if($_ =~ /^\#/ || $_ =~ /\-undef/){next;}else{my @s = split(/\t/); foreach my $k (keys(%taxa)){my $idx = $taxa{$k}; my $fh = $files{$k}; print {$fh} "$s[$idx]\n";}}} foreach my $k (keys(%files)){close $files{$k};}' < illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt

# Now for the entries
## PHYLUM ##
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_phylum.list pacbio_phylum.list
File Number 1: illumina_phylum.list
File Number 2: pacbio_phylum.list
Set     Count
1       17
1;2     46

## ORDER ##
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_order.list pacbio_order.list
File Number 1: illumina_order.list
File Number 2: pacbio_order.list
Set     Count
1       196
1;2     177

## FAMILY ##
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_family.list pacbio_family.list
File Number 1: illumina_family.list
File Number 2: pacbio_family.list
Set     Count
1       393
1;2     332

## GENUS ##
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_genus.list pacbio_genus.list
File Number 1: illumina_genus.list
File Number 2: pacbio_genus.list
Set     Count
1       1120
1;2     893
2       2

## SPECIES ##
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_species.list pacbio_species.list
File Number 1: illumina_species.list
File Number 2: pacbio_species.list
Set     Count
1       3517
1;2     2159
2       14

# I'm curious about the unique illumina entries
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o illumina_phylum.list pacbio_phylum.list
mv group_1.txt illumina_only_phylum.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt -l illumina_only_phylum.list -c 27 -i '#' -o illumina_only_phylum.lists.tab

# There are 1172 lines, only! Probably artifacts
# The only contigs above 10kb in length are from Dictyoglomales and Mollusca-undef phyla

# Now for the orders
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o illumina_order.list pacbio_order.list
mv group_1.txt illumina_only_order.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt -l illumina_only_order.list -c 31 -i '#' -o illumina_only_order.lists.tab

# That yielded only 4,627 novel lines, most of which were below 1500 bp in length

# How about the two genus entries from Pacbio that are unique?
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o illumina_genus.list pacbio_genus.list
mv group_2.txt pacbio_only_genus.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt -l pacbio_only_genus.list -c 39 -i '#' -o pacbio_only_genus.lists.tab

# Only two and they're small (8kb and 16kb) and two spurious assignments.
```

Oh well, it was good to check this. Let's see the pacbio contigs that have no illumina alignments next.

## PacBio unique sequence

Now I want to parse out the pacbio contigs that have no illumina alignments to see how they stack up against the rest of the assemblies.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools

```bash
# cum_sum less than the median and the USDA illumina coverage less than 1
perl -lane 'if($F[0] =~ /^#/){next;}elsif($F[22] <= 235 && $F[9] <=  1){print $_;}' < pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt | wc -l
75

perl -lane 'if($F[0] =~ /^#/){next;}elsif($F[22] <= 235 && $F[9] <=  1){print $_;}' < pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt | perl -lane 'print "$F[1]";' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   75
Minimum 1581
Maximum 31005
Average 7504.013333
Median  6265
Standard Deviation      5870.727140
Mode(Highest Distributed Value) 3380

# some of these are large contigs! Let's see what they look like; Allot of no-hits, but some definite taxa
# OK, let's see if any of these ended up in bins
perl -lane 'if($F[0] =~ /^#/){next;}elsif($F[22] <= 235 && $F[9] <=  1){print $_;}' < pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt | perl -lane 'print $F[0];' > pacbio_low_illumina_cov_ctgs.list

# Progressively adding stats to the tables
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-"){print "$F[1]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-" && $F[56] eq "-"){print "$F[1]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-" && $F[56] eq "-" && $F[57] eq "-"){print "$F[1]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
```
#### Pacbio Unique contig count table

| Combination | Count | Total BP | Min BP | Max BP | Avg BP | Median BP |
| :---------- | -----:| --------:|-------:|-------:|-------:|----------:|
|-Ilmn        | 796   | 4862643  | 1047   | 31005  | 6108   | 5776      |
|-Ilmn -RMG   | 682   | 4089505  | 1047   | 31005  | 5996   | 5700      |
|-Ilmn -RMG -H| 655   | 3919783  | 1047   | 31005  | 5984   | 5649      |

```bash
# I need comparisons, so let's do the same for the Illumina side
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-"){print "$F[1]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-" && $F[56] eq "-"){print "$F[1]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
perl -ne 'chomp; @F = split(/\t/); if($F[55] eq "-" && $F[56] eq "-" && $F[57] eq "-"){print "$F[1]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
``` 
#### Illumina Unique contig count table

| Combination |  Count |  Total BP | Min BP | Max BP | Avg BP | Median BP |
| :---------- | ------:| ---------:|-------:|-------:|-------:|----------:|
| -Pb         | 1474099| 2894091717| 1000   | 219586 | 1963   | 1481      |
| -Pb -RMG    | 1279779| 2458724722| 1000   | 219586 | 1921   | 1465      |
| -Pb -RMG -H | 1236195| 2363549480| 1000   | 219586 | 1911   | 1462      |

```bash
# calculating how much is not present in one previously published dataset vs another
perl -ne 'chomp; @F = split(/\t/); if($F[56] ne "-" && $F[1] > 2000){print "$F[0]\n";}' < illumina_megahit_master_table_2018_07_31.tab > illumina_megahit_present_inMICK.list
perl -ne 'chomp; @F = split(/\t/); if($F[57] ne "-" && $F[1] > 2000){print "$F[0]\n";}' < illumina_megahit_master_table_2018_07_31.tab > illumina_megahit_present_inHUN.list
perl -lane 'if($F[1] > 2000){print $F[0];}' < ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.gt2kb.fa.fai > illumina_megahit_gt2kb.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_megahit_gt2kb.list illumina_megahit_present_inMICK.list illumina_megahit_present_inHUN.list
File Number 1: illumina_megahit_gt2kb.list
File Number 2: illumina_megahit_present_inMICK.list
File Number 3: illumina_megahit_present_inHUN.list
Set     Count
1       334360
1;2     139951
1;2;3   73218
1;3     26853
2       35445
2;3     10962
3       7042

# Hmm, something fishy with that venn. Let's try my previously created list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa.fai -l illumina_megahit_no_other_db_aligns.gt15kb.ctgs.list -c 0 | perl -lane 'if($F[1] > 2000){print $F[1];}' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   327615
Sum:    1134301921
Minimum 2001
Maximum 219586
Average 3462.301546
Median  2747
Standard Deviation      2528.372692
Mode(Highest Distributed Value) 2027

# now for Pacbio
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai -l pacbio_final_no_other_db_aligns.ctgs.list -c 0 | perl -lane 'print $F[1];' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   655
Sum:    3919783
Minimum 1047
Maximum 31005
Average 5984.401527
Median  5649
Standard Deviation      3481.092459
Mode(Highest Distributed Value) 7263
```

## KAT plots for asm kmer comparisons

I should have done this far earlier, but I am going to generate some kmer comparisons for the new assembly and the reads. My hope is to pick up the same trend I noticed in the whisker plot to see what we're missing.

I already compiled 21 mer counts for each of read datasets (using the pacbio error corrected reads rather than the raw reads). 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina

```bash
module load jellyfish/2.2.3
sbatch --nodes=1 --mem=300000 --ntasks-per-node=20 -p assemble1 --wrap="jellyfish count -m 21 -s 100M -t 20 -C -o pacbio_final_asm_21mer ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa"

sbatch --nodes=1 --mem=300000 --ntasks-per-node=20 -p assemble1 --wrap="jellyfish count -m 21 -s 100M -t 20 -C -o illumina_megahit_asm_21mer ../illumina_usda_accumulated/mick_megahit_final_full.rfmt.fa"

# Let's try to generate a two way read plot first
sbatch --nodes=1 --mem=500000 --ntasks-per-node=25 -p assemble2 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/kat-2.3.4/src/kat comp -t 25 -o illumina_run3_pb_errorcorrect_21comp.kat illumina_run3_21mer pacbio_error_corrected_21mer"

# Now illumina kmer vs pacbio
export LD_LIBRARY_PATH=/mnt/nfs/nfs2/bickhart-users/binaries/anaconda3/lib:/mnt/nfs/nfs2/bickhart-users/binaries/KAT/deps/boost/build/lib:$LD_LIBRARY_PATH
sbatch --nodes=1 --mem=400000 --ntasks-per-node=20 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat comp -t 20 -o illumina_run3_pbpilon_asm 'YMPrepCannula_run3_L1_R1.fastq.gz YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz' ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa"

# Now illumina kmer vs illumina
sbatch --nodes=1 --mem=500000 --ntasks-per-node=20 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat comp -t 20 -o illumina_run3_illumina_asm 'YMPrepCannula_run3_L1_R1.fastq.gz YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz' ../illumina_usda_accumulated/mick_megahit_final_full.rfmt.fa"

# Because the default spectra plots are always undersized, let's resize them
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot spectra-cn -o illumina_run3_illumina_asm-extend.spectra-cn.png -t "Illumina Reads vs Illumina megahit ASM" -y 55000000 -x 150 illumina_run3_illumina_asm-main.mx
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot spectra-cn -o illumina_run3_pacbio_asm-extend.spectra-cn.png -t "Illumina Reads vs PacBio Pilon ASM" -y 90000000 -x 150 illumina_run3_pbpilon_asm-main.mx

# And with minimum cov filters
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot spectra-cn -o illumina_run3_pacbio_asm-mincov.spectra-cn.png -t "Illumina Reads vs PacBio Pilon ASM" -y 90000000 -x 150 -i 1 illumina_run3_pbpilon_asm-main.mx
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot spectra-cn -o illumina_run3_illumina_asm-mincov.spectra-cn.png -t "Illumina Reads vs Illumina megahit ASM" -y 55000000 -x 150 -i 1 illumina_run3_illumina_asm-main.mx

# Density analysis
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot density -o illumina_run3_illumina_asm.density.png -x 150 illumina_run3_illumina_asm-main.mx
/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat plot density -o illumina_run3_pbpilon_asm.density.png -x 150 illumina_run3_pbpilon_asm-main.mx

# Now pacbio kmer vs illumina
sbatch --nodes=1 --mem=500000 --ntasks-per-node=20 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat comp -t 20 -o pacbio_ecreads_illumina_asm ../error_corrected_reads/rumen_pacbio_corrected.fasta.gz ../illumina_usda_accumulated/mick_megahit_final_full.rfmt.fa"
# And pacbio kmer vs pacbio
sbatch --nodes=1 --mem=500000 --ntasks-per-node=20 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat comp -t 20 -o pacbio_ecreads_pbpilon_asm ../error_corrected_reads/rumen_pacbio_corrected.fasta.gz ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa"

# Finally, generating kat sect plots for everything
sbatch --nodes=1 --mem=500000 --ntasks-per-node=40 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat sect -o illumina_run3_pbpilon_asm-sect -t 40 -E ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz"
sbatch --nodes=1 --mem=500000 --ntasks-per-node=40 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat sect -o illumina_run3_illumina_asm-sect -t 40 -E ../illumina_usda_accumulated/mick_megahit_final_full.rfmt.fa YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz"
```

## TODO: Generate KAT plots for everything
## TODO: Also, combine the kmerspectra plots for Illumina and Pacbio
## TODO: Compute rand index for bins and check bin fidelity for pacbio-unique contigs with no illumina mappings
## TODO: look into viral genomes to estimate size, frequency and bias in pacbio and illumina datasets 

## Kmer complexity analysis

I want to demonstrate the complexity of our dataset once and for all. Let's compare to the Stewart et al. data, PRJEB21624. I'll use a bloom filter to remove the single frequency kmers.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/kmer

```bash
module load jellyfish2/2.2.9
grep 'PRJEB21624' $RUMEN/sequence_data/public_datasets/cattle_rumen_illumina_data_spreadsheet.revised.tab | perl -lane 'print "gunzip -c $F[0]\ngunzip -c $F[1]";' > PRJEB21624.generators

sbatch --nodes=1 --mem=350000 --ntasks-per-node=20 -p mem --wrap="jellyfish count -m 21 -s 100M -t 20 --bf-size 100G -C -g PRJEB21624.generators -o PRJEB21624_21mer_high"
# That failed! Trying higher mem
sbatch --nodes=1 --mem=650000 --ntasks-per-node=20 -p mem --wrap="jellyfish count -m 21 -s 100M -t 20 --bf-size 100G -C -g PRJEB21624.generators -o PRJEB21624_21mer_high"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="jellyfish stats -o PRJEB21624_21mer_high.stats PRJEB21624_21mer_high"

```


## Viral genomes

I want to fully characterize the viruses and see how they pan out with respect to host specificity and kmer composition.

Here are the specific types of analysis I want to perform:
1. alignment of error corrected pacbio reads to identify "junctions" with other contigs
2. the difference in kmers between the two sets of viral contigs
3. the difference in read counts between viral contigs

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools

```bash
srun --nodes=1 --ntasks-per-node=4 --mem=15000 -p short --pty bash

module load samtools bwa
grep 'Viruses' illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt | perl -lane 'system("samtools faidx ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.fa $F[0] >> illumina_megahit_viruses.fa");'

grep 'Viruses' pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt | perl -lane 'system("samtools faidx ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa $F[0] >> pacbio_pilon_viruses.fa");'

module load minimap2/2.6
sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="minimap2 -x map-pb pacbio_pilon_viruses.fa ~/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta.gz > pacbio_pilon_viruses_ecpbreads.paf"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="minimap2 -x map-pb illumina_megahit_viruses.fa ~/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta.gz > illumina_megahit_viruses_ecpbreads.paf"

# Number of error corrected reads that map to the tail ends of viral contigs and may impart some host specificity data
perl -lane 'if($F[9] > 500 && ($F[7] < 200 || $F[6] - $F[8] < 200)){print $_;}' < pacbio_pilon_viruses_ecpbreads.paf | wc -l
12406

# Far fewer, suggesting that there's either a chimerism thing with the pacbio data or that the illumina contigs are somehow off
perl -lane 'if($F[9] > 500 && ($F[7] < 200 || $F[6] - $F[8] < 200)){print $_;}' < illumina_megahit_viruses_ecpbreads.paf | wc -l
9422

# I wrote a custom script to try to filter and organize the data
### PACBIO FIRST ###
perl selectLikelyViralOverhangs.pl pacbio_pilon_viruses_ecpbreads.paf pacbio_pilon_viruses_ecpbreads.filt

# Hmm... most of the alignments were smaller than expected. Let's keep things that are larger than 150 bp so that we can get a good alignment on other contigs
perl -lane 'if($F[2] - $F[1] > 150){print $_;}' < pacbio_pilon_viruses_ecpbreads.filt.subread.bed > pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.bed
# Now to generate a new fasta and align to the whole PacBio dataset
samtools faidx ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta

perl -e '@list; while(<>){chomp; @s = split(/\t/); push(@list, "$s[0]:$s[1]-$s[2]"); if(scalar(@list) > 500){print "Printing...\n"; system("samtools faidx ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta " . join(" ", @list) . " >> pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.fa"); @list = ();}} system("samtools faidx ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta " . join(" ", @list) . " >> pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.fa");' < pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.bed

sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="minimap2 -x map-pb ../../assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.fa > pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.paf"

perl -lane 'print $F[11];' < pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.paf | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   30443
Sum:    337001
Minimum 0
Maximum 60
Average 11.069901
Median  0
Standard Deviation      20.821972
Mode(Highest Distributed Value) 0

# Note: the first contig column is the non-viral contig
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); if($s[11] == 0){next;} $s[0] =~ s/\:\d+\-\d+$//; push(@{$data{$s[0]}}, $s[5]);} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($data{$s[0]})){push(@{$data{$s[0]}}, $s[5]);}} close IN; foreach my $k (keys(%data)){print "$k\t" . join("\t", @{$data{$k}}) . "\n";}' pacbio_pilon_viruses_ecpbreads.filt.subread.gt150.paf pacbio_pilon_viruses_ecpbreads.paf > pacbio_pilon_viruses_ecpbreads.assoc.filt.tab

# Stringent filtering of association sites
perl -lane 'if(scalar(@F) > 3){next;}else{print $_;}' < pacbio_pilon_viruses_ecpbreads.assoc.filt.tab > pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.tab

# I scripted a means of condensing down all of these associations into a summary table
perl generateViralAssociationGraph.pl pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.tab pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab

# There were 183 associations! Of these, 124 associations were greater than one observation. 
perl -lane 'if($F[2] > 1){print $_;}' < pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab | perl -alne 'print $F[0];' | sort | uniq | wc -l
45  <- unique viral contigs with potential host data!

# I am now downloading this tab file for viewing in cytoscape and eventual conversion into a circos plot if I see allot of co-infectivity.


#### Illumina comparison ####
perl selectLikelyViralOverhangs.pl illumina_megahit_viruses_ecpbreads.paf illumina_megahit_viruses_ecpbreads.filt
Finished. Identified 3403 errors
perl -lane 'if($F[2] - $F[1] > 150){print $_;}' < illumina_megahit_viruses_ecpbreads.filt.subread.bed > illumina_megahit_viruses_ecpbreads.filt.subread.gt150.bed
perl -e '@list; while(<>){chomp; @s = split(/\t/); push(@list, "$s[0]:$s[1]-$s[2]"); if(scalar(@list) > 500){print "Printing...\n"; system("samtools faidx ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta " . join(" ", @list) . " >> illumina_megahit_viruses_ecpbreads.filt.subread.gt150.fa"); @list = ();}} system("samtools faidx ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta " . join(" ", @list) . " >> illumina_megahit_viruses_ecpbreads.filt.subread.gt150.fa");' < illumina_megahit_viruses_ecpbreads.filt.subread.gt150.bed
sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="minimap2 -x map-pb ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa illumina_megahit_viruses_ecpbreads.filt.subread.gt150.fa > illumina_megahit_viruses_ecpbreads.filt.subread.gt150.paf"

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); if($s[11] == 0){next;} $s[0] =~ s/\:\d+\-\d+$//; push(@{$data{$s[0]}}, $s[5]);} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($data{$s[0]})){push(@{$data{$s[0]}}, $s[5]);}} close IN; foreach my $k (keys(%data)){print "$k\t" . join("\t", @{$data{$k}}) . "\n";}' illumina_megahit_viruses_ecpbreads.filt.subread.gt150.paf illumina_megahit_viruses_ecpbreads.paf > illumina_megahit_viruses_ecpbreads.assoc.filt.tab
perl -lane 'if(scalar(@F) > 3){next;}else{print $_;}' < illumina_megahit_viruses_ecpbreads.assoc.filt.tab > illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.tab

# And now for the Cytoscape tab file
perl generateViralAssociationGraph.pl illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.tab illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab


# General stats on association:
wc -l *.cyto.tab
  223 illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab
  184 pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab

for i in *.cyto.tab; do echo $i; perl -lane 'if($F[2] > 1){print $_;}' < $i | wc -l; done
illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab
73
pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab
124

for i in *.cyto.tab; do echo $i; perl -lane 'if($F[2] > 1){print $F[0];}' < $i | sort | uniq | wc -l; done
illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab
64
pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab
45

# So the illumina data shows far fewer edges between host-virus integration


#### Preparing files for master tables ####
perl -lane 'if($F[0] eq "VirusCtg"){next;}elsif($F[2] < 2){next;} print "$F[1]\t$F[0]\t$F[3]";' < pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab > pacbio_pilon_viruses.host.assoc.tab
perl -lane 'if($F[0] eq "VirusCtg"){next;}elsif($F[2] < 2){next;} print "$F[1]\t$F[0]\t$F[3]";' < illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab > illumina_megahit_viruses.host.assoc.tab
```

## Alignment to other generated rumen microbial datasets

I want to see how much of our assembly space is already occupied by other assemblies for the rumen.

#### Hungate 1000

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate

```bash
for i in assembly_fasta/*.fasta.gz; do echo $i; gunzip -c $i >> hungate_combined_unordered_reference.fasta; done

# Damn, the hungate fastas have redundant scaffold names! I need to automate this differently
perl -e '@files = `ls ./assembly_fasta/*.gz`; chomp(@files); foreach $f (@files){@fname = split(/[\/\.\-]/, $f); my $use = $fname[3]; open(IN, "gunzip -c $f |"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//g; print ">$use\.$_\n";}else{print $_;}} close IN; }' > hungate_combined_unordered_reference.fasta

module load minimap2/2.6
sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 hungate_combined_unordered_reference.fasta ../pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > hungate_vs_pb_pilon.paf"

## Quick test to see if they assembled any prophage that we observe in our dataset
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f hungate_vs_pb_pilon.paf -c 0 -l /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools/pacbio_pilon_viruses_ecpbreads.list
# Nothing! Wow!

sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 hungate_combined_unordered_reference.fasta ../pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa > hungate_vs_ilmn_megahit.paf"

# and the association tab files for later use in the master table
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < hungate_vs_pb_pilon.paf | perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' > hungate_vs_pb_pilon.assoc.condensed.tab
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < hungate_vs_ilmn_megahit.paf | perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' > hungate_vs_ilmn_megahit.assoc.condensed.tab
```

#### Mick's RUG

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug

```bash
for i in genomes/*.fa; do echo $i; cat $i >> mick_combined_900_rug.fa; done

sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 mick_combined_900_rug.fa ../pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > mick_vs_pb_pilon.paf"

sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 mick_combined_900_rug.fa ../pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa > mick_vs_ilmn_megahit.paf"

# Now to generate simple association files for my master table
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < mick_vs_ilmn_megahit.paf > mick_vs_ilmn_megahit.assoc.tab
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < mick_vs_pb_pilon.paf > mick_vs_pb_pilon.assoc.tab

# Concatenating the multiplehit entries by semi-colon
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' < mick_vs_pb_pilon.assoc.tab > mick_vs_pb_pilon.assoc.condensed.tab
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' < mick_vs_ilmn_megahit.assoc.tab > mick_vs_ilmn_megahit.assoc.condensed.tab

### Redoing binning ###
# Because my previous bin assignments were confusing and inconsistent, I'm going to update the master tables with new bin assignments from my list

perl -lane '@b = split(/[_\.]/, $F[1]); print "$F[0]\t$b[-2]\_$b[-1]";' < illumina_dastool_analysis_binset_lt10redund.bins > illumina_dastool_analysis_binset_lt10redund.tab
perl -lane '@b = split(/[_\.]/, $F[1]); print "$F[0]\t$b[-2]\_$b[-1]";' < illumina_dastool_high_quality_dasbins.contigs > illumina_dastool_high_quality_dasbins.contigs.tab

perl -lane '@b = split(/[_\.]/, $F[1]); print "$F[0]\t$b[-2]\_$b[-1]";' < pacbio_dastool_analysis_binset_lt10redund.bins > pacbio_dastool_analysis_binset_lt10redund.tab
perl -lane '@b = split(/[_\.]/, $F[1]); print "$F[0]\t$b[-2]\_$b[-1]";' < pacbio_dastool_high_quality_dasbins.contigs > pacbio_dastool_high_quality_dasbins.contigs.tab
```


## Concatenating everything into master tables.

My goal is to generate "master tables" of all known data about each contig. For each assembly, this includes the following information.

##### Essential table data
* Contig name
* Length
* GC percent
* Alignment data
	* From 17 datasets including the YMprep run3 illumina data
	* Other historic rumen WGS data
* Taxonomic assignment
* General functional categorization (ie. FAPROTAX) 
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools/illumina_megahit_contigs_functional_assign.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools/pacbio_secpilon_contigs_functional_assign.tab*
* Binning results
	* Metabat
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/illumina_megahit_public_metabat.unsorted.bins*
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/pacbio_final_public_metabat.unsorted.bins*
	* Hi-C
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/illumina_megahit_hic.unsorted.bins*
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/pacbio_final_public_hic.unsorted.bins*
	* DAS_tool less than 10 SCG redundancy
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/illumina_dastool_analysis_binset_lt10redund.tab*
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/pacbio_dastool_analysis_binset_lt10redund.tab*
	* DAS_tool high quality draft (< 5 redundancy; > 80 completeness)
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/illumina_dastool_high_quality_dasbins.contigs.tab*
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/pacbio_dastool_high_quality_dasbins.contigs.tab*
* Number of predicted ORFs
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal/illumina_megahit_prodigal_proteins.count.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal/pacbio_final_prodigal_proteins.count.tab*
* Contains CRISPR array **TODO**
* Non-viral host with association to viral contig
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools/illumina_megahit_viruses.host.assoc.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools/pacbio_pilon_viruses.host.assoc.tab*
* Analog alignment in the opposite assembly
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_vs_pacbio_minimap2.assoc.condensed.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_vs_illumina_minimap2.assoc.condensed.tab*
* Alignment to Mick's RMGs
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_vs_ilmn_megahit.assoc.condensed.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_vs_pb_pilon.assoc.condensed.tab*
* Alignment to Hungate 1000
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_vs_ilmn_megahit.assoc.condensed.tab*
	* */home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_vs_pb_pilon.assoc.condensed.tab*

Blobtools has already generated tabular data on the rumen microbial data up to the taxonomic assignment. I want to generate association data on all of the other tabs and then join them all into the larger dataset. 


#### FAPROTAX association

OK, I need to format the data so that it is suitable for analysis by FAPROTAX. There is a format conversion from BIOM in the QIME package but my blobplots data is not directly convertible. Custom scripting to the rescue!

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools

```bash
perl convert_blob_to_otu_table.pl -b illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt -l 21 -c 23,27,31,35,39,43 -o illumina_blobplot_mock_otu_table.tab

perl convert_blob_to_otu_table.pl -b pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt -l 21 -c 23,27,31,35,39,43 -o pacbio_blobplot_mock_otu_table.tab
```

> Assembler2: /mnt/nfs/nfs2/bickhart-users/binaries/FAPROTAX_1.1

```bash
python2.7 collapse_table.py -i illumina_blobplot_mock_otu_table.tab -g FAPROTAX.txt -c '#' -d 'taxonomy' --omit_columns 0 --column_names_are_in last_comment_line -r illumina_report.txt -n columns_after_collapsing -v -o illumina_functional_table.tsv
  Reading input table..
  Loaded 8349 out of 8349 rows amongst 8350 lines, and 18 columns, from file 'illumina_blobplot_mock_otu_table.tab'
  Reading groups..
  Read 8563 lines from file 'FAPROTAX.txt', found 90 groups with 7820 members (4724 unique members)
  Assigning rows to groups..
  Collapsing table..
  Note: No entry found for group 'arsenite_oxidation_energy_yielding'
  Note: No entry found for group 'chlorate_reducers'
  Note: No entry found for group 'chloroplasts'
  Assigned 3571 records to groups, 4778 records were leftovers
  Writing report 'illumina_report.txt'..

python2.7 collapse_table.py -i pacbio_blobplot_mock_otu_table.tab -g FAPROTAX.txt -c '#' -d 'taxonomy' --omit_columns 0 --column_names_are_in last_comment_line -r pacbio_report.txt -n columns_after_collapsing -v -o pacbio_functional_table.tsv
  Reading input table..
  Loaded 2952 out of 2952 rows amongst 2953 lines, and 18 columns, from file 'pacbio_blobplot_mock_otu_table.tab'
  Reading groups..
  Read 8563 lines from file 'FAPROTAX.txt', found 90 groups with 7820 members (4724 unique members)
  Assigning rows to groups..
  Collapsing table..
  Note: No entry found for group 'arsenite_oxidation_detoxification'
  Note: No entry found for group 'arsenite_oxidation_energy_yielding'
  Note: No entry found for group 'dissimilatory_arsenite_oxidation'
  Note: No entry found for group 'aliphatic_non_methane_hydrocarbon_degradation'
  Note: No entry found for group 'chlorate_reducers'
  Note: No entry found for group 'chloroplasts'
  Note: No entry found for group 'anoxygenic_photoautotrophy_Fe_oxidizing'
  Assigned 1378 records to groups, 1574 records were leftovers
  Writing collapsed table 'pacbio_functional_table.tsv'..
  Writing report 'pacbio_report.txt'..
```

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/blobtools

```bash
# Now to convert those report.txt files into a direct association of contig with metabolic record
perl -e '$/ = "\#"; %h; while(<>){@lines = split(/\n/); $lines[0] =~ s/\s+//g; $lines[0] =~ s/\(.+\)\://;  for($x = 1; $x < scalar(@lines); $x++){$lines[$x] =~ s/\s+//g; if(length($lines[$x]) < 2){next;}push(@{$h{$lines[0]}}, $lines[$x]);}} foreach my $k (keys(%h)){foreach my $r (@{$h{$k}}){print "$r\t$k\n";}}' < illumina_report_shortened.txt > illumina_tax_functional_association.tab

perl -e '$/ = "\#"; %h; while(<>){@lines = split(/\n/); $lines[0] =~ s/\s+//g; $lines[0] =~ s/\(.+\)\://;  for($x = 1; $x < scalar(@lines); $x++){$lines[$x] =~ s/\s+//g; if(length($lines[$x]) < 2){next;}push(@{$h{$lines[0]}}, $lines[$x]);}} foreach my $k (keys(%h)){foreach my $r (@{$h{$k}}){print "$r\t$k\n";}}' < pacbio_shortened_report.txt > pacbio_tax_functional_association.tab

# Finally, process the reads the same way, but reverse associate the previous taxonomic assignments into a simple tab delimited list
perl generate_blob_tax_functional_tab.pl -b illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt -c 23,27,31,35,39,43 -o illumina_megahit_contigs_functional_assign.tab -t illumina_tax_functional_association.tab

perl generate_blob_tax_functional_tab.pl -b pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt -c 23,27,31,35,39,43 -o pacbio_secpilon_contigs_functional_assign.tab -t pacbio_tax_functional_association.tab
```

#### Final data associations for master table

##### Accumulating orf prediction totals

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal 

```bash
# Column 2 is the count of full protein orfs, column 3 is the count of partial protein orfs
gunzip -c ../prodigal/illumina_megahit_prodigal_proteins.shortform.tab.gz | perl -e '<>; %p; %h; while(<>){chomp; @s = split(/\t/); @base = split(/_/, $s[0]); pop(@base); $name = join("_", @base); if(!exists($h{$name})){$h{$name} = 0;} if(!exists($p{$name})){$p{$name} = 0;} if($s[5] =~ /00/){$h{$name} += 1;}else{$p{$name} += 1;}} foreach my $ctg (sort {$a cmp $b} keys(%h)){print "$ctg\t$h{$ctg}\t$p{$ctg}\n";}' > illumina_megahit_prodigal_proteins.count.tab

gunzip -c pacbio_final_prodigal_proteins.shortform.tab.gz | perl -e '<>; %p; %h; while(<>){chomp; @s = split(/\t/); @base = split(/_/, $s[0]); pop(@base); $name = join("_", @base); if(!exists($h{$name})){$h{$name} = 0;} if(!exists($p{$name})){$p{$name} = 0;} if($s[5] =~ /00/){$h{$name} += 1;}else{$p{$name} += 1;}} foreach my $ctg (sort {$a cmp $b} keys(%h)){print "$ctg\t$h{$ctg}\t$p{$ctg}\n";}' > pacbio_final_prodigal_proteins.count.tab

# Simple stat collection on full protein orfs
## PACBIO ##
cat pacbio_final_prodigal_proteins.count.tab | cut -f2 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   77664
Sum:    1191681
Minimum 0
Maximum 663
Average 15.344059
Median  12
Standard Deviation      16.430460
Mode(Highest Distributed Value) 8

## Illumina ##
cat illumina_megahit_prodigal_proteins.count.tab | cut -f2 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   2178467
Sum:    3082127
Minimum 0
Maximum 351
Average 1.414815
Median  1
Standard Deviation      3.011941
Mode(Highest Distributed Value) 0

```

##### Alignment to other assemblies

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project

```bash
module load minimap2/2.6

# Pacbio as reference
sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa illumina_megahit/illumina_megahit_final_contigs.perl.fa > illumina_vs_pacbio_minimap2.paf"

# Illumina as reference
sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 illumina_megahit/illumina_megahit_final_contigs.perl.fa pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > pacbio_vs_illumina_minimap2.paf"

# And the association tab file generation
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < illumina_vs_pacbio_minimap2.paf | perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' > illumina_vs_pacbio_minimap2.assoc.condensed.tab
perl -lane 'if($F[11] > 0 && $F[10] > 100){print "$F[0]\t$F[5]";}' < pacbio_vs_illumina_minimap2.paf | perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]);} foreach my $k (keys(%h)){print "$k\t" . join(";", @{$h{$k}}) . "\n";}' > pacbio_vs_illumina_minimap2.assoc.condensed.tab
```

#### Generation of master table files using JSON python script

OK, I have tabulated nearly every extraneous data file and now I can reformat the tables. First, let's preprocess the blobtools data tables before adding the data in.

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
cp ../blobtools/illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt ./illumina_preliminary_blobtools_table.tab
cp ../blobtools/pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt ./pacbio_preliminary_blobtools_table.tab

# I used vim to remove the superceeding comment lines at the top of the file

# OK, after some fumbling through Python's non-intuitive json interface, I think I managed to concatenate everything!
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j pacbio_data_files.json -t pacbio_preliminary_blobtools_table.tab -o pacbio_final_pilon_master_table_2018_07_31.tab

python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j illumina_data_files.json -t illumina_preliminary_blobtools_table.tab -o illumina_megahit_master_table_2018_07_31.tab

# Redoing master tables to get better DAS_tool bins incorporated
cp illumina_data_files.json illumina_data_files_9_2018.json
cp pacbio_data_files.json pacbio_data_files_9_2018.json

python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j pacbio_data_files_9_2018.json -t pacbio_preliminary_blobtools_table.tab -o pacbio_final_pilon_master_table_2018_09_07.tab
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j illumina_data_files_9_2018.json -t illumina_preliminary_blobtools_table.tab -o illumina_megahit_master_table_2018_09_07.tab

# And adding the final hic bins to the list
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j illumina_data_files_11_2018.json -t illumina_preliminary_blobtools_table.tab -o illumina_megahit_master_table_2018_11_15.tab
```


#### Gathering stats on major elements of the master tables.

> CERES: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```R
data <- read.delim("pacbio_final_pilon_master_table_2018_07_31.tab", header=TRUE)
summary(data)

# Here are the noteworthy stats in full tables
          name           length             GC               N
 tig00000002:    1   Min.   :  1044   Min.   :0.1544   Min.   :0.00e+00
 tig00000003:    1   1st Qu.:  6622   1st Qu.:0.4313   1st Qu.:0.00e+00
 tig00000004:    1   Median :  9916   Median :0.4898   Median :0.00e+00
 tig00000005:    1   Mean   : 13859   Mean   :0.4795   Mean   :2.57e-05
 tig00000009:    1   3rd Qu.: 15798   3rd Qu.:0.5363   3rd Qu.:0.00e+00
 tig00000012:    1   Max.   :710205   Max.   :0.9799   Max.   :1.00e+00
 (Other)    :77664

      cov6                cov7             cov_sum         superkingdom.t.24
 Min.   :   0.0000   Min.   :   0.000   Min.   :     0.8   Archaea  :  943
 1st Qu.:   0.0000   1st Qu.:   0.149   1st Qu.:    99.6   Bacteria :72099
 Median :   0.0000   Median :   0.479   Median :   235.0   Eukaryota:  757
 Mean   :   0.0344   Mean   :  14.256   Mean   :   460.8   no-hit   : 3755
 3rd Qu.:   0.0000   3rd Qu.:   2.549   3rd Qu.:   535.1   Viruses  :  116
 Max.   :1455.9300   Max.   :7501.035   Max.   :343451.0

                                                                 TaxFuncCategory
 N/A                                                                     :33914
 animal_parasites_or_symbionts                                           :16891
 fermentation;chemoheterotrophy                                          :11710
 fermentation;chemoheterotrophy;xylanolysis;animal_parasites_or_symbionts: 3583
 fermentation;chemoheterotrophy;animal_parasites_or_symbionts            : 1953
 fermentation;chemoheterotrophy;cellulolysis                             : 1267
 (Other)                                                                 : 8352
   MetabatBin        HiCBin          DASBin       CompleteORFs
 NOBIN  :34496   NOBIN  :14599   NOBIN  :68760   Min.   :  0.00
 292    : 1348   215    :  410   504    :  320   1st Qu.:  7.00
 632    :  717   173    :  312   194    :  283   Median : 12.00
 365    :  681   186    :  310   269    :  261   Mean   : 15.34
 108    :  542   86     :  310   542    :  234   3rd Qu.: 18.00
 188    :  457   128    :  305   372    :  229   Max.   :663.00
 (Other):39429   (Other):61424   (Other): 7583


# And the illumina data
data <- read.delim("illumina_megahit_master_table_2018_07_31.tab", header=TRUE)
# Noteworthy stats again:
           name             length             GC               N
 k127_100    :      1   Min.   :  1000   Min.   :0.0573   Min.   :0
 k127_1000   :      1   1st Qu.:  1216   1st Qu.:0.3108   1st Qu.:0
 k127_10000  :      1   Median :  1583   Median :0.4201   Median :0
 k127_1000000:      1   Mean   :  2342   Mean   :0.4083   Mean   :0
 k127_1000002:      1   3rd Qu.:  2415   3rd Qu.:0.5060   3rd Qu.:0
 k127_1000003:      1   Max.   :320393   Max.   :0.8146   Max.   :0
 (Other)     :2182257

      cov6               cov7             cov_sum         superkingdom.t.24
 Min.   :    0.00   Min.   :    0.00   Min.   :     2.6   Archaea  :  20620
 1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    21.5   Bacteria :1503853
 Median :    0.00   Median :    0.03   Median :    53.0   Eukaryota:  57227
 Mean   :    0.01   Mean   :    4.12   Mean   :   163.2   no-hit   : 598994
 3rd Qu.:    0.00   3rd Qu.:    0.11   3rd Qu.:   142.9   undef    :     32
 Max.   :22747.03   Max.   :38175.30   Max.   :373395.2   Viruses  :   1537

                                                                 TaxFuncCategory
 N/A                                                                     :1365778
 fermentation;chemoheterotrophy                                          : 251542
 animal_parasites_or_symbionts                                           : 249208
 fermentation;chemoheterotrophy;xylanolysis;animal_parasites_or_symbionts:  46614
 chemoheterotrophy;aerobic_chemoheterotrophy                             :  36805
 fermentation;chemoheterotrophy;xylanolysis;cellulolysis                 :  34543
 (Other)                                                                 : 197773

   MetabatBin          HiCBin            DASBin         CompleteORFs
 NOBIN  :1795096   NOBIN  :1998861   NOBIN  :2133769   Min.   :  0.000
 1393   :  32608   1145   :   1426   204    :   1392   1st Qu.:  0.000
 659    :  17345   16     :   1392   450    :   1011   Median :  1.000
 164    :  14025   1117   :   1306   1130   :    691   Mean   :  1.412
 1256   :   8064   1288   :   1299   177    :    690   3rd Qu.:  2.000
 634    :   7627   2285   :   1278   1096   :    673   Max.   :351.000
 (Other): 307498   (Other): 176701   (Other):  44037
  PartialORFs          ViralContig              ViralClass
 Min.   :0.000   -           :2182190   -            :2182190
 1st Qu.:1.000   k127_176453 :      3   Phikzvirus   :     32
 Median :2.000   k127_1792982:      3   C5virus      :      8
 Mean   :1.439   k127_2201563:      3   Cp8virus     :      6
 3rd Qu.:2.000   k127_162296 :      2   Schizot4virus:      6
 Max.   :2.000   k127_75260  :      2   N4virus      :      5
                 (Other)     :     60   (Other)      :     16

```

#### PacBio assembly Master Table Noteworthy stats

| Attribute | Value | Description |
| :--- | :--- | :--- |
| Max Len | 710,205 | Max Contig len |
| Median Len | 9,916 | Median Contig len |
| Tax No-hits | 3,755 | Contigs with no Diamond taxonomic hit assignments |
| FAPROTAX No-hits | 33,914 | Contigs with no FAPROTAX metabolic summaries|
| Metabat NoBIN | 34,496 | Contigs with no assigned Metabat bins |
| Hi-C NoBIN | 14,599 | Contigs with no assigned Hi-C bins |
| DAStool NoBin| 68,760 | Contigs with no assigned DASTools bin |
| Max ORF count | 663 | The maximum number of complete ORFs found in one contig |
| Not found in Ilmn | 796 | The count of contigs with no large alignments in the Illumina dataset |
| Not found in MICK | 25,056 | Contigs with no large alignments in the MICK RMG dataset |
| Not found in Hungate | 45,389 | Contigs with no large alignments in the Hungate 1000 dataset|

#### Illumina assembly Master Table Noteworthy stats

| Attribute | Value | Description |
| :--- | :--- | :--- |
| Max Len | 320,393 | Max Contig len |
| Median Len | 1,583 | Median Contig len |
| Tax No-hits | 598,994 | Contigs with no Diamond taxonomic hit assignments |
| FAPROTAX No-hits | 1,365,778 | Contigs with no FAPROTAX metabolic summaries|
| Metabat NoBIN | 1,998,861 | Contigs with no assigned Metabat bins |
| Hi-C NoBIN | 1,998,861 | Contigs with no assigned Hi-C bins |
| DAStool NoBin| 2,133,769 | Contigs with no assigned DASTools bin |
| Max ORF count | 351 | The maximum number of complete ORFs found in one contig |
| Not found in Pacbio | 1,474,099 | The count of contigs with no large alignments in the Pacbio dataset |
| Not found in MICK | 1,624,079 | Contigs with no large alignments in the MICK RMG dataset |
| Not found in Hungate | 1,955,529 | Contigs with no large alignments in the Hungate 1000 dataset|

#### Finding additional values for the manuscript

This is a collection of methods for me to try to find very obvious associations in the data for use in the manuscript.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/master_tables

```bash
# Simple trick to find index of columns for association
head -n 1 illumina_megahit_master_table_2018_07_31.tab | perl -lane 'for($x = 0; $x < scalar(@F); $x++){print "$F[$x]\t$x";}'

# First, DASTool bins compared to previously observed datasets
# Vs Micks
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[56]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less
# All accounted for in Mick's dataset

# Vs Hungate
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[57]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less
# again, all accounted for. 

# Vs PacBio
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[55]\n";}' < illumina_megahit_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less
# Again, all accounted for.

# Let's see this in the pacbio dataset
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[56]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[57]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less
perl -ne '@F = split(/\t/); if($F[0] eq "name"){next;}elsif($F[50] ne "NOBIN"){print "$F[55]\n";}' < pacbio_final_pilon_master_table_2018_07_31.tab | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 | less

# Nothing! Looks like all of the DAStool bins are fully accounted for in already observed data. Let's try the hic bins instead
# Nothing again. Let's try a correlation style analysis instead
# First, let's check the pacbio unique contigs

summary(data.filt[!data.filt$IlluminaCtgAligns, c(1:11)])
     length            GC             cov17                cov1
 Min.   : 1047   Min.   :0.1628   Min.   :     0.02   Min.   :    0.47
 1st Qu.: 3285   1st Qu.:0.4142   1st Qu.:     1.94   1st Qu.:    9.27
 Median : 5776   Median :0.4899   Median :     3.03   Median :   19.45
 Mean   : 6109   Mean   :0.4802   Mean   :   423.42   Mean   :  153.79
 3rd Qu.: 8163   3rd Qu.:0.5538   3rd Qu.:     4.81   3rd Qu.:   42.73
 Max.   :31005   Max.   :0.9799   Max.   :289050.96   Max.   :43629.37

 superkingdom.t.24         phylum.t.28             order.t.32
 Archaea  :  9     Firmicutes    :389   Clostridiales   :313
 Bacteria :571     no-hit        :212   no-hit          :212
 Eukaryota:  4     Bacteroidetes : 77   Bacteroidales   : 64
 no-hit   :212     Actinobacteria: 29   Firmicutes-undef: 28
 Viruses  :  0     Proteobacteria: 26   Eggerthellales  : 20
                   Spirochaetes  : 14   Spirochaetales  : 14
                   (Other)       : 49   (Other)         :145
          family.t.36                  genus.t.40
 no-hit         :212   no-hit               :212
 Lachnospiraceae:120   Lachnospiraceae-undef: 50
 Ruminococcaceae: 63   Clostridium          : 48
 Clostridiaceae : 55   Ruminococcus         : 33
 Eubacteriaceae : 34   Eubacterium          : 30
 Prevotellaceae : 33   Firmicutes-undef     : 28
 (Other)        :279   (Other)              :395
                        species.t.44 TaxFuncCategory
 no-hit                       :212   Mode :logical
 Clostridiales bacterium      : 12   FALSE:543
 Sphingobacteriaceae bacterium:  7   TRUE :253
 Denitrobacterium detoxificans:  6   NA's :0
 Ruminococcus albus           :  6
 Ruminococcus sp. YE71        :  6
 (Other)                      :547

# So, most of these have read coverage and only 212 have no taxonomic assignment
summary(data.filt[!data.filt$IlluminaCtgAligns, c(12:17, 19,20)])
 MetabatBin        HiCBin          DASBin         CompleteORFs
 Mode :logical   Mode :logical   Mode :logical   Min.   : 0.000
 FALSE:711       FALSE:546       FALSE:770       1st Qu.: 4.000
 TRUE :85        TRUE :250       TRUE :26        Median : 9.000
 NA's :0         NA's :0         NA's :0         Mean   : 9.907
                                                 3rd Qu.:14.000
                                                 Max.   :54.000
  PartialORFs    ViralContig     MICKRMGAligns   Hungate1000Aligns
 Min.   :0.000   Mode :logical   Mode :logical   Mode :logical
 1st Qu.:1.000   FALSE:796       FALSE:682       FALSE:748
 Median :1.000   NA's :0         TRUE :114       TRUE :48
 Mean   :1.041                   NA's :0         NA's :0
 3rd Qu.:2.000
 Max.   :2.000

# Most (682 - 748) haven't been seen before and have complete ORFs! And are binned! Let's how the binned population turns out
sum(data.filt[!data.filt$IlluminaCtgAligns & data.filt$DASBin & !data.filt$MICKRMGAligns,c(1)])
[1] 137992  # <- 137,992 bases of sequence that have DASbins and not found in Mick's dataset

# Most at the higher GC spectrum end. High coverage in cov1 and cov17 libraries. All bacteria. Most Firmicutes. Some of these could be plasmid dna? No viruses

# Let's see the "no-hit" super-kingdom results
# 3755 contigs. Generally, lower GC hits. high cov17 and cov1 values. 59 DASBins and 1589 Hi-CBined contigs. Max 111 complete ORFs. 212 aren't present in the Illumina dataset. 8 are viral hosts. 349 found in Mick. 25 found in Hungate 1000
```

## Hacking together a Krona diagram

I want to do one more visualization for the team.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/master_tables

```bash
# reformatting taxonomic affiliation and printing counts to file
perl -ne '@F = split(/\t/); @b = ("root"); @t = (23,27,31,35,39,43); foreach my $k (@t){push(@b, $F[$k]);} print join(";", @b) . "\n";' < illumina_megahit_master_table_2018_07_31.tab | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 > illumina_megahit_master_table_tax.labels

perl -e '<>; while(<>){chomp; @s = split(/\t/); @b = split(/\;/, $s[0]); print "$s[1]\t" . join("\t", @b); print "\n";}' < illumina_megahit_master_table_tax.labels > illumina_megahit_master_table_tax.ktinput
/mnt/nfs/nfs2/bickhart-users/binaries/KronaTools-2.7/bin/bin/ktImportText -o illumina_megahit_master_table_tax.krona.html -n 'root' illumina_megahit_master_table_tax.ktinput

# That worked, let's try the pacbio now
perl -ne '@F = split(/\t/); @b = ("root"); @t = (23,27,31,35,39,43); foreach my $k (@t){push(@b, $F[$k]);} print join(";", @b) . "\n";' < pacbio_final_pilon_master_table_2018_07_31.tab |  perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 > pacbio_final_pilon_master_table_tax.labels

perl -e '<>; while(<>){chomp; @s = split(/\t/); @b = split(/\;/, $s[0]); print "$s[1]\t" . join("\t", @b); print "\n";}' < pacbio_final_pilon_master_table_tax.labels > pacbio_final_pilon_master_table_tax.ktinput
/mnt/nfs/nfs2/bickhart-users/binaries/KronaTools-2.7/bin/bin/ktImportText -o pacbio_final_pilon_master_table_tax.krona.html -n 'root' pacbio_final_pilon_master_table_tax.ktinput
```

## Generating annotation of prodigal identified ORFS

The goal of this exercise is to generate a listing of COGs, KEGGs and Go terms for both assemblies and to see what proportion of unique entries are present in the pacbio data.

> Ceres: ~/rumen_metagenome_assembly_project/binaries/eggnogmapper

```bash
# Downloading eggnogmapper databases
python download_eggnog_data.py euk bact arch viruses
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/eggnog

```bash
# Separating out the kingdom assignments for each ORF
module load diamond/0.9.18 hmmer3/3.1b2
sbatch --nodes=1 --ntasks-per-node=1 --mem=8000 -p short --wrap="perl segregate_kingdom_fastas.pl ../master_tables/illumina_megahit_master_table_2018_07_31.tab 23 ../prodigal/illumina_megahit_prodigal_proteins.faa illumina_megahit"
sbatch --nodes=1 --ntasks-per-node=1 --mem=8000 -p short --wrap="perl segregate_kingdom_fastas.pl ../master_tables/pacbio_final_pilon_master_table_2018_07_31.tab 23 ../prodigal/pacbio_final_prodigal_proteins.faa pacbio_final"

# Hmm, the version of diamond on the server is “too new” for the diamond databases used in the search
# I'm going to use the distributed Diamond binary in the eggnog-mapper package
export PATH=~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/bin/:$PATH

# Illumina
sbatch --nodes=1 --ntasks-per-node=15 --mem=120000 -p mem --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d bact -i illumina_megahit.Bacteria.split.faa --output illumina_megahit.Bacteria.eggnog --usemem --cpu 14"
sbatch --nodes=1 --ntasks-per-node=15 --mem=120000 -p mem --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d euk -i illumina_megahit.Eukaryota.split.faa --output illumina_megahit.Eukaryota.eggnog --usemem --cpu 14"
sbatch --nodes=1 --ntasks-per-node=15 --mem=100000 -p mem --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d arch -i illumina_megahit.Archaea.split.faa --output illumina_megahit.Archaea.eggnog --usemem --cpu 14"

# PacBio
sbatch --nodes=1 --ntasks-per-node=15 --mem=120000 -p mem --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d bact -i pacbio_final.Bacteria.split.faa --output pacbio_final.Bacteria.eggnog --usemem --cpu 14"
sbatch --nodes=1 --ntasks-per-node=15 --mem=120000 -p short --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d euk -i pacbio_final.Eukaryota.split.faa --output pacbio_final.Eukaryota.eggnog --usemem --cpu 14"
sbatch --nodes=1 --ntasks-per-node=15 --mem=120000 -p short --wrap="python ~/rumen_longread_metagenome_assembly/binaries/eggnog-mapper/emapper.py -m diamond -d arch -i pacbio_final.Archaea.split.faa --output pacbio_final.Archaea.eggnog --usemem --cpu 14"

# Now to generate the tables I need for further analysis
for i in Eukaryota Archaea Bacteria; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f pacbio_final.${i}.eggnog.emapper.annotations -c 11 -i '#' -d '\t' -m -o pacbio_final.${i}.full.cog.counts.md; done
for i in Eukaryota Archaea Bacteria; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f pacbio_final.${i}.eggnog.emapper.annotations -c 11 -i '#' -d '\t' -o pacbio_final.${i}.full.cog.counts.tab; done

# Using a complex pipe to separate out the pacbio unique cogs
# Eukaryotes
perl -ne '@F = split(/\t/); @bsegs = split(/_/, $F[0]); pop(@bsegs); $F[0] = join("_", @bsegs); print join("\t", @F);' < pacbio_final.Eukaryota.eggnog.emapper.annotations | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l ../master_tables/pacbio_final_no_other_db_aligns.ctgs.list -c 0 -i '#' -d '\t' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 11 -i '#' -d '\t'


# Bacteria
perl -ne '@F = split(/\t/); @bsegs = split(/_/, $F[0]); pop(@bsegs); $F[0] = join("_", @bsegs); print join("\t", @F);' < pacbio_final.Bacteria.eggnog.emapper.annotations | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l ../master_tables/pacbio_final_no_other_db_aligns.ctgs.list -c 0 -i '#' -d '\t' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 11 -i '#' -d '\t' > pacbio_final.Bacteria.noalgn.cog.counts.tab

# Archaea
perl -ne '@F = split(/\t/); @bsegs = split(/_/, $F[0]); pop(@bsegs); $F[0] = join("_", @bsegs); print join("\t", @F);' < pacbio_final.Archaea.eggnog.emapper.annotations | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l ../master_tables/pacbio_final_no_other_db_aligns.ctgs.list -c 0 -i '#' -d '\t' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 11 -i '#' -d '\t' > pacbio_final.Archaea.noalgn.cog.counts.tab
```

Let's generate a precursor COG proportion plot to show how the proportions change in the pacbio unique ORFs

> pwd: /home/dbickhart/share/metagenomics/pilot_manuscript/figure_drafts/eggnog

```bash
perl -e '<>; %h; $c = 0; while(<>){chomp; @s = split(/\t/); @bsegs = split(/, /, $s[0]); foreach my $k (@bsegs){$h{$k} += $s[1]; $c += $s[1];}} foreach my $k (sort{$a cmp $b}keys(%h)){print "$k\t" . ($h{$k} / $c) . "\tFull\n";}' < pacbio_final.Bacteria.full.cog.counts.tab > pacbio_final.Bacteria.full.cog.rfmt.tab
perl -e '<>; %h; $c = 0; while(<>){chomp; @s = split(/\t/); @bsegs = split(/, /, $s[0]); foreach my $k (@bsegs){$h{$k} += $s[1]; $c += $s[1];}} foreach my $k (sort{$a cmp $b}keys(%h)){print "$k\t" . ($h{$k} / $c) . "\tUnique\n";}' < pacbio_final.Bacteria.noalgn.cog.counts.tab > pacbio_final.Bacteria.noalgn.cog.rfmt.tab
```

```R
library(dplyr)
library(ggplot2)
library(RColorBrewer)
full <- read.delim("pacbio_final.Bacteria.full.cog.rfmt.tab", header=FALSE)
unique <- read.delim("pacbio_final.Bacteria.noalgn.cog.rfmt.tab", header=FALSE)
colnames(full) <- c("COG", "Proportion", "Type")
colnames(unique) <- c("COG", "Proportion", "Type")

combined <- bind_rows(full, unique)
combined <- mutate(combined, Perc = percent(Proportion))
combined <- group_by(combined, Type) %>% mutate(pos = cumsum(Proportion) - (0.5 * Proportion))

ggplot(combined, aes(y=Proportion, x = Type)) + geom_bar(aes(fill = COG), stat="identity") + geom_text(aes(label = Perc, y = pos), size = 3, vjust = 0.5, hjust = 0.5) + scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(24), guide = guide_legend(nrow=2)) + theme(legend.position="bottom")

dev.copy2pdf(file="pacbio_uniq_vs_full_cog_categories.pdf", useDingbats=FALSE)
```

## Misc figure and table generation

These are my notes for making miscellaneous figures and tables for the manuscript.

#### Figure 1 alt to blobplot

I want to make a nicer alternative to a blobplot that shows the types of contigs we are assembling, and compare them to the Illumina data. I may need to muck around with formats, but it may make for a nicer plot in the end.

> pwd: /home/dbickhart/share/metagenomics/pilot_manuscript/figure_drafts/master_tables

```bash
# First, reduce the size of the tables so that I can actually load the data into R!
perl -lane 'print "$F[1]\t$F[22]\t$F[23]\tIllumina";' < illumina_megahit_master_table_2018_07_31.tab > illumina_megahit_master_table_short_form_for_plot.tab
perl -lane 'print "$F[1]\t$F[22]\t$F[23]\tPacBio";' <  pacbio_final_pilon_master_table_2018_07_31.tab > pacbio_final_pilon_master_table_short_form_for_plot.tab
```
Now to try to plot them

```R
library(dplyr)
library(ggplot2)

illumina <- read.delim("illumina_megahit_master_table_short_form_for_plot.tab")
pacbio <- read.delim("pacbio_final_pilon_master_table_short_form_for_plot.tab")


colnames(pacbio) <- c("length", "covSum", "superKingdom", "Tech")
colnames(illumina) <- c("length", "covSum", "superKingdom", "Tech")
combined <- bind_rows(illumina, pacbio)

ggplot(combined, aes(x = superKingdom, y = length, size = covSum, fill=Tech)) + geom_point(shape=21, position = "jitter") + theme_bw() + scale_y_log10()
dev.copy2pdf(file="a_mess_scatterplot.pdf", useDingbats=FALSE)

library(ggridges)
combined.filt$superKingdom <- factor(combined.filt$superKingdom, levels = c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))
ggplot(combined.filt, aes(y=superKingdom, x=length, fill=Tech)) + geom_density_ridges() + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(100, 1000, 10000, 100000, 500000), limits=c(500, 500000), labels=c("100", "1000", "10,000", "100,000", "500,000")) + xlab(label = "Log10 Contig Length (bp)") + ylab(label = "Contig SuperKingdom Taxonomic Assignment")
dev.copy2pdf(file="contig_len_ridgeplots.pdf", useDingbats=FALSE)
```

#### Figure 2 Network alignment similarity plot

I want to generate a network plot of mash profile similarities between our datasets, Mick's, the Hungate1000 and refseq. This should give a good estimate of how similar they all are and how much novelty we contribute. I may want to kick off refseq at first, given the huge number of novel contigs in the dataset. I'll also remove all contigs sub 2kb from our illumina dataset to reduce the size of the hash.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/mash

```bash
# First pacbio, as that is the easiest to deal with
sbatch --nodes=1 --mem=30000 --ntasks-per-node=20 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -i -p 20 -o pacbio_final_pilon_indv_seq_k21_s1000 ../../assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa"

# Running all of Mick's genomes individually
mkdir mick_rug
for i in ../../assemblies/mick_rug/genomes/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch --nodes=1 --mem=2000 --ntasks-per-node=1 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -p 1 -o mick_rug/${name} $i"; done
~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash paste mick_rug_k21_s1000_combined ./mick_rug/*.msh

# Now to selectively filter the Illumina assembly
samtools faidx ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa
samtools faidx ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.gt2kb.fa
wc -l ../../assemblies/pilot_project/illumina_megahit/*.fai
  2182263 ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.fa.fai
   574382 ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.gt2kb.fa.fai
  2182263 ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa.fai
  4938908 total

# OK, > 2kb it is!
sbatch --nodes=1 --mem=30000 --ntasks-per-node=20 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash sketch -i -p 20 -o illumina_megahit_gt2kb_indiv_seq_k21_s1000 ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.gt2kb.fa"

# OK, now to try the distance estimation
# Serge used a cutoff of D <= 0.05 and p <=10^-10
# Looks like I may need to estimate the distance cutoff post-facto
sbatch --nodes=1 --mem=20000 --ntasks-per-node=10 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -p 10 -v 0.00001 -t illumina_megahit_gt2kb_indiv_seq_k21_s1000.msh pacbio_final_pilon_indv_seq_k21_s1000.msh hungate_1000.msh mick_rug_k21_s1000_combined.msh > rumen_only_mash_distance_table.tab"
```

#### Figure 3 bin quality assessment margin plots

I hope that this doesn't get messy too quickly, but I'm hoping to show differences between either the two different assemblies or the different binning methods in terms of overall assembly quality.

> pwd: /home/dbickhart/share/metagenomics/pilot_manuscript/figure_drafts/das_tool

```R
library(dplyr)
library(ggplot2)

# blank plot for spacing
blankPlot <- ggplot()+geom_blank(aes(1,1))+
   theme(
     plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank())

library(gridExtra)
illumina.hic <- read.delim("illumina_megahit_dastool_HiC.eval")
illumina.metabat <- read.delim("illumina_megahit_dastool_metabat.eval")
pacbio.hic <- read.delim("pacbio_final_dastool_HiC.eval")
pacbio.metabat <- read.delim("pacbio_final_dastool_metabat.eval")

pacbio.hic <- mutate(pacbio.hic, Tech = c("PacBio"), Bin = c("HiC"))
illumina.hic <- mutate(illumina.hic, Tech = c("Illumina"), Bin = c("HiC"))
illumina.metabat <- mutate(illumina.metabat, Tech = c("Illumina"), Bin = c("MetaBat"))
pacbio.metabat <- mutate(pacbio.metabat, Tech = c("PacBio"), Bin = c("MetaBat"))

keep <- c(1,12,13,14,15)
combined <- bind_rows(illumina.hic[,keep], illumina.metabat[,keep], pacbio.hic[,keep], pacbio.metabat[,keep])
combined$Tech <- as.factor(combined$Tech)
combined$Bin <- as.factor(combined$Bin)

# Now the individual plots
scatterPlot <- ggplot(combined, aes(x=SCG_completeness, y=SCG_redundancy, color=Tech)) + geom_point() + scale_color_brewer(palette="Dark2") + theme(legend.position=c(0,1), legend.justification=c(0,1))

xdensity <- ggplot(combined, aes(x=SCG_completeness, fill=Tech)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Dark2") + theme(legend.position = "none")

ydensity <- ggplot(combined, aes(x=SCG_redundancy, fill=Tech)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Dark2") + theme(legend.position = "none") + coord_flip()

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4,1.4), heights = c(1.4,4))
dev.copy2pdf(file="tech_density_scg_bins_plot.pdf", useDingbats=FALSE)

# That worked out well, but the theme was still present and the X axis text was still present. I can try to change that later.
# Now let's try a binning strategy plot instead
scatterPlot <- ggplot(combined, aes(x=SCG_completeness, y=SCG_redundancy, color=Bin)) + geom_point() + scale_color_brewer(palette="Paired") + theme(legend.position=c(0,1), legend.justification=c(0,1))

xdensity <- ggplot(combined, aes(x=SCG_completeness, fill=Bin)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Paired") + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank())

ydensity <- ggplot(combined, aes(x=SCG_redundancy, fill=Bin)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Paired") + theme(legend.position = "none", axis.title.y=element_blank(), axis.text.y=element_blank()) + coord_flip()

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4,1.4), heights = c(1.4,4))
dev.copy2pdf(file="binner_density_scg_bins_plot.pdf", useDingbats=FALSE)

# Final attempt, generate different shapes for each technology
scatterPlot <- ggplot(combined, aes(x=SCG_completeness, y=SCG_redundancy, color=Bin)) + geom_point(aes(shape=Tech)) + scale_color_brewer(palette="Paired") + theme(legend.position=c(0,1), legend.justification=c(0,1))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4,1.4), heights = c(1.4,4))
dev.copy2pdf(file="binner_density_scg_bins_plot_plustech.pdf", useDingbats=FALSE)
```

#### Redoing the CRISPR counts figure

I will need to do this again after I recalculate the CRISPR arrays on the frozen assemblies, but let's make a good placeholder plot in the meantime.

> pwd: ~/share/metagenomics/pilot_manuscript/figure_drafts/raw_stats

```R
# Using previously saved workspace
ggplot(total, aes(x=Dataset, y=Count, color=Dataset)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.75) + theme_bw() + scale_fill_brewer(palette="Dark2") + scale_colour_brewer(palette="Dark2")

dev.copy2pdf(file="crispr_violin_jitter_plot.pdf", useDingbats=FALSE)
```

#### Network analysis 

The mash analysis is far too data intensive and is unlikely to provide us with a useful plot. I want to plot the interactions between our two datasets with each other (Illumina + Pacbio) and with the previously posted Hungate and Mick RUG datasets. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/network

```bash
# Generating large-scale binary association files
perl -lane 'print "$F[0] pp 913_RUG";' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_vs_pb_pilon.assoc.condensed.tab > pb_pilon_mick_align.sif

perl -lane 'print "$F[0] pp 913_RUG";' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_vs_ilmn_megahit.assoc.condensed.tab > ilmn_megahit_mick_align.sif

perl -lane 'print "$F[0] pp HUNGATE";' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_vs_pb_pilon.assoc.condensed.tab > pb_pilon_hungate_align.sif

perl -lane 'print "$F[0] pp HUNGATE";' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_vs_ilmn_megahit.assoc.condensed.tab > ilmn_megahit_hungate_align.sif

# Now the pairwise comparisons. Do I want to do this only once? Let's test first
perl -lane '@bsegs = split(/;/, $F[1]); foreach $k (@bsegs){print $k;}' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_vs_pacbio_minimap2.assoc.condensed.tab > pb_contigs_ilmn_ref.list
perl -lane '@bsegs = split(/;/, $F[1]); foreach $k (@bsegs){print $k;}' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_vs_illumina_minimap2.assoc.condensed.tab > ilmn_contigs_pb_ref.list
perl -lane 'print $F[0];' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_vs_pacbio_minimap2.assoc.condensed.tab > ilmn_contigs_ilmn_ref.list
perl -lane 'print $F[0];' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_vs_illumina_minimap2.assoc.condensed.tab > pb_contigs_pb_ref.list

perl -lane 'print $F[0];' < ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai > pb_contigs_total.list
perl -lane 'print $F[0];' < ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa.fai > ilmn_contigs_total.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl pb_contigs_total.list pb_contigs_pb_ref.list pb_contigs_ilmn_ref.list
File Number 1: pb_contigs_total.list
File Number 2: pb_contigs_pb_ref.list
File Number 3: pb_contigs_ilmn_ref.list
Set     Count
1       777
1;2     2851
1;2;3   74023
1;3     19

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl ilmn_contigs_total.list ilmn_contigs_ilmn_ref.list ilmn_contigs_pb_ref.list
File Number 1: ilmn_contigs_total.list
File Number 2: ilmn_contigs_ilmn_ref.list
File Number 3: ilmn_contigs_pb_ref.list
Set     Count
1       1461320
1;2     286400
1;2;3   421764
1;3     12779

# OK, so I need to include both relationships just to make sure. Also, I need to capture the unique contigs to each dataset (within reason! The Illumina junk under 2.5 kb needs to go)
cat ilmn_megahit_hungate_align.sif ilmn_megahit_mick_align.sif | perl -lane 'print $F[0];' | sort | uniq > ilmn_ctgs_mick_hungate_align.list
cat pb_pilon_hungate_align.sif pb_pilon_mick_align.sif | perl -lane 'print $F[0];' | sort | uniq > pb_ctgs_mick_hungate_align.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o ilmn_contigs_total.list ilmn_contigs_ilmn_ref.list ilmn_contigs_pb_ref.list ilmn_ctgs_mick_hungate_align.list
File Number 1: ilmn_contigs_total.list
File Number 2: ilmn_contigs_ilmn_ref.list
File Number 3: ilmn_contigs_pb_ref.list
File Number 4: ilmn_ctgs_mick_hungate_align.list
Set     Count
1       1231048
1;2     119436
1;2;3   184248
1;2;3;4 237516
1;2;4   166964
1;3     5147
1;3;4   7632
1;4     230272

mv group_1.txt ilmn_ctgs_noaligns.list
# Let's remove all contigs that aren't larger than 2.5 kb to reduce the node list on cytoscape
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa.fai -c 0 -l ilmn_ctgs_noaligns.list | perl -lane 'if($F[1] < 2500){next;}else{print $F[0];}' > ilmn_ctgs_noaligns_gt2500.sif
wc -l ilmn_ctgs_noaligns_gt2500.sif
	201853 ilmn_ctgs_noaligns_gt2500.sif

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o pb_contigs_total.list pb_contigs_pb_ref.list pb_contigs_ilmn_ref.list pb_ctgs_mick_hungate_align.list
File Number 1: pb_contigs_total.list
File Number 2: pb_contigs_pb_ref.list
File Number 3: pb_contigs_ilmn_ref.list
File Number 4: pb_ctgs_mick_hungate_align.list
Set     Count
1       642
1;2     669
1;2;3   20569
1;2;3;4 53454
1;2;4   2182
1;3     13
1;3;4   6
1;4     135

mv group_1.txt pb_ctgs_noaligns.list
# In the sake of fairness, we'll apply the same filter to the pb contigs
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai -l pb_ctgs_noaligns.list -c 0 | perl -lane 'if($F[1] < 2500){next;}else{print $F[0];}' > pb_ctgs_noaligns_gt2500.sif
wc -l pb_ctgs_noaligns_gt2500.sif
	539 pb_ctgs_noaligns_gt2500.sif

# finally, let's create the cross dataset sif files
perl -lane '@bsegs = split(/;/, $F[1]); print "$F[0] pp " . join(" ", @bsegs);' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_vs_pacbio_minimap2.assoc.condensed.tab > ilmn_ref_pb_aligns.sif
perl -lane '@bsegs = split(/;/, $F[1]); print "$F[0] pp " . join(" ", @bsegs);' < /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_vs_illumina_minimap2.assoc.condensed.tab > pb_ref_ilmn_aligns.sif

wc -l *.sif
  201853 ilmn_ctgs_noaligns_gt2500.sif
  226734 ilmn_megahit_hungate_align.sif
  558184 ilmn_megahit_mick_align.sif
  708164 ilmn_ref_pb_aligns.sif
     539 pb_ctgs_noaligns_gt2500.sif
   32281 pb_pilon_hungate_align.sif
   52614 pb_pilon_mick_align.sif
   76874 pb_ref_ilmn_aligns.sif
 1857243 total

# Now, making a master file
cat *.sif > master_comparison_network_8_9_2018.sif

# OK, my attempts at running a huge analysis on cytoscape crashed. Time to run it on Gephi instead!
perl -lane 'if(scalar(@F) == 1){print "$F[0];$F[0]";}else{push(@temp, $F[0]); for($x = 2; $x < scalar(@F);$x++){push(@temp, $F[$x]);} print join(";", @temp); @temp = ()}' < master_comparison_network_8_9_2018.sif > master_comparison_network_8_9_2018.gephi.csv

# OK, the graph is still far too big. Let's convert the contigs into bins instead
perl convert_ctgs_to_bins.pl -i master_comparison_network_8_9_2018.gephi.csv -p ilmn.metabat -b ../dastool/illumina_megahit_public_metabat.unsorted.bins -o master_comparison_network_8_9_2018.ilmn.csv
perl convert_ctgs_to_bins.pl -i master_comparison_network_8_9_2018.ilmn.csv -e 'tig' -p pb.hic -b ../dastool/pacbio_final_public_hic.unsorted.bins -o master_comparison_network_8_9_2018.ilmn.pb.csv

sort master_comparison_network_8_9_2018.ilmn.pb.csv | uniq > master_comparison_network_8_9_2018.ilmn.pb.uniq.csv

# Damn, It's still too big. Trying to remove all meaningless self-self dataset associations
perl -lane '@bsegs = split(/;/, $F[0]); if($F[0] =~ /RUG/ || $F[0] =~ /HUNGATE/){print $F[0];}elsif($bsegs[0] eq $bsegs[1]){print $F[0];}elsif(scalar(@bsegs) == 2){$fn = substr($bsegs[0], 0, 2); $sn = substr($bsegs[1], 0, 2); if($fn eq $sn){print $F[0];}}' < master_comparison_network_8_9_2018.ilmn.pb.uniq.csv > master_comparison_network_8_9_2018.ilmn.pb.uniq.supcrop.csv

# OK, so lessons learned -- stick to only the important nodes! 
# Also, I need to use a better heuristic (not a repeated frequency scan) to condense the relationships I want to display in the graph

# Let's start over, albeit in a bin-oriented fashion
perl -ne 'chomp; @F = split(/\t/); if($F[49] ne "NOBIN" && $F[0] ne "name"){$F[49] = "pb.hic.$F[49]"; $F[56] = ($F[56] ne "-")? "913_RMG" : "-"; $F[57] = ($F[57] ne "-")? "HUNGATE" : "-"; if($F[56] eq "-" && $F[57] eq "-"){$F[56] = $F[49]; $F[57] = "";} print "$F[49]\t$F[56]\t$F[57]\n";}' < ../master_tables/pacbio_final_pilon_master_table_2018_07_31.tab | sort | uniq > pb_hic_bin_outside_db_relationship.tab

perl -e '%keep = ("913_RMG" => 1, "HUNGATE" => 1); %h; while(<>){chomp; @s = split(/\t/); for($x = 1; $x < scalar(@s); $x++){if($s[$x] ne "-"){$h{$s[0]}->{$s[$x]} = 1;}}} foreach my $bin (keys(%h)){my @assoc = (keys(%{$h{$bin}})); print "$bin"; foreach my $i (@assoc){if(scalar(@assoc) > 1 && $i ne $bin){print "\t$i";}elsif(scalar(@assoc) == 1){print "\t$i";}} print "\n";}' < pb_hic_bin_outside_db_relationship.tab > pb_hic_bin_outside_db_relationship.uniq.tab
perl -lane 'print join(";", @F);' < pb_hic_bin_outside_db_relationship.uniq.tab > pb_hic_bin_outside_db_relationship.uniq.csv

# I can actually load these into cytoscape now!
perl -lane 'for($x = 1; $x < scalar(@F); $x++){print "$F[0]\tpp\t$F[$x]";}' < pb_hic_bin_outside_db_relationship.uniq.tab > pb_hic_bin_outside_db_relationship.uniq.cyto.tab

# Now to do this for the Illumina megahit dataset
perl -ne 'chomp; @F = split(/\t/); if($F[48] ne "NOBIN" && $F[0] ne "name"){$F[48] = "ilmn.metab.$F[48]"; $F[56] = ($F[56] ne "-")? "913_RMG" : "-"; $F[57] = ($F[57] ne "-")? "HUNGATE" : "-"; if($F[56] eq "-" && $F[57] eq "-"){$F[56] = $F[48]; $F[57] = "";} print "$F[48]\t$F[56]\t$F[57]\n";}' < ../master_tables/illumina_megahit_master_table_2018_07_31.tab | sort | uniq > ilmn_metab_bin_outside_db_relationship.tab

perl -e '%keep = ("913_RMG" => 1, "HUNGATE" => 1); %h; while(<>){chomp; @s = split(/\t/); for($x = 1; $x < scalar(@s); $x++){if($s[$x] ne "-"){$h{$s[0]}->{$s[$x]} = 1;}}} foreach my $bin (keys(%h)){my @assoc = (keys(%{$h{$bin}})); print "$bin"; foreach my $i (@assoc){if(scalar(@assoc) > 1 && $i ne $bin){print "\t$i";}elsif(scalar(@assoc) == 1){print "\t$i";}} print "\n";}' < ilmn_metab_bin_outside_db_relationship.tab > ilmn_metab_bin_outside_db_relationship.uniq.tab

perl -lane 'for($x = 1; $x < scalar(@F); $x++){print "$F[0]\tpp\t$F[$x]";}' < ilmn_metab_bin_outside_db_relationship.uniq.tab > ilmn_metab_bin_outside_db_relationship.uniq.cyto.tab
```

#### Quick alignment of Itzik's rumen plasmidome dataset

I saw that Itzik had previously sequenced and assembled portions of the "plasmidome." I'd like to see how much of this is represented in our dataset as well. The MGRast accessions are [here](http://www.mg-rast.org/mgmain.html?mgpage=project&project=mgp505).

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mgrast_rumen_plasmidome

```bash
# Gathering general size information
samtools faidx 19988.fna
perl -lane 'print $F[1];' < 19988.fna.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   5771
Sum:    2710501
Minimum 100
Maximum 7779
Average 469.676139
Median  170
Standard Deviation      712.562593
Mode(Highest Distributed Value) 100


sbatch --nodes=1 --mem=15000 --ntasks-per-node=3 -p short --wrap="minimap2 -x asm5 ../pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa 19988.fna > itzik_plasmid_vs_pacbio.paf"

perl -lane 'if($F[11] > 0){print $_;}' < itzik_plasmid_vs_pacbio.paf | wc -l
239 <- alignments with MAPQ > 0
perl -lane 'if($F[11] > 0){print $F[5];}' < itzik_plasmid_vs_pacbio.paf | sort | uniq | wc -l
179 <- # unique pacbio contigs

perl -lane 'if($F[11] > 0){print $F[5];}' < itzik_plasmid_vs_pacbio.paf | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0
Entry   Value
tig00122470     10	<- 4kb contig 
tig00178536     9   <- 2.7 kb contig
tig00125079     6   <- 3.2 kb contig
tig00082331     5
tig00020805     5
tig00137488     4
tig00000072     3
tig00044700     3
tig00499823     3
tig00139131     2
tig00150954     2
tig00104376     2
tig00033950     2
tig02637158     2
tig00002986     2
tig00008812     2
tig00054211     2
tig00059277     2
tig00039471     2
tig00072063     2
...

perl -lane 'if($F[11] > 0){print $F[5];}' < itzik_plasmid_vs_pacbio.paf | sort | uniq > itzik_plasmid_vs_pacbio.pbctgs.gtzero.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../../analysis/master_tables/pacbio_final_pilon_master_table_2018_07_31.tab -l itzik_plasmid_vs_pacbio.pbctgs.gtzero.list -c 0 | perl -lane 'print $F[1];' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   179
Sum:    5552721
Minimum 1351
Maximum 414595
Average 31020.787709
Median  16273
Standard Deviation      48824.597125
Mode(Highest Distributed Value) 23649

sbatch --nodes=1 --mem=15000 --ntasks-per-node=3 -p short --wrap="minimap2 -x asm5 ../pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa 19988.fna > itzik_plasmid_vs_ilmn.paf"
perl -lane 'if($F[11] > 0){print $_;}' < itzik_plasmid_vs_ilmn.paf | wc -l
1676 <- alignments with MapQ > 0
perl -lane 'if($F[11] > 0){print $F[5];}' < itzik_plasmid_vs_ilmn.paf | sort | uniq | wc -l
832 <- unique illumina contigs

perl -lane 'if($F[11] > 0){print $F[5];}' < itzik_plasmid_vs_ilmn.paf | sort | uniq > itzik_plasmid_vs_ilmn.ilmctgs.gtzero.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../../analysis/master_tables/illumina_megahit_master_table_2018_07_31.tab -l  itzik_plasmid_vs_ilmn.ilmctgs.gtzero.list -c 0 | perl -lane 'print $F[1];' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   832
Sum:    3017981
Minimum 1002
Maximum 180507
Average 3627.381010
Median  2010
Standard Deviation      7700.645672
Mode(Highest Distributed Value) 1287
```

OK, I think that Itzik's plasmid sequence was contaminated or we've got some plasmid integration/chimerism in our assembly. What is interesting is that the pacbio assembly hits are fewer but tend to be larger contigs than the illumina. 

#### Generating contigs per bin plots

> pwd: share/metagenomics/pilot_manuscript/figure_drafts/das_tool/

```bash
library(dplyr)
library(ggplot2)

bins.ilmn.meta <- read.delim("illumina_megahit_dastool_metabat.eval")
bins.pac.meta <- read.delim("pacbio_final_dastool_metabat.eval")
bins.ilmn.hic <- read.delim("illumina_megahit_dastool_HiC.eval")
bins.pac.hic <- read.delim("pacbio_final_dastool_HiC.eval")

bins.ilmn.meta <- mutate(bins.ilmn.meta, Tech = "Illumina")
bins.ilmn.hic <- mutate(bins.ilmn.hic, Tech = "Illumina")
bins.pac.meta <- mutate(bins.pac.meta, Tech = "PacBio")
bins.pac.hic <- mutate(bins.pac.hic, Tech = "PacBio")
bins.ilmn.meta <- mutate(bins.ilmn.meta, Bin = "MetaBat")
bins.ilmn.hic <- mutate(bins.ilmn.hic, Bin = "HiC")
bins.pac.meta <- mutate(bins.pac.meta, Bin = "MetaBat")
bins.pac.hic <- mutate(bins.pac.hic, Bin = "HiC")
combined <- bind_rows(bins.ilmn.meta, bins.ilmn.hic, bins.pac.meta, bins.pac.hic)

combined <- mutate(combined, TechBin = paste0(Tech, Bin))
# Ridgeline plot didn't look as good here
ggplot(combined, aes(y=Bin, x=contigs, fill=Tech)) + geom_density_ridges() + theme_bw() + scale_fill_brewer(palette="Blues") + scale_x_log10()
library(gridExtra)

lower <- ggplot(combined, aes(x=TechBin, y=contigs, fill=Bin)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + theme_bw() + scale_fill_brewer(palette="Blues") + ylim(c(0,1000))
upper <- ggplot(combined, aes(x=TechBin, y=contigs, fill=Bin)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + theme_bw() + scale_fill_brewer(palette="Blues") + ylim(c(1000, 35000)) + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank())
grid.arrange(upper, lower, ncol=1)
dev.copy2pdf(file="bin_scores_by_contig_count.pdf", useDingbats=FALSE)
```

## Running DESMAN on the contigs

I want to address Adam's concern with strain-level identification using the longer-read contigs. Let's try to run portions of the DESMAN pipeline on the reads to see if we can pull the appropriate information from them.

First, making the DBs

> Assembler2:

```bash
wget https://desmandatabases.s3.climb.ac.uk/rpsblast_cog_db.tar.gz
ls rpsblast_cog_db
mv rpsblast_cog_db ./Databases/

wget https://desmandatabases.s3.climb.ac.uk/nr.faa
```

> Assembler2:/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon

```bash
module load samtools bwa
module load bedtools/2.27.1
export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin:$PATH

# I need the prodigal gff file first
sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p assemble2 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/Prodigal/prodigal -a usda_pacbio_second_pilon_indelsonly.prod.faa -c -d usda_pacbio_second_pilon_indelsonly.prod.fna -f gff -p meta -i usda_pacbio_second_pilon_indelsonly.fa -o usda_pacbio_second_pilon_indelsonly.prod.gff"

export COGSDB_DIR=/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/Databases/rpsblast_cog_db/
/mnt/nfs/nfs2/bickhart-users/binaries/CONCOCT/scripts/RPSBLAST.sh -f usda_pacbio_second_pilon_indelsonly.prod.faa -p -c 8 -r 1

# OK, I ran the first script, but it needs a cog to cddid file. Here's the comment line in the python script:
# Read a simple tsv file with two columns, cddid and cogid respectively.

```

#### attempts at running desman automated

This may be a disaster, but here goes!

```bash
source activate python3
module load samtools bwa java/jdk1.8.0_92 hmmer/3.1b2 
export PATH=$PATH:/mnt/nfs/nfs2/bickhart-users/binaries/bin
conda install cython numpy biopython

export DESMANHOME=/mnt/nfs/nfs2/bickhart-users/binaries/DESMAN/

bwa index usda_pacbio_second_pilon_indelsonly.fa
# Creating the contig lists
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); push(@{$h{$s[1]}}, $s[0]);} close IN; foreach my $k (keys(%h)){open(OUT, "> hic_bin_$k.list"); foreach my $c (@{$h{$k}}){print OUT "$c\n";} close OUT;}' pacbio_final_public_hic.unsorted.bins
# I rewrote the DESMAN nextflow script to avoid the lengthy alignment to the reference.
mv hic_*.list ./hic_bins/

# OK, I'm going to start with one of the high performing (ie. DAS-tool selected) bins
# NOTE: DAS-Tool fudges the numbers on the bins. You must search using contig assignment
grep 'tig00002711' ../pacbio_final_pilon/pacbio_final_dastool_DASTool_scaffolds2bin.txt
# That corresponded to hic_bin_4.list

/mnt/nfs/nfs2/bickhart-users/binaries/DESMAN/scripts/desmanflow_bams.nf --assembly=usda_pacbio_second_pilon_indelsonly.fa --speciescontigs=hic_bins/hic_bin_4.list --inputbams=./ --output=./hic_bin_4

# NOTE: The phylosift program is screwed up -- it tries to download files from a defunct database URL
# Worse: the desman script originaly overwrote my changes in the phylosiftrc file!
/mnt/nfs/nfs2/bickhart-users/binaries/DESMAN/scripts/desmanflow_bams.nf --assembly=usda_pacbio_second_pilon_indelsonly.fa --speciescontigs=hic_bins/hic_bin_4.list --inputbams=./ --output=./hic_bin_4 -resume --straincount=4

# It crashed again due to the weird pileup logic. I need to script this in another way.
# I gave up on the nextflow script. Let’s try this directly

python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/usda_pacbio_second_pilon_indelsonly.fa -c /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/hic_bins/hic_bin_345.list -b /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/USDA.sorted.merged.bam -o hic_bin_345 -d /mnt/nfs/nfs2/bickhart-users/binaries/DESMAN


# OK, my test worked. let's try processing all of the Hi-c bins with my strain inference script
for i in hic_bins/*.list; do name=`basename $i | cut -d'.' -f1`; echo $name; python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/usda_pacbio_second_pilon_indelsonly.fa -c /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/$i -b /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/USDA.sorted.merged.bam -o $name -d /mnt/nfs/nfs2/bickhart-users/binaries/DESMAN ; done

# Allot of the bins worked. If there wasn't enough information or the number of contigs per bin were small, the program failed to run.
# The PDF files are a good means of assessing program success though:
cp ./*/*.pdf ./strain_pdfs/

# I wrote a strain prediction script to automate the process. Let's see how it performs
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/usda_pacbio_second_pilon_indelsonly.fa -c /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/$i -b /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon/USDA.sorted.merged.bam -o hic_bin_4 -d /mnt/nfs/nfs2/bickhart-users/binaries/DESMAN

# ACK! It looks like the core genes I'm using to train the model are Ecoli centric! I've got the core of the model working, but I need to refine the gene selection and variant calling
```

I have an idea: let's use the DAS_tool output and the rest of the pipeline from DESMAN to calculate SCG frequency from there. I need to grep out ALL bins that have more than 0% completion. Just to flesh things out, I also want to add more samples to the analysis. Let's add all of Micks' samples and the high quality stuff (> 20 Gigbaytes total bam) to this list:

* PRJEB10338	cov0
* PRJEB21624	cov1
* PRJEB8939		cov2
* PRJNA214227	cov3
* PRJNA291523	cov7
* PRJNA60251	cov16
* USDA	<- ours	cov17

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics/pilot_project/pacbio_final_pilon

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i pacbio_final_public_hic.unsorted.bins,pacbio_final_public_metabat.unsorted.bins -c usda_pacbio_second_pilon_indelsonly.fa -o pacbio_final_dastool -l HiC,metabat --search_engine diamond -t 10 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db --write_bins 1 --proteins pacbio_final_prodigal_proteins.faa --score_threshold 0

# Converting them into the "analysis" bin set
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 10){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_dastool_analysis_binset_lt10redund.bins

# and the high quality bins for strain discovery
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 5 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_dastool_hq_binset_lt5redun_gt50comp.bins

wc -l *.bins
  55002 pacbio_dastool_analysis_binset_lt10redund.bins
   4929 pacbio_dastool_hq_binset_lt5redun_gt50comp.bins

cat pacbio_dastool_hq_binset_lt5redun_gt50comp.bins | cut -f2 | sort | uniq | wc -l
78  <- that's a healthy start

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 10){next;}else{ print "$s[7]\n";}}' < pacbio_final_dastool_DASTool_summary.txt | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
Sum:    798,163,918

# Now to get the SCG gene locations that DAS_tool used
cat pacbio_final_prodigal_proteins.faa.bacteria.scg pacbio_final_prodigal_proteins.faa.archaea.scg | perl -lane 'print $F[0];' > pacbio_final_prodigal_proteins.scg.cat.list
source activate python3

# Printing them to BED format
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_final_prodigal_proteins.shortform.tab -c 0 -l pacbio_final_prodigal_proteins.scg.cat.list | perl -lane '$r = $F[0]; $r =~ s/_\d{1,3}//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > pacbio_final_prodigal_proteins.scg.loc.bed

# Creating the coglists that desman needs
cat pacbio_final_dastool_proteins.faa.archaea.scg pacbio_final_dastool_proteins.faa.bacteria.scg > pacbio_final_dastool_proteins.faa.combined.scg
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$r = $s[0]; $r =~ s/_\d{1,3}//; print "$h{$s[0]},$r,$s[1],$s[2],$s[0],$s[3]\n";}} close IN;' pacbio_final_dastool_proteins.faa.combined.scg pacbio_final_prodigal_proteins.shortform.tab > pacbio_final_prodigal_master_cogs.csv
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/pacbio

```bash
export RUMEN=/scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly
module load gcc/8.1.0 prodigalorffinder/2.6.3 samtools/1.9 r/3.4.3

# I need to generate "list" files of all of the contigs belonging to das-tool bins that I will analyze
perl -lane 'open(OUT, ">> $F[1].hqdas.bin.list"); print OUT "$F[0]"; close OUT;' < ../../dastool/pacbio_dastool_hq_binset_lt5redun_gt50comp.bins

mkdir bin_lists
mv *.list ./bin_lists/

mkdir bed_lists
# Getting gene locations from the bin lists
for i in bin_lists/*.list; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ../../dastool/pacbio_final_prodigal_proteins.scg.loc.bed -c 0 -l $i > bed_lists/${name}.scg.bed; done

sbatch -p debug $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $RUMEN/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -c $RUMEN/analysis/desman/pacbio/bin_lists/pacbio_final_public_hic.110.hqdas.bin.list -g $RUMEN/analysis/desman/pacbio/bed_lists/pacbio_final_public_hic.110.scg.bed -d $RUMEN/binaries/DESMAN -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB10338/PRJEB10338.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB21624/PRJEB21624.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB8939/PRJEB8939.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA214227/PRJNA214227.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA291523/PRJNA291523.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA60251/PRJNA60251.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/USDA/USDA.sorted.merged.bam -o pacbio_hic_110

# My test worked out well enough. Let's queue the entire shebang
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch -p short $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $RUMEN/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -c $RUMEN/analysis/desman/pacbio/$i -g $RUMEN/analysis/desman/pacbio/bed_lists/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB10338/PRJEB10338.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB21624/PRJEB21624.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJEB8939/PRJEB8939.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA214227/PRJNA214227.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA291523/PRJNA291523.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/PRJNA60251/PRJNA60251.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/USDA/USDA.sorted.merged.bam -o $name; done

# Since the desman R plot script needs a ~/.Rlibs directory, it failed in each case
for i in */desman_dic.fits; do name=`echo $i | cut -d'/' -f1`; echo $name; /scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly/binaries/DESMAN/scripts/PlotDev.R -l $i -o ${name}.pdf; done

mkdir fit_pdfs
mv *.pdf ./fit_pdfs/

# Generating cog tables
for i in bin_lists/*.list; do num=`basename $i | cut -d '.' -f1,2`; echo $num; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} close IN; open(OUT, "> $ARGV[2]"); open(IN, "< $ARGV[1]"); while($l = <IN>){chomp($l); @s = split(/,/, $l); if(exists($h{$s[1]})){print OUT "$l\n";}} close OUT; close IN;' $i pacbio_final_prodigal_master_cogs.csv $num.scg.cogs.csv; done

mkdir cog_lists
mv *.cogs.csv ./cog_lists/

# Testing the strain prediction script
sbatch -p debug $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a $RUMEN/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -c $RUMEN/analysis/desman/pacbio/bin_lists/pacbio_final_public_hic.110.hqdas.bin.list -g $RUMEN/analysis/desman/pacbio/bed_lists/pacbio_final_public_hic.110.scg.bed -t $RUMEN/analysis/desman/pacbio/cog_lists/pacbio_final_public_hic.110.scg.cogs.csv -d $RUMEN/binaries/DESMAN -o pacbio_hic_110

# Queueing them all up
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch -p short $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a $RUMEN/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -c $RUMEN/analysis/desman/pacbio/$i -g $RUMEN/analysis/desman/pacbio/bed_lists/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -o $name -t $RUMEN/analysis/desman/pacbio/cog_lists/${bed}.scg.cogs.csv; done
```

I'm a little concerned that the strain selection results are a bit too monomorphic. I mean, the haplotypes are always around a 3 or 4 in the dataset. Let's try to rerun with only the USDA sequence data to see how that parses out.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/pacbio_sole

```bash
cp -r ../pacbio/bed_lists ./
cp -r ../pacbio/bin_lists ./
cp -r ../pacbio/cog_lists ./

for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch -p short $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $RUMEN/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa -c $RUMEN/analysis/desman/pacbio_sole/$i -g $RUMEN/analysis/desman/pacbio_sole/bed_lists/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -b $RUMEN/assemblies/pilot_project/pacbio_final_pilon/publicdb/USDA/USDA.sorted.merged.bam -o $name; done

# Queuing things up with my script:
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_pacbio_sole.sh $name $bed; done
```

I also want to see how bin strain assignment varies if I exclude the USDA dataset and only use the other datasets.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/pacbio_rest

```bash
# Just to keep things moving along, I've generated a queuing script for both sections of the pipeline
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_pacbio_rest.sh $name $bed; done
```


#### And the Illumina data

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin:$PATH; /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c mick_megahit_final_full.rfmt.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 10 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db --write_bins 1 --proteins illumina_megahit_prodigal_proteins.faa --score_threshold 0

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 10){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_megahit_dastool_DASTool_bins/$s[0].contigs.fa");while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < illumina_megahit_dastool_DASTool_summary.txt > illumina_dastool_analysis_binset_lt10redund.bins

cat illumina_dastool_analysis_binset_lt10redund.bins | cut -f2 | sort | uniq | wc -l
1630 <- more bins than the pacbio data!

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 5 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_megahit_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < illumina_megahit_dastool_DASTool_summary.txt > illumina_dastool_analysis_binset_lt5redun_gt50comp.bins

cat illumina_dastool_analysis_binset_lt5redun_gt50comp.bins | cut -f2 | sort | uniq | wc -l
119 <- OK, this is fair. More bins than pacbio but far less than expected given the coverage differential

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 10){next;}else{ print "$s[7]\n";}}' < illumina_megahit_dastool_DASTool_summary.txt | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
Sum:    1,475,288,770

# Generating the cog scg bins and ancillary files for desman
source activate python3
cat illumina_megahit_prodigal_proteins.faa.archaea.scg illumina_megahit_prodigal_proteins.faa.bacteria.scg > illumina_megahit_prodigal_proteins.faa.combined.scg
cat illumina_megahit_prodigal_proteins.faa.archaea.scg illumina_megahit_prodigal_proteins.faa.bacteria.scg | perl -lane 'print $F[0];' > illumina_megahit_prodigal_proteins.scg.cat.list
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f illumina_megahit_prodigal_proteins.shortform.tab -c 0 -l illumina_megahit_prodigal_proteins.scg.cat.list | perl -lane '$r = $F[0]; $r  =~ s/_\d{1,3}$//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > illumina_megahit_prodigal_proteins.scg.loc.bed
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$r = $s[0]; $r =~ s/_\d{1,3}$//; print "$h{$s[0]},$r,$s[1],$s[2],$s[0],$s[3]\n";}} close IN;' illumina_megahit_prodigal_proteins.faa.combined.scg illumina_megahit_prodigal_proteins.shortform.tab > illumina_megahit_prodigal_master_cogs.csv
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/illumina

```bash
export RUMEN=/scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly
module load gcc/8.1.0 prodigalorffinder/2.6.3 samtools/1.9 r/3.4.3

perl -lane 'open(OUT, ">> $F[1].hqdas.bin.list"); print OUT "$F[0]"; close OUT;' < ../../dastool/illumina_dastool_analysis_binset_lt5redun_gt50comp.bins

mkdir bin_lists
mv *.list ./bin_lists/

mkdir bed_lists
for i in bin_lists/*.list; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ../../dastool/illumina_megahit_prodigal_proteins.scg.loc.bed -c 0 -l $i > bed_lists/${name}.scg.bed; done

# Queuing up the major jobs
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch -p short $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainInference.py -a $RUMEN/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa -c $RUMEN/analysis/desman/illumina/$i -g $RUMEN/analysis/desman/illumina/bed_lists/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJEB10338/PRJEB10338.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJEB21624/PRJEB21624.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJEB8939/PRJEB8939.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJNA214227/PRJNA214227.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJNA291523/PRJNA291523.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/PRJNA60251/PRJNA60251.sorted.merged.bam -b $RUMEN/assemblies/pilot_project/illumina_megahit/publicdb/USDA/USDA.sorted.merged.bam -o $name; done

mkdir cog_lists
for i in bin_lists/*.list; do num=`basename $i | cut -d '.' -f1,2`; echo $num; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} close IN; open(OUT, "> $ARGV[2]"); open(IN, "< $ARGV[1]"); while($l = <IN>){chomp($l); @s = split(/,/, $l); if(exists($h{$s[1]})){print OUT "$l\n";}} close OUT; close IN;' $i illumina_megahit_prodigal_master_cogs.csv $num.scg.cogs.csv; done

mkdir fit_pdfs
cp ./*/*.pdf ./fit_pdfs/

for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch -p short $RUMEN/binaries/python_toolchain/metagenomics/desmanStrainPrediction.py -a $RUMEN/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa -c $RUMEN/analysis/desman/illumina/$i -g $RUMEN/analysis/desman/illumina/bed_lists/${bed}.scg.bed -d $RUMEN/binaries/DESMAN -t $RUMEN/analysis/desman/illumina/cog_lists/${bed}.scg.cogs.csv -o $name; done


#### And the limited sources of data ####
mkdir illumina_sole
mkdir illumina_rest
for i in bed_lists bin_lists cog_lists; do echo $i; cp -r illumina/$i illumina_sole; cp -r illumina/$i illumina_rest/; done

cd illumina_sole
export RUMEN=/scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly; module load gcc/8.1.0 prodigalorffinder/2.6.3 samtools/1.9 r/3.4.3; for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_illumina_sole.sh $name $bed; done

cd ../illumina_rest/
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_illumina_rest.sh $name $bed; done
```


And using Hansel and Gretel

```bash
module load samtools bcftools
bcftools mpileup -Ou --threads 5 -f usda_pacbio_second_pilon_indelsonly.fa USDA.sorted.merged.bam | bcftools call -vmO z -o USDA.sorted.pacbio.pilon.vcf.gz

bcftools index USDA.sorted.pacbio.pilon.vcf.gz
samtools index USDA.sorted.merged.bam

# I ran into problem with a single contig test. There was a divide by zero error

```

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
module load samtools bcftools
bcftools mpileup -Ou --threads 5 -f mick_megahit_final_full.rfmt.fa USDA.sorted.merged.bam | bcftools call -vmO z -o USDA.sorted.illumina.megahit.vcf.gz
```

## Finalizing bins

These are my notes on the selection of final DAS_tool bins from the dataset. 

Let's first check some distributions and see how the data stacks up.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_final_pilon/

```R
pacbio <- read.delim("pacbio_final_dastool_DASTool_summary.txt")

pacbio.filt <- pacbio[,c("binScore", "SCG_completeness", "SCG_redundancy")]
pacbio.filt$binScore <- pacbio.filt$binScore * 100

summary(pacbio.filt)
 SCG_completeness SCG_redundancy
 Min.   : 50.98   Min.   : 0.000
 1st Qu.: 62.75   1st Qu.: 0.000
 Median : 74.51   Median : 2.630
 Mean   : 74.58   Mean   : 5.231
 3rd Qu.: 87.25   3rd Qu.: 5.880
 Max.   :100.00   Max.   :37.250

library(ggplot2)
library(reshape)

pacbio.melt <- melt(pacbio.filt)
pdf(file="pacbio_dastool_select_bins.pdf", useDingbats=FALSE)
ggplot(pacbio.melt, aes(x=variable, y=value, fill=variable)) + geom_violin(trim=TRUE) + geom_boxplot(width=0.1, fill="white") + theme_minimal()

# That showed that only about 7 bins were outliers in the redundancy category
nrow(pacbio.filt[pacbio.filt$SCG_redundancy < 10,])
[1] 90

# Let's see that plot
pdf(file="pacbio_dastool_select_bins_lt10redundant.pdf", useDingbats=FALSE)
pacbio.melt <- melt(pacbio.filt[pacbio.filt$SCG_redundancy < 10,])
ggplot(pacbio.melt, aes(x=variable, y=value, fill=variable)) + geom_violin(trim=TRUE) + geom_boxplot(width=0.1, fill="white") + theme_minimal()
dev.off()

cor(pacbio.filt)
                   binScore SCG_completeness SCG_redundancy
binScore          1.0000000        0.7958932     -0.1421622
SCG_completeness  0.7958932        1.0000000      0.4738371
SCG_redundancy   -0.1421622        0.4738371      1.0000000

# OK, lets do this with the illumina data then.
illumina <- read.delim(file="../illumina_usda_accumulated/illumina_megahit_dastool_DASTool_summary.txt")
illumina.filt <- illumina[,c("binScore", "SCG_completeness", "SCG_redundancy")]
illumina.filt$binScore <- illumina.filt$binScore * 100

summary(illumina)
 SCG_completeness SCG_redundancy
 Min.   : 50.98   Min.   : 0.000
 1st Qu.: 66.67   1st Qu.: 1.960
 Median : 82.35   Median : 3.920
 Mean   : 79.09   Mean   : 4.969
 3rd Qu.: 90.20   3rd Qu.: 5.880
 Max.   :100.00   Max.   :23.530

illumina.melt <- melt(illumina.filt)
pdf(file="illumina_dastool_select_bins.pdf", useDingbats=FALSE)
ggplot(illumina.melt, aes(x=variable, y=value, fill=variable)) + geom_violin(trim=TRUE) + geom_boxplot(width=0.1, fill="white") + theme_minimal()
dev.off()

# Again, only a few that are greater than 10% SCG redundancy
nrow(illumina.filt[illumina.filt$SCG_redundancy < 10,])
[1] 179

cor(illumina.filt)
                   binScore SCG_completeness SCG_redundancy
binScore          1.0000000        0.9227543     -0.1079953
SCG_completeness  0.9227543        1.0000000      0.2664311
SCG_redundancy   -0.1079953        0.2664311      1.0000000

# Now, let's print a density plot of SCG completeness and redundancy across groups
illumina.melt <- melt(illumina.filt[illumina.filt$SCG_redundancy < 10,])

library(dplyr)
pacbio.melt <- mutate(pacbio.melt, Tech = c("PacBio"))
illumina.melt <- mutate(illumina.melt, Tech = c("Illumina"))
combined <- bind_rows(pacbio.melt, illumina.melt)

pdf(file="illumina_pacbio_bins_lt10redundancy.pdf", useDingbats=FALSE)
ggplot(combined, aes(x=variable, y=value, fill=variable)) + geom_violin(fill= "lightgrey") + geom_jitter(position = position_jitter(0.2)) + theme_bw() + facet_grid(Tech ~ .)
dev.off()

# Raw stats
sum(pacbio[pacbio$SCG_redundancy < 10,"contigs"])  =  6,496
sum(pacbio[pacbio$SCG_redundancy < 10,"size"])   =  153,755,142
mean(pacbio[pacbio$SCG_redundancy < 10,"N50"])   =  55,733.68

sum(illumina[illumina$SCG_redundancy < 10,"contigs"])  = 40,894
sum(illumina[illumina$SCG_redundancy < 10,"size"])   = 279,703,145
mean(illumina[illumina$SCG_redundancy < 10,"N50"])  =  13,458.51

# Printing bin names for easier association later:
write.table(pacbio[pacbio$SCG_redundancy < 10, c("bin", "size", "contigs", "N50", "binScore", "SCG_completeness", "SCG_redundancy")], file="pacbio_dastool_bins_lt10_redundancy.tab", quote=FALSE)
write.table(illumina[illumina$SCG_redundancy < 10, c("bin", "size", "contigs", "N50", "binScore", "SCG_completeness", "SCG_redundancy")], file="illumina_dastool_bins_lt10_redundancy.tab", quote=FALSE)
```

Now let's find out the statistics of the high quality bins vs these arbitrary bins vs the entire dataset for each method. Let's separate them into two categories

```bash
# Columns: 1. dastools bin, 2. bin length(bp), 3. # contigs, 4. SCG_completeness, 5. SCG_redundancy
perl -e '$outbase = "pacbio_dastool"; open(HQ, "> $outbase\_high_quality_dasbins.tab"); open(REST, "> $outbase\_analysis_dasbins.tab"); <STDIN>; while(<STDIN>){chomp; @s = split(/\s+/); if($s[6] > 80 && $s[7] < 5){print HQ "$s[1]\t$s[2]\t$s[3]\t$s[6]\t$s[7]\n";} print REST "$s[1]\t$s[2]\t$s[3]\t$s[6]\t$s[7]\n";}' < pacbio_dastool_bins_lt10_redundancy.tab
wc -l pacbio_dastool_analysis_dasbins.tab pacbio_dastool_high_quality_dasbins.tab
  90 pacbio_dastool_analysis_dasbins.tab
  22 pacbio_dastool_high_quality_dasbins.tab
 112 total

perl -e '$outbase = "illumina_dastool"; open(HQ, "> $outbase\_high_quality_dasbins.tab"); open(REST, "> $outbase\_analysis_dasbins.tab"); <STDIN>; while(<STDIN>){chomp; @s = split(/\s+/); if($s[6] > 80 && $s[7] < 5){print HQ "$s[1]\t$s[2]\t$s[3]\t$s[6]\t$s[7]\n";} print REST "$s[1]\t$s[2]\t$s[3]\t$s[6]\t$s[7]\n";}' < illumina_dastool_bins_lt10_redundancy.tab
wc -l illumina_dastool_high_quality_dasbins.tab illumina_dastool_analysis_dasbins.tab
   47 illumina_dastool_high_quality_dasbins.tab
  179 illumina_dastool_analysis_dasbins.tab
  226 total

source activate python3
perl -lane 'print $F[0];' < illumina_dastool_analysis_dasbins.tab > illumina_dastool_analysis_dasbins.list
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ../illumina_usda_accumulated/illumina_megahit_dastool_DASTool_scaffolds2bin.txt -c 1 -l illumina_dastool_analysis_dasbins.list > illumina_dastool_analysis_dasbins.contigs

perl -lane 'print $F[0];' < illumina_dastool_high_quality_dasbins.tab > illumina_dastool_high_quality_dasbins.list
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f ../illumina_usda_accumulated/illumina_megahit_dastool_DASTool_scaffolds2bin.txt -c 1 -l illumina_dastool_high_quality_dasbins.list > illumina_dastool_high_quality_dasbins.contigs

perl -lane 'print $F[0];' < pacbio_dastool_analysis_dasbins.tab > pacbio_dastool_analysis_dasbins.list
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_final_dastool_DASTool_scaffolds2bin.txt -c 1 -l pacbio_dastool_analysis_dasbins.list > pacbio_dastool_analysis_dasbins.contigs

perl -lane 'print $F[0];' < pacbio_dastool_high_quality_dasbins.tab > pacbio_dastool_high_quality_dasbins.list
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_final_dastool_DASTool_scaffolds2bin.txt -c 1 -l pacbio_dastool_high_quality_dasbins.list > pacbio_dastool_high_quality_dasbins.contigs
```

Alternative strategy: use DAS_tools SCG redundancy to filter the best sets of bins from each method.

```R
ilmetabat <- read.delim("../illumina_usda_accumulated/illumina_megahit_dastool_metabat.eval")
pbhic <- read.delim("pacbio_final_dastool_HiC.eval")

pbhic.filt <- mutate(pbhic[, c("contigs", "SCG_completeness", "SCG_redundancy")], Tech = "PacBio")
ilmetabat.filt <- mutate(ilmetabat[, c("contigs", "SCG_completeness", "SCG_redundancy")], Tech = "Illumina")
bin.filt <- bind_rows(pbhic.filt, ilmetabat.filt)

pdf(file="full_dastool_bin_SCG_stats.pdf", useDingbats=FALSE)

```