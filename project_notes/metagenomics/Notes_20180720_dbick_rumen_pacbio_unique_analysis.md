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
```


## TODO: Generate KAT plots for everything
## TODO: Also, combine the kmerspectra plots for Illumina and Pacbio
## TODO: Compute rand index for bins and check bin fidelity for pacbio-unique contigs with no illumina mappings
## TODO: look into viral genomes to estimate size, frequency and bias in pacbio and illumina datasets 