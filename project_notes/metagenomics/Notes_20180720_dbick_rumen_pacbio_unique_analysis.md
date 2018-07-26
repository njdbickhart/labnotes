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


# Finally, generating kat sect plots for everything
sbatch --nodes=1 --mem=500000 --ntasks-per-node=40 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat sect -o illumina_run3_pbpilon_asm-sect -t 40 -E ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz"
sbatch --nodes=1 --mem=500000 --ntasks-per-node=40 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/KAT/src/kat sect -o illumina_run3_illumina_asm-sect -t 40 -E ../illumina_usda_accumulated/mick_megahit_final_full.rfmt.fa YMPrepCannula_run3_L2_R1.fastq.gz YMPrepCannula_run3_L2_R2.fastq.gz YMPrepCannula_run3_L3_R1.fastq.gz YMPrepCannula_run3_L3_R2.fastq.gz YMPrepCannula_run3_L4_R1.fastq.gz YMPrepCannula_run3_L4_R2.fastq.gz"
```

## TODO: Generate KAT plots for everything
## TODO: Also, combine the kmerspectra plots for Illumina and Pacbio
## TODO: Compute rand index for bins and check bin fidelity for pacbio-unique contigs with no illumina mappings
## TODO: look into viral genomes to estimate size, frequency and bias in pacbio and illumina datasets 

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
```