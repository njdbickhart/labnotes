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
	* DAS_tool concatenation
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/illumina_megahit_dastool_DASTool_bins.tab*
		* */home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool/pacbio_final_dastool_DASTool_bins.tab*
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

##### Alignment to other assembly

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

```