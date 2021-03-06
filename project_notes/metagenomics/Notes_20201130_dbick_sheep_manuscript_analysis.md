# Sheep manuscript analysis
---
*11/30/2020*

## Table of Contents

## Checklist and analysis

OK, let's organize this and prepare for the main manuscript push. I need to run some basic analysis and generate some useful statistics on all of the assemblies before we can prepare the first full draft.

Let's list them first:

### Necessary datasets

* Taxonomy 
	* List of classified taxa by bin -- GTDB-TK??
	* Comparisons of taxa by dataset -- Contigs vs length vs completeness in CLR vs HiFi
* Virus and plasmid analysis
	* Develop workflow for network analysis -- follow [tutorial](https://programminghistorian.org/en/lessons/exploring-and-analyzing-network-data-with-python#metrics-available-in-networkx)?
	* Get plasmid data from Itzik's lab
* Mechanisms of contig breaks
	* Align comparable taxa and identify sources of misassemblies and/or contig breaks
	* Compare against BGC data? 


## Taxonomic classification

I have some simple methods for taxonomic classification, but let's try a more complex method like GTDB-Tk first. What could go wrong?

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/gtdbtk

# make sure that GTDBTK_DATA_PATH=/KEEP/rumen_longread_metagenome_assembly/gtdbtk/share/gtdbtk-1.3.0/db/

mkdir gtdbtk_bins

# NOTE: It's infuriating, but I need to convert the '.' in the bin names to '_' for the program's purposes
module load samtools; for i in clr1 clr2 clr3 flye4; do echo $i; mkdir gtdbtk_bins/${i}.contigs; sbatch create_bins.pl b3c_${i}_dastool/${i}.das_bin3c.eval binning/bin3c/${i}.contigs/bin3c.full.clusters.tab ${i}.contigs.fasta gtdbtk_bins/${i}.contigs;  done

# ANOTHER NOTE: pplacer is apparently a nightmare memory hog. If it uses more than 1 thread, it explodes! 
mkdir gtdbtk_output
for i in clr1 clr2 clr3 flye4; do mkdir gtdbtk_output/${i}.contigs; sbatch -N 1 -n 30 --mem=700000 -p priority-mem -q msn-mem --wrap="gtdbtk classify_wf --genome_dir gtdbtk_bins/${i}.contigs --out_dir gtdbtk_output/${i}.contigs --cpus 30 --pplacer_cpus 1"; done

# OK let's concatenate the files and try to make it fit through Krona
cat gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/flye4.contigs/gtdbtk.ar122.summary.tsv | perl -ne 'chomp; @F = split(/\t/); if($F[0] =~ /^user_genome/){next;} print "$F[0]\t$F[1]\n";' > gtdbtk_output/flye4.contigs/flye4.concat_taxonomy.tab


perl -e 'chomp(@ARGV); %conv; %sum; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[1] =~ s/\./_/; $conv{$F[0]} = $F[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; if($_ =~ /^#/){next;} @F = split(/\t/); if(exists($conv{$F[0]})){$sum{$conv{$F[0]}} += $F[28];}} close IN; open(IN, "< $ARGV[2]"); while(<IN>){chomp; @F = split(/\t/); print "$_\t" . $sum{$F[0]} . "\n";} close IN;' b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt blobtools/table.flye4.contigs.blobDB.table.txt gtdbtk_output/flye4.contigs/flye4.concat_taxonomy.tab > gtdbtk_output/flye4.contigs/flye4.concat_taxonomy.pluscov.tab

# OK, now to package this up into something that Krona can understand
perl -ne 'chomp; @F = split(/\t/); $F[1] =~ s/\s+/_/g; @s = split(/;/, $F[1]); for($i = 0; $i < scalar(@s); $i++){$s[$i] =~ s/^.__//;} $F[2] = int($F[2]); print "$F[2]\t" . join("\t", @s) . "\n";' < gtdbtk_output/flye4.contigs/flye4.concat_taxonomy.pluscov.tab > gtdbtk_output/flye4.contigs/flye4.prekrona.tab

/lustre/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportText -o gtdbtk_output/flye4.contigs/flye4.krona.html gtdbtk_output/flye4.contigs/flye4.prekrona.tab


# That workflow made it happen! Let's run the CLR datasets through it and try to compare the output
for i in clr1 clr2 clr3; do echo $i; cat gtdbtk_output/${i}.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/${i}.contigs/gtdbtk.ar122.summary.tsv | perl -ne 'chomp; @F = split(/\t/); if($F[0] =~ /^user_genome/){next;} print "$F[0]\t$F[1]\n";' > gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.tab; perl -e 'chomp(@ARGV); %conv; %sum; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @F = split(/\t/); $F[1] =~ s/\./_/; $conv{$F[0]} = $F[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; if($_ =~ /^#/){next;} @F = split(/\t/); if(exists($conv{$F[0]})){$sum{$conv{$F[0]}} += $F[28];}} close IN; open(IN, "< $ARGV[2]"); while(<IN>){chomp; @F = split(/\t/); print "$_\t" . $sum{$F[0]} . "\n";} close IN;' b3c_${i}_dastool/${i}.das_DASTool_scaffolds2bin.txt blobtools/table.${i}.contigs.blobDB.table.txt gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.tab > gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.pluscov.tab; perl -ne 'chomp; @F = split(/\t/); $F[1] =~ s/\s+/_/g; @s = split(/;/, $F[1]); for($i = 0; $i < scalar(@s); $i++){$s[$i] =~ s/^.__//;} $F[2] = int($F[2]); print "$F[2]\t" . join("\t", @s) . "\n";' < gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.pluscov.tab > gtdbtk_output/${i}.contigs/${i}.prekrona.tab; /lustre/project/rumen_longread_metagenome_assembly/binaries/Krona/bin/ktImportText -o gtdbtk_output/${i}.contigs/${i}.krona.html gtdbtk_output/${i}.contigs/${i}.prekrona.tab; done

# Let's do a quick set analysis to see how things look from the raw taxonomy side
for i in flye4 clr1 clr2 clr3; do echo $i; perl -ne 'chomp; @F = split(/\t/); $F[1] =~ s/\s/_/g; print "$F[1]\n";' < gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.tab > gtdbtk_output/${i}.contigs/${i}.tax.list; done

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl gtdbtk_output/flye4.contigs/flye4.tax.list gtdbtk_output/clr1.contigs/clr1.tax.list gtdbtk_output/clr2.contigs/clr2.tax.list gtdbtk_output/clr3.contigs/clr3.tax.list
File Number 1: gtdbtk_output/flye4.contigs/flye4.tax.list
File Number 2: gtdbtk_output/clr1.contigs/clr1.tax.list
File Number 3: gtdbtk_output/clr2.contigs/clr2.tax.list
File Number 4: gtdbtk_output/clr3.contigs/clr3.tax.list
Set     Count
1       36
1;2     21
1;2;3   86
1;2;3;4 107
1;2;4   10
1;3     14
1;3;4   8
1;4     7
2       9
2;3     25
2;3;4   12
2;4     2
3       10
3;4     2
4       8
```

#### Graphlan trees

I'm going to do a quick Graphlan tree to show the diversity of the bacteria sequenced

> ubuntu: 

```bash
conda activate /mnt/c/SharedFolders/metagenomics/tim_sheep/graphlan

graphlan_annotate.py --annot ../initial_annotation.tab flye4/gtdbtk.bac120.classify.tree flye4.bac120.initial.xml

graphlan.py flye4.bac120.initial.xml  flye4.bac120.initial.png 
```

## Mechanisms of contig breaks

I think that this will be one of the big features. We can try to get the graph-guys to confirm, but let's make their jobs easier by identifying obvious regions that were broken within similar genome bins.

First, I'd like to do pair-wise bin alignments and then use minimap2 alignments to identify breaks in contigs (if they exist!). We'll use the gtdbk bin folders first.

In order to identify the most likely orthologous bin, I will use mash distance estimation on the bin3c bins. 

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
for i in clr1 clr2 clr3 flye4; do echo $i; sbatch mash_sketch_bins.pl gtdbtk_bins/${i}.contigs $i gtdbtk_bins/${i}_k21_s100000_combined; done

# Now, the pairwise comparisons between clr and hifi datasets
for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 2 --mem=35000 -p priority -q msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -v 0.1 gtdbtk_bins/flye4_k21_s100000_combined.msh gtdbtk_bins/${i}_k21_s100000_combined.msh > gtdbtk_bins/flye4_${i}_distance.tab"; done

# close to one-to-one alignments
for i in clr1 clr2 clr3; do echo $i; perl -lane 'if($F[2] < 0.1){print $_;}' < gtdbtk_bins/flye4_${i}_distance.tab | wc -l; done
clr1
421
clr2
403
clr3
186

# Bins that have > 1 alignments with distance < 0.1
for i in clr1 clr2 clr3; do echo $i; perl -e 'use File::Basename; %h; %d; while(<>){chomp; @s = split(/\t/); $r = basename($s[0]); $t = basename($s[1]); $r =~ s/\.fna//; $t =~ s/\.fna//; if($s[2] < 0.1){$h{$r} += 1; $d{$t} += 1}} $c = 0; foreach $k (keys(%h)){if($h{$k} > 1){$c++; print "$k,";}} print "\n"; $j = 0; foreach $k (keys(%d)){if($d{$k} > 1){$j++; print "$k,";}} print "\n$c\t$j\n";' < gtdbtk_bins/flye4_${i}_distance.tab; done
clr1
bin3c_12481,bin3c_2808,bin3c_237,bin3c_105,bin3c_180,bin3c_373,bin3c_186,bin3c_270,bin3c_459,
bin3c_335,bin3c_574,bin3c_497,bin3c_490,bin3c_35348,bin3c_617,bin3c_13,bin3c_216,bin3c_451,bin3c_258,bin3c_384,bin3c_303,bin3c_328,bin3c_424,bin3c_161,bin3c_35103,bin3c_457,bin3c_137,
9       18
clr2
bin3c_67,bin3c_304,bin3c_374,bin3c_166,bin3c_9757,bin3c_19505,bin3c_155,bin3c_110,bin3c_186,bin3c_236,
bin3c_122,bin3c_125,bin3c_4884,bin3c_171,bin3c_201,bin3c_327,bin3c_224,bin3c_14,bin3c_134,bin3c_89,bin3c_320,bin3c_12595,bin3c_581,bin3c_487,bin3c_35073,bin3c_74,bin3c_4880,
10      17
clr3
bin3c_118,bin3c_285,
bin3c_61,bin3c_59,bin3c_210,bin3c_167,bin3c_296,bin3c_13743,bin3c_139,bin3c_70,bin3c_236,bin3c_69,bin3c_346,bin3c_152,bin3c_60,bin3c_329,bin3c_142,
2       15

# So there are several possibilities here. A. the HiFi assembly delineated strains better (first number), B. the CLR assemblies had errors (second number) or C. one or the other assembly is wrong
# Let's make a graph and then try to associate taxonomy to each bin
for i in clr1 clr2 clr3; do echo $i; perl -e 'use File::Basename; while(<>){chomp; @s = split(/\t/); $r = basename($s[0]); $t = basename($s[1]); $r =~ s/\.fna//; $t =~ s/\.fna//; if($s[2] < 0.1){print "$r\t$t\t$s[2]\n";}}' < gtdbtk_bins/flye4_${i}_distance.tab > gtdbtk_bins/flye4_${i}_associations.tab; done

for i in clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); %bac = {}; %arch = {}; foreach $i (0,1){open(IN, $ARGV[$i]); while(<IN>){chomp; @s = split(/\t/);

for i in clr1 clr2 clr3; do echo $i; perl combine_bin_tax_and_pairing.pl gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/flye4.contigs/gtdbtk.ar122.summary.tsv gtdbtk_output/${i}.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/${i}.contigs/gtdbtk.ar122.summary.tsv gtdbtk_bins/flye4_${i}_associations.tab gtdbtk_output/flye4_${i}_assoc_tax.tab; done

# Let's do a quick count to see how many bins line up
for i in clr1 clr2 clr3; do echo $i; perl -ne 'chomp; @s = split(/\t/); if($s[-2] ne $s[-1]){print $_;}' < gtdbtk_output/flye4_${i}_assoc_tax.tab | wc -l ; done
clr1
0
clr2
0
clr3
0

# OK, let's queue up a series of "HiFi vs all" alignments based on the pairs of top hits
module load minimap2
mkdir gtdbtk_output/clr1_aligns; perl -lane 'system("sbatch minimap2_align.sh gtdbtk_bins/flye4.contigs/$F[0].fna gtdbtk_bins/clr1.contigs/$F[1].fna gtdbtk_output/clr1_aligns/$F[0].$F[1].algn.paf");' < gtdbtk_output/flye4_clr1_assoc_tax.tab

# Now printing out dotplots for quick viewing
conda activate /KEEP/rumen_longread_metagenome_assembly/r
mkdir gtdbtk_output/clr1_aligns/plots; for i in gtdbtk_output/clr1_aligns/*.paf; do name=`basename $i | cut -d'.' -f1,2`; echo $name;  sbatch -N 1 -n 2 --mem=9000 -p priority -q msn --wrap="Rscript /lustre/project/rumen_longread_metagenome_assembly/binaries/Themis-ASM/scripts/pafDotPlotly.R -i $i -o gtdbtk_output/clr1_aligns/plots/$name -q 5000 -m 5000 -l -s"; done


tar -czvf gtdbtk_output/clr1_aligns/clr1.plots.tar.gz gtdbtk_output/clr1_aligns/plots/

# I wonder what stats I can draw from my previous paf variant analysis
mkdir gtdbtk_output/clr1_aligns/vartabs; for i in gtdbtk_output/clr1_aligns/*.paf; do name=`basename $i | cut -d'.' -f1,2`; echo $name;  sbatch -N 1 -n 2 --mem=9000 -p priority -q msn --wrap="python3 /lustre/project/rumen_longread_metagenome_assembly/binaries/Themis-ASM/scripts/betweenAlignmentVariants.py -f $i -a 1000000 -q 1000000 -o gtdbtk_output/clr1_aligns/vartabs/$name.var.tab"; done

cat gtdbtk_output/clr1_aligns/vartabs/*.tab | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 6 -d '\t' -m
|Entry              | Value|
|:------------------|-----:|
|Repeat_contraction |  9963|
|Interchromosomal   |  4503|
|Repeat_expansion   |  3161|
|Inversion          |   998|
|Deletion           |   958|
|Tandem_contraction |   737|
|Tandem_expansion   |   491|
|type               |   421|
|Longrange          |   304|
|Insertion          |   247|

## TODO: Let's queue up the rest of the dataset and then plot distributions by size. I'm curious if there's a size bias that could be linked back to any particular gene families.
```

## Master table creation

I need to get most of these statistics under control so I will create major tables using my previous methodology. Starting with Blobtools tables and adding in GTDB-tk, and other tables.

```bash
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} print join("\t", @F[0..2, 8..11, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.flye4.contigs.blobDB.table.txt > flye4.mastertable.template.tab
```

##### Generating new information for the table

```bash
### Checkm
module load miniconda/3.6

# NOTE: If I don't remove the contigs under a certain size, they will clog up checkm and stop
module load samtools; for i in clr1 clr2 clr3; do echo $i; mkdir sub_ctgs_${i}; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %ctg; while(<IN>){chomp; @s = split(/\t/); if($s[1] > 60000){$ctg{$s[0]} = $s[1];}} $count = 0; @l; foreach $k (sort{$ctg{$b} <=> $ctg{$a}} keys(%ctg)){system("samtools faidx $ARGV[1] $k > $ARGV[2]/$k.fasta"); $count++; if($count > 1200){last;}}' $i.contigs.fasta.fai $i.contigs.fasta sub_ctgs_${i}; done

conda activate /KEEP/rumen_longread_metagenome_assembly/checkm
for i in clr1 clr2 clr3; do echo $i; sbatch --nodes=1 --mem=45000 --ntasks-per-node=8 -p priority -q msn --wrap="checkm lineage_wf -f $i.contigs.checkm.txt -t 8 -x fasta sub_ctgs_${i} checkm_contigs_${i}"; done


## bin-level associations
for i in flye4 clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[-2], $s[-1]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $d = join("\t", @{$data{$s[1]}}); print "$s[0]\t$d\n";} close IN;' b3c_${i}_dastool/${i}.das_DASTool_summary.txt b3c_${i}_dastool/${i}.das_DASTool_scaffolds2bin.txt > ${i}_binlevel_dastoolscores.tab; done

## GTDBK tax association
for i in clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $k = shift(@s); $data{$k} = [$s[0], $s[1], $s[2]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $d = join("\t", @{$data{$s[1]}}); print "$s[0]\t$d\n";} close IN;' gtdbtk_output/flye4_${i}_assoc_tax.tab b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt > flye4_${i}_gtdbk.tab; done

## And for flye4 GTDBK bin to contig conversion
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %data; while(<IN>){chomp; @s = split(/\t/); $k = shift(@s); $data{$k} = $s[0] . "";} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $d =  $data{$s[1]}; print "$s[0]\t$d\n";} close IN;' gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt > flye4_gtdbk_tax.tab

## Prodigal
module load prodigalorffinder/2.6.3
for i in clr1 clr2 clr2 flye4; do echo $i; sbatch -N 1 -n 2 --mem=45000 -p priority -q msn --wrap="prodigal -i $i.contigs.fasta -a $i.prodigal.proteins.faa -o $i.prodigal.stdout -s $i.prodigal.scores"; done

for i in flye4 clr1 clr2 clr3; do echo $i; grep '>' $i.prodigal.proteins.faa | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > $i.prodigal.shortform.tab; done

# Converting to short form for table
for i in *.prodigal.shortform.tab; do echo $i; perl -e '%tot; %comp; %part; while(<>){ chomp; @s = split(/\t/); if($_ =~ /^ContigID/){next;} @nsegs = split(/_/, $s[0]); $key = "$nsegs[0]\_$nsegs[1]"; if($s[5] eq "10" || $s[5] eq "01"){$part{$key} += 1;}else{$comp{$key} += 1;} $tot{$key} = 1;} foreach my $k (sort {$a cmp $b} keys(%tot)){$p = exists($part{$k})? $part{$k} : 0; $c = exists($comp{$k})? $comp{$k} : 0; print "$k\t$c\t$p\n";}' < $i > $i.short; done
```

OK, that's enough for now. Let's see how the table looks.

```bash
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j flye4_data_files_12_2020.json -o flye4_master_table_12_2020.tab -t flye4.mastertable.template.tab
```

##NOTE: todo use phylophlan analysis for species heatmap
##NOTE: todo use lefse to plot dendrogram of species occurrence and compare CLR propensity vs HiFI
##NOTE: todo use alignment depth comparison plots to check against low abundant strain assembly