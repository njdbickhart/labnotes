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

# for R heatmap viewing later
for i in clr1 clr2 clr3; do echo $i; perl -MFile::Basename -lane 'if($F[2] < 0.1){$s = basename($F[0]); $r = basename($F[1]); $s =~ s/\..+$//; $r =~ s/\..+$//; print "$s\t$r\t$F[2]";}' < gtdbtk_bins/flye4_${i}_distance.tab > gtdbtk_bins/flye4_${i}_filtreformat_distance.tab; done 

###TODO: Download this data and run it through a heatmap program ####

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

## heatmap creation
for i in flye4 clr1 clr2 clr3; do echo $i; perl -lane '$F[1] =~ s/\./_/; print "$F[1].fna";' < b3c_${i}_dastool/$i.das_DASTool_scaffolds2bin.txt |sort | uniq > b3c_${i}_dastool/$i.scaffolds.list; done

# Had to do this for flye4:
cat binning/bin3c/flye4.contigs/*.tab | perl -lane '$F[1] =~ s/\./_/; print "$F[1].fna";' |sort | uniq > b3c_flye4_dastool/flye4.scaffolds.list
# And the rest apparently!
for i in clr1 clr2 clr3; do echo $i; cat binning/bin3c/$i.contigs/*.tab | perl -lane '$F[1] =~ s/\./_/; print "$F[1].fna";' |sort | uniq > b3c_${i}_dastool/$i.scaffolds.list; done

conda activate /KEEP/rumen_longread_metagenome_assembly/seaborn
###NOTE: I NEED to rewrite the code to eliminate combinations that don't happen in the file
sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python3 ~/python_toolchain/metagenomics/mashSparseToSortMatrix.py -r b3c_flye4_dastool/flye4.scaffolds.list -q b3c_clr1_dastool/clr1.scaffolds.list -d gtdbtk_bins/flye4_clr1_distance.tab -o gtdbtk_bins/flye4_clr1_distance.matrix"

for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python3 ~/python_toolchain/metagenomics/mashSparseToSortMatrix.py -r b3c_flye4_dastool/flye4.scaffolds.list -q b3c_${i}_dastool/${i}.scaffolds.list -d gtdbtk_bins/flye4_${i}_distance.tab -o gtdbtk_bins/flye4_${i}_distance.matrix"; done
```

## Master table creation

I need to get most of these statistics under control so I will create major tables using my previous methodology. Starting with Blobtools tables and adding in GTDB-tk, and other tables.

```bash
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 8..11, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.flye4.contigs.blobDB.table.txt > flye4.mastertable.template.tab

for i in clr1 clr2 clr3; do echo $i; perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 8..11, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.$i.contigs.blobDB.table.txt > $i.mastertable.template.tab; done

# NOTE: the coverage files are different! Don't do the above!
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 16..19, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.clr1.contigs.blobDB.table.txt > clr1.mastertable.template.tab
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 20..23, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.clr2.contigs.blobDB.table.txt > clr2.mastertable.template.tab
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 24..27, 29, 32, 35, 38, 41, 44]) . "\n";' < blobtools/table.clr3.contigs.blobDB.table.txt > clr3.mastertable.template.tab
```

##### Generating new information for the table

```bash
### Checkm
module load miniconda/3.6

# NOTE: If I don't remove the contigs under a certain size, they will clog up checkm and stop
module load samtools; for i in clr1 clr2 clr3; do echo $i; mkdir sub_ctgs_${i}; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %ctg; while(<IN>){chomp; @s = split(/\t/); if($s[1] > 90000){$ctg{$s[0]} = $s[1];}} $count = 0; @l; foreach $k (sort{$ctg{$b} <=> $ctg{$a}} keys(%ctg)){system("samtools faidx $ARGV[1] $k > $ARGV[2]/$k.fasta"); $count++; if($count > 1200){last;}}' $i.contigs.fasta.fai $i.contigs.fasta sub_ctgs_${i}; done

conda activate /KEEP/rumen_longread_metagenome_assembly/checkm
for i in clr1 clr2 clr3; do echo $i; sbatch --nodes=1 --mem=45000 --ntasks-per-node=8 -p priority -q msn --wrap="checkm lineage_wf -f $i.contigs.checkm.txt -t 8 -x fasta sub_ctgs_${i} checkm_contigs_${i}"; done

for i in clr1 clr2 clr3; do echo $i; perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); print "$s[0]\t$s[-3]\t$s[-2]\n";}' < $i.contigs.checkm.txt > ${i}_checkm_contigs.tab; done


## bin-level associations
for i in flye4 clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[-2], $s[-1]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $d = join("\t", @{$data{$s[1]}}); print "$s[0]\t$d\n";} close IN;' b3c_${i}_dastool/${i}.das_DASTool_summary.txt b3c_${i}_dastool/${i}.das_DASTool_scaffolds2bin.txt > ${i}_binlevel_dastoolscores.tab; done

## GTDBK tax association
for i in clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $k = shift(@s); for($x = 0; $x < scalar(@s); $x++){$s[$x] =~ s/\s+/_/g;} $data{$k} = [$s[0], $s[1], $s[2]];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $d = join("\t", @{$data{$s[1]}}); print "$s[0]\t$d\n";} close IN;' gtdbtk_output/flye4_${i}_assoc_tax.tab b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt > flye4_${i}_gtdbk.tab; done

## And for flye4 GTDBK bin to contig conversion
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %data; while(<IN>){chomp; @s = split(/\t/); $k = shift(@s); for($x = 0; $x < scalar(@s); $x++){$s[$x] =~ s/\s+/_/g;} $data{$k} = $s[0] . "";} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $s[1] =~ s/\./_/g; $d =  $data{$s[1]}; print "$s[0]\t$d\n";} close IN;' gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt > flye4_gtdbk_tax.tab

## Prodigal
module load prodigalorffinder/2.6.3
for i in clr1 clr2 clr2 flye4; do echo $i; sbatch -N 1 -n 2 --mem=45000 -p priority -q msn --wrap="prodigal -i $i.contigs.fasta -a $i.prodigal.proteins.faa -o $i.prodigal.stdout -s $i.prodigal.scores"; done

for i in flye4 clr1 clr2 clr3; do echo $i; grep '>' $i.prodigal.proteins.faa | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > $i.prodigal.shortform.tab; done

# Converting to short form for table
for i in *.prodigal.shortform.tab; do echo $i; perl -e '%tot; %comp; %part; while(<>){ chomp; @s = split(/\t/); if($_ =~ /^ContigID/){next;} @nsegs = split(/_/, $s[0]); $key = "$nsegs[0]\_$nsegs[1]"; if($s[5] eq "10" || $s[5] eq "01"){$part{$key} += 1;}else{$comp{$key} += 1;} $tot{$key} = 1;} foreach my $k (sort {$a cmp $b} keys(%tot)){$p = exists($part{$k})? $part{$k} : 0; $c = exists($comp{$k})? $comp{$k} : 0; print "$k\t$c\t$p\n";}' < $i > $i.short; done

# Viruses
for i in clr1 clr2 clr3; do echo $i; perl -lane 'if($F[0] =~ /^Virus/){next;}else{print "$F[1]\t$F[0]\t$F[3]";}' < $i.contigs.vassoc.final.tab > $i.contigs.vassoc.assoc.tab; done
```


OK, that's enough for now. Let's see how the table looks.

```bash
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j flye4_data_files_12_2020.json -o flye4_master_table_12_2020.tab -t flye4.mastertable.template.tab

for i in clr1 clr2 clr3; do echo $i; python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j ${i}_data_files_12_2020.json -o ${i}_master_table_12_2020.tab -t ${i}.mastertable.template.tab; done
```

## Strain characterization

I rewrote portions of the cDNA_cupcake API to allow me to run a separate MAG-level strain phaser. Let's see how that goes!

```bash
module load miniconda/3.6 samtools
conda activate /KEEP/rumen_longread_metagenome_assembly/desman

# Testing it out first
sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g desman/bed_lists/flye4/bin3c.522.scg.bed -o test.pbhap"

mkdir strain_phase
mkdir strain_phase/flye4
# OK, that works but let's try to queue it up for all bins
for i in desman/bed_lists/flye4/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g $i -o strain_phase/flye4/$name.strain"; done

# Let's try this on the other samples
for j in clr1 clr2 clr3; do echo $j; mkdir strain_phase/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/mag_phaser.py -a $j.contigs.fasta -b $j.contigs.fasta.ccs.bam -g $i -o strain_phase/$j/$name.strain"; done; done

for i in flye4 clr1 clr2 clr3; do echo $i; ls strain_phase/$i/*.human_readable.txt | wc -l ; done
flye4
169
clr1
192
clr2
204
clr3
176

for i in flye4 clr1 clr2 clr3; do echo $i; wc -l strain_phase/$i/*.human_readable.txt | perl -lane 'print "$F[0]";' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl ; done
flye4
total   170
Sum:    506282
Minimum 3
Maximum 253141
Average 2978.129412
Median  327
Standard Deviation      19474.550808
Mode(Highest Distributed Value) 3

clr1
total   193
Sum:    1349252
Minimum 1
Maximum 674626
Average 6990.943005
Median  583
Standard Deviation      49332.758647
Mode(Highest Distributed Value) 3

clr2
total   205
Sum:    1350248
Minimum 1
Maximum 675124
Average 6586.575610
Median  785
Standard Deviation      47450.792098
Mode(Highest Distributed Value) 3

clr3
total   177
Sum:    4542594
Minimum 1
Maximum 2271297
Average 25664.372881
Median  801
Standard Deviation      184853.213062
Mode(Highest Distributed Value) 5

# I think this is promising! Let me first test the results if I get a full bin analysis going

# Testing it out on a large bin 3.9 megabases.
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -c 0 -f flye4.contigs.fasta.fai  -l contig_15947,contig_1616,contig_1934,contig_1935,contig_22736,contig_23300,contig_26510,contig_27544,contig_27546,contig_27548,contig_27551,contig_38385,contig_38393,contig_39828,contig_4661,contig_49279,contig_54841,contig_54854,contig_58423,contig_59232,contig_5930,contig_63482,contig_65750,contig_65874,contig_74799,contig_7633 -d '\t' | perl -lane 'print "$F[0]\t1\t$F[1]";' > temp_bin3c_79_flye4.bed

sbatch -N 1 -n 2 -p priority -q msn --mem=55000 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g temp_bin3c_79_flye4.bed -o test.bin3c_79"

# Nope, I don't think that it's worth it. 

# Instead, let's plot out the differences
python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/plotPhaseCount.py -f b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.115.contigs.fa.fai -b flye4.contigs.fasta.ccs.bam -u strain_phase/flye4/bin3c.115.strain.human_readable.txt -o test.bin3c.115.flye4

# OK, queuing up the whole shebang
mkdir strain_phase/plots
for i in flye4 clr1 clr2 clr3; do echo $i; mkdir strain_phase/plots/$i; for j in strain_phase/$i/*.human_readable.txt; do name=`basename $j | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 --mem=55000 -p priority -q msn --wrap=" samtools faidx b3c_${i}_dastool/${i}.das_DASTool_bins/${name}.contigs.fa ; python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/plotPhaseCount.py -f b3c_${i}_dastool/${i}.das_DASTool_bins/${name}.contigs.fa.fai -b $i.contigs.fasta.ccs.bam -u $j -o strain_phase/plots/$i/$name.strains"; done; done

# Now, associating flye4 bins with their nearest neighbors in the clr datasets for strain counts
for i in flye4 clr1 clr2 clr3; do for j in `ls strain_phase/$i/* | xargs -I {} basename {} | cut -d'.' -f1,2 | sort | uniq`; do perl -e 'chomp(@ARGV); if( -s $ARGV[0].strain.NO_SNPS_FOUND ){print "$ARGV[0]\t0\n";}else{ open(IN, "< strain_phase/$ARGV[1]/$ARGV[0].strain.human_readable.txt"); %count; <IN>; while(<IN>){chomp; @s = split(/\t/); $count{$s[4]} += 1;} close IN; print "$ARGV[0]\t" . scalar(keys(%count)) . "\n";}' $j $i; done > strain_phase/$i.bin_straincount.tab; echo $i; done

for i in strain_phase/*_straincount.tab; do echo $i; cat $i | cut -f2 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
strain_phase/clr1.bin_straincount.tab
total   557
Sum:    414
Minimum 0
Maximum 4
Average 0.743268
Median  0
Standard Deviation      1.064711
Mode(Highest Distributed Value) 0

strain_phase/clr2.bin_straincount.tab
total   552
Sum:    448
Minimum 0
Maximum 4
Average 0.811594
Median  0
Standard Deviation      1.094787
Mode(Highest Distributed Value) 0

strain_phase/clr3.bin_straincount.tab
total   279
Sum:    409
Minimum 0
Maximum 4
Average 1.465950
Median  2
Standard Deviation      1.216533
Mode(Highest Distributed Value) 2

strain_phase/flye4.bin_straincount.tab
total   695
Sum:    352
Minimum 0
Maximum 3
Average 0.506475
Median  0
Standard Deviation      0.904457
Mode(Highest Distributed Value) 0

# Adding completeness and contamination estimates for bins
for i in flye4 clr1 clr2 clr3; do echo $i; perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if(!exists($h{$s[0]})){print STDERR "Error in $s[0]\n!";} print "$s[0]\t$h{$s[0]}\t$s[11]\t$s[12]\n"} close IN;' strain_phase/$i.bin_straincount.tab b3c_${i}_dastool/$i.das_DASTool_summary.txt > strain_phase/$i.bin_pluscomp.tab; done

python3 ~/python_toolchain/metagenomics/createMAGstrainAssoc.py -r strain_phase/flye4.bin_pluscomp.tab -o strain_phase/flye4_strain_associations -n clr1 -a gtdbtk_bins/flye4_clr1_associations.tab -s strain_phase/clr1.bin_pluscomp.tab -n clr2 -a gtdbtk_bins/flye4_clr2_associations.tab -s strain_phase/clr2.bin_pluscomp.tab -n clr3 -a gtdbtk_bins/flye4_clr3_associations.tab -s strain_phase/clr3.bin_pluscomp.tab
```

Now to try to plot it in R so we can highlight strain differences.

```R
library(dplyr)
library(ggplot2)
library(ggExtra)

setwd("C:/SharedFolders/metagenomics/tim_sheep/strain_determination/")


data <- read.delim("flye4_strain_associations.melt", header=TRUE, sep="\t")
data <- data[!data$AssocComp == -1,]

data$de <- "Nominal"
data$de[data$StrainDelta >= 2 | data$StrainDelta <= -2] <- "Different"
data$de <- as.factor(data$de)

p <- ggplot(data=data, aes(y=AssocCont, x=as.factor(StrainDelta), col=de)) + geom_boxplot() + geom_jitter() + scale_fill_brewer(palette="Dark2") + theme_bw()
p

ggplot(data=data, aes(x=StrainDelta, fill="Count")) +geom_density() + theme_bw() + xlab(label="Strain Delta (<- more in HIFI; -> more in CLR)")

## TODO: find a better way to visualize this dataset
```

##NOTE: todo use phylophlan analysis for species heatmap
##NOTE: todo use lefse to plot dendrogram of species occurrence and compare CLR propensity vs HiFI
##NOTE: todo use alignment depth comparison plots to check against low abundant strain assembly

## Coverage depth vs completeness and bin completeness plots

```bash
for i in clr1 clr2 clr2; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[19] < 5.0 && $s[18] ne "-"){$sum = $s[3] + $s[4] + $s[5] + $s[6]; print "$s[0]\t$sum\t$s[18]\t$ARGV[1]\n";}} close IN;' ${i}_master_table_12_2020.tab $i > ${i}_ctg_cov_by_comp.tab; done

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[18] ne "-"){$sum = $s[3] + $s[4] + $s[5] + $s[6]; print "$s[0]\t$sum\t$s[18]\t$ARGV[1]\n";}} close IN;' flye4_master_table_12_2020.tab hifi > flye4_ctg_cov_by_comp.tab

```
 Now to plot it in R

```R
library(dplyr)
library(ggplot2)
library(ggridges)

flye4 <- read.delim("flye4_ctg_cov_by_comp.tab", header=FALSE)
clr1 <- read.delim("clr1_ctg_cov_by_comp.tab", header=FALSE)
clr2 <- read.delim("clr2_ctg_cov_by_comp.tab", header=FALSE)
clr3 <- read.delim("clr3_ctg_cov_by_comp.tab", header=FALSE)

colnames(flye4) <- c("Contig", "Coverage", "Completeness", "Assembly")
colnames(clr1) <- c("Contig", "Coverage", "Completeness", "Assembly")
colnames(clr2) <- c("Contig", "Coverage", "Completeness", "Assembly")
colnames(clr3) <- c("Contig", "Coverage", "Completeness", "Assembly")

combined <- bind_rows(flye4, clr1, clr2, clr3)

combined$Assembly <- as.factor(combined$Assembly)
ggplot(data = combined, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10()

ggplot(data = combined, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.copy2pdf(file="facet_cov_by_comp_scatter.pdf", useDingbats=FALSE)

# Let's view this for just the top 80% complete contigs
top80 <- combined[combined$Completeness >= 80,]
summary(top80)
    Contig             Coverage         Completeness    Assembly  
 Length:333         Min.   :   2.614   Min.   : 80.01   clr1: 68  
 Class :character   1st Qu.:  13.101   1st Qu.: 91.55   clr2: 73  
 Mode  :character   Median :  31.237   Median : 94.50   clr3: 37  
                    Mean   : 170.323   Mean   : 93.58   hifi:155  
                    3rd Qu.:  95.643   3rd Qu.: 97.37             
                    Max.   :3370.070   Max.   :100.00

ggplot(data = top80, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.copy2pdf(file="facet_cov_by_comp_top80.pdf", useDingbats=FALSE)

# OK, clear improvement. Let's see though if we can show this in a ridgeline distribution
ggplot(data=top80[top80$Completeness >=95,], aes(y=Assembly, x=Coverage, fill=Assembly)) + geom_density_ridges() + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
dev.copy2pdf(file="ridge_cov_at_95_log.pdf", useDingbats=FALSE)

# SO, looks like there are bimodal (trimodal?) distribution of contig coverage in the CLR datasets but not the HIFI
# Strain assembly related? Or artifacts?

# Let's try to prepare this for the figure by combining CLR datasets and setting colors
top80 <- top80 %>% mutate(ASM = ifelse(Assembly == "hifi", "HIFI", "CLR"))
ggplot(data=top80[top80$Completeness >=95,], aes(y=Assembly, x=Coverage, fill=ASM)) + geom_density_ridges() + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
dev.copy2pdf(file="ridge_covcol_95_log.pdf", useDingbats=FALSE)

ggplot(data = top80, aes(x=Completeness, y=Coverage, colour=ASM)) + geom_point() + scale_colour_brewer(palette='Dark2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.copy2pdf(file="facet_cov_by_comp_top80.pdf", useDingbats=FALSE)

```


## Viral statistics

```bash
for i in flye4 clr1 clr2 clr3; do echo $i; perl -ne '@F=split(/\t/); if($F[7] == "Viruses" && $F[1] > 1000){print "$F[1]\n";}' < ${i}_master_table_12_2020.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl ; done

flye4
total   57256
Sum:    3424465269
Minimum 1001
Maximum 5546585
Average 59809.718964
Median  15287
Standard Deviation      192075.704295
Mode(Highest Distributed Value) 1489

clr1
total   48334
Sum:    2984842227
Minimum 1001
Maximum 7141023
Average 61754.504634
Median  23646.5
Standard Deviation      158810.782120
Mode(Highest Distributed Value) 1831

clr2
total   48790
Sum:    3008459117
Minimum 1001
Maximum 6067615
Average 61661.387928
Median  23121.5
Standard Deviation      159437.994227
Mode(Highest Distributed Value) 3449

clr3
total   24315
Sum:    1675674347
Minimum 1001
Maximum 5583265
Average 68915.251779
Median  34266
Standard Deviation      163346.658411
Mode(Highest Distributed Value) 3038

```

### Strain file preparation

I need to prepare files for Liz to test. I will select from the following MAGs so that she can see differences between the assemblies as well.

* Different HIFI - bin3c_44 and bin3c_42 
	* CLR1 bin3c_13
	* CLR2 bin3c_14
* Normal HIFI - bin3c_146
	* CLR1 bin3c_171
	* CLR2 bin3c_154
* Normal HIFI - bin3c_120
	* CLR1 bin3c_94
	* CLR2 bin3c_198
* Different HIFI - bin3c_319
	* CLR1 bin3c_346
	* CLR2 bin3c_337

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
# Here is my template for files needed
sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/mag_phase/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g desman/bed_lists/flye4/bin3c.522.scg.bed -o test.pbhap"

# Let's first subsection the bams
mkdir liz_test
vim liz_test/hifi_bin_lists.list
vim liz_test/clr1_bin_lists.list
vim liz_test/clr2_bin_lists.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f b3c_clr1_dastool/clr1.das_DASTool_scaffolds2bin.txt -l liz_test/clr1_bin_lists.list -c 1 -d '\t' | perl -lane 'print $F[0];' |sort | uniq > liz_test/clr1_ctg_lists.list
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f b3c_clr2_dastool/clr2.das_DASTool_scaffolds2bin.txt -l liz_test/clr2_bin_lists.list -c 1 -d '\t' | perl -lane 'print $F[0];' |sort | uniq > liz_test/clr2_ctg_lists.list
cat b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l liz_test/hifi_bin_lists.list -c 1 -d '\t' | perl -lane 'print $F[0];' |sort | uniq > liz_test/hifi_ctg_lists.list

perl -ne 'chomp; print "$_ ";' < liz_test/hifi_ctg_lists.list; echo

samtools faidx flye4.contigs.fasta contig_1416 contig_1837 contig_19687 contig_19927 contig_20371 contig_2504 contig_25423 contig_25426 contig_25428 contig_25435 contig_25437 contig_30414 contig_38194 contig_40897 contig_46424 contig_49499 contig_49519 contig_50421 contig_52351 contig_58524 contig_59402 contig_59409 contig_61593 contig_67291 contig_67293 contig_77971 contig_796 > liz_test/hifi.bincontigs.only.fa

samtools faidx liz_test/hifi.bincontigs.only.fa

samtools view flye4.contigs.fasta.ccs.bam contig_1416 contig_1837 contig_19687 contig_19927 contig_20371 contig_2504 contig_25423 contig_25426 contig_25428 contig_25435 contig_25437 contig_30414 contig_38194 contig_40897 contig_46424 contig_49499 contig_49519 contig_50421 contig_52351 contig_58524 contig_59402 contig_59409 contig_61593 contig_67291 contig_67293 contig_77971 contig_796 > liz_test/hifi.ccs.aligns.sam 
samtools view -bt liz_test/hifi.bincontigs.only.fa.fai liz_test/hifi.ccs.aligns.sam > liz_test/hifi.ccs.aligns.bam

# Let's see if I can automate this for CLR1 and CLR2 now
for i in clr1 clr2; do echo $i; contigs=`perl -ne 'chomp; print "$_ ";' < liz_test/${i}_ctg_lists.list`; echo $contigs; samtools faidx $i.contigs.fasta $contigs > liz_test/$i.bincontigs.only.fa; samtools faidx liz_test/$i.bincontigs.only.fa; samtools view $i.contigs.fasta.ccs.bam $contigs > liz_test/$i.ccs.aligns.sam; samtools view -bt liz_test/$i.bincontigs.only.fa.fai liz_test/$i.ccs.aligns.sam > liz_test/$i.ccs.aligns.bam; done

# Now to just copy the bin beds over
for i in bin3c.44 bin3c.42 bin3c.146 bin3c.120 bin3c.319; do echo $i; cp desman/bed_lists/flye4/$i.scg.bed liz_test/hifi.$i.scg.bed; done
for i in bin3c.13 bin3c.171 bin3c.94 bin3c.346; do echo $i; cp desman/bed_lists/clr1/$i.scg.bed liz_test/clr1.$i.scg.bed; done
for i in bin3c.14 bin3c.154 bin3c.198 bin3c.337; do echo $i; cp desman/bed_lists/clr2/$i.scg.bed liz_test/clr2.$i.scg.bed; done

# Removing the sam files to save space
tar -czvf liz_test.tar.gz liz_test
```

## CheckV analysis

> Ceres:

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/checkv

# Does not need to be done again
checkv download_database ./

# Must be done each time:
export CHECKVDB=/lustre/project/forage_assemblies/sheep_project/complete_flye/checkv-db-v0.6

# Separating contigs for checkv analysis from each assembly
module load samtools; for i in clr1 clr2 clr3 flye4; do echo $i; list=`perl -ne 'chomp; @s = split(/\t/); print "$s[0] ";' < $i.viral.contigs.list`; samtools faidx $i.contigs.fasta $list > $i.viral.contigs.fa; done

for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 25 --mem=75000 -p priority -q msn --wrap="checkv end_to_end $i.viral.contigs.fa ${i}_checkv -t 25"; done

# It worked! Now to tabulate the results and compare them
for i in clr1 clr2 clr3 flye4; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f ${i}_checkv/quality_summary.tsv -c 7 -d '\t' -m; done

```

|Entry          |  HIFI|  CLR1|  CLR2|  CLR3|
|:--------------|-----:|-----:|-----:|-----:|
|Low-quality    |   183|   104|   115|    83|
|High-quality   |    91|    78|    78|    45|
|Medium-quality |    76|    35|    40|    18|
|Not-determined |    22|    30|    28|    13|
|Complete       |    11|    12|     9|     6|


## Split-read analysis

#### TODO: Finish this analysis!

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
# Creating paf files
module load samtools minimap2 bedtools/2.25.0
for i in flye4.contigs.fasta clr1.contigs.fasta clr2.contigs.fasta clr3.contigs.fasta; do echo $i; sbatch -N 1 -n 3 -p priority -q msn --mem=55000 -t 3-0 --wrap="minimap2 -x asm20 $i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel1_CCS.fasta /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta > $i.ccs.paf"; done

# Finding multi-mappers
for i in flye4 clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 2 --mem=55000 -p priority -q msn --wrap="perl ../identify_multimapped_ccs_reads.pl $i.contigs.fasta.ccs.paf $i.ccs.mult.algn.bed $i.ccs.mult.algn.graph $i.ccs.mult.algn.links"; done

# OK, so I just need the graph file for now, but it must be modified to fit within my R paradigm
for i in flye4 clr1 clr2 clr3; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); print "$s[2]\t$ARGV[1]\n";}' $i.ccs.mult.algn.graph $i; done > comparison.ccs.mult.algn.tab

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f comparison.ccs.mult.algn.tab -c 1 -m
|Entry | Value|
|:-----|-----:|
|flye4 | 71372|
|clr1  | 52720|
|clr2  | 52513|
|clr3  | 28531|

# It's clear that there are more multimappers in the flye4 assembly. Is this due to strain separation of contigs?
# Let's check on a per read basis how many times they multimap in one assembly vs another
sbatch -N 1 -n 2 --mem=8000 -p priority -q msn --wrap="perl consolidate_multimapper_table.pl consolidate.multimapper.ccs flye4 flye4.ccs.mult.algn.bed clr1 clr1.ccs.mult.algn.bed clr2 clr2.ccs.mult.algn.bed clr3 clr3.ccs.mult.algn.bed"

# done quickly in R
summary(reads)
                                READ             flye4             clr1
 m54033_180919_161442/10027080/ccs:      1   Min.   : 0.000   Min.   : 0.0000
 m54033_180919_161442/10027085/ccs:      1   1st Qu.: 0.000   1st Qu.: 0.0000
 m54033_180919_161442/10027086/ccs:      1   Median : 0.000   Median : 0.0000
 m54033_180919_161442/10027102/ccs:      1   Mean   : 1.012   Mean   : 0.6996
 m54033_180919_161442/10027123/ccs:      1   3rd Qu.: 2.000   3rd Qu.: 2.0000
 m54033_180919_161442/10027136/ccs:      1   Max.   :33.000   Max.   :31.0000
 (Other)                          :3979872
      clr2             clr3
 Min.   : 0.000   Min.   : 0.0000
 1st Qu.: 0.000   1st Qu.: 0.0000
 Median : 0.000   Median : 0.0000
 Mean   : 0.739   Mean   : 0.6147
 3rd Qu.: 2.000   3rd Qu.: 2.0000
 Max.   :30.000   Max.   :29.0000

# So let's calculate the proportion of reads from flye4 that have few multimapping regions in other assemblies
perl -lane 'if($F[0] eq "READ" || $F[1] == 0 || $F[1] eq ""){next;} if($F[2] + $F[3] + $F[4] < 2 && $F[1] > 2){print "$F[0]";}' < consolidate.multimapper.ccs.reads | wc -l
255732

perl -lane 'if($F[0] eq "READ" || $F[1] == 0 || $F[1] eq ""){next;} if($F[2] + $F[3] + $F[4] < 2 && $F[1] > 2){print "$F[0]\t$F[1]";}' < consolidate.multimapper.ccs.reads > consolidate.multimapper.ccs.reads.list

sbatch -N 1 -n 2 --mem=45000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel1_CCS.fasta -o consolidate.multimapper.sequel1_CCS.fasta -l consolidate.multimapper.ccs.reads.list"
sbatch -N 1 -n 2 --mem=45000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta -o consolidate.multimapper.sequelII_CCS.fasta -l consolidate.multimapper.ccs.reads.list"

cat consolidate.multimapper.sequelII_CCS.fasta consolidate.multimapper.sequel1_CCS.fasta > consolidate.multimapper.all_CCS.fasta

grep 'flye4' consolidate.multimapper.ccs.coords.bed | bedtools sort -i stdin | bedtools merge -i stdin -c 3 -o count > consolidate.multimapper.flye4.multcoords.bed


# Some of these are whole contigs or near the beginnings. I need to assess how many reads are near the start and ends
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $end = $h{$s[0]}; $isstart = 0; $isend = 0; if($s[1] < 100){$isstart = 1;} if($s[2] + 100 > $end){$isend = 1;} $len = $s[2] - $s[1]; print join("\t", @s) . "\t$len\t$isstart\t$isend\n";} close IN;' flye4.contigs.fasta.fai consolidate.multimapper.flye4.multcoords.bed > consolidate.multimapper.flye4.multcoords.table

# Total	AtStart	AtEnd	Both	Neither
# 95893	38485	37269	23485	43624

# Let's work with the Neither category for now and count overlaps with genes
perl -lane 'if(!$F[-1] && !$F[-2]){print $_;}' < consolidate.multimapper.flye4.multcoords.table > consolidate.multimapper.flye4.multcoords.notEnd.bed

# Overlapping with SCGs
grep '>' b3c_flye4_dastool/flye4.das_proteins.faa | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > b3c_flye4_dastool/flye4.das_proteins.bed

cat b3c_flye4_dastool/flye4.das_proteins.faa.archaea.scg b3c_flye4_dastool/flye4.das_proteins.faa.bacteria.scg > b3c_flye4_dastool/flye4.das_proteins.faa.both.scg
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f b3c_flye4_dastool/flye4.das_proteins.bed -c 0 -l b3c_flye4_dastool/flye4.das_proteins.faa.both.scg -d '\t' | perl -lane '@bsegs = split(/_/, $F[0]); $val = pop(@bsegs); $F[0] = join("_", @bsegs); print "$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]";' > b3c_flye4_dastool/flye4.das_proteins.scg.bed

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f b3c_flye4_dastool/flye4.das_proteins.bed -c 0 -l b3c_flye4_dastool/flye4.das_proteins.faa.both.scg -d '\t' -v | perl -lane '@bsegs = split(/_/, $F[0]); $val = pop(@bsegs); $F[0] = join("_", @bsegs); print "$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]";' > b3c_flye4_dastool/flye4.das_proteins.nonscg.bed

# Let's now compare against just the SCGS and the normal genes to see which has better representation
# Non SCGs
bedtools intersect -a consolidate.multimapper.flye4.multcoords.notEnd.bed -b b3c_flye4_dastool/flye4.das_proteins.nonscg.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       545411
        Total Length:           485423528
bedtools intersect -a consolidate.multimapper.flye4.multcoords.notEnd.bed -b b3c_flye4_dastool/flye4.das_proteins.nonscg.bed -v | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       279
        Total Length:           1109710

# SCGs
bedtools intersect -a consolidate.multimapper.flye4.multcoords.notEnd.bed -b b3c_flye4_dastool/flye4.das_proteins.scg.bed -v | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       38425
        Total Length:           430136476
cat consolidate.multimapper.flye4.multcoords.notEnd.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       43624
        Total Length:           563728756

# It's actually not SCGs! Interesting!
# Both sets
cat b3c_flye4_dastool/flye4.das_proteins.scg.bed b3c_flye4_dastool/flye4.das_proteins.nonscg.bed | bedtools intersect -a consolidate.multimapper.flye4.multcoords.notEnd.bed -b stdin | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       560459
        Total Length:           496049014
cat b3c_flye4_dastool/flye4.das_proteins.scg.bed b3c_flye4_dastool/flye4.das_proteins.nonscg.bed | bedtools intersect -a consolidate.multimapper.flye4.multcoords.notEnd.bed -b stdin -v | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       243
        Total Length:           1030552

# testing contig_4157 or bin3c_21828
grep -P 'contig_4157\t' consolidate.multimapper.flye4.multcoords.notEnd.bed > flye4.contig_4157.multimaps.bed
perl -lane 'system("samtools faidx flye4.contigs.fasta $F[0]:$F[1]-$F[2] >> flye4.contig_4157.multimaps.fasta");' < flye4.contig_4157.multimaps.bed

minimap2 -x asm5 clr1.contigs.fasta flye4.contig_4157.multimaps.fasta > clr1.contig_4157.maps.paf

# Not very convincing. Let's find the bins that actually have mash distance correspondance. Like this:
### Refbin  RefComp RefCont RefStrain       clr2_bin        clr2_comp       clr2_cont       clr2_dist       clr2_strain     clr3_bin        clr3_compclr3_cont        clr3_dist       clr3_strain     clr1_bin        clr1_comp       clr1_cont       clr1_dist       clr1_strain
### bin3c_236       100.0   0.0     0       bin3c_293;bin3c_89      94.12;96.08     0.00;7.84       0.01;0.10       0;7     NA      NA      NA      NANA      bin3c_261       82.35   0.00    0.01    0

mkdir split_reads; 
minimap2 -x asm5 b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.236.contigs.fa b3c_clr2_dastool/clr2.das_DASTool_bins/bin3c.293.contigs.fa > split_reads/flye4_b3c_236_vs_clr2_b3c_293.paf


# Actually, I found a better example:
# HiFi assembly bins 451, 452 and 471 map to clr1 bin 451, clr2 327 and clr3 236. All classified as: d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__;s__
# Let's plot the strain haplotype counts

for i in 451 452 471; do echo $i; gunzip b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b flye4.contigs.fasta.ccs.bam -u mag_phase_hr/flye4/bin3c.$i.strain.human_readable_by_pos.txt -i 5000 -o mag_phase_hr/flye4.bin3c.$i.plot -e desman/bed_lists/flye4/bin3c.$i.scg.bed"; done

j=clr1; i=451; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 5000 -o mag_phase_hr/$j.bin3c.451.plot -e desman/bed_lists/$j/bin3c.$i.scg.bed"
j=clr2; i=327; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 5000 -o mag_phase_hr/$j.bin3c.327.plot -e desman/bed_lists/$j/bin3c.$i.scg.bed"
j=clr3; i=236; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 5000 -o mag_phase_hr/$j.bin3c.236.plot -e desman/bed_lists/$j/bin3c.$i.scg.bed"

# And without the bed break lines
j=clr1; i=451; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 1000 -o mag_phase_hr/$j.bin3c.451.nolines"
j=clr2; i=327; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 1000 -o mag_phase_hr/$j.bin3c.327.nolines"
j=clr3; i=236; gunzip b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.gz; samtools faidx b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_${j}_dastool/$j.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b $j.contigs.fasta.ccs.bam -u mag_phase_hr/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 1000 -o mag_phase_hr/$j.bin3c.236.nolines"

for i in 451 452 471; do echo $i; samtools faidx b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.$i.contigs.fa; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.$i.contigs.fa.fai -b flye4.contigs.fasta.ccs.bam -u mag_phase_hr/flye4/bin3c.$i.strain.human_readable_by_pos.txt -i 1000 -o mag_phase_hr/flye4.bin3c.$i.nolines"; done
```


```R
library(dplyr)
mult <- read.delim("comparison.ccs.mult.algn.tab", header=FALSE)
colnames(mult) <- c("Count", "ASM")
mult <- within(mult, quantile <- as.integer(cut(Count, quantile(Count, probs=c(0.1, 0.5, 0.75, 1.0)), include.lowest=TRUE)))

mult %>% group_by(quantile, ASM) %>% summarize(num = n(), sum=sum(Count))
   quantile ASM     num     sum
      <int> <fct> <int>   <int>
 1        1 clr1  30545   37054
 2        1 clr2  30321   36817
 3        1 clr3  15185   18513
 4        1 flye4 35955   44721
 5        2 clr1  10782   47160
 6        2 clr2  10689   46868
 7        2 clr3   5690   25298
 8        2 flye4 16020   71205
 9        3 clr1  11393  812094
10        3 clr2  11503  810891
11        3 clr3   7656  761660
12        3 flye4 19397 1539472

mult <- mutate(mult, Quant = cut(mult$Count, breaks=c(0, 1, 7, 21, 62891)))
mult %>% group_by(Quant, ASM) %>% summarize(num = n(), sum=sum(Count))
`summarise()` regrouping output by 'Quant' (override with `.groups` argument)
# A tibble: 16 x 4
# Groups:   Quant [4]
   Quant         ASM     num     sum
   <fct>         <fct> <int>   <int>
 1 (0,1]         clr1  24036   24036
 2 (0,1]         clr2  23825   23825
 3 (0,1]         clr3  11857   11857
 4 (0,1]         flye4 27189   27189
 5 (1,7]         clr1  17291   60178
 6 (1,7]         clr2  17185   59860
 7 (1,7]         clr3   9018   31954
 8 (1,7]         flye4 24786   88737
 9 (7,21]        clr1   6055   76787
10 (7,21]        clr2   6145   77784
11 (7,21]        clr3   3663   46877
12 (7,21]        flye4  9921  126393
13 (21,6.29e+04] clr1   5338  735307
14 (21,6.29e+04] clr2   5358  733107
15 (21,6.29e+04] clr3   3993  714783
16 (21,6.29e+04] flye4  9476 1413079
```

## Figure 1

I am going to generate a general plot following my summary plot style from my snakemake pipeline. I need to consolidate the CLR datasets and add checkM data to the mix instead of the GC content plot.

Preparing the data first

```bash
echo -e 'LEN\tGC\tKING\tCOMP' > blobtools/table.flye4.extended.lens.tab; perl -ne 'if($_ =~ /^name/){next;} chomp; @s = split(/\t/); print "$s[1]\t$s[2]\t$s[7]\t$s[18]\n";' < flye4_master_table_12_2020.tab >> blobtools/table.flye4.extended.lens.tab

for i in clr1 clr2 clr3; do echo $i; echo -e 'LEN\tGC\tKING\tCOMP' > blobtools/table.$i.extended.lens.tab; perl -ne 'if($_ =~ /^#/){next;} chomp; @s = split(/\t/); $s[18] = ($s[19] <= 5)? $s[18] : "-"; print "$s[1]\t$s[2]\t$s[7]\t$s[18]\n";' < ${i}_master_table_12_2020.tab >> blobtools/table.$i.extended.lens.tab; done


### CLR3 rerun
# /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun
echo -e 'LEN\tGC\tKING\tCOMP' > blobtools/table.clr3.extended.lens.tab; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); if($s[-1] <= 5){$data{$s[0]} = $s[1];}} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; if($_ =~ /^#/){next;} @s = split(/\t/); $val = "-"; if(exists($data{$s[0]})){$val = $data{$s[0]};} print "$s[1]\t$s[2]\t$s[9]\t$val\n";} close IN;' /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun/tables/clr3.contigs.ctg_cov_by_comp.simp.tab blobtools/table.clr3.contigs.blobDB.table.txt > blobtools/table.clr3.extended.lens.tab
```

Now to plot the extended data

```R
library(ggridges)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(doMC)
library(data.table)
library(RColorBrewer)

sumdata <- data.frame(LEN=NA, GC=NA, KING=NA, COMP=NA, ASM=NA)[-1,]
sumdata$KING <- as.factor(sumdata$KING)
sumdata$COMP <- as.factor(sumdata$COMP)
sumdata$ASM <- as.factor(sumdata$ASM)

tp <- read.delim("../blobtools/table.flye4.extended.lens.tab", header=TRUE)
tp <- tp %>% mutate(ASM = c("HIFI"))
tp[tp$KING == "undef", "KING"] <- "no-hit"
sumdata <- bind_rows(sumdata, tp)

tp <- read.delim("../blobtools/table.clr1.extended.lens.tab", header=TRUE)
tp <- tp %>% mutate(ASM = c("CLR1"))
tp[tp$KING == "undef", "KING"] <- "no-hit"
sumdata <- bind_rows(sumdata, tp)
tp <- read.delim("../blobtools/table.clr2.extended.lens.tab", header=TRUE)
tp <- tp %>% mutate(ASM = c("CLR2"))
tp[tp$KING == "undef", "KING"] <- "no-hit"
sumdata <- bind_rows(sumdata, tp)
tp <- read.delim("blobtools/table.clr3.extended.lens.tab", header=TRUE)
tp <- tp %>% mutate(ASM = c("CLR3"))
tp[tp$KING == "undef", "KING"] <- "no-hit"
sumdata <- bind_rows(sumdata, tp)

sumdata$ASM <- as.factor(sumdata$ASM)
sumdata <- sumdata[sumdata$LEN > 1000,]

sumdata$COMP <- as.character(sumdata$COMP)
tmp_wide <- sumdata[sumdata$COMP != "-",] %>% group_by(ASM) %>% summarize(`>90%` = sum(as.numeric(COMP) > 90), `>80%` = sum(as.numeric(COMP) > 80 & as.numeric(COMP) <= 90), `>70%` = sum(as.numeric(COMP) > 70 & as.numeric(COMP) <=80), `>60%` = sum(as.numeric(COMP) > 60 & as.numeric(COMP) <=70))
plot_table <- melt(tmp_wide,id.vars = 'ASM', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table$title <- "CheckM Contig Completeness"
colors <- c(brewer.pal(4, "Oranges"), brewer.pal(4, "Greens"), brewer.pal(4, "Greens"), brewer.pal(4, "Greens"))

plot_table$Completeness <- ordered(plot_table$Completeness, levels=c(">60%", ">70%", ">80%", ">90%"))
sumdata$KING <- ordered(sumdata$KING, levels=c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))


tsumdata <- sumdata %>% group_by(ASM, KING) %>% summarize(SUM=sum(as.numeric(LEN)))
tsumdata$KING <- ordered(tsumdata$KING, levels=c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))

## IT worked but the checkm plot was scrunched and didn't use  my color scheme

p2 <- ggplot(plot_table, aes(ASM, Bins, fill=Completeness)) + geom_bar(stat="identity", position="stack", show.legend = FALSE) + theme_bw() +theme(axis.text.x = element_text(size=11, angle = 45, hjust = 1), axis.title=element_text(size=13, face="bold"), axis.text.y = element_text(size=11), legend.position = "none") + labs(title="Contig Completeness") + ylab("Contigs") + guides(fill=FALSE)

p3 <- ggplot(tsumdata, aes(x=ASM, y=SUM, fill=ASM)) + geom_bar(stat="identity") + scale_color_brewer(palette = "Dark2") + theme_bw() + scale_fill_brewer(palette="Dark2") + theme(axis.title=element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + theme(legend.position = "none") + labs(title="Nucleotides per Kingdom") + facet_wrap(~KING, scales = "free")

p1 <- ggplot(sumdata, aes(y=KING, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=11, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(1000, 10000, 100000, 1000000, 6000000), limits=c(100, 6000000), labels=c("1000", "10,000", "100,000", "1,000,000", "6,000,000")) + xlab(label = "Log10 Contig Lengths (bp)") + ylab(label= "Superkingdom Taxonomy")

pdf(file="figures/hifi_clr_contig_comp_plot.newclr3.pdf", useDingbats=FALSE)
grid.arrange(grobs=list(p3, p2, p1), layout_matrix=rbind(c(1,1,2), c(1,1,2), c(3,3,3)))
dev.off()
```

#### Coverage by cumulative sum plot

> /mnt/c/SharedFolders/metagenomics/tim_sheep/master_tables 

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[1] < 1000){next;} $c = $s[3] + $s[4] + $s[5] + $s[6]; $s[18] = ($s[18] eq "-")? "NA" : $s[18]; print "$s[0]\t$s[1]\t$c\t$s[18]\t$ARGV[1]\n";}' flye4_master_table_12_2020.tab  HIFI | sort -k3n > hifi_ctg_by_cov.tab

for i in clr1 clr2 clr3; do echo $i;  perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $ARGV[1] = uc($ARGV[1]); <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[1] < 1000){next;} $c = $s[3] + $s[4] + $s[5] + $s[6]; $s[18] = ($s[18] eq "-")? "NA" : $s[18]; print "$s[0]\t$s[1]\t$c\t$s[18]\t$ARGV[1]\n";}' ${i}_master_table_12_2020.tab $i | sort -k3n > ${i}_ctg_by_cov.tab; done

for i in clr1_ctg_by_cov.tab clr2_ctg_by_cov.tab clr3_ctg_by_cov.tab hifi_ctg_by_cov.tab; do echo $i; perl -lane 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]";' <$i >> combined_ctg_by_cov.tab ; done

##### NOTE THIS WAS DONE ON CERES
perl -e '%l = ("CLR1" => "clr1", "CLR2" => "clr2", "CLR3" => "clr3", "HIFI" => "flye4"); $h = <>; chomp($h); $h =~ s/\r//g; print "$h\tcomp\tcont\n"; while($c = <>){chomp $c; $c =~ s/\r//g; @s = split(/\t/, $c); $f = $l{$s[0]}; print "$c"; open(IN, "< b3c_" . $f . "_dastool/" . $f . ".das_DASTool_summary.txt"); while($t = <IN>){chomp $t; @n = split(/\t/, $t); if ($n[0] eq $s[1]){print "\t$n[11]\t$n[12]\n"; last;}} close IN; }' < combined_per_bin_averagecount.tab > combined_per_bin_averagecount.pluscomp.tab
```

```R
library(ggplot2)
library(dplyr)

data <- read.delim("combined_ctg_by_cov.tab", header=TRUE)
data.csum <- data %>% group_by(ASM) %>% mutate(CSUM= cumsum(as.numeric(LEN)))
data.csum %>% group_by(ASM) %>% summarise(max=max(CSUM))
`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 4 x 2
  ASM          max
  <chr>      <dbl>
1 CLR1  2984846227
2 CLR2  3008462117
3 CLR3  1675674347
4 HIFI  3424468269


ggplot(data=data.csum, aes(x=COV, y=CSUM, color=ASM)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + xlim(c(0,50)) + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(1000000000, 2000000000, 3000000000), labels=c("1.0", "2.0", "3.0")) + ylab("Cumulative length (Gbp)") + xlab("X Coverage")

# Plot saved as pdf: cumulative_length_by_xcov.pdf

# Now for all contigs that had completeness estimates
data.comp <- na.omit(data[data$COMP > 60,]) %>% group_by(ASM) %>% mutate(CSUM = as.numeric(cumsum(LEN)))
data.comp %>% summarize(cmax = max(CSUM), count = n())
`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 4 x 3
  ASM        cmax count
  <chr>     <dbl> <int>
1 CLR1  246569508   133
2 CLR2  247760547   132
3 CLR3  127249599    70
4 HIFI  477237589   239

ggplot(data=data.comp, aes(x=COV, y=CSUM, color=ASM)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(100000000, 200000000, 300000000, 400000000, 500000000), labels=c("100", "200", "300", "400", "500")) + ylab("Cumulative length (Mbp)") + xlab("X Coverage") + xlim(c(0,100)) + ggtitle("Contigs greater than 60% Complete")

# Plot saved as pdf: cumlength_vs_xcov_ctgsgt60_comp.pdf

# Now, without the 60% complete filter. Just less than 5% redundancy
data.all <-na.omit(data) %>% group_by(ASM) %>% mutate(CSUM = as.numeric(cumsum(LEN)))
data.all %>% summarize(cmax = max(CSUM), count = n())
`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 4 x 3
  ASM         cmax count
  <chr>      <dbl> <int>
1 CLR1   982610924  1201
2 CLR2   993465163  1201
3 CLR3   699592799  1201
4 HIFI  1242372220  1112

ggplot(data=data.all, aes(x=COV, y=CSUM, color=ASM)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(0, 400000000, 800000000, 1200000000), labels=c("0", "400", "800", "1200")) + ylab("Cumulative length (Mbp)") + xlab("X Coverage") + xlim(c(0,100)) + ggtitle("Contigs less than 5% redundant")
# Plot saved as pdf: cumlength_vs_xcov_ctgslt5_redund.pdf

# Now to see if we can reverse by Xcov and plot cumsum that way instead
data.rev <- data %>% arrange(desc(COV)) %>% group_by(ASM) %>% mutate(CSUM = cumsum(as.numeric(LEN)))

ggplot(data=data.rev, aes(x=COV, y=CSUM, color=ASM)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(1000000000, 2000000000, 3000000000), labels=c("1.0", "2.0", "3.0")) + ylab("Cumulative length (Gbp)") + xlab("X Coverage") + scale_x_reverse(lim=c(100, 0)

# Pdf was saved: cumsum_reverse_xcov_plot.pdf
```

Now to try to do this for the bins. This can be done with the information from the master tables, but I will have to group by the common bins to generate the statistics I need. I generated the file with aggregate bins in an ipython notebook and now I will run it through my R code.

```R
library(ggplot2)
library(dplyr)
library(ggridges)
setwd("C:/SharedFolders/metagenomics/tim_sheep/master_tables/")

data <- read.delim("combined_per_bin_averagecount.tab", header=TRUE)
data.rev <- data %>% arrange(desc(avgCov)) %>% group_by(Assembly) %>% mutate(CSUM = cumsum(as.numeric(totlen)))

data.rev %>% summarize(cmax = max(CSUM), count = n())
`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 4 x 3
  Assembly       cmax count
  <chr>         <dbl> <int>
1 CLR1     1292130032   557
2 CLR2     1315962075   552
3 CLR3      561776173   279
4 HIFI     1505258010   695

ggplot(data=data.rev, aes(x=avgCov, y=CSUM, color=Assembly)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(500000000, 1000000000, 1500000000), labels=c("0.5", "1.0", "1.5")) + ylab("Cumulative length (Gbp)") + xlab("average X Coverage") + scale_x_reverse(lim=c(100, 0))

# Pdf was saved as: cumulative_sum_bin_avgcov.pdf

data <- read.delim("combined_per_bin_averagecount.pluscomp.tab", header=TRUE)
data.fix <- data[data$cont < 5,]

data.fix$Assembly <- as.factor(data.fix$Assembly)
data.top <- data.fix[data.fix$comp >= 80,]
ggplot(data=data.top, aes(x=comp, y=avgCov, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Dark2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
# File was saved as bin_top80_cov_by_comp.pdf

data.ninety <- data.top[data.top$comp >=95,]
ggplot(data=data.ninety, aes(y=Assembly, x=avgCov, fill=Assembly)) + geom_density_ridges() + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
# File was saved as bin_top95_cov_ridgeline.pdf

data.ninety[data.ninety$Assembly == "HIFI" & data.ninety$avgCov < 10 & data.ninety$count == 1,]
     Assembly    Bin3cBin count  totlen  avggc avgCov stdCov  comp cont
1397     HIFI bin3c.10610     1 1525910 0.3021 6.5267     NA 96.08    0
1539     HIFI bin3c.21828     1 1769107 0.2890 5.5499     NA 98.04    0
1549     HIFI bin3c.22331     1 1901240 0.2990 3.9807     NA 98.04    0

data.ninety[data.ninety$Assembly == "HIFI" & data.ninety$avgCov < 10 & data.ninety$comp > 99.99 & data.ninety$cont == 0.0,]
     Assembly  Bin3cBin count  totlen     avggc   avgCov    stdCov comp
1449     HIFI bin3c.152     4 3291527 0.5027250 4.687100 1.6213751  100
1454     HIFI bin3c.156     7 3242186 0.5695000 3.433200 0.3331985  100
1530     HIFI bin3c.211     5 2981387 0.5483400 6.264920 0.7890848  100
1536     HIFI bin3c.217    18 2969005 0.5064889 5.268000 1.2365500  100
1563     HIFI bin3c.236     3 2886616 0.5652000 6.497267 1.0095671  100
1607     HIFI bin3c.271     3 2765004 0.5570000 4.715433 2.8267415  100
1698     HIFI  bin3c.35     8 5011542 0.4541000 4.196600 1.1767137  100
1738     HIFI bin3c.387     3 2317290 0.4104333 7.068800 0.4549975  100
1877     HIFI bin3c.506     5 1880153 0.4752400 6.177640 1.0878781  100
1882     HIFI bin3c.510     4 1861749 0.4220500 9.541050 1.3099928  100
     cont
1449    0
1454    0
1530    0
1536    0
1563    0
1607    0
1698    0
1738    0
1877    0
1882    0


### New CLR3 bins
tp <-read.delim("clr1.contigs.bin_cov_by_comp.len.tab")
data <- tp
tp <-read.delim("clr2.contigs.bin_cov_by_comp.len.tab")
data <- bind_rows(data, tp)
tp <-read.delim("clr3.contigs.bin_cov_by_comp.len.tab")
data <- bind_rows(data, tp)
tp <-read.delim("flye4.contigs.bin_cov_by_comp.len.tab")
data <- bind_rows(data, tp)

data$Assembly <- as.factor(data$Assembly)

data.rev <- data %>% arrange(desc(Coverage)) %>% group_by(Assembly) %>% mutate(CSUM = cumsum(as.numeric(Length)))
data.rev %>% summarize(cmax = max(CSUM), count = n())
# A tibble: 4 x 3
  Assembly            cmax count
* <fct>              <dbl> <int>
1 clr1.contigs  1299194391   563
2 clr2.contigs  1311659869   551
3 clr3.contigs  1253996764   534
4 flye4.contigs 1504965682   694

ggplot(data=data.rev, aes(x=avgCov, y=CSUM, color=Assembly)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.y = element_blank()) + scale_y_continuous(breaks=c(500000000, 1000000000, 1500000000), labels=c("0.5", "1.0", "1.5")) + ylab("Cumulative length (Gbp)") + xlab("average X Coverage") + scale_x_reverse(lim=c(100, 0))

# plot was saved as: assembly_cumulative_sum_update.pdf

data.fix <- data[data$Contamination < 5,]
data.top <- data.fix[data.fix$Completeness >= 80,]

ggplot(data=data.top, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Dark2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
# File was saved as bin_top80_cov_by_comp.nocont.pdf

data.ninety <- data.top[data.top$Completeness >=95,]
ggplot(data=data.ninety, aes(y=Assembly, x=Coverage, fill=Assembly)) + geom_density_ridges() + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
# File was saved as bin_top95_cov_ridgeline.nocont.pdf

data.ninety <- data.top[data.top$Completeness >=90,]
ggplot(data=data.ninety, aes(y=Assembly, x=Coverage, color=Assembly)) + geom_violin(aes(color="black", fill=Assembly), alpha=0.4) + geom_jitter(height = 0.3) + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
# Figure was saved as top_90_jitter_violin.pdf
```

#####TODO: Find a good example for associating with the multimappers and showing the divergence between CLR and HIFI assemblies

HiFi assembly bins 451, 452 and 471 map to clr1 bin 451, clr2 327 and clr3 236. All classified as: d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__;s__

#### Modification of contig coverage plot


OK, I'm going to revise the coverage estimates based on the total number of bases mapped to each contig and add that to our plot.

```bash
for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 2 --mem=14000 -p priority -q msn --wrap="perl paf_base_count.pl $i.contigs.fasta.ccs.paf $i.contigs.fasta.ccs.alignedbases.tab"; done

# On second thought, this is probably super-biased. Let's try running that JGI estimation program
module load metabat/2.12.1
for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 2 --mem=14000 -p priority -q msn --wrap="jgi_summarize_bam_contig_depths --outputDepth $i.contigs.fasta.ccs.alignedbases.tab $i.contigs.fasta.ccs.bam"; done

# No, that just gave me the average. I need to do this as a percentage of the total CCS sequence
# Getting total bases of CCS reads
cat /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel1_CCS.fasta.fai /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta.fai | perl -e '$c = 0; while(<STDIN>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
255708235755

for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 2 --mem=20000 -p priority -q msn --wrap="perl generate_percentage_total_coverage.pl $i.contigs.fasta.ccs.bam 255708235755 $i.contigs.ccs.percdataset.tab"; done

# Now to convert this to something I can plot in R
for i in clr1.contigs.ccs.percdataset.tab clr3.contigs.ccs.percdataset.tab clr2.contigs.ccs.percdataset.tab flye4.contigs.ccs.percdataset.tab; do echo $i; perl -e '<>; while(<>){chomp; @s = split(/\t/); $s[1] *= 100; print join("\t", @s) . "\n";}' < $i > $i.rfmt; done

echo -e "ASM\tLEN\tPERC" > combined_percdataset.tab; for i in clr1 clr2 clr3 flye4; do echo $i; perl -e 'chomp(@ARGV); if($ARGV[0] eq "flye4"){$ARGV[0] = "HIFI";} open(IN, "< $ARGV[2]"); %c; while(<IN>){chomp; @s = split(/\t/); $c{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $l = $c{$s[0]}; print "$ARGV[0]\t$l\t$s[1]\n";} close IN;' $i $i.contigs.ccs.percdataset.tab.rfmt $i.contigs.fasta.fai>> combined_percdataset.tab; done
```

```R
library(ggplot2)
library(dplyr)
data <- read.delim("combined_percdataset.tab", header=TRUE)
data.rev <- data %>% arrange(desc(LEN)) %>% group_by(ASM) %>% mutate(CSUM = cumsum(as.numeric(PERC)), CLEN = cumsum(as.numeric(LEN)))

data.rev %>% summarize(cmax = max(CSUM), count = n())
`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 4 x 3
  ASM    cmax count
  <fct> <dbl> <int>
1 clr1  105.  48408
2 clr2  105.  48748
3 clr3   97.2 24415
4 HIFI  103.  55711

pdf(file="combined_percdataset_ccs_reads.pdf", useDingbats=FALSE)
ggplot(data=data.rev, aes(x=CLEN, y=CSUM, color=ASM)) + geom_line() + theme_bw() + scale_color_brewer(palette="Dark2") + theme(axis.ticks.x = element_blank()) + scale_x_continuous(breaks=c(500000000, 1000000000, 1500000000, 2000000000, 3000000000), labels=c("0.5", "1.0", "1.5", "2.0", "3.0")) + ylab("Cumulative percentage of mapped CCS reads") + xlab("Cumulative length of assembly (Gbp)")
dev.off()
```

#### Contig sequence comparisons

I am going to use SASinspector on the circular HIFI contigs compared to the CLR datasets to identify mapped and missing portions. The goal is to see if there are common elements that cause the breakdown of contigs from the CLR assemblies against comparable HIFI assemblies. 

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/saspector

# Now to identify the nearest hit for all of the HIFI HQ contigs
mkdir saspector
mkdir saspector/hq_hifi_contigs; perl -lane 'if($F[1] >= 80){print($_); system("samtools faidx flye4.contigs.fasta $F[0] > saspector/hq_hifi_contigs/$F[0].fna");}' < "/lustre/project/forage_assemblies/sheep_project/complete_flye/flye4_checkm_contigs.tab"

# I had to reduce the size of the fasta files to avoid hitting a mash error!
module load samtools; for i in clr1.contigs.fasta clr2.contigs.fasta clr3.contigs.fasta; do echo $i; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0].fai"); while(<IN>){chomp; @s = split(/\t/); if($s[1] > 100000){system("samtools faidx $ARGV[0] $s[0] >> $ARGV[0].gt100k.fa");}} close IN;' $i; done

for i in clr1.contigs.fasta.gt100k.fa clr2.contigs.fasta.gt100k.fa clr3.contigs.fasta.gt100k.fa; do echo $i; sbatch -N 1 -n 2 -p priority -q msn --mem=35000 --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.2/mash sketch -i -k 21 -s 100000 $i -o saspector/$i"; done

mkdir saspector/hq_hifi_mash; 
for i in saspector/hq_hifi_contigs/*.fna; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=15000 --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.2/mash sketch -i -k 21 -s 100000 $i -o saspector/hq_hifi_mash/$name"; done

mkdir saspector/clr1; mkdir saspector/clr2; mkdir saspector/clr3; 
for i in saspector/hq_hifi_mash/*.msh; do name=`basename $i | cut -d'.' -f1`; echo $name; for j in clr1 clr2 clr3; do echo $j;  sbatch -N 1 -n 2 --mem=15000 -p priority -q msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.2/mash dist -d 0.1 saspector/${j}.contigs.fasta.gt100k.fa.msh $i > saspector/$j/$name.$j.dist.tab"; done; done

# Now to try to queue it all up
for i in clr1 clr2 clr3; do echo $i; for j in saspector/$i/*.tab; do name=`basename $j | cut -d'.' -f1`; echo $name; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); system("samtools faidx $ARGV[1] $s[0] >> saspector/$ARGV[2].fa");}' $j $i.contigs.fasta ${i}_${name}; sbatch -N 1 -n 2 -p priority -q msn --mem=35000 --wrap="SASpector -dir saspector/clr1/$name saspector/hq_hifi_contigs/$name.fna saspector/${i}_${name}.fa"; done ; done

# That didn't work because the path for the reference genome was incorporated in file names!
# I'll have to queue in each directory specifically
cd saspector/clr1
mv ../*.fa ./
cp ../hq_hifi_contigs/*.fna ./

for j in *.dist.tab; do name=`echo $j | cut -d'.' -f1`; echo $name; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); system("samtools faidx $ARGV[1] $s[0] >> $ARGV[2].fa");}' $j ../../clr1.contigs.fasta clr1_${name}; sbatch -N 1 -n 2 -p priority -q msn --mem=35000 --wrap="SASpector -dir $name -p $name $name.fna clr1_${name}.fa"; done

cd ../

for i in clr2 clr3; do echo $i; cd $i; cp ../hq_hifi_contigs/*.fna ./; for j in *.dist.tab; do name=`echo $j | cut -d'.' -f1`; echo $name; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); system("samtools faidx $ARGV[1] $s[0] >> $ARGV[2].fa");}' $j ../../$i.contigs.fasta ${i}_${name}; sbatch -N 1 -n 2 -p priority -q msn --mem=35000 --wrap="SASpector -dir $name -p $name $name.fna ${i}_${name}.fa"; done; cd ../; done

# Now to turn this all into something that's actionable
for i in clr1 clr2 clr3; do echo $i; echo -e "ASM\tHIFIContig\tGCContent\tLength\tNumberMappedRegions\tNumberUnmappedRegions\tFilteredUnmappedRegions\tFractionMapped\tFractionUnmapped" > ${i}_reference_summary.tab; for j in $i/*/*_referencesummary.tsv; do name=`basename $j | cut -d '_' -f1,2`; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; $s = <IN>; print "$ARGV[1]\t$ARGV[2]\t$s"; close IN;' $j $i $name; done >> ${i}_reference_summary.tab; done

# And a bed file for all of the unmapped regions
for i in clr1 clr2 clr3; do echo $i; for j in $i/*/*_unmapsummary.tsv; do name=`basename $j | cut -d '_' -f1,2`; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; while($s = <IN>){chomp $s; @F = split(/\t/, $s); @bsegs = split(/[_:]/, $F[0]); $chr = "$bsegs[0]\_$bsegs[1]"; print "$chr\t$bsegs[2]\t$bsegs[3]\t$ARGV[1]\n";} close IN;' $j $i; done > ${i}_unmapped_regions.bed; done

# Now to see if anything matches up
module load bedtools/2.25.0
# Missing in just clr1 & clr2
cat *.bed | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o distinct -delim ';' | grep 'clr1;clr2' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       584
        Total Length:           32263237
        Length Average:         55245.2688356164
        Length Median:          2187.5
        Length Stdev:           176440.311961825
        Smallest Length:        102
        Largest Length:         1775719

# Missing in all three CLR
cat *.bed | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o distinct -delim ';' | grep 'clr1;clr2;clr3' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       299
        Total Length:           26401940
        Length Average:         88300.8026755853
        Length Median:          2592
        Length Stdev:           231826.181477718
        Smallest Length:        102
        Largest Length:         1775719

# Wow! What is the region that is missing 1775719?
cat *.bed | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o distinct -delim ';' | grep 'clr1;clr2;clr3' | perl -lane 'if($F[2] - $F[1] > 1000000){print $_;}'
contig_1656     517199  1873256 clr1;clr2;clr3
contig_22264    295902  1564788 clr1;clr2;clr3
contig_3316     551730  2327449 clr1;clr2;clr3
contig_3340     902674  2514109 clr1;clr2;clr3
contig_5786     185504  1255027 clr1;clr2;clr3
```

## Read classification rates

I need to calculate how well each assembly "covers" the short-read and CCS data. 

```bash
module load kraken2/2.1.1
# download NCBI taxonomy
sbatch -N 1 -n 2 --mem=20000 -p priority -q msn -t 1-0 --wrap="kraken2-build --download-taxonomy --db NCBI"

# Create latest NCBI Kraken database
sbatch -N 1 -n 32 --mem=300000 -p priority -q msn -t 1-0 --wrap="kraken2-build --standard --threads 32 --db REFSEQ"

# Hmm, to make custom DBs, I need to add the NCBI taxonomy link to the fasta header of each bin
# I can do this by generating consensus taxID assignment using the blobtools diamond taxify output
python3 ~/python_toolchain/metagenomics/ncbiTaxIDFromBlobtools.py -c /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.taxids -t blobtools/taxify.flye4.contigs.diamondout.tsv.taxified.out -d b3c_flye4_dastool/flye4.das_DASTool_bins -o flye4.taxconsensus.tab

perl -lane 'print $F[2];' < flye4.taxconsensus.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   695
Sum:    3877
Minimum 1
Maximum 65
Average 5.578417
Median  4
Standard Deviation      5.874710
Mode(Highest Distributed Value) 1

# OK, most bins are able to be classified using this method, but the median is that at least four contigs differ from the consensus
for i in clr1 clr2 clr3; do echo $i; python3 ~/python_toolchain/metagenomics/ncbiTaxIDFromBlobtools.py -c /project/rumen_longread_metagenome_assembly/assemblies/protists/uniprot_ref_proteomes.taxids -t blobtools/taxify.$i.contigs.diamondout.tsv.taxified.out -d b3c_${i}_dastool/$i.das_DASTool_bins -o $i.taxconsensus.tab; done

for i in clr1 clr2 clr3; do echo $i; perl -lane 'print $F[2];' < $i.taxconsensus.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl; done
clr1
total   557
Sum:    4139
Minimum 1
Maximum 57
Average 7.430880
Median  6
Standard Deviation      6.290708
Mode(Highest Distributed Value) 5

clr2
total   552
Sum:    4456
Minimum 1
Maximum 65
Average 8.072464
Median  6
Standard Deviation      7.745158
Mode(Highest Distributed Value) 1

clr3
total   279
Sum:    2317
Minimum 1
Maximum 87
Average 8.304659
Median  6
Standard Deviation      8.890846
Mode(Highest Distributed Value) 1

# OK, now to edit the fastas so that they can be incorporated into kraken databases
for i in flye4 clr1 clr2 clr3; do echo $i; mkdir kraken_${i}_bins; perl -e 'use File::Basename; chomp(@ARGV); open(IN, "<$ARGV[0]"); %ctg; while(<IN>){chomp; @s = split(/\t/); $c = 3; $t = "NA"; while($t eq "NA"){@csegs = split(/:/, $s[$c]); $c++; $t = $csegs[0];} if($t eq "NA"){$t = 1; print "$s[0]\t$t\n";} $ctg{$s[0]} = $t;} close IN;   @f = `ls $ARGV[1];`; foreach $j (@f){open(IN, "< $j"); open(OUT, "> $ARGV[2]/" . basename($j)); $tid = $ctg{basename($j)}; while(<IN>){if($_ =~ /^>/){chomp; $_ .= "\|kraken:taxid\|$tid\n";} print {OUT} $_;} close IN;}' $i.taxconsensus.tab b3c_${i}_dastool/$i.das_DASTool_bins kraken_${i}_bins; done

# Now I should be able to generate Kraken databases!
for i in flye4 clr1 clr2 clr3; do echo $i; sbatch -t 0-8 create_kraken_customdb.sh kraken_${i}_db kraken_${i}_bins; done

# OK, now to start classifying reads
# I'm going to feed each lane of the illumina data in separately
mkdir kraken_classifications
for i in clr1 clr2 clr3 flye4; do echo $i; dbname=kraken_${i}_db; sbatch -N 1 -n 2 --mem=45000 -p priority -q msn -t 1-0 --wrap="kraken2 --quick --db $dbname --output kraken_classifications/$i.L1ILMN.kraken.out --paired /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L001_R1_001.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L001_R2_001.fastq.gz; kraken2 --quick --db $dbname --output kraken_classifications/$i.L2ILMN.kraken.out --paired /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L002_R1_001.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L002_R2_001.fastq.gz; kraken2 --quick --db $dbname --output kraken_classifications/$i.L3ILMN.kraken.out --paired /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L003_R1_001.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L003_R2_001.fastq.gz; kraken2 --quick --db $dbname --output kraken_classifications/$i.L4ILMN.kraken.out --paired /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L004_R1_001.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_S4_L004_R2_001.fastq.gz;"; done

for i in kraken_classifications/*.out; do echo $i; sbatch -N 1 -n 2 --mem=8000 -p priority -q msn -t 1-0 --wrap="python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 0 -d '\t' > $i.count"; done

for i in kraken_classifications/*.count; do name=`basename $iname; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1] * 1;} close IN; $ratio = $h{"C"} / ($h{"C"} + $h{"U"}); print "$ARGV[1]\t$ratio\n";' $i $name; done

# Output was pretty trash. 50% classification rate, but probably due to read issues in the illumina data
```

## MMGenome2 plots

Just to keep Mads happy, I want to compare coverage profiles between the CLR datasets, the CCS reads and short reads for each of our assembled bins.

I already have CCS read alignment depth of coverage:
* *.contigs.fasta.ccs.alignedbases.tab

I also have average depth of coverage from short reads:
* blobtools/*.contigs.Lib101_L1.bam.cov

But I have not aligned individual clr datasets to each assembly. Let's do that now

```bash
module load minimap2 samtools
mkdir minimap_clr_bams

for i in clr.1 clr.2 clr.3; do echo $i; for j in clr1 clr2 clr3 flye4; do echo $j; sbatch -N 1 -n 3 -p priority -q msn --mem=45000 -t 2-0 --wrap="minimap2 -ax map-pb -R '@RG\tID:$i\tSM:$i' $j.contigs.fasta /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/*_${i}.fastq.gz | samtools sort -T $j.$i -o minimap_clr_bams/$j.$i.sorted.bam -"; done; done

# I renamed the files to remove the extraneous '.' in the clr names
module load metabat/2.12.1
for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 2 --mem=24000 -p priority -q msn -t 1-0 --wrap="jgi_summarize_bam_contig_depths --percentIdentity 85 --noIntraDepthVariance --minMapQual 1 --outputDepth $i.jgidepth.alllr.error.tab $i.contigs.fasta.ccs.bam minimap_clr_bams/$i.clr1.sorted.bam minimap_clr_bams/$i.clr2.sorted.bam minimap_clr_bams/$i.clr3.sorted.bam"; done

for i in *.jgidepth.alllr.error.tab; do echo $i; name=`basename $i | cut -d'.' -f1`; echo $name; perl -lane 'if($F[0] eq "contigName"){for($x = 3; $x < scalar(@F); $x++){@bsegs = split(/\./, $F[$x]); $bsegs[1] = ($bsegs[1] eq "contigs")? "hifi" : $bsegs[1]; $F[$x] = $bsegs[1];}} print join("\t", @F[0,3,4,5,6]);' < $i > $name.mmgenome2.input.tab; done

# Generating the illumina coverage
for i in clr1 clr2 clr3 flye4; do echo $i; perl -lane 'if($F[0] =~ /##/){next;}elsif($F[0] =~ /#{1}/){print "contigName\tilmn";}else{print join("\t", @F[0,2]);}' < blobtools/$i.contigs.Lib101_L1.bam.cov > $i.mmgenome.ilumina.tab; done

# Now to just plot out everything
```

> C:/SharedFolders/metagenomics/tim_sheep/mmgenome/

```R
library(mmgenome2)
library(dplyr)

ccsdata <- read.delim("flye4.mmgenome2.input.tab", header=TRUE)
ilmdata <- read.delim("flye4.mmgenome.ilumina.tab", header=TRUE)

combined <- left_join(ccsdata, ilmdata, by="contigName")
mm <- mmload("flye4.contigs.fasta", coverage=combined)

for (c in c("clr1", "clr2", "clr3", "ilmn")){
filename <- paste0("hifi.", c, ".covplot.pdf")
pdf(file=filename, useDingbats = FALSE)
mt <- mmplot(mm, x = "cov_hifi", y = paste0("cov_", c), color_by = "gc", x_scale="log10", y_scale = "log10")
print(mt)
dev.off()
}
```

## MagPhase gene-level run

OK, Liz came through on the improvements! I need to install the newest version of the tool and try to run through all of my bins to compare.

Most of the tools are meant to be run with only path level access. I will try that first.

```bash
module load miniconda/3.6 samtools
conda activate /KEEP/rumen_longread_metagenome_assembly/desman
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/sequence/
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/rarefaction/

# Testing it out first
sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g desman/bed_lists/flye4/bin3c.522.scg.bed -o test.pbhap.magphased"

# It works! Let's queue up the rest
for j in flye4 clr1 clr2 clr3; do echo $j; mkdir mag_phase/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a $j.contigs.fasta -b $j.contigs.fasta.ccs.bam -g $i -o mag_phase/$j/$name.strain"; done; done

# Checking for bins that are totally resolved
for i in clr1 clr2 clr3 flye4; do echo $i; ls mag_phase/$i/*.NO_SNPS_FOUND | wc -l; done
clr1
202
clr2
186
clr3
39
flye4
311

# Trying to extract some means of estimating strain heterogeneity
for i in clr1 clr2 clr3 flye4; do grep 'contig' mag_phase/$i/*.human_readable.txt | perl -e '%d; while(<>){chomp; @s = split(/\t/); ($h) = $s[0] =~ m/(bin3c.\d{1,4})./; $d{$h}->{$s[5]} += 1;} foreach $bin (keys(%d)){foreach $alt (keys(%{$d{$bin}})){if($alt eq "base"){next;} print "$bin\t$alt\t" . $d{$bin}->{$alt} . "\n"}}' > mag_phase/$i.alt_counts.tab; done

for i in clr1 clr2 clr3 flye4; do echo $i; cat mag_phase/$i.alt_counts.tab | cut -f2 | sort | uniq -c; done
clr1
    352 ALT0
     94 ALT1
      7 ALT2
    352 REF
clr2
    365 ALT0
    122 ALT1
      9 ALT2
    365 REF
clr3
    240 ALT0
    146 ALT1
     31 ALT2
    240 REF
flye4
    379 ALT0
     47 ALT1
    379 REF

## New version -- only major change is the output file format
mkdir mag_phase_hr
for j in flye4 clr1 clr2 clr3; do echo $j; mkdir mag_phase_hr/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a $j.contigs.fasta -b $j.contigs.fasta.ccs.bam -g $i -o mag_phase_hr/$j/$name.strain"; done; done

# Now to try to consolidate output and generate usable stats for plotting in R
mkdir mag_phase_hr/consolidated
for j in flye4 clr1 clr2 clr3; do echo $j; mkdir mag_phase_hr/consolidated/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/python_toolchain/metagenomics/calcMagPhaseOutputVals.py -f mag_phase_hr/$j -p $name -d b3c_${j}_dastool/${j}.das_DASTool_summary.txt -o mag_phase_hr/consolidated/$j/$name; done; done

# Let's try to pull some base stats out of this run
for i in clr1 clr2 clr3 flye4; do echo $i; cat mag_phase_hr/consolidated/$i/*.short | perl -lane 'print "$F[1]";' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl ; done
clr1
total   557
Sum:    2306
Minimum 0
Maximum 60
Average 4.140036
Median  3
Standard Deviation      6.148211
Mode(Highest Distributed Value) 0

clr2
total   552
Sum:    2445
Minimum 0
Maximum 63
Average 4.429348
Median  4
Standard Deviation      5.881967
Mode(Highest Distributed Value) 0

clr3
total   534
Sum:    2223
Minimum 0
Maximum 68
Average 4.162921
Median  3
Standard Deviation      6.798469
Mode(Highest Distributed Value) 0

flye4
total   695
Sum:    1975
Minimum 0
Maximum 29
Average 2.841727
Median  0
Standard Deviation      4.014489
Mode(Highest Distributed Value) 0

# Calculating stats from long form and skipping smaller haplotypes
# Number of rows, Max phased distance, max haplotype size
# I found out why these are so long -- the SCG coordinates overlap! 
for i in clr1 clr2 clr3 flye4; do echo $i; cat mag_phase_hr/consolidated/$i/*.long | perl -e '$c = 0; $m = 0; $l = 0; while(<STDIN>){chomp; @s = split(/\t/); if($s[2] < 4){next;} $c++; if($s[6] - $s[5] > $m){$m = $s[6] - $s[5];}  if($s[2] > $l){$l = $s[2];} } print "$c\t$m\t$l\n";' ; done
clr1
22826   463082  438
clr2
24018   792893  378
clr3
23583   493333  345
flye4
12751   336899  309

# OK, let's combine and see if we can generate the same plot for the strain phasing.
for i in clr1 clr2 clr3 flye4; do echo $i; cat mag_phase_hr/consolidated/$i/*.short > mag_phase_hr/consolidated/$i.consolidated.short.tab; done

# OK, let's queue it up!
python3 ~/python_toolchain/metagenomics/createMAGstrainAssoc.py -r mag_phase_hr/consolidated/flye4.consolidated.short.tab -o mag_phase_hr/consolidated/flye4_strain_associations -n clr1 -a gtdbtk_bins/flye4_clr1_associations.tab -s mag_phase_hr/consolidated/clr1.consolidated.short.tab -n clr2 -a gtdbtk_bins/flye4_clr2_associations.tab -s mag_phase_hr/consolidated/clr2.consolidated.short.tab -n clr3 -a gtdbtk_bins/flye4_clr3_associations.tab -s mag_phase_hr/consolidated/clr3.consolidated.short.tab


# Averages
for i in clr1 clr2 clr3 flye4; do echo $i; cat mag_phase_hr/consolidated/$i/*.long | perl -e '$c = 0; $m = 0; $l = 0; while(<STDIN>){chomp; @s = split(/\t/); if($s[2] < 4){next;} $c++; $m += $s[6] - $s[5];  $l += $s[2];} $m /= $c; $l /= $c; print "$c\t$m\t$l\n";' ; done
clr1
22826   1092.63970910365        28.1955226496101
clr2
24018   906.666208676826        29.1952285785661
clr3
30864   761.680728356661        52.0766264904095
flye4
12751   1174.2172378637 19.7565681123049

for i in clr1 clr2 flye4; do echo $i; cat ../mag_phase_hr/consolidated/$i/*.long | perl -e '%data; $count=0; while(<>){chomp; @F = split(/\t/); if(length($F[1]) > 1){$count++; $t = $F[0] . $F[4] . $F[5]; $data{$t} += 1;}} $sum = 0; foreach $v (values(%data)){$sum += $v;} print "$count\t$sum\t" . ($sum / scalar(keys(%data))) . "\n"; '; done
clr1
27558   27558   5.06674020959735
clr2
28902   28902   4.84688914975683
flye4
17481   17481   4.38780120481928
cat magphase/consolidated/clr3/*.long | perl -e '%data; $count=0; while(<>){chomp; @F = split(/\t/); if(length($F[1]) > 1){$count++; $t = $F[0] . $F[4] . $F[5]; $data{$t} += 1;}} $sum = 0; foreach $v (values(%data)){$sum += $v;} print "$count\t$sum\t" . ($sum / scalar(keys(%data))) . "\n"; '
28109   28109   5.35307560464673


# New version again
mkdir magphase

for j in flye4 clr1 clr2 ; do echo $j; mkdir magphase/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a $j.contigs.fasta -b $j.contigs.fasta.ccs.filt.bam -g $i --bhFDR 0.01 -o magphase/$j/$name.strain"; done; done

mkdir magphase/consolidated
for j in flye4 clr1 clr2; do echo $j; mkdir magphase/consolidated/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/python_toolchain/metagenomics/calcMagPhaseOutputVals.py -f magphase/$j -p $name -d b3c_${j}_dastool/${j}.das_DASTool_summary.txt -o magphase/consolidated/$j/$name; done; done

for i in clr1 clr2 flye4; do echo $i; cat magphase/consolidated/$i/*.short > magphase/consolidated/$i.consolidated.short.tab; done


### Newest version
# averages
for i in clr1 clr2 flye4; do echo $i; cat magphase/consolidated/$i/*.long | perl -e '$c = 0; $m = 0; $l = 0; $max = 0; while(<STDIN>){chomp; @s = split(/\t/); if($s[2] < 4){next;} $c++; $m += $s[6] - $s[5];  $max = ($s[6] - $s[5] > $max)? $s[6] - $s[5] : $max; $l += $s[2];} $m /= $c; $l /= $c; print "$c\t$m\t$max\t$l\n";' ; done
clr1
16139   1008.12714542413        463082  25.4517628105831
clr2
16160   992.585829207921        480257  27.1715965346535
flye4
9919    1100.73727190241        336899  19.6123601169473

for i in clr3; do echo $i; cat clr3_rerun/magphase/consolidated/$i.filt/*.long | perl -e '$c = 0; $m = 0; $l = 0; $max = 0; while(<STDIN>){chomp; @s = split(/\t/); if($s[2] < 4){next;} $c++; $m += $s[6] - $s[5];  $max = ($s[6] - $s[5] > $max)? $s[6] - $s[5] : $max; $l += $s[2];} $m /= $c; $l /= $c; print "$c\t$m\t$max\t$l\n";' ; done
clr3
16313   893.287132961442        493333  26.3890761968982

# Counts
for i in clr1 clr2 flye4; do echo $i; cat magphase/consolidated/$i/*.long | perl -e '%data; $count=0; while(<>){chomp; @F = split(/\t/); if(length($F[1]) > 1){$count++; $t = $F[0] . $F[4] . $F[5]; $data{$t} += 1;}} $sum = 0; $max = 0; foreach $v (values(%data)){$sum += $v; if($v > $max){$max = $v;}} print "$count\t$sum\t$max\t" . ($sum / scalar(keys(%data))) . "\n"; '; done
clr1
19845   19845   59      4.78538702676634
clr2
19954   19954   54      4.61577608142494
flye4
13894   13894   25      4.23210478221139

for i in clr3; do echo $i; cat clr3_rerun/magphase/consolidated/$i.filt/*.long | perl -e '%data; $count=0; while(<>){chomp; @F = split(/\t/); if(length($F[1]) > 1){$count++; $t = $F[0] . $F[4] . $F[5]; $data{$t} += 1;}} $sum = 0; $max = 0; foreach $v (values(%data)){$sum += $v; if($v > $max){$max = $v;}} print "$count\t$sum\t$max\t" . ($sum / scalar(keys(%data))) . "\n"; '; done
clr3
19914   19914   64      4.95619711299154

# New stats in CLR3_RERUN
# NOTE THESE ARE FOR THE NEW DEFITION OF COMPLETE MAGS
for i in clr1 clr2 clr3 flye4; do echo $i; cat magphase/consolidated/$i.consolidated.long.tab | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l  binning/DASTool/$i.contigs.full_DASTool_summary.txt.complete.list -c 0 -d '\t' | perl -e '$c = 0; $m = 0; $l = 0; $max = 0; while(<STDIN>){chomp; @s = split(/\t/); if($s[2] < 4){next;} $c++; $m += $s[6] - $s[5];  $max = ($s[6] - $s[5] > $max)? $s[6] - $s[5] : $max; $l += $s[2];} $m /= $c; $l /= $c; print "$c\t$m\t$max\t$l\n";' ; done
clr1
8516    1106.05648191639        463082  27.0534288398309
clr2
10725   1057.83188811189        480257  24.0620046620047
clr3
16294   961.246164232233        493333  33.9976678531975
flye4
6571    1151.29584538122        336899  20.1491401613149

for i in clr1 clr2 clr3 flye4; do echo $i; cat magphase/consolidated/$i.consolidated.long.tab | python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f stdin -l  binning/DASTool/$i.contigs.full_DASTool_summary.txt.complete.list -c 0 -d '\t' | perl -e '%data; $count=0; while(<>){chomp; @F = split(/\t/); if(length($F[1]) > 1){$count++; $t = $F[0] . $F[4] . $F[5]; $data{$t} += 1;}} $sum = 0; $max = 0; foreach $v (values(%data)){$sum += $v; if($v > $max){$max = $v;}} print "$count\t$sum\t$max\t" . ($sum / scalar(keys(%data))) . "\n"; '; done

```

Now to try to plot the melted data to show the differences between strains

```R
library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggridges)

data <- read.delim("mag_phase_hr/consolidated/flye4_strain_associations.melt", header=TRUE, sep="\t")
data <- data[!data$AssocComp == -1,]

#data$de <- "Nominal"
#data$de[data$StrainDelta >= 2 | data$StrainDelta <= -2] <- "Different"
#data$de <- as.factor(data$de)

pdf(file="mag_phase_hr/consolidated_scg_magphase_strain_delta.pdf", useDingbats=FALSE)
ggplot(data=data, aes(x=StrainDelta, fill=Dataset)) + geom_density(alpha=0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab(label="Strain Delta (<- more in HIFI; -> more in CLR)")
dev.off()

pdf(file="mag_phase_hr/consolidated_scg_magphase_strain_delta.ridges.pdf", useDingbats=FALSE)
ggplot(data=data, aes(x=StrainDelta, y=Dataset, fill=Dataset)) + geom_density_ridges(alpha=0.5) + theme_bw() + theme(text = element_text(size=20)) + xlab(label="Strain Delta (<- more in HIFI; -> more in CLR)")
dev.off()

## TODO: find a better way to visualize this dataset
```

## MagPhase ANI estimates

I can't just pull all SNVs so easily out of MagPhase output, but perhaps I can calculate estimates of ANI% based on the expected length of all genes and the haplotyped SNPs, and I can calculate linked SNPs via the length of haplotypes.

> Ceres:  /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
mkdir mag_phase_hr/ani_estimates;
echo -e "bin\tscglen\tSNPs\tMaxPhase\tANI" > mag_phase_hr/ani_estimates/flye4.scgani.tab; for i in desman/bed_lists/flye4/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; j=mag_phase_hr/consolidated/flye4/$name.long; perl -e 'chomp(@ARGV); $l = 0; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[2] - $s[1] > 0){$l += $s[2] - $s[1];}} close IN; $ls = 0; $ll = 0; $ln = 0; open(IN, "< $ARGV[1]"); %d; while(<IN>){chomp; @s = split(/\t/); if(!exists($d{$s[4]}->{$s[5]}) && $s[6] - $s[5] > 1){$d{$s[4]}->{$s[5]} = 1; $ln += length($s[1]); if(length($s[1]) > $ll){$ll = length($s[1]);}}} close IN; $ls = $ln / $l; $ls = 1 - $ls; print "$ARGV[2]\t$l\t$ln\t$ll\t$ls\n";' $i $j $name >> mag_phase_hr/ani_estimates/flye4.scgani.tab; done

# Now to draw this in with inStrain and then Desman
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; <IN>; while(<IN>){chomp; @s = split(/\t/); $s[17] = ($s[17] eq "")? 0 : $s[17]; $h{$s[0]} = [$s[13],$s[14],$s[17],$s[28]];} close IN; open(IN, "< $ARGV[1]"); <IN>; print "bin\tscglen\tSNPs\tMaxPhase\tANI\tconANI\tpopANI\tLinkSNPs\tdivergents\n"; while($l = <IN>){chomp $l; @s = split(/\t/, $l); if(!exists($h{$s[0]})){next;} $t = join("\t", @{$h{$s[0]}}); print "$l\t$t\n";} close IN;' instrain/subset_flye4/flye4.minmapone/output/flye4.minmapone_genome_info.tsv mag_phase_hr/ani_estimates/flye4.scgani.tab > mag_phase_hr/ani_estimates/flye4.scgani.instrain.tab

# Now strain estimates
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; <IN>; while(<IN>){chomp; @s = split(/\t/); if($s[2] eq "HiFi"){@bsegs = split(/_/, $s[0]); $h{"bin3c.$bsegs[1]"} = $s[1];}} close IN; open(IN, "< $ARGV[1]"); $h = <IN>; chomp($h); $h = $h . "\tDesmanS\n"; print "$h"; while($l = <IN>){chomp $l; @s = split(/\t/, $l); if(!exists($h{$s[0]})){print "$l\tNA\n";}else{print "$l\t" . $h{$s[0]} . "\n";}} close IN;' desman/initial_strain_counts.tab mag_phase_hr/ani_estimates/flye4.scgani.instrain.tab > mag_phase_hr/ani_estimates/flye4.scgani.instrain.desman.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = [$s[1], $s[2], $s[3]];} close IN; open(IN, "< $ARGV[1]"); $h = <IN>; chomp($h); $h = "$h\tMagStrain\tComp\tCont\n"; print "$h"; while($l = <IN>){chomp $l; @s = split(/\t/, $l); if(!exists($h{$s[0]})){next;} $t = join("\t", @{$h{$s[0]}}); print "$l\t$t\n";} close IN;' mag_phase_hr/consolidated/flye4.consolidated.short.tab mag_phase_hr/ani_estimates/flye4.scgani.instrain.desman.tab > mag_phase_hr/ani_estimates/flye4.scgani.instrain.desman.full.tab

# For the CLR datasets
for z in clr1 clr2; do echo $z; echo -e "bin\tscglen\tSNPs\tMaxPhase\tANI" > mag_phase_hr/ani_estimates/$z.scgani.tab; for i in desman/bed_lists/$z/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; j=mag_phase_hr/consolidated/$z/$name.long; perl -e 'chomp(@ARGV); $l = 0; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[2] - $s[1] > 0){$l += $s[2] - $s[1];}} close IN; $ls = 0; $ll = 0; $ln = 0; open(IN, "< $ARGV[1]"); %d; while(<IN>){chomp; @s = split(/\t/); if(!exists($d{$s[4]}->{$s[5]}) && $s[6] - $s[5] > 1){$d{$s[4]}->{$s[5]} = 1; $ln += length($s[1]); if(length($s[1]) > $ll){$ll = length($s[1]);}}} close IN; $ls = $ln / $l; $ls = 1 - $ls; print "$ARGV[2]\t$l\t$ln\t$ll\t$ls\n";' $i $j $name >> mag_phase_hr/ani_estimates/$z.scgani.tab; done; done

for z in clr1 clr2; do echo $z; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = [$s[1], $s[2], $s[3]];} close IN; open(IN, "< $ARGV[1]"); $h = <IN>; chomp($h); $h = "$h\tMagStrain\tComp\tCont\n"; print "$h"; while($l = <IN>){chomp $l; @s = split(/\t/, $l); if(!exists($h{$s[0]})){next;} $t = join("\t", @{$h{$s[0]}}); print "$l\t$t\n";} close IN;' mag_phase_hr/consolidated/$z.consolidated.short.tab mag_phase_hr/ani_estimates/$z.scgani.tab > mag_phase_hr/ani_estimates/$z.scgani.magphaseonly.full.tab; done

for z in clr3; do echo $z; echo -e "bin\tscglen\tSNPs\tMaxPhase\tANI" > mag_phase_hr/ani_estimates/$z.scgani.tab; for i in clr3_rerun/desman/bed_lists/$z/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; j=clr3_rerun/magphase/consolidated/$z/$name.long; perl -e 'chomp(@ARGV); $l = 0; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); if($s[2] - $s[1] > 0){$l += $s[2] - $s[1];}} close IN; $ls = 0; $ll = 0; $ln = 0; open(IN, "< $ARGV[1]"); %d; while(<IN>){chomp; @s = split(/\t/); if(!exists($d{$s[4]}->{$s[5]}) && $s[6] - $s[5] > 1){$d{$s[4]}->{$s[5]} = 1; $ln += length($s[1]); if(length($s[1]) > $ll){$ll = length($s[1]);}}} close IN; $ls = $ln / $l; $ls = 1 - $ls; print "$ARGV[2]\t$l\t$ln\t$ll\t$ls\n";' $i $j $name >> mag_phase_hr/ani_estimates/$z.scgani.tab; done; done

for z in clr3; do echo $z; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = [$s[1], $s[2], $s[3]];} close IN; open(IN, "< $ARGV[1]"); $h = <IN>; chomp($h); $h = "$h\tMagStrain\tComp\tCont\n"; print "$h"; while($l = <IN>){chomp $l; @s = split(/\t/, $l); if(!exists($h{$s[0]})){next;} $t = join("\t", @{$h{$s[0]}}); print "$l\t$t\n";} close IN;' clr3_rerun/magphase/consolidated//$z.consolidated.short.tab mag_phase_hr/ani_estimates/$z.scgani.tab > mag_phase_hr/ani_estimates/$z.scgani.magphaseonly.full.tab; done
```

Comparing some stats now

```R
setwd("C:/SharedFolders/metagenomics/tim_sheep/strain_determination/")
table <- read.delim("flye4.scgani.instrain.desman.full.tab", header=TRUE)
table.dstrains <- table %>% filter(!is.na(DesmanS))

cor(table[,c(2,3,4,5,6,7,8,9,11,12,13)])
                scglen        SNPs    MaxPhase          ANI      conANI
scglen      1.00000000  0.21286291  0.21373817 -0.091186038  0.13300744
SNPs        0.21286291  1.00000000  0.88094913 -0.964819925  0.04482685
MaxPhase    0.21373817  0.88094913  1.00000000 -0.871305609  0.03639086
ANI        -0.09118604 -0.96481992 -0.87130561  1.000000000 -0.02706479
conANI      0.13300744  0.04482685  0.03639086 -0.027064789  1.00000000
popANI      0.08599783  0.08072495  0.07724750 -0.078623754  0.92049504
LinkSNPs    0.04626894  0.39559045  0.36205337 -0.414818581  0.07981850
divergents  0.22667959  0.57029141  0.52884238 -0.557730133  0.05386158
MagStrain   0.11800586  0.52972483  0.56217740 -0.593569546  0.04988282
Comp        0.64663229  0.08828972  0.10022013  0.009963343  0.29808300
Cont        0.47430273  0.11912303  0.16638492 -0.057965589  0.03612167
                popANI    LinkSNPs  divergents   MagStrain         Comp
scglen      0.08599783  0.04626894  0.22667959  0.11800586  0.646632287
SNPs        0.08072495  0.39559045  0.57029141  0.52972483  0.088289716
MaxPhase    0.07724750  0.36205337  0.52884238  0.56217740  0.100220131
ANI        -0.07862375 -0.41481858 -0.55773013 -0.59356955  0.009963343
conANI      0.92049504  0.07981850  0.05386158  0.04988282  0.298082997
popANI      1.00000000  0.06191319  0.09797255  0.10758429  0.174743582
LinkSNPs    0.06191319  1.00000000  0.51223517  0.57680948  0.040239953
divergents  0.09797255  0.51223517  1.00000000  0.60216897  0.110693593
MagStrain   0.10758429  0.57680948  0.60216897  1.00000000 -0.020318304
Comp        0.17474358  0.04023995  0.11069359 -0.02031830  1.000000000
Cont        0.00908554  0.04972587  0.12628234  0.10285933  0.229883113
                  Cont
scglen      0.47430273
SNPs        0.11912303
MaxPhase    0.16638492
ANI        -0.05796559
conANI      0.03612167
popANI      0.00908554
LinkSNPs    0.04972587
divergents  0.12628234
MagStrain   0.10285933
Comp        0.22988311
Cont        1.00000000

# I find it really interesting that the ANI estimates for InStrain are not correlated with the number of SNPs identified.


clr1 <- read.delim("clr1.consolidated.short.tab", header=FALSE)
colnames(clr1) <- c("bin", "MaxHap", "Comp", "Cont")
clr1 <- clr1 %>% mutate(asm="CLR1")

clr2 <- read.delim("clr2.consolidated.short.tab", header=FALSE)
colnames(clr2) <- c("bin", "MaxHap", "Comp", "Cont")
clr2 <- clr2 %>% mutate(asm="CLR2")

clr3 <- read.delim("clr1.consolidated.short.tab", header=FALSE)
colnames(clr3) <- c("bin", "MaxHap", "Comp", "Cont")
clr3 <- clr3 %>% mutate(asm="CLR3")

flye4 <- read.delim("flye4.consolidated.short.tab", header=FALSE)
colnames(flye4) <- c("bin", "MaxHap", "Comp", "Cont")
flye4 <- flye4 %>% mutate(asm="HIFI")

consolidate <- bind_rows(clr1, clr2, clr3, flye4)
ggplot(consolidate, aes(x=MaxHap, fill=asm)) + geom_histogram(aes(y= ..count..), colour="black") + theme_bw() + facet_wrap(~ asm) + xlab("MagPhase Maximum Haplotype Count") + ylab("Number of Bins") + scale_fill_brewer(palette = "Dark2")

# Saved as magphase_strain_count_histo.pdf
```

## MagPhase contig bubble correspondence

I am going to try to take Misha's identified bubbles and correlate them with magphase phased haplotypes. I have to reduce the number of regions compared, but I think that the program can handle it. Misha generated paf files that show correspondence with the edge maps of the contigs so hopefully those will be easy to associate after the pipeline finishes.

```bash
module load miniconda/3.6 samtools bedtools/2.25.0
conda activate /KEEP/rumen_longread_metagenome_assembly/desman
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/sequence/
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/rarefaction/


# Converting the paf files into useable bed files
for i in 1307 1377 2269 5786 68292; do echo $i; perl -lane 'print "$F[5]\t$F[7]\t$F[8]\t$F[0]";' < $i/mapping.paf | bedtools sort -i stdin | bedtools merge -i stdin -delim ';' -c 4 -o 'distinct' -d -1000 > $i/$i.edgemaps.bed; done

mv * magphase-bubbles

# I also should only pull the portions of the beds that represent the smaller bubble regions. 
for i in 1307 1377 2269 5786 68292; do echo $i; perl -lane 'if($F[1] < 350000){print "$F[5]\t$F[7]\t$F[8]\t$F[0]";}' < magphase-bubbles/$i/mapping.paf | bedtools sort -i stdin | bedtools merge -i stdin -delim ';' -c 4 -o 'distinct'  > magphase-bubbles/$i/$i.shortmaps.bed; done

# Now running individual magphase runs on these initial beds 
mkdir magphase-bubbles/phasing
mkdir magphase-bubbles/bubbles

# Non-bubble run is in magphase-bubbles/phasing
for i in magphase-bubbles/*/*.edgemaps.bed; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 10 --mem=300000 -p priority -q msn -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g $i -o magphase-bubbles/phasing/$name.strain"; done

# Bubble run is in magphase-bubbles/bubbles
for i in magphase-bubbles/*/*.shortmaps.bed; do name=`basename $i | cut -d'.' -f1`; echo $name; sbatch -N 1 -n 10 --mem=100000 -p priority -q msn -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a flye4.contigs.fasta -b flye4.contigs.fasta.ccs.bam -g $i -o magphase-bubbles/bubbles/$name.strain"; done


# Consolidating the bubble data
for i in 1307 1377 2269 5786 68292; do echo $i; mkdir magphase-bubbles/bubbles/$i; python3 ~/python_toolchain/metagenomics/calcMagPhaseOutputVals.py -f magphase-bubbles/bubbles -p $i -o magphase-bubbles/bubbles/$i/$i; done

cat magphase-bubbles/bubbles/*/*.short > magphase-bubbles/bubbles_consolidated.short.tab
cat magphase-bubbles/bubbles/*/*.long > magphase-bubbles/bubbles_consolidated.long.tab

# Trying to plot variant positions
grep -P 'contig_1307\t' flye4.contigs.fasta.fai > magphase-bubbles/bubbles/1307.contigs.fasta.fai
python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f magphase-bubbles/bubbles/1307.contigs.fasta.fai -b flye4.contigs.fasta.ccs.bam -u magphase-bubbles/bubbles/1307.strain.human_readable_by_pos.txt -i 15000 -o magphase-bubbles/bubbles/1307.strain.plot -e magphase-bubbles/1307/1307.shortmaps.bed

for i in 1307 1377 2269 68292; do echo $i; grep -P 'contig_'$i'\t' flye4.contigs.fasta.fai > magphase-bubbles/bubbles/$i.contigs.fasta.fai; python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f magphase-bubbles/bubbles/$i.contigs.fasta.fai -b flye4.contigs.fasta.ccs.bam -u magphase-bubbles/bubbles/$i.strain.human_readable_by_pos.txt -i 15000 -o magphase-bubbles/bubbles/$i.strain.plot -e magphase-bubbles/$i/$i.shortmaps.bed; done

# extracting data for Liz to use in her analysis
samtools faidx flye4.contigs.fasta contig_1307 contig_1377 contig_2269 contig_5786 contig_68292 > magphase-bubbles/all_contigs.fasta
samtools faidx magphase-bubbles/all_contigs.fasta

samtools view -t magphase-bubbles/all_contigs.fasta.fai flye4.contigs.fasta.ccs.bam contig_1307 contig_1377 contig_2269 contig_5786 contig_68292 > magphase-bubbles/all_contigs.bubbles.ccs.sam
samtools view -bt magphase-bubbles/all_contigs.fasta.fai magphase-bubbles/all_contigs.bubbles.ccs.sam > magphase-bubbles/all_contigs.bubbles.ccs.bam
```

## InStrain comparison run

I want to compare results from the latest tool, InStrain, against our MagPhase results. I may need to do this on a contig-contig basis because of the phasing, but this should be an interesting comparison.

```bash
module load miniconda/3.6
conda activate /KEEP/rumen_longread_metagenome_assembly/instrain

perl -lane 'print $F[0];' < b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt | sort | uniq > b3c_flye4_dastool/flye4.das_DASTool_scaffolds.list

sbatch -N 1 -n 3 --mem=45000 -p priority -q msn --wrap="prodigal -d b3c_flye4_dastool/flye4.prodigal.fna -o b3c_flye4_dastool/flye4.prodigal.gbk -p meta -i flye4.contigs.fasta -c "

# Trying to speed this up by only listing the contigs that were das-binned
samtools faidx flye4.contigs.fasta -r b3c_flye4_dastool/flye4.das_DASTool_scaffolds.list > flye4.dasbin.only.fasta

sbatch -N 1 -n 3 --mem=45000 -p priority -q msn --wrap="prodigal -d b3c_flye4_dastool/flye4.dasbin.prodigal.fna -o b3c_flye4_dastool/flye4.dasbin.prodigal.gbk -p meta -i flye4.dasbin.only.fasta -c "

# The authors of this program make a big deal about how much memory this takes, so I am going to start with a subset of the short-reads first
# I'm also limiting analysis to the contigs that were part of our bins. Makes sense for future comparisons
sbatch -n 70 -N 1 --mem=300000 -p priority -q msn --wrap="mkdir instrain; mkdir instrain/subset_flye4; inStrain profile -p 70 -o instrain/subset_flye4/flye4.nominmap -g b3c_flye4_dastool/flye4.dasbin.prodigal.fna -s b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt --database_mode --scaffolds_to_profile b3c_flye4_dastool/flye4.das_DASTool_scaffolds.list mapping/flye4.contigs/Lib101_L1.bam flye4.contigs.fasta"

# Now for minquality score filtration
sbatch -n 70 -N 1 --mem=300000 -p priority -q msn --wrap="inStrain profile -p 70 -o instrain/subset_flye4/flye4.minmapone -g b3c_flye4_dastool/flye4.dasbin.prodigal.fna -s b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt --database_mode --min_mapq 1 --scaffolds_to_profile b3c_flye4_dastool/flye4.das_DASTool_scaffolds.list mapping/flye4.contigs/Lib101_L1.bam flye4.contigs.fasta"

# Instrain only calls strains in compare mode, so let's test it out on the duplicate objects

sbatch -n 70 -N 1 --mem=300000 -p priority -q msn --wrap="inStrain compare -p 70 -o instrain/subset_flye4/flye4.minmapone.comp -i instrain/subset_flye4/flye4.minmapone instrain/subset_flye4/flye4.nominmap -s b3c_flye4_dastool/flye4.das_DASTool_scaffolds2bin.txt --database_mode"

# Nope, it has a distinct requirement for more than one sample
# Maybe we can compare estimates for linked SNPs vs estimated ANI to Mag_phase?

# Mini notes for interesting regions to compare via IGV plots
samtools mpileup -r contig_15470:274592-275904 -f flye4.contigs.fasta mapping/flye4.contigs/merged.bam

```

### Plasmid analysis

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye/plasmids

```bash
module load minimap2 samtools
# OK, let's first assess the level of self-hits
sbatch -N 1 -n 3 -p priority -q msn --mem=43000 --wrap="minimap2 -cx asm5 ../flye4.contigs.fasta hifi/hifi_assembly_converted.cycs.fasta > hifi/hifi_assembly.cycs.flye4.map.paf"

# Note: this was a mistaken alignment off of the wrong file in the output
sbatch -N 1 -n 3 -p priority -q msn --mem=43000 --wrap="minimap2 -cx asm20 ../flye4.contigs.fasta hifi/hifi_assembly_converted.cycs.paths.txt > hifi/hifi_assembly.cycs.flye4.asm20.map.paf"

# Checking the distribution of read mappings
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f hifi/hifi_assembly.cycs.flye4.map.paf -d '\t' -c 0 | perl -lane 'if($F[0] eq "Entry"){next;}else{print $F[1];}' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   5479
Sum:    12704
Minimum 1
Maximum 51
Average 2.318671
Median  1
Standard Deviation      4.359142
Mode(Highest Distributed Value) 1

# OK, so I need to filter the data anyways, let's try to filter the data first before alignment
# Contig length
perl -lane '@bsegs = split(/_/, $F[0]); print $bsegs[3];' < hifi/assembly_graph.nodes.scores | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   158066
Sum:    6856596184
Minimum 463
Maximum 5537745
Average 43378.058431
Median  6908
Standard Deviation      165162.885107
Mode(Highest Distributed Value) 503

# Contig coverage
perl -lane '@bsegs = split(/_/, $F[0]); if($bsegs[5] =~ /\.\d{1}.{1}$/){next;}else{print $bsegs[5];}' < hifi/assembly_graph.nodes.scores | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   79033
Sum:    83089054
Minimum 0.0
Maximum 20349898.0
Average 1051.321018
Median  13.0
Standard Deviation      76512.029165
Mode(Highest Distributed Value) 3.0

# And the plasmid scores
perl -lane 'print $F[1];' < hifi/assembly_graph.nodes.scores | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   145680
Sum:    49351.5126359785
Minimum 0.0
Maximum 1.0
Average 0.338767
Median  0.382403994096126
Standard Deviation      0.254138
Mode(Highest Distributed Value) 0.0

# Let's see how many are left if we apply a filter of length 1000000> x > 10 kb, coverage > 50 and plasmid score > 0.5
perl -lane '@bsegs = split(/_/, $F[0]); if($bsegs[5] =~ /\.\d{1}.{1}$/){next;}else{if($bsegs[3] > 10000 && $bsegs[3] < 1000000 && $bsegs[5] > 50.0 && $F[1] >= 0.5){print $F[0];}}' < hifi/assembly_graph.nodes.scores | wc -l
614
perl -lane '@bsegs = split(/_/, $F[0]); if($bsegs[5] =~ /\.\d{1}.{1}$/){next;}else{if($bsegs[3] > 10000 && $bsegs[3] < 1000000 && $bsegs[5] > 50.0 && $F[1] >= 0.5){print $F[0];}}' < hifi/assembly_graph.nodes.scores > hifi/assembly_graph.names.filtered.list

# OK, I misinterpreted the data files. The "scores" are from edges (overlaps between genome and nodes) and the program automatically outputs the supposed plasmids into the fasta file. The fastas are only 5,000 entries, so perhaps I should filter on size and sequence overlaps?
perl -lane 'print $F[1];' < plasmids/hifi/hifi_assembly_converted.cycs.fasta.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   5528
Sum:    267998427
Minimum 1000
Maximum 3587075
Average 48480.178546
Median  5153
Standard Deviation      268859.323068
Mode(Highest Distributed Value) 1573

# It's pretty clear that several whole genomes are getting caught up here because of the algorithm. Let's remove plasmids above 1 mbp and also clear out plasmids under 5kb in size. That would leave 2731 in total
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while($n = <IN>){chomp($n); $n =~ s/>//g; $s = <IN>; if(exists($h{$n})){chomp($s); $s =~ s/(.{1,60})/$1\n/g; print ">$n\n$s\n";}} close IN;' hifi/hifi_assembly_converted.cycs.filtered.list hifi/hifi_assembly_converted.cycs.fasta > hifi/hifi_assembly.filtered.cycs.fasta

# Then, for my algorithm to work, I need to remove fasta entries from the assembly that might overlap with the plasmid sequence
sbatch -N 1 -n 3 -p priority -q msn --mem=43000 --wrap="minimap2 -cx asm20 ../flye4.contigs.fasta hifi/hifi_assembly.filtered.cycs.fasta > hifi/hifi_assembly.filtered.flye4.asm20.map.paf"

# OK, there are some partial overlaps pulled from larger contigs. These could be integration plasmids or misidentified viral sequence. To avoid issues, I will pull only the contigs that map perfectly with the identified plasmid nodes (about 2,191 in total)
perl -lane 'if($F[3] == $F[1] && $F[8] == $F[6] && $F[7] == $F[2]){print "$F[0]\t$F[5]";}' < hifi/hifi_assembly.cycs.flye4.asm20.map.paf > hifi/hifi_assembly.cycs.flye4.map.perfects.tab

# There are 8 nodes that are redundant with 3-2 other contigs each
# printing out the list without the overlap contigs
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[1]} = 0;} foreach $k (keys(%h)){if($h{$k}){print "$k\n"; system("samtools faidx ../flye4.contigs.fasta $k >> $ARGV[2]");}}' ../flye4.contigs.fasta.fai hifi/hifi_assembly.cycs.flye4.map.perfects.tab hifi/flye4.plasmid.contigs.fasta

# Now to add the plasmid nodes
samtools faidx hifi/hifi_assembly.filtered.cycs.fasta
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 0;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = 1;} foreach $k (keys(%h)){if($h{$k}){print "$k\n"; system("samtools faidx hifi/hifi_assembly.filtered.cycs.fasta $k >> $ARGV[2]");}}' hifi/hifi_assembly.filtered.cycs.fasta.fai hifi/hifi_assembly.cycs.flye4.map.perfects.tab hifi/flye4.plasmid.contigs.fasta

# now to align the Hi-C reads
sbatch -N 1 -n 8 -p priority -q msn --mem=45000 --wrap="module load bwa; bwa index hifi/flye4.plasmid.contigs.fasta; bwa mem hifi/flye4.plasmid.contigs.fasta /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R1_001.fastq.gz /project/rumen_longread_metagenome_assembly/sheep_poop/hic_data/Smith_Sheep_63_HC_S2_L001_R2_001.fastq.gz | samtools sort -T temp -o hifi/flye4.plasmid.hiclinks.bam -"

# Now to create the blobtools table that I will use in the script
perl -lane 'print $F[1];' < hifi/hifi_assembly.cycs.flye4.map.perfects.tab > hifi/hifi_assembly.cycs.flye4.map.perfects.ctgs.list

# Reverse grepping a blobtools file
# 29 41
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../blobtools/table.flye4.contigs.blobDB.table.txt -i '#' -d '\t' -c 0 -l hifi/hifi_assembly.cycs.flye4.map.perfects.ctgs.list -v | perl -e 'while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[29]\t$s[41]\n";}' > plasmid_blobtools_file_flye4.tab
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../blobtools/table.flye4.contigs.blobDB.table.txt -i '#' -d '\t' -c 0 -l hifi/hifi_assembly.cycs.flye4.map.perfects.ctgs.list | perl -e 'while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[29]\t$s[41]\n";}' > replacements_blobtools_flye4.tab

# Now to steal the annotation of the analog contigs
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[1]} = $s[0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s =split(/\t/); $s[0] = $h{$s[0]}; print join("\t", @s); print "\n";} close IN;' hifi/hifi_assembly.cycs.flye4.map.perfects.tab replacements_blobtools_flye4.tab > replacements_blobtools_flye4.converted.tab

cat plasmid_blobtools_file_flye4.tab replacements_blobtools_flye4.converted.tab > final_plasmid_blobtools_flye4.tab

# Now to create the viral contig fasta file
python3 ~/python_toolchain/sequenceData/filterFastaqFromList.py -f hifi/flye4.plasmid.contigs.fasta -l hifi/hifi_assembly.cycs.flye4.map.perfects.tab -o plasmid_contigs.fasta

perl -e 'while($h =<>){$s = <>; chomp $s; $s =~ s/(.{1,60})/$1\n/g; print "$h$s";}' < plasmid_contigs.fasta > temp.fasta

mv temp.fasta plasmid_contigs.fasta
samtools faidx plasmid_contigs.fasta

# Now to try running the pipeline!
sbatch -N 1 -n 3 --mem=100000 -p priority -q msn --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a hifi/flye4.plasmid.contigs.fasta -g plasmid_contigs.fasta -b final_plasmid_blobtools_flye4.tab -i hifi/flye4.plasmid.hiclinks.bam -v plasmid_contigs.fasta.fai -m minimap2 -s samtools -o plasmid_contigs.association"

```

#### Create CLR graph bubbles

These are my notes on how to generate bubble graphs for the CLR assemblies using Misha's instructions

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep

```bash
module load miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

cp -r sheep_clr1 sheep_bubbles_clr1
sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.1.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.1.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.1.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.1.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.1.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.1.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.1.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.1.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.1.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.1.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.1.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.1.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.1.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.1.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.1.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.1.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.1.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.1.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.1.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.1.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.1.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.1.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.1.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.1.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.1.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.1.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.1.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.1.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.1.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.1.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.1.fastq.gz -o sheep_bubbles_clr1 --meta -t 70 --keep-haplotypes --resume-from repeat

cp -r sheep_clr2 sheep_bubbles_clr2
sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.2.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.2.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.2.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.2.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.2.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.2.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.2.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.2.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.2.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.2.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.2.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.2.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.2.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.2.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.2.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.2.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.2.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.2.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.2.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.2.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.2.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.2.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.2.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.2.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.2.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.2.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.2.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.2.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.2.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.2.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.2.fastq.gz -o sheep_bubbles_clr2 --meta -t 70 --keep-haplotypes --resume-from repeat

cp -r sheep_clr3 sheep_bubbles_clr3
sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.3.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.3.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.3.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.3.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.3.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.3.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.3.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.3.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.3.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.3.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.3.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.3.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.3.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.3.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.3.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.3.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.3.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.3.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.3.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.3.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.3.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.3.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.3.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.3.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.3.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.3.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.3.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.3.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.3.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.3.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.3.fastq.gz -o sheep_bubbles_clr3 --meta -t 70 --keep-haplotypes --resume-from repeat

mv sheep_bubbles_clr1/assembly_graph.gfa sheep_bubbles_clr1/clr1_bubbles_graph.gfa
mv sheep_bubbles_clr1/assembly_info.txt sheep_bubbles_clr1/clr1_assembly_info.txt
sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gzip sheep_bubbles_clr1/clr1_bubbles_graph.gfa"

mv sheep_bubbles_clr2/assembly_graph.gfa sheep_bubbles_clr2/clr2_bubbles_graph.gfa
mv sheep_bubbles_clr2/assembly_info.txt sheep_bubbles_clr2/clr2_assembly_info.txt
sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gzip sheep_bubbles_clr2/clr2_bubbles_graph.gfa"

mv sheep_bubbles_clr3/assembly_graph.gfa sheep_bubbles_clr3/clr3_bubbles_graph.gfa
mv sheep_bubbles_clr3/assembly_info.txt sheep_bubbles_clr3/clr3_assembly_info.txt
sbatch -N 1 -n 1 --mem=5000 -p priority -q msn --wrap="gzip sheep_bubbles_clr3/clr3_bubbles_graph.gfa"
```

#### Second sheep assembly

Just in case we need it for review:

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep

```bash
module load miniconda/3.6
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-hifi /lustre/project/gaur_genome_assembly/metagenome_assembly/sheep_614208/PacBio/sheep_201614208.CCS.Q20.fastq -o sheep_sample2_metaflye --meta -t 70
```


#### Clostridia example alignments and plots

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye

```bash
module load mauve/2015-02-13 samtools bedtools/2.25.0
progressiveMauve --output=clostridia_aligns/clostridia_mauve.xmfa gtdbtk_bins/flye4.contigs/bin3c_451.fna gtdbtk_bins/flye4.contigs/bin3c_452.fna gtdbtk_bins/flye4.contigs/bin3c_471.fna gtdbtk_bins/clr1.contigs/bin3c_451.fna gtdbtk_bins/clr2.contigs/bin3c_327.fna gtdbtk_bins/clr3.contigs/bin3c_236.fna

# Making average coverage files
for i in gtdbtk_bins/flye4.contigs/bin3c_451.fna gtdbtk_bins/flye4.contigs/bin3c_452.fna gtdbtk_bins/flye4.contigs/bin3c_471.fna; do name=`basename $i | cut -d '.' -f1`; echo $name; samtools faidx $i; cp $i.fai clostridia_aligns/flye4.$name.fai; bedtools makewindows -g $i.fai -w 10000 > clostridia_aligns/flye4.$name.wins; samtools bedcov clostridia_aligns/flye4.$name.wins flye4.contigs.fasta.ccs.bam | perl -lane '$l = $F[2] - $F[1]; $F[3] /= $l; print join("\t", @F);' > clostridia_aligns/flye4.$name.covs.bed; done

j=clr1; i=gtdbtk_bins/clr1.contigs/bin3c_451.fna; name=`basename $i | cut -d '.' -f1`; echo $name; samtools faidx $i;  cp $i.fai clostridia_aligns/$j.$name.fai; bedtools makewindows -g $i.fai -w 10000 > clostridia_aligns/$j.$name.wins; samtools bedcov clostridia_aligns/$j.$name.wins $j.contigs.fasta.ccs.bam | perl -lane '$l = $F[2] - $F[1]; $F[3] /= $l; print join("\t", @F);' > clostridia_aligns/$j.$name.covs.bed;
j=clr2; i=gtdbtk_bins/clr2.contigs/bin3c_327.fna; name=`basename $i | cut -d '.' -f1`; echo $name; samtools faidx $i;  cp $i.fai clostridia_aligns/$j.$name.fai; bedtools makewindows -g $i.fai -w 10000 > clostridia_aligns/$j.$name.wins; samtools bedcov clostridia_aligns/$j.$name.wins $j.contigs.fasta.ccs.bam | perl -lane '$l = $F[2] - $F[1]; $F[3] /= $l; print join("\t", @F);' > clostridia_aligns/$j.$name.covs.bed;
j=clr3; i=gtdbtk_bins/clr3.contigs/bin3c_236.fna; name=`basename $i | cut -d '.' -f1`; echo $name; samtools faidx $i;  cp $i.fai clostridia_aligns/$j.$name.fai; bedtools makewindows -g $i.fai -w 10000 > clostridia_aligns/$j.$name.wins; samtools bedcov clostridia_aligns/$j.$name.wins $j.contigs.fasta.ccs.bam | perl -lane '$l = $F[2] - $F[1]; $F[3] /= $l; print join("\t", @F);' > clostridia_aligns/$j.$name.covs.bed;

# Let's prepare some files for Liz to create IGV plots from
mkdir clostridia_aligns/bams
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b flye4.contigs.fasta.ccs.bam contig_15470 contig_23521 contig_23530 contig_23535 contig_47300 contig_54438 contig_55624 >  clostridia_aligns/bams/hifi.451.ccs.bam"
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b flye4.contigs.fasta.ccs.bam contig_23527 contig_23533 contig_23537 contig_4105 >  clostridia_aligns/bams/hifi.452.ccs.bam "
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b flye4.contigs.fasta.ccs.bam contig_23522 contig_23534 contig_23536 contig_23538 contig_4092 >  clostridia_aligns/bams/hifi.471.ccs.bam"

sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b mapping/flye4.contigs/merged.bam contig_15470 contig_23521 contig_23530 contig_23535 contig_47300 contig_54438 contig_55624 >  clostridia_aligns/bams/hifi.451.sr.bam"
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b mapping/flye4.contigs/merged.bam contig_23527 contig_23533 contig_23537 contig_4105 >  clostridia_aligns/bams/hifi.452.sr.bam"
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b mapping/flye4.contigs/merged.bam contig_23522 contig_23534 contig_23536 contig_23538 contig_4092 >  clostridia_aligns/bams/hifi.471.sr.bam"

cp desman/bed_lists/flye4/bin3c.451.scg.bed clostridia_aligns/bams/hifi.451.scg.bed
cp desman/bed_lists/flye4/bin3c.452.scg.bed clostridia_aligns/bams/hifi.452.scg.bed
cp desman/bed_lists/flye4/bin3c.471.scg.bed clostridia_aligns/bams/hifi.471.scg.bed

cp b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.451.contigs.fa clostridia_aligns/bams/hifi.451.contigs.fa
cp b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.452.contigs.fa clostridia_aligns/bams/hifi.452.contigs.fa
cp b3c_flye4_dastool/flye4.das_DASTool_bins/bin3c.471.contigs.fa clostridia_aligns/bams/hifi.471.contigs.fa

for i in clr1 clr2 clr3; do echo $i; 
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b minimap_clr_bams/flye4.$i.sorted.bam contig_15470 contig_23521 contig_23530 contig_23535 contig_47300 contig_54438 contig_55624 >  clostridia_aligns/bams/hifi.451.$i.bam"
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b minimap_clr_bams/flye4.$i.sorted.bam contig_23527 contig_23533 contig_23537 contig_4105 >  clostridia_aligns/bams/hifi.452.$i.bam "
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b minimap_clr_bams/flye4.$i.sorted.bam contig_23522 contig_23534 contig_23536 contig_23538 contig_4092 >  clostridia_aligns/bams/hifi.471.$i.bam"; done

for i in 451 452 471; do echo $i; sbatch -N 1 -n 1 --mem=8000 -p priority -q msn --wrap="samtools merge clostridia_aligns/bams/hifi.$i.clr.bam clostridia_aligns/bams/hifi.$i.clr1.bam clostridia_aligns/bams/hifi.$i.clr2.bam clostridia_aligns/bams/hifi.$i.clr3.bam"; done

# And let's repeat the whole shebang for CLR1 bin 451
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b clr1.contigs.fasta.ccs.bam contig_13955 contig_18132 contig_24073 contig_25426 contig_25429 contig_25431 contig_27381 contig_27384 contig_27387 contig_28324 contig_28933 contig_28939 contig_28969 contig_35732 contig_37225 contig_37668 contig_43105 contig_43107 contig_43111 contig_43291 contig_43293 contig_45819 contig_46371 contig_49116 contig_62537 contig_62539 contig_62540 >  clostridia_aligns/bams/clr1.451.ccs.bam"

sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b mapping/clr1.contigs/merged.bam contig_13955 contig_18132 contig_24073 contig_25426 contig_25429 contig_25431 contig_27381 contig_27384 contig_27387 contig_28324 contig_28933 contig_28939 contig_28969 contig_35732 contig_37225 contig_37668 contig_43105 contig_43107 contig_43111 contig_43291 contig_43293 contig_45819 contig_46371 contig_49116 contig_62537 contig_62539 contig_62540 >  clostridia_aligns/bams/clr1.451.sr.bam"

cp desman/bed_lists/clr1/bin3c.451.scg.bed clostridia_aligns/bams/clr1.451.scg.bed

cp b3c_clr1_dastool/clr1.das_DASTool_bins/bin3c.451.contigs.fa clostridia_aligns/bams/clr1.451.contigs.fa

for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b minimap_clr_bams/clr1.$i.sorted.bam contig_13955 contig_18132 contig_24073 contig_25426 contig_25429 contig_25431 contig_27381 contig_27384 contig_27387 contig_28324 contig_28933 contig_28939 contig_28969 contig_35732 contig_37225 contig_37668 contig_43105 contig_43107 contig_43111 contig_43291 contig_43293 contig_45819 contig_46371 contig_49116 contig_62537 contig_62539 contig_62540 >  clostridia_aligns/bams/clr1.451.$i.bam"; done

sbatch -N 1 -n 1 --mem=8000 -p priority -q msn --wrap="samtools merge clostridia_aligns/bams/clr1.451.clr.bam clostridia_aligns/bams/clr1.451.clr1.bam clostridia_aligns/bams/clr1.451.clr2.bam clostridia_aligns/bams/clr1.451.clr3.bam"

```

#### Testing magphase hap edit distance

```bash
module load clustalomega/1.2.4
mkdir mag_phase_hr/hapaligns; for i in mag_phase_hr/flye4/*.strain.human_readable_by_hap.txt; do name=`basename $i | cut -d'.' -f1,2`; echo $name; perl -e '@rows; chomp(@ARGV); open(IN, "< $ARGV[1]"); while($l = <IN>){chomp $l; @s = split(/\t/, $l); if($s[1] == 0){ if(scalar(@rows) > 10){open(OUT, "> mag_phase_hr/hapaligns/$ARGV[0]." . $rows[0]->[2] . ".fasta"); foreach $r (@rows){ print {OUT} ">" . $r->[2] . "." . $r->[1] . "." . $r->[3] . "\n" . $r->[0] . "\n";} close OUT;} @rows = (); } if($s[0] =~ /\?/){next;}else{push(@rows, [$s[0], $s[1], $s[2], $s[3]]);}} if(scalar(@rows) > 10){open(OUT, "> mag_phase_hr/hapaligns/$ARGV[0]." . $rows[0]->[2] . ".fasta"); foreach $r (@rows){ print {OUT} ">" . $r->[2] . "." . $r->[1] . "." . $r->[3] . "\n" . $r->[0] . "\n";} close OUT;}' $name $i; done

for i in mag_phase_hr/hapaligns/*.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; clustalo -i $i -o $i.align; done

# Getting the lengths:
for i in mag_phase_hr/hapaligns/*.fasta; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $h = <IN>; $s = <IN>; chomp($s); print "$ARGV[0]\t" . length($s) . "\n";' $i; done

# Full distribution of values where hap counts > 10
for i in mag_phase_hr/hapaligns/*.faperl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $h = <IN>; $s = <IN>; chomp($s); print length($s) . "\n";' $i; done | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   48
Sum:    1577
Minimum 4
Maximum 168
Average 32.854167
Median  12.5
Standard Deviation      39.837964
Mode(Highest Distributed Value) 4

# values that are 4
for i in mag_phase_hr/hapaligns/*.fasta; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $h = <IN>; $s = <IN>; chomp($s); if(length($s) < 5){print length($s) . "\n";}' $i; done | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   14
Sum:    56
Minimum 4
Maximum 4
Average 4.000000
Median  4
Standard Deviation      0.000000
Mode(Highest Distributed Value) 4

# Calc hamming distances between adjacent elements
for i in mag_phase_hr/hapaligns/*.fasta; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $h = <IN>; $s = <IN>; chomp($s); $last = ""; if(length($s) > 10){$last = $s; while($h = <IN>){ @asegs = split(//, $last); $s = <IN>; chomp($s); @bsegs = split(//, $s); $count = 0; for($x = 0; $x < scalar(@asegs); $x++){if($asegs[$x] ne $bsegs[$x]){$count++;}} chomp($h); print "$ARGV[0]\t$h\t$count\n"; $last = $s;} }' $i; done

for i in mag_phase_hr/hapaligns/*.faperl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $h = <IN>; $s = <IN>; chomp($s); $last = ""; if(length($s) > 10){$last = $s; while($h = <IN>){ @asegs = split(//, $last); $s = <IN>; chomp($s); @bsegs = split(//, $s); $count = 0; for($x = 0; $x < scalar(@asegs); $x++){if($asegs[$x] ne $bsegs[$x]){$count++;}} chomp($h); print "$ARGV[0]\t$h\t$count\n"; $last = $s;} }' $i; done | perl -lane 'print "$F[2]";' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   305			<- 122 are less than 4 bases different
Sum:    5681
Minimum 1
Maximum 131
Average 18.626230
Median  8
Standard Deviation      24.262096
Mode(Highest Distributed Value) 2

# Creating a table
echo -e "file\tcontig\tHapNum\tCoverage\tHapSeq\tcontig\tHapNum\tCoverage\tHapSeq\tDistance" > mag_phase_hr/hapaligns/haps_gt10_distances.tab; for i in mag_phase_hr/hapaligns/*.fasta; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $fh = <IN>; chomp($fh); $fh =~ s/>//g; $s = <IN>; chomp($s); $last = ""; if(length($s) > 10){$last = $s; while($h = <IN>){ @asegs = split(//, $last); $s = <IN>; chomp($s); @bsegs = split(//, $s); $count = 0; for($x = 0; $x < scalar(@asegs); $x++){if($asegs[$x] ne $bsegs[$x]){$count++;}} chomp($h); $h =~ s/>//g; @fsegs = split(/\//, $ARGV[0]); @fhsegs = split(/\./, $fh); @nhsegs = split(/\./, $h); print "$fsegs[-1]\t" . join("\t", @fhsegs) . "\t$last\t" . join("\t", @nhsegs) . "\t$s\t$count\n";} }' $i; done >> mag_phase_hr/hapaligns/haps_gt10_distances.tab
```

## intra-contig Mash distances and mapq0 wins

We think that better resolution of bins could be accomplished if we can compare Mash distances, and contig read mapping can be better estimated by counting Mapq0 bins. Let's try that.

```bash
mkdir mash_intra_bin

for j in flye4 clr1 clr2 clr3; do echo $j; sbatch -N 1 -n 2 --mem=10000 -p priority -q msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist gtdbtk_bins/${j}_k21_s100000_combined.msh gtdbtk_bins/${j}_k21_s100000_combined.msh > mash_intra_bin/${j}.intrabin.pairwise"; done

# Fixing with new CLR3 assembly
sbatch -N 1 -n 2 --mem=10000 -p priority -q msn --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist clr3_rerun/gtdbtk_bins/clr3_k21_s100000_combined.msh clr3_rerun/gtdbtk_bins/clr3_k21_s100000_combined.msh > mash_intra_bin/clr3.intrabin.pairwise"

# I wrote a script to pull out only the nearest neighbors
python3 mash_intra_bin/create_table.py mash_intra_bin/flye4.intrabin.pairwise HIFI gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv mash_intra_bin/flye4.intrabin.neighbors.tab
for i in 1 2; do echo $i; python3 mash_intra_bin/create_table.py mash_intra_bin/clr${i}.intrabin.pairwise CLR${i} gtdbtk_output/clr${i}.contigs/gtdbtk.bac120.summary.tsv mash_intra_bin/clr${i}.intrabin.neighbors.tab; done
for i in 3; do echo $i; python3 mash_intra_bin/create_table.py mash_intra_bin/clr${i}.intrabin.pairwise CLR${i} clr3_rerun/gtdbtk_output/clr3.contigs/gtdbtk.bac120.summary.tsv mash_intra_bin/clr${i}.intrabin.neighbors.tab; done


for i in flye4 clr1 clr2 clr3; do echo $i; perl -lane 'if($F[0] eq "Bin"){next;}else{print "$F[2]";}' < mash_intra_bin/${i}.intrabin.neighbors.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl ; done
flye4
total   654
Minimum 0.0277154
Maximum 0.389565
Average 0.198642
Median  0.199592

clr1
total   525
Minimum 0.0224163
Maximum 0.396904
Average 0.214028
Median  0.212502

clr2
total   518
Minimum 0.0390759
Maximum 0.393093
Average 0.212392
Median  0.2101145

clr3
total   241
Minimum 0.00398972
Maximum 0.410602
Average 0.235563
Median  0.232418

### And the Mapq0 bins
mkdir mapq_checks
sbatch -N 1 -n 2 --mem=30000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/getBAMMapQ0Ratios.py -b mapping/flye4.contigs/merged.bam -o mapq_checks/flye4.shortreads.mapq.bed"

sbatch -N 1 -n 2 --mem=30000 -p priority -q msn --wrap="python3 ~/python_toolchain/sequenceData/getBAMMapQ0Ratios.py -b flye4.contigs.fasta.ccs.bam -o mapq_checks/flye4.hifireads.mapq.bed"

file=mapq_checks/flye4.shortreads.mapq; echo -e "TECH\tRATIO\tZERO" > $file.tab;  perl -lane '$tech = "SHORT"; $ratio = ($F[-2] != 0)? $F[-1] / $F[-2] : 0; $zero = ($F[-2] != 0)? "FALSE" : "TRUE"; print "$tech\t$ratio\t$zero";' < $file.bed >> $file.tab;
file=mapq_checks/flye4.hifireads.mapq; echo -e "TECH\tRATIO\tZERO" > $file.tab;  perl -lane '$tech = "HIFI"; $ratio = ($F[-2] != 0)? $F[-1] / $F[-2] : 0; $zero = ($F[-2] != 0)? "FALSE" : "TRUE"; print "$tech\t$ratio\t$zero";' < $file.bed >> $file.tab;
```

```R
setwd("C:/SharedFolders/metagenomics/tim_sheep/mapq_checks/")
library(dplyr)
library(ggplot2)

temp <- read.delim("flye4.hifireads.mapq.tab", header = TRUE)
temp2 <- read.delim("flye4.shortreads.mapq.tab", header = TRUE)
combined <- bind_rows(temp, temp2)
combined$TECH <- as.factor(combined$TECH)

ggplot(data=combined, aes(x=TECH, y=RATIO, fill=TECH)) + geom_violin() + scale_fill_brewer(palette="Dark2") + theme_bw()
# File was saved as hifi_assembly_short_vs_hifi_mapq0.pdf
```

## Rerunning CLR3

Misha found out that CLR3 looks kinda odd and should not have the issues we identified previously. He suggested that I restart from the consensus stage to see if that resolves the assembly. 

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/

```bash
# First, backing up data
sbatch -N 1 -n 1 -p priority -q msn --mem=5000 --wrap="cp -r sheep_clr3 sheep_clr3_backup"

module load miniconda/3.6 minimap2
source activate /KEEP/rumen_longread_metagenome_assembly/flye

sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.3.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.3.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.3.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.3.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.3.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.3.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.3.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.3.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.3.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.3.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.3.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.3.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.3.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.3.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.3.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.3.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.3.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.3.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.3.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.3.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.3.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.3.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.3.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.3.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.3.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.3.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.3.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.3.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.3.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.3.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.3.fastq.gz -o sheep_clr3 --meta -t 70 --resume-from consensus

cp sheep_clr3/assembly.fasta sheep_clr3/asm1/sheep_clr3.fasta

conda activate /KEEP/rumen_longread_metagenome_assembly/python3/
cd sheep_clr3/

sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/pilonCorrection/Snakefile --cluster "sbatch --nodes=1 --ntasks-per-node={threads} --mem=33000 -p priority -q msn -o logs/{rule}.stdout" --latency-wait 40 --jobs 50

for i in pilon_round1/*.fasta; do echo $i; perl -ne 'if($_ =~ /^>/){chomp; $_ =~ s/_pilon//; $_ .= "\n";} print $_;' < $i >> sheep_clr3_pilon1.fasta; done
sbatch -N 1 -n 2 --mem=15000 -p priority -q msn --wrap="bwa index sheep_clr3_pilon1.fasta"
sbatch -N 1 -n 70 --mem=300000 --dependency=afterok:5653656 -p priority -q msn --wrap="bwa mem -t 34 sheep_clr3_pilon1.fasta /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz | samtools sort -T sheep_pilon2 -@ 34 -o sheep_clr3_pilon1.bam -"

sbatch -N 1 -n 2 --mem=15000 -p priority -q msn --dependency=afterok:5653661 --wrap="samtools index sheep_clr3_pilon1.bam; python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f sheep_clr3_pilon1.bam -g sheep_clr3_pilon1.fasta -o sheep_clr3_pilon2 -p priority -q msn -e 30000 -m"


# And the bubbles..
conda activate /KEEP/rumen_longread_metagenome_assembly/flye/
cp -r sheep_clr3 sheep_bubbles_clr3
sbatch --nodes=1 --mem=1000000 --ntasks-per-node=70 -p priority-mem -q msn-mem flye -g 1.0g --pacbio-raw ../../sheep_poop/m54033_180919_161442_clr.3.fastq.gz ../../sheep_poop/m54033_180921_182853_clr.3.fastq.gz ../../sheep_poop/m54033_181128_171116_clr.3.fastq.gz ../../sheep_poop/m54033_181130_213801_clr.3.fastq.gz ../../sheep_poop/m54033_181203_205641_clr.3.fastq.gz ../../sheep_poop/m54033_181204_171640_clr.3.fastq.gz ../../sheep_poop/m54033_181217_171601_clr.3.fastq.gz ../../sheep_poop/m54033_181219_095558_clr.3.fastq.gz ../../sheep_poop/m54033_181220_061614_clr.3.fastq.gz ../../sheep_poop/m54033_190206_141429_clr.3.fastq.gz ../../sheep_poop/m54033_190207_103359_clr.3.fastq.gz ../../sheep_poop/m54337_181123_182744_clr.3.fastq.gz ../../sheep_poop/m54337_181124_144902_clr.3.fastq.gz ../../sheep_poop/m54337_181125_110904_clr.3.fastq.gz ../../sheep_poop/m54337_181126_072918_clr.3.fastq.gz ../../sheep_poop/m54337_181210_204959_clr.3.fastq.gz ../../sheep_poop/m54337_181212_153618_clr.3.fastq.gz ../../sheep_poop/m54337_181214_195839_clr.3.fastq.gz ../../sheep_poop/m54337_181215_161949_clr.3.fastq.gz ../../sheep_poop/m54337_181217_200214_clr.3.fastq.gz ../../sheep_poop/m54337_181218_162309_clr.3.fastq.gz ../../sheep_poop/m54337_181219_124328_clr.3.fastq.gz ../../sheep_poop/m54337_181220_090342_clr.3.fastq.gz ../../sheep_poop/m54337U_200203_184822_clr.3.fastq.gz ../../sheep_poop/m54337U_200211_145859_clr.3.fastq.gz ../../sheep_poop/m54337U_200213_024013_clr.3.fastq.gz ../../sheep_poop/m54337U_200214_212611_clr.3.fastq.gz ../../sheep_poop/m54337U_200220_195233_clr.3.fastq.gz ../../sheep_poop/m54337U_200222_075656_clr.3.fastq.gz ../../sheep_poop/m54337U_200223_200916_clr.3.fastq.gz ../../sheep_poop/m54337U_200227_201603_clr.3.fastq.gz -o sheep_bubbles_clr3 --meta -t 70 --keep-haplotypes --resume-from repeat
```

## Rerunning on fixed CLR3 assembly

This is going to be a disaster, but here goes nothing!

> Ceres: /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun

```bash
cp -r ../assembly/ ./assembly
cp -r ../mapping/ ./mapping
cp -r ../binning ./binning

# Note: only works with module load python_3/3.6.6
sbatch -N 1 -n 2 --mem=10000 -p priority -q msn snakemake -s ~/python_toolchain/snakeMake/hifiMAGManuscript/Snakefile --cluster-config ~/python_toolchain/snakeMake/hifiMAGManuscript/cluster.json --cluster "sbatch -N 1 --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} -p priority -q msn -o {cluster.stdout}" -p --use-conda --jobs 250 --verbose --latency-wait 40

for i in clr1.contigs clr2.contigs clr3.contigs flye4.contigs; do echo $i; mkdir assembly/$i; perl -e '%data; chomp(@ARGV); open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = $s[1];} close IN; $count = 0; foreach $k (sort {$data{$b} <=> $data{$a}} keys(%data)){system("samtools faidx assembly/$ARGV[1].fa $k > assembly/$ARGV[1]/$k.fasta"); $count++; if($count >= 1000){last;}}' assembly/$i.fa.fai $i; done

for i in clr1 clr2 clr3 flye4; do echo $i; sbatch -N 1 -n 30 --mem=150000 -p priority -q msn --wrap="checkm lineage_wf -f tables/$i.contigs.contigs.checkm.txt -t 30 -x fasta {output.directory}"; done

### NOTE: Redo ridgeline plots and distributions after removing bins with > 5% contamination


### GTDB-TK recalculation
mkdir gtdbtk_bins
for i in clr1 clr2 clr3 flye4; do echo $i; mkdir gtdbtk_bins/${i}.contigs; sbatch ../create_bins.pl binning/DASTool/${i}.contigs.full_bin3c.eval binning/DASTool/${i}.contigs.full_cluster_attribution.tsv assembly/${i}.contigs.fa gtdbtk_bins/${i}.contigs; done

for i in gtdbtk_bins/clr3.contigs/*.fna; do echo -n -e "$i\t"; samtools faidx $i; perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < $i.fai; done > gtdbtk_bins/clr3.contigs.binsizes.tab
perl -lane '@bsegs = split(/[\/\.]/, $F[0]); print "$bsegs[-2]\t$F[1]";' < gtdbtk_bins/clr3.contigs.binsizes.tab > gtdbtk_bins/clr3.contigs.binsizes.fmt.tab


conda activate /KEEP/rumen_longread_metagenome_assembly/gtdbtk

mkdir gtdbtk_output
for i in clr3; do mkdir gtdbtk_output/${i}.contigs; sbatch -N 1 -n 30 --mem=700000 -p priority-mem -q msn-mem --wrap="gtdbtk classify_wf --genome_dir gtdbtk_bins/${i}.contigs --out_dir gtdbtk_output/${i}.contigs --cpus 30 --pplacer_cpus 1"; done

cat gtdbtk_output/clr3.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/clr3.contigs/gtdbtk.ar122.summary.tsv | perl -ne 'chomp; @F = split(/\t/); if($F[0] =~ /^user_genome/){next;} print "$F[0]\t$F[1]\n";' > gtdbtk_output/clr3.contigs/clr3.concat_taxonomy.tab

for i in clr3; do echo $i; perl -ne 'chomp; @F = split(/\t/); $F[1] =~ s/\s/_/g; print "$F[1]\n";' < gtdbtk_output/${i}.contigs/${i}.concat_taxonomy.tab > gtdbtk_output/${i}.contigs/${i}.tax.list; done

### MASH recalculation
for i in clr3; do echo $i; sbatch ../mash_sketch_bins.pl gtdbtk_bins/${i}.contigs $i gtdbtk_bins/${i}_k21_s100000_combined; done

cp ../gtdbtk_bins/flye4_k21_s100000_combined.msh ./gtdbtk_bins/
for i in clr3; do echo $i; sbatch -N 1 -n 2 --mem=35000 -p priority -q msn --dependency=afterok:5661415 --wrap="~/rumen_longread_metagenome_assembly/binaries/mash-Linux64-v2.0/mash dist -v 0.1 gtdbtk_bins/flye4_k21_s100000_combined.msh gtdbtk_bins/${i}_k21_s100000_combined.msh > gtdbtk_bins/flye4_${i}_distance.tab"; done

for i in clr1 clr2; do echo $i; cp ../gtdbtk_bins/flye4_${i}_* ./gtdbtk_bins/; cp ../gtdbtk_bins/${i}_*.msh ./gtdbtk_bins/; done

for i in clr3; do echo $i; perl -MFile::Basename -lane 'if($F[2] < 0.1){$s = basename($F[0]); $r = basename($F[1]); $s =~ s/\..+$//; $r =~ s/\..+$//; print "$s\t$r\t$F[2]";}' < gtdbtk_bins/flye4_${i}_distance.tab > gtdbtk_bins/flye4_${i}_filtreformat_distance.tab; done

for i in `ls gtdbtk_bins/clr3.contigs/*.fna`; do name=`basename $i`; echo $name; done > gtdbtk_bins/clr3.scaffolds.list

conda activate /KEEP/rumen_longread_metagenome_assembly/seaborn
sbatch -N 1 -n 2 -p priority -q msn --mem=45000 --wrap="python3 ~/python_toolchain/metagenomics/mashSparseToSortMatrix.py -r ../b3c_flye4_dastool/flye4.scaffolds.list -q gtdbtk_bins/clr3.scaffolds.list -d gtdbtk_bins/flye4_clr3_distance.tab -o gtdbtk_bins/flye4_clr3_distance.matrix"


for i in assembly/clr3.contigs.fa; do echo $i; sbatch -N 1 -n 70 -p priority -q msn --mem=155000 -t 3-0 --wrap="minimap2 -ax asm20 -t 35 -R '@RG\tID:CCS\tSM:CCS' $i /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequel1_CCS.fasta.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta.gz | samtools sort -@ 35 -T $i.tmp -o $i.ccs.bam -"; done

for i in clr3; do echo $i; perl -e 'use File::Basename; while(<>){chomp; @s = split(/\t/); $r = basename($s[0]); $t = basename($s[1]); $r =~ s/\.fna//; $t =~ s/\.fna//; if($s[2] < 0.1){print "$r\t$t\t$s[2]\n";}}' < gtdbtk_bins/flye4_${i}_distance.tab > gtdbtk_bins/flye4_${i}_associations.tab; done

for i in clr3; do echo $i; perl combine_bin_tax_and_pairing.pl gtdbtk_output/flye4.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/flye4.contigs/gtdbtk.ar122.summary.tsv gtdbtk_output/${i}.contigs/gtdbtk.bac120.summary.tsv gtdbtk_output/${i}.contigs/gtdbtk.ar122.summary.tsv gtdbtk_bins/flye4_${i}_associations.tab gtdbtk_output/flye4_${i}_assoc_tax.tab; done

### Magphase
# making the scg beds
cat binning/DASTool/clr3.contigs.full_proteins.faa.archaea.scg binning/DASTool/clr3.contigs.full_proteins.faa.bacteria.scg > binning/DASTool/clr3.contigs.full_proteins.faa.combined.scg; cat binning/DASTool/clr3.contigs.full_proteins.faa | grep '>' | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' >  desman/clr3.das_prodigal.shortform.tab; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$r = $s[0]; $r =~ s/_\d{1,3}//; print "$h{$s[0]},$r,$s[1],$s[2],$s[0],$s[3]\n";}} close IN;' binning/DASTool/clr3.contigs.full_proteins.faa.combined.scg desman/clr3.das_prodigal.shortform.tab > desman/clr3.das_prodigal.master_cogs.csv;

for i in clr3; do echo $i; cat binning/DASTool/$i.contigs.full_proteins.faa.combined.scg | perl -lane 'print $F[0];' > desman/$i.scg.list; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f desman/$i.das_prodigal.shortform.tab -c 0 -l desman/$i.scg.list | perl -lane '$r = $F[0]; $r  =~ s/_\d{1,3}$//; print "$r\t$F[1]\t$F[2]\t$F[0]";' > desman/$i.scg.loc.bed; perl -lane 'open(OUT, ">> desman/bin_lists/$F[1].hqdas.bin.list"); print OUT "$F[0]"; close OUT;' < binning/DASTool/$i.contigs.full_cluster_attribution.tsv; mv desman/bin_lists/*.list desman/bin_lists/$i/; for j in desman/bin_lists/$i/*.list; do name=`basename $j | cut -d'.' -f1,2`; echo " $name"; python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f desman/$i.scg.loc.bed -c 0 -l $j > desman/bed_lists/$i/${name}.scg.bed; done; done

mkdir magphase
conda activate /KEEP/rumen_longread_metagenome_assembly/desman
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/sequence/
export PATH=$PATH:/lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/rarefaction/

for j in clr3; do echo $j; mkdir magphase/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a assembly/$j.contigs.fa -b assembly/$j.contigs.fa.ccs.check.bam -g $i --bhFDR 0.01 -o magphase/$j/$name.strain"; done; done

sbatch -N 1 -n 1 --mem=45000 -p priority -q msn --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/sequence/filter_bam_by_coverage.py -c 0.9 --filter_secondary --filter_supp assembly/clr3.contigs.fa.ccs.check.bam assembly/clr3.contigs.fa.ccs.filt.bam"
samtools index assembly/clr3.contigs.fa.ccs.filt.bam

for j in clr3; do echo $j; mkdir magphase/$j.filt; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; sbatch -N 1 -n 2 -p priority -q msn --mem=45000 -t 1-0 --wrap="python /lustre/project/rumen_longread_metagenome_assembly/binaries/cDNA_Cupcake/phasing/mag_phaser.py -a assembly/$j.contigs.fa -b assembly/$j.contigs.fa.ccs.filt.bam -g $i --bhFDR 0.01 -o magphase/$j.filt/$name.strain"; done; done

mkdir magphase/consolidated
for j in clr3; do echo $j; mkdir magphase/consolidated/$j; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/python_toolchain/metagenomics/calcMagPhaseOutputVals.py -f magphase/$j -p $name -d binning/DASTool/clr3.contigs.full_DASTool_summary.txt -o magphase/consolidated/$j/$name; done; done

for j in clr3; do echo $j; mkdir magphase/consolidated/$j.filt; for i in desman/bed_lists/$j/*.scg.bed; do name=`basename $i | cut -d'.' -f1,2`; echo $name; python3 ~/python_toolchain/metagenomics/calcMagPhaseOutputVals.py -f magphase/$j.filt -p $name -d binning/DASTool/clr3.contigs.full_DASTool_summary.txt -o magphase/consolidated/$j.filt/$name; done; done

cat magphase/consolidated/clr3/*.short > magphase/consolidated/clr3.consolidated.short.tab
cat magphase/consolidated/clr3.filt/*.short > magphase/consolidated/clr3.consolidated.filt.short.tab

python3 ~/python_toolchain/metagenomics/createMAGstrainAssoc.py -r ../mag_phase_hr/consolidated/flye4.consolidated.short.tab -o magphase/consolidated/flye4_strain_associations -n clr1 -a ../gtdbtk_bins/flye4_clr1_associations.tab -s ../mag_phase_hr/consolidated/clr1.consolidated.short.tab -n clr2 -a ../gtdbtk_bins/flye4_clr2_associations.tab -s ../mag_phase_hr/consolidated/clr2.consolidated.short.tab -n clr3 -a gtdbtk_bins/flye4_clr3_associations.tab -s magphase/consolidated/clr3.consolidated.short.tab

python3 ~/python_toolchain/metagenomics/createMAGstrainAssoc.py -r ../magphase/consolidated/flye4.consolidated.short.tab -o magphase/consolidated/flye4_strain_associations_filt -n clr1 -a ../gtdbtk_bins/flye4_clr1_associations.tab -s ../magphase/consolidated/clr1.consolidated.short.tab -n clr2 -a ../gtdbtk_bins/flye4_clr2_associations.tab -s ../magphase/consolidated/clr2.consolidated.short.tab -n clr3 -a gtdbtk_bins/flye4_clr3_associations.tab -s magphase/consolidated/clr3.consolidated.filt.short.tab

### TODO: copy over the new output file and update the table in the manuscript.

# NOTE: example clostridia clr3 bin is bin3c_367 now with mash distance 0.06 to all three flye4 bins
mkdir example_clostridia
conda activate /KEEP/rumen_longread_metagenome_assembly/desman
j=clr3; i=367; samtools faidx gtdbtk_bins/clr3.contigs/bin3c_367.fna; sbatch -N 1 -n 2 --mem=75000 -p priority -q msn --wrap="python3 ~/python_toolchain/metagenomics/plotMagPhaseOutput.py -f gtdbtk_bins/clr3.contigs/bin3c_367.fna.fai -b assembly/clr3.contigs.fa.ccs.check.bam -u magphase/$j/bin3c.$i.strain.human_readable_by_pos.txt -i 1000 -o magphase/$j.bin3c.367.nolines"

grep 'bin3c.367' binning/DASTool/clr3.contigs.full_cluster_attribution.tsv | perl -ne 'chomp; @s = split(/\t/); print "$s[0] ";' ; echo
contig_33698 contig_37094 contig_38905 contig_39380 contig_39384 contig_39657 contig_39905 contig_40140 contig_44383 contig_48574 contig_48576 contig_48599 contig_48603 contig_57059 contig_59619 contig_59629 contig_59631 contig_59632 contig_63024 contig_67956 scaffold_49001 scaffold_59627 contig_18119 contig_18123 contig_18973 contig_28140 contig_28142 contig_29347 contig_29348 contig_31958 contig_18654
grep 'bin3c.367' binning/DASTool/clr3.contigs.full_cluster_attribution.tsv | perl -ne 'chomp; @s = split(/\t/); print "$s[0]\n";' > example_clostridia/clr3.bin3c_367.contig.list

sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b assembly/clr3.contigs.fa.ccs.check.bam contig_33698 contig_37094 contig_38905 contig_39380 contig_39384 contig_39657 contig_39905 contig_40140 contig_44383 contig_48574 contig_48576 contig_48599 contig_48603 contig_57059 contig_59619 contig_59629 contig_59631 contig_59632 contig_63024 contig_67956 scaffold_49001 scaffold_59627 contig_18119 contig_18123 contig_18973 contig_28140 contig_28142 contig_29347 contig_29348 contig_31958 contig_18654 >  example_clostridia/clr3.367.ccs.bam"

sbatch -N 1 -n 30 --mem=35000 -p priority -q msn --wrap="samtools merge -@ 30 mapping/clr3.contigs/merged.bam mapping/clr3.contigs/Lib101_L1.bam mapping/clr3.contigs/Lib101_L2.bam mapping/clr3.contigs/Lib101_L3.bam mapping/clr3.contigs/Lib101_L4.bam"
sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --dependency=afterok:5668143 --wrap="samtools index mapping/clr3.contigs/merged.bam; samtools view -b mapping/clr3.contigs/merged.bam contig_33698 contig_37094 contig_38905 contig_39380 contig_39384 contig_39657 contig_39905 contig_40140 contig_44383 contig_48574 contig_48576 contig_48599 contig_48603 contig_57059 contig_59619 contig_59629 contig_59631 contig_59632 contig_63024 contig_67956 scaffold_49001 scaffold_59627 contig_18119 contig_18123 contig_18973 contig_28140 contig_28142 contig_29347 contig_29348 contig_31958 contig_18654 >  example_clostridia/clr3.367.sr.bam"

cp desman/bed_lists/clr3/bin3c.367.scg.bed example_clostridia/clr3.367.scg.bed

cp gtdbtk_bins/clr3.contigs/bin3c_367.fna example_clostridia/clr3.367.contigs.fa

for i in clr.1 clr.2 clr.3; do echo $i; for j in clr3; do echo $j; sbatch -N 1 -n 70 -p priority -q msn --mem=65000 -t 2-0 --wrap="minimap2 -ax map-pb -t 35 -R '@RG\tID:$i\tSM:$i' assembly/$j.contigs.fa /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/*_${i}.fastq.gz | samtools sort -@ 35 -T $j.$i -o minimap_clr_bams/$j.$i.sorted.bam -"; done; done
for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools index minimap_clr_bams/clr3.$i.sorted.bam; samtools view -b minimap_clr_bams/clr3.$i.sorted.bam contig_33698 contig_37094 contig_38905 contig_39380 contig_39384 contig_39657 contig_39905 contig_40140 contig_44383 contig_48574 contig_48576 contig_48599 contig_48603 contig_57059 contig_59619 contig_59629 contig_59631 contig_59632 contig_63024 contig_67956 scaffold_49001 scaffold_59627 contig_18119 contig_18123 contig_18973 contig_28140 contig_28142 contig_29347 contig_29348 contig_31958 contig_18654 > clr3.367.$i.bam"; done

samtools merge example_clostridia/clr3.367.clr.bam clr3.367.clr1.bam clr3.367.clr2.bam clr3.367.clr3.bam


##### getting stats


#####
# Note: these are the remaining contigs from clr2 for liz
grep 'bin3c.327' binning/DASTool/clr2.contigs.full_cluster_attribution.tsv | perl -ne 'chomp; @s = split(/\t/); print "$s[0] ";' ; echo
contig_11551 contig_11554 contig_20468 contig_22549 contig_27175 contig_28326 contig_28328 contig_33627 contig_34719 contig_34730 contig_34732 contig_34827 contig_37523 contig_38166 contig_43350 contig_45847 contig_48452 contig_50830 contig_51331 contig_51614 contig_51624 contig_51627 contig_52726 contig_59129 contig_64508 contig_65919

samtools view -b ../clr2.contigs.fasta.ccs.bam contig_11551 contig_11554 contig_20468 contig_22549 contig_27175 contig_28326 contig_28328 contig_33627 contig_34719 contig_34730 contig_34732 contig_34827 contig_37523 contig_38166 contig_43350 contig_45847 contig_48452 contig_50830 contig_51331 contig_51614 contig_51624 contig_51627 contig_52726 contig_59129 contig_64508 contig_65919 > example_clostridia/clr2.327.ccs.bam

samtools index ../mapping/clr2.contigs/merged.bam; samtools view -b ../mapping/clr2.contigs/merged.bam contig_11551 contig_11554 contig_20468 contig_22549 contig_27175 contig_28326 contig_28328 contig_33627 contig_34719 contig_34730 contig_34732 contig_34827 contig_37523 contig_38166 contig_43350 contig_45847 contig_48452 contig_50830 contig_51331 contig_51614 contig_51624 contig_51627 contig_52726 contig_59129 contig_64508 contig_65919 > example_clostridia/clr2.327.sr.bam

for i in clr1 clr2 clr3; do echo $i; sbatch -N 1 -n 1 --mem=9000 -p priority -q msn --wrap="samtools view -b ../minimap_clr_bams/clr2.$i.sorted.bam contig_11551 contig_11554 contig_20468 contig_22549 contig_27175 contig_28326 contig_28328 contig_33627 contig_34719 contig_34730 contig_34732 contig_34827 contig_37523 contig_38166 contig_43350 contig_45847 contig_48452 contig_50830 contig_51331 contig_51614 contig_51624 contig_51627 contig_52726 contig_59129 contig_64508 contig_65919 > clr2.327.$i.bam"; done

samtools merge example_clostridia/clr2.327.clr.bam clr2.327.clr1.bam clr2.327.clr2.bam clr2.327.clr3.bam
#####

# viruses
source activate /KEEP/rumen_longread_metagenome_assembly/seaborn/
mkdir viruses
for i in clr3; do echo $i; perl -ne '@s = split(/\t/); if($s[9] eq "Viruses"){print "$s[0]\t$s[1]\n";}' < blobtools/table.${i}.contigs.blobDB.table.txt > assembly/${i}.viral.contigs.list; done

for i in clr3; do echo $i; sbatch --nodes=1 --mem=30000 --ntasks-per-node=4 -p priority -q msn -J ${i}_virus --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/RumenLongReadASM/viralAssociationPipeline.py -a assembly/${i}.contigs.fa -b blobtools/table.${i}.contigs.blobDB.table.txt -i /lustre/project/forage_assemblies/sheep_project/complete_flye/clr3_rerun/mapping/${i}.contigs/hic/hic_Sau3AI.bam -v assembly/${i}.viral.contigs.list -l /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/sheep_poop_sequelII_CCS.fasta -m minimap2 -o viruses/${i}.contigs.vassoc"; done

perl -lane 'if($F[0] eq 'VirusCtg'){next;}else{print "$F[1]\t$F[0]\t$F[3]";}' < viruses/clr3.contigs.vassoc.final.tab > viruses/clr3.contigs.vassoc.assoc.tab

# Master table
perl -ne 'chomp; @F = split(/\t/); if($_ =~ /^##/){next;} for($x = 0; $x < scalar(@F); $x++){$F[$x] =~ s/\s+/_/g;} print join("\t", @F[0..2, 4..9, 12,15,18,21,24]) . "\n";' <blobtools/table.clr3.contigs.blobDB.table.txt > clr3.mastertable.template.tab
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j clr3_data_files_03_2021.json -o clr3_master_table_03_2021.tab -t clr3.mastertable.template.tab
```


### Basic, basic stats

```bash
for i in /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz /lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz; do echo $i; gunzip -c $i | perl -e '$c = 0; $bp = 0; while(<>){$s = <>; <>; <>; chomp($s); $c++; $bp += length($s);} print "$c\t$bp\n";'; done

/lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R1.fastq.gz
512187895       76972169332
/lustre/project/rumen_longread_metagenome_assembly/sheep_poop/illumina_wgs/LIB101996_R2.fastq.gz
512187895       76975222468

for i in *.fai; do echo $i; perl -lane 'if($F[1] < 1000){next;}else{print "$_";}' < $i > temp.gt1kb.fai; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i temp.gt1kb.fai; done
clr1.contigs.fasta.fai
Total length:   2984846227
Num contigs:    48338
N50 length:     1492428347
N50 value:      185361
L50 value:      3134
Max:    7141023
Min:    1000
clr2.contigs.fasta.fai
Total length:   3008462117
Num contigs:    48793
N50 length:     1504349034
N50 value:      187989
L50 value:      3112
Max:    6067615
Min:    1000
flye4.contigs.fasta.fai
Total length:   3424468269
Num contigs:    57259
N50 length:     1712285155
N50 value:      279880
L50 value:      2274
Max:    5546585
Min:    1000

for i in clr3_rerun/assembly/clr3.contigs.fa.fai; do echo $i; perl -lane 'if($F[1] < 1000){next;}else{print "$_";}' < $i > temp.gt1kb.fai; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/calculateContigN50.pl -i temp.gt1kb.fai; done
clr3_rerun/assembly/clr3.contigs.fa.fai
Total length:   2977556825
Num contigs:    52607
N50 length:     1488853997
N50 value:      181150
L50 value:      3100
Max:    5235210
Min:    1000

perl -lane 'if($F[1] > 90 && $F[2] < 5){print $_;}' < clr3_rerun/tables/clr3.contigs.ctg_cov_by_comp.simp.tab | wc -l
64

perl -lane 'if($F[1] > 90 && $F[2] < 5){print "$F[0]";}' < clr3_rerun/tables/clr3.contigs.ctg_cov_by_comp.simp.tab > clr3_rerun/tables/clr3.contigs.ctg_cov_by_comp.simp.tab.hqctgs.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/sheep_clr3/assembly_info.txt -c 0 -l clr3_rerun/tables/clr3.contigs.ctg_cov_by_comp.simp.tab.hqctgs.list -d '\t' | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 3 -d '\t' -m
|Entry | Value|
|:-----|-----:|
|N     |    42|
|Y     |    22|

perl -lane 'if($F[1] > 1000000){print $_;}' < /lustre/project/rumen_longread_metagenome_assembly/assemblies/sheep/sheep_clr3/assembly_info.txt | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 3 -d '\t' -m
|Entry | Value|
|:-----|-----:|
|N     |   258|
|Y     |    26|
```

#### IGV coordinates for figures

HiFi 451:   contig_23521:224,065-231,684
HiFi 451:   contig_15470:274,833-275,900

clr1 451: contig_37225:2,983-10,287  (is the analog of the first HiFi 451 region)



> /mnt/c/Sharedfolders

```bash
for i in clr1.451.ccs.bam clr1.451.sr.bam hifi.451.sr.bam hifi.451.ccs.bam; do echo $i; name=`echo $i | cut -d'.' -f1,2,3`; python3 filter_bam_by_coverage.py $i $name.filt.bam -c 0.9 --filter_secondary --filter_supp; done

for i in hifi*filt.bam; do name=`echo $i | cut -d'.' -f1,2,3,4`; base=`echo $i | cut -d'.' -f1,2`; echo $name; ctg="contig_23521"; start=227243; end=230008; python3 paint_bam_post_phaser.py $i $name.$ctg.$start.$end.bam magphase_v23_bhr/$base.strain.human_readable_by_read.txt -c $ctg -s $start -e $end; done 


for i in hifi*filt.bam; do name=`echo $i | cut -d'.' -f1,2,3,4`; base=`echo $i | cut -d'.' -f1,2`; echo $name; ctg="contig_15470"; start=274833; end=275900; python3 paint_bam_post_phaser.py $i $name.$ctg.$start.$end.bam magphase_v23_bhr/$base.strain.human_readable_by_read.txt -c $ctg -s $start -e $end; samtools index $name.$ctg.$start.$end.bam; done
```

#### Assigning empty taxa

```bash
grep -v 'NOBIN' flye4_master_table_12_2020.tab | perl -lane 'if(scalar(@F) < 22){print $_;}' | perl -e '%data; while(<>){chomp; @s = split(/\t/); $data{$s[13]}->{"d__$s[7];p__$s[8];o__$s[9];f__$s[10];g__$s[11];s__$s[12]"} += 1;} foreach my $k (keys(%data)){ print "$k"; foreach my $v (sort {$data{$k}->{$b} <=> $data{$k}->{$a}} keys(%{$data{$k}})){print "\t$v\t" . $data{$k}->{$v};} print "\n";}' > gtdbtk_output/flye4_nongtdbtk_bins.tab

```