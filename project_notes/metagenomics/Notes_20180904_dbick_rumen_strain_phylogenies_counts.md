# Strain haplotyping
---
*9/4/2018*



> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman

```bash
for i in illumina illumina_rest illumina_sole pacbio pacbio_rest pacbio_sole; do echo $i; echo -e "experiment\tbin\tstrains" >> ${i}.combined.strain.counts; for j in $i/*/*.strain.count; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\s+/); print "$ARGV[0]\t$s[0]\t$s[1]\n";} close IN;' $i $j >> ${i}.combined.strain.counts; done; done
```


#### Testing strain counts


```R
library(ggplot2)
library(dplyr)

combined <- data.frame(experiment=character(), bin=character(), strains=numeric(), stringsAsFactors = FALSE); for(i in c("illumina_rest.combined.strain.counts", "illumina.combined.strain.counts", "pacbio_sole.combined.strain.counts", "illumina_sole.combined.strain.counts", "pacbio_rest.combined.strain.counts", "pacbio.combined.strain.counts")){print(i); temp <- read.delim(i, header=TRUE); temp$strains <- as.numeric(temp$strains); combined <- bind_rows(combined, temp);}
ggplot(combined, aes(x=experiment, y=strains, fill=experiment)) + geom_violin()
dev.copy2pdf(file="first_round_strain_counts.pdf", useDingbats=FALSE)
```


I'm a little concerned that the "illumina" and "pacbio" runs were initiated with a different version of the scripts and that may be a cause of discrepancy in the strain counts I see. Let me rerun them to be sure.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/illumina

```bash
mkdir old
mv illumina_* ./old
cp ../illumina_rest/run_desman_pipeline_illumina_rest.sh ./

for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "illumina_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_illumina.sh $name $bed; done
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman/pacbio

```bash
mkdir old
mv pacbio_* ./old/

cp ../pacbio_rest/run_desman_pipeline_pacbio_rest.sh ./
for i in bin_lists/*.list; do name=`echo $i | perl -e '$h = <STDIN>; chomp($h); @bsegs = split(/[\._]/, $h); print "pacbio_$bsegs[4]\_$bsegs[5]";'`; echo $name; bed=`basename $i | cut -d'.' -f1,2`; echo $bed; sbatch run_desman_pipeline_pacbio.sh $name $bed; done
```

It looks like the values are the same. It's reproducible at least!

#### Preparing for alignments

I want to determine the minimum gene set for the alignment, and then determine the minimum alignment size for each SCG. This will be a crazy exercise, but that's why they pay me the small bucks!

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/desman

```bash
# I think that my venn list program can handle the large number of comparisons! 
for i in pacbio/pacbio_*; do file=`basename $i`; ls $i/scg_haps/*.fa | perl -lane '@dsegs = split(/\//, $F[0]); @bsegs = split(/[_\.]/, $dsegs[-1]); shift(@bsegs); shift(@bsegs); shift(@bsegs); pop(@bsegs); pop(@bsegs); print join("_", @bsegs);' > pacbio/$file.scg.list; done

# Removing zero value entries
rm pacbio/pacbio_metabat_345.scg.list pacbio/pacbio_hic_716.scg.list pacbio/pacbio_hic_630.scg.list pacbio/pacbio_hic_456.scg.list pacbio/pacbio_hic_438.scg.list pacbio/pacbio_hic_375.scg.list

# Removing single digit lists
rm pacbio/pacbio_metabat_374.scg.list pacbio/pacbio_metabat_263.scg.list pacbio/pacbio_metabat_233.scg.list pacbio/pacbio_hic_872.scg.list pacbio/pacbio_hic_732.scg.list pacbio/pacbio_hic_728.scg.list pacbio/pacbio_hic_585.scg.list pacbio/pacbio_hic_446.scg.list pacbio/pacbio_hic_395.scg.list pacbio/pacbio_hic_217.scg.list pacbio/pacbio_hic_200.scg.list

# Removing lists under 20 genes
rm pacbio/pacbio_metabat_573.scg.list pacbio/pacbio_metabat_506.scg.list pacbio/pacbio_metabat_498.scg.list pacbio/pacbio_metabat_172.scg.list pacbio/pacbio_metabat_082.scg.list pacbio/pacbio_hic_972.scg.list pacbio/pacbio_hic_772.scg.list pacbio/pacbio_hic_361.scg.list pacbio/pacbio_hic_358.scg.list pacbio/pacbio_hic_323.scg.list pacbio/pacbio_hic_283.scg.list pacbio/pacbio_hic_222.scg.list pacbio/pacbio_hic_1148.scg.list

# OK, that got me nowhere. Let's start from the >=70% complete lists from das_tool
for i in `perl -lane 'if($F[-1] < 5 && $F[-2] > 70){$F[0] =~ s/_final_public//; $F[0] =~ s/\./_/g; print $F[0];}' < ../dastool/pacbio_final_dastool_DASTool_summary.txt`; do file=`basename $i`; ls pacbio/$i/scg_haps/*.fa | perl -lane '@dsegs = split(/\//, $F[0]); @bsegs = split(/[_\.]/, $dsegs[-1]); shift(@bsegs); shift(@bsegs); shift(@bsegs); pop(@bsegs); pop(@bsegs); print join("_", @bsegs);' > pacbio/$file.scg.list; done
for i in pacbio/*.list; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $c = 0; while(<IN>){$c++;} close IN; if($c < 30){system("rm $ARGV[0]");}' $i ; done

for i in pacbio/*.list; do echo -n "$i "; done; echo
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl pacbio/pacbio_hic_1061.scg.list pacbio/pacbio_hic_110.scg.list pacbio/pacbio_hic_228.scg.list pacbio/pacbio_hic_239.scg.list pacbio/pacbio_hic_245.scg.list pacbio/pacbio_hic_272.scg.list pacbio/pacbio_hic_339.scg.list pacbio/pacbio_hic_506.scg.list pacbio/pacbio_hic_517.scg.list pacbio/pacbio_hic_617.scg.list pacbio/pacbio_hic_683.scg.list pacbio/pacbio_hic_839.scg.list pacbio/pacbio_metabat_066.scg.list pacbio/pacbio_metabat_162.scg.list pacbio/pacbio_metabat_265.scg.list pacbio/pacbio_metabat_309.scg.list pacbio/pacbio_metabat_434.scg.list pacbio/pacbio_metabat_493.scg.list pacbio/pacbio_metabat_576.scg.list pacbio/pacbio_metabat_636.scg.list

# That revealed about 9 of the genes common to a large number of datasets.
# Let's flip the lists though
perl -e '@files = `ls pacbio/*.scg.list`; chomp; %h; foreach my $f (@files){open(IN, "< $f"); while(<IN>){chomp; $h{$_}->{$f} = 1;} close IN;} foreach my $k (sort {scalar(keys($h{$b})) <=> scalar(keys($h{$a}))} keys(%h)){ print "$k\t" . scalar(keys($h{$k})) . "\n";}'
alanyl_tRNA_synthetase  19
ribosomal_protein_L19   18
arginyl_tRNA_synthetase 17
Phenylalanyl-tRNA_synthetase_alpha      17
ribosomal_protein_L14   16
gyrA    16
tRNA_N6-adenosine_threonylcarbamoyltransferase  16


# OK, let's try this with Illumina now
for i in `perl -lane 'if($F[-1] < 5 && $F[-2] > 70){$F[0] =~ s/_megahit_public//; if($F[0] =~ /hic/){$F[0] =~ s/_megahit_hic//; $F[0] .= "_hqdas";} $F[0] =~ s/\./_/g; print $F[0];}' < ../dastool/illumina_megahit_dastool_DASTool_summary.txt`; do file=`basename $i`; ls illumina/$i/scg_haps/*.fa | perl -lane '@dsegs = split(/\//, $F[0]); @bsegs = split(/[_\.]/, $dsegs[-1]); shift(@bsegs); shift(@bsegs); shift(@bsegs); shift(@bsegs); pop(@bsegs); pop(@bsegs); print join("_", @bsegs);' > illumina/$file.scg.list; done
for i in illumina/*.list; do perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $c = 0; while(<IN>){$c++;} close IN; if($c < 30){system("rm $ARGV[0]");}' $i ; done

for i in illumina/*.list; do echo -n "$i "; done; echo
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina/illumina_047_hqdas.scg.list illumina/illumina_072_hqdas.scg.list illumina/illumina_100_hqdas.scg.list illumina/illumina_381_hqdas.scg.list illumina/illumina_397_hqdas.scg.list illumina/illumina_metabat_018.scg.list illumina/illumina_metabat_020.scg.list illumina/illumina_metabat_039.scg.list illumina/illumina_metabat_1030.scg.list illumina/illumina_metabat_1031.scg.list illumina/illumina_metabat_1040.scg.list illumina/illumina_metabat_1056.scg.list illumina/illumina_metabat_1080.scg.list illumina/illumina_metabat_1085.scg.list illumina/illumina_metabat_1095.scg.list illumina/illumina_metabat_1097.scg.list illumina/illumina_metabat_1112.scg.list illumina/illumina_metabat_1119.scg.list illumina/illumina_metabat_1134.scg.list illumina/illumina_metabat_1148.scg.list illumina/illumina_metabat_117.scg.list illumina/illumina_metabat_1196.scg.list illumina/illumina_metabat_1226.scg.list illumina/illumina_metabat_1234.scg.list illumina/illumina_metabat_1239.scg.list illumina/illumina_metabat_1250.scg.list illumina/illumina_metabat_1278.scg.list illumina/illumina_metabat_1281.scg.list illumina/illumina_metabat_1305.scg.list illumina/illumina_metabat_1316.scg.list illumina/illumina_metabat_172.scg.list illumina/illumina_metabat_193.scg.list illumina/illumina_metabat_203.scg.list illumina/illumina_metabat_225.scg.list illumina/illumina_metabat_244.scg.list illumina/illumina_metabat_275.scg.list illumina/illumina_metabat_280.scg.list illumina/illumina_metabat_288.scg.list illumina/illumina_metabat_336.scg.list illumina/illumina_metabat_406.scg.list illumina/illumina_metabat_421.scg.list illumina/illumina_metabat_485.scg.list illumina/illumina_metabat_510.scg.list illumina/illumina_metabat_572.scg.list illumina/illumina_metabat_613.scg.list illumina/illumina_metabat_655.scg.list illumina/illumina_metabat_665.scg.list illumina/illumina_metabat_734.scg.list illumina/illumina_metabat_743.scg.list illumina/illumina_metabat_804.scg.list illumina/illumina_metabat_889.scg.list illumina/illumina_metabat_896.scg.list illumina/illumina_metabat_981.scg.list illumina/illumina_metabat_997.scg.list

# 9 of the genes are unique to illumina/illumina_metabat_1040.scg.list!!
perl -e '@files = `ls illumina/*.scg.list`; chomp; %h; foreach my $f (@files){open(IN, "< $f"); while(<IN>){chomp; $h{$_}->{$f} = 1;} close IN;} foreach my $k (sort {scalar(keys($h{$b})) <=> scalar(keys($h{$a}))} keys(%h)){ print "$k\t" . scalar(keys($h{$k})) . "\n";}'
Preprotein_translocase_subunit_SecY     49
ribosomal_protein_L5    48
gyrA    46
ribosomal_protein_L6P-L9E       46
alanyl_tRNA_synthetase  46
ribosomal_protein_S3    46
ribosomal_protein_L18   45
ribosomal_protein_L14   45
recA    45
arginyl_tRNA_synthetase 45
ribosomal_protein_L16-L10E      45
Methionyl-tRNA_synthetase       45
ribosomal_protein_S8    45
ribosomal_protein_L2    45


# What the hell, let's just try to pick the haplotypes that belong to gyr A and make the phylogenies from there
grep 'gyrA' pacbio/*.scg.list | perl -e '%h; while(<>){chomp; @s = split(/[\/\.]/); $h{$s[1]} = 1;} foreach my $k (keys(%h)){system("cat pacbio/$k/scg_haps/*_gyrA_hap.fa >> pacbio/pacbio_scg_gyra_haps.fa");}'
# Removing the segments that were smaller than 500 bp
vim pacbio/pacbio_scg_gyra_haps.fa
# Now to align the fasta and try to find the minimum length

```

## Counting ORFs near the peripheries of contigs

I need to see if the partial ORFs are more frequently found near the edges of contigs, which would explain why the illumina assembly has more of them!

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[0] ne "name"){print "$F[0]\t$F[1]\n";}' < ../master_tables/illumina_megahit_master_table_2018_09_07.tab > illumina_an_bins_contiglens.tab
perl count_partial_ORFs_near_contig_ends.pl -s illumina_megahit_prodigal_proteins.shortform.tab -f illumina_an_bins_contiglens.tab -o illumina_an_bins_orfs

perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[0] ne "name"){print "$F[0]\t$F[1]\n";}' < ../master_tables/pacbio_final_pilon_master_table_2018_09_07.tab > pacbio_an_bins_contiglens.tab
perl count_partial_ORFs_near_contig_ends.pl -s pacbio_final_prodigal_proteins.shortform.tab -f pacbio_an_bins_contiglens.tab -o pacbio_an_bins_orfs

```

```R
illumina <- read.delim("illumina_an_bins_orfs.partial.cats", header=FALSE)
colnames(illumina) <- c("Location", "Completion")

illumina$CompCat[illumina$Completion == 0] <- "Full"
illumina$CompCat[illumina$Completion == 1] <- "MissEnd"
illumina$CompCat[illumina$Completion == 10] <- "MissStart"
illumina$CompCat[illumina$Completion == 11] <- "MissBoth"

illumina$CompCat <- as.factor(illumina$CompCat)
library(MASS)
chisq.test(table(illumina$Location, illumina$CompCat))

        Pearson's Chi-squared test

data:  table(illumina$Location, illumina$CompCat)
X-squared = 2826800, df = 6, p-value < 2.2e-16

library(graphics)
library(ggplot2)
pdf(file="illumina_orf_categories_mosaic_plot.pdf", useDingbats=FALSE)
mosaicplot(table(illumina$Location, illumina$CompCat), shade=TRUE, las=2, main="Illumina ORF Positions on Contigs vs Completeness")
dev.off()

pacbio <- read.delim("pacbio_an_bins_orfs.partial.cats", header=FALSE)
colnames(pacbio) <- c("Location", "Completion")
pacbio$CompCat[pacbio$Completion == 0] <- "Full"
pacbio$CompCat[pacbio$Completion == 1] <- "MissEnd"
pacbio$CompCat[pacbio$Completion == 10] <- "MissStart"
pacbio$CompCat[pacbio$Completion == 11] <- "MissBoth"

chisq.test(table(pacbio$Location, pacbio$CompCat))

        Pearson's Chi-squared test

data:  table(pacbio$Location, pacbio$CompCat)
X-squared = 1680200, df = 6, p-value < 2.2e-16

pdf(file="pacbio_orf_categories_mosaic_plot.pdf", useDingbats=FALSE)
mosaicplot(table(pacbio$Location, pacbio$CompCat), shade=TRUE, las=2, main="Pacbio ORF Positions on Contigs vs Completeness")
dev.off()

table(pacbio$Location, pacbio$CompCat)

          Full MissBoth MissEnd MissStart
  End     3905        0   34730         0
  Mid   860711        0       0         0
  Start   3813       56       6     34615

table(illumina$Location, illumina$CompCat)

          Full MissBoth MissEnd MissStart
  End    12510        0  249262         0
  Mid   959925        0       0         0
  Start  12188    23458    1395    248553
```

## Grabbing only the high quality bins

I need to pull out the high quality bins from both datasets and characterize the contents.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Illumina first 
perl -lane 'print $F[1];' < ../dastool/illumina_dastool_high_quality_dasbins.contigs.tab | sort | uniq | perl -ne 'chomp; @s = split(/_/); print "$s[1]\n";' > illumina_dastool_high_quality_dasbins.nums.list

python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnGrep.py -f illumina_megahit_master_table_2018_09_07.tab -c 51 -l illumina_dastool_high_quality_dasbins.nums.list > illumina_megahit_master_table_2018_09_07.HQbins.only.tab

# OK, that didn't work out because of the bin naming conventions I used.
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print "$_\n";}' < illumina_megahit_master_table_2018_09_07.tab > illumina_megahit_master_table_2018_09_07.HQbins.only.tab

perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,51,52,53,57,58]) . "\n";' < illumina_megahit_master_table_2018_09_07.HQbins.only.tab > illumina_megahit_master_table_2018_09_07.HQbins.only.short.tab

# Pacbio next
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print "$_\n";}' < pacbio_final_pilon_master_table_2018_09_07.tab > pacbio_final_pilon_master_table_2018_09_07.HQbins.only.tab

perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,51,52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_2018_09_07.HQbins.only.tab > pacbio_final_pilon_master_table_2018_09_07.HQbins.only.short.tab

perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[12]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < illumina_megahit_master_table_2018_09_07.HQbins.only.short.tab | perl -lane 'my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d);' > illumina_megahit_master_table_2018_09_07.HQbins.taxtable.tab

perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[12]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < pacbio_final_pilon_master_table_2018_09_07.HQbins.only.short.tab | perl -lane 'my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d);' > pacbio_final_pilon_master_table_2018_09_07.HQbins.taxtable.tab

perl -lane 'my @counts; for($x = 1; $x < scalar(@F); $x++){@j = split(/;/, $F[$x]); push(@counts, scalar(@j));} print "$F[0]\t" . join("\t", @counts);' < pacbio_final_pilon_master_table_2018_09_07.HQbins.taxtable.tab > pacbio_final_pilon_master_table_2018_09_07.HQbins.taxcounts
perl -lane 'my @counts; for($x = 1; $x < scalar(@F); $x++){@j = split(/;/, $F[$x]); push(@counts, scalar(@j));} print "$F[0]\t" . join("\t", @counts);' < illumina_megahit_master_table_2018_09_07.HQbins.taxtable.tab > illumina_megahit_master_table_2018_09_07.HQbins.taxcount

### the AN bin datasets ###
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print "$_\n";}' < illumina_megahit_master_table_2018_09_07.tab > illumina_megahit_master_table_2018_09_07.ANbins.tab
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50, 51, 52,53,57,58]) . "\n";' < illumina_megahit_master_table_2018_09_07.ANbins.tab > illumina_megahit_master_table_2018_09_07.ANbins.short.tab
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[12]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < illumina_megahit_master_table_2018_09_07.ANbins.short.tab |perl -lane 'my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d);' > illumina_megahit_master_table_2018_09_07.ANbins.taxtable.tab
perl -lane 'my @counts; for($x = 1; $x < scalar(@F); $x++){@j = split(/;/, $F[$x]); push(@counts, scalar(@j));} print "$F[0]\t" . join("\t", @counts);' < illumina_megahit_master_table_2018_09_07.ANbins.taxtable.tab > illumina_megahit_master_table_2018_09_07.ANbins.taxcounts

perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print "$_\n";}' < pacbio_final_pilon_master_table_2018_09_07.tab > pacbio_final_pilon_master_table_2018_09_07.ANbins.tab
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50, 51, 52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_2018_09_07.ANbins.tab > pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[12]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab |perl -lane 'my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d);' > pacbio_final_pilon_master_table_2018_09_07.ANbins.taxtable.tab
perl -lane 'my @counts; for($x = 1; $x < scalar(@F); $x++){@j = split(/;/, $F[$x]); push(@counts, scalar(@j));} print "$F[0]\t" . join("\t", @counts);' < pacbio_final_pilon_master_table_2018_09_07.ANbins.taxtable.tab > pacbio_final_pilon_master_table_2018_09_07.ANbins.taxcounts
```

Getting brief summaries of the data for the manuscript.

```R
ilmn <- read.delim("illumina_megahit_master_table_2018_09_07.HQbins.only.short.tab",header=TRUE)
           name          length             GC            cov_sum
 k127_1000532:   1   Min.   :  1000   Min.   :0.1815   Min.   :   8.208
 k127_1001057:   1   1st Qu.:  3445   1st Qu.:0.4155   1st Qu.:  45.613
 k127_1002355:   1   Median :  5361   Median :0.4837   Median :  90.898
 k127_1003354:   1   Mean   :  8793   Mean   :0.4758   Mean   : 173.495
 k127_1003632:   1   3rd Qu.:  9485   3rd Qu.:0.5412   3rd Qu.: 240.206
 k127_1003884:   1   Max.   :320393   Max.   :0.6813   Max.   :3425.903
 (Other)     :9227
 superkingdom.t.24         phylum.t.28                 order.t.32
 Archaea  : 199    Bacteroidetes :4360   Bacteroidales      :3860
 Bacteria :8748    Firmicutes    :3017   Clostridiales      :2485
 Eukaryota:  26    Bacteria-undef: 423   Bacteria-undef     : 423
 no-hit   : 258    Spirochaetes  : 390   Spirochaetales     : 370
 Viruses  :   2    no-hit        : 258   Bacteroidetes-undef: 307
                   Synergistetes : 245   Firmicutes-undef   : 288
                   (Other)       : 540   (Other)            :1500
   MetabatBin       HiCBin      CompleteORFs      PartialORFs
 NOBIN  : 570   NOBIN  :7382   Min.   :  0.000   Min.   :0.000
 82     : 393   25     : 471   1st Qu.:  2.000   1st Qu.:1.000
 752    : 362   114    : 437   Median :  4.000   Median :2.000
 1302   : 346   238    : 179   Mean   :  6.835   Mean   :1.485
 601    : 337   253    : 141   3rd Qu.:  8.000   3rd Qu.:2.000
 121    : 319   1389   :  83   Max.   :238.000   Max.   :2.000
 (Other):6906   (Other): 540
      MICKRMGAligns              Hungate1000Aligns
 -           :4263   -                    :7019
 scaffold_93 :  27   3960633.scaffold00002:  54
 scaffold_15 :  23   3960633.scaffold00001:  53
 scaffold_60 :  23   3960793.scaffold00003:  34
 scaffold_211:  22   3960633.scaffold00004:  30
 scaffold_53 :  21   3960793.scaffold00001:  30
 (Other)     :4854   (Other)              :2013

ilan <- read.delim("illumina_megahit_master_table_2018_09_07.ANbins.short.tab", header=TRUE)

# Now for the pacbio
pb <- read.delim("pacbio_final_pilon_master_table_2018_09_07.HQbins.only.short.tab", header=TRUE)
          name          length             GC            cov_sum
 tig00002456:   1   Min.   :  1759   Min.   :0.2743   Min.   :  18.74
 tig00002711:   1   1st Qu.: 11038   1st Qu.:0.4526   1st Qu.: 155.46
 tig00002903:   1   Median : 22108   Median :0.4855   Median : 300.64
 tig00002934:   1   Mean   : 38997   Mean   :0.4693   Mean   : 562.54
 tig00003179:   1   3rd Qu.: 48016   3rd Qu.:0.5051   3rd Qu.: 691.33
 tig00003214:   1   Max.   :710205   Max.   :0.6193   Max.   :4648.39
 (Other)    :1216
 superkingdom.t.24         phylum.t.28                order.t.32
 Archaea  :  20    Bacteroidetes :946   Bacteroidales      :893
 Bacteria :1188    Firmicutes    :180   Veillonellales     : 82
 Eukaryota:   8    Bacteria-undef: 41   Clostridiales      : 81
 no-hit   :   6    Euryarchaeota : 20   Bacteroidetes-undef: 46
                   Proteobacteria: 13   Bacteria-undef     : 41
                   no-hit        :  6   Methanobacteriales : 19
                   (Other)       : 16   (Other)            : 60

   MetabatBin      HiCBin     CompleteORFs     PartialORFs
 410    :185   55     :125   Min.   :  0.00   Min.   :0.000
 NOBIN  :116   10     : 96   1st Qu.: 11.00   1st Qu.:1.000
 246    :113   18     : 92   Median : 21.00   Median :1.000
 619    : 94   NOBIN  : 61   Mean   : 36.08   Mean   :1.164
 633    : 79   20     : 48   3rd Qu.: 43.75   3rd Qu.:2.000
 491    : 76   3      : 48   Max.   :663.00   Max.   :2.000
 (Other):559   (Other):752
     MICKRMGAligns                        Hungate1000Aligns
 -          : 125   -                              :425
 scaffold_55:  11   3960801.scaffold00001          :  9
 scaffold_22:   8   3960801.scaffold00004          :  9
 scaffold_36:   8   4304392.scaffold00001          :  7
 scaffold_59:   8   3960633.scaffold00002          :  5
 scaffold_23:   6   3353505.scf7180000000012|quiver:  4
 (Other)    :1056   (Other)                        :763

pban <- read.delim("pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab", header=TRUE)

sum(pb[pb$MICKRMGAligns == "-" & pb$Hungate1000Aligns == "-", 2])
[1] 1703508 <- length of contigs with no alignments to either dataset
sum(ilmn[ilmn$MICKRMGAligns == "-" & ilmn$Hungate1000Aligns == "-", 2])
[1] 23553844

# Now for the ANbin counts
sum(ilan[ilan$MICKRMGAligns == "-" & ilan$Hungate1000Aligns == "-", 2])
[1] 669890409
sum(pban[pban$MICKRMGAligns == "-" & pban$Hungate1000Aligns == "-", 2])
[1] 137209071

nrow(ilan[ilan$MICKRMGAligns == "-" & ilan$Hungate1000Aligns == "-",])
[1] 207599
nrow(pban[pban$MICKRMGAligns == "-" & pban$Hungate1000Aligns == "-",])
[1] 12421

# Calculating data novel to each dataset from the AN bins
ilall <- read.delim("illumina_megahit_master_table_2018_09_07.ANbins.tab", header=TRUE)
pball <- read.delim('pacbio_final_pilon_master_table_2018_09_07.tab', header=TRUE)

nrow(ilall[ilall$MICKRMGAligns == "-" & ilall$Hungate1000Aligns == "-" & ilall$PacBioCtgAligns == "-" & ilall$DASBin_AN != "NOBIN",])
[1] 152793
sum(ilall[ilall$MICKRMGAligns == "-" & ilall$Hungate1000Aligns == "-" & ilall$PacBioCtgAligns == "-" & ilall$DASBin_AN != "NOBIN",2])
[1] 435645025

nrow(pball[pball$IlluminaCtgAligns == "-" & pball$MICKRMGAligns == "-" & pball$Hungate1000Aligns == "-" & pball$DASBin_AN != "NOBIN",])
[1] 186
sum(pball[pball$IlluminaCtgAligns == "-" & pball$MICKRMGAligns == "-" & pball$Hungate1000Aligns == "-" & pball$DASBin_AN != "NOBIN",2])
[1] 1186026


# Testing GC
summary(pball[pball$IlluminaCtgAligns == "-" & pball$MICKRMGAligns == "-" & pball$Hungate1000Aligns == "-" & pball$DASBin_AN != "NOBIN",])
Mean   :0.5010
Median :0.5074

summary(pball$GC)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.1544  0.4313  0.4898  0.4795  0.5363  0.9799

ks.test(pball$GC, pball[pball$IlluminaCtgAligns == "-" & pball$MICKRMGAligns == "-" & pball$Hungate1000Aligns == "-" & pball$DASBin_AN != "NOBIN",3], alternative="greater")
        Two-sample Kolmogorov-Smirnov test

data:  pball$GC and pball[pball$IlluminaCtgAligns == "-" & pball$MICKRMGAligns == pball$GC and     "-" & pball$Hungate1000Aligns == "-" & pball$DASBin_AN != pball$GC and     "NOBIN", 3]
D^+ = 0.16495, p-value = 4.12e-05
alternative hypothesis: the CDF of x lies above that of y



icount <- read.delim("illumina_megahit_master_table_2018_09_07.HQbins.taxcount", header=FALSE)
summary(icount)
           V1           V2              V3           V4               V5
 1          : 1   Min.   :1.000   Min.   : 1   Min.   : 1.000   Min.   : 1.00
 DASBin_HQ  : 1   1st Qu.:2.000   1st Qu.: 3   1st Qu.: 5.000   1st Qu.: 5.00
 hic_047    : 1   Median :2.000   Median : 4   Median : 7.000   Median : 9.00
 hic_397    : 1   Mean   :2.327   Mean   : 5   Mean   : 7.959   Mean   :10.21
 metabat_018: 1   3rd Qu.:3.000   3rd Qu.: 6   3rd Qu.:11.000   3rd Qu.:15.00
 metabat_020: 1   Max.   :4.000   Max.   :26   Max.   :20.000   Max.   :27.00
 (Other)    :43                                                 NA's   :1
       V6              V7       
 Min.   : 1.00   Min.   : 1.000 
 1st Qu.: 5.75   1st Qu.: 5.000  
 Median :13.50   Median : 8.500  
 Mean   :13.50   Mean   : 9.271  
 3rd Qu.:21.00   3rd Qu.:12.000  
 Max.   :42.00   Max.   :43.000  
 NA's   :1       NA's   :1      

summary(pcount)
         V1           V2             V3              V4
 7        : 1   Min.   :1.00   Min.   :1.000   Min.   : 1.000
 8        : 1   1st Qu.:1.00   1st Qu.:1.000   1st Qu.: 1.000
 DASBin_HQ: 1   Median :1.00   Median :2.500   Median : 3.000
 hic_110  : 1   Mean   :1.75   Mean   :2.667   Mean   : 3.333
 hic_195  : 1   3rd Qu.:2.00   3rd Qu.:3.000   3rd Qu.: 4.000
 hic_217  : 1   Max.   :9.00   Max.   :8.000   Max.   :11.000
 (Other)  :19   NA's   :1      NA's   :1       NA's   :1
       V5               V6               V7
 Min.   : 1.000   Min.   : 1.000   Min.   : 1.000
 1st Qu.: 3.750   1st Qu.: 4.000   1st Qu.: 5.000
 Median : 5.000   Median : 7.000   Median : 6.500
 Mean   : 5.583   Mean   : 7.333   Mean   : 6.458
 3rd Qu.: 8.000   3rd Qu.:10.000   3rd Qu.: 8.000
 Max.   :15.000   Max.   :17.000   Max.   :13.000
 NA's   :1        NA's   :1        NA's   :1
```

## Coverage metrics to indicate difference in sample composition

Let's use coverage averages/information to try to distinguish our samples from previously sequenced datasets

> Ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,8,9,14,15,16,17,21,23,39,47,50,57,58]) . "\n";' < illumina_megahit_master_table_2018_09_07.ANbins.tab > illumina_megahit_master_table_2018_09_07.ANbins.rd.tab
perl -ne 'chomp; @F = split(/\t/); print join("\t", @F[0,1,2,8,9,14,15,16,17,21,23,39,47,50,57,58]) . "\n";' < pacbio_final_pilon_master_table_2018_09_07.ANbins.tab > pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.tab
```

```R
library(dplyr)
library(vegan)

ilmn <- read.delim("illumina_megahit_master_table_2018_09_07.ANbins.rd.tab", header=TRUE)
rddat <- ilmn[,c(4,5,6,7,8,9,10)]
rownames(rddat) <- ilmn$name
rddat <- t(rddat)

zscores <- t(scale(rddat, center = TRUE, scale = TRUE))
summary(zscores)
     cov16             cov17              cov0              cov1
 Min.   :-1.5398   Min.   :-0.9301   Min.   :-1.7249   Min.   :-0.7967
 1st Qu.:-0.5592   1st Qu.:-0.1705   1st Qu.:-0.5753   1st Qu.: 0.6235
 Median :-0.4520   Median : 0.4655   Median :-0.4790   Median : 1.6736
 Mean   :-0.2732   Mean   : 0.7223   Mean   :-0.4016   Mean   : 1.3755
 3rd Qu.:-0.2525   3rd Qu.: 1.7243   3rd Qu.:-0.3525   3rd Qu.: 2.1579
 Max.   : 2.2678   Max.   : 2.2678   Max.   : 2.2670   Max.   : 2.2678
      cov2              cov3              cov7
 Min.   :-1.5211   Min.   :-1.3831   Min.   :-1.8757
 1st Qu.:-0.6422   1st Qu.:-0.5932   1st Qu.:-0.7615
 Median :-0.5396   Median :-0.4551   Median :-0.6384
 Mean   :-0.5185   Mean   :-0.2646   Mean   :-0.6399
 3rd Qu.:-0.4429   3rd Qu.:-0.2357   3rd Qu.:-0.5373
 Max.   : 2.2526   Max.   : 2.2676   Max.   : 2.2677

zscores.rda <- rda(zscores)
pdf(file="test_total_zscore_pca.pdf", useDingbats=FALSE)
biplot(zscores.rda, display = c("sites", "species"), type= c("text", "points"))
dev.off()

rddat <- ilmn[,c(4,5,6,7,8,9,10)]
rddat.norm <- wisconsin(rddat)
rddat.norm.rda
pdf(file="test_total_zscore_pca.pdf", useDingbats=FALSE)
biplot(rddat.norm.rda, display = c("sites", "species"), type= c("text", "points"))
ordihull(rddat.norm.rda, group=ilmn$superkingdom.t.24)
dev.off()

```

Let's try a python approach, because the data manipulation in R is incredibly tedious and unintuitive. I am using the hypergeomtric test to assess enrichment of specific taxonomic groupings in each dataset.

```bash
# Genus illumina associations
python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/metagenomics/readdepthHyperGeomEnrichmentTest.py -f illumina_megahit_master_table_2018_09_07.ANbins.rd.tab -c 3,5,6,7,8,9 -s 4 -g 11 -o illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo -m 10,12,13,14
wc -l illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo
1942 illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo

# Since there were 1942 tests, and I want to use an alpha of 0.05, the Benjamini correction cutoff is 2.57 x 10-5
perl -ne 'chomp; @F=split(/\t/); if($F[4] < 0.000026){print "$_\n";}' < illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo | wc -l
205

# Genus pacbio associations
python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/metagenomics/readdepthHyperGeomEnrichmentTest.py -f pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.tab -c 3,5,6,7,8,9 -s 4 -g 11 -o pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo -m 10,12,13,14

wc -l pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo
958 pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo

# So the Benjamini correction is 5 x 10-5
perl -ne 'chomp; @F=split(/\t/); if($F[4] < 0.000052){print "$_\n";}' < pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo | wc -l
44

# Let's see how many geni are common between the two datasets
perl -ne 'chomp; @F=split(/\t/); if($F[10] eq "True" && $F[4] ne "Group"){print "$F[4]\n";}' < pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo > pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo.list
perl -ne 'chomp; @F=split(/\t/); if($F[10] eq "True" && $F[4] ne "Group"){print "$F[4]\n";}' < illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo > illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo.list

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo.list pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo.list
File Number 1: illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo.list
File Number 2: pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo.list
Set     Count
1       120
1;2     29
2       27

# 29! Impressive! And 27 unique to pacbio! Let's see what they are
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_2 illumina_megahit_master_table_2018_09_07.ANbins.rd.hypgeo.list pacbio_final_pilon_master_table_2018_09_07.ANbins.rd.hypgeo.list
Phi29virus
Anaeromyces
Ichthyophthirius
Candidatus Izimaplasma
no-hit
Elusimicrobium
Microscilla
Stentor
Enterococcus
Candidatus Pelagibacter
Clostridium
Tetrahymena
Eubacterium
Bacillus
Acetobacter
Paramecium
Coprobacillus
Azospirillum
Neocallimastix
Bacteroidetes-undef
Proteobacteria-undef
Gossypium
Stylonychia
Dysgonomonas
Oxytricha
Phikzvirus
Thermoplasmatales-undef
Euryarchaeota-undef
Pseudocohnilembus
```

## Generating ARG figure

I want to summarize the ARG genes in a way that shows the different classes of genes present, while showing the quality of the PacBio calls.

> pwd: F:/SharedFolders/metagenomics/pilot_manuscript/figure_drafts/amr_genes

```R
library(ggplot2)
library(gridExtra)
library(scales)
data <- read.delim("amr_gene_columns.txt")

# The following plot was too dense
ggplot(data, aes(x=qcoverage, color=ARG, fill=ARG)) + geom_histogram(position="identity", alpha=0.5) + scale_color_brewer(palette = "dark2") + theme_bw() + facet_grid(~Tech)

# This created a servicable but boring plot
ggplot(data, aes(x=ARGClass, color=ARGClass, fill=ARGClass)) + geom_bar(position="stack") + scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + facet_grid(~Tech)

# OK, I think that this is the only way to show this data in a reasonable way
upper <- ggplot(data, aes(x=ARGClass, fill=ARGClass)) + geom_bar(position="stack") + scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits=c(20,175), oob=rescale_none) + facet_grid(~Tech)
lower <- ggplot(data, aes(x=ARGClass, fill=ARGClass)) + geom_bar(position="stack") + scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + scale_y_continuous(limits=c(0,15), oob= rescale_none) + facet_grid(~Tech)
pdf(file="arg_gene_class_counts.pdf", useDingbats = FALSE)
grid.arrange(upper, lower, ncol=1)

### Bradd asked me to change the plot so that it displays the high confidence AMR genes
data <- read.delim("ResultsGreaterThan90QC.txt")
upper <- ggplot(data, aes(x=Class, fill=Class)) + geom_bar(position="stack") + scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits=c(20,90), oob=rescale_none) + facet_grid(~Tech)
lower <- ggplot(data, aes(x=Class, fill=Class)) + geom_bar(position="stack") + scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + scale_y_continuous(limits=c(0,15), oob= rescale_none) + facet_grid(~Tech)
pdf(file="arg_gene_class_counts_revised.pdf", useDingbats = FALSE)
grid.arrange(upper, lower, ncol=1)
```

## Generating a hive plot

This may be the death of me, but I want to convey the relatedness of our datasets to Mick's and the Hungate1000's datasets. I had a crazy (crazy!) idea to use a hive plot instead of a network to demonstrate this. I plan on making two hives, one for each assembly, and plotting the proportional counts of bases onto each.

First, I need to condense the data into a non-overlapping set of edges and nodes with weights that correspond to a proportional count.

OK, there are no good tools for modifying edge weight scale to make the sort of "proportional" graph that I saw in the hiveplot promotional webpage. I may have to generate a "hack" to make this work. Basically, I divide the plot into nodes of 100 for each super-set and use those as proxies for the proportions. I would then have to sort the nodes by numerical order to create the illusion of a propotional comparison.

Thankfully, the data is only going to be a "does it align?" or a "it aligns" binary comparison, and I will use the filtered PAF files from my previous alignments to make things easier.

Let's start by queueing the Hungate vs Mick alignment and proceed from there.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Hungate vs Mick
module load minimap2/2.6
module load bedtools

# Because of the way that the fasta headers were named, I need to add unique IDs to each file
for i in /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/genomes/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; perl -e 'chomp @ARGV; open(IN, "< $ARGV[0]"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; $_ = ">" . $ARGV[1] . "_$_\n"; print $_;}else{print $_;}}' $i $name >> /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug_reformat.fasta; done

sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug_reformat.fasta > hungate_vs_mick.paf"
sbatch --nodes=1 --ntasks-per-node=5 --mem=20000 -p short --wrap="minimap2 -x asm5 -t 5 /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug_reformat.fasta /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa > mickreformat_vs_ilmn.paf"

# Now, how to process the PAF files? Keep it simple?
# Let's use my previous threshold. > 100 bp alignment length and >0 mapping quality
# Actually, let's incorporate the AN bin filter as well
sbatch -p short convert_paf_to_beds.pl -p /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_vs_ilmn_megahit.paf -t hungate_ilmn_paf_filt.bed -q ilmn_hungate_paf_filt.bed -f ../dastool/illumina_dastool_analysis_binset_lt10redund.bins
sbatch -p short convert_paf_to_beds.pl -p mickreformat_vs_ilmn.paf -t mick_ilmn_paf_filt.bed -q ilmn_mick_paf_filt.bed -f ../dastool/illumina_dastool_analysis_binset_lt10redund.bins

bedtools sort -i ilmn_hungate_paf_filt.bed | bedtools merge -i stdin > ilmn_hungate_paf_filt.sort.merge.bed
bedtools sort -i hungate_ilmn_paf_filt.bed | bedtools merge -i stdin > hungate_ilmn_paf_filt.sort.merge.bed

cat hungate_ilmn_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           118637195	# Hungate asm in ilmn AN bins
cat ilmn_hungate_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           172100230	# Ilmn AN bins in Hungate asm

bedtools sort -i ilmn_mick_paf_filt.bed | bedtools merge -i stdin > ilmn_mick_paf_filt.sort.merge.bed
bedtools sort -i mick_ilmn_paf_filt.bed | bedtools merge -i stdin > mick_ilmn_paf_filt.sort.merge.bed

cat ilmn_mick_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           411825315	# Ilmn AN bins in Mick asm
cat mick_ilmn_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           342782266	# Mick asm in Ilmn AN bins

# Now to calculate the universal hungate vs mick data
sbatch -p short convert_paf_to_beds.pl -p hungate_vs_mick.paf -t hungate_vs_mick_paf_filt.bed -q mick_vs_hungate_paf_filt.bed

bedtools sort -i hungate_vs_mick_paf_filt.bed | bedtools merge -i stdin > hungate_vs_mick_paf_filt.sort.merge.bed
bedtools sort -i mick_vs_hungate_paf_filt.bed | bedtools merge -i stdin > mick_vs_hungate_paf_filt.sort.merge.bed

cat hungate_vs_mick_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           143129710	# hungate asm on mick's bins
cat mick_vs_hungate_paf_filt.sort.merge.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
Total Length:           357991392	# micks bins on hungate asm

# Now to calculate total asm size for each one
perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";' < ~/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta.fai
1317539873	# Hungate total asm size
for i in /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/genomes/*.fa; do echo $i; sbatch --nodes=1 --ntasks-per-node=1 --mem=1000 -p short --wrap="samtools faidx $i"; done
cat /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/mick_rug/genomes/*.fai | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
2176657067	# Mick total asm size
perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); if($F[0] =~ /bin/){next;} if($F[-1] < 10){$c += $F[7];}} print "$c\n";' < ../dastool/illumina_megahit_dastool_DASTool_summary.txt
1475288770	# Illumina AN bin total asm size

# And my script to calculate the proportions for the network analysis:
perl calculate_network_for_hive.pl -f 1475288770 -s 2176657067 -t 1317539873 -j 411825315 -k 172100230 -l 342782266 -h 118637195 -y 357991392 -b 143129710 -o test.network.tab
```

```R
library("HiveR")
network <- read.delim("test.network.tab", header=FALSE)
colnames(network) <- c("source", "target", "weight")
network$weight <- as.numeric(network$weight)
hpd <- edge2HPD(network, axis.cols = c("red", "green", "blue"))
hpd$nodes[grep("first", hpd$nodes$lab), "axis"] = as.integer(1)
hpd$nodes[grep("second", hpd$nodes$lab), "axis"] = as.integer(2
hpd$nodes[grep("third", hpd$nodes$lab), "axis"] = as.integer(3)
hpd$edges$color[hpd$edges$weight == 2] <- "green"
hpd$edges$color[hpd$edges$weight == 1] <- "red"
hpd$nodes$radius <- as.numeric(sapply(str_split(hpd$nodes$lab, "_"), '[', 2))

# Inverting the nodes so that the matching proportion is on the outside:
hpd$nodes$radius <- 100 - as.numeric(sapply(str_split(hpd$nodes$lab, "_"), '[', 2)) + 1

plotHive(hpd, axLabs = c("Illumina", "Stewart et al.", "Hungate"), bkgnd = "black", dr.nodes = FALSE, ch = 0)
dev.copy2pdf(file="First_hive_plot_illumina.pdf", useDingbats=FALSE)
```

**OK, This didn't work out as well as I had hoped! Let's recalculate based on my master tables and generate an association matrix instead for a chord diagram**

## Chord diagram

The Chord diagram is going to be a bit messy, but let's see if we can calculate the proportions of the assembly better this way.

```bash
perl -e '<>; $m = 0; $h = 0; $j = 0; $t = 0; $a = 0; open(MI, "> illumina_ANbins_mickaligns.list"); open(HU, "> illumina_ANbins_hunaligns.list"); open(BO, "> illumina_ANbins_both.list"); while(<>){chomp; @s = split(/\t/); if($s[-2] eq "-" && $s[-1] eq "-"){$a += $s[1];} if($s[-2] ne "-" && $s[-1] ne "-"){$j += $s[1]; @b = split(/;/, $s[-2]); push(@b, split(/;/, $s[-1])); print BO join("\n", @b) . "\n";}elsif($s[-2] ne "-"){$m += $s[1]; @b = split(/;/, $s[-2]); print MI join("\n", @b) . "\n";}elsif($s[-1] ne "-"){$h += $s[1]; @b = split(/;/, $s[-1]); print HU join("\n", @b) . "\n";} $t += $s[1];} close MI; close HU; print "mick\thun\tboth\tsole\ttot\n"; print "$m\t$h\t$j\t$a\t$t\n";' < illumina_megahit_master_table_2018_09_07.ANbins.short.tab
mick    hun     both    sole    tot
392273034       84158897        328966430       669890409       1475288770

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug.fa.fai -c 0 -l illumina_ANbins_mickaligns.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
305642558	# Aligns only to MICK	Remaining: 1871014509 -  248619611 (both) =  1622394898
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta.fai -c 0 -l illumina_ANbins_hunaligns.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
665055395	# Aligns only to Hungate Remaining: 454890449

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta.fai -c 0 -l illumina_ANbins_both.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
759461029	# Both from Hungate
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug.fa.fai -c 0 -l illumina_ANbins_both.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
248619611	# Both from Mick


# Now for Pacbio
perl -e '<>; $m = 0; $h = 0; $j = 0; $t = 0; $a = 0; open(MI, "> pacbio_ANbins_mickaligns.list"); open(HU, "> pacbio_ANbins_hunaligns.list"); open(BO, "> pacbio_ANbins_both.list"); while(<>){chomp; @s = split(/\t/); if($s[-2] eq "-" && $s[-1] eq "-"){$a += $s[1];} if($s[-2] ne "-" && $s[-1] ne "-"){$j += $s[1]; @b = split(/;/, $s[-2]); push(@b, split(/;/, $s[-1])); print BO join("\n", @b) . "\n";}elsif($s[-2] ne "-"){$m += $s[1]; @b = split(/;/, $s[-2]); print MI join("\n", @b) . "\n";}elsif($s[-1] ne "-"){$h += $s[1]; @b = split(/;/, $s[-1]); print HU join("\n", @b) . "\n";} $t += $s[1];} close MI; close HU; print "mick\thun\tboth\tsole\ttot\n"; print "$m\t$h\t$j\t$a\t$t\n";' < pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab
mick    hun     both    sole    tot
208361045       28424376        424169426       137209071       798163918

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug.fa.fai -c 0 -l pacbio_ANbins_mickaligns.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
188401211	# Aligns only to MICK	Remaining: 1988255856 - (both) = 1766328688
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta.fai -c 0 -l pacbio_ANbins_hunaligns.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
357677323	# Aligns only to Hungate	Remaining: (after greping all nonredundant scaffolds): 594380222

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/hungate/hungate_combined_unordered_reference.fasta.fai -c 0 -l pacbio_ANbins_both.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
669683165	#Both from hungate
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ~/rumen_longread_metagenome_assembly/assemblies/mick_rug/mick_combined_900_rug.fa.fai -c 0 -l pacbio_ANbins_both.list | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); $c += $s[1];} print "$c\n";'
221927168	# both from MICK	
```

```R
library(chorddiag)
# Illumina chord graph
m <- matrix(c(669890409, 392273034, 84158897, 328966430, 305642558, 1622394898, 0, 248619611, 665055395, 0, 454890449, 197594029, 100000000, 100000000, 100000000, 0), byrow = TRUE, nrow=4, ncol=4)
assemblies <- c("Illumina", "Stewart et. al", "Hungate1000", "Both")
dimnames(m) <- list(origin = assemblies, destination = assemblies)
colors <- brewer.pal(4, "Dark2")
chorddiag(m, groupColors=colors, groupnamePadding = 20)
# I then opened in chrome, saved as a pdf and used Inkscape to grep out the diagram and modify it


#PacBio chord graph
m <- matrix(c(137209071, 208361045, 28424376, 424169426, 188401211, 1766328688, 0, 221927168, 357677323, 0, 594380222, 365482328, 100000000, 100000000, 100000000, 0), byrow=TRUE, nrow=4, ncol=4)
assemblies <- c("PacBio", "Stewart et. al", "Hungate1000", "Both")
dimnames(m) <- list(origin=assemblies, destination=assemblies)
colors <- brewer.pal(4, "Dark2")
chorddiag(m, groupColors=colors, groupnamePadding = 20)
```

## Generating supplementary tables

I want to subsection the master tables to generate our draft supplementary info. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

#### Supplementary tables 2 and 3

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,48,49,50,51]) . "\n";' < illumina_megahit_master_table_2018_09_07.tab > supplementary_table_2_shortread_bins.tab
# REDO:
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,48,49,50,51]) . "\n";' < illumina_megahit_master_table_2018_11_28.tab > supplementary_table_1_shortread_bins.tab
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,48,49,50,51]) . "\n";' < pacbio_final_pilon_master_table_2018_09_07.tab > supplementary_table_3_longread_bins.tab
```

#### Supplementary tables 4 and 5

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-" && $F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51]) . "\n";}' < illumina_megahit_master_table_2018_09_07.tab > supplementary_table_4_short_read_unique.tab
# REDO:
perl -ne 'chomp; @F = split(/\t/); if($F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-" && $F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51]) . "\n";}' < illumina_megahit_master_table_2018_11_28.tab > supplementary_table_6_short_read_unique.tab
perl -ne 'chomp; @F = split(/\t/); if($F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-" && $F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51]) . "\n";}' < pacbio_final_pilon_master_table_2018_09_07.tab > supplementary_table_5_long_read_unique.tab
```

#### Supplementary tables 6 and 7

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51,56,57,58]) . "\n";}' < illumina_megahit_master_table_2018_09_07.tab > supplementary_table_6_short_read_taxonomy.tab
# REDO:
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51,56,57,58]) . "\n";}' < illumina_megahit_master_table_2018_11_28.tab > supplementary_table_3_short_read_taxonomy.tab
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51,56,57,58]) . "\n";}' < pacbio_final_pilon_master_table_2018_09_07.tab > supplementary_table_7_long_read_taxonomy.tab
```

#### Supplementary tables 8 and 9

```bash
cat ../prodigal/illumina_megahit_prodigal_proteins.shortform.tab > supplementary_table_8_short_read_prodigal.tab
cat ../prodigal/pacbio_final_prodigal_proteins.shortform.tab > supplementary_table_9_long_read_prodigal.tab
```

## Testing hypotheses: partial proteins and assembled contigs

I want to test a few final hypotheses and generate a few more talking points for the manuscript.

#### Partial ORF dataset.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
module load samtools
cp grep_out_fasta_sequence.pl grep_out_partial_sequence.pl
sbatch grep_out_partial_sequence.pl illumina_megahit_prodigal_proteins.shortform.tab ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa illumina_megahit_prodigal_partial_proteins.fna
sbatch grep_out_partial_sequence.pl pacbio_final_prodigal_proteins.shortform.tab ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa pacbio_final_prodigal_partial_proteins.fna

perl -lane 'print $F[0];' < ../dastool/illumina_dastool_analysis_binset_lt10redund.bins > illumina_an_bin_contigs.list
sbatch grep_partial_prodigal_aaseq.pl illumina_megahit_prodigal_proteins.faa illumina_an_bin_contigs.list illumina_megahit_prodigal_partial_proteins.faa

perl -lane 'print $F[0];' < ../dastool/pacbio_dastool_analysis_binset_lt10redund.bins > pacbio_an_bin_contigs.list
sbatch grep_partial_prodigal_aaseq.pl pacbio_final_prodigal_proteins.faa pacbio_an_bin_contigs.list pacbio_final_prodigal_partial_proteins.faa

# Lets calculate the overall length stats for the dataset
perl -lane 'if($F[0] =~ /ContigID/){next;} $len = $F[2] - $F[1]; print "$F[5]\t$len";' < pacbio_final_prodigal_proteins.shortform.tab > pacbio_final_prodigal_proteins.lenXpart.tab
perl -lane 'if($F[0] =~ /ContigID/){next;} $len = $F[2] - $F[1]; print "$F[5]\t$len";' < illumina_megahit_prodigal_proteins.shortform.tab > illumina_megahit_prodigal_proteins.lenXpart.tab
```

```R
library(dplyr)
library(ggplot2)

pb <- read.delim("pacbio_final_prodigal_proteins.lenXpart.tab", header=FALSE)
il <- read.delim("illumina_megahit_prodigal_proteins.lenXpart.tab", header=FALSE)

colnames(pb) <- c("Partial", "Len")
colnames(il) <- c("Partial", "Len")
pb <- pb %>% mutate(Tech = c("PacBio"))
il <- il %>% mutate(Tech = c("Illumina"))
total <- bind_rows(pb, il)
total$Partial <- as.factor(total$Partial)
total <- total %>% mutate(outlier = Len > 10000)

pdf(file="partial_by_len_plots.pdf", useDingbats=FALSE)
ggplot(total, aes(x=Partial, y=Len, color=Partial)) + geom_violin(trim=FALSE) + geom_point(data = function(x) dplyr::filter_(x, ~ outlier), position = 'jitter') + theme_bw() + scale_color_brewer(palette="Dark2") + facet_grid(. ~ Tech)

```

#### Read overlaps vs GC percent

This is a test of the read depth bias. Basically, is there a huge tranche of low GC reads from the PacBio end that are just not getting assembled?

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/pacbio_aligns

```bash
module load minimap2
sbatch --nodes=1 --ntasks-per-node=5 --mem=16000 -p short --wrap="minimap2 -t 5 -x map-pb ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa ~/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta > pacbioasm_vs_errcorrected_reads.paf"

sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnCounter.py -f pacbioasm_vs_errcorrected_reads.paf -d '\t' -c 0 > pacbioasm_vs_errcorrected_reads.align.counts"

perl -lane 'print "$F[0]\t$F[2]";' < rumen_pacbio_corrected.gc.tab > rumen_pacbio_corrected.gc.corrected.tab
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); <IN>; print "overlaps\tgc\n"; %d; while(<IN>){chomp; @s = split(/\t/); $d{$s[0]} = 1; if(exists($h{$s[0]})){ print "$s[1]\t$h{$s[0]}\n";}} close IN; foreach my $k (keys(%h)){if(!exists($d{$k})){print "0\t$h{$k}\n";}}' rumen_pacbio_corrected.gc.corrected.tab pacbioasm_vs_errcorrected_reads.align.counts > pacbio_ovlp_by_gc.tab

wc -l pacbioasm_vs_errcorrected_reads.align.counts
4405455 pacbioasm_vs_errcorrected_reads.align.counts
wc -l pacbio_ovlp_by_gc.tab
5931683 pacbio_ovlp_by_gc.tab
# There are allot of missing reads with no alignment overlap!

# I'm curious: how does it look when I align pacbio reads to Illumina contigs?
sbatch --nodes=1 --ntasks-per-node=5 --mem=16000 -p short --wrap="minimap2 -t 5 -x map-pb ~/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa ~/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta > illuminaasm_vs_errcorrected_reads.paf"

sbatch --nodes=1 --mem=20000 --ntasks-per-node=3 -p short --wrap="python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/utils/tabFileColumnCounter.py -f illuminaasm_vs_errcorrected_reads.paf -c 0 -d '\t' > illuminaasm_vs_errcorrected_reads.align.counts"
```

Now to try to plot it all out in R

```R
library(dplyr)
library(ggplot2)
library(MASS)
library(viridis)

data <- read.delim("pacbio_ovlp_by_gc.tab")
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data$density <- get_density(data$overlaps, data$gc)

pdf(file="pacbio_gc_by_aln_overlaps.pdf", useDingbats=FALSE)
ggplot(data) + geom_point(aes(gc, overlaps, color=density)) + scale_color_viridis()
dev.off()

bitmap(file = "pacbio_gc_by_aln_overlaps.png", res=300, type="jpeg")

# Now for a violin plot
data$ovcat <- ifelse(data$overlaps > 0, "Overlaps", "NoOverlap")
pdf(file="pacbio_overlaps_gc_violin.pdf", useDingbats=FALSE)
ggplot(data, aes(x=ovcat, y=gc, color=ovcat)) + geom_violin(trim=FALSE) + theme_bw() + scale_colour_brewer(palette="Dark2")
dev.off()

# Now to generate tables
data %>% group_by(ovcat) %>% summarize('25%'= quantile(gc, probs=0.25), '50%' = quantile(gc, probs=0.50), '75%' = quantile(gc, probs=0.75), avg=# A tibble: 2 x 7= median(gc), n=n())
      ovcat     `25%`     `50%`     `75%`       avg    median       n
      <chr>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>   <int>
1 NoOverlap 0.2927553 0.4047442 0.4960352 0.4013871 0.4047442 1526228
2  Overlaps 0.4100586 0.4850455 0.5349027 0.4623783 0.4850455 4405454
```


#### Validation of bins pre- and post-DAS_Tool

I just want to run checkm on the pre-and post-DAS_Tool bins to estimate stats using a second tool. Then we compare the rankings to see which binning method was best, overall.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool

```bash
module load checkm
module load prodigalorffinder

sbatch --nodes=1 --mem=240000 --ntasks-per-node=40 -p mem --wrap="checkm lineage_wf -t 20 -x fa illumina_megahit_dastool_DASTool_bins illumina_megahit_dastool_checkm"

sbatch --nodes=1 --mem=100000 --ntasks-per-node=20 -p short --wrap="checkm lineage_wf -t 20 -x fa pacbio_final_dastool_DASTool_bins pacbio_final_dastool_checkm"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p short --wrap="checkm qa --tab_table -o 1 -f pacbio_final_dastool_checkm.results.tab pacbio_final_dastool_checkm/lineage.ms pacbio_final_dastool_checkm"
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
sbatch --nodes=1 --mem=100000 --ntasks-per-node=20 -p short --wrap="checkm lineage_wf -t 20 -x fa public_metabat illumina_megahit_checkm"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p short --wrap="checkm qa --tab_table -o 1 -f illumina_megahit_checkm.results.tab illumina_megahit_checkm/lineage.ms illumina_megahit_checkm"
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon

```bash
sbatch --nodes=1 --mem=100000 --ntasks-per-node=20 -p short --wrap="checkm lineage_wf -t 20 -x fa metabat2 pacbio_megahit_checkm"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p short --wrap="checkm qa --tab_table -o 1 -f pacbio_megahit_checkm.results.tab pacbio_megahit_checkm/lineage.ms pacbio_megahit_checkm"

sbatch --nodes=1 --mem=100000 --ntasks-per-node=20 -p short --wrap="checkm lineage_wf -t 20 -x fasta hic_clusters/best_genome_clusters pacbio_hic_checkm"
sbatch --nodes=1 --mem=12000 --ntasks-per-node=2 -p short --wrap="checkm qa --tab_table -o 1 -f pacbio_hic_checkm.results.tab pacbio_hic_checkm/lineage.ms pacbio_hic_checkm"
```

#### Counting the effects of DAS_Tool dereplication

I wrote a brute-force association script to try to determine how frequently the DASbins break up the original bins after dereplication.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool

```bash
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p short --wrap="perl backtrace_bin_affiliation.pl illumina_megahit_public_metabat.unsorted.bins illumina_megahit_hic.unsorted.bins illumina_dastool_analysis_binset_lt10redund.bins illumina_dastool_backtrace"

sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p short --wrap="perl backtrace_bin_affiliation.pl pacbio_final_public_metabat.unsorted.bins pacbio_final_public_hic.unsorted.bins pacbio_dastool_analysis_binset_lt10redund.bins pacbio_dastool_backtrace"

# And now a test to see just how frequently the bins change stats
perl compare_bin_stats.pl pacbio_dastool_backtrace.key ../../assemblies/pilot_project/pacbio_final_pilon/pacbio_hic_checkm.results.tab ../../assemblies/pilot_project/pacbio_final_pilon/pacbio_megahit_checkm.results.tab pacbio_final_dastool_checkm.results.tab pacbio_final_dastool_checkm.hic_mega_comp.tab
```

```R
library(dplyr)
library(ggplot2)

data <- read.delim("pacbio_final_dastool_checkm.hic_mega_comp.tab")

cor(as.matrix(data[,c(4:9)]))
            DasComp   PreComp    DasCon    PreCon DasStrain PreStrain
DasComp   1.0000000 0.7977572 0.5420134 0.2479558 1.0000000 0.3322989
PreComp   0.7977572 1.0000000 0.4458649 0.4269740 0.7977572 0.4307700
DasCon    0.5420134 0.4458649 1.0000000 0.5183714 0.5420134 0.3356574
PreCon    0.2479558 0.4269740 0.5183714 1.0000000 0.2479558 0.2802857
DasStrain 1.0000000 0.7977572 0.5420134 0.2479558 1.0000000 0.3322989
PreStrain 0.3322989 0.4307700 0.3356574 0.2802857 0.3322989 1.0000000

pdf(file="dastool_pre_post_completion.pdf", useDingbats=FALSE)
ggplot(data, aes(x = DasComp, y =PreComp, fill = Tech)) + geom_point(size=3) + geom_smooth(method="lm") + theme_bw() + xlab(label="Das_Tool bin completeness") + ylab(label="Pre-Das bin completeness") + facet_wrap( ~ Tech)

pdf(file="dastool_pre_post_contamination.pdf", useDingbats=FALSE)
ggplot(data, aes(x = DasCon, y =PreCon, fill = Tech)) + geom_point(size=3) + geom_smooth(method="lm") + theme_bw() + xlab(label="Das_Tool bin contamination") + ylab(label="Pre-Das bin contamination") + facet_wrap( ~ Tech)

# Now for the counts of changed bins
count <- read.delim("pacbio_dastool_backtrace.count")
pdf(file="dastool_contig_counts_per_bin.pdf", useDingbats=FALSE)
ggplot(count, aes(x = DasBinCtgNum, y =OrigBinCtgNum, fill = Tech)) + geom_point(size=3) + geom_smooth(method="lm") + theme_bw() + xlab(label="Das_Tool contigs per bin") + ylab(label="Pre-Das contigs per bin") + facet_wrap( ~ Tech)

pdf(file="dastool_contig_counts_violon.pdf", useDingbats=FALSE)
ggplot(count, aes(x = Tech, y = Diff, color=Tech)) + geom_violin(trim=FALSE) + ylab(label="Difference in Orig contigs DAS_Tool")


data <- data %>% mutate(DiffCon = PreCon - DasCon)

```
#### Plotting taxonomic information based on VIR and AMR gene hi-c links

I want to try to replot Max's heatmaps but with the superkingdom affiliation of each bin.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/amr_vir_heatmaps

```bash
# Creating tax tables for the association
for i in *.tsv; do dos2unix $i; done
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[10]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < ../master_tables/illumina_megahit_master_table_2018_09_07.ANbins.short.tab |perl -ne 'chomp; @F = split(/\t/); my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d) . "\n";' > illumina_megahit_master_table_09_07.hicbin.taxtable.tab

perl -ne 'chomp; @F = split(/\t/); my @values; push(@values, $F[0]); for($x = 1; $x < scalar(@F); $x++){my %h; @bsegs = split(/;/, $F[$x]); foreach my $row (@bsegs){@hsegs = split(/:/, $row); $h{$hsegs[0]} = $hsegs[1];} @sorted = sort{$h{$b} <=> $h{$a}} keys(%h); push(@values, $sorted[0] . "");} print join("\t", @values); print "\n";' < illumina_megahit_master_table_09_07.hicbin.taxtable.tab > illumina_megahit_master_table_09_07.hicbin.taxconsensus.tab 


# Pacbio
perl -e '%h; while(<>){chomp; @s = split(/\t/); push(@{$h{$s[10]}}, [$s[4], $s[5], $s[6], $s[7], $s[8], $s[9]]);} foreach my $bin (keys(%h)){%t = (); foreach my $row (@{$h{$bin}}){for($x = 0; $x < scalar(@{$row}); $x++){push(@{$t{$x}}, $row->[$x]);}} @j = (); foreach $idx (sort {$a <=> $b}keys(%t)){push(@j, join(";", @{$t{$idx}}));} print "$bin\t" . join("\t", @j) . "\n";}' < ../master_tables/pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab | perl -ne 'chomp; @F = split(/\t/); my @d; for($x = 1; $x < scalar(@F); $x++){my %h; @jsegs = split(/;/, $F[$x]); foreach $k (@jsegs){$h{$k} += 1;} my @tmp; foreach $k (sort{$a cmp $b}keys(%h)){push(@tmp,"$k:$h{$k}");} push(@d, join(";", @tmp));} print "$F[0]\t" . join("\t", @d) . "\n";' > pacbio_final_pilon_master_table_2018_09_07.hicbin.taxtable.tab
perl -ne 'chomp; @F = split(/\t/); my @values; push(@values, $F[0]); for($x = 1; $x < scalar(@F); $x++){my %h; @bsegs = split(/;/, $F[$x]); foreach my $row (@bsegs){@hsegs = split(/:/, $row); $h{$hsegs[0]} = $hsegs[1];} @sorted = sort{$h{$b} <=> $h{$a}} keys(%h); push(@values, $sorted[0] . "");} print join("\t", @values); print "\n";' < pacbio_final_pilon_master_table_2018_09_07.hicbin.taxtable.tab > pacbio_final_pilon_master_table_2018_09_07.hicbin.taxconsensus.tab


perl -e '$h = <>; print "bin\t$h"; while(<>){chomp; @s = split(/\t/); @csegs = split(/\./, $s[0]); $csegs[1] = $csegs[1] * 1; $s[0] = $csegs[1]; print join("\t", @s); print "\n";}' < ilmn_arg_hic_links_filt.tsv > ilmn_arg_hic_links_filt.mod.tsv
perl -e '$h = <>; print "bin\t$h"; while(<>){chomp; @s = split(/\t/); @csegs = split(/\./, $s[0]); $csegs[1] = $csegs[1] * 1; $s[0] = $csegs[1]; print join("\t", @s); print "\n";}' < pb_vir_hic_links_filt.tsv > pb_vir_hic_links_filt.mod.tsv
perl -e '$h = <>; print "bin\t$h"; while(<>){chomp; @s = split(/\t/); @csegs = split(/\./, $s[0]); $csegs[1] = $csegs[1] * 1; $s[0] = $csegs[1]; print join("\t", @s); print "\n";}' < pb_arg_hic_links_filt.tsv > pb_arg_hic_links_filt.mod.tsv


# Now for the virus stats
python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f ../blobtools/pacbio_pilon_viruses.host.assoc.tab -c 2 -d "\t"
Entry   Value
Phikzvirus      49
N4virus 32
Bpp1virus       23
Cp8virus        8
C5virus 4
T5virus 3
Cjw1virus       2
Schizot4virus   1
Spbetavirus     1
Sk1virus        1

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); <IN>; %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[1]} = "NA";} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$h{$s[0]} = $s[1];}} close IN; foreach my $k (keys(%h)){print "$k\t$h{$k}\n";}' ../blobtools/pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab ../dastool/pacbio_final_public_hic.unsorted.bins > pacbio_hostcontig_assoc.tab
grep -v NA pacbio_hostcontig_assoc.tab | grep -v NOBIN > pacbio_hostcontig_assoc.filt.tab
```

First let's try the AMR genes.

```R
library(dplyr)
library(RColorBrewer)
arg.ilhic <- read.delim("ilmn_arg_hic_links_filt.mod.tsv", header=TRUE)
tax.ilhic <- read.delim("illumina_megahit_master_table_09_07.hicbin.taxconsensus.tab", header=FALSE)
colnames(tax.ilhic) <- c("bin", "Kingdom", "Phylum", "Class", "Family", "Genus", "FAPRO")
arg.ilhic$bin <- as.factor(arg.ilhic$bin)

comp.ilhic <- inner_join(arg.ilhic, tax.ilhic, by="bin")
comp.ilhic <- mutate(comp.ilhic, Name = paste0(bin, "_", Kingdom, "_", Genus))
rownames(comp.ilhic) <- comp.ilhic$Name
my.col <- brewer.pal(3, "Set1")[comp.ilhic$Kingdom]

pdf(file="ilmn_hic_arg_heatmap.pdf", useDingbats=FALSE)
heatmap(as.matrix(comp.ilhic[,2:12]), cexRow=1.5, labRow=paste(comp.ilhic$Kingdom, comp.ilhic$Genus, sep=" "), Colv = NA, Rowv = NA, RowSideColors =my.col)
dev.off()


# Now for the Pacbio data
arg.pbhic <- read.delim("pb_arg_hic_links_filt.mod.tsv", header=TRUE)
tax.pbhic <- read.delim("pacbio_final_pilon_master_table_2018_09_07.hicbin.taxconsensus.tab", header=FALSE)
colnames(tax.pbhic) <- c("bin", "Kingdom", "Phylum", "Class", "Family", "Genus", "FAPRO")
arg.pbhic$bin <- as.factor(arg.pbhic$bin)
tax.pbhic$bin <- as.factor(tax.pbhic$bin)

comp.pbhic <- inner_join(arg.pbhic, tax.pbhic, by="bin")
my.col <- brewer.pal(3, "Set1")[comp.pbhic$Kingdom]

pdf(file="pb_hic_arg_heatmap.pdf", useDingbats=FALSE)
heatmap(as.matrix(comp.pbhic[,2:139]), cexRow=1.5, labRow=paste(comp.pbhic$Kingdom, comp.pbhic$Genus, sep=" "), Colv=NA, Rowv=NA, RowSideColors=my.col)
dev.off()

# It was pretty ugly and needs some clustering to look ok. I'm going to do a hierarchical clustering first and reorder the data
rownames(comp.pbhic) <- comp.pbhic$bin
pbhic.arg.clust <- hclust(dist(comp.pbhic[2:139]))
pbhic.arg.clust.r <- hclust(dist(t(comp.pbhic[2:139])))

comp.pbhic.order <- comp.pbhic[pbhic.arg.clust$order, ]
my.col <- brewer.pal(3, "Set1")[comp.pbhic.order$Kingdom]

pdf(file="pb_hic_arg_heatmap.pdf", useDingbats=FALSE)
heatmap(as.matrix(comp.pbhic.order[,2:139]), cexRow=1.5, labRow=paste(comp.pbhic.order$Kingdom, comp.pbhic.order$Genus, sep=" "), Colv=NA, Rowv=NA, RowSideColors=my.col)
dev.off()
```

Now the viruses

```R
library(RColorBrewer)
library(dplyr)
library(gplots)

pbvir.names <- read.delim("../blobtools/pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab", header=TRUE)
pbvir.names.f <- unique(pbvir.names[,c(1,4)])

pbvir.links <- read.delim("pb_vir_hic_links_filt.mod.tsv", header=TRUE)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pbvir.aligns <- read.delim("pacbio_hostcontig_assoc.filt.tab", header=FALSE)

pbvir.names.sort <- pbvir.names.f[pbvir.names.f$VirusCtg %in% colnames(pbvir.links),]
cols <- palette(brewer.pal(8, "Dark2"))[pbvir.names.sort$VirusGenus]

pdf("pacbio_virus_hic_assoc_heatmap.pdf", useDingbats=FALSE)
heatmap.2(as.matrix(pbvir.links[,2:48]), labCol=pbvir.names.sort$VirusGenus, labRow=pbvir.links$bin, trace="none", ColSideColors=cols, col=hmcol)
dev.off()
```

#### Full length 16S

I need to pull out the 16S (> 1500 bp) that are full length from the pacbio dataset and then see if their categorization fits what the general consensus of the bin is from Blobtools annotation.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
perl -e '$h = <>; while(<>){chomp; @s = split(/\t/); if($s[6] < 1500){next;}else{$s[1] =~ s/\_pilon.*$//; print join("\t", @s); print "\n";}}' < mick_16s_pacbio_summary.tsv > mick_16s_pacbio_summary.reformat.tab

perl -lane 'print $F[1]' < mick_16s_pacbio_summary.reformat.tab > mick_16s_pacbio_summary.contig.list

python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab -c 0 -l mick_16s_pacbio_summary.contig.list | perl -ne 'chomp; @F = split(/\t/); print "$F[11]\t$F[0]\t$F[4]\t$F[10]\t$F[12]\n";' > mick_16s_pacbio_summary.reformat.mastertable.tab

python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f mick_16s_pacbio_summary.reformat.mastertable.tab -f ../amr_vir_heatmaps/pacbio_final_pilon_master_table_2018_09_07.hicbin.taxconsensus.tab -c 0 -o test
perl -lane '$F[2] = lc($F[2]); $F[4] = lc($F[4]); print "$F[0]\t$F[2]\t$F[4]\t$F[9]";' < merged_mick_ssu_collection.tab > merged_mick_ssu_collection.tax.tab

perl -lane 'if($F[3] eq $F[2] || $F[3] eq $F[1]){print "Match";}else{print "Nope";}' < merged_mick_ssu_collection.tax.tab | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0
Entry   Value
Match   176
Nope    5

```

## Rerunning cluster analysis

Now that I have the full set of proximeta, Illumina bins. I need to rerun das_tool and regenerate the statistics for publication.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated/

```bash
export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin/:$PATH

# First, generate the bin assignments using the fasta files from Max
perl -e '@f = `ls ../final_clusters/*.fasta`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-2]\n";}} close IN;}' > illumina_megahit_hic.final.bins

# now to run DAS_tool
# I'm going to run it with protein prediction just to see if it changes anything
/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i illumina_megahit_hic.final.bins,illumina_megahit_public_metabat.unsorted.bins -c mick_megahit_final_full.rfmt.fa -o illumina_final_dastool -l HiC,metabat --search_engine diamond -t 20 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db --write_bins 1 --score_threshold 0

tar -czvf illumina_final_dastool_DASTool_bins.tar.gz illumina_final_dastool_DASTool_bins

perl -e 'use File::Basename; @f = `ls illumina_megahit_dastool_DASTool_bins/*.fa`; chomp(@f); foreach $h (@f){@hsegs = split(/[\.\_]/, basename($h)); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-4]\_$hsegs[-3]\n";}} close IN;}' > illumina_final_dastool_DASTool_bins.tab

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 10){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_analysis_binset_lt10redund.bins

perl -lane '@b = split(/[_\.]/, $F[1]); print "$F[0]\t$b[-2]\_$b[-1]";' < illumina_final_dastool_analysis_binset_lt10redund.bins > illumina_final_dastool_analysis_binset_lt10redund.tab

perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] > 5 || $s[-2] < 80){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; print "$l\t$bsegs[-2].$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_HQbins.bins


```

Running checkm to generate summary stats independent of DAS_Tool. Also generating the appropriate master tables and other stats that I need for the manuscript.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool

```bash
module load checkm prodigalorffinder
sbatch --nodes=1 --mem=300000 --ntasks-per-node=20 -p mem --wrap="checkm lineage_wf -t 20 -x fa illumina_final_dastool_DASTool_bins illumina_final_dastool_checkm"
```



## Bin testing for significance

```R
library(dplyr)
library(ggplot2)
library(MASS)
ctgtabs <- read.delim("illumina_megahit_master_table_2018_11.contig.bin.stats.tab")

ctgtabs <- mutate(ctgtabs, MetabatBin = ifelse(Metabat == "NOBIN", "NONE", "BIN"), HiCBin = ifelse(HiC == "NOBIN", "NONE", "BIN"))
ctgtabs$MetabatBin <- as.factor(ctgtabs$MetabatBin)
ctgtabs$HiCBin <- as.factor(ctgtabs$HiCBin)
table(ctgtabs$len, ctgtabs$MetabatBin)

ctgtabs <- mutate(ctgtabs, Length = ifelse(len <= 2500, "Small", "Large"), GCProp = ifelse(GC < 0.42, "Low", "High"))
ctgtabs$Length <- as.factor(ctgtabs$Length)
ctgtabs$GCProp <- as.factor(ctgtabs$GCProp)

ctgtabs <- mutate(ctgtabs, BinStatus = ifelse(MetabatBin == "BIN" | HiCBin == "BIN", ifelse(MetabatBin == "BIN" & HiCBin == "BIN", "Both", ifelse(MetabatBin == "BIN", "META", "HIC")), "NONE"))
ctgtabs$BinStatus <- as.factor(ctgtabs$BinStatus)

bitmap("illumina_gc_vs_len_binning.png", res=300)
ggplot(ctgtabs, aes(x=len, y=GC, color=BinStatus)) + geom_point() + theme_bw() + stat_ellipse() + scale_x_log10() + xlab("Log10 Contig Length (bp)") + ylab("Contig GC proportion")
dev.off()

pdf("illumina_gc_vs_len_binning.pdf", useDingbats=FALSE)
ggplot(ctgtabs, aes(x=len, y=GC, color=BinStatus)) + stat_ellipse() + theme_bw() + xlab("Contig Length (bp)") + ylab("Contig GC proportion")
dev.off()

as.matrix(table(ctgtabs$BinStatus, ctgtabs$Length))

        Large  Small
  Both 185960  46064
  HIC  163840 685213
  META  91295  63848
  NONE  73273 872770

as.matrix(table(ctgtabs$BinStatus, ctgtabs$GCProp))

         High    Low
  Both 142407  89617
  HIC  598866 250187
  META  29891 125252
  NONE 320646 625397


ctgtabs.pb <- read.delim("pacbio_final_pilon_master_table_2018_09.contig.bin.stats.tab")
ctgtabs.pb <- mutate(ctgtabs.pb, MetabatBin = ifelse(Metabat == "NOBIN", "NONE", "BIN"), HiCBin = ifelse(HiC == "NOBIN", "NONE", "BIN"))
ctgtabs.pb$MetabatBin <- as.factor(ctgtabs.pb$MetabatBin)
ctgtabs.pb$HiCBin <- as.factor(ctgtabs.pb$HiCBin)

ctgtabs.pb <- mutate(ctgtabs.pb, Length = ifelse(len <= 10000, "Small", "Large"), GCProp = ifelse(GC < 0.48, "Low", "High"))
ctgtabs.pb$Length <- as.factor(ctgtabs.pb$Length)
ctgtabs.pb$GCProp <- as.factor(ctgtabs.pb$GCProp)

ctgtabs.pb <- mutate(ctgtabs.pb, BinStatus = ifelse(MetabatBin == "BIN" | HiCBin == "BIN", ifelse(MetabatBin == "BIN" & HiCBin == "BIN", "Both", ifelse(MetabatBin == "BIN", "META", "HIC")), "NONE"))
ctgtabs.pb$BinStatus <- as.factor(ctgtabs.pb$BinStatus)

table(ctgtabs.pb$BinStatus, ctgtabs.pb$Length)

       Large Small
  Both 23601 13942
  HIC  10530 14998
  META  1960  3671
  NONE  2328  6640

table(ctgtabs.pb$BinStatus, ctgtabs.pb$GCProp)

        High   Low
  Both 23622 13921
  HIC  14496 11032
  META  1997  3634
  NONE  2916  6052
```

## Solden Viral contigs and resurrecting my network plot

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/amr_vir_heatmaps

```bash
cat */*.fna > solden_viral_contigs.fa

grep '>' solden_viral_contigs.fa | wc -l
1497

sbatch --nodes=1 --ntasks-per-node=6 --mem=40000 -p short --wrap="minimap2 -x asm5 -t 6 solden_viral_contigs.fa /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > solden_vs_pacbio_asm5_align.paf"
sbatch --nodes=1 --ntasks-per-node=6 --mem=40000 -p short --wrap="minimap2 -x asm5 -t 6 solden_viral_contigs.fa /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.fa > solden_vs_illumina_asm5_align.paf"

#OK let's filter the paf file
cat solden_vs_pacbio_asm5_align.paf | cut -f11 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   11843
Sum:    16275828
Minimum 40
Maximum 62031
Average 1374.299417
Median  714
Standard Deviation      2063.410233
Mode(Highest Distributed Value) 49

# Looks like 1kb would be a good cutoff
perl -lane 'if($F[10] > 1000 && $F[11] > 0){print $_;}' < solden_vs_pacbio_asm5_align.paf | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 5 -d '\t' | wc -l
200 <- viruses from their data that match ours

perl -lane 'if($F[10] > 1000 && $F[11] > 0){print "$F[0]";}' < solden_vs_pacbio_asm5_align.paf | sort | uniq > solden_vs_pacbio_asm5_align.pbuniqmaps.list
wc -l solden_vs_pacbio_asm5_align.pbuniqmaps.list
1600 solden_vs_pacbio_asm5_align.pbuniqmaps.list 	<- perhaps there are prophage contaminants or other ORFs in these contigs that are found in their viruses?

cat solden_vs_illumina_asm5_align.paf | cut -f11 | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   31430
Sum:    27030078
Minimum 40
Maximum 31453
Average 860.008845
Median  476
Standard Deviation      1272.475682
Mode(Highest Distributed Value) 40

# Again 1kb and no zero mapq maps
perl -lane 'if($F[10] > 1000 && $F[11] > 0){print $_;}' < solden_vs_illumina_asm5_align.paf | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 5 -d '\t' | wc -l
381	<- viruses from their dataset that match the illumina data
perl -lane 'if($F[10] > 1000 && $F[11] > 0){print $F[0];}' < solden_vs_illumina_asm5_align.paf > iluniqmaps.list
wc -l iluniqmaps.list
3915 iluniqmaps.list	<- again, prophage contaminants?

# getting the hic links
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f new_rumen_pacbio.counts -l pacbio_pilon_viruses.list -c 0 > pacbio_pilon_viruses.hiclinks.tab
perl -lane 'if($F[2] > 20){print $_;}' < pacbio_pilon_viruses.hiclinks.tab > pacbio_pilon_viruses.hiclinks.filt.tab

perl -lane 'print "$F[0]";' < ../blobtools/pacbio_pilon_viruses.fa.fai > pacbio_pilon_viruses.list
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl solden_vs_pacbio_asm5_align.pbuniqmaps.list pacbio_pilon_viruses.list
File Number 1: solden_vs_pacbio_asm5_align.pbuniqmaps.list
File Number 2: pacbio_pilon_viruses.list
Set     Count
1       1600
2       116

# No overlap??
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_2018_09_07.ANbins.short.tab -l solden_vs_pacbio_asm5_align.pbuniqmaps.list -c 0 > solden_vs_pacbio_asm5_align.pbuniqmaps.master.tab


# I'm really confused -- there's little overlap between their contigs and ours. The overlaps that I do find tend to map to contigs in HQ bins and the genes in the contigs are defined in eggnogmapper
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f new_rumen_pacbio.counts -l solden_vs_pacbio_asm5_align.pbuniqmaps.list -c 0 > solden_pacbio_hiclinks.tab
perl -lane 'if($F[2] > 20){print $_;}' < solden_pacbio_hiclinks.tab > solden_pacbio_hiclinks.filt.tab

# I'm giving up on this crap. Our methods aren't compatible at all
perl generateViralAssociationGraph.test.pl ../blobtools/pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt pacbio_pilon_viruses.hiclinks.filt.tab pacbio_pilon_viruses.hiclinks.filt.cyto.tab

perl -lane 'if($F[2] > 1){print $_;}' < ../blobtools/pacbio_pilon_viruses_ecpbreads.assoc.filt.stringent.cyto.tab > pacbio_pilon_viruses_ecpbreads.precombine.cyto.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]}->{$s[1]} = [$s[3], $s[4], $s[5]];} close IN; print "VirusCtg\tHostCtg\tCategory\tVirusGenus\tHostKingdom\tHostGenus\n"; open(IN, "< $ARGV[1]"); %seen; while(<IN>){chomp; @s = split(/\t/); if(exists($data{$s[0]}->{$s[1]})){print "$s[0]\t$s[1]\tBOTH\t$s[3]\t$s[4]\t$s[5]\n";}else{print "$s[0]\t$s[1]\tPACB\t$s[3]\t$s[4]\t$s[5]\n";} $seen{$s[0]}->{$s[1]} = 1;} foreach my $v (keys(%data)){foreach my $c (keys(%{$data{$v}})){ if(!exists($seen{$v}->{$c})){print "$v\t$c\tHIC\t" . join("\t", @{$data{$v}->{$c}}) . "\n";}}}' pacbio_pilon_viruses.hiclinks.filt.cyto.tab pacbio_pilon_viruses_ecpbreads.precombine.cyto.tab > pacbio_pilon_viruses.combined.cyto.tab

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f pacbio_pilon_viruses.combined.cyto.tab -c 2 -d '\t'
Entry   Value
PACB    105
HIC     64
BOTH    19


## Now to repeat the same thing for the Illumina data
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f new_rumen_illumina.counts -l illumina_megahit_viruses.list -c 0 > illumina_viruses.hiclinks.tab
perl -lane 'if($F[0] ne $F[1] && $F[2] > 10){print $_;}' < illumina_viruses.hiclinks.tab > illumina_viruses.hiclinks.filt.tab
perl -lane 'if($F[2] > 1){print $_;}' < ../blobtools/illumina_megahit_viruses_ecpbreads.assoc.filt.stringent.cyto.tab > illumina_megahit_viruses_ecpbreads.precombine.cyto.tab

perl generateViralAssociationGraph.test.pl ../blobtools/illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt illumina_viruses.hiclinks.filt.tab illumina_viruses.hiclinks.filt.cyto.tab

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f illumina_megahit_viruses.combined.cyto.tab -c 2 -d '\t'
Entry   Value
PACB    67
HIC     36
BOTH    6
```

## Plotting the CheckM stats

This is for **Fig 2c + d**. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool


```R
library(doMC)
library(data.table)
library(ggplot2)
library(dplyr)

pbmethods <- list(HiC=read.delim("../../assemblies/pilot_project/pacbio_final_pilon/pacbio_hic_checkm.results.tab", sep="\t", header=TRUE), MetaBat=read.delim("../../assemblies/pilot_project/pacbio_final_pilon/pacbio_megahit_checkm.results.tab", sep="\t", header=TRUE), DASTool=read.delim("pacbio_final_dastool_checkm.results.tab", sep="\t", header=TRUE))

methods <- names(pbmethods)
for(i in 1:length(pbmethods)){
pbmethods[[i]]$Method <- methods[i]
}
result_table <- do.call(rbind.data.frame, pbmethods)
tmp_wide <- result_table %>% group_by(Method) %>% filter(Contamination < 5) %>% summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90), `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70))

melt(tmp_wide,id.vars = 'Method', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table <- melt(tmp_wide,id.vars = 'Method', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table$title <- "CheckM Bin Statistics (< 5% Contamination)"
packageVersion("ggplot2")
plot_table$Completeness <- factor(plot_table$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%')))
colors <- rev(c("#08306B","#1664AB","#4A97C9","#93C4DE"))
pdf("check_dastool_stats_pacbio_plot_cont5.pdf", useDingbats=FALSE)
ggplot(plot_table, aes(Method, Bins, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)
dev.off()


ilmethods <- list(HiC=read.delim("illumina_final_hic_checkm.tab", sep="\t", header=TRUE), MetaBat=read.delim("../../assemblies/pilot_project/illumina_megahit/illumina_megahit_checkm.results.tab", sep="\t", header=TRUE), DASTool=read.delim("illumina_final_dastool_checkm.results.tab", sep="\t", header=TRUE))

ilhicshift <- ilmethods$HiC[,c("cluster_id", "Completeness", "Contamination")]
colnames(ilhicshift) <- c("Bin.Id", "Completeness", "Contamination")
ilmetashift <- ilmethods$MetaBat[,c("Bin.Id", "Completeness", "Contamination")]
ildasshift <- ilmethods$DASTool[,c("Bin.Id", "Completeness", "Contamination")]
ilmethods <- list(HiC=ilhicshift, MetaBat=ilmetashift, DASTool=ildasshift)

methods <- names(ilmethods)
for(i in 1:length(ilmethods)){
ilmethods[[i]]$Method <- methods[i]
}
result_table <- do.call(rbind.data.frame, ilmethods)
tmp_wide <- result_table %>% group_by(Method) %>% filter(Contamination < 5) %>% summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90), `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70))

melt(tmp_wide,id.vars = 'Method', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table <- melt(tmp_wide,id.vars = 'Method', measure.vars = c('>90%','>80%','>70%', '>60%'), value.name = 'Bins',variable.name ="Completeness")
plot_table$title <- "CheckM Bin Statistics (< 5% Contamination)"

plot_table$Completeness <- factor(plot_table$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%')))
colors <- rev(c("#08306B","#1664AB","#4A97C9","#93C4DE"))
pdf("check_dastool_stats_illumina_plot_cont5.pdf", useDingbats=FALSE)
ggplot(plot_table, aes(Method, Bins, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)
dev.off()
```

## Replotting the DAS_tool eval stats

This is for **Fig 2e**

> pwd: F:/SharedFolders/metagenomics/pilot_manuscript/figure_drafts/das_tool/

```R
library(dplyr)
library(ggplot2)
library(gridExtra)

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

illumina.hic <- read.delim("illumina_final_dastool_HiC.eval")
illumina.metabat <- read.delim("illumina_final_dastool_metabat.eval")
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

scatterPlot <- ggplot(combined, aes(x=SCG_completeness, y=SCG_redundancy, color=Bin)) + geom_point(aes(shape=Tech), position="jitter") + scale_color_brewer(palette="Paired") + theme(legend.position=c(0,1), legend.justification=c(0,1))
xdensity <- ggplot(combined, aes(x=SCG_completeness, fill=Bin)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Paired") + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank())
ydensity <- ggplot(combined, aes(x=SCG_redundancy, fill=Bin)) + geom_density(alpha=0.5) + scale_fill_brewer(palette="Paired") + theme(legend.position = "none", axis.title.y=element_blank(), axis.text.y=element_blank()) + coord_flip()
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4,1.4), heights = c(1.4,4))
```

## ARG allele hic analysis

I just want to glean the important features of the Hi-C ARG analysis that Max prepared.

> pwd: /cygdrive/f/SharedFolders/metagenomics/pilot_manuscript/figure_drafts/amr_genes

```bash
head -n 1 pb_arg_hic_links_filt.tsv | perl -lane 'for($x = 0; $x < scalar(@F); $x++){print "$x\t$F[$x]";}' | grep '20805'  <- the contig with a big smear in the middle of the plot

# Identifying the clusters associated with this contig
perl -lane 'if($F[82] > 4){print "$F[0]\t$F[82]";}' < pb_arg_hic_links_filt.tsv

```

| Hi-C Bin | MarkerGene tax | Mash profile |
| :--- | :--- | :--- |
| pacbio_final_public_hic.194.contigs.fa*(original bin)| p__Bacteroidetes |	Bacteroidales_bacterium_ph8.fna |
| pacbio_final_public_hic.217.contigs.fa | k__Bacteria	| Francisella_tularensis_subsp._tularensis_TI0902.fna |
| pacbio_final_public_hic.588.contigs.fa | k__Bacteria  |	Prevotella_ruminicola_23.fna |
| pacbio_final_public_hic.617.contigs.fa | k__Bacteria	| Ruminobacter_sp._RM87.fna |
| pacbio_final_public_hic.907.contigs.fa | k__Bacteria	| Methanothermobacter_thermautotrophicus_str._Delta_H.fna |
| pacbio_final_public_hic.1061.contigs.fa | root	| Clostridiales_bacterium_oral_taxon_876_str._F0540.fna

```bash
head -n 1 pb_arg_hic_links_filt.tsv | perl -lane 'for($x = 0; $x < scalar(@F); $x++){print "$x\t$F[$x]";}' | grep '28353' <- the first contig in the list

perl -lane 'if($F[111] > 4){print "$F[0]\t$F[111]";}' < pb_arg_hic_links_filt.tsv
# I found allot of associations with fibrolytic strains here
```


## Checking DAS_tool dereplication for binning comparison

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); %hic; while(<IN>){chomp; @s = split(/\t/); $hic{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[2]"); %meta; while(<IN>){chomp; @s = split(/\t/); $meta{$s[0]} = $s[1];} close IN; my %total; foreach my $k (keys(%data), keys(%hic), keys(%meta)){$total{$k} = 1;} foreach my $k (sort {$a cmp $b} keys(%total)){$dbin = (exists($data{$k}))? $data{$k} : "NONE"; $hbin = (exists($hic{$k}))? $hic{$k} : "NONE"; $mbin = (exists($meta{$k}))? $meta{$k} : "NONE"; print "$k\t$dbin\t$hbin\t$mbin\n";}' illumina_final_dastool_DASTool_scaffolds2bin.txt illumina_megahit_hic.final.bins illumina_megahit_public_metabat.unsorted.bins > illumina_final_dastool_DASTool_scaffolds2bin.table.tab

perl -e '%data; while(<>){chomp; @s = split(/\t/); if($s[1] eq "NONE" || $s[3] eq "NONE"){next;} $data{$s[1]}->{$s[3]} += 1;} foreach my $d (keys(%data)){print "$d"; foreach my $m (sort{$data{$d}->{$b} <=> $data{$d}->{$a}} keys(%{$data{$d}})){ print "\t$m\t" . $data{$d}->{$m} . "\n"; last;}}' < illumina_final_dastool_DASTool_scaffolds2bin.table.tab > illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count
grep -v 'hic' illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count > illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count.filt

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f illumina_megahit_public_metabat.unsorted.bins -c 1 | grep -v 'Entry' > illumina_final_dastool_DASTool_scaffolds2bin.metabat.prederep.count

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); push(@s, $data{$s[1]}, $data{$s[1]} - $s[2]); print join("\t", @s) . "\n";} close IN;' illumina_final_dastool_DASTool_scaffolds2bin.metabat.prederep.count illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count.filt > illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count.comp
cat illumina_final_dastool_DASTool_scaffolds2bin.metabat.best.count.comp | cut -f5 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
Sum:    64485
total   928
Minimum 0
Maximum 3781
Average 69.488147
Median  28
Standard Deviation      204.996319
Mode(Highest Distributed Value) 0


perl -e '%data; while(<>){chomp; @s = split(/\t/); if($s[1] eq "NONE" || $s[2] eq "NONE"){next;} $data{$s[1]}->{$s[2]} += 1;} foreach my $d (keys(%data)){print "$d"; foreach my $m (sort{$data{$d}->{$b} <=> $data{$d}->{$a}} keys(%{$data{$d}})){ print "\t$m\t" . $data{$d}->{$m} . "\n"; last;}}' < illumina_final_dastool_DASTool_scaffolds2bin.table.tab > illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count
grep -v metabat illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count > illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count.filt

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f illumina_megahit_hic.final.bins -c 1 | grep -v 'Entry' > illumina_final_dastool_DASTool_scaffolds2bin.hic.prederep.count
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); push(@s, $data{$s[1]}, $data{$s[1]} - $s[2]); print join("\t", @s) . "\n";} close IN;' illumina_final_dastool_DASTool_scaffolds2bin.hic.prederep.count illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count.filt > illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count.comp

cat illumina_final_dastool_DASTool_scaffolds2bin.hic.best.count.comp | cut -f5 | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
Sum:    72306
total   3056
Minimum 0
Maximum 356
Average 23.660340
Median  14
Standard Deviation      30.951519
Mode(Highest Distributed Value) 0
```


## Reviewer responses

Here are my notes for generating statistics and data to respond to the reviewer comments on our manuscript.

> Ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/dastool

```bash
# Getting the Bowers standard counts:
# High-quality draft (HQ)
perl -lane 'if($F[-2] > 90 && $F[-1] < 5){print $F_;}' < pacbio_final_dastool_DASTool_summary.txt | wc -l
10
# Medium-quality draft (MQ)
perl -lane 'if($F[-2] > 50 && $F[-1] < 10){print $F_;}' < pacbio_final_dastool_DASTool_summary.txt | wc -l
103

# High-quality (HQ)
perl -lane 'if($F[-2] > 90 && $F[-1] < 5){print $F_;}' < illumina_final_dastool_DASTool_summary.txt | wc -l
42
# Medium-quality (MQ)
perl -lane 'if($F[-2] > 50 && $F[-1] < 10){print $F_;}' < illumina_final_dastool_DASTool_summary.txt | wc -l
325

# Now, printing them out so that I can generate new master tables
# Illumina
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 5 || $s[-2] < 90){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa");while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_DASTool_HQD.tab
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 10 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< illumina_final_dastool_DASTool_bins/$s[0].contigs.fa");while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < illumina_final_dastool_DASTool_summary.txt > illumina_final_dastool_DASTool_MQD.tab

# Pacbio
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 5 || $s[-2] < 90){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_final_dastool_DASTool_HQD.tab
perl -e 'while(<>){chomp; @s = split(/\t/); if($s[0] eq "bin" || $s[-1] >= 10 || $s[-2] < 50){next;}else{@bsegs = split(/\./, $s[0]); open(IN, "< pacbio_final_dastool_DASTool_bins/$s[0].contigs.fa"); while($l = <IN>){if($l =~ /^>/){chomp $l; $l =~ s/>//g; @jsegs = split(/[_\.]/, $bsegs[-2]); print "$l\t$jsegs[-1]\_$bsegs[-1]\n";}} close IN;}}' < pacbio_final_dastool_DASTool_summary.txt > pacbio_final_dastool_DASTool_MQD.tab
```

#### Generating new master tables and calculating stats

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Editting the json files
cp illumina_data_files_11_2018.json illumina_data_files_4_2019.json
cp pacbio_data_files_9_2018.json pacbio_data_files_4_2019.json

python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j illumina_data_files_4_2019.json -t illumina_preliminary_blobtools_table.tab -o illumina_megahit_master_table_4_1_2019.tab
python3 ~/python_toolchain/metagenomics/addJSONColumnsToTable.py -j pacbio_data_files_4_2019.json -t pacbio_preliminary_blobtools_table.tab -o pacbio_final_pilon_master_table_4_1_2019.tab

# Subsectioning them for ease of use in my other analysis pipelines
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < illumina_megahit_master_table_4_1_2019.tab > illumina_megahit_master_table_4_1_2019.MQbins.short.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] eq "NOBIN"){next;} print join("\t", @F[0,1,2,22,23,27,31,35,39,47,48,49,50,51,52,53,57,58]) . "\n";' < illumina_megahit_master_table_4_1_2019.tab > illumina_megahit_master_table_4_1_2019.HQbins.short.tab
```

#### Calculating new stats for the manuscript

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/prodigal

```bash
# 16S stats
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab -c 0 -l mick_16s_pacbio_summary.contig.list | perl -ne 'chomp; @F = split(/\t/); print "$F[11]\t$F[0]\t$F[4]\t$F[10]\t$F[12]\n";' > mick_16s_pacbio_summary.reformat.MQbins.tab
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f ../master_tables/pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab -c 0 -l mick_16s_pacbio_summary.contig.list | perl -ne 'chomp; @F = split(/\t/); print "$F[11]\t$F[0]\t$F[4]\t$F[10]\t$F[12]\n";' > mick_16s_pacbio_summary.reformat.HQbins.tab

wc -l mick_16s_pacbio_summary.reformat.MQbins.tab mick_16s_pacbio_summary.reformat.HQbins.tab               64 mick_16s_pacbio_summary.reformat.MQbins.tab
  14 mick_16s_pacbio_summary.reformat.HQbins.tab
  78 total

python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f mick_16s_pacbio_summary.reformat.MQbins.tab -c 2 -d '\t'
Entry   Value
Bacteria        59
Eukaryota       4
Archaea 1
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Dataset comparison counts
perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   22
Sum:    317319
Minimum 3552
Maximum 30215
Average 14423.590909
Median  13725
Standard Deviation      7221.517982

perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < illumina_megahit_master_table_4_1_2019.HQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   3650
Sum:    25491578
Minimum 1000
Maximum 90246
Average 6983.993973
Median  4899.5
Standard Deviation      6453.684480

perl -lane 'if($F[-2] eq "-" && $F[-1] eq "-"){print $F[1];}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   1254
Sum:    16151563
Minimum 1508
Maximum 367662
Average 12880.034290
Median  9669
Standard Deviation      14777.483832

# The "does not align to any other dataset" stats
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[1]\n";}' < illumina_megahit_master_table_4_1_2019.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   27120
Sum:    87817564
Minimum 1000
Maximum 71737
Average 3238.110767
Median  2455
Standard Deviation      2838.136951

perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[1]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   20
Sum:    137081
Minimum 1762
Maximum 16216
Average 6854.050000
Median  7242.5
Standard Deviation      4053.889237

# Now to find out if we still have firmicutes in the PACBIO dataset
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[0]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_MQ_unique_ctgs.list
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN" && $F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-"){print "$F[0]\n";}' < pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_final_pilon_MQ_unique_ctgs.list

# Still Firmicutes and the statistic still stands
python3 ~/python_toolchain/utils/tabFileColumnGrep.py -f pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab -c 0 -d '\t' -l pacbio_final_pilon_MQ_unique_ctgs.list | perl -lane 'print $F[2];' > pacbio_final_pilon_MQ_unique_ctgs.list.gc
perl -lane 'if($F[0] eq "name"){next;}else{print $F[2];}' < pacbio_final_pilon_master_table_4_1_2019.MQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.MQbins.short.gc

##### Generating new supplementary table for HQ bins for Reviewer #2
# First, let's calculate the total coverage for all datasets from the short-reads
# The coverage values are averages, so I need to multiply by the lengths to get an estimate
perl -e '$sum = 0; $sa = 0; $n = 0; while(<>){chomp; @F = split(/\t/); if($F[0] eq "name"){next;} $sum += ($F[9] * $F[1]); $sa += $F[9]; $n++}; $avg = $sa / $n; print "$sum\t$sa\t$n\t$avg\n";' < pacbio_final_pilon_master_table_4_1_2019.tab
66888601387.3526        4652765.74600004        77670   59.9042840993954
perl -e '$sum = 0; $sa = 0; $n = 0; while(<>){chomp; @F = split(/\t/); if($F[0] eq "name"){next;} $sum += ($F[9] * $F[1]); $sa += $F[9]; $n++}; $avg = $sa / $n; print "$sum\t$sa\t$n\t$avg\n";' < illumina_megahit_master_table_4_1_2019.tab
129991919911.382        54111883.2900027        2182263 24.7962245109791

# OK, it's a slight overestimate, but still in the ballpark
# Let's generate a "proportion of aligned reads" statistic per HQ bin

perl consolidate_HQ_bin_table.pl illumina_megahit_master_table_4_1_2019.tab 129991919911 illumina_megahit_HQBIN_supplementary_summary.tab
perl consolidate_HQ_bin_table.pl pacbio_final_pilon_master_table_4_1_2019.tab 66888601387 pacbio_final_HQBIN_supplementary_summary.tab

#OK so most of the data worked out great, but the proportional coverage is bad. Damn. I'm going to have to calculate this by hand from the bams!
```


#### Creating read subsamples for Phase team

We're testing the hypothesis that PE links are more informative than Hi-C data in binning metagenome assembly contigs. 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/pilot_project/illumina

```bash
# I think that the files need to be unzipped for this, unfortunately!
for i in YMPrepCannula_S1_*.fastq.gz; do sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p msn --wrap="gunzip $i"; done

sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq.gz -s YMPrepCannula_S1_L001_R2_001.fastq.gz -f YMPrepCannula_S1_L002_R1_001.fastq.gz -s YMPrepCannula_S1_L002_R2_001.fastq.gz -f YMPrepCannula_S1_L003_R1_001.fastq.gz -s YMPrepCannula_S1_L003_R2_001.fastq.gz -f YMPrepCannula_S1_L004_R1_001.fastq.gz -s YMPrepCannula_S1_L004_R2_001.fastq.gz -o illumina_pe_reads_m1_subsample -l 22487509"

# Ooops! Looks like gziping doesn't work for file line counts
sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq -s YMPrepCannula_S1_L001_R2_001.fastq -f YMPrepCannula_S1_L002_R1_001.fastq -s YMPrepCannula_S1_L002_R2_001.fastq -f YMPrepCannula_S1_L003_R1_001.fastq -s YMPrepCannula_S1_L003_R2_001.fastq -f YMPrepCannula_S1_L004_R1_001.fastq -s YMPrepCannula_S1_L004_R2_001.fastq -o illumina_pe_reads_m1_subsample -l 22487509 -i subsample_m1.log -b 2280780244"
sbatch --nodes=1 --mem=25000 --ntasks-per-node=5 -p msn --wrap="python3 ~/python_toolchain/sequenceData/randomSubsampleFastqs.py -f YMPrepCannula_S1_L001_R1_001.fastq -s YMPrepCannula_S1_L001_R2_001.fastq -f YMPrepCannula_S1_L002_R1_001.fastq -s YMPrepCannula_S1_L002_R2_001.fastq -f YMPrepCannula_S1_L003_R1_001.fastq -s YMPrepCannula_S1_L003_R2_001.fastq -f YMPrepCannula_S1_L004_R1_001.fastq -s YMPrepCannula_S1_L004_R2_001.fastq -o illumina_pe_reads_s3_subsample -l 40889499 -i subsample_s3.log -b 2280780244"
```

I sent the files off to Max for analysis. We'll see what he finds!

#### Generating log differential coverage plots

The first reviewer wants to see log differential coverage plots. I had generated these with Blobplots, but I think he wants us to use his pipeline. First, let's align the Hi-C reads and error corrected reads to the HQ bins I just prepared.

Given the fact that MMGenome2 can plot a grid of differential coverage, let's see if we can plot ALL datasets (Hi-C, PacbioEC, normal short-read, 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/diff_cov

```bash
module load bwa samtools metabat/2.12.1
# Generating the lists I need
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print join("\t", @F[(0,1,2,9,14,15,16,17,21,8,23,27,39,51,52,53)]) . "\n";}' < ../master_tables/illumina_megahit_master_table_4_1_2019.tab > illumina_HQBIN_only_shortcov.tab
perl -ne 'chomp; @F = split(/\t/); if($F[51] ne "NOBIN"){print join("\t", @F[(0,1,2,9,14,15,16,17,21,8,23,27,39,51,52,53)]) . "\n";}' < ../master_tables/pacbio_final_pilon_master_table_4_1_2019.tab > pacbio_HQBIN_only_shortcov.tab

# Generating the HQ bin fastas
perl -lane 'if($F[0] eq "name"){next;}else{print $F[0];}' < illumina_HQBIN_only_shortcov.tab > illumina_HQBIN_only_shortcov.ctg.list
samtools faidx -r illumina_HQBIN_only_shortcov.ctg.list ../../assemblies/pilot_project/illumina_megahit/illumina_megahit_final_contigs.perl.fa > illumina_megahit_final_onlyHQBIN.fa
samtools faidx illumina_megahit_final_onlyHQBIN.fa
sbatch --nodes=1 --mem=10000 --ntasks-per-node=1 -p msn --wrap="bwa index illumina_megahit_final_onlyHQBIN.fa"

perl -lane 'if($F[0] eq "name"){next;}else{print $F[0];}' < pacbio_HQBIN_only_shortcov.tab > pacbio_HQBIN_only_shortcov.ctg.list
samtools faidx -r pacbio_HQBIN_only_shortcov.ctg.list ../../assemblies/pilot_project/pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > pacbio_final_pilon_onlyHQBIN.fa
samtools faidx pacbio_final_pilon_onlyHQBIN.fa
sbatch --nodes=1 --ntasks-per-node=1 --mem=10000 -p msn --wrap="bwa index pacbio_final_pilon_onlyHQBIN.fa"

# OK, I've got our illumina coverage down, but I haven't aligned the pacbio EC reads using bwa mem yet. Let's queue that up
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem illumina_megahit_final_onlyHQBIN.fa ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta | samtools sort -o illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam -T illumina_pbtemp - ; samtools index illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam"
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem pacbio_final_pilon_onlyHQBIN.fa ../../sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta | samtools sort -o pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam -T pacbio_pbtemp - ; samtools index pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam"

# Let's also queue up the hi-c read pairs
ls /beegfs/project/rumen_longread_metagenome_assembly/sequence_data/pilot_project/hic/*.fastq > hic_reads_spreadsheet.tab
# Editted the file to fit my pipeline's preferred format
vim hic_reads_spreadsheet.tab

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b illumina_hic -t hic_reads_spreadsheet.tab -f /beegfs/project/rumen_longread_metagenome_assembly/analysis/diff_cov/illumina_megahit_final_onlyHQBIN.fa -m -p msn
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b pacbio_hic -t hic_reads_spreadsheet.tab -f /beegfs/project/rumen_longread_metagenome_assembly/analysis/diff_cov/pacbio_final_pilon_onlyHQBIN.fa -m -p msn

sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p msn --wrap="jgi_summarize_bam_contig_depths --outputDepth illumina_hic_jgi_depth.tab illumina_hic/M1HIC/M1HIC.sorted.merged.bam illumina_hic/S3HIC/S3HIC.sorted.merged.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p msn --wrap="jgi_summarize_bam_contig_depths --outputDepth pacbio_hic_jgi_depth.tab pacbio_hic/M1HIC/M1HIC.sorted.merged.bam pacbio_hic/S3HIC/S3HIC.sorted.merged.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="jgi_summarize_bam_contig_depths --outputDepth pacbio_final_pilon_onlyHQBIN.pbecreads.depth.tab pacbio_final_pilon_onlyHQBIN.pbecreads.sorted.bam"
sbatch --nodes=1 --mem=20000 --ntasks-per-node=2 -p short --wrap="jgi_summarize_bam_contig_depths --outputDepth illumina_megahit_final_onlyHQBIN.pbecreads.depth.tab illumina_megahit_final_onlyHQBIN.pbecreads.sorted.bam"


# Now to start consolidating everything together
perl -lane 'print "$F[0]\t$F[3]";' < pacbio_HQBIN_only_shortcov.tab > pacbio_HQBIN_only_shortcov.temp.tab
perl -lane 'print "$F[0]\t$F[3]\t$F[5]";' < pacbio_hic_jgi_depth.tab > pacbio_hic_jgi_depth.temp
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' pacbio_hic_jgi_depth.temp pacbio_HQBIN_only_shortcov.temp.tab > pacbio_HQBIN_shortread_avgcov.tab

perl -lane 'print "$F[0]\t$F[3]";' < illumina_HQBIN_only_shortcov.tab > illumina_HQBIN_only_shortcov.temp.tab
perl -lane 'print "$F[0]\t$F[3]\t$F[5]";' < illumina_hic_jgi_depth.tab > illumina_hic_jgi_depth.temp.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' illumina_hic_jgi_depth.temp.tab illumina_HQBIN_only_shortcov.temp.tab > illumina_HQBIN_shortread_avgcov.tab

# Now for the PB reads vs Illumina reads on both assemblies
perl -lane 'print "$F[0]\t$F[3]";' < illumina_megahit_final_onlyHQBIN.pbecreads.depth.tab > illumina_megahit_final_onlyHQBIN.pbecreads.depth.temp
perl -lane 'print "$F[0]\t$F[3]";' < pacbio_final_pilon_onlyHQBIN.pbecreads.depth.tab > pacbio_final_pilon_onlyHQBIN.pbecreads.depth.temp


perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' illumina_megahit_final_onlyHQBIN.pbecreads.depth.temp illumina_HQBIN_only_shortcov.temp.tab > illumina_HQBIN_sr_vs_lr_avgcov.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %data; while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); $data{$name} = [@s];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $name = shift(@s); print "$name\t" . join("\t", @s) . "\t" . join("\t", @{$data{$name}}) . "\n";} close IN;' pacbio_final_pilon_onlyHQBIN.pbecreads.depth.temp pacbio_HQBIN_only_shortcov.temp.tab > pacbio_HQBIN_sr_vs_lr_avgcov.tab
```

##### New ideas: 
1. Mapping efficiency plays a big role in Pacbio read alignment, so subsection the reads?
2. Also check to see if there's a phylum bias in the bins. This will be easier. Not sure if the pacbio alignments will change much, but we can do them to address any questions about alignment bias


```bash
module load bwa samtools metabat/2.12.1
# Generating new PacBio short read fastqs
sbatch --nodes=1 --ntasks-per-node=6 --mem=20000 -p msn --wrap="python3 ~/python_toolchain/sequenceData/sampleLongReadToShortRead.py -f /beegfs/project/rumen_longread_metagenome_assembly/sequence_data/pilot_project/pacbio/rumen_pacbio_corrected.fasta -s rumen_pacbio_corrected.sr.R1.fq -r rumen_pacbio_corrected.sr.R2.fq -t 6"

sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem illumina_megahit_final_onlyHQBIN.fa rumen_pacbio_corrected.sr.R1.fq rumen_pacbio_corrected.sr.R2.fq | samtools sort -o illumina_HQBIN.pbecreads.leapfrog.sorted.bam -T illumina_pbtemp -; samtools index illumina_HQBIN.pbecreads.leapfrog.sorted.bam; jgi_summarize_bam_contig_depths --outputDepth illumina_pb_ecreads_leapfrog_depth.tab illumina_HQBIN.pbecreads.leapfrog.sorted.bam"
sbatch --nodes=1 --ntasks-per-node=1 --mem=12000 -p msn --wrap="bwa mem pacbio_final_pilon_onlyHQBIN.fa rumen_pacbio_corrected.sr.R1.fq rumen_pacbio_corrected.sr.R2.fq | samtools sort -o pacbio_HQBIN.pbecreads.leapfrog.sorted.bam -T pacbio_pbtemp -; samtools index pacbio_HQBIN.pbecreads.leapfrog.sorted.bam; jgi_summarize_bam_contig_depths --outputDepth pacbio_pb_ecreads_leapfrog_depth.tab pacbio_HQBIN.pbecreads.leapfrog.sorted.bam"

wc -l rumen_pacbio_corrected.sr.R1.fq
399818640 rumen_pacbio_corrected.sr.R1.fq
# That's approximately 99 million read pairs

# Now to prep them the same way as before and generate the plots for comparison
perl -lane 'if($F[0] eq "contigName"){next;} print "$F[0]\t$F[3]";' < pacbio_pb_ecreads_leapfrog_depth.tab > pacbio_pb_ecreads_leapfrog_depth.temp
perl -lane 'if($F[0] eq "contigName"){next;} print "$F[0]\t$F[3]";' < illumina_pb_ecreads_leapfrog_depth.tab > illumina_pb_ecreads_leapfrog_depth.temp

 perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq "contigName"){print join("\t", @s) . "\n";next;} $long = $h{$s[0]}; print "$s[0]\t$s[1]\t$long\n";} close IN;' illumina_pb_ecreads_leapfrog_depth.temp illumina_HQBIN_sr_vs_lr_avgcov.tab > illumina_HQBIN_sr_vs_lr_avgcov.leapfrog.tab
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq "contigName"){print join("\t", @s) . "\n";next;} if(!exists($h{$s[0]})){print STDERR "Error!\n";} $long = $h{$s[0]}; print "$s[0]\t$s[1]\t$long\n";} close IN;' pacbio_pb_ecreads_leapfrog_depth.temp pacbio_HQBIN_sr_vs_lr_avgcov.tab > pacbio_HQBIN_sr_vs_lr_avgcov.leapfrog.tab
```

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/analysis/master_tables

```bash
# Generating the tax consensus for the contigs in the HQ bins
perl generate_tax_table.pl pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab > pacbio_final_pilon_master_table_4_1_2019.HQbins.taxtable.tab
perl generate_tax_table.pl illumina_megahit_master_table_4_1_2019.HQbins.short.tab > illumina_megahit_master_table_4_1_2019.HQbins.taxtable.tab

# I am just pulling the BlobTools phyla information (most stable taxonomic category in the bins it appears)
perl -ne 'chomp; @F = split(/\t/); print "$F[0]\t$F[5]\n";' < pacbio_final_pilon_master_table_4_1_2019.HQbins.short.tab > ../diff_cov/pacbio_HQbins_taxassignment.tab
perl -ne 'chomp; @F = split(/\t/); print "$F[0]\t$F[5]\n";' < illumina_megahit_master_table_4_1_2019.HQbins.short.tab > ../diff_cov/illumina_HQbins_taxassignment.tab
```


Now to plot the coverage differential plots in mmgenome2.

> F:/SharedFolders/metagenomics/pilot_manuscript/diff_cov/

```R
# PacbioHIC
library(mmgenome2)
covdataframe <- read.delim("pacbio_HQBIN_shortread_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe)

pdf(file="pacbio_short_read_M1HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_M1HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

pdf(file="pacbio_short_read_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_S3HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

pdf(file="pacbio_M1HIC_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_M1HIC", y = "cov_S3HIC", color_by = "gc", x_scale= "log10", y_scale = "log10")
dev.off()

#IlluminaHIC
covdataframe <- read.delim("illumina_HQBIN_shortread_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe)

pdf(file="illumina_short_read_M1HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_M1HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_short_read_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_USDA", y = "cov_S3HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_M1HIC_S3HIC_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_M1HIC", y = "cov_S3HIC", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()


# Illumina sr vs lr
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe)

pdf(file="illumina_sr_vs_lr_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

# Generating coverage by len plots to see if there are biases
mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="illumina_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="illumina_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()

# I just want to make sure that the short read and long read data doesn't have an obvious bias. Let's get the coefficients from a linear model to make sure
ilsrlm <- lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)
summary(ilsrlm)

Call:
lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)

Residuals:
   Min     1Q Median     3Q    Max 
-33.22 -14.55  -7.16   4.38 807.70 

Coefficients:
                  Estimate  Std. Error t value Pr(>|t|)    
(Intercept)    45.26095831  1.85465448  24.404  < 2e-16 ***
mm.trim$length -0.00019230  0.00002624  -7.328 2.61e-13 ***
mm.trim$gc     -0.28991344  0.03881071  -7.470 9.02e-14 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 28.81 on 6855 degrees of freedom
Multiple R-squared:  0.01357,	Adjusted R-squared:  0.01328 
F-statistic: 47.14 on 2 and 6855 DF,  p-value: < 2.2e-16

illrlm <- lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)
summary(illrlm)

Call:
lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
 -4.246  -2.629  -1.719   1.120 202.217 

Coefficients:
                   Estimate   Std. Error t value Pr(>|t|)    
(Intercept)    -0.356198367  0.339818384  -1.048   0.2946    
mm.trim$length -0.000010549  0.000004808  -2.194   0.0283 *  
mm.trim$gc      0.067749436  0.007111078   9.527   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 5.278 on 6855 degrees of freedom
Multiple R-squared:  0.0151,	Adjusted R-squared:  0.01481 
F-statistic: 52.54 on 2 and 6855 DF,  p-value: < 2.2e-1


# Pacbio sr vs lr
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe)

pdf(file = "pacbio_sr_vs_lr_covplot.pdf", useDingbats = FALSE)
mmplot(mm, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

mm.trim <- mm %>% filter(cov_short != "NA")
pdf(file="pacbio_sr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_short)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()
pdf(file="pacbio_lr_vs_length_binplot.pdf", useDingbats = FALSE)
ggplot(mm.trim, aes(x=length, y=cov_long)) + geom_bin2d() + theme_bw() + scale_x_log10() + scale_y_log10()
dev.off()

pbsrlm <- lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)
summary(pbsrlm)

Call:
lm(formula = mm.trim$cov_short ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
-88.386 -18.075  -5.624  11.008 174.108 

Coefficients:
                   Estimate   Std. Error t value    Pr(>|t|)    
(Intercept)    109.37309258  11.39537216   9.598     < 2e-16 ***
mm.trim$length   0.00024681   0.00002309  10.691     < 2e-16 ***
mm.trim$gc      -1.22724830   0.24173592  -5.077 0.000000573 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 29.84 on 430 degrees of freedom
Multiple R-squared:  0.2597,	Adjusted R-squared:  0.2562 
F-statistic: 75.42 on 2 and 430 DF,  p-value: < 2.2e-16

pblrlm <- lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)
summary(pblrlm)

Call:
lm(formula = mm.trim$cov_long ~ mm.trim$length + mm.trim$gc)

Residuals:
    Min      1Q  Median      3Q     Max 
-13.517  -2.942  -0.721   2.172  41.477 

Coefficients:
                   Estimate   Std. Error t value Pr(>|t|)    
(Intercept)    -0.093548111  1.781133609  -0.053  0.95814    
mm.trim$length  0.000033695  0.000003608   9.338  < 2e-16 ***
mm.trim$gc      0.110666364  0.037784107   2.929  0.00358 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.664 on 430 degrees of freedom
Multiple R-squared:  0.1757,	Adjusted R-squared:  0.1719 
F-statistic: 45.83 on 2 and 430 DF,  p-value: < 2.2e-16


# Let's try a plot with the contigs colored by phyla
# Pacbio by phyla
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
phyladatafram <- read.delim("pacbio_HQbins_taxassignment.tab", header=TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe, taxonomy = phyladatafram)

# One of the contigs was misassigned to the Arthropoda, let's change it to no-hit
mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda", "no-hit"))

pdf(file = "pacbio_sr_vs_lr_covplot.byphylum.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Illumina by phyla
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header =TRUE)
phyladatafram <- read.delim("illumina_HQbins_taxassignment.tab", header = TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe, taxonomy = phyladatafram)

mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda" | phylum == "Chordata" | phylum == "Mollusca" | phylum == "Nematoda", "no-hit"))
mm.trim <- mm.trim %>% mutate(phylum=replace(phylum, phylum == "Acidobacteria" | phylum == "Candidatus Desantisbacteria" | phylum == "Candidatus Melainabacteria" | phylum == "Candidatus Moranbacteria" | phylum == "Candidatus Nomurabacteria" | phylum == "Candidatus Saccharibacteria" | phylum == "Chlorobi" | phylum == "Chlorophyta" | phylum == "Elusimicrobia" | phylum == "Gemmatimonadetes" | phylum == "Kiritimatiellaeota" | phylum == "Thermotogae", "other"))

pdf(file="illumina_sr_vs_lr_covplot.byphylum.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Now trying the Leapfrog read alignments
covdataframe <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.leapfrog.tab", header=TRUE)
phyladatafram <- read.delim("illumina_HQbins_taxassignment.tab", header = TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/illumina_megahit_final_contigs.perl.fa", coverage=covdataframe, taxonomy = phyladatafram)

mm.trim <- mm %>% mutate(phylum=replace(phylum, phylum == "Arthropoda" | phylum == "Chordata" | phylum == "Mollusca" | phylum == "Nematoda", "no-hit"))
mm.trim <- mm.trim %>% mutate(phylum=replace(phylum, phylum == "Acidobacteria" | phylum == "Candidatus Desantisbacteria" | phylum == "Candidatus Melainabacteria" | phylum == "Candidatus Moranbacteria" | phylum == "Candidatus Nomurabacteria" | phylum == "Candidatus Saccharibacteria" | phylum == "Chlorobi" | phylum == "Chlorophyta" | phylum == "Elusimicrobia" | phylum == "Gemmatimonadetes" | phylum == "Kiritimatiellaeota" | phylum == "Thermotogae", "other"))

pdf(file="illumina_sr_vs_lr_covplot.bygc.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file="illumina_sr_vs_lr_covplot.byphylum.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y = "cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

# Testing zero coverage contigs based on mapping bias
prevcovdata <- read.delim("illumina_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
count(prevcovdata[prevcovdata$long == 0,])
1  2836
count(covdataframe[covdataframe$long == 0,])
1   178  <- less mapping bias, but it's still there!

# pacbio
covdataframe <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.leapfrog.tab", header = TRUE)
phyladatafram <- read.delim("pacbio_HQbins_taxassignment.tab", header=TRUE)
colnames(phyladatafram) <- c("scaffold", "phylum")
mm <- mmload("F:/Globus/usda_pacbio_second_pilon_indelsonly.fa", coverage=covdataframe, taxonomy = phyladatafram)

pdf(file = "pacbio_sr_vs_lr_covplot.bygc.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "gc", x_scale = "log10", y_scale = "log10")
dev.off()

pdf(file = "pacbio_sr_vs_lr_covplot.byphylum.leap.pdf", useDingbats = FALSE)
mmplot(mm.trim, x = "cov_short", y="cov_long", color_by = "phylum", x_scale = "log10", y_scale = "log10")
dev.off()

prevcovdata <- read.delim("pacbio_HQBIN_sr_vs_lr_avgcov.tab", header=TRUE)
cor.test(covdataframe$short, covdataframe$long)

	Pearson's product-moment correlation

data:  covdataframe$short and covdataframe$long
t = 20.767, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6567436 0.7513733
sample estimates:
      cor 
0.7072118

count(prevcovdata[prevcovdata$long == 0,])
1     8
count(covdataframe[covdataframe$long == 0,])
1     0
```