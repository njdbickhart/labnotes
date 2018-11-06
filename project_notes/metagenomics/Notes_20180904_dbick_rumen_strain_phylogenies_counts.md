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
perl -ne 'chomp; @F = split(/\t/); if($F[50] eq "NOBIN"){next;} print join("\t", @F[0,1,2,48,49,50,51]) . "\n";' < pacbio_final_pilon_master_table_2018_09_07.tab > supplementary_table_3_longread_bins.tab
```

#### Supplementary tables 4 and 5

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-" && $F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51]) . "\n";}' < illumina_megahit_master_table_2018_09_07.tab > supplementary_table_4_short_read_unique.tab
perl -ne 'chomp; @F = split(/\t/); if($F[56] eq "-" && $F[57] eq "-" && $F[58] eq "-" && $F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51]) . "\n";}' < pacbio_final_pilon_master_table_2018_09_07.tab > supplementary_table_5_long_read_unique.tab
```

#### Supplementary tables 6 and 7

```bash
perl -ne 'chomp; @F = split(/\t/); if($F[50] ne "NOBIN"){print join("\t", @F[0,1,2,8,9,14,15,16,17,21,52,53,23,27,31,35,39,43,48,49,50,51,56,57,58]) . "\n";}' < illumina_megahit_master_table_2018_09_07.tab > supplementary_table_6_short_read_taxonomy.tab
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