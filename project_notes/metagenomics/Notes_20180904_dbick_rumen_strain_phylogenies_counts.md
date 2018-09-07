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