# Comparing marker placement on the new Pig assemblies
---
*6/22/2018*

These are my notes of the alignment and ordering of SNP markers on the new Pig assemblies for comparative analysis.

## Table of Contents


## Preparing the files

> Assembler2: /mnt/nfs/nfs2/bickhart-users/pig_genomes

```bash
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Sus_scrofa/ARCHIVE/ANNOTATION_RELEASE.104/Assembled_chromosomes/seq/*.fa.gz'

for i in ssc_ref*.fa.gz; do echo $i; unpigz $i; done
# Reformating fasta names and printing to external file
for i in `seq 1 18` MT X Y; do echo $i; perl -ne 'if($_ =~ /^>/){@bsegs = split(/\|/); print ">$bsegs[3]\n";}else{print $_;}' < ssc_ref_Sscrofa10.2_chr${i}.fa >> ssc_10.2_ncbi.reference.fa; done

perl -ne 'if($_ =~ /^>/){chomp; @bsegs = split(/\s+/); print ">$bsegs[-1]\n";}else{print $_;}' < ssc_ref_Sscrofa10.2_unplaced.fa >> ssc_10.2_ncbi.reference.fa

# Indexing the assemblies that will be used for the alignment
sbatch --nodes=1 --ntasks-per-node=1 --mem=20000 -p assemble2 --wrap="bwa index ssc_10.2_ncbi.reference.fa;"
sbatch --nodes=1 --ntasks-per-node=1 --mem=20000 -p assemble2 --wrap="bwa index GCA_002844635.1_USMARCv1.0_genomic.fna;"

# Converting SNP file tsvs into marker fasta files
dos2unix ggp_marker_sites.tab
dos2unix pig_60k_marker_sites.tab
dos2unix pig_80k_marker_sites.tab

perl -ne '@s = split(/\t/); print ">$s[0].$s[9].$s[10]\n$s[5]\n";' < pig_60k_marker_sites.tab > pig_60k_marker_sites.fa
perl -ne '@s = split(/\t/); print ">$s[0].$s[9].$s[10]\n$s[5]\n";' < ggp_marker_sites.tab > ggp_marker_sites.fa
perl -ne '@s = split(/\t/); print ">$s[0].$s[9].$s[10]\n$s[5]\n";' < pig_80k_marker_sites.tab > pig_80k_marker_sites.fa

# Preparing the markers from the Axiom array
perl -lane 'if($F[0] =~ /^#/ || $F[0] =~ /^\"Probe/){next;} $F[0] =~ s/\"//g; @segs = split(/[,\[\/]/, $F[0]); print ">$segs[0].$segs[3].$segs[4]\n$segs[6]";' < Axiom_PigHD_v1_Annotation.r3.csv > Axiom_PigHD_v1_Annotation.r3.forprobe.fa
```

## Alignment

```bash
# Aligning to USMarc
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCA_002844635.1_USMARCv1.0_genomic.fna -p $i.fa -o $i.USMARC.snps; done

# Aligning to ss11.1
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCF_000003025.6_Sscrofa11.1_genomic.fna -p $i.fa -o $i.ROSLIN.snps; done

# Finally, aligning to ss10.2
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a ssc_10.2_ncbi.reference.fa -p $i.fa -o $i.SS10.snps; done

module load samtools bwa; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCA_002844635.1_USMARCv1.0_genomic.fna -p Axiom_PigHD_v1_Annotation.r3.forprobe.fa -o axiom_pigHD.USMARC.snps
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCF_000003025.6_Sscrofa11.1_genomic.fna -p Axiom_PigHD_v1_Annotation.r3.forprobe.fa -o axiom_pigHD.ROSLIN.snps
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a ssc_10.2_ncbi.reference.fa -p Axiom_PigHD_v1_Annotation.r3.forprobe.fa -o axiom_pigHD.SS10.snps
```

Gah! The 80k snps have HORRIBLE names! Gotta fix the file

> Ceres: /home/derek.bickharhth/cattle_genome_assemblies/Pig/pig_genomes

```bash
perl -ne 'if($_ =~ /^>/){$_ =~ s/LD-Porcine80K_WU_10\.//; if($_ =~ /LD-Porcine80K_M1GA0025150-10.2/){$_ = ">imasillymarker_M1GA0025150.10.2056617";}} print $_;' < pig_80k_marker_sites.fa > pig_80k_marker_sites.fixed.fa

module load samtools bwa
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a GCA_002844635.1_USMARCv1.0_genomic.fna -p pig_80k_marker_sites.fixed.fa -o pig_80k_marker_sites.fixed.USMARC.snps
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a GCF_000003025.6_Sscrofa11.1_genomic.fna -p pig_80k_marker_sites.fixed.fa -o pig_80k_marker_sites.fixed.ROSLIN.snps
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/alignAndOrderSnpProbes.pl -a ssc_10.2_ncbi.reference.fa -p pig_80k_marker_sites.fixed.fa -o pig_80k_marker_sites.fixed.SS10.snps
```

## Correlation analysis

```bash
# Generating relative marker order for Spearman rank correlation
grep '>' pig_60k_marker_sites.fa | perl -e '%h; while(<>){chomp; $_ =~ s/^>//; @b = split(/\./); $h{$b[1]}->{$b[2]} = $_;} foreach my $chr (sort {$a <=> $b} keys(%h)){$c = 1; foreach my $pos (sort{$a <=> $b} keys(%{$h{$chr}})){print "$h{$chr}->{$pos}\t$c\n"; $c++;}}' > pig_60k_marker_order.tab
cat pig_60k_marker_sites.ROSLIN.snps.tab | perl -e '%h; while(<>){chomp; @b = split(/\t/); $h{$b[1]}->{$b[2]} = "$b[0].$b[4].$b[5]";} $c = 1; foreach my $chr (sort {$a <=> $b} keys(%h)){foreach my $pos (sort{$a <=> $b} keys(%{$h{$chr}})){print "$h{$chr}->{$pos}\t$c\n"; $c++;}}' > pig_60k_marker_ROSLIN_snporder.tab
```

```R
library(dplyr)
original <- read.delim("pig_60k_marker_order.tab", header=FALSE)
colnames(original) <- c("PROBE", "ORank")
roslin <- read.delim("pig_60k_marker_ROSLIN_snporder.tab", header=FALSE)
colnames(roslin) <- c("PROBE", "RRank")

combined <- left_join(original, roslin, by="PROBE")
combined.ranks <- select(combined, -PROBE)
cor(combined.ranks, method="kendall", use="pairwise")
          ORank     RRank
ORank 1.0000000 0.7762156
RRank 0.7762156 1.0000000
```

I need to generate a better comparison input data between methods first. I wrote a python script to join together the output for better comparisons.

```bash
source activate python3
python3 ../binaries/python_toolchain/utils/tabFileLeftJoinTable.py -f pig_60k_marker_sites.ROSLIN.snps.tab -f pig_60k_marker_sites.SS10.snps.tab -f pig_60k_marker_sites.USMARC.snps.tab -c 0 -c 4 -c 5 -m 1 -m 2 -o pig_60k_marker_sites.left_join.snps.tab
for i in axiom_pigHD pig_60k_marker_sites; do echo $i; python ../binaries/python_toolchain/utils/tabFileLeftJoinTable.py -f $i.ROSLIN.snps.tab -f $i.SS10.snps.tab -f $i.USMARC.snps.tab -c 0 -c 4 -c 5 -m 1 -m 2 -o $i.left_join.snps.tab; done
python ../binaries/python_toolchain/utils/tabFileLeftJoinTable.py -f pig_80k_marker_sites.ROSLIN.snps.tab -f pig_80k_marker_sites.SS10.snps.tab -f pig_80k_marker_sites.USMARC.snps.tab -c 0 -c 4 -c 5 -m 1 -m 2 -o pig_80k_marker_sites.left_join.snps.tab

#### Redo of 80k markers!
python3 ~/python_toolchain/utils/tabFileLeftJoinTable.py -f pig_80k_marker_sites.fixed.ROSLIN.snps.tab -f pig_80k_marker_sites.fixed.SS10.snps.tab -f pig_80k_marker_sites.fixed.USMARC.snps.tab -c 0 -c 4 -c 5 -m 1 -m 2 -o pig_80k_marker_sites.fixed.left_join.snps.tab

cat *left_join*.tab | sort -k2,2n -k3,3n > all_pig_snps.left_join.tab

### Redo post 80k marker mappings!
cat axiom_pigHD.left_join.snps.tab pig_60k_marker_sites.left_join.snps.tab pig_80k_marker_sites.fixed.left_join.snps.tab | sort -k2,2n -k3,3n > all_pig_snps.fixed.left_join.tab

# Removing all previously unplaced chromosomes from the chips -- we will take care of them later
perl -lane 'if($F[1] == 0 || $F[1] =~ /^GL/){next;}else{print $_;}' < all_pig_snps.left_join.tab > all_pig_snps.left_join.filt.tab

# Now to use a custom script to rank things
perl generate_probe_rank_order.pl all_pig_snps.left_join.filt.tab all_pig_snps.left_join.rankings

### Redo post 80k marker mappings!
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/pigGenomeSNPSortRankOrder.pl all_pig_snps.fixed.left_join.filt.tab all_pig_snps.fixed.left_join.rankings
```

```R
data <- read.delim("pig_60k_marker_sites.left_join.snps.tab", header=FALSE)
colnames(data) <- c("Probe", "chr", "pos", "Rchr", "Rpos", "Schr", "Spos", "Mchr", "Mpos")
summary(data[data$chr == 0 & (data$Rchr == "*" | data$Schr == "*" | data$Mchr == "*"),])
          Rchr          Rpos                    Schr          Spos
 NC_010448.4:151   Min.   :        1   *          :993   Min.   :        1
 *          :139   1st Qu.:  4922196   NC_010445.3: 21   1st Qu.:        1
 NC_010447.5: 90   Median : 32875492   NC_010454.3: 19   Median :        1
 NC_010455.5: 85   Mean   : 54318145   NC_010456.4: 19   Mean   : 10424795
 NC_010445.4: 71   3rd Qu.: 89301410   NC_010455.4: 16   3rd Qu.:        1
 NC_010456.5: 63   Max.   :272967198   NC_010447.4: 15   Max.   :314046274
 (Other)    :621                       (Other)    :137
         Mchr          Mpos
 *         :243   Min.   :        1
 CM009091.1:125   1st Qu.:    38301
 CM009090.1: 74   Median : 20546664
 CM009097.1: 61   Mean   : 46837785
 CM009098.1: 60   3rd Qu.: 81869888
 CM009088.1: 50   Max.   :280093874
 (Other)   :607

# So out of 12,627 original unmapped probes (chr0 from the probe list), only 139 are still unmapped on Roslin, 993 on SS10 and 243 on USMARC

# now calculating the rank correlation
colnames(fullranks) <- c("probe", "ref", "Roslin", "SS10", "MARC")
summary(fullranks)
                          probe             ref             Roslin
 ALGA0000009-0_B_R_2031571545:     1   Min.   :     1   Min.   :    -1
 ALGA0000014-0_B_F_1895350382:     1   1st Qu.:199375   1st Qu.:148967
 ALGA0000021-0_T_F_2031571546:     1   Median :374102   Median :330088
 ALGA0000022-0_T_F_2031571547:     1   Mean   :371185   Mean   :330630
 ALGA0000046-0_B_R_2031571548:     1   3rd Qu.:548830   3rd Qu.:511052
 ALGA0000047-0_B_R_2031571549:     1   Max.   :723558   Max.   :693837
 (Other)                     :698906
      SS10             MARC
 Min.   :    -1   Min.   :    -1
 1st Qu.:153291   1st Qu.:139751
 Median :334182   Median :320540
 Mean   :334685   Mean   :321980
 3rd Qu.:515404   3rd Qu.:502325
 Max.   :698121   Max.   :684107

# Now the rank order tests
cor.test( ~ ref + Roslin, data=fullranks[,c(-1)], method="spearman")
data:  ref and Roslin
S = 5.9408e+15, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.8955929

cor.test( ~ ref + SS10, data=fullranks[,c(-1)], method="spearman")
data:  ref and SS10
S = 6.605e+15, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
    rho
0.88392

cor.test( ~ ref + MARC, data=fullranks[,c(-1)], method="spearman")
data:  ref and MARC
S = 9.6871e+15, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
     rho
0.829754

# Now let's try to plot them as a linear increasing plot
library(ggplot2)
library(dplyr)
library(reshape)

fullranks.melt <- melt(fullranks, id = c("probe", "ref"))
png(file="combined_pig_snp_rank_order.png", width=1200, height=1200)
ggplot(fullranks.melt, aes(x=ref, y=value, color=variable)) + geom_point(shape=1) + scale_color_brewer(palette="Dark2") + theme_bw() + theme(text = element_text(size=20)) + labs(title = "Assembly SNP rank concordance vs reported positions", x = "Reported positions (rank order)", y = "Mapped positions (rank order)")
dev.off()
```

### Bugfixes

It turns out that the 80k chip data had an extra column that was mucking up my analysis! I've written a script to correct this and now the numbers look far better.

```R
fullranks <- read.delim("all_pig_snps.fixed.left_join.rankings", header=FALSE)
colnames(fullranks) <- c("probe", "ref", "Roslin", "SS10", "MARC")

library(ggplot2)
library(dplyr)
library(reshape)

cor.test( ~ ref + Roslin, data=fullranks[,c(-1)], method="spearman")
sample estimates:
      rho 
0.8888981 

cor.test( ~ ref + SS10, data=fullranks[,c(-1)], method="spearman")
sample estimates:
      rho 
0.8846437

cor.test( ~ ref + MARC, data=fullranks[,c(-1)], method="spearman")
sample estimates:
      rho 
0.8126023


fullranks.melt <- melt(fullranks, id = c("probe", "ref"))
png(file="combined_pig_snp_rank_order_fixed.png", width=1200, height=1200)
ggplot(fullranks.melt, aes(x=ref, y=value, color=variable)) + geom_point(shape=1) + scale_color_brewer(palette="Dark2") + theme_bw() + theme(text = element_text(size=20)) + labs(title = "Assembly SNP rank concordance vs reported positions", x = "Reported positions (rank order)", y = "Mapped positions (rank order)")
dev.off()
```

Let's try a "leave-one-out" approach to see if we can correct the rank order artifacts that show a "v" in the plot.

```bash
cat pig_60k_marker_sites.left_join.snps.tab pig_80k_marker_sites.fixed.left_join.snps.tab | sort -k2,2n -k3,3n > 60k_80k_snps.fixed.left_join.tab
cat axiom_pigHD.left_join.snps.tab pig_80k_marker_sites.fixed.left_join.snps.tab | sort -k2,2n -k3,3n > axiom_80k_snps.fixed.left_join.tab
cat axiom_pigHD.left_join.snps.tab pig_60k_marker_sites.left_join.snps.tab | sort -k2,2n -k3,3n > axiom_60k_snps.fixed.left_join.tab

for i in 60k_80k_snps.fixed axiom_80k_snps.fixed axiom_60k_snps.fixed; do echo $i; perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/pigGenomeSNPSortRankOrder.pl $i.left_join.tab $i.left_join.rankings; done
```

```R
testRankings <- function(tab, name){
ranks <- read.delim(tab, header=FALSE)
colnames(ranks) <- c("probe", "ref", "Roslin", "SS10", "MARC")
print(cor.test( ~ ref + Roslin, data=ranks[,c(-1)], method="spearman"))
print(cor.test( ~ ref + SS10, data=ranks[,c(-1)], method="spearman"))
print(cor.test( ~ ref + MARC, data=ranks[,c(-1)], method="spearman"))
ranks.melt <- melt(ranks, id=c("probe", "ref"))
png(file=name, width=1200, height=1200)
print(ggplot(ranks.melt, aes(x=ref,y=value, color=variable)) + geom_point(shape=1) + scale_color_brewer(palette="Dark2") + theme_bw() + theme(text = element_text(size=20)) + labs(title = "Assembly SNP rank concordance vs reported positions", x = "Reported positions (rank order)", y = "Mapped positions (rank order)"))
dev.off()
}

# Just a way to automate all of the tests
for(x in c("axiom_60k_snps.fixed.left_join", "axiom_80k_snps.fixed.left_join", "60k_80k_snps.fixed.left_join")){tab <- paste0(x, ".rankings"); name <- paste0(x, ".png"); testRankings(tab, name);}


# And a color-blind appropriate palette
fullranks.melt <- melt(fullranks, id=c("probe", "ref"))
png(file="combined_pig_snp_rank_order_fixed.cb.png", width=1200, height=1200)
ggplot(fullranks.melt, aes(x=ref, y=value, color=variable)) + geom_point(shape=1, alpha = 0.4) + scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9")) + theme_bw() + theme(text = element_text(size=20)) + labs(title = "Assembly SNP rank concordance vs reported positions", x = "Reported positions (rank order)", y = "Mapped positions (rank order)")
dev.off()
```

| Combination | Roslin | SS10 | USMARC |
| :--- | ---: | ---: | ---:|
|axiom-60k | 0.6452 | 0.5157 | 0.5844 |
|axiom-80k | 0.6515 | 0.5270 | 0.5901 |
|60k-80k   | 0.5762 | 0.4241 | 0.5283 |

The "v-shaped" artifact is caused by the 60k SNP marker reported positions.

#### Trying to get reordered USMarc chromosomes

```bash
minimap2 -x asm10 GCF_000003025.6_Sscrofa11.1_genomic.fna GCA_002844635.1_USMARCv1.0_genomic.fna > roslin_to_usmarc_minimap.paf
perl -lane 'if($F[0] =~ /^NPJO/){next;} print "$F[0];$F[5];$F[4]";' < roslin_to_usmarc_minimap.paf | python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f stdin -c 0 > roslin_to_usmarc_orientation.tab

```

#### Chromosomes that are inverted in USMarc

| chr name | chr num |
|:---        |---:|
| CM009086.1 |  1 |
| CM009090.1 |  5 |
| CM009091.1 |  6 |
| CM009092.1 |  7 |
| CM009093.1 |  8 |
| CM009094.1 |  9 |
| CM009095.1 |  10|
| CM009096.1 |  11|
| CM009098.1 |  13|
| CM009099.1 |  14|
| CM009100.1 |  15|
| CM009101.1 |  16|
 

```bash
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/assembly_scripts/pigGenomeSNPSortRankOrder.pl all_pig_snps.fixed.left_join.filt.tab all_pig_snps.fixed.left_join.correct.rankings correct

# SNP spacing
## Roslin
perl -e '%keys; while(<>){chomp; @s = split(/\t/); if($s[5] eq "*"){next;}else{push(@{$keys{$s[5]}}, $s[6]);}} $sum = 0; foreach my $c (keys(%keys)){my @d = @{$keys{$c}}; @j = sort{$a <=> $b} @d; for($x = 1; $x < scalar(@j); $x++){print ($j[$x] - $j[$x - 1]); print "\n";}}' < all_pig_snps.fixed.left_join.filt.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   692920
Sum:    2711433846
Minimum 0
Maximum 2563404
Average 3913.054676
Median  2722
Standard Deviation      14990.315604
Mode(Highest Distributed Value) 0

## SS10
perl -e '%keys; while(<>){chomp; @s = split(/\t/); if($s[3] eq "*"){next;}else{push(@{$keys{$s[3]}}, $s[4]);}} $sum = 0; foreach my $c (keys(%keys)){my @d = @{$keys{$c}}; @j = sort{$a <=> $b} @d; for($x = 1; $x < scalar(@j); $x++){print ($j[$x] - $j[$x - 1]); print "\n";}}' < all_pig_snps.fixed.left_join.filt.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   690761
Sum:    2469665425
Minimum 0
Maximum 14310303
Average 3575.282080
Median  2647
Standard Deviation      26587.544585
Mode(Highest Distributed Value) 0

## USMARC
perl -e '%keys; while(<>){chomp; @s = split(/\t/); if($s[7] eq "*"){next;}else{push(@{$keys{$s[7]}}, $s[8]);}} $sum = 0; foreach my $c (keys(%keys)){my @d = @{$keys{$c}}; @j = sort{$a <=> $b} @d; for($x = 1; $x < scalar(@j); $x++){print ($j[$x] - $j[$x - 1]); print "\n";}}' < all_pig_snps.fixed.left_join.filt.tab | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   677172
Sum:    2492657220
Minimum 0
Maximum 29053888
Average 3680.980932
Median  2657
Standard Deviation      39304.300111
Mode(Highest Distributed Value) 0
```

#### SNP spacing on major chromosomes

| Assembly | Mapped SNPs | Max Distance | Avg Dis. | Median Dis. | Stdev Dis. | 
| :------- | ----------: | -----------: | --------:| ----------: | ---------: |
| Ss10     | 690,761     | 14,310,303   | 3,575.28 | 2,647       | 26,587.54  |
| Ss11(Ros)| 692,920     | 2,563,404    | 3,913.05 | 2,722       | 14,990.32  |
| MARCv1   | 677,172     | 29,053,888   | 3,680.98 | 2,657       | 39,304.30  |


```R
fullranks <- read.delim("all_pig_snps.fixed.left_join.correct.rankings", header=FALSE)
colnames(fullranks) <- c("probe", "ref", "Sscrofa11.1", "Sscrofa10.2", "USMARCv1.0")
cor.test(fullranks$ref, fullranks$USMARCv1.0, method="spearman")
sample estimates:
      rho 
0.8188725
fullranks.melt <- melt(fullranks, id=c("probe", "ref"))
png(file="combined_pig_snp_rank_order_fixed.corrected.cb.png", width=1200, height=1200)
ggplot(fullranks.melt, aes(x=ref, y=value, color=variable)) + geom_point(shape=1, alpha = 0.4) + scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9")) + theme_bw() + theme(text = element_text(size=20)) + labs(title = "Assembly SNP rank concordance vs reported positions", x = "Reported positions (rank order)", y = "Mapped positions (rank order)")
dev.off()
```