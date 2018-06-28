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
```

## Alignment

```bash
# Aligning to USMarc
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCA_002844635.1_USMARCv1.0_genomic.fna -p $i.fa -o $i.USMARC.snps; done

# Aligning to ss11.1
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a GCF_000003025.6_Sscrofa11.1_genomic.fna -p $i.fa -o $i.ROSLIN.snps; done

# Finally, aligning to ss10.2
for i in ggp_marker_sites pig_60k_marker_sites pig_80k_marker_sites; do echo $i; perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a ssc_10.2_ncbi.reference.fa -p $i.fa -o $i.SS10.snps; done
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