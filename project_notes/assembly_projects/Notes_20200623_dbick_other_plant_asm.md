# Plant assembly repeatmasking
---
*6/23/2020*


## Table of contents


## Preparing the assemblies

Just some quick processing to get things ready for downstream analysis.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/plants

```bash
module load samtools bwa

for i in *.fasta; do echo $i; sbatch --nodes=1 --mem=15000 --ntasks-per-node=2 -p priority -q msn --wrap="samtools faidx $i; bwa index $i;"; done

wc -l *.fai
 14808 HiFi_canu.contigs.fasta.fai
  6598 two_cell_CO46_hifiasm.a.fasta.fai
  5040 two_cell_CO46_hifiasm.p.fasta.fai
 26446 total

# Looks like hifiasm separated out two haplotypes
```

## Repeat analysis

OK, I'm going to queue up Repeatmasker to try to call repeats on each assembly. I want to make sure that consecutive RM tasks don't collide in the same directory as well! I have no information on which species of plant this is, so I will use the generic repbase libraries.

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/plants

```bash
# First, making sure that I have the species in the database
/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/util/queryTaxonomyDatabase.pl -species "Camelina sativa" -taxDBFile /lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/Libraries/taxonomy.dat

RepeatMasker Taxonomy Database Utility
======================================
Species = Camelina sativa
Lineage = Camelina sativa
...
# It's there

# Now to see what's in the repeat database
/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/util/queryRepeatDatabase.pl -species "Camelina sativa" -stat
...
176 ancestral and ubiquitous sequence(s) with a total length of 51081 bp
0 lineage specific sequence(s) with a total length of 0 bp
# Did not see satellite repeats -- not sure if they're in there! Will have to test

# Checking arabidopsis instead
/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/util/queryRepeatDatabase.pl -species "arabidopsis" -stat
...
176 ancestral and ubiquitous sequence(s) with a total length of 51081 bp
558 lineage specific sequence(s) with a total length of 1701864 bp
# That actually had the satellite sequences marked. Might be better to use arabidopsis then?

# HiCanu first
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 70 -species arabidopsis -no_is -gff HiFi_canu.contigs.fasta"

# Now the other two
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 70 -species arabidopsis -no_is -gff two_cell_CO46_hifiasm.a.fasta"
sbatch --nodes=1 --mem=100000 --ntasks-per-node=70 -p priority -q msn --wrap="/lustre/project/rumen_longread_metagenome_assembly/binaries/RepeatMasker/RepeatMasker -pa 70 -species arabidopsis -no_is -gff two_cell_CO46_hifiasm.p.fasta"

# Quick table comparison
python3 summarizeRepeatMaskerTables.py rm_table_summary.tab HiFi_canu.contigs.fasta.tbl two_cell_CO46_hifiasm.a.fasta.tbl two_cell_CO46_hifiasm.p.fasta.tbl

# hacking a basepair summary of these elements compared against canu
# "small" == "small rna"
perl -lane 'if($_ =~ /^Class/){print "Category\tCanu\tHifiASM\tPCanu"; next;} @fsegs = split(/,/, $F[1]); @hsegs = split(/,/, $F[2]); @gsegs = split(/,/, $F[3]); $c = $fsegs[1]; $h = $hsegs[1] + $gsegs[1]; $p = sprintf("%0.2f", ($h > 0)? $c / $h : 0.0); print "$F[0]\t$c\t$h\t$p";' < rm_table_summary.tab
```

#### Centromere analysis

Now I need to calculate how frequently the centromeric repeat is present and then see how many contigs are mostly derived from centromeres

> Ceres: /lustre/project/rumen_longread_metagenome_assembly/plants

```bash
 module load bedtools/2.25.0
# Transforming them all to bed format
for i in *.out; do echo $i; perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < $i > $i.bed; done

# Now to generate tables
for i in *.bed; do echo $i; python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 5 -o $i.tab -d '\t' -m; done

# And checking on centromeric repeats
grep Satellite *.bed.tab
HiFi_canu.contigs.fasta.out.bed.tab:|Satellite          |   7468|
HiFi_canu.contigs.fasta.out.bed.tab:|Satellite/centr    |      1|
two_cell_CO46_hifiasm.a.fasta.out.bed.tab:|Satellite          |   1527|
two_cell_CO46_hifiasm.p.fasta.out.bed.tab:|Satellite          |   6123|
two_cell_CO46_hifiasm.p.fasta.out.bed.tab:|Satellite/centr    |      1|

#Only one in each?? Let's see
grep 'Satellite/centr' *.out.bed
HiFi_canu.contigs.fasta.out.bed:tig00000701     6049515 6049575 -       COLAR12 Satellite/centr 178
two_cell_CO46_hifiasm.p.fasta.out.bed:utg000002l        11316286        11316346        -       COLAR12 Satellite/centr       178

# Both are smack-dab in the center of the contigs!
# From looking at two examples of inserts that Tim gave (utg000061 and utg0000431), it looks like a repeating span of LTRs and Satellites in both regions. Maybe heterochromatin involved in the Centromere? Let's grep things out and make a bed track
cat HiFi_canu.contigs.fasta.out.bed |perl -lane 'if($F[5] =~ /Satellite/){print "$F[0]\t$F[1]\t$F[2]";}' | bedtools merge -i stdin > HiFi_canu.contigs.fasta.out.bed.sats
cat two_cell_CO46_hifiasm.p.fasta.out.bed | perl -lane 'if($F[5] =~ /Satellite/){print "$F[0]\t$F[1]\t$F[2]";}' | bedtools merge -i stdin > two_cell_CO46_hifiasm.p.fasta.out.bed.sats

# Getting 1 mb window counts of satellite regions
bedtools makewindows -g HiFi_canu.contigs.fasta.fai -w 1000000 | bedtools intersect -a stdin -b HiFi_canu.contigs.fasta.out.bed.sats -wa -c | perl -lane 'if($F[-1] == 0){next;}else{print $_;}'
bedtools makewindows -g two_cell_CO46_hifiasm.p.fasta.fai -w 1000000 | bedtools intersect -a stdin -b two_cell_CO46_hifiasm.p.fasta.out.bed.sats -wa -c | perl -lane 'if($F[-1] == 0){next;}else{print $_;}' | less
```