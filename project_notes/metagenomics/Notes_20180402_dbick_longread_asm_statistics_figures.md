# USDA Pacbio assembly statistics comparison
---
*4/2/2018*

These are my notes on generating figures, tables and other statistics for our assembly manuscript.

## Table of Contents


## Analysis ideas

#### Cluster assessment

We have so many different clustering and binning methods that it would make sense to generate some metric by which the clustering algorithms compare to one another. One method that I found online was the Rand Index, which rates the group membership of clusters. If I calculate the Rand Index between tetra-nucleotide-based binning methods and then by Hi-C clustering, and then run an ANOVA on the clusters based on their GC + read depth (from illumina reads), then that might result in interesting data points


## Full-assembly comparison alignments

I want to generate alignment comparisons between the long-read and short read assemblies from our sample. I think that this comparison will show the unique regions to further interrogate and may reveal interesting differences between the technologies. My inspiration is from [this paper](http://genome.cshlp.org/content/25/4/534.full).

I want to separate our "good consensus microbes" from the other assemblies. There were about 16-20 clusters from the phase proxi-meta analysis that were "good," so let's start there. I want to separate segments based on > 5000 bp (chosen because of mashmap  and > 97% sequence identity at first to see how many bases overlapped.

> assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/assembly_align_comp

```bash
# Illumina as reference
/mnt/nfs/nfs2/bickhart-users/binaries/mashmap-Linux64-v2.0/mashmap -r ../illumina_usda_accumulated/mick_megahit_final_full.fasta -q ../pacbio_pilon_accumulated/total_contigs.fasta --perc_identity 97 -t 10 -o usda_ilmn_vs_pacbio_mashmap_5k_97.out
# PacBio as reference
/mnt/nfs/nfs2/bickhart-users/binaries/mashmap-Linux64-v2.0/mashmap -q ../illumina_usda_accumulated/mick_megahit_final_full.fasta -r ../pacbio_pilon_accumulated/total_contigs.fasta --perc_identity 97 -t 10 -o usda_pacbio_vs_ilmn_mashmap_5k_97.out

python3 /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/metagenomics/mashmapCalculationCondenser.py -i usda_ilmn_vs_pacbio_mashmap_5k_97.out -o usda_ilmn_vs_pacbio_mashmap_5k_97.c
python3 /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/metagenomics/mashmapCalculationCondenser.py -i usda_pacbio_vs_ilmn_mashmap_5k_97.out -o usda_pacbio_vs_ilmn_mashmap_5k_97.c

# Now to separate out the queries that had few mappings
perl -lane 'if($F[2] / $F[1] > 0.75){print $F[0];}' < usda_ilmn_vs_pacbio_mashmap_5k_97.c.segs.tab > pacbio_query_gt_75_unmapped.list
perl -lane 'if($F[2] / $F[1] > 0.75){print $F[0];}' < usda_pacbio_vs_ilmn_mashmap_5k_97.c.segs.tab > ilmn_query_gt_75_unmapped.list

perl -lane 'print $F[0]' < usda_ilmn_vs_pacbio_mashmap_5k_97.c.segs.tab > pacbio_query_full_map.list
perl -lane 'print $F[0]' < usda_pacbio_vs_ilmn_mashmap_5k_97.c.segs.tab > ilmn_query_full_map.list

perl -lane 'print $F[0]' < ../pacbio_pilon_accumulated/total_contigs.fasta.fai > pacbio_fullfasta_contigs.list
perl -lane 'print $F[0]' < ../illumina_usda_accumulated/mick_megahit_final_full.fasta.fai > ilmn_fullfasta_contigs.list

wc -l *.list
 2182263 ilmn_fullfasta_contigs.list
   32352 ilmn_query_full_map.list
     603 ilmn_query_gt_75_unmapped.list
   53098 pacbio_fullfasta_contigs.list
   37622 pacbio_query_full_map.list
     910 pacbio_query_gt_75_unmapped.list

# Pacbio stats first
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl pacbio_fullfasta_contigs.list pacbio_query_full_map.list pacbio_query_gt_75_unmapped.list
File Number 1: pacbio_fullfasta_contigs.list
File Number 2: pacbio_query_full_map.list
File Number 3: pacbio_query_gt_75_unmapped.list
Set     Count
1       15476
1;2     36712
1;2;3   910

# so 16386 pacbio contigs with few or no Ilumina contig mappings

# Illumina stats next
perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl ilmn_fullfasta_contigs.list ilmn_query_full_map.list ilmn_query_gt_75_unmapped.list
File Number 1: ilmn_fullfasta_contigs.list
File Number 2: ilmn_query_full_map.list
File Number 3: ilmn_query_gt_75_unmapped.list
Set     Count
1       2149911
1;2     31749
1;2;3   603

# Let's calculate it based on Mbp
# Aligned pacbio read Mbp
perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); if($F[2] < 500){$F[2] = 0;} $c += $F[1] - $F[2];}print "$c\n";' < usda_ilmn_vs_pacbio_mashmap_5k_97.c.segs.tab
506976006
perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); if($F[2] < 500){$F[2] = 0;} $c += $F[1] - $F[2];}print "$c\n";' < usda_pacbio_vs_ilmn_mashmap_5k_97.c.segs.tab
277258885

perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); $c += $F[1];} print "$c\n";' < ../illumina_usda_accumulated/mick_megahit_final_full.fasta.fai
5111042186
perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); $c += $F[1];} print "$c\n";' < ../pacbio_pilon_accumulated/total_contigs.fasta.fai
845929544

# So the pacbio contigs are more repetitive apparently. They capture double the amount of bases found in the Illumina contigs
# Here is how the venn would look:
# PBunique -> shared -> ILunique
# 338953538	-> 277258885 -> 4833783301


# Now to calculate the number of covered bases in each assembly
# Illumina first
samtools depth -aa aligns/USDA/USDA.sorted.merged.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); $h{$s[2]} += 1;} foreach my $k (sort {$a <=> $b} keys(%h)){print "$k\t$h{$k}\n";}' > ilmn_depth_histogram.tab
# now pacbio (pre-pilon)
cd ../pacbio_usda_retry_pilon
samtools depth -aa aligns/USDA/USDA.sorted.merged.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); $h{$s[2]} += 1;} foreach my $k (sort {$a <=> $b} keys(%h)){print "$k\t$h{$k}\n";}' > pacbio_prepilon_depth_histogram.tab
```

#### Contig lengths by SuperKingdom

I want to try to separate out some information based on taxonomic assignment of contigs. Namely, what is the distribution of lengths of contigs generated by each technology.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated

```bash
pwd
	/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated
perl -lane 'if($F[0] =~ /\#/){next;} else{print "$F[0]\t$F[1]\t$F[5]";}' < usda_pacbio_supkingdom.usda_pacbio_pilon_blobplot.blobDB.table.txt > usda_pacbio_supkingdom.contig.lens.tab

perl -lane 'if($F[0] =~ /\#/){next;} else{print "$F[0]\t$F[1]\t$F[5]";}' < ../illumina_usda_accumulated/usda_illumina_dontuse_supkingdom.mick_megahit_illumina_blobplot.blobDB.table.txt > usda_illumina_dontuse_supkingdom.contig.lens.tab
```

```R
library(ggplot2)
data <- read.delim("usda_pacbio_supkingdom.contig.lens.tab", header=FALSE)

colnames(data) <- c("Contig", "Len", "Kingdom")
pdf(file="usda_pacbio_supkingdom_contiglens.pdf", useDingbats=FALSE)
ggplot(data, aes(x = Len, fill = Kingdom)) + geom_density(alpha = 0.3) + scale_x_log10(breaks = c(1000, 10000, 25000, 100000, 250000))
dev.off()

# Let's see if separate boxplots make the data easier to interrogate
# Boxplots didn't work, let's try faceting the density plots
pdf(file="usda_pacbio_supkingdom_contiglens_boxplot.pdf", useDingbats=FALSE)
ggplot(data, aes(x= Len, fill=Kingdom)) + geom_density() + facet_grid(Kingdom ~ .) + xlim(c(0,200000))
dev.off()

# This gives me an idea: Let's plot both the pacbio and illumina contig lengths per kingdom instead, and facet the data
library(dplyr)
illumina <- read.delim("usda_illumina_dontuse_supkingdom.contig.lens.tab", header=FALSE)
colnames(illumina) <- c("Contig", "Len", "Kingdom")

data <- mutate(data, Tech = c("PacBio"))
illumina <- mutate(illumina, Tech = c("Illumina"))
total <- bind_rows(illumina, data)

total$Kingdom <- as.factor(total$Kingdom)
total$Tech <- as.factor(total$Tech)

pdf(file="usda_ilmn_pacbio_supkingdom_contiglens_density.pdf", useDingbats=FALSE)
p <- ggplot(total, aes(x= Len)) + geom_density(aes(fill=Tech), alpha=0.5) + xlim(c(0,150000)) + facet_grid(Kingdom ~ .) + ylim(c(0, 0.0006))
p
dev.off()
```

Now to generate the merged file. 

## CRISPR and defense system detection

In this section I hope to identify and annotate elements of the CRISPR-CAS9 and other viral defense systems, and I hope to also identify spacer sequence that is suitable for analysis.

Important things to check:
* Number of CRISPR systems in pacbio vs illumina
* Number of repeats per array in pacbio vs illumina
* Bias of array elements in pacbio vs illumina


Let's start by searching for spacers from CRISPR.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/crispr_spacers

```bash
# Using Piler-cr
/mnt/nfs/nfs2/bickhart-users/binaries/pilercr1.06/pilercr -in ../illumina_usda_accumulated/mick_megahit_final_full.fasta -out mick_megahit_final_full.crispr.spacer.out
   0 s       3 Mb (  0%)  Reading sequences from ../illumina_usda_accumulated/mick_megahit_final_full.fasta

*** FATAL ERROR ***  ReadMFA: buffer too small
# It's probably because Mick's asm doesn't have normal newline delimits in the fasta sequence. I can reformat but it may take some time

samtools faidx mick_megahit_final_full.fasta k127_21 k127_37 k127_39 k127_42 k127_43 k127_46 k127_48 k127_49 k127_51 k127_60 k127_62 k127_63 k127_68 k127_72 k127_73 k127_75 k127_81 k127_83 k127_86 k127_93 k127_94 k127_95 k127_98 k127_99 k127_100 k127_103 k127_106 k127_115 k127_121 k127_122 > test.fasta

/mnt/nfs/nfs2/bickhart-users/binaries/pilercr1.06/pilercr -in test.fasta -out test.crispr.out
# That was it. The fasta file just needed to be properly formatted

# Testing it with the pacbio data
/mnt/nfs/nfs2/bickhart-users/binaries/pilercr1.06/pilercr -in ../pacbio_pilon_accumulated/pacbio_pilon_unique_contigs.fasta -out pacbio_pilon_unique_contigs.crispr.spacer.out
# Could not create parameters errors. I think that the program was designed to work under 32 bit constraints
# Rather than rewrite it, let's try to parse it out on the pre-accumulated pilon corrected pacbio data
mkdir pacbio
for i in ../pacbio_usda_pilon/*.fasta; do name=`basename $i | cut -d'.' -f1,2`; echo $name; /mnt/nfs/nfs2/bickhart-users/binaries/pilercr1.06/pilercr -in $i -out pacbio/$name.crispr.spacer.out -minarray 3 -minrepeat 23 -maxrepeat 50 -minspacer 26 -maxspacer 50 -seq pacbio/$name.crispr.fasta -trimseqs; done

# It's annoying, but the software was developed long ago for another type of computer.
find -maxdepth 1 -size +0 -print -type f | grep 'fasta' | wc -l
70 <- fastas that contained CRISPR repeats according to pilecr

# The file format is not easy to parse. Trying to parse it with a one-liner, but I'll have to conver this to a script later
for i in `find -maxdepth 1 -size +0 -print -type f | grep 'fasta'`; do crispr=`basename $i | cut -d'.' -f1,2,3`; echo $crispr; perl -e '$name; @spacers; @pos; $start = 0; $int = 0; while(<>){chomp; if($start == 0){if($_ =~ DETAIL){$start = 1; next;}}else{if($_ =~ /^>(.+)_pilon/){$name = $1; next;} if($_ =~ /^=/){$int = ($int == 1)? 0 : 1; next;} if($int){$_ =~ s/^\s+//; @s = split(/\s+/); push(@pos, $s[0]); push(@spacers, $s[-1]);}}} for($x = 0; $x < scalar(@spacers); $x++){print ">$name\_$pos[$x]\n$spacers[$x]\n";}' < $crispr.spacer.out > $crispr.spacers.fasta; done

wc -l *.spacers.fasta
3902 total <- 1951 spacers in total so far

# Let's combine them all
cat *.spacers.fasta > ../pacbio_pilon_unique_first.crispr.spacers.fasta
# And let's gather basic statistics (ie. contig and count of spacers)
grep '>' pacbio_pilon_unique_first.crispr.spacers.fasta | perl -e '%h; while(<>){chomp; $_ =~ /^>(.+)\_/; $h{$1} += 1;} foreach my $k (sort{$a cmp $b} keys(%h)){print "$k\t$h{$k}\n";}' > pacbio_pilon_unique_first.crispr.spacer.contig.cnt

wc -l pacbio_pilon_unique_first.crispr.spacer.contig.cnt
70 pacbio_pilon_unique_first.crispr.spacer.contig.cnt

# OOPs! There was a flaw in my script! Creating a formal script to process the data instead
for i in `find -maxdepth 1 -size +0 -print -type f | grep 'fasta'`; do crispr=`basename $i | cut -d'.' -f1,2,3`; echo $crispr; perl ~/sperl/metagenomics_scripts/pilecr_processing_script.pl -i $crispr.spacer.out -o $crispr.r; done
cat *.r*.fasta > ../pacbio_pilon_unique_first.crispr.spacers.fasta
cat *.r*.stats > ../pacbio_pilon_unique_first.crispr.spacer.contig.cnt
wc -l ../pacbio_pilon_unique_first.crispr.spacer.contig.cnt
90 ../pacbio_pilon_unique_first.crispr.spacer.contig.cnt  <- that's more like it

# Now to do this on the illumina data using my pipeline script
mkdir illumina
cd illumina
perl ~/sperl/metagenomics_scripts/pilecr_processing_script.pl -f ../../illumina_usda_accumulated/mick_megahit_final_full.fasta -o mick_megahit_crispr

# Now let's clean up the junk
rm *temp*.fa
find . -empty -type f -delete

### NOTE: this was on my local Virtualbox desktop
# I am having some success with minced, but it processes the fasta in serial and takes about one second per fasta entry -- which would be over 1 million seconds!
# Actually, I was wrong! It processed everything in an hour and found 5,000+ unique spacers
crisprviz.sh -f rumen_illuminaR3PCRFree_megahit.final.contigs.fa -m 15 -n 60
grep 'Repeats:' rumen_illuminaR3PCRFree_megahit.final.contigs.fa.crisprs | perl -e '$c = 0; while(<STDIN>){chomp; $_ =~ /Repeats: (\d+)/; $c += $1 - 1;} print "$c\n";'
5202

# And for the first-round pilon-corrected asm
crisprviz.sh -f pacbio_pilon_unique_contigs.fasta -m 15 -n 60
grep 'Repeats:' pacbio_pilon_unique_contigs.fasta.crisprs | perl -e '$c = 0; while(<STDIN>){chomp; $_ =~ /Repeats: (\d+)/; $c += $1 - 1;} print "$c\n";'
1516

# Damn, it's fewer! But these are only the clustered pacbio pilon CRISPR arrays and they're not validated
# Creating stat tables for plotting in the /figure_drafts/raw_stats folder
grep '>' ../../../pilot_project/assemblies/rumen_illuminaR3PCRFree_megahit.final.contigs.fa_spacers.fa | perl -e '%h; while(<>){chomp; $_ =~ s/>//g; @s = split(/_/); $h{"$s[0]_$s[1]"} += 1;} foreach my $k (keys(%h)){print "$k\t$h{$k}\n";}' > mick_megahit_crisprviz_spacer_counts.tab
grep '>' ../../../pilot_project/assemblies/pacbio_pilon/pacbio_pilon_unique_contigs.fasta_spacers.fa | perl -e '%h; while(<>){chomp; $_ =~ s/>//g; @s = split(/_/); $h{"$s[0]_$s[1]"} += 1;} foreach my $k (keys(%h)){print "$k\t$h{$k}\n";}' > pacbio_pilon_firsttry_spacers.fa
```

#### Alignment and association of spacer elements

I want to get a sense of what types of contigs the spacers align to.

```bash
# Pacbio first
bwa mem ../pacbio_pilon_accumulated/total_contigs.fasta pacbio_pilon_unique_first.crispr.spacers.fasta > pacbio_pilon_unique_first.crispr.spacers.sam

# Finding non-self hits -- note: there are lots of alternative alignments that I'm not checking!
perl -lane 'if($F[0] =~ /^@/){next;} @nsegs = split(/_/, $F[0]); if($nsegs[0] ne $F[2] && $F[2] ne "*"){print "$F[0]\t$F[2]";}' < pacbio_pilon_unique_first.crispr.spacers.sam > pacbio_pilon_unique_first.crispr.spacers.assignable.tab

cat pacbio_pilon_unique_first.crispr.spacers.assignable.tab | cut -f2 | sort | uniq | xargs -I {} grep {} ../pacbio_pilon_accumulated/usda_pacbio_supkingdom.usda_pacbio_pilon_blobplot.blobDB.table.txt | perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 5
Entry   Count
Bacteria        48
Eukaryota       1
Viruses 3
no-hit  2

# looks like many could be other CRISPR arrays on other bacteria, or may be due to false positive CRISPR assignment

```


## ORF detection and comparison

I will be running Prodigal as a means of finding ORFs in our contigs.

#### Illumina first

> Assemble2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
sbatch --nodes=1 --mem=100000 --ntasks-per-node=2 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/Prodigal/prodigal -a mick_megahit_final_full.prod.prottrans -c -d mick_megahit_final_full.prod.genenuc -f gff -i mick_megahit_final_full.fasta -o mick_megahit_final_full.prod.out -p meta"

```

## Preliminary stat drawings

Let's make some small figures and tables to show how the data looks.

#### Venn of bp overlap between datasets

```R
library(VennDiagram)
pdf(file="pacbio_illumina_megabase_overlap.pdf", useDingbats=FALSE)
# Using MBp overlap values drawn from my first mashmap -- LIKELY TO CHANGE!!!!!!!!!
draw.pairwise.venn(area2=5111, area1=845, cross.area=277, category=c("PacBio", "Illumina"), fill=c("blue", "red"), cat.pos=c(-10,10), cex=c(2,2,2), cat.cex=c(3,3))
dev.copy2pdf(file="pacbio_illumina_megabase_overlap.pdf", useDingbats=FALSE)

```

#### GC content and contig length histograms

I need basic stat drawings to fill in some of our multi-part figures in the schema.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated

```bash
python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/sequenceData/calcGCcontentFasta.py  -f pacbio_pilon_unique_contigs.fasta -o pacbio_pilon_unique_contigs.gc -t 10

python /mnt/nfs/nfs2/bickhart-users/binaries/python_toolchain/sequenceData/calcGCcontentFasta.py  -f ../illumina_usda_accumulated/mick_megahit_final_full.fasta -o mick_megahit_final_full.gc -t 10
perl -lane 'print "$F[0]\t$F[-1]";' < mick_megahit_final_full.gc > mick_megahit_final_full.simpnames.gc
```

And now for the plots. I'll keep it simple so that we have a placeholder for future refinement.

```R
library(ggplot2)
pacbio.gc <- read.delim("pacbio_pilon_unique_contigs.gc", header = FALSE)
colnames(pacbio.gc) <- c("Contig", "GC Perc")

illumina.gc <- read.delim("mick_megahit_final_full.simpnames.gc", header = FALSE)
colnames(illumina.gc) <- c("Contig", "GC Perc")

ks.test(pacbio.gc$GC_Perc, illumina.gc$GC_Perc)

        Two-sample Kolmogorov-Smirnov test

data:  pacbio.gc$GC_Perc and illumina.gc$GC_Perc
D = 0.35114, p-value < 2.2e-16
alternative hypothesis: two-sided

Warning message:
In ks.test(pacbio.gc$GC_Perc, illumina.gc$GC_Perc) :
  p-value will be approximate in the presence of ties

library(dplyr)
pacbio.gc <- mutate(pacbio.gc, Tech = c("PacBio"))
pacbio.gc$Tech <- as.factor(pacbio.gc$Tech)
illumina.gc <- mutate(illumina.gc, Tech = c("Illumina"))
illumina.gc$Tech <- as.factor(illumina.gc$Tech)

total.gc <- bind_rows(illumina.gc, pacbio.gc)
total.gc$Tech <- as.factor(total.gc$Tech)

pdf(file="ilmn_pacbio_usda_gcperc_boxplot.pdf", useDingbats=FALSE)
ggplot(total.gc, aes(x=Tech, y=GC_Perc, fill=Tech)) + geom_boxplot() + ylab("Average GC percentage per Contig")
dev.off()

# In this folder I have a preloaded data from my taxonomic facet_grid plots
ks.test(total[total$Tech == "PacBio","Len"], total[total$Tech == "Illumina","Len"])

        Two-sample Kolmogorov-Smirnov test

data:  total[total$Tech == "PacBio", "Len"] and total[total$Tech == "Illumina", "Len"]
D = 0.83783, p-value < 2.2e-16
alternative hypothesis: two-sided

Warning message:
In ks.test(total[total$Tech == "PacBio", "Len"], total[total$Tech ==  :
  p-value will be approximate in the presence of ties


pdf(file="ilmn_pacbio_usda_ctglen_boxplot.pdf", useDingbats=FALSE)
ggplot(total, aes(x=Tech, y=Len, fill=Tech)) + geom_boxplot() + ylab("Contig length (bp)") + scale_y_continuous(trans = 'log10')
dev.off()

p <- ggplot(total, aes(x=Tech, y=Len, fill=Tech)) + geom_boxplot() + ylab("Log10 Contig length (bp)") + scale_y_continuous(trans = 'log10')
ggsave("ilmn_pacbio_usda_ctglen_boxplot.png", plot = p)
```

#### AlignQC charts and graphs

I am going to generate some AlignQC data on the pilon-corrected reads and data.


#### ANVIO statistics

After some incredibly painful workarounds, I was finally able to create an anvio virtual environment by using anaconda3. The virtual environment is here:

> /mnt/nfs/nfs2/bickhart-users/binaries/virtual_env/anvio

Now to try to generate some anvio data on our pacbio contigs.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated

```bash
source activate /mnt/nfs/nfs2/bickhart-users/binaries/virtual_env/anvio

anvi-gen-contigs-database -f pacbio_pilon_unique_contigs.fasta -o pacbio_pilon_unique_contigs.anvio.db -n 'USDA PacBio Pilon Database'

#HMM checking
anvi-run-hmms -c pacbio_pilon_unique_contigs.anvio.db

# I need to run this once to set up the NCBI cog database
anvi-setup-ncbi-cogs --cog-data-dir /mnt/nfs/nfs2/bickhart-users/binaries/anvio_ncbi_cogs -T 10

# COG annotation
anvi-run-ncbi-cogs -c pacbio_pilon_unique_contigs.anvio.db -T 10 --cog-data-dir /mnt/nfs/nfs2/bickhart-users/binaries/anvio_ncbi_cogs

# Bam read profiling
anvi-profile -c pacbio_pilon_unique_contigs.anvio.db -i aligns/USDA/USDA.sorted.merged.bam --output-dir anvio_pacbio_profile --sample-name 'PacBioPilon' -T 10
```

And for the Illumina data.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
source activate /mnt/nfs/nfs2/bickhart-users/binaries/virtual_env/anvio

perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/metagenomics_scripts/pilecr_processing_script.pl -f mick_megahit_final_full.fasta -d 2000 -o mick_megahit_final_full

anvi-gen-contigs-database -f mick_megahit_final_full.rfmt.fa -o mick_megahit_final_full.rfmt.anvio.db -n 'USDA Illumina Database'

anvi-run-ncbi-cogs -c mick_megahit_final_full.rfmt.anvio.db -T 10 --cog-data-dir /mnt/nfs/nfs2/bickhart-users/binaries/anvio_ncbi_cogs

anvi-profile -c mick_megahit_final_full.rfmt.anvio.db -i aligns/USDA/USDA.sorted.merged.bam --output-dir anvio_illumina_profile --sample-name 'IlluminaMega' -T 10
```

#### Crisprviz spacer length plots

```R
library(dplyr)
# This will need to be adjusted with the pilecr data, but let's generate some stats now
pilon <- read.delim("pacbio_pilon_firsttry_spacers.fa", header=FALSE)
illumina <- read.delim("mick_megahit_crisprviz_spacer_counts.tab", header=FALSE)

colnames(pilon) <- c("Contig", "Count")
colnames(illumina) <- c("Contig", "Count")

pilon <- mutate(pilon, Dataset=c("PacBio"))
illumina <- mutate(illumina, Dataset=c("Illumina"))

total <- bind_rows(pilon, illumina)
total$Dataset <- as.factor(total$Dataset)

ks.test(pilon$Count, illumina$Count)

	Two-sample Kolmogorov-Smirnov test

data:  pilon$Count and illumina$Count
D = 0.25786, p-value = 2.822e-06
alternative hypothesis: two-sided

Warning message:
In ks.test(pilon$Count, illumina$Count) :
  p-value will be approximate in the presence of ties

> summary(pilon$Count)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   3.00    5.00    7.00   12.43   16.00  115.00 
> summary(illumina$Count)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  3.000   3.000   5.000   8.714   9.000 112.000 


plot(ecdf(x = pilon$Count), main = "ECDF of PacBio and Illumina CRISPR spacer counts per array", col = "blue")
lines(ecdf(x = illumina$Count), col = "red")
legend("center", c("Illumina", "PacBio"), col = c("red", "blue"), pch = c(20))
dev.copy2pdf(file="pacbio_illumina_crispr_ecdf.pdf", useDingbats=FALSE)

library(ggplot2)
ggplot(total, aes(x=Dataset, y=Count, fill=Dataset)) + geom_boxplot() + ylab("Count of CRISPR spacers per array")
dev.copy2pdf(file="pacbio_illumina_crispr_spacers_boxplot.pdf", useDingbats=FALSE)
```