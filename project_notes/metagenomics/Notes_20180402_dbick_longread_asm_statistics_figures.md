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

# Violin plot
pdf("rumen_contig_gc_violin_plot.pdf", useDingbats=FALSE)
p <- ggplot(total.gc, aes(x=Tech, y=GC_Perc, fill=Tech)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white")
p + labs(x="Technology", y="Avg GC ratio") + scale_fill_brewer(palette="Dark2") + theme_classic()
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

# Read depth profile merger
anvi-merge -o SAMPLES-MERGED -c pacbio_pilon_unique_contigs.anvio.db --sample-name pacbio_pilon /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated/anvio_pacbio_profile/PROFILE.db
# NOTE: if you try to merge only one RD profile, the program exits. Apparently it wants more data!

# I can import the phase cluster results as separate bins now
# Generating the simple tab delimited format they want
for i in ../pacbio_usda_clusters/*.fai; do echo $i; fname=`basename $i`; perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); $ARGV[1] =~ s/\.fasta\.fai//; while(<IN>){chomp; @s = split(/\t/); print "$s[0]\t$ARGV[1]\n";} close IN;' $i $fname >> pacbio_hic_clusters.tab; done

# OK, I got an error about the duplicate assignment of bins. Let me keep only the first observed bin in the file for each contig
## NOTE: Anvio hates non-underscore characters. I had to replace them in the cluster names ##
perl -e '%h; while(<>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){next;} else{$s[1] =~ s/[.-]/_/g; print join("\t", @s); print "\n"; $h{$s[0]} = 1;}}' < pacbio_hic_clusters.tab > pacbio_hic_clusters.firstobs.tab

wc -l pacbio_hic_clusters.tab pacbio_hic_clusters.firstobs.tab
  53222 pacbio_hic_clusters.tab
  53098 pacbio_hic_clusters.firstobs.tab

# There weren't too many. Oh well, continuing!
anvi-import-collection -c pacbio_pilon_unique_contigs.anvio.db -p /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated/anvio_pacbio_profile/PROFILE.db -C "ProxiMeta"  --contigs-mode pacbio_hic_clusters.firstobs.tab

# Finally, summarizing before using interative output
anvi-summarize -p /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated/anvio_pacbio_profile/PROFILE.db -c pacbio_pilon_unique_contigs.anvio.db -C "ProxiMeta" -o PACBIO-PRELIM

# That crashed with an error. The problem is that even though I have clustering via ProxiMeta, it was done only on a contig-by-contig basis. I would need more alignment files on my contigs to proceed.
# Maybe I can run alignments using Hess et al and some others?
ls /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/datasources/SRR094* > ../hess_rumen_fiber_sequence.tab
perl -lane 'print "$_\tHESS\tHESS";' < ../hess_rumen_fiber_sequence.tab > ../temp
mv ../temp ../hess_rumen_fiber_sequence.tab

perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b hess -t ../hess_rumen_fiber_sequence.tab -f pacbio_pilon_unique_contigs.fasta -m
# ARG! SRA screwed up with the read names! I need to reformat
perl -lane 'system("sbatch reformat_sra_files.pl $F[0] $F[1] /mnt/nfs/nfs1/derek.bickhart/metagenomics/reformated_sra/");' < hess_rumen_fiber_sequence.tab

ls /mnt/nfs/nfs1/derek.bickhart/metagenomics/reformated_sra/*.fq* > hess_rumen_fiber_sequence.tab
perl -lane 'print "$_\tHESS\tHESS";' < hess_rumen_fiber_sequence.tab > temp
mv temp hess_rumen_fiber_sequence.tab

sleep 1h; for i in SRR094926_1.fq SRR094926_2.fq SRR094437_1.fq SRR094437_2.fq; do echo $i; pigz /mnt/nfs/nfs1/derek.bickhart/metagenomics/reformated_sra/$i; done; perl ~/sperl/sequence_data_pipeline/generateAlignSlurmScripts.pl -b hess -t ../hess_rumen_fiber_sequence.tab -f pacbio_pilon_unique_contigs.fasta -m
rm hess/HESS/*sorted.bam*

## Now, I have an idea to make this complete in a reasonable amount of time: let's reprofile without indvidual base composition metrics ##
# Hess first
anvi-profile -i hess/HESS/HESS.sorted.merged.bam -c pacbio_pilon_unique_contigs.anvio.db -T 25 --skip-SNV-profiling --sample-name PbHessNoSNV -o PbHessNoSNVProfile
# Now the illumina data
anvi-profile -i aligns/USDA/USDA.sorted.merged.bam -c pacbio_pilon_unique_contigs.anvio.db -T 25 --skip-SNV-profiling --sample-name PbIlmnNoSNV -o PbIlmnNoSNVProfile
# And merging:
anvi-merge -o NoSNV-merged -c pacbio_pilon_unique_contigs.anvio.db PbHessNoSNVProfile/PROFILE.db PbIlmnNoSNVProfile/PROFILE.db
# Adding Proximeta bins
anvi-import-collection -p NoSNV-merged/PROFILE.db -c pacbio_pilon_unique_contigs.anvio.db -C "ProxiMeta" --contigs-mode pacbio_hic_clusters.firstobs.tab
# Now summarizing...
anvi-summarize -p NoSNV-merged/PROFILE.db -c pacbio_pilon_unique_contigs.anvio.db -C "CONCOCT" -o PbCombNoSNV
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

## Methylation check

I think that DNA methylation is an important metric to use to bin our samples. There has been some work done on detecting Methylated sites in both Nanopore and Pacbio reads using HMMs and Neural networks. I want to see if I can run some software on our corrected sites to detect methylation signatures.

#### Nanopore methylation detection

While several groups have published theory and datasets, only a handful have tried to detect methyl-adenine reads from their data. Let's use the existing Nanopore data from my sequencing of the sample to try to classify reads. Here's my conceptualized order of operations:

1. Determine methyl adenine sites on assembled contigs from Serge's assembly of our data
2. Map those contigs back to the PacBio data and determine meythlyation status
3. Eventually, determine Pacbio contig methylation status and compare back to Nanopore calls

I am first going to try to use [signalAlign](https://github.com/ArtRand/signalAlign). 

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/nanopore

```bash
# Downloading data
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 0BxbRXPCzWa5-Q1VMa1RLOTNhWm8 nanopore_asm_mhap.fasta.gz

sbatch /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 0BxbRXPCzWa5-SlpWWHVpNFp1TkE 20170802_YuAndMorrison6.tar.gz

# I need to basecall the fast5 files using the latest version of albacore. I setup a virtual env on the Agil cluster
# to do this.
source activate albacore

read_fast5_basecaller.py --flowcell FLO-MIN106 --kit SQK-LSK108 --barcoding --output_format fast5 --input /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/nanopore/fast5 --save_path basecalled --worker_threads 30 -r

virtualenv --python=python2.7 virtual_env/signalAlign
source virtual_env/signalAlign/bin/activate
# Apparently there's a big issue with signalAlign, I may need to try a different solution

# Trying tombo
conda install -c bioconda ont-tombo
# Conda couldn't install on the cluster due to the outdated version of R we have
pip install ont-tombo[full]

tombo resquiggle basecalled/ ../pacbio_pilon_accumulated/pacbio_pilon_unique_contigs.fasta --processes 20

tombo test_significance --fast5-basedirs basecalled --alternate-bases 5mC 6mA --statistics-file-basename YuAndMorrison6 --processes 30
```

#### PacBio methylation detection

I will try to use [mbin](https://github.com/fanglab/mbin) to classify some CCS reads from Cheryl's dataset and use them similar to my strategy with the Nanopore methylated reads.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_methylation

```bash
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1FmsVqMhUmA66TcyybbX2gxF-N6dYDWbR chery_rumen.css.tar.gz

# Damn! That tarball did not have the bam files I needed.
sbatch /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 0BxbRXPCzWa5-OGhweUVzMXdSZW8 01Aug17_USMARC_Sequel.tar.gz
tar -xvf 01Aug17_USMARC_Sequel.tar.gz

# OK, this one did!
## INSTALLING MBIN ##
# Now I need to set up a virtual environment
virtualenv --python=python2.7 virtual_env/mbin
source virtual_env/mbin/bin/activate
pip install mbin

# Adding bhtsne to PATH
pwd
/mnt/nfs/nfs2/bickhart-users/binaries/bhtsne

export PATH=$PATH:/mnt/nfs/nfs2/bickhart-users/binaries/bhtsne
# OK, the solution was to check the stacktrace error, pip uninstall the offending package, and then reinstall without the use of the cache dir (--no-cache-dir)


# I need to align non-adenine methylated bases to my dataset using the SMRT tools pipeline so that the IPD length gets carried over.
# This is part of the in initial "buildcontrols" step of the pipeline that is needed to interrogate the signal. It's also super painful and I may need to return to it later. 
```


#### Testing new CCS data

I want to see how similar our pacbio CCS data is to our previous sample.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics/pilot_project/new_ccs_pacbio

```bash
# Let's sketch this and check against the error corrected pacbio reads
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 10 -s 100000 -o ccs.fastq.21.100k ccs.fastq.gz

# Error corrected reads
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash sketch -p 10 -s 100000 -o rumen_pacbio_corrected.fasta.10.100k ../error_corrected_reads/rumen_pacbio_corrected.fasta.gz

# Mash gave me a warning about the kmer size for the pacbio error corrected reads. Still, the similarity is quite close
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash dist ccs.fastq.21.100k.msh rumen_pacbio_corrected.fasta.10.100k.msh
ccs.fastq.gz    ../error_corrected_reads/rumen_pacbio_corrected.fasta.gz        0.185048        0       1037/100000
```

## RePilon and finalization of PacBio dataset

I need to run two rounds of Pilon to check and confirm that the Pacbio assembly is correct. Then we'll decide on the final analysis steps and hand off the assembly to our working teams to generate some information.

I will be running this on Ceres for the first time. Since Ceres does not have Mouse installed, and because I'm worried about asking the team to install Perl modules, I rewrote my pipeline script for read alignment in python just for this purpose.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_first_assembly

```bash
# Preparing the input fasta file list
/home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_first_assembly

# Just checking to see if the path is recognized on external nodes
sbatch --nodes=1 --mem=100 --ntasks-per-node=1 --partition short --wrap="ls /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/pilot_project/illumina/*.gz"

# Looks like they're recognized. Let's test out the python script to see if it generates files correctly
module load python_3/gcc/64/3.6.2
python3 ~/rumen_longread_metagenome_assembly/binaries/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b align -t ../ymprep_run3_illumina_reads.tab -f usda_pacbio_multiround.asm.fasta -p medium -m

# Note, due to a bug at the time, the dependencies for the merger script didn't queue up. I need to run that separately after the alignment jobs finish
# Now to run the pilon correction
python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f USDA.sorted.merged.bam -g usda_pacbio_multiround.asm.fasta -o pilon -p short

# That didn't work out so well, since the scheduler has a job submission limit of 1000 jobs per user!
# Doing a manual edit to queue them all up
ls pilon/scripts/*.sh | perl -e '$c = 0; while(<STDIN>){chomp; open(IN, "< $_"); open(OUT, "> pilon/newscripts/pilon_$c\.sh"); $n = 0; whil
e(<IN>){if($n < 14){print OUT $_;} $n++;} close IN; for($x = 0; $x < 6; $x++){$d = <STDIN>; chomp $d; $n = 0; open(IN, "< $d"); while(<IN>){if($n
>= 12 && $n < 14){ print OUT $_;} $n++;} close IN;} $c++; close OUT;}'

# still too many scripts. Queueing up the last 100
for i in `seq 901 999` `seq 90 99` 9; do echo $i; sbatch pilon/newscripts/pilon_${i}.sh; done

# I have a simple script for combining all of the entries and removing the pilon tags
sbatch ../../consolidate_pilon_fastas.sh usda_pacbio_pilon_round1.fa
cp pacbio_first_assembly/pilon/usda_pacbio_pilon_round1.fa ./pacbio_first_pilon/
```

OK, now I can run the second round of pilon correction on the data.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_first_pilon

```bash
module load bwa/0.7.12
sbatch --nodes=1 --ntasks-per-node=2 --mem=6000 -p short --wrap="bwa index usda_pacbio_pilon_round1.fa; samtools faidx usda_pacbio_pilon_round1.fa"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b align -t ../ymprep_run3_illumina_reads.tab -f usda_pacbio_pilon_round1.fa -p medium -m

# The dependencies still didn't load for the merger step! Grr!
sbatch --dependency=afterany:225158:225159:225160:225161 align/USDA/scripts/samMerger_22929245490.sh

# I was finally able to fix the bugs in the pilon pipeline script
python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f align/USDA/USDA.sorted.merged.bam -g usda_pacbio_pilon_round1.fa
 -o pilon -p short -m

sbatch ../../consolidate_pilon_fastas.sh usda_pacbio_second_pilon_round.fa
samtools faidx usda_pacbio_second_pilon_round.fa
perl -e 'chomp(@ARGV); %h; foreach $a (@ARGV){open(IN, "< $a"); while(<IN>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]); } close IN; } foreach $c (keys(%h)){for($x = 1; $x < scalar(@{$h{$c}}); $x++){ $d = abs($h{$c}->[$x] - $h{$c}->[$x - 1]); if($d != 0){print "$c\t" . join("\t", @{$h{$c}}) . "\t$d\n";}}}' usda_pacbio_second_pilon_round.fa.fai ../usda_pacbio_pilon_round1.fa.fai > ../usda_pacbio_pilon_round2.changes.tab

# Damn! I corrected all bases instead of indels! Going to redo the correction
python3 ~/python_toolchain/sequenceData/slurmPilonCorrectionPipeline.py -f align/USDA/USDA.sorted.merged.bam -g usda_pacbio_pilon_round1.fa -o pilon2 -p short -m
sbatch ../../consolidate_pilon_fastas.sh usda_pacbio_second_pilon_indelsonly.fa

perl -e 'chomp(@ARGV); %h; foreach $a (@ARGV){open(IN, "< $a"); while(<IN>){chomp; @s = split(/\t/); push(@{$h{$s[0]}}, $s[1]); } close IN; } foreach $c (keys(%h)){for($x = 1; $x < scalar(@{$h{$c}}); $x++){ $d = abs($h{$c}->[$x] - $h{$c}->[$x - 1]); if($d != 0){print "$c\t" . join("\t", @{$h{$c}}) . "\t$d\n";}}}' usda_pacbio_second_pilon_indelsonly.fa.fai ../usda_pacbio_pilon_round1.fa.fai > ../usda_pacbio_pilon_round2_indelsonly.changes.tab

# Tabulating the results of the pilon correction
for i in pilon2/outLog/*.out; do echo $i; perl -e '%h; while($l = <>){chomp $l; if($l =~ /^Confirmed/){$b = <>; $s = <>; if($s =~ /^Large/){while(1){$s = <>; if($s =~ /^Large/){}else{last;}}} chomp $b; chomp $s; ($samp, $len) = $s =~ /(.+):\d+-(\d+)/; ($perc) = $l =~ /\((.+)\%\)/; ($snp, $amb, $ins, $del) = $b =~ /Found (\d+) snps; (\d+) ambiguous .* corrected (\d+) .* (\d+) small deletions .*/; $h{$samp} = [$len, $perc, $snp, $amb, $ins, $del];}} foreach $k (sort {$a cmp $b} keys(%h)){print "$k\t" . join("\t", @{$h{$k}}) . "\n";}' < $i >> pilon_round2_correction_stats.tab ; done
```

#### Anvio run on pacbio2 dataset

I want to quantify the improvement of the second round of pilon correction on the pacbio dataset. I will attempt to run Anvio on it with alignment of the Hess dataset and our dataset

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_usda_round2pilon

```bash
python /mnt/nfs/nfs2/bickhart-users/binaries/download_from_gdrive.py 1Ky0CWfmc6dwrrRUaAC8_9bWaFFqLqOuJ usda_pacbio_second_pilon_indelsonly.fa.gz
unpigz usda_pacbio_second_pilon_indelsonly.fa.gz

sbatch --mem=15000 --ntasks-per-node=1 --nodes=1 -p assemble3 -wrap="bwa index usda_pacbio_second_pilon_indelsonly.fa; module load samtools; samtools faidx usda_pacbio_second_pilon_indelsonly.fa"

```


#### Kmeans clustering to determine clustering averages

> Assembler2: 

```R
library(dplyr)
library(broom)
library(ggplot2)
data <- read.delim("usda_clusters_rg_counts.ext.full.tab", header=TRUE)
data <- data[data$ReadGroup == "USDA", c(1,2,4,5,6,7,8)]

data.filt <- mutate(data, Cov = Count * 150 / Len, Name = paste(Cluster, Contig, sep="_"))

# Removing high coverage outliers
data.final.filt <- filter(data.final, Cov < 137944)

data.matrix <- cbind(GC = data.final.filt$GC,Cov =  data.final.filt$Cov)
kclusts <- data.frame(k=2:10) %>% group_by(k) %>% do(kclust=kmeans(data.matrix, .$k))
clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))
assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], data.matrix))

png("usda_pacbio_kmeans_clustering.png")
p1 <- ggplot(assignments, aes(GC, Cov)) + geom_point(aes(color=.cluster)) + facet_wrap(~ k)
p1
dev.off()

png("usda_pacbio_kmeans_variance.png")
ggplot(clusterings, aes(k, tot.withinss)) + geom_line()
dev.off()

# And now, how to draw out summary stat information from the kmeans
five <- assignments[assignments$k == 5,]
five.means <- group_by(five, .cluster)

summarize(five.means, count= n(), gc = mean(GC), cov = mean(Cov), stdv = sd(Cov))
# A tibble: 5 x 5
  .cluster count    gc   cov  stdv
  <chr>    <int> <dbl> <dbl> <dbl>
1 1           46 0.422 25838  4843
2 2          119 0.476 11361  3081
3 3           12 0.297 59433  9749
4 4          670 0.508  2953  1375
5 5        42682 0.491   107   174

ten <- assignments[assignments$k == 10,]
ten.means <- group_by(ten, .cluster)
summarize(ten.means, count = n(), gc = mean(GC), cov = mean(Cov), stdv = sd(Cov))
# A tibble: 10 x 5
   .cluster count    gc     cov   stdv
   <chr>    <int> <dbl>   <dbl>  <dbl>
 1 1           52 0.472 14370   1816
 2 10          17 0.391 30979   3440
 3 2           29 0.440 22824   2342
 4 3          841 0.504  1316    351
 5 4        37364 0.491    56.4   49.9
 6 5           66 0.477  9053   1314
 7 6          317 0.508  2867    559
 8 7           12 0.297 59433   9749
 9 8         4712 0.491   373    148
10 9          119 0.496  5508    839

```

## Public dataset alignment and binning

I want to download and analyze publically accessible rumen WGS reads for further binning and classification of the rumen microbes in this project. I wanted to download everything possible in the SRA. Here was my command for finding the illumina-only, WGS-only datasets: **((cattle rumen) NOT "pcr"[Selection]) NOT "ls454"[Platform] **. Unfortunately, that wasn't enough, and I had to manually filter away the ION torrent, 454 and amplicon reads that remained. I also had to add the HESS dataset in (it was missing). Now I need to download the files sequentially on Ceres.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/public_datasets

```bash
# Downloading the files
cat cattle_rumen_illumina_run_sheet_sra.txt | perl -e '$h = <STDIN>; while(<STDIN>){chomp; @s = split(/\t/); print "$s[0]\n";}' | xargs -I {} sbatch --nodes=1 --tasks-per-node=1 --mem=500 -p short --wrap="fastq-dump.2 -I --gzip --split-files {}"

# Quite a few worked, but quite a few also failed! Perhaps it was a simultaneous query issue? I will have to ensure that all of the files are properly accounted for in the folder.
# Also, many of the files were corrupt

# An online git error shows that most of the modern SRA tools need to prefetch data and convert them inplace. It's a pain, but let me see what I can do. 
# In true NCBI fashion, it tries to download stuff to home even if you don't want it to!
rm ~/ncbi/public/sra/*
prefetch --max-size 50000000 ERR2282092
fastq-dump.2 -I --gzip --split-files ERR2282092

# This worked, but it took a long time! Let's script it so that I have 3 running processes to download all of the files
head -n 41 cattle_rumen_illumina_run_sheet_sra.txt | tail -n 40 | perl -lane 'print $F[0];' > download_set_one.list
head -n 81 cattle_rumen_illumina_run_sheet_sra.txt | tail -n 40 | perl -lane 'print $F[0];' > download_set_two.list
tail -n 40 cattle_rumen_illumina_run_sheet_sra.txt | perl -lane 'print $F[0];' > download_set_three.list

sbatch download_sra_dataset.sh download_set_one.list
sbatch download_sra_dataset.sh download_set_two.list
sbatch download_sra_dataset.sh download_set_three.list

cat remaining_set.list | xargs -I {} sbatch --nodes=1 --ntasks-per-node=2 --mem=2000 --time=12:00:00 --wrap="prefetch --max-size 45000000 {}; fastq-dump.2 -I --gzip --split-files {}"

# OK, now to generate a spreadsheet tab file for the projects
perl -ne '$d = "/scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly/sequence_data/public_datasets"; chomp; @F = split(/\t/); if($F[0] =~ /^Run/){next;} @files = `ls $F[0]*`; chomp(@files); if(scalar(@files) < 2){next;} print "$d/$files[0]\t$d/$files[1]\t$F[0]\t$F[21]\n";' < cattle_rumen_illumina_run_sheet_sra.txt > cattle_rumen_illumina_data_spreadsheet.tab

# We should be set for alignment... fingers crossed!
```

I am going to run alignments of all 109 of these datasets against each respective assembly for comparisons.

#### Pacbio

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon

```bash
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p short --wrap="module load bwa/gcc/64/0.7.12; bwa index usda_pacbio_second_pilon_indelsonly.fa"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b publicdb -t ../../../sequence_data/public_datasets/cattle_rumen_illumina_data_spreadsheet.tab -f usda_pacbio_second_pilon_indelsonly.fa -p short -m

# Dammit all! I forgot that the SRA creates read names with BWA errors!
# Remaking all of the fastq files:

perl -lane 'system("sbatch --nodes=1 --mem=6000 --ntasks-per-node=2 -p short --wrap=\"python3 ~/python_toolchain/sequenceData/fixSRAFastqFiles.py -f $F[0] -r $F[1] -o $F[2]\_r -l $F[2]\_r.log\"");' < cattle_rumen_illumina_data_spreadsheet.tab

# Some of the files were too big. Rather than restarting and wasting more time, I'm going to truncate them
for i in SRR094437_r SRR094926_r SRR094424_r ERR2282092_r ERR2530126_r ERR2027896_r; do echo $i; sbatch --nodes=1 -p short --mem=2000 --ntasks-per-node=4 --wrap="for b in ${i}.*; do echo $b; unpigz -c $b | wc -l; done"; done

# Damn, the compression corrupted the files. Let's reprocess them instead
for i in SRR094437 SRR094926 SRR094424 ERR2282092 ERR2530126 ERR2027896; do echo $i; sbatch --nodes=1 --mem=6000 --ntasks-per-node=2 -p long --wrap="python3 ~/python_toolchain/sequenceData/fixSRAFastqFiles.py -f ${i}_1.fastq.gz -r ${i}_2.fastq.gz -o ${i}_r -l ${i}_r.log"; done

# Updating the spreadsheet
perl -lane '$F[0] =~ s/1.fastq.gz/r.1.fq.gz/; $F[1] =~ s/2.fastq.gz/r.2.fq.gz/; if( -s $F[0] ||  -s $F[1]){print join("\t", @F);}' < cattle_rumen_illumina_data_spreadsheet.tab > cattle_rumen_illumina_data_spreadsheet.revised.tab

# Finger's crossed! Let's try again
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b publicdb -t ../../../sequence_data/public_datasets/cattle_rumen_illumina_data_spreadsheet.revised.tab -f usda_pacbio_second_pilon_indelsonly.fa -p short -m

# OK, now to generate anvio fasta profiles
sbatch --nodes=1 --ntasks-per-node=2 --mem=18000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-gen-contigs-database -f usda_pacbio_second_pilon_indelsonly.fa -o usda_pacbio_second_pilon_indelsonly.anvio.db"

# Generating NCBI cogs
sbatch --nodes=1 --ntasks-per-node=2 --mem=18000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-setup-ncbi-cogs --cog-data-dir ~/rumen_longread_metagenome_assembly/binaries/anvi_ncbi_cogs;"

#NOTE: the program didn't see ncbi tools installed, so it did not create a blast search DB for the files!

sbatch --nodes=1 --ntasks-per-node=6 --mem=20000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-run-hmms -c usda_pacbio_second_pilon_indelsonly.anvio.db --num-threads 6;"

sbatch --nodes=1 --ntasks-per-node=6 --mem=20000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2 blast/gcc/64/2.2.26; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-run-ncbi-cogs -c usda_pacbio_second_pilon_indelsonly.anvio.db --num-threads 6 --cog-data-dir ~/rumen_longread_metagenome_assembly/binaries/anvi_ncbi_cogs"
```

#### Illumina

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
sbatch --nodes=1 --mem=10000 --ntasks-per-node=2 -p short --wrap="module load bwa/gcc/64/0.7.12; bwa index illumina_megahit_final_contigs.fa"

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b publicdb -t ../../../sequence_data/public_datasets/cattle_rumen_illumina_data_spreadsheet.revised.tab -f illumina_megahit_final_contigs.fa -p short -m

# Now to reformat the fasta and generate anvio input data
sbatch --nodes=1 --ntasks-per-node=2 --mem=18000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-script-reformat-fasta -o illumina_megahit_final_contigs.reformated.fa --simplify-names illumina_megahit_final_contigs.fa"

# Ah, so anvio strips ALL contig names! That would be a huge pain to deal with later
# Going to do it the safe way with perl
perl -ne 'if($_ =~ /^>/){@b = split(/\s+/); print "$b[0]\n";}else{print $_;}' < illumina_megahit_final_contigs.fa > illumina_megahit_final_contigs.perl.fa

sbatch --nodes=1 --ntasks-per-node=2 --mem=18000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-gen-contigs-database -f illumina_megahit_final_contigs.perl.fa -o illumina_megahit_final_contigs.perl.anvio.db"

# It ran out of memory! Let's increase the memory size
sbatch --nodes=1 --ntasks-per-node=3 --mem=30000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-gen-contigs-database -f illumina_megahit_final_contigs.perl.fa -o illumina_megahit_final_contigs.perl.anvio.db"

# I think it failed again. Let's try this again, but this time with only the contigs that were above 2342 bp in length (the average; the median contig size was 1583 bp).
# I'm taking advantage of the fact that the megahit fasta is two lines per fasta entry
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while($n = <IN>){$s = <IN>; chomp $n; chomp $s; $b = substr($n, 1); if(exists($h{$b})){print "$n\n$s\n";}}' illumina_megahit_contigs.tosave.list illumina_megahit_final_contigs.perl.fa > illumina_megahit_final_contigs.gt2kb.fa

sbatch --nodes=1 --ntasks-per-node=3 --mem=30000 -p short --wrap="module load prodigalorffinder/gcc/64/2.6.3 samtools/gcc/64/1.4.1 hmmer3/gcc/64/3.1b2; source ~/rumen_longread_metagenome_assembly/binaries/virtual_envs/anvio-4/bin/activate; anvi-gen-contigs-database -f illumina_megahit_final_contigs.gt2kb.fa -o illumina_megahit_final_contigs.gt2kb.anvio.db"
```

#### Anvio processing


#### Metabat binning

I am going to generate Metabat bin files for each assembly based on the public db bams.

#### Pacbio
> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon

```bash
module load metabat/2.12.1

# I need to add in the YMprep3 data
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b publicdb -t ../ymprep_run3_illumina_reads.tab -p short -m -f usda_pacbio_second_pilon_indelsonly.fa

sbatch --nodes=1 --mem=25000 --ntasks-per-node=2 -p medium --wrap="jgi_summarize_bam_contig_depths --outputDepth usda_second_pilon_publicdb.depths.tab --pairedContigs usda_second_pilon_publicdb.paired.txt publicdb/*/*.merged.bam"

sbatch --nodes=1 --dependency=afterany:208209 --mem=50000 --ntasks-per-node=8 --wrap="metabat2 -i usda_pacbio_second_pilon_indelsonly.fa -a usda_second_pilon_publicdb.depths.tab -o usda_second_pilon_publicdb_metabat -t 8 -v"

# Now to organize the data and run checkM
mv usda_second_pilon_publicdb* ./metabat/
module load pplacer/v1.1.alpha19 hmmer3/gcc/64/3.1b2 prodigalorffinder/gcc/64/2.6.3
sbatch --nodes=1 --mem=45000 --ntasks-per-node=8 -p short --wrap="checkm lineage_wf -f metabat/CheckM.txt -t 8 -x fa metabat/ metabat/SCG"

```
#### Illumina
> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b publicdb -t ../ymprep_run3_illumina_reads.tab -p short -m -f illumina_megahit_final_contigs.fa

sbatch --nodes=1 --mem=25000 --ntasks-per-node=2 -p medium --wrap="jgi_summarize_bam_contig_depths --outputDepth illumina_megahit_contigs_publicdb.depths.tab --pairedContigs illumina_megahit_contigs_publicdb.paired.txt publicdb/*/*.merged.bam"

sbatch --nodes=1 --mem=70000 --ntasks-per-node=10 -p short --wrap="metabat2 -i illumina_megahit_final_contigs.perl.fa -a illumina_megahit_contigs_publicdb.depths.tab -o usda_illumina_publicdb_metabat -t 10 -v"

mkdir public_metabat
mv usda_illumina_publicdb_metabat*.fa ./public_metabat/
module load pplacer/v1.1.alpha19 hmmer3/gcc/64/3.1b2 prodigalorffinder/gcc/64/2.6.3 checkm/v1.0.11
sbatch --nodes=1 --mem=45000 --ntasks-per-node=8 -p short --wrap="checkm lineage_wf -f metabat/CheckM.txt -t 8 -x fa public_metabat/ public_metabat/SCG"
```


## Hybrid assembly comparison

Serge generated opera scaffolds from our Illumina megahit assembly and the raw pacbio reads. We need to validate the assembly, check to see what's been incorporated and how it all looks. First, let's start with an alignment comparison. Here are the assumptions:

* All of the contigs in the Illumina assembly 

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/opera_scaffolds

```bash
module load minimap2/2.6
sbatch --nodes=1 --ntasks-per-node=6 --mem=40000 -p short --wrap="minimap2 -x asm5 -t 6 rumen_opera_pacbio.scaffoldSeq.filled.fasta ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > opera_vs_pacbio_asm5_align.paf"
sbatch --nodes=1 --ntasks-per-node=6 --mem=40000 -p short --wrap="minimap2 -x asm5 -t 6 rumen_opera_pacbio.scaffoldSeq.filled.fasta ../illumina_megahit/illumina_megahit_final_contigs.fa > opera_vs_illumina_asm5_align.paf"

# Let's find out how many Illumina contigs are missing
srun --nodes=1 --ntasks-per-node=8 --mem=45000 --pty bash
perl -lane 'print $F[0];' < opera_vs_illumina_asm5_align.paf > opera_vs_illumina_asm5_align.illumina.ctgs.list
perl -lane 'print "$F[0]";' < ../illumina_megahit/illumina_megahit_final_contigs.fa.fai > illumina_megahit_contigs.list
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl opera_vs_illumina_asm5_align.illumina.ctgs.list illumina_megahit_contigs.list
	File Number 1: opera_vs_illumina_asm5_align.illumina.ctgs.list
	File Number 2: illumina_megahit_contigs.list
	Set     Count
	1;2     1712719
	2       469544
# Storing the results for later comparison
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o opera_vs_illumina_asm5_align.illumina.ctgs.list illumina_megahit_contigs.list
wc -l opera_vs_illumina_asm5_align.missingilmn.ctgs.list


# Now let's find out how many pacbio contigs are missing
perl -lane 'print "$F[0]";' < opera_vs_pacbio_asm5_align.paf > opera_vs_pacbio_asm5_align.pacbio.ctgs.list
perl -lane 'print "$F[0]";' < ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai > pacbio_pilon_contigs.list
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl opera_vs_pacbio_asm5_align.pacbio.ctgs.list pacbio_pilon_contigs.list
	File Number 1: opera_vs_pacbio_asm5_align.pacbio.ctgs.list
	File Number 2: pacbio_pilon_contigs.list
	Set     Count
	1;2     76401
	2       1269	<- interesting!!
# Storing for later comparisons again
mv group_2.txt opera_vs_pacbio_asm5_align.missingpacb.ctgs.list

# Finally, let's see what opera scaffolds are completely unique and are unmapped in each dataset
perl -lane 'print "$F[0]";' < rumen_opera_pacbio.scaffoldSeq.filled.fasta.fai > opera_contigs.list
perl -lane 'print "$F[5]";' < opera_vs_illumina_asm5_align.paf > opera_in_illumina_contigs.list
perl -lane 'print "$F[5]";' < opera_vs_pacbio_asm5_align.paf > opera_in_pacbio_contigs.list
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl opera_contigs.list opera_in_illumina_contigs.list opera_in_pacbio_contigs.list
	File Number 1: opera_contigs.list
	File Number 2: opera_in_illumina_contigs.list
	File Number 3: opera_in_pacbio_contigs.list
	Set     Count
	1;2     857761
	1;2;3   281175

# OK, that makes sense, they're contained entirely in the illumina contig list

###### Checking Pacbio-unique contigs ######
# Let's see what's unique to the pacbio dataset
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){print "$s[1]\n";}}' opera_vs_pacbio_asm5_align.missingpacb.ctgs.list ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   1269
Minimum 1047
Maximum 31005
Average 6299.130024
Median  6149
Standard Deviation      3620.589719
Mode(Highest Distributed Value) 1298

# OK, these tend to be the smaller contigs, but there half of this list is over 6kb in size
# To confirm, let's test this on the illumina contigs (just to make sure that the smaller contigs are overrepresented
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){print "$s[1]\n";}}' opera_vs_illumina_asm5_align.missingilmn.ctgs.list ../illumina_megahit/illumina_megahit_final_contigs.fa.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   469544
Minimum 1000
Maximum 1499
Average 1200.532195
Median  1178
Standard Deviation      140.280609
Mode(Highest Distributed Value) 1012

# Oh yeah, everything under 1500 as Serge predicted. I wonder if those illumina contigs have matches to the pacbio missing contigs though...
# Creating the fasta files
perl -e '@d; while(<>){chomp; push(@d, $_); if(scalar(@d > 100)){$wl = join(" ", @d); system("samtools faidx ../illumina_megahit/illumina_megahit_final_contigs.fa $wl >> opera_vs_illumina_asm5_align.missingilmn.ctgs.fa"); print "joining " . scalar(@d) . " contigs\n"; @d = ();}}' < opera_vs_illumina_asm5_align.missingilmn.ctgs.list

# That's taking a long, long time! Let's try aligning to the whole dataset in the meantime
sbatch --ntasks-per-node=6 --nodes=1 --mem=45000 -p short --wrap="minimap2 -x asm5 -t 6 ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa ../illumina_megahit/illumina_megahit_final_contigs.fa > pacbio_pilon_to_illumina_minimap2.paf"
sbatch --ntasks-per-node=6 --nodes=1 --mem=45000 -p short --wrap="minimap2 -x asm5 -t 6 ../illumina_megahit/illumina_megahit_final_contigs.fa ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa > illumina_to_pacbio_pilon_minimap2.paf"

perl -lane 'print $F[5];' < pacbio_pilon_to_illumina_minimap2.paf > pacbio_contigs_with_illumina_aligns.list
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl opera_vs_pacbio_asm5_align.missingpacb.ctgs.list pacbio_contigs_with_illumina_aligns.list pacbio_pilon_contigs.list
File Number 1: opera_vs_pacbio_asm5_align.missingpacb.ctgs.list
File Number 2: pacbio_contigs_with_illumina_aligns.list
File Number 3: pacbio_pilon_contigs.list
Set     Count
1;2;3   562
1;3     707
2;3     76070
3       331

# OK, so it looks like there are 707 contigs unaccounted for. Let's see what they look like
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o opera_vs_pacbio_asm5_align.missingpacb.ctgs.list pacbio_contigs_with_illumina_aligns.list pacbio_pilon_contigs.list
mv group_1_3.txt opera_vs_pacbio_no_ilm_aligns.ctgs.list

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){print "$s[1]\n";}}' opera_vs_pacbio_no_ilm_aligns.ctgs.list ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   707
Minimum 1047
Maximum 31005
Average 6094.702970
Median  5774
Standard Deviation      3455.965324
Mode(Highest Distributed Value) 3401

# Still some larger contigs present! Now let's look at it from the Illumina side
perl -lane 'print $F[5];' < illumina_to_pacbio_pilon_minimap2.paf > illumina_contigs_with_pacbio_aligns.list
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl opera_vs_illumina_asm5_align.missingilmn.ctgs.list illumina_contigs_with_pacbio_aligns.list illumina_megahit_contigs.list
File Number 1: opera_vs_illumina_asm5_align.missingilmn.ctgs.list
File Number 2: illumina_contigs_with_pacbio_aligns.list
File Number 3: illumina_megahit_contigs.list
Set     Count
1;2;3   26977
1;3     442567
2;3     488110
3       1224609

# There are still a ton of smaller contigs with no mapping sites in the pacbio dataset!
perl /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o opera_vs_illumina_asm5_align.missingilmn.ctgs.list illumina_contigs_with_pacbio_aligns.list illumina_megahit_contigs.list
mv group_1_3.txt opera_vs_illumina_no_pacb_aligns.ctgs.list
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){print "$s[1]\n";}}' opera_vs_illumina_no_pacb_aligns.ctgs.list ../illumina_megahit/illumina_megahit_final_contigs.fa.fai | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/statStd.pl
total   442567
Minimum 1000
Maximum 1499
Average 1200.050101
Median  1177
Standard Deviation      140.181060
Mode(Highest Distributed Value) 1012

# I want to generate a venn of all shared nucleotides, but it's difficult to use shared coordinates here. Let's use the Illumina megahit assembly here (it's a subset of all assemblies) and subtract downwards
perl -lane 'print "$F[0]\t$F[2]\t$F[3]";' < pacbio_pilon_to_illumina_minimap2.paf > pacbio_pilon_ilmn_coords.unsorted.bed
perl -lane 'print "$F[0]\t$F[2]\t$F[3]";' < opera_vs_illumina_asm5_align.paf > opera_ilmn_coords.unsorted.bed
perl -lane 'print "$F[0]\t1\t$F[1]";' < ../illumina_megahit/illumina_megahit_final_contigs.fa.fai > illmn_coords.unsorted.bed

# Pacbio bases not in the Opera or Illumina assemblies
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$c += $s[1];}} print "$c\n";' opera_vs_pacbio_no_ilm_aligns.ctgs.list ../pacbio_final_pilon/usda_pacbio_second_pilon_indelsonly.fa.fai
4308955

##### Are there any regions of the opera assembly that have no counterparts in the Illumina assembly? ####
module load bedtools/2.25.0
perl -lane 'print "$F[5]\t$F[7]\t$F[8]";' < opera_vs_illumina_asm5_align.paf | bedtools sort -i stdin | bedtools merge -i stdin > opera_vs_illumina_asm5_operaonly.bed
perl -lane 'print "$F[0]\t1\t$F[1]";' < rumen_opera_pacbio.scaffoldSeq.filled.fasta.fai | bedtools sort -i stdin > rumen_opera_pacbio.scaffoldSeq.filled.lens.bed

# Raw overlap and count
bedtools subtract -a rumen_opera_pacbio.scaffoldSeq.filled.lens.bed -b opera_vs_illumina_asm5_operaonly.bed | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       2149505
        Total Length:           55309422
        Length Average:         25.7312367265952
        Length Median:          10
        Length Stdev:           224.486278327532
        Smallest Length:        1
        Largest Length:         25809

# Most of these are tiny little junk segments though, what happens if we take only the 1kb+ fraction?
bedtools subtract -a rumen_opera_pacbio.scaffoldSeq.filled.lens.bed -b opera_vs_illumina_asm5_operaonly.bed | perl -lane 'if($F[2] - $F[1] > 1000){print $_;}' | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       11345
        Total Length:           28682425
        Length Average:         2528.19964742177
        Length Median:          1961
        Length Stdev:           1724.66393169366
        Smallest Length:        1001
        Largest Length:         25809

# OK, that's more convincing. Let's see if these regions have NO alignment to the pacbio assembly at all.
bedtools subtract -a rumen_opera_pacbio.scaffoldSeq.filled.lens.bed -b opera_vs_illumina_asm5_operaonly.bed | perl -lane 'if($F[2] - $F[1] > 1000){print $_;}' > opera_vs_illumina_asm5_operaonly.gt1kb.bed
perl -lane 'print "$F[5]\t$F[7]\t$F[8]";' < opera_vs_pacbio_asm5_align.paf | bedtools sort -i stdin | bedtools merge -i stdin > opera_vs_pacbio_asm5_operaonly.bed
bedtools subtract -a opera_vs_illumina_asm5_operaonly.gt1kb.bed -b opera_vs_pacbio_asm5_operaonly.bed |perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       7705
        Total Length:           21215715
        Length Average:         2753.49967553537
        Length Median:          2079
        Length Stdev:           1957.79685571838
        Smallest Length:        1001
        Largest Length:         25809

#### VENN information ####
# Interesting, so the assembled pacbio reads don't map very well to the pacbio scaffold sites
# I think that I have enough info to generate a venn
# Pacbio only: 4,308,955
# Illumina only: 
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$c += $s[1];}} print "$c\n";' opera_vs_illumina_asm5_align.missingilmn.ctgs.list ../illumina_megahit/illumina_megahit_final_contigs.fa.fai
	563,702,689
# Opera only: 21,215,715

# Pacbio area: 1076426244
# Illumina area: 5111042186
# Opera area: 3942687706

# Illumina + Pacbio area: 
bedtools sort -i pacbio_pilon_ilmn_coords.unsorted.bed | bedtools merge -i stdin | perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
1386185112  <- longer than the total pacbio assembly area. Why? Because of indels in the alignment

#### Digging into the Unique PacBio contigs ####
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o pacbio_pilon_contigs.list pacbio_contigs_no_illumina_aligns.list opera_vs_pacbio_asm5_align.missingpacb.ctgs.list
Set     Count
1       76070
1;2     331	<- no illumina aligns at all
1;2;3   707	<- no illumina aligns but used in scaffolding
1;3     562 <- extra contigs used in scaffolding

mv group_1_2.txt pacbio_contigs_no_illumina_no_opera.list
mv group_1_3.txt pacbio_contigs_with_illumina_and_opera.list
mv group_1_2_3.txt pacbio_contigs_no_illumina_and_opera.list
```

Because I could not install pybedtools on CERES, I transfered everything over to the AGIL server.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/opera_scaffolds

```bash
scp -pr derek.bickharhth@ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/opera_scaffolds/*.bed ./

bedtools sort -i illmn_coords.unsorted.bed | bedtools merge -i stdin > illmn_coords.sorted.bed
bedtools sort -i opera_ilmn_coords.unsorted.bed | bedtools merge -i stdin > opera_ilmn_coords.sorted.bed
bedtools sort -i pacbio_pilon_ilmn_coords.unsorted.bed | bedtools merge -i stdin > pacbio_pilon_ilmn_coords.sorted.bed
``` 

```python
import pybedtools
ilmn = pybedtools.BedTool('illmn_coords.sorted.bed')
opera = pybedtools.BedTool('opera_ilmn_coords.sorted.bed')
pacbio = pybedtools.BedTool('pacbio_pilon_ilmn_coords.sorted.bed')

# Total bp in Illumina
sum(map(len, ilmn))
	# 5108859923
# Total bp in Opera
sum(map(len, opera))
	# 4212340985
# Total bp in Pacbio
sum(map(len, pacbio))
	# 1386185112 (+ 4308955 = 1390494067)

# Illumina + opera
sum(map(len, (ilmn + opera)))
	# 4545626778
# Illumina + pacbio
sum(map(len, (ilmn + pacbio)))
	# 2406188564
# opera + pacbio
sum(map(len, (opera + pacbio)))
	# 

# Bases only in Illumina
sum(map(len, (ilmn - opera - pacbio)))
	# 520752537
# Bases in Illumina + opera, but not in pacbio
sum(map(len, (ilmn + opera - pacbio)))
	# 2181918822
# Bases in Illumina + pacbio, but not in opera
sum(map(len, (ilmn - opera + pacbio)))
	# 42480608
# Bases in all three
sum(map(len, (ilmn + opera + pacbio)))
	# 2363707956
# Bases in Pacbio + opera but not in Illumina
sum(map(len, (opera + pacbio - ilmn)))
	# 0
# Bases only in opera
sum(map(len, (opera - pacbio - ilmn)))
	# 0
# Bases only in Pacbio
sum(map(len, (pacbio - opera - ilmn)))
	# 0 (here, but not in reality), it should be 4308955

# damn, this is not a good utility because the intersection is dependent on the order of input
```

#### Detecting gaps in the assembly

Mick identified some absurd number of gaps in our Opera assembly! Let's see what's there.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/opera_scaffolds

```bash
module load java/64/1.8.0_121
sbatch --nodes=1 --ntasks-per-node=4 --mem=35000 -p short --wrap="java -jar ~/rumen_longread_metagenome_assembly/binaries/GetMaskBedFasta/store/GetMaskBedFasta.jar -f rumen_opera_pacbio.scaffoldSeq.filled.fasta -o rumen_opera_pacbio.scaffoldSeq.filled.gaps.bed -s rumen_opera_pacbio.scaffoldSeq.filled.gaps.stats"

```

## Final blobplots

I am going to include all of the public data in our final blobplots

#### Pacbio

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon

```bash
module load samtools/gcc/64/1.4.1
sbatch --nodes=1 --ntasks-per-node=5 --mem=35000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools map2cov -i usda_pacbio_second_pilon_indelsonly.fa -b publicdb/PRJEB10338/PRJEB10338.sorted.merged.bam -b publicdb/PRJEB21624/PRJEB21624.sorted.merged.bam -b publicdb/PRJEB8939/PRJEB8939.sorted.merged.bam -b publicdb/PRJNA214227/PRJNA214227.sorted.merged.bam -b publicdb/PRJNA255688/PRJNA255688.sorted.merged.bam -b publicdb/PRJNA270714/PRJNA270714.sorted.merged.bam -b publicdb/PRJNA280381/PRJNA280381.sorted.merged.bam -b publicdb/PRJNA291523/PRJNA291523.sorted.merged.bam -b publicdb/PRJNA366460/PRJNA366460.sorted.merged.bam -b publicdb/PRJNA366463/PRJNA366463.sorted.merged.bam -b publicdb/PRJNA366471/PRJNA366471.sorted.merged.bam -b publicdb/PRJNA366487/PRJNA366487.sorted.merged.bam -b publicdb/PRJNA366591/PRJNA366591.sorted.merged.bam -b publicdb/PRJNA366667/PRJNA366667.sorted.merged.bam -b publicdb/PRJNA366681/PRJNA366681.sorted.merged.bam -b publicdb/PRJNA398239/PRJNA398239.sorted.merged.bam -b publicdb/PRJNA60251/PRJNA60251.sorted.merged.bam -b publicdb/USDA/USDA.sorted.merged.bam -o usda_pacbio_secpilon"

sbatch --nodes=1 --ntasks-per-node=2 --mem=15000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools taxify -f usda_pacbio_second_pilon_indelsonly.diamond.hits -m ../uniprot_ref_proteomes.taxids -s 0 -t 2 -o usda_pacbio_second_pilon_indelsonly.uniprot"

mv usda_pacbio_second_pilon_indelsonly.uniprot.usda_pacbio_second_pilon_indelsonly.diamond.hits.taxified.out usda_pacbio_second_pilon_indelsonly.diamond.hits.uniprot.taxified.out
for i in *.cov; do echo -n "-c $i "; done; echo
sbatch --nodes=1 --ntasks-per-node=2 --mem=25000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools create -i usda_pacbio_second_pilon_indelsonly.fa -c usda_pacbio_secpilon.PRJEB10338.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJEB21624.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJEB8939.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA214227.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA255688.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA270714.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA280381.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA291523.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366460.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366463.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366471.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366487.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366591.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366667.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA366681.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA398239.sorted.merged.bam.cov -c usda_pacbio_secpilon.PRJNA60251.sorted.merged.bam.cov -c usda_pacbio_secpilon.USDA.sorted.merged.bam.cov -t usda_pacbio_second_pilon_indelsonly.diamond.hits.uniprot.taxified.out -o pacbio_secpilon_blobplot --db /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/blobtools/data/nodesDB.txt"

# I am using multiplots to try to separate out the different components here
# Phylum plot
sbatch --nodes=1 --ntasks-per-node=2 --mem=25000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools plot -i pacbio_secpilon_blobplot.blobDB.json --notitle --sort_first "no-hit,other,undef" -p 14 -o pacbio_secpilon_blobplot_phylum -m"

# Superkingdom plot
sbatch --nodes=1 --ntasks-per-node=2 --mem=25000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools plot -i pacbio_secpilon_blobplot.blobDB.json --notitle --sort_first "no-hit,other,undef" -p 14 -o pacbio_secpilon_blobplot_supkingdom -m -r superkingdom"

# Generating the full table
sbatch --nodes=1 --ntasks-per-node=2 --mem=5000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools view -i pacbio_secpilon_blobplot.blobDB.json --out pacbio_secpilon_blobplot_all -r all -b"

# reducing the verbosity of the filenames
for i in *.png; do echo $i; perl -e 'chomp(@ARGV); @asegs = split(/\./, $ARGV[0]); shift(@asegs);$temp = join(".", @asegs); $temp =~ s/\.blobDB\.json\.bestsum//g; $temp =~ s/\.p14\.span\.100//g; print "$temp\n"; system("mv $ARGV[0] $temp");' $i; done
mv pacbio_secpilon_blobplot_supkingdom.pacbio_secpilon_blobplot.blobDB.json.bestsum.superkingdom.p14.span.100.blobplot.multiplot.stats.txt pacbio_secpilon_blobplot.superkingdom.blobplot.multiplot.stats.txt
mv pacbio_secpilon_blobplot_phylum.pacbio_secpilon_blobplot.blobDB.json.bestsum.phylum.p14.span.100.blobplot.multiplot.stats.txt pacbio_secpilon_blobplot.phylum.blobplot.multiplot.stats.txt

# compressing and storing the files
mkdir pacbio_blobplot_superkingdom; ls *superkingdom.* | xargs -I {} mv {} ./pacbio_blobplot_superkingdom/; tar -czvf pacbio_blobplot_superkingdom.tar.gz ./pacbio_blobplot_superkingdom
mkdir pacbio_blobplot_phylum; ls *phylum.* | xargs -I {} mv {} ./pacbio_blobplot_phylum; tar -czvf pacbio_blobplot_phylum.tar.gz ./pacbio_blobplot_phylum

# I created a script to automate the steps above: 
sbatch ../run_blobplot_generation.sh genus pacbio pacbio_secpilon_blobplot.blobDB.json
```

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_final_pilon

```bash
sbatch --nodes=1 --mem=45000 --ntasks-per-node=3 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/bin/diamond blastx --query usda_pacbio_second_pilon_indelsonly.fa --db ../../diamond/uniprot_ref_proteosomes.diamond.dmnd --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 > usda_pacbio_second_pilon_indelsonly.diamond.hits"
```

#### Illumina

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
sbatch --nodes=1 --ntasks-per-node=8 --mem=70000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools map2cov -i illumina_megahit_final_contigs.fa -b publicdb/PRJEB10338/PRJEB10338.sorted.merged.bam -b publicdb/PRJEB21624/PRJEB21624.sorted.merged.bam -b publicdb/PRJEB8939/PRJEB8939.sorted.merged.bam -b publicdb/PRJNA214227/PRJNA214227.sorted.merged.bam -b publicdb/PRJNA255688/PRJNA255688.sorted.merged.bam -b publicdb/PRJNA270714/PRJNA270714.sorted.merged.bam -b publicdb/PRJNA280381/PRJNA280381.sorted.merged.bam -b publicdb/PRJNA291523/PRJNA291523.sorted.merged.bam -b publicdb/PRJNA366460/PRJNA366460.sorted.merged.bam -b publicdb/PRJNA366463/PRJNA366463.sorted.merged.bam -b publicdb/PRJNA366471/PRJNA366471.sorted.merged.bam -b publicdb/PRJNA366487/PRJNA366487.sorted.merged.bam -b publicdb/PRJNA366591/PRJNA366591.sorted.merged.bam -b publicdb/PRJNA366667/PRJNA366667.sorted.merged.bam -b publicdb/PRJNA366681/PRJNA366681.sorted.merged.bam -b publicdb/PRJNA398239/PRJNA398239.sorted.merged.bam -b publicdb/PRJNA60251/PRJNA60251.sorted.merged.bam -b publicdb/USDA/USDA.sorted.merged.bam -o usda_illum_megahit"

scp -pr illumina_accum_diamond_uniprot.ilmn_accum_diamond_uniprot.tsv.taxified.out derek.bickharhth@ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/ilmn_accum_diamond_uniprot.tsv.taxified.out

for i in *.cov; do echo -n "-c $i "; done; echo
sbatch --nodes=1 --ntasks-per-node=2 --mem=25000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools create -i illumina_megahit_final_contigs.perl.fa -t ilmn_accum_diamond_uniprot.tsv.taxified.out --db /home/derek.bickharhth/rumen_longread_metagenome_assembly/binaries/blobtools/data/nodesDB.txt -o illumina_megahit_blobplot -c usda_illum_megahit.PRJEB10338.sorted.merged.bam.cov -c usda_illum_megahit.PRJEB21624.sorted.merged.bam.cov -c usda_illum_megahit.PRJEB8939.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA214227.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA255688.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA270714.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA280381.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA291523.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366460.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366463.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366471.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366487.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366591.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366667.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA366681.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA398239.sorted.merged.bam.cov -c usda_illum_megahit.PRJNA60251.sorted.merged.bam.cov -c usda_illum_megahit.USDA.sorted.merged.bam.cov"

# Generating the table
sbatch --nodes=1 --ntasks-per-node=4 --mem=50000 -p short --wrap="~/rumen_longread_metagenome_assembly/binaries/blobtools/blobtools view -i illumina_megahit_blobplot.blobDB.json --out illumina_blobplot_all -r all -b"

# Generating the plots
for i in phylum superkingdom genus; do echo $i; sbatch ../run_blobplot_generation.sh $i illumina illumina_megahit_blobplot.blobDB.json; done

for i in *.png; do echo $i; perl -e 'chomp(@ARGV); @asegs = split(/\./, $ARGV[0]); shift(@asegs);$temp = join(".", @asegs); $temp =~ s/\.blobDB\.json\.bestsum//g; $temp =~ s/\.p14\.span\.100//g; print "$temp\n"; system("mv $ARGV[0] $temp");' $i; done

```

## Mash contig assignment

> Assembler2: 

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash screen -p 10 -i 0.85 -w  ../hungate.msh usda_pacbio_second_pilon_indelsonly.fa > pacbio_second_pilon_hungate_mashscreen.tab
/mnt/nfs/nfs2/bickhart-users/binaries/mash-Linux64-v2.0/mash screen -p 10 -i 0.85 -w  ../refseq.nucl_plas.k21s1000.msh usda_pacbio_second_pilon_indelsonly.fa > pacbio_second_pilon_refseq_mashscreen.tab
```

## Rarefaction using Nonpareil and kmerspectrum analyzer

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina

```bash
# First, testing nonpareil on one of the reads from the YMPrep_run3 Illumina data
sbatch --nodes=1 --mem=150000 --ntasks-per-node=60 -p assemble3 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/nonpareil/nonpareil -s YMPrepCannula_run3_L3_R1.fastq.gz -T kmer -b YMPrepCannula_nonpareil_test_L3_R1 -t 20 -R 150000"

module load jellyfish/2.2.3
cp /mnt/nfs/nfs2/bickhart-users/binaries/kmerspectrumanalyzer/scripts/countkmers21.sh ./
# I hard-linked the kmer-tool2 invocation in the script
vim countkmers21.sh

sbatch --nodes=1 --mem=60000 --ntasks-per-node=4 -p assemble1 --wrap="zcat YMPrepCannula_run3_L*.fastq.gz | ./countkmers21.sh"
# Bah, it's really outdated and unsuitable for analysis. Let me generate the jellyfish db the old-fashioned way

ls *.fastq.gz | xargs -n 1 echo gunzip -c > generators
sbatch --nodes=1 --mem=200000 --ntasks-per-node=20 -p assemble1 --wrap="jellyfish count -m 21 -s 100M -t 20 -C -g generators -o illumina_run3_21mer"

jellyfish histo -o illumina_run3_21mer.histo -t 10 illumina_run3_21mer

# Damn, the jellyfish2 output didn't work! Let's try jellyfish1
sbatch --nodes=1 --mem=200000 --ntasks-per-node=20 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/jellyfish-1.1.11/bin/jellyfish count -m 21 -s 20000M -t 20 -C -o illumina_run3_21mer.j1 YMPrepCannula_run3_L3_R1.fastq"

# Oh, pssh! I should have read the script more carefully. It interprets kmer histograms instead of the jellyfish dump. Let's do a jellyfish screen with a huge spectrum plot and include the pacbio reads too
sbatch --nodes=1 --mem=15000 --ntasks-per-node=2 -p assemble3 --wrap="jellyfish histo -o illumina_run3_21mer.high.histo -t 10 --high=10000000 illumina_run3_21mer"

echo gunzip -c ../error_corrected_reads/rumen_pacbio_corrected.fasta.gz > pbgenerators
sbatch --nodes=1 --mem=300000 --ntasks-per-node=20 -p assemble1 --wrap="jellyfish count -m 21 -s 100M -t 20 -C -o pacbio_error_corrected_21mer -g pbgenerators"

# Plotting illumina jf reads
python /mnt/nfs/nfs2/bickhart-users/binaries/kmerspectrumanalyzer/scripts/plotkmerspectrum.py illumina_run3_21mer.high.histo -w png -g 1
python /mnt/nfs/nfs2/bickhart-users/binaries/kmerspectrumanalyzer/scripts/plotkmerspectrum.py illumina_run3_21mer.high.histo -w png -g 5
python /mnt/nfs/nfs2/bickhart-users/binaries/kmerspectrumanalyzer/scripts/plotkmerspectrum.py illumina_run3_21mer.high.histo -w png -g 6

# Now for the Pacbio jf reads
sbatch --nodes=1 --mem=15000 --ntasks-per-node=2 -p assemble3 --wrap="jellyfish histo -o pacbio_error_corrected_21mer.high.histo -t 10 --high=10000000 pacbio_error_corrected_21mer"

ls illumina_run3_21mer.high.histo pacbio_error_corrected_21mer.high.histo > histos_list.tab
for i in `seq 1 6`; do python /mnt/nfs/nfs2/bickhart-users/binaries/kmerspectrumanalyzer/scripts/plotkmerspectrum.py -l histos_list.tab -w png -g $i; done

```


Attempting to use another tool for rarefaction. This one is bbmap and it has a "saturation" curve generation program. It needs 100 bytes per read pair, and I think that we have close to 500 million read pairs to analyze (ie. 50 gbytes). Let's confirm first before I dump a big memory job on the cluster.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/pilot_project/illumina

```bash
module load java_sdk/64/1.8.0_121
# testing read pair counts
sbatch --nodes=1 --mem=5000 --ntasks-per-node=1 -p short --wrap="gunzip -c YMPrepCannula_S1_L001_R1_001.fastq.gz | wc -l"

# I think that it is actually going to be unnecessary because the script only processes a pair at a time!
# Using lane 3 because that had the fewest issues on the Nextseq
sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --mail-type=ALL --mail-user=derek.bickhart@ars.usda.gov ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_S1_L003_R1_001.fastq.gz in2=YMPrepCannula_S1_L003_R2_001.fastq.gz out=YMPrepCannula_S1_L003.bbunique.stats interval=100000 cumulative=t count=t printlastbin=t minprob=10 -Xmx20g

# Hmm, that only gave me quality estimates! Trying something else
sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_S1_L003_R1_001.fastq.gz in2=YMPrepCannula_S1_L003_R2_001.fastq.gz out=YMPrepCannula_S1_L003.bbunique.nominprob.stats interval=100000 cumulative=t count=t printlastbin=t -Xmx20g

# That worked, but let's try to generate the histogram without the cumulative function
sbatch --nodes=1 --ntasks-per-node=4 --mem=31000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_S1_L003_R1_001.fastq.gz in2=YMPrepCannula_S1_L003_R2_001.fastq.gz out=YMPrepCannula_S1_L003.bbunique.nocumulative.stats interval=100000 count=t printlastbin=t -Xmx30g

# Now to combine all the data together for one larger estimate
sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p short --wrap='for i in YMPrepCannula_S1_L001_R1*; do gunzip -c $i; done | gzip > YMPrepCannula_combined_R1.fastq.gz'
sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p short --wrap='for i in YMPrepCannula_S1_L001_R2*; do gunzip -c $i; done | gzip > YMPrepCannula_combined_R2.fastq.gz'

# In the meantime, let's estimate the kmer rarity of our error corrected pacbio reads
sbatch --nodes=1 --ntasks-per-node=4 --mem=31000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=../pacbio/rumen_pacbio_corrected.fasta.gz out=rumen_pacbio_corrected.bbunique.stats count=t printlastbin=t -Xmx30g

sbatch --nodes=1 --ntasks-per-node=4 --mem=61000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_combined_R1.fastq.gz in2=YMPrepCannula_combined_R2.fastq.gz out=YMPrepCannula_combined.bbunique.stats count=t printlastbin=t -Xmx60g

### The GC problem ###
# I think that low GC protists are screwing up the rarefaction estimates from BB unique. Let's separate out the components and see if that changes rarefaction estimates
sbatch extractLowGCReads.pl -f YMPrepCannula_combined_R1.fastq.gz -r YMPrepCannula_combined_R2.fastq.gz -c 0.3 -s YMPrepCannula_combined_30_gccutoff.stats
cat YMPrepCannula_combined_30_gccutoff.stats
Total pairs: 142952897
GCCutoff: 0.3
LowGC: 51622141
HighGC: 90795373
LowQualFilt: 535383

for i in YMPrepCannula_combined_R1.f.hgc.fq YMPrepCannula_combined_R1.f.lgc.fq YMPrepCannula_combined_R2.f.hgc.fq YMPrepCannula_combined_R2.f.lgc.fq; do echo $i; sbatch --nodes=1 --ntasks-per-node=4 --mem=5000 -p short --wrap="pigz $i"; done

sbatch --nodes=1 --ntasks-per-node=4 --mem=61000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_combined_R1.f.hgc.fq.gz in2=YMPrepCannula_combined_R2.f.hgc.fq.gz out=YMPrepCannula_combined.hgc.bbunique.stats count=t printlastbin=t -Xmx60g
sbatch --nodes=1 --ntasks-per-node=4 --mem=61000 -p short --mail-type=ALL --mail-user=derek.bickharhth ~/rumen_longread_metagenome_assembly/binaries/bbmap/bbcountunique.sh in=YMPrepCannula_combined_R1.f.lgc.fq.gz in2=YMPrepCannula_combined_R2.f.lgc.fq.gz out=YMPrepCannula_combined.lgc.bbunique.stats count=t printlastbin=t -Xmx60g
```

## DAS_tool concatenation

#### Illumina megahit

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit

```bash
# OK, I need to grep out the BIN ids for each file from the metabat and Hi-C data
perl -e '@f = `ls public_metabat/*.fa`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-2]\n";}} close IN;}' > illumina_megahit_public_metabat.unsorted.bins

mkdir hicbins
mv illumina_megahit_hicbins.zip ./hicbins/

perl -e '@f = `ls hicbins/best_genome_clusters/*.fasta`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-2]\n";}} close IN;}' > illumina_megahit_hic.unsorted.bins

wc -l *.bins
 183402 illumina_megahit_hic.unsorted.bins
 387167 illumina_megahit_public_metabat.unsorted.bins

# I think that the hi-c bins are missing quite a few fasta entries! Google must be screwing up  here. I'll check later
module load dastool/1.1.0 diamond prodigalorffinder/2.6.3
sbatch --nodes=1 --ntasks-per-node=10 --mem=45000 -p short --wrap="DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c illumina_megahit_final_contigs.perl.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 10"

# OK, the Iowa team screwed up the DAS_tool installation and did not extract the db zip into the proper directory. Doing that now
cd ..
cp /software/7/apps/dastool/1.1.0/db.zip ./
unzip db.zip
mkdir db
mv *.faa ./db
mv *.lookup ./db
cd illumina_megahit/


sbatch --nodes=1 --ntasks-per-node=10 --mem=45000 -p short --wrap="DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c illumina_megahit_final_contigs.perl.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 10 --db_directory /scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly/assemblies/pilot_project/db"

```

OK, they didn't even install ruby!! I have to transition this to Steve's server.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/illumina_usda_accumulated

```bash
scp -pr derek.bickharhth@ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/illumina_megahit/*.bins ./

export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin/:$PATH
sbatch --nodes=1 --ntasks-per-node=20 --mem=100000 -p assemble1 --wrap="/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c mick_megahit_final_full.rfmt.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 20 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db"

# Dammit! It crashed because it's missing the "doMC" package! It's an R package that needed to be installed
# Ah, I see, Rscript and R interpreter use different installation pathways for local libs
/usr/bin/Rscript -e "install.packages('doMC', repos='http://cran.us.r-project.org')"
/usr/bin/Rscript -e "install.packages('./package/DASTool_1.1.0.tar.gz', type='source')"
/usr/bin/Rscript -e "install.packages('ggplot2', repos='http://cran.us.r-project.org')"

# Ah, I know why! It's because Steve doesn't show the home directory to the nodes, and all of my local R libs are being packaged in my home directory
/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c mick_megahit_final_full.rfmt.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 20 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db

# Preparing prodigal output for shipment to Bradd
gunzip -c illumina_megahit_prodigal_proteins.faa.gz | grep '>' | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > illumina_megahit_prodigal_proteins.shortform.tab

/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i illumina_megahit_hic.unsorted.bins,illumina_megahit_public_metabat.unsorted.bins -c mick_megahit_final_full.rfmt.fa -o illumina_megahit_dastool -l HiC,metabat --search_engine diamond -t 30 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db --write_bins 1

perl -e '@f = `ls illumina_megahit_dastool_DASTool_bins/*.fa`;chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-3]\n";}} close IN;}'> illumina_megahit_dastool_DASTool_bins.tab
```

#### Pacbio pilon

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon

```bash
# Gathering the bins
perl -e '@f = `ls metabat2/*.fa`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-2]\n";}} close IN;}'> pacbio_final_public_metabat.unsorted.bins

mkdir hic_clusters
mv new_rumen_pacbio_best_genome_clusters.zip ./hic_clusters/
perl -e '@f = `ls hic_clusters/best_genome_clusters/*.fasta`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-2]\n";}} close IN;}' > pacbio_final_public_hic.unsorted.bins
```

Now to do the DAS_tool concatenation

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_final_pilon

```bash
scp -pr derek.bickharhth@ceres:/home/derek.bickharhth/rumen_longread_metagenome_assembly/assemblies/pilot_project/pacbio_final_pilon/*.bins ./

export PATH=/mnt/nfs/nfs2/bickhart-users/binaries/bin/:$PATH
/mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/DAS_Tool -i pacbio_final_public_hic.unsorted.bins,pacbio_final_public_metabat.unsorted.bins -c usda_pacbio_second_pilon_indelsonly.fa -o pacbio_final_dastool -l HiC,metabat --search_engine diamond -t 30 --db_directory /mnt/nfs/nfs2/bickhart-users/binaries/DAS_Tool/db --write_bins 1

# Generating list of master bins
perl -e '@f = `ls pacbio_final_dastool_DASTool_bins/*.fa`; chomp(@f); foreach $h (@f){@hsegs = split(/\./, $h); open(IN, "< $h"); while(<IN>){if($_ =~ /^>/){chomp; $_ =~ s/>//; print "$_\t$hsegs[-3]\n";}} close IN;}'> pacbio_final_dastool_DASTool_bins.tab

# Transforming prodigal predictions into the shortform tab delimited file
grep  '>' pacbio_final_dastool_proteins.faa | grep '>' | perl -e '@ids = ("ID", "partial", "start_type", "rbs_motif", "rbs_spacer", "gc_cont"); print "ContigID\tStart\tEnd\tOrient\t" . join("\t", @ids) . "\n"; while(<STDIN>){chomp; $_ =~ s/\>//g; @lsegs = split(/\s*\#\s*/); @bsegs = split(/\;/, $lsegs[-1]); print "$lsegs[0]\t$lsegs[1]\t$lsegs[2]\t$lsegs[3]"; foreach my $k (@bsegs){$k =~ s/^.+\=//; print "\t$k";} print "\n";}' > pacbio_final_prodigal_proteins.shortform.tab


```

## Tallying all data into a larger table

We have lots of bin information and metrics on each assembly. I want to do cross tool comparisons, but it's difficult with all of the data separated in different compartments. Let's gather it all together into one larger table.

Here's the information that I'm dealing with:
* Bin assignment (DAStool concatenation?)
* Read depth
* GC percentage
* Diamond tax classification
* Keggs/cogs present
* Present in Hungate and/or refseq above X cutoff?

It turns out that Blobtools did this for me! Now to try to sort out the data into something more graphical (ie. an NMDS).

> Assembler2: /mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/blob_data_tables

```bash
grep -v '##' illumina_blobplot_all.illumina_megahit_blobplot.blobDB.table.txt > illumina_blobplot_all.reformat.tab
grep -v '##' pacbio_secpilon_blobplot_all.pacbio_secpilon_blobplot.blobDB.table.txt > pacbio_secpilon_blobplot_all.reformat.tab
```

```R
library(dplyr)
library(vegan)

illumina <- read.delim("illumina_blobplot_all.reformat.tab", header=TRUE)
illumina.reformat <- illumina[,c(2:24)]
rownames(illumina.reformat) <- illumina[,1]
illumina.reformat$Tech <- c("Illumina")

pacbio <- read.delim("pacbio_secpilon_blobplot_all.reformat.tab", header=TRUE)
pacbio.reformat <- pacbio[,c(2:24)]
rownames(pacbio.reformat) <- pacbio[,1]
pacbio.reformat$Tech <- c("Pacbio")

combined <- bind_rows(illumina.reformat, pacbio.reformat)

# OK, this probably won't end well! Let's explode things anyways
nmds <- metaMDS(combined[,c(1:23)], k=2, trymax=1000)

# OK, it can only work on numeric values. Let's also remove the coverage sum as that's a derived variable
nmds <- metaMDS(combined[,c(1:21)], k=2, trymax=1000)
Error: cannot allocate vector of size 19026.2 Gb

# Hah! OK, let's try it with fewer values
nmds <- metaMDS(combined[,c(1:3, 22)], k=2, trymax=1000)
Error: cannot allocate vector of size 19026.2 Gb

# OK, it's because of the rows, not the columns. I need to remove some of the "chaff" from the Illumina assembly
combined.filtered <- filter(combined, length > 2572)
nmds <- metaMDS(combined.filtered[,c(1:3, 22)], k=2, trymax=1000)
Error: cannot allocate vector of size 1188.5 Gb

# Getting closer! Let's do it by 5000
combined.filtered <- filter(combined, length > 5000)
nmds <- metaMDS(combined.filtered[,c(1:3, 22)], k=2, trymax=1000)
Error in distfun(comm, method = distance, ...) :
  long vectors (argument 4) are not supported in .Fortran

# damn! Let's try a principal component analysis instead
pca <- rda(combined.filtered[,c(1:3, 22)])
png("rumen_combined_pca_plot.first.png", width=2000, height=1500)
biplot(pca, display = c("sites", "species"), type= c("text", "points"))
def.off()

# OK, not the best! Let's try again
pca <- rda(combined.filtered[,c(1:21)])
biplot(pca, display = c("sites", "species"), type= c("text", "points"), cex=3)
ordihull(pca, group=combined.filtered$Tech)
dev.off()

# Still not good. Perhaps I need to include something like # of orfs and the like?
```

#### Hypothesis: GC bias is not present in the reads

I'm testing a hypothesis that the GC bias identified in the pacbio dataset is not present in the original reads.

> Ceres: /home/derek.bickharhth/rumen_longread_metagenome_assembly/sequence_data/pilot_project/illumina

```bash
for i in YMPrepCannula_S1*.fastq.gz; do echo $i; sbatch -p short ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/metagenomics_scripts/calculateReadGCHisto.pl -f $i -o $i.gcbin; done
sbatch -p short ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/metagenomics_scripts/calculateReadGCHisto.pl -f ../pacbio/rumen_pacbio_corrected.fasta.gz -o rumen_pacbio_corrected.fasta.gz.gcbin

# now to combine all of the Illumina bins
perl -e '@files = `ls YMPrep*.gcbin`; chomp(@files); %gcbin; foreach my $f (@files){open(IN, "< $f"); <IN>; while(<IN>){chomp; @s = split(/\t/); $gcbin{$s[0]} += $s[1];} close IN;} open(OUT, "> YMPrepCombined_allfastas.gcbin"); print OUT "GC\tCount\n"; foreach my $gc (sort{$a <=> $b}keys(%gcbin)){print OUT "$gc\t$gcbin{$gc}\n";} close OUT;'

# The Nextseq kinda hosed us here. The GC percentages of individual reads is quite high
module load r
srun --nodes=1 --mem=8000 --ntasks-per-node=2 --pty bash

# Rerunning with revised program and output
for i in YMPrepCannula_S1*.fastq.gz; do echo $i; sbatch -p short ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/metagenomics_scripts/calculateReadGCHisto.pl -f $i -o $i.gclist; done
sbatch -p short ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/metagenomics_scripts/calculateReadGCHisto.pl -f ../pacbio/rumen_pacbio_corrected.fasta.gz -o rumen_pacbio_corrected.fasta.gz.gclist

# Cleaning up the illumina reads to decrease the verbosity
sbatch average_gc_lists.pl YMPrepCannula_S1_L001_R1_001.fastq.gz.gclist YMPrepCannula_S1_L001_R2_001.fastq.gz.gclist YMPrepCannula_S1_L001_avg.gclist
sbatch average_gc_lists.pl YMPrepCannula_S1_L002_R1_001.fastq.gz.gclist YMPrepCannula_S1_L002_R2_001.fastq.gz.gclist YMPrepCannula_S1_L002_avg.gclist
sbatch average_gc_lists.pl YMPrepCannula_S1_L003_R1_001.fastq.gz.gclist YMPrepCannula_S1_L003_R2_001.fastq.gz.gclist YMPrepCannula_S1_L003_avg.gclist
sbatch average_gc_lists.pl YMPrepCannula_S1_L004_R1_001.fastq.gz.gclist YMPrepCannula_S1_L004_R2_001.fastq.gz.gclist YMPrepCannula_S1_L004_avg.gclist

sbatch average_gc_lists.pl YMPrepCannula_S1_L001_avg.gclist YMPrepCannula_S1_L002_avg.gclist YMPrepCannula_S1_L003_avg.gclist YMPrepCannula_S1_L004_avg.gclist YMPrepCannula_S1_combined_avg.gclist
srun --nodes=1 --mem=50000 --ntasks-per-node=3 -p short --pty bash
```

```R
illumina <- read.delim("YMPrepCombined_allfastas.gcbin", header=TRUE)
pacbio <- read.delim("rumen_pacbio_corrected.fasta.gz.gcbin", header=TRUE)

# Removing the "All N's" and "All G's" artifacts from the count
illumina.filt <- illumina[seq(from=3, to=nrow(illumina)-2),]

colnames(illumina.filt) <- c("GC", "Illumina")
colnames(pacbio) <- c("GC", "PacBio")

combined <- left_join(illumina.filt, pacbio, by="GC")
combined$PacBio[is.na(combined$PacBio)] <- 0
combined.format <- combined[,c(2,3)]
rownames(combined.format) <- combined[,1]

library(plotly)
pdf("rumen_read_gc_perc.pdf", useDingbats=FALSE)
p <- ggplot(combined, aes(x=GC)) + geom_point(aes(y=Illumina, color="red")) + geom_point(aes(y=PacBio, color="blue")) + stat_smooth(aes(y=Illumi> p color="darkred")) + stat_smooth(aes(y=PacBio, color="darkblue"))
dev.off()

# That did not go the way I had hoped! I need to create a file with GC counts per read. 
# Trying with the flat list
files <- c("YMPrepCannula_S1_L001_avg.gclist", "YMPrepCannula_S1_L002_avg.gclist", "YMPrepCannula_S1_L003_avg.gclist", "YMPrepCannula_S1_L004_avg.gclist")

library(data.table)
#illumina <- do.call(rbind, lapply(files, fread))
illumina <- read.delim("YMPrepCannula_S1_combined_avg.gclist", header=FALSE)
pacbio <- read.delim("rumen_pacbio_corrected.fasta.gz.gclist", header=FALSE)
colnames(illumina) <- c("GC")
colnames(pacbio) <- c("GC")

# Adding categories
illumina <- mutate(illumina, Tech=c("Illumina"))
pacbio <- mutate(pacbio, Tech=c("PacBio"))

combined <- bind_rows(illumina, pacbio)
pdf("rumen_read_gc_violin_plot.pdf", useDingbats=FALSE)
p <- ggplot(combined, aes(x=Tech, y=GC, fill=Tech)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white")
p + labs(x="Technology", y="Avg GC ratio") + scale_fill_brewer(palette="Dark2") + theme_classic()

# I think that I may have "smoothed out" too much in the combined averaging
illumina <- read.delim("YMPrepCannula_S1_L003_avg.gclist", header=FALSE)

```

##### Conclusions: The Bias is present in both datasets, but the low GC data isn't being assembled by either dataset. The Protists are too complex and have reference genomes that are too large. 

##### Idea: Remove the low GC reads before doing kmer rarefaction, as we are interested in teasing out the coverage of those datasets as well.