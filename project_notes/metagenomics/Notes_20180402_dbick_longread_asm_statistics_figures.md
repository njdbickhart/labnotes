# USDA Pacbio assembly statistics comparison
---
*4/2/2018*

These are my notes on generating figures, tables and other statistics for our assembly manuscript.

## Table of Contents


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

```bash
pwd
	/mnt/nfs/nfs2/bickhart-users/metagenomics_projects/pilot_project/pacbio_pilon_accumulated
perl -lane 'if($F[0] =~ /\#/){next;} else{print "$F[0]\t$F[1]\t$F[5]";}' < usda_pacbio_supkingdom.usda_pacbio_pilon_blobplot.blobDB.table.txt > usda_pacbio_supkingdom.contig.lens.tab
```

```R
library(ggplot2)
data <- read.delim("usda_pacbio_supkingdom.contig.lens.tab", header=FALSE)

colnames(data) <- c("Contig", "Len", "Kingdom")
pdf(file="usda_pacbio_supkingdom_contiglens.pdf", useDingbats=FALSE)
ggplot(data, aes(x = Len, fill = Kingdom)) + geom_density(alpha = 0.3) + scale_x_log10(breaks = c(1000, 10000, 25000, 100000, 250000))
dev.off()
```