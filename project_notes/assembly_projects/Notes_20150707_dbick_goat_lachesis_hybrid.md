# Lachesis reverse mapping for hybrid assembly
---

*7/7/2015*

These are my notes on mapping the pacbio contigs back to the Lachesis scaffolds in order to generate comparisons for the final hybrid assembly.

## Table of Contents

* [Preparation](#preparation)
* [Order File Alignment](#order)
* [Summary Table of Results](#summary)
* [Oriented RH map assignments](#orient)
* [Gap information scan](#gap)

<a name="preparation"></a>
## Preparation

I'm going to do this locally because I don't want to wait on the slow file transfer to Blade14. 

> pwd: /home/dbickhart/share/goat_assembly_paper/super_scaffolding/lachesis

```bash
# It's inefficient, but I'm going to use BWA mem for quick contig alignments
bwa index goat_sds_31scaffolds.fasta

# BWA mem has an upper limit of 1 Mbp for alignments. I'm going to have to split some of the larger PacBio contigs to make this work
# I wrote a script to split each PacBio contig into chunks under 9900000
perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaByLength.pl -f goat_split_36ctg_assembly.fa -o lachesis/goat_sub_1mb_assembly.fa

# Just checking my work
samtools faidx lachesis/goat_sub_1mb_assembly.fa
wc -l lachesis/goat_sub_1mb_assembly.fa.fai
	5080 lachesis/goat_sub_1mb_assembly.fa.fai	# far fewer than I was expecting!
perl -lane 'if($F[1] > 990000){print $_;}' < lachesis/goat_sub_1mb_assembly.fa.fai
# Returned nothing
```

OK, things are prepared for alignment. I'm going to try to leave this running all night and hope it doesn't suck up all of the memory on my machine.

```bash
bwa mem goat_sds_31scaffolds.fasta goat_sub_1mb_assembly.fa > goat_sub_1mb_assembly_lachesis_mappings.sam
```

*7/8/2015*

--

OK, this took longer than I expected. I'm going to have to move this to the server because of HD IO issues. Let's first ensure that I'm not redoing the same reads:

> pwd: /home/dbickhart/share/goat_assembly_paper/super_scaffolding/lachesis

```bash
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] < 2000){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < goat_sub_1mb_assembly_lachesis_mappings.sam | cut -f1 | sort | uniq > finished_bwa_mem_first_run.txt

# Now to remove the segments that finished last night
perl ../../../programs_source/Perl/perl_toolchain/sequence_data_scripts/filterFastaByName.pl -f goat_sub_1mb_assembly.fa -l finished_bwa_mem_first_run.txt -o goat_sub_1mb_assembly_restart2.fa
```

I then gzipped both files and shipped them to Blade14 via the nfs mount.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/lachesis/

```bash
# I couldn't transfer the bwa indicies, so I will regenerate them
bwa index goat_sds_31scaffolds.fasta
bwa mem -t 10 goat_sds_31scaffolds.fasta goat_sub_1mb_assembly_restart2.fa > goat_sub_1mb_assembly_lachesis_restart_mappings.sam
```

*7/10/2015*

--

It turns out that I just needed to ask Shawn for the ordering list, rather than brute-forcing contig assignments! I'm going to try to do direct comparisons between the RH mappings and the ordering files.

<a name="order"></a>
## Order file comparisons
Let's first pull the exact contig names from Brian's RH map, then I'm going to do a modified Smith-Waterman alignment in order to find the optimal order. 

> pwd: /home/dbickhart/share/goat_assembly_paper/super_scaffolding/lachesis

```bash
perl -lane 'if($F[0] eq "RH" || !($F[5] =~ /^utg.+/)){next;}else{$F[5] =~ /(utg\d+)\_.+/; print "$F[1]\t$1";}' < ../RHmap_to_PacBio-v3_contigs.txt | uniq > RH_map_pacbio_contig_order.tab

# To do alignments, I really need to segregate these by chromosome
for i in `seq 1 29`; do echo $i; grep "^$i\s" RH_map_pacbio_contig_order.tab | perl -lane 'print $F[1];' > chr${i}.rh_map.order; done

grep "^X" RH_map_pacbio_contig_order.tab | perl -lane 'print $F[1];' > chrX.rh_map.order

mkdir rh_map_order
mv chr*.order ./rh_map_order/
```

Now, I want to make a modification of the Smith-Waterman alignment to generate contig comparisons. 

*7/11/2015*

--

I made the modification program, but I must be very careful about the alignment gap scoring method. I did not assign a gap opening penalty in order to be very conscious of the affline penalty (as numerous single contig gaps are likely between the methods). I may need to revisit this after consideration.

I may also need to do better space formatting on the output. 

Let's start with the input formatting on Shawn's order files.

> pwd: /home/dbickhart/share/goat_assembly_paper/super_scaffolding/lachesis

```bash
for i in Order_files-2015-07-10/Order_files/*.ordering; do echo $i; outname=`basename $i | cut -d'.' -f1`; perl -lane 'if($F[0] =~ /^#/){next;}else{$F[1] =~ s/\.\d+//g; print $F[1];}' < $i > Order_files-2015-07-10/${outname}.list; done

# Just making sure that the line numbers add up
grep -v "^#" Order_files-2015-07-10/Order_files/group0.ordering | wc -l
	71

wc -l Order_files-2015-07-10/group0.list
	71 Order_files-2015-07-10/group0.list

# Now for the test run of the SW program
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr1.rh_map.order -b Order_files-2015-07-10/group0.list -o test.out

# The first set of conditions (match = 5, gap = -4, mismatch = -1) wasn't permissive enough for gapping the A region. Let's reduce it to a more simplistic setting: (match = 1, gap = -1, mismatch = 0)
# It worked! I also set the max values to the lower diagonal of the matrix to force all values to be accounted in the matrix
```

It's not perfect, and needs manual editing, but it will be far faster than scoring the matches by hand!

The script is located on my perl_toolchain [github repository](https://github.com/njdbickhart/perl_toolchain/blob/master/assembly_scripts/assemblyScaffoldSWAlign.pl).

Let's automate all of this and generate all of the contig alignments.

```bash
mkdir alignments
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr1.rh_map.order -b Order_files-2015-07-10/group0.list -o alignments/order0_chr1.align
...

# Oops! Ran into the first problem! The group orders do not frequently match
# Let's try to get my venn program in action
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group1.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr13.rh_map.order 36
	rh_map_order/chr19.rh_map.order 26
	rh_map_order/chr22.rh_map.order 7
	rh_map_order/chr2.rh_map.order 8

# OK, so group1 is split between more than 4 chromosomes on the RH map
# Let's see if we can get alignments for them
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr13.rh_map.order -b Order_files-2015-07-10/group1.list -o alignments/order1_chr13.align

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr19.rh_map.order -b Order_files-2015-07-10/group1.list -o alignments/order1_chr19.align

# This one was a clear forward alignment:
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr22.rh_map.order -b Order_files-2015-07-10/group1.list -o alignments/order1_chr22.align

# Some contig order issues:
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr2.rh_map.order -b Order_files-2015-07-10/group1.list -o alignments/order1_chr2.align

# It is actually quite fragmented! There are clear inversions and segregated fragments here! Perhaps BioNano superscaffolds could resolve this?
```

Using my [venn comparison script](https://github.com/njdbickhart/perl_toolchain/blob/master/bed_cnv_fig_table_pipeline/nameListVennCount.pl) was very handy here.

My thoughts: the venn components are clearly important for identifying which chromosome "cluster" the Lachesis data should be compared against. Also, the inter-chromosome associations are clearly biasing the clustering data. This is clearly important biological data, but it interferes with the scaffolding. 

There may be a way to associate superscaffolds in a node-based way to resolve these issues. 

Let's finalize the data and send it out to the group.

```bash
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group2.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr2.rh_map.order 43
	rh_map_order/chr3.rh_map.order 2
	rh_map_order/chr5.rh_map.order 2

# We'll just do chr2 alignment
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr2.rh_map.order -b Order_files-2015-07-10/group2.list -o alignments/order2_chr2.align -d
	# The end of the RH chr2 is missing in this cluster

# group3
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group3.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr3.rh_map.order 63
	# and a few "1's" on other chrs

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr3.rh_map.order -b Order_files-2015-07-10/group3.list -o alignments/order3_chr3.align
	# There appears to be an inversion in the middle of the chromosome

# group4
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group4.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr4.rh_map.order 32
	# and a few "1's"

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr4.rh_map.order -b Order_files-2015-07-10/group4.list -o alignments/order4_chr4.align

# group5
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group5.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr5.rh_map.order 49
	rh_map_order/chr9.rh_map.order 3

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr5.rh_map.order -b Order_files-2015-07-10/group5.list -o alignments/order5_chr5.align
	# Inversion at the end of the chromosome

# group6
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group6.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr6.rh_map.order 51
	# from now on, "1's" will be noted if absent in the association

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr6.rh_map.order -b Order_files-2015-07-10/group6.list -o alignments/order6_chr6.align
	# Another inversion, but most of it is here

# group7
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group7.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr8.rh_map.order 56

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr8.rh_map.order -b Order_files-2015-07-10/group7.list -o alignments/order7_chr8.align
	# Another inversion

# group8
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group8.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr7.rh_map.order 43

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr7.rh_map.order -b Order_files-2015-07-10/group8.list -o alignments/order8_chr7.align
	# Most intact assembly (apart from chr1) so far

# group9
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group9.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr11.rh_map.order 52

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr11.rh_map.order -b Order_files-2015-07-10/group9.list -o alignments/order9_chr11.align
	# Inversion

# group10
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group10.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr10.rh_map.order 45

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr10.rh_map.order -b Order_files-2015-07-10/group10.list -o alignments/order10_chr10.align
	# Very good forward alignment

# group11
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group11.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr14.rh_map.order 46

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr14.rh_map.order -b Order_files-2015-07-10/group11.list -o alignments/order11_chr14.align
	# Inversion

# group12
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group12.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr9.rh_map.order 32

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr9.rh_map.order -b Order_files-2015-07-10/group12.list -o alignments/order12_chr9.align
	# Inversion at the beginning

# group13
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group13.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr12.rh_map.order 21

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr12.rh_map.order -b Order_files-2015-07-10/group13.list -o alignments/order13_chr12.align
	# Lots of genomic segment islands that are not represented in the alignment

# group14
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group14.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr15.rh_map.order 43

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr15.rh_map.order -b Order_files-2015-07-10/group14.list -o alignments/order14_chr15.align

# group15
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group15.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr16.rh_map.order 45

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr16.rh_map.order -b Order_files-2015-07-10/group15.list -o alignments/order15_chr16.align
	# inversion 

# group16
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group16.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chrX.rh_map.order 222
	# Very few "1's"

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chrX.rh_map.order -b Order_files-2015-07-10/group16.list -o alignments/order16_chrX.align
	# The alignments are very small -- lots of misarranged regions
	# Likely a problem with both the RH map and the Lachesis? maybe Lachesis is better here?

# group17
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group17.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr17.rh_map.order 45
	# Very few "1's"

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr17.rh_map.order -b Order_files-2015-07-10/group17.list -o alignments/order17_chr17.align
	# Inversion

# group18
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group18.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr20.rh_map.order 17

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr20.rh_map.order -b Order_files-2015-07-10/group18.list -o alignments/order18_chr20.align
	# Lots of rearranged chromosome segments

# group 19
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group19.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr21.rh_map.order 43

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr21.rh_map.order -b Order_files-2015-07-10/group19.list -o alignments/order19_chr21.align
	# Inversion

# group20
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group20.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr18.rh_map.order 46

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr18.rh_map.order -b Order_files-2015-07-10/group20.list -o alignments/order20_chr18.align
	# Lots of rearranged segments

# group21
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group21.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr24.rh_map.order 23

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr24.rh_map.order -b Order_files-2015-07-10/group21.list -o alignments/order21_chr24.align
	# Inversion

# group22
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group22.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr23.rh_map.order 26

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr23.rh_map.order -b Order_files-2015-07-10/group22.list -o alignments/order22_chr23.align
	# Inversion

# group23
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group23.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr29.rh_map.order

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr29.rh_map.order -b Order_files-2015-07-10/group23.list -o alignments/order23_chr29.align
	# Inversion

# group24
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group24.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr26.rh_map.order 21

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr26.rh_map.order -b Order_files-2015-07-10/group24.list -o alignments/order24_chr26.align
	# Inversion

# group25
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group25.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chrX.rh_map.order 196  <- maybe this is part of the X and Y??

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chrX.rh_map.order -b Order_files-2015-07-10/group25.list -o alignments/order25_chrX.align
	# Again, more disorder -- perhaps the PAR has confused the clustering?

# group26
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group26.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr22.rh_map.order 20

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr22.rh_map.order -b Order_files-2015-07-10/group26.list -o alignments/order26_chr22.align
	# slight gap towards the end

# group27
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group27.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr28.rh_map.order

# Note, I fixed a bug with my assembly scaffold alignment script here
perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr28.rh_map.order -b Order_files-2015-07-10/group27.list -o alignments/order27_chr28.align
	# Inversion

# group28
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group28.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr27.rh_map.order 23

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr27.rh_map.order -b Order_files-2015-07-10/group28.list -o alignments/order28_chr27.align -d
	# Inversion

# group29
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group29.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr25.rh_map.order 31

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr25.rh_map.order -b Order_files-2015-07-10/group29.list -o alignments/order29_chr25.align
	#Inversion

# group30
for i in rh_map_order/chr*.order; do echo -n "$i "; value=`perl ../../../programs_source/Perl/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl $i Order_files-2015-07-10/group30.list | grep "1;2" | cut -f2`; echo $value; done
	rh_map_order/chr19.rh_map.order 2
	rh_map_order/chr2.rh_map.order 1
	# Novel scaffolds?!

perl ../../../programs_source/Perl/perl_toolchain/assembly_scripts/assemblyScaffoldSWAlign.pl -a rh_map_order/chr19.rh_map.order -b Order_files-2015-07-10/group30.list -o alignments/order30_chr19.align
	# Only two contigs
```
<a name="summary"></a>
## Alignment Summary Table

Now to summarize all of that in a table:

| Cluster # | RHmap chrs | Counts | Comments |
| :--- | :--- | :--- | :--- |
group0 | Chr1 | Chr1:63 | Most map to Chr1 -- some single contig mappings to other chrs |
group1 | Chr13,Chr19,Chr22,Chr2 | Chr13:36; Chr19:26; Chr22:7, Chr2:8| The most fragmented cluster
group2 |  Chr2,Chr3,Chr5 | Chr2:43; Chr3:2, Chr5:2 | a bit more reasonable 
group3 | Chr3 | Chr3:63 | |
group4 | Chr4 | Chr4:32 | |
group5 | Chr5,Chr9 | Chr5:43; Chr9:3 | Clear split + Inversion
group6 | Chr6 | Chr6:51 | Inversion
group7 | Chr8 | Chr8:56 | Inversion
group8 | Chr7 | Chr7:43 | Intact -- no inversion
group9 | Chr11 | Chr11:52 | Inversion
group10 | Chr10 | Chr10:45 | Intact -- no inversion
group11 | Chr14 | Chr14:46 | Inversion
group12 | Chr9 | Chr9:42 | Inversion at the beginning
group13 | Chr12 | Chr12:21 | Most contigs were not in contiguous alignment
group14 | Chr15 | Chr15:43 | 
group15 | Chr16 | Chr16:45 | Inversion
group16 | ChrX | ChrX:222 | Lots of small sub-alignments -- very discontinuous
group17 | Chr17 | Chr17:45 | Inversion
group18 | Chr20 | Chr20:17 | Most contigs were not in contiguous alignment
group19 | Chr21 | Chr21:43 | Inversion
group20 | Chr18 | Chr18:46 | Most contigs were not in contiguous alignment
group21 | Chr24 | Chr24:23 | Inversion
group22 | Chr23 | Chr23:26 | Inversion
group23 | Chr29 | Chr29:30 | Inversion
group24 | Chr26 | Chr26:21 | Inversion
group25 | ChrX | ChrX:196 | Lots of small sub-alignments -- very discontinuous
group26 | Chr22 | Chr22:20 | Gap of contigs at the end
group27 | Chr28 | Chr28:11 | Inversion
group28 | Chr27 | Chr27:23 | Inversion
group29 | Chr25 | Chr25:31 | Inversion
group30 | Chr19,Chr2 | Chr19:2; Chr2:1 | Hard to tell!


In light of this information, it's clear that group1 is a chimeric scaffold derived from inter-chromosome interactions from several autosomes. Also, the two groups that have the X chromosome contigs COULD be the two sex chromosomes, but it is difficult to scaffold them from there! Another possibility is that the scaffolding "forked" once it reached the PAR and that one cluster is part of the X and the other is a chimeric Y + X? 

Ivan and Shawn will know best how to reanalyze the data.

<a name="orient"></a>
## Oriented RH map assignments

OK, I have simple contig order lists, but I don't have the contigs assigned in terms of their orientation on the RH map. I will write a quick one-shot perl script to do this.

Here's the script (it's so singular in purpose that it's not worth archiving in my github repository):

```perl
#!/usr/bin/perl
# This is a one shot script designed to quickly interrogate and order the pacbio contigs on the RH map
# Output format:
#	chr	contigname	start	end	orientation["+", "-", "?"]

use strict;

my $output = "ordered_rh_map_contigs.tab";
my $input = "RHmap_to_PacBio-v3_contigs.txt";

open(my $IN, "< $input") || die "Could not find input file!\n";
open(my $OUT, "> $output");
my $header = <$IN>;
my @current;
my %data; # {chr}->[]->[contigname, min, max, orient]
my $chrcontig = "0"; my $curchr = "0";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[5] =~ /^\*/){next;}
	if($segs[5] =~ /\#N\/A/){next;}
	
	# get the simple contig name
	my @ctgsegs = split(/\_/, $segs[5]);
	$segs[5] = $ctgsegs[0];
	
	if($segs[5] ne $chrcontig && $chrcontig ne "0" && scalar(@current) > 0){
		# new contig, and we have data to sort through
		push(@{$data{$curchr}}, determineOrder(\@current));
		@current = ();
	}elsif($segs[5] ne $chrcontig && scalar(@current) == 0){
		# new contig and there's only one entry for it!
		push(@{$data{$curchr}}, [$chrcontig, $segs[6], $segs[6], "?"]);
	}
	
	push(@current, \@segs);
	$chrcontig = $segs[5];
	$curchr = $segs[1];
}
if(scalar(@current) >= 1){
	push(@{$data{$curchr}}, determineOrder(\@current));
}

foreach my $chr (sort {$a <=> $b} keys(%data)){
	foreach my $row (@{$data{$chr}}){
		print {$OUT} "$chr\t" . $row->[0] . "\t" . $row->[1] . "\t" . $row->[2] . "\t" . $row->[3] . "\n";
	}
}

close $IN;
close $OUT;
	
exit;

sub determineOrder{
	my ($array) = @_;
	my $min = 20000000;
	my $max = 0;
	my $orient = ($array->[0]->[6] - $array->[1]->[6] < 0)? "+" : "-";
	# I'm using empty string concatenation and redundant math to ensure that the Perl interpreter
	# Doesn't just store the address of the array in memory
	my $ctgname = $array->[0]->[5] . "";
	foreach my $v (@{$array}){
		if($v->[6] < $min){
			$min = $v->[6] + 1 - 1;
		}
		if($v->[6] > $max){
			$max = $v->[6] + 1 - 1;
		}
	}
	return [$ctgname, $min, $max, $orient];
}
```

I ran this on my local virtualbox installation:

> pwd: /home/dbickhart/share/goat_assembly_paper/super_scaffolding

```bash
perl generate_oriented_rhmap_contigs.pl

# Gzipping original RHmap file just in case
gzip RHmap_to_PacBio-v3_contigs.txt
```

The output file has the following columns:

1. Chromosome
2. PacBio contig abbreviation
3. minimum alignment position (start coordinate)
4. maximum alignment position (end coordinate)
5. Orientation (can be "+", "-" or "?" if there is only one point of reference)

<a name="gap"></a>
## Gap information scan

*8/24/2015*

I need to pull information from the assemblies so that gap regions and lengths can be easily determined. I'm going to write a Java program to pull this from the fasta quickly.

--

*8/25/2015*

OK, I wrote the program [GetMaskBedFasta](https://github.com/njdbickhart/GetMaskBedFasta) and ran it from my IDE in console mode. It produced 10,712 gap regions, and 10,520 of which were greater than 100bp in length.

Example settings for running the program:

> pwd: /home/dbickhart/share/goat_assembly_paper/bng_scaffolds

```bash
java -jar GetMaskBedFasta.jar -f Goat333HybScaffolds1242contigs0723.fasta -o bng_gaps.bed -s bng_stats.tab
```

I'm going to try splitting the BNG superscaffolds to help Shawn with the Lachesis alignment in the meantime. I'll select only the gap regions that are larger than 100 bp (larger than the read length?) and I've updated my fasta splitting script, [splitFastaWBreakpointBed.pl](https://github.com/njdbickhart/perl_toolchain/blob/master/sequence_data_scripts/splitFastaWBreakpointBed.pl), so that it does not pull gap sequence from the beginning and ends of scaffold sequence.

> pwd: /home/dbickhart/share/goat_assembly_paper/bng_scaffolds

```bash
# Going to extract only the gap regions above 100 bp
perl -lane 'if($F[2] - $F[1] > 100){print $_;}' < bng_gaps.bed > bng_gaps_above100bp.bed

perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f Goat333HybScaffolds1242contigs0723.fasta -o bng_split_gap.fa -b bng_gaps_above100bp.bed

# Using samtools to get summary statistics from the split fasta
samtools faidx bng_split_gap.fa

wc -l *.fai
  2349 bng_split_gap.fa.fai			<-774 more segments than the original BNG super scaffolds
  1575 Goat333HybScaffolds1242contigs0723.fasta.fai

tail *.fai
==> bng_split_gap.fa.fai <==
...
utg40992        23779   2,657,212,952      23779   23780	<- a loss of 39 megabases

==> Goat333HybScaffolds1242contigs0723.fasta.fai <==
...
utg40992        23779   2,696,524,662      23779   23780
```

OK, there were far fewer split regions than I expected, but this can be attributed to close gaps in the scaffolds that probably corresponded to nickase binding sites? I'll check to confirm:

```bash
# From bng_gaps_above100bp.bed
	Super-Scaffold_1        2673354 2677638
	Super-Scaffold_1        2677646 2682610

# So the extract region would be Super-Scaffold_1:2677638-2677646
samtools faidx Goat333HybScaffolds1242contigs0723.fasta Super-Scaffold_1:2677638-2677646
	>Super-Scaffold_1:2677638-2677646
	NGCTCTTCN

# From Alex's email, the nickase site is: GCTCTTC, so that's a site!
```

I think that I'm ready to gzip this and send it to Shawn.

Wait... getting a bed file with the exact split coordinates would help substantially, I think. Let's modify the perl script and run it again to get the right coordinates.

```bash
perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f Goat333HybScaffolds1242contigs0723.fasta -o bng_split_gap.fa -b bng_gaps_above100bp.bed -s bng_split_gap_coords.bed

# The bng_split_gap_coords.bed file will contain the coordinates that the split fasta entry has on the original scaffold

# Now to set the old fai file in a new, sorted location
mv bng_split_gap.fa.fai bng_split_gap.fa.fai.bak
sort -k1,1 bng_split_gap.fa.fai.bak | cut -f1,2 > bng_split_gap.fa.fai.bak.sorted

samtools faidx bng_split_gap.fa
sort -k1,1 bng_split_gap.fa.fai | cut -f1,2 > bng_retest.fa.fai.sorted

diff bng_retest.fa.fai.sorted bng_split_gap.fa.fai.bak.sorted
# Nothing -- so we're good!
```

--

*8/26/2015*

Just generating some quick summary statistics on the split assembly.

> pwd: /home/dbickhart/share/goat_assembly_paper/bng_scaffolds

```bash
# Calculating contig N50
perl ../../programs_source/Perl/perl_toolchain/assembly_scripts/calculateContigN50.pl bng_split_gap.fa
	N50 length:     1309271228
	N50 value:      12,149,468
	L50 value:      58
```

#### Initial check for 2kb gap split cutoff

Alex mentioned that a 100bp cutoff would be too stringent. Let's try a 2kb gap cutoff and see if that appreciably improves the statistics.

> pwd: /home/dbickhart/share/goat_assembly_paper/bng_scaffolds

```bash
perl -lane 'if($F[2] - $F[1] > 2000){print $_;}' < bng_gaps.bed > bng_gaps_above2kb.bed

wc -l bng_gaps*
  10520 bng_gaps_above100bp.bed
   9329 bng_gaps_above2kb.bed
  10712 bng_gaps.bed

# OK, that was a substantial decrease. Let's make the split fasta now
perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f Goat333HybScaffolds1242contigs0723.fasta -o bng_split_2kb_gap.fa -b bng_gaps_above2kb.bed -s bng_split_2kb_gap_coords.bed

samtools faidx bng_split_2kb_gap.fa

wc -l *.fai
  2583 bng_split_2kb_gap.fa.fai
  2349 bng_split_gap.fa.fai
  1242 bng_unscaffolded_pacbio.fai
  1575 Goat333HybScaffolds1242contigs0723.fasta.fai

# Hmm... that made more segments! How is that possible?
# Let's look at the data
grep -P 'Super-Scaffold_336\.' *.fai | perl -lane '$F[0] =~ s/\:/\|/g; pop(@F); pop(@F); print join("|", @F);'
```

| File | Scaffold | ScaffLen | total Bases |
| :--- | :--- | :--- | ---: |
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.1|**177701**|22
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.2|**155403**|180707
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.3|869399|338723
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.4|1523|1222634
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.5|80668|1224205
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.6|3279|1306240
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.7|72515|1309596
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.8|103326|1383342
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.9|1765|1488413
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.10|1773|1490231
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.11|**326432**|1492057
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.12|1856|1823953
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.13|**105876**|1825863
bng_split_2kb_gap.fa.fai|Super-Scaffold_336.14|198|1933527
bng_split_gap.fa.fai|Super-Scaffold_336.1|**177701**|293574308
bng_split_gap.fa.fai|Super-Scaffold_336.2|**154708**|293754993
bng_split_gap.fa.fai|Super-Scaffold_336.3|96051|293912302
bng_split_gap.fa.fai|Super-Scaffold_336.4|554331|294009976
bng_split_gap.fa.fai|Super-Scaffold_336.5|148739|294573568
bng_split_gap.fa.fai|Super-Scaffold_336.6|69431|294724808
bng_split_gap.fa.fai|Super-Scaffold_336.7|78872|294795419
bng_split_gap.fa.fai|Super-Scaffold_336.8|72515|294875628
bng_split_gap.fa.fai|Super-Scaffold_336.9|102867|294949374
bng_split_gap.fa.fai|Super-Scaffold_336.10|**326432**|295053979
bng_split_gap.fa.fai|Super-Scaffold_336.11|**105876**|295385875

Bold values are the same as in the previous file. It seems obvious that smaller distances between nickase sites allowed my scripts's internal "100bp" length filter for splitting to kick in. Hell, the end of the chromosome filter failed for the 2kb gap file because there was 200bp uninterrupted!

Calculating the N50:

```bash
perl ../../programs_source/Perl/perl_toolchain/assembly_scripts/calculateContigN50.pl bng_split_2kb_gap.fa
N50 length:     1316999836
N50 value:      14,451,513
L50 value:      47
```

#### Implementation of minimum size filter for scaffolds

Hmm, despite the increase in number of fragments, the N50 went up by 2Mb! Let's try improving the Perl script so that the minimum size filter is customizable.

```bash
# 2kb minimum size filter
perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f Goat333HybScaffolds1242contigs0723.fasta -o bng_split_2kb_2kbf_gap.fa -b bng_gaps_above2kb.bed -s bng_split_2kb_2kbf_gap_coords.bed -m 2000

samtools faidx bng_split_2kb_2kbf_gap.fa

wc -l *.fai
  2146 bng_split_2kb_2kbf_gap.fa.fai
  2583 bng_split_2kb_gap.fa.fai
  2349 bng_split_gap.fa.fai
  1242 bng_unscaffolded_pacbio.fai
  1575 Goat333HybScaffolds1242contigs0723.fasta.fai
  9895 total

perl -lane 'print $F[1];' < bng_split_2kb_2kbf_gap.fa.fai | statStd.pl
	total   2146
	Minimum 435
	Maximum 66863457
	Average 1220070.247903
	Median  58891.5
	Standard Deviation      4774317.508145
	Mode(Highest Distributed Value) 21821

# I wanted to see if the smaller size contigs contained a copious amount of N's
samtools faidx Goat333HybScaffolds1242contigs0723.fasta Super-Scaffold_880:11738277-11741706
>Super-Scaffold_880:11738277-11741706
	# all N's apart from a nickase site that breaks up the distance

# N50 size calculation
perl ../../programs_source/Perl/perl_toolchain/assembly_scripts/calculateContigN50.pl bng_split_2kb_2kbf_gap.fa
	N50 length:     1316999836
	N50 value:      14,451,513
	L50 value:      47
```

So, the 2kb split contains far more contigs with only "N's" due to nickase site interruptions. Here is a table with the scaffold N50 values:

| filter criteria | N50 | 
| :--- | ---: |
100bp gap, 100bp min | 12,149,468
2kb gap, 100bp min | 14,451,513
2kb gap, 2kb min | 14,451,513

#### Implementation of 'N' ratio filter for smaller split contigs

Let's see if we can filter the 2kb gap split file by removing segments that have mostly N's. I rewrote the [splitFastaWBreakpointBed.pl](https://github.com/njdbickhart/perl_toolchain/blob/master/sequence_data_scripts/splitFastaWBreakpointBed.pl) script again to accommodate an N base filter.

> pwd: /home/dbickhart/share/goat_assembly_paper/bng_scaffolds

```bash
# 2kb gaps, 1000bp min size filter, 50% N base filter
perl ../../programs_source/Perl/perl_toolchain/sequence_data_scripts/splitFastaWBreakpointBed.pl -f Goat333HybScaffolds1242contigs0723.fasta -o bng_split_2kb_1kbf_50n_gap.fa -b bng_gaps_above2kb.bed -s bng_split_2kb_1kbf_50n_gap_coords.bed -m 1000 -n 0.5

samtools faidx bng_split_2kb_1kbf_50n_gap.fa

wc -l *.fai
  2075 bng_split_2kb_1kbf_50n_gap.fa.fai
  2146 bng_split_2kb_2kbf_gap.fa.fai
  2583 bng_split_2kb_gap.fa.fai
  2349 bng_split_gap.fa.fai
  1242 bng_unscaffolded_pacbio.fai
  1575 Goat333HybScaffolds1242contigs0723.fasta.fai

perl ../../programs_source/Perl/perl_toolchain/assembly_scripts/calculateContigN50.pl bng_split_2kb_1kbf_50n_gap.fa
	N50 length:     1316999836
	N50 value:      14451513
	L50 value:      47

perl -lane 'print $F[1];' < bng_split_2kb_1kbf_50n_gap.fa.fai | statStd.pl
	total   2075
	Minimum 1501
	Maximum 66863457
	Average 1261708.550843
	Median  61520
	Standard Deviation      4849948.272566
	Mode(Highest Distributed Value) 51884

```

OK, so ultimately this looks to be the best split fasta. Let's package it up and send it to Alex and Shawn.