# Goat assembly analysis
---
*4/12/2016*

These are my notes on the conclusion of the analysis of the goat reference genome. I wanted to segregate these from the main file, which was becoming quite lengthy!

## Table of Contents
* [Spearman rank correlation of Lachesis](#correlation)
* [Preparing data for NCBI](#ncbi)
* [Centromere repeat check](#centromere)
* [Gap check update](#gapupdate)

<a name="correlation"></a>
## Spearman rank correlation of Lachesis

I want to see if I can create a correlation statistic to describe the RH map concordance of our Lachesis assemblies. I'm going to try to use the [Spearman rank correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) to accomplish this.

Some notes on how I'm going to do this:
* I'm going to do this on the three Lachesis assemblies
	* Pre-BNG
	* Post-BNG
	* Pre-Polishing
* All comparisons will be against the RH map indicated order of the contigs/scaffolds
* Ranks will be discrete integers, but will continue through the chromosome segments (ie. chr2 start will be chr1 end + 1)

I'll write a script to parse the RH mappings into a tab delimited file with two columns. Missing contigs and scaffolds will be removed as they may cause unintended behavior with the covariance estimate. 

Here's the script:

```perl
#!/usr/bin/perl
# This is a one-shot script designed to translate Shawn's RH map order file into a simple tab delimited file

use strict;
chomp @ARGV;
my $usage = "perl $0 <input rh text file> <output tab file>\n";
unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open input text file! $ARGV[0]\n";
open(my $OUT, "> $ARGV[1]");
my ($begin, $eqcount, $chrenter, $chr);
$begin = 0; $eqcount = 0; $chrenter = 0; $chr = "NA";
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	if($line =~ /^=/){
		$eqcount++;
		if($eqcount == 2){
			$begin = 1;
		}
		next;
	}
	if($begin){
		if($line =~ /^LACHESIS/){
			$chrenter = 1;
			my $head = <$IN>;
			my @hsegs = split(/\s+/, $line);
			$hsegs[1] =~ m/CHR(\d+)/;
			$chr = $1;
			next;
		}elsif($chrenter && $line =~ /^Summary/){
			$chrenter = 0;
			next;
		}elsif($chrenter){
			my @segs = split(/\s+/, $line);
			if($segs[0] =~ /^\*/){
				next; # skip missed contigs
			}
			unshift(@segs, $chr); # Adding lachesis chr to the front of the array
			print $OUT join("\t", @segs) . "\n";
		}
	}
}
close $IN;
close $OUT;
```

Now I'm going to create another one-shot script to condense these tab files into vectors of ranks that are discrete.

And here is that script:

```perl
#!/usr/bin/perl
# This is a one-shot script designed to process the tabular rh order files into discrete rankings
# 6/20/2016: modified to take into account interchromosome assignments in pre-bng data

use strict;
chomp(@ARGV);
my $usage = "perl $0 <input rh tab file> <output vector file>\n";
unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open input text file! $ARGV[0]\n";
open(my $OUT, "> $ARGV[1]");
my ($chr, $maxl, $maxr, $addl, $addr, @data);
$chr = 0; $maxl = 0; $maxr = 0; $addl = 0; $addr = 0;
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	if($segs[0] != $chr){
		my @values = processdata($addl, $addr, @data);
		foreach my $v (@values){
			print {$OUT} join("\t", @{$v}) . "\n";
		}
		@data = ();
		# update rank variables
		$addl += $maxl + 1;
		$addr += $maxr + 1;
		$maxl = 0; $maxr = 0; $chr = $segs[0];		
	}
	if($segs[2] > $maxl){
		$maxl = $segs[2];
	}
	if($segs[4] > $maxr){
		$maxr = $segs[4];
	}
	push(@data, [@segs]);
}
my @values = processdata($addl, $addr, @data);
foreach my $v (@values){
	print {$OUT} join("\t", @{$v}) . "\n";
}
close $IN;
close $OUT;
	
exit;

sub processdata{
	my ($addl, $addr, @data) = @_;
	my @processed;
	# estimate most frequent rhchr
	my %chrs;
	foreach my $d (@data){
		$chrs{$d->[3]} += 1;
	}
	my $mostchr; my $maxval = 0;
	foreach my $k (keys(%chrs)){
		if($k eq '.'){next;}
		if($chrs{$k} > $maxval){
			$maxval = $chrs{$k};
			$mostchr = $k;
		}
	}
	# process data
	foreach my $d (@data){
		# remove bad entries
		if($d->[4] eq '.'){
			next;
		}
		if($d->[3] != $mostchr){
			push(@processed, [$d->[1], $d->[2] + $addl, -1]);
		}else{
			push(@processed, [$d->[1], $d->[2] + $addl, $d->[4] + $addr]);
		}
	}
	return @processed;
}
```

With those scripts, let's test out how well this correlation metric works to see if I need to tweak how I'm processing the data:

> pwd: /home/dbickhart/share/goat_assembly_paper/lachesis/rh_map_correlations

```bash
# pre BNG
perl convert_rh_order_to_tab.pl pre_bng_lachesis_rhorder.txt pre_bng_lachesis_rhorder.tab
perl generate_rh_order_vectors.pl pre_bng_lachesis_rhorder.tab pre_bng_lachesis_rhorder.vecs

# post BNG
perl convert_rh_order_to_tab.pl post_bng_lachesis_rhorder.txt post_bng_lachesis_rhorder.tab
perl generate_rh_order_vectors.pl post_bng_lachesis_rhorder.tab post_bng_lachesis_rhorder.vecs

```

The "vecs" files should contain a simple list that is easy to import into R. Now I just need to run the correlation plots to see how the values stack up. I'm expecting the post-bng results to be far better than the pre-bng results. We'll see!

```r
data <- read.delim("pre_bng_lachesis_rhorder.vecs", header=FALSE)
cor.test(data$V2, data$V3, alternative = "two.sided", method = "spearman", exact = TRUE)

        Spearmans rank correlation rho

	data:  data$V2 and data$V3
	S = 1017700, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	      rho
	0.9950017

post <- read.delim("post_bng_lachesis_rhorder.vecs", header=FALSE)
cor.test(post$V2, post$V3, alternative = "two.sided", method = "spearman", exact = TRUE)

        Spearmans rank correlation rho

	data:  post$V2 and post$V3
	S = 1027600, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	      rho
	0.9812317

# I plotted the data using the base "plot" function, and it looks like the X chromosome is screwing things up for the "post" BNG rho

```



<a name="ncbi"></a>
## Preparing data for NCBI

I kicked the assembly to the genome submission portal, but NCBI found some contamination issues. I'm going to extract a contig and see how it aligns to the reference for one of the contaminants, then I'm going to prepare the assembly for final submission again (without the problem sequence).

Ah, Serge beat me to the punch here. Here's his email:

>Yeah, I’m always suspicious when they turn up contaminant, I presume whenever the cow gets submitted it will be flagged as contaminated. We recently had a bacterial genome kicked back for having too many psedugenes and therefore mis-assembled when all its near neighbors had the same number of pseudogenes (but a relatively high % compared to an average bacteria). Sometimes these checks are overly ambitious.

>I mapped the contaminant contigs to both campy and delftia. The campy hit is definitely real, the >whole scaffold maps over most of its length at 97% identity. The delftia seems real too, the contigs >only cover about 10% of the genome but they are spread equally across the genome which I wouldn’t >expect if it was an issue with the bacterial assembly so likely just really low coverage. I’d suggest >tossing all the contigs the flagged as trim too (except the adaptor one) because i wouldn’t trust a >contig that’s half bacteria and half goat (which is what their result is implying). We also have a >small hit to a tree genome which also makes no sense. 

Let's drop the sequence that we don't need.

> pwd: /home/dbickhart/share/goat_assembly_paper/final

```bash
# Removing everything
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %skip; while(<IN>){chomp; $_ =~ s/\r//g; @s = split(/\t/); $skip{$s[0]} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($skip{$s[0]})){next;} system("samtools faidx papadum-v13.full.fa $s[0] >> papadum-v13.trimmed.fa"); print STDERR "$s[0] .. ";} print STDERR "\n";' ncbi_contaminant_drop.tab papadum-v13.full.fa.fai

# Checking parity
samtools faidx papadum-v13.trimmed.fa

# Looks good. Now to trim that adaptor sequence
rm papadum-v13.trimmed.fa.fai
samtools faidx papadum-v13.full.fa unplaced_10982:1-10843 | perl -e 'open(OUT, ">> papadum-v13.trimmed.fa"); $h = <>; $h =~ m/(>unplaced_10982):/; print OUT "$1\n"; while(<>){print OUT $_;} close OUT;'

# Final parity check
samtools faidx papadum-v13.trimmed.fa
```

I made sure that the last contig was in place.

## NCBI genome format

In order to submit this to NCBI I need to reformat the fasta so that it can be processed by them.

Here are the facets of the genome that need to be fixed so that we can submit it:

1. Please remove any N nucleotides from the beginning or end of the sequence 
2. Any run of exactly 5 N's is an unknown length gap (could be 1 could be 10000).  All other N's are estimated length gaps.  Is this correct? Are the 55 runs of exactly 5 N's completely unknown length gaps?  
3. It looks like you have a single scaffold for chromosomes 1-29 but chromosome X is in two different scaffolds.  Is this correct? 
4. I noticed there isn't a chromosome Y but you indicated this is a male goat. 
 Is this correct?
5. Are any of the runs of N's in the fasta file scaffold-breaking gaps? 

I think that the only things we need to really address is (1.) above.

Here's my perl script to do this:

```perl
#!/usr/bin/perl
# This is a one-shot script that is designed to trim the beginning and ending "N's" of a fasta

use strict;

chomp(@ARGV);
my $usage = "perl $0 <input fasta> <output fasta> <changes bed>\n";

unless(scalar(@ARGV) == 3){
	print $usage;
	exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open fasta file!\n$usage";
open(my $OUT, "> $ARGV[1]");
open(my $BED, "> $ARGV[2]");
my $header = "NA"; my $seq;
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^>/){
		if($header eq "NA"){
			$header = $line;
			next;
		}
		my ($start) = $seq =~ /^(N+)/;
		my ($end) = $seq =~ /(N+)$/;
		my $nseq = substr($seq, length($start), length($seq) - length($end) - length($start)); 
		$nseq =~ s/(.{60})/$1\n/g;
		chomp($nseq);
		print {$OUT} "$header\n$nseq\n";
		chomp($header);
		$header =~ s/>//g;
		
		if(length($start) > 0){
			my $s = length($start);
			print {$BED} "$header\t1\t$s\ttrim\n";
		}
		if(length($end) > 0){
			my $fs = length($seq) - length($end);
			my $fe = length($seq);
			print {$BED} "$header\t$fs\t$fe\ttrim\n";
		}
		$header = $line;
		$seq = "";
	}else{
		chomp $line;
		$seq .= $line;
	}
}
close $IN;
close $OUT;
close $BED;
```

> pwd: /home/dbickhart/share/goat_assembly_paper/final

```bash
perl remove_beginning_and_end_ns.pl papadum-v13.trimmed.fa papadum-v13.nstrim.fa papadum-v13.nstrim.bed
samtools faidx papadum-v13.nstrim.fa

```

## Identifying filled gaps

I have written a script to try to interrogate gap regions and to see if we filled them in our assembly. My strategy is to take flanking sequence from the CHI_1.0 reference assembly, then align it to our assembly and count the non-N bases in between that we have filled.

I will also generate plots related to gap counts, gap size disparity and other metrics that arise in the data. 

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000317765.1_CHIR_1.0/GCF_000317765.1_CHIR_1.0_genomic.fna.gz
gunzip GCF_000317765.1_CHIR_1.0_genomic.fna.gz
samtools faidx GCF_000317765.1_CHIR_1.0_genomic.fna

# Gotta fix the NCBI naming scheme
perl -ne 'if($_ =~ /^>/){chomp; @s = split(/\s+/); print "$s[0]\n";}else{print $_;}' < GCF_000317765.1_CHIR_1.0_genomic.fna > CHIR_1.0_fixed.fa
samtools faidx CHIR_1.0_fixed.fa

perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o CHIR_1.0_fixed.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d papadumv13_gap_fills.tab

# I had to change the trans-chr naming scheme just to ensure that I capture all info
mv papadumv13_gap_fills.tab papadumv13_gap_fills.old.tab
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o CHIR_1.0_fixed.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d papadumv13_gap_fills.tab
```

Some interesting things about the transchr regions: I found that some of the discrepancies in the CHI_1.0 reference may be due to their permissive scaffolding of smaller contigs of repetitive regions. Example:

Type | CHI chr | CHI start | CHI end | gap size | v13 r1 coords | v13 r2 coords|
:--- | :--- | :--- | :--- | :--- | :--- | :--- |
Trans  | NC_022293.1   |  100025065   |    100025197   |    132  |         cluster_1:101396029-101396530 |  cluster_16:59202103-59202526   
Trans  | NC_022293.1   |  100027515   |    100027852   |    337   |       cluster_1:101398497-101398992  | cluster_6:16544269-16544441    
Trans  | NC_022293.1   |  100028024    |   100028916   |    892   |       cluster_6:16544268-16544440  |   cluster_1:101400723-101401154

The RH map reveals this to be a 400kb region between two sheep probes, whereas our assembly is only about 200kb on the same large contig. This suggests that they merged some stuff together that they shouldn't have!

Let's start collecting some stats on the data.

```bash
# Some bad apples
cat /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/papadumv13_gap_fills.old.tab | grep -v 'Trans' | perl -lane 'if($F[-1] > 100000){print ($F[-1] - $F[4]);}' | wc -l
	1081
# 0.5% of gaps correspond to large differences in assembly sequence placement

# count of deviations
cat /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/papadumv13_gap_fills.old.tab | grep -v 'Trans' | perl -lane 'if($F[-1] < 100000){print ($F[-1] - $F[4]);}' | statStd.pl
	total   241830
	Minimum -22675
	Maximum 99295
	Average 138.423165
	Median  16
	Standard Deviation      2444.962691
	Mode(Highest Distributed Value) 11

# On average, we filled 138 more bases for each gap

# Count of different gap types
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/papadumv13_gap_fills.old.tab -c 0 -m
```
#### Gap counts
|Entry  |  Count|
|:------|------:|
|Closed | 242888|
|Open   |     23|
|Trans  |  13853|

So, 5% of the gaps were trans-chr events, and < 0.1% were unable to be filled by our assembly. 95% of the gaps (- the 0.5% that were large differences on the same cluster) were successfully filled.


## Working on CHI_2.0 gaps

I received the sequence from Serge and I want to rerun the analysis on CHI_2.0. Let's prepare the fasta and work on the data.

> Blade14: /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/

```bash
perl -ne 'if($_ =~ /^>/){chomp; @s = split(/\|/); print ">$s[3]\n";}else{print $_;}' < GCA_000317765.2_CHIR_2.0_genomic.fna > CHIR_2.0_fixed.fa
samtools faidx CHIR_2.0_fixed.fa

wc -l CHIR_2.0_fixed.fa.fai
	102896 CHIR_2.0_fixed.fa.fai # CHI 1.0 was 77,021 ... not sure if they decided to just add more??
```

OK, now to work in a directory off of the shared drive.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/chi2

```bash
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d papadumv13_gap_fills.chi2.tab

# Just a quick comparison of gap regions
wc -l temp.gap.bed
	221091 temp.gap.bed # vs 260474 in CHI_1.0

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/table_bed_length_sum.pl temp.gap.bed
FName   IntNum  TotLen  LenAvg  LenStdev        LenMedian       SmallestL       LargestL
temp.gap.bed    221091  84973034        384.335110881945        1409.90000053502        0       0       22812
# CHI_1.0 had ~137 megabases of gaps, this one has 84.97 mbps

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f papadumv13_gap_fills.chi2.tab -c 0 -m
```

#### CHI_2.0 gap closures

|Entry  | Count|
|:------|-----:|
|Closed | 59737|
|Open   |     4|
|Trans  |  9252|

OK, not as impressive as above. Discounting transchr and open segments, we have 86.58% closure rate (n = 68993 for the tab file).

Let's dig further to get more data.

```bash
# Checking average gap size deviations
cat papadumv13_gap_fills.chi2.tab | grep -v 'Trans' | perl -lane 'if($F[-1] < 100000){print ($F[-1] - $F[4]);}' | statStd.pl
total   59380
Minimum -22653
Maximum 99230
Average -202.261485
Median  20
Standard Deviation      3476.456333
Mode(Highest Distributed Value) 10

# So on average, the BGI assembly was predicting far larger gaps and we were filling them with smaller counts of bases

```

Serge brought up the great point that I should try to validate our gap closures with contigs (ie. anchoring and spanning PacBio reads) just to ensure that we have a good set here. I'm going to map the version 3 contigs to our v13 assembly and then see how many gaps are internal to the contig maps.

**Strategy:** pull 500bp from each end of the contigs and map to the v13 assembly to generate a bed coordinate file. Then intersect with the gap bed file.

> blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/remap

```bash
# I'm using papadum v4, which contains my split contigs that had errors
perl -lane 'if($F[0] =~ /^utg/){$e = $F[1] - 500; print "$F[0]\t1\t500\n$F[0]\t$e\t$F[1]";}' < /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v4/papadum-v4s.fa.fai > papadum_v4s.ends.bed

# Generating the fasta file for BWA realignment
perl -lane '$ucsc = "$F[0]:$F[1]-$F[2]"; open(IN, "samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v4/papadum-v4s.fa $ucsc |"); $h = <IN>; print ">$F[0]_$F[1]_$F[2]"; while(<IN>){ print $_;} close IN;' < papadum_v4s.ends.bed > papadum_v4s.ends.fa

# BWA alignment
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz papadum_v4s.ends.fa > papadum_v4s.ends.sam

# Just checking alignment stats and manipulating the files
samtools view -bS papadum_v4s.ends.sam | samtools sort -T papadum -o papadum_v4s.ends.bam -
samtools index papadum_v4s.ends.bam
samtools idxstats papadum_v4s.ends.bam | tail
	*       0       0       5

# 5 didn't align! Let's take a closer look
samtools view papadum_v4s.ends.bam | tail | perl -lane 'print "$F[0]\t$F[2]";'
	utg1292_730853_731353   *
	utg3966_1_500   *
	utg3966_835613_836113   *
	utg7363_218102_218602   *
	utg8805_1_500   *

# I'll fight over them later, let's work on the ones that did
samtools view papadum_v4s.ends.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); if($s[1] & 2048){next;} @b = split(/_/, $s[0]); push(@{$h{$b[0]}}, [$s[2], $s[3]]);} foreach my $k (keys(%h)){if(scalar(@{$h{$k}}) > 2){print STDERR "Error with $k!\n";}foreach my $r (@{$h{$k}}){print $r->[0] . "\t" . $r->[1] . "\t";} print "$k\n";}' > first_tier_mappings.notbed
perl -lane 'if($F[0] ne $F[2]){print $_;}' < first_tier_mappings.notbed | wc -l
	541

# ok, new strategy, that's too many contigs to lose. Let's go internal here
perl -lane 'if($F[0] =~ /^utg/){$e = $F[1] - 500; $e1 = $F[1] - 1500; print "$F[0]\t500\t1500\n$F[0]\t$e1\t$e";}' < /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v4/papadum-v4s.fa.fai > papadum_v4s.ends.int.bed

perl -lane '$ucsc = "$F[0]:$F[1]-$F[2]"; open(IN, "samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v4/papadum-v4s.fa $ucsc |"); $h = <IN>; print ">$F[0]_$F[1]_$F[2]"; while(<IN>){ print $_;} close IN;' < papadum_v4s.ends.int.bed > papadum_v4s.ends.int.fa

bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz papadum_v4s.ends.int.fa > papadum_v4s.ends.int.sam
samtools view -bS papadum_v4s.ends.int.sam | samtools sort -T papapdum -o papadum_v4s.ends.int.bam -

samtools view papadum_v4s.ends.int.bam | perl -e '%h; while(<>){chomp; @s = split(/\t/); if($s[1] & 2048){next;} @b = split(/_/, $s[0]); push(@{$h{$b[0]}}, [$s[2], $s[3]]);} foreach my $k (keys(%h)){if(scalar(@{$h{$k}}) > 2){print STDERR "Error with $k!\n";}foreach my $r (@{$h{$k}}){print $r->[0] . "\t" . $r->[1] . "\t";} print "$k\n";}' > second_tier_mappings.notbed

perl -lane 'if($F[0] ne $F[2]){print $_;}' < second_tier_mappings.notbed | wc -l
670

# That's even worse!
# working with the first tier mappings for now...
perl -lane 'if($F[0] ne $F[2]){next;}else{print "$F[0]\t$F[1]\t$F[3]\t$F[4]";}' < first_tier_mappings.notbed | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > first_tier_mappings.sort.bed

# Let's try the intersection with the cluster aligning gaps
grep 'Closed' ../chi2/papadumv13_gap_fills.chi2.tab | perl -lane 'if($F[8] > 50000){next;} unless($F[5] =~ /^cluster/){next;} ($fs, $fe) = $F[6] =~ /.+:(\d+)-(\d+)/; ($ss, $se) = $F[7] =~ /.+:(\d+)-(\d+)/; push(@n, ($fs, $fe, $ss, $se)); @n = sort{$a <=> $b}(@n); print "$F[5]\t$n[1]\t$n[2]"; @n = ();' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > cluster_gaps_closed.bed
```

OK, these strategies failed, but Serge was able to generate a nucmer plot of his v3 assembly against v13 so that I can create coordinate bed files for contig assignments to the new assembly. 
The **out.1coords** file is the filtered, high confidence alignments here. Let's condense it down to contig/sequence alignments. Like with everything, I'm going to try the easy way first and then try a more nuanced approach if that gives unsatisfactory results.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check

```bash
perl -lane '@b = split(/_/, $F[12]); print "$F[11]\t$F[0]\t$F[1]\t$b[0]";' < out.1coords > goat_v3_vs_v13_aligns_unrefined.bed

# Now to merge them with a good base "slop" to ensure that close ends are joined in the final bed file.
perl condense_chr_cluster_seq.pl goat_v3_vs_v13_aligns_unrefined.bed goat_v3_vs_v13_aligns_merged.bed

# OK, that seemed to work!
```

Now, I will attempt to find gaps that are closed internal to the contigs as evidence that the closure was completed.

```bash
mkdir closed_confirm

# all gap cluster intersections
intersectBed -a goat_v3_vs_v13_aligns_merged.bed -b remap/cluster_gaps_closed.bed -wa -wb | wc -l
	48,403
# all gaps that are not covered by the cluster bed file
intersectBed -a remap/cluster_gaps_closed.bed -b goat_v3_vs_v13_aligns_merged.bed -v | wc -l
	8,984
# gaps that are not fully within contig boundaries
intersectBed -a goat_v3_vs_v13_aligns_merged.bed -b remap/cluster_gaps_closed.bed -wa -wb | perl -lane 'if($F[5] >= $F[1] && $F[6] <= $F[2]){next;}else{print $_;}' | wc -l
	219
# gaps that ARE within contig boundaries
intersectBed -a goat_v3_vs_v13_aligns_merged.bed -b remap/cluster_gaps_closed.bed -wa -wb | perl -lane 'if($F[5] >= $F[1] && $F[6] <= $F[2]){print $_;}' | wc -l
	48,184
```

OK, let's tabulate the stats. Here are the numbers:

#### CHI 2.0 gap fills

Entry | Count
:--- | ---:
Total gaps in chi2	| 68,993
Closed gaps in chi2 | 59,737
Confirmed closed gaps in chi2 | 48,184
Unconfirmed closed gaps in chi2 | 9,203
Ambiguous gap closures in chi2 (trans) | 9,252
Failed gap closures in chi2 | 4

So **86.6%** of gaps were closed, and **69.83%** were confirmed by PacBio consensus contigs. 

## Identifying CHI 1.0 to 2.0 closures

Now, my goal is to take the CHI 1.0 and CHI 2.0 gaps, identify the improvements that CHI 2.0 made and then see how many of those were confirmed in our assembly.

> Blade14:

```bash
# Generating CHI 1.0 bed file
mkdir chi1
cat /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/papadumv13_gap_fills.old.tab | grep 'Closed' | perl -lane 'if($F[-1] < 10000){($s1, $e1) = $F[6] =~ /.+:(\d+)-(\d+)/; ($s2, $e2) = $F[7] =~ /.+:(\d+)-(\d+)/; my @u; push(@u, ($s1, $s2, $e1, $e2)); @u = sort {$a <=> $b} @u; print "$F[5]\t$u[1]\t$u[2]";}' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > chi1/chi1_closed_gaps.bed
wc -l chi1/chi1_closed_gaps.bed
	241,128 chi1/chi1_closed_gaps.bed

# OK, now to intersect with chi2 to see how many gaps were left open
intersectBed -a chi1/chi1_closed_gaps.bed -b remap/cluster_gaps_closed.bed -v | wc -l
	199,199

# Gaps putatively closed by the assembly
intersectBed -a chi1/chi1_closed_gaps.bed -b remap/cluster_gaps_closed.bed -v | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       199199
        Total Length:           88,459,287

# Total CHI1 gaps
cat chi1/chi1_closed_gaps.bed | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       241128
        Total Length:           133,597,731

# exclusive CHI 2.0 gaps (for some reason!)
intersectBed -b chi1/chi1_closed_gaps.bed -a remap/cluster_gaps_closed.bed -v | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/bed_length_sum.pl
        Interval Numbers:       15548
        Total Length:           14,423,862

intersectBed -a chi1/chi1_closed_gaps.bed -b remap/cluster_gaps_closed.bed -v > chi1/chi1_gaps_closed_by_chi2.bed


# Gaps closed within contig boundaries
intersectBed -a goat_v3_vs_v13_aligns_merged.bed -b chi1/chi1_gaps_closed_by_chi2.bed -wa -wb | perl -lane 'if($F[5] >= $F[1] && $F[6] <= $F[2]){print $_;}' | wc -l
	168,826

# Gaps not covered by contig alignments
intersectBed -a chi1/chi1_gaps_closed_by_chi2.bed -b goat_v3_vs_v13_aligns_merged.bed -v | wc -l
	30,279
```

That makes **70%** of gaps unambiguously closed in this version.

I'm not confident that I'm identifying the gaps correctly. Let's try remapping chi1.0 gaps onto chi2.0 using my alignment block strategy.

```bash
bwa index /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o CHIR_1.0_fixed.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d chi1_to_chi2_gap_summary.tab

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi1_to_chi2_gap_summary.tab -c 0 -m
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f papadumv13_gap_fills.tab -c 0 -m

```

#### Raw count improvement of v13 over Chi_2.0

|Chi_2.0 Entry  |  Count|v13 Entry  |  Count| Perc improvement |
|:------|------:|:------|------:| ----: 
|Closed | 230587|Closed | 242888| 5.3%
|Open   |     51|Open   |     23| 121%
|Trans  |  26126|Trans  |  13853| 88.6%

Chi_2.0 closes **89.8%** of gaps, whereas v13 closes **94.6%**.

Now I'm going to "merge" the documents to generate a different type of count.

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); push(@d, $s[0]);} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); push(@n, ($d[$c], $s[0])); print join(";", @n); @n = (); print "\n"; $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -m
```
256764

#### Direct comparison of gap closure status in v13 (left) and chi_2.0 (right)

|Entry         |  Count| Percentage
|:-------------|------:| ----:
|Closed;Closed | 226783| 88.3%
|Closed;Open   |     51| < 0.0001%
|Closed;Trans  |  16054| 6.3 %
|Open;Closed   |     23| < 0.0001%
|Trans;Closed  |   3781| 1.47%
|Trans;Trans   |  10072| 3.92%

Let's do our due dilligence and identify the Open;Closed category gaps.

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); push(@d, $s[0]);} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); push(@n, ($d[$c], $s[0])); $y = join(";", @n); if($y eq "Open;Closed"){open(TEMP, "grep -P \"$s[1]\t$s[2]\" $ARGV[0] |"); $l = <TEMP>; close TEMP; @h = split(/\t/, $l); print join("\t", @s); print "\t$h[6]\t$h[7]\t$h[8]\t$h[9]\t$h[10]\n";} @n = (); $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab

# All of the Open;Closed regions had neither pair map to our assembly
# Let's grab them and see if there are actually reads that align to these regions

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); push(@d, $s[0]);} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); push(@n, ($d[$c], $s[0])); $y = join(";", @n); if($y eq "Open;Closed"){@u = split("[:-]", $s[6]); @r = split("[:-]", $s[7]); push(@m, ($u[1], $u[2], $r[1], $r[2])); @m = sort{$a <=> $b} @m; print "$s[5]\t$m[0]\t$m[3]\n"; @m = ();} @n = (); $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab > chi2_closed_open_regions.bed

perl -lane '$F[1] -= 100; $F[2] += 100; system("samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa $F[0]:$F[1]-$F[2]");' < chi2_closed_open_regions.bed > chi2_closed_open_regions.fa
bwa index chi2_closed_open_regions.fa
# Doing a small test first to see if this works out
bwa mem chi2_closed_open_regions.fa /mnt/nfs/nfs2/GoatData/sequence_data/C0LYPACXX/Sample_Goat250/Goat250_NoIndex_L001_R1_001.fastq.gz /mnt/nfs/nfs2/GoatData/sequence_data/C0LYPACXX/Sample_Goat250/Goat250_NoIndex_L001_R2_001.fastq.gz > chi2_closed_open_regions.1.1.sam
samtools view -bS chi2_closed_open_regions.1.1.sam | samtools sort -T chi2 -o chi2_closed_open_regions.1.1.bam -

samtools index chi2_closed_open_regions.1.1.bam
samtools idxstats chi2_closed_open_regions.1.1.bam | perl -lane 'if($F[2]){print $_;}'
	CM001714_2:110755006-110756284  1279    6       0
	CM001720_2:92504225-92505523    1299    51      51		<	<
	CM001731_2:61530410-61532110    1701    2       0
	CM001731_2:61534256-61535509    1254    11      1		<
	CM001731_2:61537861-61539039    1179    4       0			<
	CM001731_2:62979583-62980535    953     3       1		<
	AJPT02042923.1:239614-241686    2073    9       0
	AJPT02042923.1:238569-240007    1439    1       0
	CM001739_2:59089786-59091045    1260    4       0
	CM001739_2:75243202-75244649    1448    1       1		<
	CM001734_2:44465535-44466979    1445    2       0
	CM001734_2:44466632-44468626    1995    1       0
	CM001734_2:44468602-44492165    23564   86      55		<	<

# Doing another to confirm
bwa mem chi2_closed_open_regions.fa /mnt/nfs/nfs2/GoatData/sequence_data/C0LYPACXX/Sample_Goat250/Goat250_NoIndex_L001_R1_002.fastq.gz /mnt/nfs/nfs2/GoatData/sequence_data/C0LYPACXX/Sample_Goat250/Goat250_NoIndex_L001_R2_002.fastq.gz | samtools view -bS -F 8 -f 4 - | samtools sort -T temp -o chi2_closed_open_regions.1.2.bam -

```

Six of the regions have high degrees of one-end mapped reads, suggesting that the region is difficult to resolve, and is not representative of live animals.

Let's generate a list of gap regions from our Closed;Open and Closed;Trans list in a bed file to send to Ben for gene intersection.

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); push(@d, $s[0]);} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); push(@n, ($d[$c], $s[0])); $y = join(";", @n); if($y eq "Closed;Open" || $y eq "Closed;Trans"){print "$s[1]\t$s[2]\t$s[3]\n";} @n = (); $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab >  chi1_closed_gaps_unresolved_in_chi2.bed

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]}->{$s[1]}->{$s[2]} = 1;} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[1]}->{$s[2]}->{$s[3]})){@u = split("[:-]", $s[6]); @r = split("[:-]", $s[7]); push(@m, ($u[1], $u[2], $r[1], $r[2])); @m = sort{$a <=> $b} @m; print "$s[5]\t$m[0]\t$m[3]\t$s[1]_$s[2]_$s[3]\n"; @m = ();}} close IN;' chi1_closed_gaps_unresolved_in_chi2.bed papadumv13_gap_fills.tab | perl -lane 'if($F[2] - $F[1] < 50000){print $_;}' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > v13_closed_gaps_unresolved_in_chi2.bed

```

Now I need to cherry pick a good example of a gap that we've completely closed, but the CHI_2.0 and CHI_1.0 assemblies failed to close, and a gene intersects with the region of the gap.

```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); $d = $s[10] - $s[8]; push(@d, "$s[0];$s[6];$s[7];$d");} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); if($s[4] < 100){$c++; next;} $d = $s[10] - $s[8]; if($d < 25){$c++; next;} push(@n, ($d[$c], "$s[0];$s[6];$s[7];$d")); print join(";", @n); @n = (); print "\n"; $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab | head

# That gave me a good list to start. Let's see how many of these events there are.

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); $d = $s[10] - $s[8]; push(@d, "$s[0];$s[6];$s[7];$d");} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); if($s[4] < 100){$c++; next;} $d = $s[10] - $s[8]; if($d < 25){$c++; next;} push(@n, ($d[$c], "$s[0];$s[6];$s[7];$d")); print join(";", @n); @n = (); print "\n"; $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab | wc -l
	38194

# That's a good proportion of gaps that still have ambiguous bases in them in CHI 2!
# Let's take the gap regions in this grouping, filter them further and then intersect them with Ben's gene list

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); $d = $s[10] - $s[8]; push(@d, "$s[0];$s[6];$s[7];$d");} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); if($s[4] < 100){$c++; next;} $d = $s[10] - $s[8]; if($d < 25){$c++; next;} push(@n, ($d[$c], "$s[0];$s[6];$s[7];$d")); print join(";", @n); @n = (); print "\n"; $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab | perl -lane '@d = split(";", $F[0]); if($d[-1] > 10000){next;}else{($c1, $s1, $e1) = $d[1] =~ /(.+):(\d+)-(\d+)/; ($c2, $s2, $e2) = $d[2] =~ /(.+):(\d+)-(\d+)/; print "$c1\t$s1\t$e1\t$c2\t$s2\t$e2\t$d[0];$d[4];$d[3];$d[7];$d[5];$d[6]";}' | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > candidate_gap_regions_forfigure.bedpe

# There were 36,954 candidate gaps after removing the large anomalies. 
# Intersecting them with the exon coordinates
bedtools bed12tobed6 -i Papadum_v13_EVM5.bed12 | intersectBed -a candidate_gap_regions_forfigure.bedpe -b stdin -wb > candidate_gap_regions_forfigure.exon_intersect.bedpe
# there were 14,181 candidates for the whole gene intersection.
# there were 2087 candidates for just the exons (command shown above)

# Here's a nice example:
scaffold_4      704058  704559  scaffold_4      704866  705368  Closed;Closed;0;5064;CM001726_2:75227873-75228374;CM001726_2:75233718-75234153  scaffold_4      703407  704807  evm.model.Scaffold_1370.15      0       -

# NOTE: I confirmed via bam alignments that this gap was correctly closed by our assembly.
# 5kb predicted by CHI_2.0. 300 bp closed gap in our assembly.
samtools view /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-HiSeq-Goat400-Freezev13.bam scaffold_4:704058-705368 | less

# extracting coordinates from the fastas (+/- 10kb)
# CM001726_2:75217873-75244153
# scaffold_4:694058-715368
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz scaffold_4:694058-715368 > scaffold_4-694058_715368.fa
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa CM001726_2:75217873-75244153 > CM001726_2-75217873_75244153.fa

# Now to do a quick nucmer plot
sh ../../john_assembled_contigs/run_nucmer_plot_automation_script.sh scaffold_4-694058_715368.fa CM001726_2-75217873_75244153.fa

# I just need to adjust the coordinates to generate the files for Jana
# 6 kb on the end of the gene
grep 'evm.model.Scaffold_1370.15' Papadum_v13_EVM5.bed12 | bedtools bed12tobed6 -i stdin > chosen_gene_exons.bed
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz scaffold_4:696058-721368 > scaffold_4-696058_721368.fa
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa CM001726_2:75217873-75250153 > CM001726_2-75217873_75250153.fa

# The plots don't look bad, but I need a really clear example
# Let's select a wide range of regions and ask Jana to create plots for each
# I'd like to get big gaps that are > 5kb, but let's see what the distribution looks like first
perl -lane '@b = split(";", $F[6]); print "$b[3]";' <candidate_gap_regions_forfigure.exon_intersect.bedpe | statStd.pl
total   2087
Minimum 25
Maximum 9920
Average 1134.033062
Median  541
Standard Deviation      1566.453545
Mode(Highest Distributed Value) 32

# OK, so we're likely to get 300 regions this way based on back of the napkin calculations
# Let's still test it.
perl -lane '@b = split(";", $F[6]); if($b[3] > 5000){print $_;}' <candidate_gap_regions_forfigure.exon_intersect.bedpe | wc -l
95

# Lower than that, but close!
# Let's automate the fasta generation
mkdir dotter_data
# NOTE: added 100kb to each side of the region -- will have to expand it out to the size of the gene later
perl -lane '@b = split(";", $F[6]); if($b[3] > 5000){print $_;}' <candidate_gap_regions_forfigure.exon_intersect.bedpe | perl -lane 'if($F[0] =~ /unplaced/){next;} my @n; push(@n, ($F[1], $F[2], $F[4], $F[5])); @n = sort{$a <=> $b} @n; my $start = $n[0] - 100000; my $end = $n[3] + 100000; my $name = "$F[0]_$start\_$end"; system("samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz $F[0]:$start-$end > dotter_data/$name.v13.fa"); my @b = split(";", $F[6]); my ($chr, $s1, $e1) = $b[4] =~ /(.+):(\d+)-(\d+)/; my ($s2, $e2) = $b[5] =~ /.+:(\d+)-(\d+)/; my @j; push(@j, ($s1, $e1, $s2, $e2)); @j = sort{$a <=> $b} @j; my $cstart = $j[0] - 100000; my $cend = $j[3] + 100000; system("samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa $chr:$cstart-$cend > dotter_data/$name.chi2.fa");'

# Had to remove 3 sets of groups that had abnormal coordinates
# That should be enough to work with now

# Now I am going to take an example and try to fit the whole gene model over it
perl -lane '@b = split(";", $F[6]); if($b[3] > 5000){print $_;}' <candidate_gap_regions_forfigure.exon_intersect.bedpe | grep AJPT02103338.1
grep evm.model.Scaffold_397.85 Papadum_v13_EVM5.bed12 | head -n 1 | bedtools bed12tobed6 -i stdin > evm.scaffold_397.85.gene.bed

# evm.model.Scaffold_397.89 is unknown
# evm.model.Scaffold_397.85 is mucin-5b-like
# evm.model.Scaffold_397.91 is mucin-5ac
# evm.model.Scaffold_397.87 is TOLLIP
# evm.model.Scaffold_397.81 is mucin-5ac

# HOLD the phone! My percent filled column on the tab file wasn't accurately reporting the number of N's correctly! Here it is:
# Gaps filled from CHI1.0 to CHI2.0 without any N bases:
grep 'Closed' chi1_to_chi2_gap_summary.tab | perl -lane 'if($F[10] - $F[8] < 1){print $_;}' | wc -l
160299

# So here is the revised table
```
|Chi_2.0 Entry  |  Count|v13 Entry  |  Count| Perc improvement |
|:------|------:|:------|------:| ----: 
|Closed |160,299|Closed |242,888| 36%
|Open   | 70,339|Open   |     23| 121%
|Trans  | 26,126|Trans  |  13853| 88.6%


```bash
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); @d; while(<IN>){chomp; @s = split(/\t/); push(@d, $s[0]);} close IN; open(IN, "< $ARGV[1]"); $c = 0; while(<IN>){chomp; @s = split(/\t/); my $state = $s[0]; if($s[0] eq "Closed" && $s[10] - $s[8] > 1){$state = "Open";} push(@n, ($d[$c], $state)); print join(";", @n); @n = (); print "\n"; $c++;} close IN;' papadumv13_gap_fills.tab chi1_to_chi2_gap_summary.tab | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -m
```
256764 total gaps in chi1

#### Final stats of V13 (right) to CHI2.0 (left)
|Entry         |  Count|
|:-------------|------:|
|Closed;Closed | 178117|
|Closed;Open   |  48717|
|Closed;Trans  |  16054|
|Open;Closed   |     20|
|Open;Open     |      3|
|Trans;Closed  |   2699|
|Trans;Open    |   1082|
|Trans;Trans   |  10072|

Now to just retrieve some small stats on CHIR1 for a supplementary table.

```bash
# ARS1 closed gaps in gene models
cat chi1/chi1_closed_gaps.bed | sortBedFileSTDIN.pl | intersectBed -a stdin -b Papadum_v13_EVM5.bed12 | wc -l
101610

# ARS1 closed gaps in exons
bedtools bed12tobed6 -i Papadum_v13_EVM5.bed12 | intersectBed -a chi1/chi1_closed_gaps.bed -b stdin | wc -l
21852

# ARS1 closed gaps upstream of genes
perl -lane 'if($F[5] eq "-"){$ns = $F[2]; $ne = $F[2] + 2000;}else{$ns = $F[1] - 2000; $ne = $F[1];} if($ns < 1){$ns = 1;} print "$F[0]\t$ns\t$ne";' < Papadum_v13_EVM5.bed12 | intersectBed -a chi1/chi1_closed_gaps.bed -b stdin | wc -l
13279

```

<a name="centromere"></a>
## Centromere repeat check

Serge identified several low complexity regions that may be centromeric. I'm interested to see if we've found it so I downloaded the consensus sequence from [an article](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-r10) that claims to have identified it (note: it's in the online supplementary notes). I'm going to pull the sequence from clusters 7 and 9 to see if it matches the centromere.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/centromere_check

```bash
# Just because it's easy, let's try bwa
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz centromere.fa > centromere.sam

samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz cluster_7:112500444-112506199
```

There were just some unplaced chromosomes that appear to harbor the centromeric repeat. I'm guessing that most of the heterochromatin is missing from our assembly or is in the degenerate contigs.

Serge has done some extra work on the centromeres and telomeres, just checking the data here.

> pwd:

```bash
perl -e '<>; while(<>){chomp; @s = split(/\t/); $s[0] =~ s/^>//; my $ss = ($s[3] > $s[4])? $s[4] : $s[3]; my $se = ($s[3] > $s[4])? $s[3] : $s[4]; print "$s[0]\t$ss\t$se\n";}' < v13.telomere > v13.telomere.bed

cat v13.telomere.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -d 500 | perl -lane 'if($F[2] - $F[1] > 2000){print $_;}' > v13.filtered.2kb.telomere.bed

cat centromere_aligns.bed | sortBedFileSTDIN.pl | mergeBed -i stdin -d 9000 > centromere_aligns.filtered.bed
```
<a name="nucmerplot"></a>
## Figure 2 nucmer plot generation

In addition to the FRC plot and the dotplot showing assembly discrepancies, I want to create a nucmer image that shows the entirety of our chr 20 scaffold in a nucmer comparison to bgi's CHI_1.0 chr20. I'll generate the nucmer plot myself here. I will juxtapose the BioNano chr20 scaffold image and generate a CHI_1.0 scaffold tiling map as a comparison.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/figure_2

```bash
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz cluster_17 > ars1_cluster17.fa
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa CM001729_2 > chi2_chr20.fa

# I think that CM001729_2 is the chr20 scaffold. We'll see!
sh ../../john_assembled_contigs/run_nucmer_plot_automation_script.sh chi2_chr20.fa ars1_cluster17.fa
sh ../../john_assembled_contigs/run_nucmer_plot_automation_script.sh ars1_cluster17.fa chi2_chr20.fa
```

In order to generate a gap figure, I had to use R.

```R
library(Sushi)
beddata <- read.delim(file="chi2_chr20.gaps.bed", sep="\t", header=FALSE)
colnames(beddata) <- c("chrom", "start", "end")
plotBed(beddata=beddata, chrom="CM001729_2", chromstart=1, chromend=74161552, row="supplied", palettes=list(SushiColors(7)), type="density")
labelgenome("CM001729_2", 1, 74161552, n=4,scale="Mb",edgeblankfraction=0.10)
dev.copy2pdf(file="chi_2_20_gap_density.pdf", useDingbats=FALSE)

```

## Final stats mash-up

These are just some odd stat pulls to try to finalize the supplementary data. 

#### Y scaffold BNG percentage

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/goat_y

```bash
# Alignments of scaffolds
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz unplaced_bng_scafs.fa > unplaced_bng_scafs.sam

# Conversion to tab delimited format
perl -lane 'if($F[0] =~ /^@/){next;}else{if($F[1] > 1000){next;}print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]";}' < unplaced_bng_scafs.sam > alignments.tab

# Getting just the scaffold names
cat alignments.tab | cut -f3 | sort | uniq

# I then marked the scaffold names with a "1" in the tab file
# Total Y scaffold counts from Serge's filtered list
perl -e 'chomp @ARGV; open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; $c = 0; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]})){$c += $h{$s[0]};}} close IN; print "$c\n";' /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz.fai y_scaffolds_ided.tab
10,446,801

# BNG scaffolds only
perl -e 'chomp @ARGV; open(IN, "< $ARGV[0]"); %h; while(<IN>){chomp; @s = split(/\t/); $h{$s[0]} = $s[1];} close IN; $c = 0; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); if(exists($h{$s[0]}) && $s[1]){$c += $h{$s[0]};}} close IN; print "$c\n";' /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz.fai y_scaffolds_ided.tab
5,243,236
```

#### Lachesis RH map correlations

I am going to redo the correlations by setting misplaced contigs (wrong chromosome determined by RH map) to "-1" in order to better represent order.

> pwd: /home/dbickhart/share/goat_assembly_paper/lachesis/rh_map_correlations

```bash
# Using rewritten script
perl generate_rh_order_vectors.pl post_bng_lachesis_rhorder.tab post_bng_lachesis_rhorder.new.vecs
perl generate_rh_order_vectors.pl pre_bng_lachesis_rhorder.tab pre_bng_lachesis_rhorder.new.vecs
```

```R
data <- read.delim("pre_bng_lachesis_rhorder.new.vecs", header=FALSE)
cor.test(data$V2, data$V3, alternative = "two.sided", method = "spearman")

        Spearmans rank correlation rho

	data:  data$V2 and data$V3
	S = 41518000, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	    rho
	0.9101831
	
	Warning message:
	In cor.test.default(data$V2, data$V3, alternative = "two.sided",  :
	  Cannot compute exact p-value with ties
data <- read.delim("post_bng_lachesis_rhorder.new.vecs", header=FALSE)
cor.test(data$V2, data$V3, alternative = "two.sided", method = "spearman")

        Spearmans rank correlation rho

	data:  data$V2 and data$V3
	S = 1174400, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	     rho
	0.978919
	
	Warning message:
	In cor.test.default(data$V2, data$V3, alternative = "two.sided",  :
	  Cannot compute exact p-value with ties
```


#### X chromosome tandem repeat collapse visualization

I'm going to take the region of the X chromosome that is collapsed in CHIR_2.0 and attempt to make a figure out of it! The data is derived from the BioNano repeat arrays. Here are the coordinates that are collapsed:

> cluster_21:51966063-60909199

I am going to gather the following information to try to make this work:

* 
*  information
* CHIR_2.0 region information (through alignments)

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/x_chr_tandem

```bash
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz cluster_21:51966063-60909199 > xchr_tandem_region.fa

# Aligning segments using my 1kb window alignment consensus map script
perl ~/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f xchr_tandem_region.fa -r /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -o xchr_tandem_region.chir2align.tab

# Because of the fragmentary nature of the alignment, I kinda had to do this by eye:
# block 1:	CM001739_2:128721766-128869308	cluster_21:51966063-52109063
# block 2:  CM001739_2:120306877-122248687	cluster_21:52109063-53786063
# block 3:  CM001739_2:119370571-119702052	cluster_21:53786063-53917063
# block 4:  CM001739_2:122361963-125080251	cluster_21:53917063-56690063
# block 5:  CM001739_2:86072248-86369469(rev) cluster_21:56690063-57145063
# block 6:  

# It looks like CHIR_2.0 may have corrected this gap. Let's check CHIR_1.0
perl ~/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f xchr_tandem_region.fa -r ../gap_check/CHIR_1.0_fixed.fa -o xchr_tandem_region.chir1align.tab
# NC_022322.1 is the X chr

# Interestingly, CHIR_2.0 has done a major rewrite here! CHIR_1.0 was a huge problem region
grep NC_022322.1 xchr_tandem_region.chir1align.tab | perl -lane 'my $s = ($F[4] > $F[5])? $F[5]: $F[4]; my $e = ($F[4] > $F[5])? $F[4]: $F[5]; my @g = split(";", $F[-1]); my $c = 0; foreach $i (@g){$c += $i;} my $avg = $c / scalar(@g); if($avg > 55){print "$F[3]\t$s\t$e";}' | sortBedFileSTDIN.pl | mergeBed -d 10000 -i stdin | perl -lane '$l = $F[2] - $F[1]; if($l > 400000){print join("\t", @F) . "\t$l";}'
NC_022322.1     83806841        84238004        431163
NC_022322.1     113525452       113968114       442662
NC_022322.1     114347970       114824125       476155
NC_022322.1     115946048       118132746       2186698
NC_022322.1     118785435       119973419       1187984
NC_022322.1     119986859       120521860       535001
NC_022322.1     120593321       121132114       538793

# These are the largest segments that exist here. Let's piece this together for a circos diagram
```
51966063

ARS1c | start | end | CHI1c | start | end | orient | actual AStart | actual AEnd | 
:--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
cluster_21 | 4740000 | 5037000 | NC_022322.1 | 83806841 | 84098991 | + | 56706063 | 57003063
cluster_21 | 5040000 | 5177000 | NC_022322.1 | 84105588 | 84238004 | - | 57006063 | 57143063
cluster_21 | 192000 | 728000 | NC_022322.1 | 113524958 | 113968114 | + | 52158063 | 52694063
cluster_21 | 927000 | 1409000 | NC_022322.1 | 114351879 | 114824204 | + | 52893063 | 53375063
cluster_21 | 2409000 | 4682000 | NC_022322.1 | 115946048 | 118131598 | + | 54375063 | 56648063
cluster_21 | 5434000 | 7734000 | NC_022322.1 | 118785435 | 121142948 | + | 57400063 | 59700063
cluster_21 | 7812000 | 8943136 | NW_005101056.1 | 42670 | 1162819 | - | 59778063 | 60909199
cluster_21 | 1 | 189000 | NC_022322.1 | 121312110 | 121524693 | -

60556417-60638863 

816952-899398

60035465-60104062

 296000-364597

```bash
~/RepeatMasker/RepeatMasker -pa 10 -qq -species goat -no_is -gff xchr_tandem_region.fa

# I found that there were 900 charlie TE's in the region!
# They also seem to correlate well with the breakpoints of CHIR_1.0
perl -e '$e = 8942890; for($x = 0; $x < $e; $x += 10000){my $te = $x + 10000; print "cluster_21\t$x\t$te\n";}' > cluster_21_region_windows.10kb.bed
grep 'DNA/hAT-Charlie' xchr_tandem_region_rmask_out.bed > xchr_tandem_region_rmask_out.hat_charlie.bed
bedtools coverage -b cluster_21_region_windows.10kb.bed -a xchr_tandem_region_rmask_out.hat_charlie.bed | sortBedFileSTDIN.pl | perl -lane 'if($F[3]){print $_;}' | bedtools merge -i stdin -d 10000 | wc -l
83
bedtools coverage -b cluster_21_region_windows.10kb.bed -a xchr_tandem_region_rmask_out.hat_charlie.bed | sortBedFileSTDIN.pl | perl -lane 'if($F[3]){print $_;}' | bedtools merge -i stdin -d 10000 | perl -lane 'if($F[2] - $F[1] > 10000){print $_;}' | wc -l
59
bedtools coverage -b cluster_21_region_windows.10kb.bed -a xchr_tandem_region_rmask_out.hat_charlie.bed | sortBedFileSTDIN.pl | perl -lane 'if($F[3]){print $_;}' | bedtools merge -i stdin -d 10000 | perl -lane 'if($F[2] - $F[1] > 10000){print $_;}' > cluster_21_region_charlie_tpases_gt10kb.bed

# I was able to find an interesting pattern of orientation for the Charlie-Tpases
# (+, -, -, -, +) 
# This coincides nicely with their hotspot regions and the breakpoints in the assembly
perl get_charlie_pattern.pl xchr_tandem_region_rmask_out.hat_charlie.bed > cluster_21_region_charlie_tpases_20kb_pattern.bed

```

```bash
Super-Scaffold_1637     2002349 2003349 109     997     1000    0.9970000
Super-Scaffold_1637     2003349 2004349 38      831     1000    0.8310000
Super-Scaffold_1637     2004349 2005349 44      939     1000    0.9390000
Super-Scaffold_1637     2005349 2006349 38      885     1000    0.8850000
Super-Scaffold_1637     2006349 2007349 25      524     1000    0.5240000
Super-Scaffold_1637     2007349 2008349 47      788     1000    0.7880000
Super-Scaffold_1637     2008349 2009349 52      830     1000    0.8300000
Super-Scaffold_1637     2009349 2010349 102     984     1000    0.9840000
Super-Scaffold_1637     2010349 2011349 70      993     1000    0.9930000
Super-Scaffold_1637     2011349 2012349 58      772     1000    0.7720000

cluster_20      1984828 1985828 107     937     1000    0.9370000
cluster_20      1985828 1986828 23      507     1000    0.5070000
cluster_20      1986828 1987828 73      982     1000    0.9820000
cluster_20      1987828 1988828 146     1000    1000    1.0000000

```

#### Repeat analysis

OK, so I need to identify repetitive element regions in both regions and then try to correlate gaps with these regions (wilcoxon rank sum test?). Perhaps by counting bovA/bovB repeats.

First, to generate the repeatmasker output.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/repeat_analysis

```bash
~/RepeatMasker/RepeatMasker -pa 10 -no_is -species goat -q /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa

~/RepeatMasker/RepeatMasker -pa 10 -no_is -species goat -q papadum-v13.full.fa

# Now to convert to bed format for intersections
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < papadum-v13.full.fa.out > papadum-v13.full.repeatmask.out.bed
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\n";}' < CHIR_2.0_fixed.fa.out > CHIR_2.0_fixed.repeatmask.out.bed

# The repeat counts are similar
wc -l *.bed
  5385732 CHIR_2.0_fixed.repeatmask.out.bed
  5391929 papadum-v13.full.repeatmask.out.bed

# Let's see if the lengths are similar
cat papadum-v13.full.repeatmask.out.bed | cut -f7 | statStd.pl
	total   5391929
	Minimum -1
	Maximum 44050
	Average 960.958183
	Median  207
	Standard Deviation      1730.468934
	Mode(Highest Distributed Value) 66

cat CHIR_2.0_fixed.repeatmask.out.bed | cut -f7 | statStd.pl
	total   5385732
	Minimum 0
	Maximum 30408
	Average 920.310183
	Median  206
	Standard Deviation      1661.982519
	Mode(Highest Distributed Value) 66

# And the actual mapping lengths
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/table_bed_length_sum.pl -c papadum-v13.full.repeatmask.out.bed CHIR_2.0_fixed.repeatmask.out.bed
```

| FName | IntNum | TotLen | LenAvg | LenStdev | LenMedian | SmallestL | LargestL |
| :---- | -----: | -----: | -----: | -------: | --------: | --------: | -------: |
| papadum-v13 | 5391929 | 1,403,055,172 | 260.2139 | 499.8013 | 149 | 5 | 33064 |
| CHIR_2.0 | 5385732 | 1,177,568,174 | 218.6458 | 290.2337 | 142 | 5 | 30069 |

```bash
perl -lane 'print "$F[1]\t$F[2]\t$F[3]\t$F[0]";' < ../gap_check/papadumv13_gap_fills.tab | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > chir_2_gaps_withars1_stats.bed

intersectBed -a chir_2_gaps_withars1_stats.bed -b CHIR_2.0_fixed.repeatmask.out.bed | uniq | wc -l
4406
# That's far less than I expected! That's a pain!
# Let's instead relate it back to ARS1

perl -lane 'if($F[0] ne "Closed"){next;} my($fc, $fs, $fe) = $F[6] =~ /(.+):(\d+)-(\d+)/; my ($sc, $ss, $se) = $F[7] =~ /(.+):(\d+)-(\d+)/; my @d; push(@d, ($fs, $fe, $ss, $se)); @d = sort {$a <=> $b} @d; if($d[2] - $d[1] > 10000){next;} print "$sc\t$d[1]\t$d[2]\t$F[1]:$F[2]-$F[3]";' < ../gap_check/chi2/papadumv13_gap_fills.chi2.tab | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/sortBedFileSTDIN.pl > chi_2_gaps_in_ars1.bed

intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed | cut -f4 | uniq | wc -l
52,000

# That is a strangely rounded number
# Stayed the same after a sort and then a uniq
# Checking how many gaps had more than one repeat within them
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed | cut -f4 | perl -e '%h; while(<>){chomp; $h{$_} += 1;} $c = 0; $sum = 0; foreach $k (keys(%h)){if($h{$k} > 1){$c++; $sum += $h{$k};}} $avg = $sum / $c; print "$c\t$avg\n";'
20406   4.07179261001666

# So, almost half of the gaps we closed had an average of 4 repeats in them
# Checking the sizes of repeats that spanned chir2 gaps
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | cut -f11 | statStd.pl
total   114683
Minimum 6
Maximum 28916
Average 1430.226886
Median  576
Standard Deviation      1785.936157
Mode(Highest Distributed Value) 3835

# Counting repeat classes
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 8 > repeat_class_gap_tab_count.tab

# BovB repeats made up 43.58% of all identified repeats in closed gap regions!
# Let's see if it's due to clustering
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | grep 'BovB' | cut -f11 | statStd.pl
total   54196
Minimum 10
Maximum 3847
Average 2251.342424
Median  1670
Standard Deviation      1344.969222
Mode(Highest Distributed Value) 3835

# The lengths indicate that these are full lenght elements
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | grep 'BovB' | cut -f4 | sort | uniq | wc -l
38,923

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f papadum-v13.full.repeatmask.out.bed -c 4 > ars1_full_repeat_class_count.tab

# Getting repeat counts from super class estimates
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f CHIR_2.0_fixed.repeatmask.out.bed -c 5 > chi2_repeat_super_class_count.tab
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f papadum-v13.full.repeatmask.out.bed -c 5 > ars1_repeat_super_class_count.tab

# Including repeat lengths
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[6];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' ars1_repeat_super_class_count.tab papadum-v13.full.repeatmask.out.bed > ars1_repeat_super_class_count.extend.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[6];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' chi2_repeat_super_class_count.tab CHIR_2.0_fixed.repeatmask.out.bed > chi2_repeat_super_class_count.extend.tab

# Just checking to see if I can infer more about the repeat sizes from the gap regions from a graph
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | cut -f11 > repeat_lengths.vec
intersectBed -a chi_2_gaps_in_ars1.bed -b papadum-v13.full.repeatmask.out.bed -wb | grep 'BovB' | cut -f11 > repeat_bovb_lengths.vec

# there may be more repeats in the unplaced regions. Checking
grep 'CM0017' CHIR_2.0_fixed.repeatmask.out.bed | wc -l
	5,140,444
	245,288	<- unplaced

grep 'cluster' papadum-v13.full.repeatmask.out.bed | wc -l
	4,980,501
	411,428 <- unplaced
```

#### GGplot of gap repeat length histogram

```R
rlens <- read.delim("repeat_lengths.vec", header=FALSE)
rbovb <- read.delim("repeat_bovb_lengths.vec", header=FALSE)
data <- data.frame(Class=c("Total"), Lengths=rlens)
data <- rbind(data, data.frame(Class=c("BovB"), Lengths=rbovb))
colnames(data) <- c("Class", "Lengths")

ggplot(data, aes(x=Lengths, fill=Class)) + geom_histogram(binwidth=10, alpha=0.5, position="identity") + xlim(0,4000) + theme_set(theme_gray(base_size = 16))
dev.copy2pdf(file="gap_repeat_lengths.pdf", useDingbats=FALSE)


# Making the full figure repeat plot
data <- read.delim("subtable_repeat_counts.formatted", sep="\t", header=TRUE)
data.ggplot <- data.frame(Class=data$X, Count=data$ARS1Count, Assembly=c("ARS1"))
data.ggplot <- rbind(data.ggplot, data.frame(Class=data$X, Count=data$CHIR2Count, Assembly=c("CHIR2.0")))

fancy_scientific <- function(l){
l <- format(l, scientific = TRUE)
l <- gsub("^(.*)e", "'\\1'e", l)
l <- gsub("e", "%*%10^", l)
 parse(text=l)
}

ggplot(data=data.ggplot, aes(x=Class, y=Count, fill=Assembly)) + geom_bar(stat="identity", position=position_dodge(), colour="black") + scale_y_log10(breaks = c(10,100,1000,10000,100000,1000000), labels = fancy_scientific) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Count (log10)") + xlab("Repeat Class (Average Length (bp))")
```

#### Repeat length completion analysis

OK, new plan, let's see what happens when I pull out only repeats that have > 75% match length to the database.

```bash
# Adding the percent divergence and the length match to the query
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[13] =~ s/[()]//g; $s[11] =~ s/[()]//g; $s[12] =~ s/[()]//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[11] + $s[13]; $totsize = $s[12] + $s[13];}else{$unmapped = $s[13] + $s[11]; $totsize = $s[11] + $s[12];} $perc = ($totsize - $unmapped)/$totsize; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\t$s[1]\t$perc\n";}' < papadum-v13.full.fa.out > papadum-v13.full.repeatmask.extend.out.bed

perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[13] =~ s/[()]//g; $s[11] =~ s/[()]//g; $s[12] =~ s/[()]//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[11] + $s[13]; $totsize = $s[12] + $s[13];}else{$unmapped = $s[13] + $s[11]; $totsize = $s[11] + $s[12];} $perc = ($totsize - $unmapped)/$totsize; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\t$s[1]\t$perc\n";}' < CHIR_2.0_fixed.fa.out > CHIR_2.0_fixed.repeatmask.extend.out.bed

cat papadum-v13.full.repeatmask.extend.out.bed | cut -f9 | statStd.pl
	Maximum 0.999959951942331
	Average 0.458806
	Median  0.320754716981132
	Standard Deviation      0.380204

cat CHIR_2.0_fixed.repeatmask.extend.out.bed | cut -f9 | statStd.pl
	Maximum 1.01666666666667
	Average 0.449307
	Median  0.311320754716981
	Standard Deviation      0.380812

cat papadum-v13.full.repeatmask.extend.out.bed | cut -f8 | statStd.pl
	Maximum 59.2
	Average 17.857126

cat CHIR_2.0_fixed.repeatmask.extend.out.bed | cut -f8 | statStd.pl
	Maximum 66.1
	Average 17.744743

# Time to start filtering and counting!
# > 80% similarity and > 75% length (OUR ASSEMBLY IS A 3% IMPROVEMENT)
perl -lane 'if($F[8] > 0.75 && $F[7] < 20){print $_;}' < papadum-v13.full.repeatmask.extend.out.bed | wc -l
	1,305,013
perl -lane 'if($F[8] > 0.75 && $F[7] < 20){print $_;}' < CHIR_2.0_fixed.repeatmask.extend.out.bed | wc -l
	1,279,437

# Let's make a table of this data
perl -e 'chomp(@ARGV); open(ARS, "< $ARGV[0]"); open(CHI, "< $ARGV[1]"); my $lenthresh = 0.99; my $idthresh = 0; print "LengthThresh\tDisimilarityThresh\tARS1Count\tCHI2Count\n"; for(my $x = 0; $x < 15; $x++){ print "$lenthresh\t$idthresh"; my $ars = 0; my $chi = 0; while(<ARS>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$ars++;}} seek(ARS, 0, 1); while(<CHI>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$chi++;}} seek(CHI, 0, 1); print "\t$ars\t$chi\n"; $lenthresh -= 0.05; $idthresh += 5;}' papadum-v13.full.repeatmask.extend.out.bed CHIR_2.0_fixed.repeatmask.extend.out.bed

# Hmm... We're beating them in all these categories! Let's go lower
perl -e 'chomp(@ARGV); open(ARS, "< $ARGV[0]"); open(CHI, "< $ARGV[1]"); my $lenthresh = 0.29; my $idthresh = 70; print "LengthThresh\tDisimilarityThresh\tARS1Count\tCHI2Count\n"; for(my $x = 0; $x < 6; $x++){ print "$lenthresh\t$idthresh"; my $ars = 0; my $chi = 0; while(<ARS>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$ars++;}} seek(ARS, 0, 0); while(<CHI>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$chi++;}} seek(CHI, 0, 0); print "\t$ars\t$chi\n"; $lenthresh -= 0.05; $idthresh += 5;}' papadum-v13.full.repeatmask.extend.out.bed CHIR_2.0_fixed.repeatmask.extend.out.bed
```

LengthThresh  |  DisimilarityThresh  |    ARS1Count   |    CHI2Count
---: | ---: | ---: | ---: 
0.99  |  0  |     236    | 127
0.94  |  5  |     223177 | 216165
0.89  |  10 |     430240 | 422684
0.84   | 15 |     886145  |867140
0.79   | 20 |     1294524 |1268993
0.74   | 25 |     1539233 |1512908
0.69   | 30 |     1758000 |1709417
0.64   | 35 |     1901918 |1850280
0.59  |  40 |     1987695 |1938407
0.54   | 45 |     2138697 |2092243
0.49  |  50 |     2250482 |2208138
0.44  |  55 |     2366002 |2328662
0.39   | 60 |     2466942 |2429621
0.34   | 65 |     2625083 |2563249
0.29   | 70 |     3022302 |2954212
0.24  |  75 |     3241224| 3149584
0.19  |  80 |     3455405| 3353230
0.14  |  85 |     3651978| 3563806
0.09  |  90 |     3970625| 3901507
0.04  |  95 |     4629640| 4572496

```bash
# I'm curious. Let's do BovB at all thresholds
grep BovB papadum-v13.full.repeatmask.extend.out.bed > papadum-v13.full.repeatmask.extend.out.bovb.bed
grep BovB CHIR_2.0_fixed.repeatmask.extend.out.bed > CHIR_2.0_fixed.repeatmask.extend.out.bovb.bed

perl -e 'chomp(@ARGV); open(ARS, "< $ARGV[0]"); open(CHI, "< $ARGV[1]"); my $lenthresh = 0.99; my $idthresh = 0; print "LengthThresh\tDisimilarityThresh\tARS1Count\tCHI2Count\n"; for(my $x = 0; $x < 20; $x++){ print "$lenthresh\t$idthresh"; my $ars = 0; my $chi = 0; while(<ARS>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$ars++;}} seek(ARS, 0, 0); while(<CHI>){chomp; @s = split(/\t/); if($s[8] >= $lenthresh && $s[7] <= $idthresh){$chi++;}} seek(CHI, 0, 0); print "\t$ars\t$chi\n"; $lenthresh -= 0.05; $idthresh += 5;}' papadum-v13.full.repeatmask.extend.out.bovb.bed CHIR_2.0_fixed.repeatmask.extend.out.bovb.bed

# The comparsion points are very good and follow our expected trend!
# OK, we've agreed to a cutoff value: 75% repeat length and 60% ID (40% dissimilarity)
perl -lane 'if($F[7] <= 40 && $F[8] >= 0.75){print $_;}' < papadum-v13.full.repeatmask.extend.out.bed > papadum-v13.full.repeatmask.extend.75thresh.out.bed
perl -lane 'if($F[7] <= 40 && $F[8] >= 0.75){print $_;}' < CHIR_2.0_fixed.repeatmask.extend.out.bed > CHIR_2.0_fixed.repeatmask.extend.75thresh.out.bed

# Now for the pipeline to generate the tables
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f CHIR_2.0_fixed.repeatmask.extend.75thresh.out.bed -c 5 > chi2_repeat_super_class_count.75thresh.tab
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f papadum-v13.full.repeatmask.extend.75thresh.out.bed -c 5 > ars1_repeat_super_class_count.75thresh.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' ars1_repeat_super_class_count.75thresh.tab papadum-v13.full.repeatmask.extend.75thresh.out.bed > ars1_repeat_super_class_count.extend.75thresh.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' chi2_repeat_super_class_count.75thresh.tab CHIR_2.0_fixed.repeatmask.extend.75thresh.out.bed > chi2_repeat_super_class_count.extend.75thresh.tab


#### On further review, my original repeat length estimate was off! 
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' ars1_repeat_super_class_count.tab papadum-v13.full.repeatmask.out.bed > ars1_repeat_super_class_count.extend.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' chi2_repeat_super_class_count.tab CHIR_2.0_fixed.repeatmask.out.bed > chi2_repeat_super_class_count.extend.tab


# formatting data for R
perl -lane 'if($F[0] eq "Type"){print $_;}else{printf("%s (%d bp)\t%d\t%d\t%d\t%d\t%d\n", $F[0], $F[5], $F[1], $F[2], $F[3], $F[4], $F[5])}' < highqual_repeat_summary.tab > highqual_repeat_summary.rformat.tab
```

> Blade14: /mnt/iscsi/vnx_gliu_7/reference/

```bash
# Generating cattle repeatmasker file
~/RepeatMasker/RepeatMasker -q -pa 10 -species cow -no_is umd3_kary_unmask_ngap.fa

~/RepeatMasker/RepeatMasker -q -pa 20 -species sheep -no_is oviAri3.fa
```

#### Checking immune regions to ensure continuity of CHIR_1 and 2

John raised a point about figure 5 that the mapping to CHIR_1.0 may have been a problem. We can wave this as being a comparison of new de novo assemblies, but let's see if there is an actual improvement in CHIR_2.0.

> NKC: NC_022297.1:91233093-91817092
> LRC: cluster_20:64067497-64886393

OK, so the NKC is a single scaffold, so let's align that to CHIR_2.0. The LRC is difficult, so we'll align that section from our assembly to CHIR_2.0.

```bash
samtools faidx ../CHIR_1.0_fixed.fa NC_022297.1:91233093-91817092 > nkc_region.fa
perl ~/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f nkc_region.fa -r /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -o nkc_region_chi2.tab
Longest aligments:      chr     start   end     length
                        CM001717_2      36755150        100123118       63367968
                        CM001739_2      22431548        43798176        21366628
                        CM001723_2      24456073        27127775        2671702
                        CM001714_2      100384654       101008884       624230

# Checking the number of gaps here
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa CM001714_2:100384654-101008884 > nkc_chir2_aligned_seg.fa
~/jdk1.8.0_05/bin/java -jar ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -f nkc_chir2_aligned_seg.fa -o nkc_chir2_aligned_seg.gap.bed -s nkc_chir2_aligned_seg.gap.stats
wc -l nkc_chir2_aligned_seg.gap.bed
	30 nkc_chir2_aligned_seg.gap.bed

# Now for the LRC region
samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz cluster_20:64067497-64886393 > lrc_papadum_region.fa
perl ~/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f lrc_papadum_region.fa -r /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -o lrc_region_chi2.tab
perl ~/perl_toolchain/assembly_scripts/alignUnitigSectionsToRef.pl -f lrc_papadum_region.fa -r /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -o lrc_region_chi2.tab
Longest aligments:      chr     start   end     length
                        CM001710_2      20394930        143358691       122963761
                        CM001739_2      7626536 126954893       119328357
                        CM001712_2      9127115 117202535       108075420
                        CM001716_2      6375053 105453689       99078636
                        CM001718_2      4454532 79322517        74867985
                        CM001727_2      480173  64940967        64460794
                        CM001723_2      6813625 69742374        62928749
                        CM001714_2      58025634        103872347       45846713
                        CM001726_2      15433754        49419276        33985522
                        CM001719_2      43668247        74739486        31071239
                        CM001711_2      91119110        117871455       26752345
                        CM001725_2      14646374        40569406        25923032
                        CM001717_2      3214536 28828586        25614050
                        CM001732_2      16150671        34999784        18849113

# OK! Allot more issues! Let's check it out
# AJPT02103297.1 is in the middle. It's the same situation
# Cluster_20:259000-691000 is missing in CHIR_2.0
```

## Additional repetitive analysis

Serge wants me to check the Sheep and Cattle assembly repetitive contents using the same 75% length threshold. 

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/repeat_analysis/other_ruminants

```bash
perl -e 'while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[9] eq "+")? "+" : "-"; $s[13] =~ s/-//g; $qlen = $s[14] - $s[13]; $s[15] =~ s/-//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[13] + $s[15]; $totsize = $s[14] + $s[15];}else{$unmapped = $s[15] + $s[13]; $totsize = $s[13] + $s[14];} my $mapdisc = $s[2] + $s[3] + $s[4]; $mapdisc /= 1000; if($totsize == 0 || $unmapped == 0 || $totsize - $unmapped == 0){next;} my $superclass = ($s[11] eq $s[12])? $s[12] : "$s[11]/$s[12]"; $perc = ($totsize - $unmapped)/$totsize; print "$s[5]\t$s[6]\t$s[7]\t$orient\t$s[10]\t$superclass\t$qlen\t$mapdisc\t$perc\n";}' < umd3_cattle_rmsk.out > umd3_cattle_rmsk.repmask.bed

perl -e 'while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[9] eq "+")? "+" : "-"; $s[13] =~ s/-//g; $qlen = $s[14] - $s[13]; $s[15] =~ s/-//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[13] + $s[15]; $totsize = $s[14] + $s[15];}else{$unmapped = $s[15] + $s[13]; $totsize = $s[13] + $s[14];} my $mapdisc = $s[2] + $s[3] + $s[4]; $mapdisc /= 1000; if($totsize == 0 || $unmapped == 0 || $totsize - $unmapped == 0){next;} my $superclass = ($s[11] eq $s[12])? $s[12] : "$s[11]/$s[12]"; $perc = ($totsize - $unmapped)/$totsize; print "$s[5]\t$s[6]\t$s[7]\t$orient\t$s[10]\t$superclass\t$qlen\t$mapdisc\t$perc\n";}' < oari3_sheep_rmsk.out > oari3_sheep_rmsk.repmask.bed

perl -lane 'if($F[7] <= 0.40 && $F[8] >= 0.75){print $_;}' < umd3_cattle_rmsk.repmask.bed > umd3_cattle_rmsk.repmask.75thresh.bed
perl -lane 'if($F[7] <= 0.40 && $F[8] >= 0.75){print $_;}' < oari3_sheep_rmsk.repmask.bed > oari3_sheep_rmsk.repmask.75thresh.bed

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f umd3_cattle_rmsk.repmask.75thresh.bed -c 5 > umd3_cattle_rmsk.repeatsuperclass.75thresh.tab
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f oari3_sheep_rmsk.repmask.75thresh.bed -c 5 > oari3_sheep_rmsk.repeatsuperclass.75thresh.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' umd3_cattle_rmsk.repeatsuperclass.75thresh.tab umd3_cattle_rmsk.repmask.75thresh.bed > umd3_cattle_rmsk.repeatsuperclass.extend.75thresh.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' oari3_sheep_rmsk.repeatsuperclass.75thresh.tab oari3_sheep_rmsk.repmask.75thresh.bed > oari3_sheep_rmsk.repeatsuperclass.extend.75thresh.tab


perl -e 'while(<>){chomp; @F = split(/\t/); if                                                          ($F[0] eq "Type"){print $_;}else{printf("%s (%d bp)\t%d\t%d\t%d\t%d\t%d\n", $F[0], $F[5], $F[1], $F[2], $F[3], $F[                                                          4], $F[5])}}' < highqual_repeat_summary.tab > highqual_repeat_summary.rformat.tab
```

Printing comparison plot to file:

```R
data <-read.delim("combined_ruminant_data.tab", sep="\t", header=TRUE)
data.ggplot <- data.frame(Type = data$Class, Count = data$ARS1Count, Assembly = c("ARS1"))
data.ggplot <- rbind(data.ggplot, data.frame(Type = data$Class, Count = data$CHIR2Count, Assembly = c("CHIR2.0")))
data.ggplot <- rbind(data.ggplot, data.frame(Type = data$Class, Count = data$CattleCount, Assembly = c("UMD3.1")))
data.ggplot <- rbind(data.ggplot, data.frame(Type = data$Class, Count = data$SheepCount, Assembly = c("OARI3.0")))

data.ggplot$Type <- factor(data.ggplot$Type, levels = data.ggplot$Type)
library(ggplot2)

ggplot(data=data.ggplot, aes(x=Type, y=Count, fill=Assembly)) + geom_bar(stat="identity", position=position_dodge(), colour="black") + scale_y_log10(breaks = c(10,100,1000,10000,100000,1000000), labels = fancy_scientific) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Count (log10)") + xlab("Repeat Class (Average Length (bp))")

dev.copy2pdf(file="ruminant_comparison.pdf", useDingbats=FALSE)
```

#### Checking cattle and sheep run outside of UCSC

> Blade14: /mnt/iscsi/vnx_gliu_7/reference

```bash
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[13] =~ s/[()]//g; $s[11] =~ s/[()]//g; $s[12] =~ s/[()]//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[11] + $s[13]; $totsize = $s[12] + $s[13];}else{$unmapped = $s[13] + $s[11]; $totsize = $s[11] + $s[12];} $perc = ($totsize - $unmapped)/$totsize; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\t$s[1]\t$perc\n";}' < umd3_kary_unmask_ngap.fa.out > umd3_kary_unmask_ngap.fa.repeat.extend.bed

perl -lane 'if($F[7] <= 40 && $F[8] >= 0.75){print $_;}' < umd3_kary_unmask_ngap.fa.repeat.extend.bed > umd3_kary_unmask_ngap.fa.repeat.extend.75thresh.bed

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f umd3_kary_unmask_ngap.fa.repeat.extend.75thresh.bed -c 5 > ../goat_assembly/repeat_analysis/other_ruminants/umd3_quick_cattle_rmsk.repeatsuperclass.75thresh.tab
perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' ../goat_assembly/repeat_analysis/other_ruminants/umd3_quick_cattle_rmsk.repeatsuperclass.75thresh.tab umd3_kary_unmask_ngap.fa.repeat.extend.75thresh.bed > ../goat_assembly/repeat_analysis/other_ruminants/umd3_quick_cattle_rmsk.repeatsuperclass.extend.75thresh.tab
```

OK, it looks like the data is totally different -- this warrants a new analysis for Sheep too just to be consistent.

```bash
perl -e '<>; <>; <>; while(<>){ $_ =~ s/^\s+//; @s = split(/\s+/); $orient = ($s[8] eq "+")? "+" : "-"; $qlen = $s[12] - $s[11]; $s[13] =~ s/[()]//g; $s[11] =~ s/[()]//g; $s[12] =~ s/[()]//g; my $unmapped; my $totsize; my $perc; if($orient eq "+"){$unmapped = $s[11] + $s[13]; $totsize = $s[12] + $s[13];}else{$unmapped = $s[13] + $s[11]; $totsize = $s[11] + $s[12];} $perc = ($totsize - $unmapped)/$totsize; print "$s[4]\t$s[5]\t$s[6]\t$orient\t$s[9]\t$s[10]\t$qlen\t$s[1]\t$perc\n";}' < oviAri3.fa.out > oviAri3_quick_sheep_rmsk.repmask.bed

perl -lane 'if($F[7] <= 40 && $F[8] >= 0.75){print $_;}' < oviAri3_quick_sheep_rmsk.repmask.bed > oviAri3_quick_sheep_rmsk.repmask.75thresh.bed

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f oviAri3_quick_sheep_rmsk.repmask.75thresh.bed -c 5 > oviAri3_quick_sheep_rmsk.repeatsuperclass.75thresh.tab

perl -e 'chomp(@ARGV); open(IN, "< $ARGV[0]"); my %data; <IN>; while(<IN>){chomp; @s = split(/\t/); $data{$s[0]} = [$s[1], 0];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/\t/); $data{$s[5]}->[1] += $s[2] - $s[1];} close IN; print "Entry\tCount\tTotLen\tAvgLen\n"; foreach my $k (sort {$a cmp $b} keys(%data)){$count = $data{$k}->[0]; $len = $data{$k}->[1]; $avg = $len / $count; print "$k\t$count\t$len\t$avg\n";}' oviAri3_quick_sheep_rmsk.repeatsuperclass.75thresh.tab oviAri3_quick_sheep_rmsk.repmask.75thresh.bed > oviAri3_quick_sheep_rmsk.repeatsuperclass.extend.75thresh.tab
```
<a name="gapupdate"></a>
## Gap check update 
I need to run a couple of tests on the filled gap analysis in order to satisfy reviewer concerns. Here are the tests I need to run:

* Check that CHIR_1.0 and CHIR_2.0 don't close any gaps in ARS1
* Check WGS alignments in predicted gap closure regions for hard clipping that indicates incorrect gap closure
* Find out how many gaps are in ARS1 gene models

#### Checking ARS1 gap closures in CHIR_2.0

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check

```bash
# First, let's compare ARS1 gaps to CHIR_2.0
mkdir ars1_gaps
samtools faidx papadum-v13.full.fa
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o papadum-v13.full.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d ars1_gaps/ars1_gaps_on_chi2.tab

grep 'Closed' ars1_gaps/ars1_gaps_on_chi2.tab | wc -l
373 <- not a good sign off the bat, lets look at this closer

# Some of these regions are huge, like on the order of 100kb or larger
# I suspect that these are just misaligned reads instead of being true closures
grep 'Closed' ars1_gaps/ars1_gaps_on_chi2.tab | perl -lane 'if($F[8] < 100000){print $_;}' | wc -l
255

# Time to check them out!
```
> Blade14 : /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/ars1_gaps/

```bash
# Starting a bwa run
/mnt/nfs/nfs2/GoatData/sequence_data/C0LYPACXX/Sample_Goat400/*.fastq.gz | perl -e 'use File::Basename; %h; while(<>){chomp; $fname = basename($_); $file = $_; @bname = split(/_/, $fname); push(@{$h{$bname[2]}->{$bname[4]}}, $file); } foreach $lane (keys(%h)){foreach $num (keys(%{$h{$lane}})){print join("\t", @{$h{$lane}->{$num}}); print "\tgoat400\tpapadum\n";}}' > papadum_400_file_list.tab

cp  ~/.mergedpipeline.cnfg chi2_quick_align.cnfg
# altered to just do alignments
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs papadum_400_file_list.tab --output chi2_aligns --config chi2_quick_align.cnfg --reference /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa --threads 20

# Had to merge them due to errors
perl ~/perl_toolchain/sequence_data_scripts/bamMergeUtility.pl -d ./chi2_aligns/papadum -o ./chi2_aligns/chi2_papadum_merged.bam

perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi2.tab -t chi2_aligns/chi2_papadum_merged.bam -o ars1_gaps_on_chi2.depth.checked.tab


### Testing CHIR_1.0
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs papadum_400_file_list.tab --output chi1_aligns --reference ../CHIR_1.0_fixed.fa --config chi2_quick_align.cnfg --threads 15

perl ~/perl_toolchain/sequence_data_scripts/bamMergeUtility.pl -d chi1_aligns/papadum/ -o chi1_aligns/ars1_chi1_merged.bam

perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o ../papadum-v13.full.fa -s ../CHIR_1.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d ars1_gaps_on_chi1.newlogic.tab

perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi1.newlogic.tab -t chi1_aligns/ars1_chi1_merged.bam -o ars1_gaps_on_chi1.newlogic.checked.tab
# other checking criteria
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi1.newlogic.tab -t chi1_aligns/ars1_chi1_merged.bam -o ars1_gaps_on_chi1.newlogic.checked.2.tab -f ../CHIR_1.0_fixed.fa -v /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi2.newlogic.tab -t chi2_aligns/chi2_papadum_merged.bam -o ars1_gaps_on_chi2.newlogic.checked.2.tab -f /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -v /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz

# There were still some large gaps that were supposedly filled in the BGI assembly but were present in ours
# This needs to be checked
grep 'FullClose' ars1_gaps_on_chi1.newlogic.checked.tab | perl -lane 'if($F[4] > 30){system("samtools faidx /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/CHIR_1.0_fixed.fa $F[3]");}' > chi1_large_gaps_doublecheck.fa
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz chi1_large_gaps_doublecheck.fa > chi1_large_gaps_doublecheck.sam


```

#### Checking ARS1 filled gaps from CHIR_2.0

I wrote a script to process the data and identify misassemblies within closed gaps.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check

```bash
# CHI2 gaps closed
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g chi2/papadumv13_gap_fills.chi2.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o chi2/papadumv13_gap_fills.chi2.checked.depth.tab
Tested: 59737
Skipped: 9256

# CHI1 gaps closed
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ../papadumv13_gap_fills.old.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o ars1_chi1_filled_gaps.check.tab
Tested: 242888
Skipped: 13876

# CHI2.0 gaps fully closed
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi2/papadumv13_gap_fills.chi2.checked.depth.tab -c 0 -m

# CHI1.0 gaps fully closed
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f ars1_chi1_filled_gaps.check.tab -c 0

# New logic
# New logic script
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -o ../papadum-v13.full.fa -s /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d ars1_gaps_on_chi2.newlogic.tab

perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi2.newlogic.tab -t chi2_aligns/chi2_papadum_merged.bam -o ars1_gaps_on_chi2.newlogic.checked.tab
```

|Entry     | CHIR_1|CHIR_2|
|:---------|------:|-----:|
|FullClose | 208401| 48298|
|GAP       |  33406| 11078|
|Large     |   1081|   361|


It looks really good, but there were no filled gaps with less than 5 reads worth of coverage? I find that a bit hard to believe... Let's test this.

Test intervals:
**cluster_1:1711-1967**
**cluster_1:6324-6507**

The first should have several "bad areas" and the second should have good coverage.

```bash
echo -e "Closed\tYup\t1000\t2000\t1000\tcluster_1\tcluster_1:1711-1967\tcluster_1:1711-1967\t200\t1\t200\nClosed\tYup\t2000\t3000\t1000\tcluster_1\tcluster_1:6324-6507\tcluster_1:6324-6507\t200\t1\t200" > test_dataset.out.tab

perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g test_dataset.out.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o test_dataset.out.checked.tab
```

OK, my second set of criteria include the following:
* Gap flanking sequence must map unambiguously within 100kb of each other
	* Gap flanking sequence aligning > 100kb are "too **Large**"
	* Gap flanking sequence that has missing alignments, or 10% of the read as soft/hard clipped are "**Unmapped**"
* Any number of ambiguous bases ("N") are considered "**Gaps**"
* All else are considered "Closed" and are subject to the following checks
* Next, I check sequence-based coverage to confirm closure
	* If the samtools depth estimate for any bases across the region are < 5, read is considered a "**GAP**" (cryptic misassembly)
	* If the gap region has any ambiguous bases, it is considered a "**Gap**"
	* Next I check the full gap region by pulling the full sequence of the region in the target assembly and aligning it back to the original assembly
		* If you map the full gap closed region from the target assembly (ie. CHIR_2.0) back to the original assembly (ie. ARS1) and the entire region aligns to the assembly without ambiguity (no hard clipping, full length alignment), then the gap region is considered "**Ambiguous**" (most likely due to repetitive structure)
		* All remaining gap regions are labeled as "**FullClose**" entries and likely represent closure of the gap


|Entry     |CHIR_2|CHIR_1|
|:---------|-----:|-----:|
|Ambiguous |    12|    14|
|FullClose |     4|     4|
|GAP       |   141|    94|
|Large     |    85|    49|
|Trans     |   100|    98|
|Unmapped  |   254|   337|


## Rerunning the gap closure pipeline on CHIR1 and CHIR2

Now that I think that I have fair criteria, let's see how frequently ARS1 closes gaps in CHIR2 and CHIR1.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/chi1

```bash
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -s ../papadum-v13.full.fa -o ../CHIR_1.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d chi1_gaps_on_ars1.newlogic.tab ; echo "done first"; perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g chi1_gaps_on_ars1.newlogic.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o chi1_gaps_on_ars1.newlogic.checked.tab -v ../CHIR_1.0_fixed.fa -f /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz
```

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/chi2

```bash
perl ~/perl_toolchain/assembly_scripts/identifyFilledGaps.pl -s ../papadum-v13.full.fa -o /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -g ~/GetMaskBedFasta/store/GetMaskBedFasta.jar -j ~/jdk1.8.0_05/bin/java -d chi2_gaps_on_ars1.newlogic.tab ; echo "done first"; perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g chi2_gaps_on_ars1.newlogic.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o chi2_gaps_on_ars1.newlogic.checked.tab -v /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa -f /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz
# New criteria: only unmapped gap regions are ambiguous
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi2_gaps_on_ars1.newlogic.checked.2.tab -c 0 -m
```
|Entry     | Count|
|:---------|-----:|
|Ambiguous |    23|
|FullClose | 41306|
|GAP       |  9718|
|Gap       |   598|
|Large     |   264|

```bash
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g chi1_gaps_on_ars1.newlogic.tab -t /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam -o chi1_gaps_on_ars1.newlogic.checked.2.tab -v ../CHIR_1.0_fixed.fa -f /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi1_gaps_on_ars1.newlogic.checked.2.tab -c 0 -m
```
|Entry     |  Count|
|:---------|------:|
|FullClose | 153975|
|GAP       |  26090|
|Gap       |   2282|
|Large     |    568|

#### Testing ARS1 gaps closed based on the same logic

```bash
perl ~/perl_toolchain/assembly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi2.newlogic.tab -t chi2
_aligns/chi2_papadum_merged.bam -o chi2_gaps_on_ars1.newlogic.checked.3.tab -f /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixe
d.fa -v /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz; echo "done 1"; rm temp.check.fa; perl ~/perl_toolchain/ass
embly_scripts/checkIdentifiedGapPairedEndDisc.pl -g ars1_gaps_on_chi1.newlogic.tab -t chi1_aligns/ars1_chi1_merged.bam -o chi1_gaps_on_ars1.newlog
ic.checked.3.tab -f ../CHIR_1.0_fixed.fa -v /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi2_gaps_on_ars1.newlogic.checked.3.tab -c 0 -m
```

|Entry     | Count|
|:---------|-----:|
|FullClose |    16|
|GAP       |    49|
|Gap       |    92|
|Large     |    85|

```bash
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chi1_gaps_on_ars1.newlogic.checked.3.tab -c 0 -m
```

|Entry     | Count|
|:---------|-----:|
|FullClose |    18|
|GAP       |    19|
|Gap       |    75|
|Large     |    49|

```bash
# Checking alignments in order to resolve fullclosed regions 
#CHIR1
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz temp.check.fa > chi1_fullclose_check.sam
perl -lane 'if($F[0] =~ /^@/){next;}else{print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]";}' < chi1_fullclose_check.sam
# 8 alignments were to unplaced contigs and should be ambiguous
# 1 alignment was to a scaffold -- may be worth a follow-up later
# NC_022298.1:45519413-45520211   16      scaffold_324    12445   38      799M  Orig:cluster_6:47855013-47855650     637

grep 'FullClose' chi2_gaps_on_ars1.newlogic.checked.3.tab | perl -lane 'system("samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa $F[3] >> chir2_ars1_gapregions.fa");'
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz chir2_ars1_gapregions.fa > chir2_ars1_gapregions.sam
perl -lane 'if($F[0] =~ /^@/){next;}else{print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]";}' < chir2_ars1_gapregions.sam
# 7 alignments were to unplaced contigs or were unmapped and should be ambiguous
# Again, another alignment to the same scaffold that should be followed up
# CM001715_2:49090466-49091264    16      scaffold_324    12445   38      799M Orig: cluster_6:47855013-47855650     637

# Used this one-liner to distinguish between actual gap regions and ambiguous small regions that were not able to be checked
perl -lane 'if($F[0] eq "Closed" && $F[10] > $F[8] && $F[8] != -1 && $F[10] < 100000){print $_;}' < ../chi1/chi1_gaps_on_ars1.newlogic.tab
```

#### Complete gap summary for ARS1 gap fills

|Entry     | CHIR2|  CHIR1|
|:---------|-----:|------:|
|Ambiguous |    23|      0|
|FullClose | 41306| 153975|
|Unconfirm | 18113|  88293|
|Mutual Gap|    29|     52|
|Large     |   264|    568|
|Trans     |  9252|  13853|
|Unmapped  |     6|     23|

#### Complete gap summary for BGI filling of ARS1 gaps

|Entry     | CHIR2| CHIR1|
|:---------|-----:|-----:|
|Ambiguous |     5|     8|
|FullClose |    11|    10|
|Unknown   |    63|    30|
|Mutual Gap|    78|    64|
|Large     |    85|    49|

#### Common ARS1 gaps filled by both assemblies

cluster_15:47423320-47423344
cluster_15:73919205-73919210
cluster_25:12174627-12174632
cluster_4:101582232-101582256
cluster_6:47855013-47855650
cluster_7:71064191-71065249
cluster_7:88496162-88496186
cluster_8:106614827-106614851
cluster_9:87048108-87048132

#### Identifying ARS1 + CHI2, CHI1 filled gaps

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/gap_check/chi1

```bash
# CHI1 gaps filled by CHI2
grep 'Closed' ../chi1_to_chi2_gap_summary.tab | perl -lane 'if($F[10] - $F[8] < 1){($s1, $e1) = $F[6] =~ /.+:(\d+)-(\d+)/; ($s2, $e2) = $F[7] =~ /.+:(\d+)-(\d+)/; @h = ($s1, $s2, $e1, $e2); @h = sort{$a <=> $b}@h; print "$F[5]\t$h[1]\t$h[2]\t$F[1]:$F[2]-$F[3]";}' | sortBedFileSTDIN.pl > chi1_gaps_filled_chi2_without_nbases.bed

grep 'Closed' ../papadumv13_gap_fills.old.tab | perl -lane 'if($F[10] - $F[8] < 1){($s1, $e1) = $F[6] =~ /.+:(\d+)-(\d+)/; ($s2, $e2) = $F[7] =~ /.+:(\d+)-(\d+)/; @h = ($s1, $s2, $e1, $e2); @h = sort{$a <=> $b}@h; print "$F[5]\t$h[1]\t$h[2]\t$F[1]:$F[2]-$F[3]";}' | sortBedFileSTDIN.pl > chi1_gaps_filled_ars1_without_nbases.bed

# Simple venn comparison
cat chi1_gaps_filled_ars1_without_nbases.bed | cut -f4 > chi1_gaps_filled_ars1_without_nbases.list
cat chi1_gaps_filled_chi2_without_nbases.bed | cut -f4 > chi1_gaps_filled_chi2_without_nbases.list

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl chi1_gaps_filled_ars1_without_nbases.list chi1_gaps_filled_chi2_without_nbases.list
	File Number 1: chi1_gaps_filled_ars1_without_nbases.list
	File Number 2: chi1_gaps_filled_chi2_without_nbases.list
	Set     Count
	1       84294   <- unique to our assembly/trans in chi2
	1;2     157884  <- shared
	2       2415	<- Trans in our assembly

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_2 chi1_gaps_filled_ars1_without_nbases.list chi1_gaps_filled_chi2_without_nbases.list > chi1_gaps_filled_both_without_nbases.list

# Retrieving closure sequence for confirmation
perl -lane 'if($F[2] - $F[1] < 10000){system("samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/BGI_chi_2/CHIR_2.0_fixed.fa $F[0]:$F[1]-$F[2] >> chi1_gaps_filled_chi2_without_nbases.fa");}' < chi1_gaps_filled_chi2_without_nbases.bed
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz chi1_gaps_filled_chi2_without_nbases.fa > chi1_gaps_filled_chi2_without_nbases.sam

# Checking mapped alignments without Mapq filter
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] != 2064 && $F[1] != 2048 && $F[2] ne "*"){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < chi1_gaps_filled_chi2_without_nbases.sam | wc -l
143,135 (89.29%)

# Now with a mapq filter of 10
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] != 2064 && $F[1] != 2048 && $F[2] ne "*" && $F[4] > 10){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]";}' < chi1_gaps_filled_chi2_without_nbases.sam | wc -l
113,219 (70.63%)

# Now to do a venn comparison
perl -lane 'if($F[0] =~ /^@/){next;} if($F[1] != 2064 && $F[1] != 2048 && $F[2] ne "*"){print "$F[0]";}' < chi1_gaps_filled_chi2_without_nbases.sam > chi1_gaps_filled_chi2_without_nbases.crossalign.pass.list
# Converting back to CHI1 coords
perl -e 'chomp(@ARGV); %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; @s = split(/\t/); $h{$s[0]}->{$s[1]}->{$s[2]} = $s[3];} close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; ($chr, $s, $e) = $_ =~ /(.+):(\d+)-(\d+)/; if(exists($h{$chr}->{$s}->{$e})){print $h{$chr}->{$s}->{$e}; print "\n";}else{print STDERR "Could not find $_!\n";}} close IN;' chi1_gaps_filled_chi2_without_nbases.bed chi1_gaps_filled_chi2_without_nbases.crossalign.pass.list > chi1_gaps_filled_chi2_without_nbases.crossalign.pass.converted.list

# Now to grep out the fullclose and partial closed gaps from CHI1 on ARS1
grep 'FullClose' chi1_gaps_on_ars1.newlogic.checked.2.tab | perl -lane 'print $F[1];' > chi1_gaps_on_ars1.newlogic.fullclose.wgs.list
grep -v 'Large' chi1_gaps_on_ars1.newlogic.checked.2.tab | perl -lane 'print $F[1];' > chi1_gaps_on_ars1.newlogic.allbutlarge.wgs.list

# OK, we're ready for venn comparisons
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl chi1_gaps_filled_chi2_without_nbases.list chi1_gaps_filled_ars1_without_nbases.list chi1_gaps_filled_chi2_without_nbases.crossalign.pass.converted.list chi1_gaps_on_ars1.newlogic.fullclose.wgs.list chi1_gaps_on_ars1.newlogic.allbutlarge.wgs.list
File Number 1: chi1_gaps_filled_chi2_without_nbases.list
File Number 2: chi1_gaps_filled_ars1_without_nbases.list	(n = 242178)
File Number 3: chi1_gaps_filled_chi2_without_nbases.crossalign.pass.converted.list
File Number 4: chi1_gaps_on_ars1.newlogic.fullclose.wgs.list
File Number 5: chi1_gaps_on_ars1.newlogic.allbutlarge.wgs.list

Set     Count
1       163
1;2     3101	<- same with chi2, but didn't pass either test (1.28%)
1;2;3   35306	<- alignment pass filter only (14.578%)
1;2;3;4;5       94212	<- confirmed in both assemblies; both methods (38.90%)
1;2;3;5 11343   <- confirmed in both by only cross-align (4.684%)
1;2;4;5 11956   <- confirmed in both by only wgs filter (4.937%)
1;2;5   1966    <- unconfirmed in both by tentative wgs filter (0.811%)
1;3     2237    <- shouldn't exist -- will check (they're Trans chr aligns)
1;3;5   13		<- shouldn't exist -- will check (they're gaps in our assembly)
1;5     2
2       21476	<- failed all tests (8.867%)
2;4;5   47807   <- confirmed in ARS1 by WGS gap fill (19.74%)
2;5     15011   <- unconfirmed in ARS1 by wgs gap fill (6.20%)
5       37
```
That makes 82.84% confirmed in our assembly (200624 / 242178 total closed gaps). 94,212 (38.90%) were confirmed by both cross assembly alignment and WGS alignment. 153,975 were confirmed by WGS methods. 143,135 were confirmed by cross-assembly alignment.
 
<a name="fosmid"></a>
## Fosmid Check
Steve has done the alignments but there are some issues. I need to remove duplicates and then merge together fosmid reports. First, lets format and fix the data.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/fosmid_check

```bash
for i in CHIR-1 CHIR-2 papadum-lachesis papadum-v03 papadum-v04 papadum-v05 papadum-v13; do echo $i; mkdir $i; wkdir=/mnt/nfs/nfs2/GoatData/Fosmid/${i}; for b in $wkdir/*.bam; do bname=`basename $b`; echo $bname; ~/jdk1.8.0_05/bin/java -jar ~/picard-tools-1.85/picard.jar MarkDuplicates I=$b O=${i}/${bname} M=${i}/${bname}.mkdups REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200; done; done

for i in ./*/*.bam; do echo $i; samtools index $i; done

# I'm getting allot of seg faults. Maybe I should run the libraries separately?
perl reformat_bams.pl

for i in CHIR-1 CHIR-2 papadum-lachesis papadum-v03 papadum-v04 papadum-v05 papadum-v13; do bfai=$i/${i}_lumpy_sv_BfaI_output.vcf; rsai=$i/${i}_lumpy_sv_RsaI_output.vcf; out=$i/${i}_sep_libs_nonoverlap.bed; for b in $bfai $rsai; do perl -lane 'if($F[0] =~ /^#/){next;}else{($type) = $F[7] =~ /SVTYPE=(.+)\;/; $F[0] =~ s/\|/_/g; $e = $F[1] + 1; print "$F[0]\t$F[1]\t$e\t$type\n";}' < $b; done > $out; done

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f ./CHIR-1/CHIR-1_sep_libs_nonoverlap.bed,./CHIR-2/CHIR-2_sep_libs_nonoverlap.bed,./papadum-lachesis/papadum-lachesis_sep_libs_nonoverlap.bed,./papadum-v03/papadum-v03_sep_libs_nonoverlap.bed,./papadum-v04/papadum-v04_sep_libs_nonoverlap.bed,./papadum-v05/papadum-v05_sep_libs_nonoverlap.bed,./papadum-v13/papadum-v13_sep_libs_nonoverlap.bed -c 3 -m
File1:  ./CHIR-1/CHIR-1_sep_libs_nonoverlap.bed
File2:  ./CHIR-2/CHIR-2_sep_libs_nonoverlap.bed
File3:  ./papadum-lachesis/papadum-lachesis_sep_libs_nonoverlap.bed
File4:  ./papadum-v03/papadum-v03_sep_libs_nonoverlap.bed
File5:  ./papadum-v04/papadum-v04_sep_libs_nonoverlap.bed
File6:  ./papadum-v05/papadum-v05_sep_libs_nonoverlap.bed
File7:  ./papadum-v13/papadum-v13_sep_libs_nonoverlap.bed
```

|CHIR1 | Count|CHIR2 | Count|LachOn| Count|PBctg | Count|PBBNG | Count|PBnLaC| Count|ARS1  | Count|
|:-----|-----:|:-----|-----:|:-----|-----:|:-----|-----:|:-----|-----:|:-----|-----:|:-----|-----:|
|BND   |  4342|BND   |   840|BND   |   870|BND   |  1158|BND   |   920|BND   |   664|BND   |   456|
|DEL   |   197|DEL   |   144|DEL   |   162|DEL   |   130|DEL   |   134|DEL   |   143|DEL   |   144|
|DUP   |    14|DUP   |    12|DUP   |    31|DUP   |     4|DUP   |     5|DUP   |    33|DUP   |     8|
|INV   |     5|INV   |    10|INV   |    22|INV   |     0|INV   |     2|INV   |    13|INV   |    10|
|100MB | 155.8|100MB |  34.4|100MB |  37.1|100MB |  44.2|100MB |  36.3|100MB |  29.2|100MB |  21.1|


OK, the lumpy stats are done. Let's calculate the number of reads that mapped.

```bash
for i in ./*/*.bam; do unmapped=`samtools idxstats ./CHIR-1/CHIR-1-goat-fosmid-BfaI.bam | grep "*" | cut -f4`; echo -e "$i\t$unmapped"; done

# they were all the same: 88594
```

OK, let's turn this into a consolidating table with an "error per mb" estimate. Done (above). Now to check repetitive regions as a response to reviewer #2.

I need to pull out the Del and Dup sequences and check them against repetitive regions. All of the BND signal appears to be from fosmids that bind to unplaced contigs, suggesting that they should be scaffolded where their anchor reads are attached. 

```bash
# Pulling only events that have more than one PE support and are less than 1 mb
for i in ./papadum-v13/*I_output.vcf; do perl -lane 'if($F[0] =~ /^#/){next;}else{($type, $end, $pe) = $F[7] =~ /SVTYPE=(.{3,5});.*END=(\d+);.*PE=(\d{1,4});/; if($pe > 0 && $type ne "BND" && $end - $F[1] < 1000000){print "$F[0]\t$F[1]\t$end\t$type\t$pe";}}' < $i; done | sortBedFileSTDIN.pl > papadum-v13/del_dup_inv_events.bed

# Intersecting with repeats
intersectBed -a papadum-v13/del_dup_inv_events.bed -b ../repeat_analysis/papadum-v13.full.repeatmask.extend.75thresh.out.bed -wb | perl -lane 'if($F[11] > 1000){print "$F[5]\t$F[6]\t$F[7]\t$F[10]\t$F[11]\t$F[3]";}' | wc -l
21

intersectBed -a papadum-v13/del_dup_inv_events.bed -b ../repeat_analysis/papadum-v13.full.repeatmask.extend.75thresh.out.bed -wb | perl -lane 'if($F[11] > 1000){print "$F[5]\t$F[6]\t$F[7]\t$F[10]\t$F[11]\t$F[3]";}'
cluster_20      62619994        62626357        LINE/L1 6888    DEL
cluster_14      62117816        62124142        LINE/L1 8390    INV
cluster_10      77586122        77592441        LINE/L1 8154    DUP
cluster_10      77794417        77796334        Satellite/centr 1965    DUP
cluster_10      77795261        77796535        Satellite/centr 1262    DUP
cluster_10      77797410        77799859        Satellite/centr 2443    DUP
cluster_10      77801227        77802588        Satellite/centr 1360    DUP
cluster_10      77803542        77805200        Satellite/centr 1686    DUP
cluster_10      77808725        77811352        Satellite/centr 2602    DUP
cluster_10      77827930        77829604        Satellite/centr 1684    DUP
cluster_10      77830878        77835736        Satellite/centr 4895    DUP
cluster_10      77834661        77835739        Satellite/centr 1402    DUP
cluster_10      77837035        77838709        Satellite/centr 1398    DUP
cluster_10      77839663        77841023        Satellite/centr 1386    DUP
cluster_10      77720101        77726629        LINE/L1 8381    DUP
cluster_19      68865813        68867308        Simple_repeat   1524    DUP
cluster_19      69035001        69038800        Simple_repeat   3888    DUP
cluster_19      69197571        69199724        Simple_repeat   2182    DUP
cluster_19      68865813        68867308        Simple_repeat   1524    DUP
cluster_19      69035001        69038800        Simple_repeat   3888    DUP
cluster_19      69197571        69199724        Simple_repeat   2182    DUP

# There were only 5 dels/dups/invs that had repetitive content
cluster_10      77479823        77856293        DUP     5
cluster_19      68374696        68587239        DEL     5
cluster_19      68410541        69257800        DUP     4
cluster_19      68410663        69244694        DUP     4
cluster_20      62517172        62761341        DEL     5
cluster_14      62101583        62183931        INV     4
```


#### Comparison of fosmid end discrepencies to gaps

I wanted to see how many of our gaps were predicted deletions/inversions/duplications based on the fosmid data.

```bash
intersectBed -a papadum-v13/del_dup_inv_events.bed -b ../gap_check/papadum-v13.full.gaps.bed -wa | uniq | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 3 -m
```
#### Gaps compared to del, inv and dup calls

|Entry | Count|
|:-----|-----:|
|DEL   |     6|
|DUP   |     3|
|INV   |     5|

That's 14 gaps (sum - 1 for a duplicate duplication on cluster_19).

```bash
intersectBed -a papadum-v13/del_dup_inv_events.bed -b ../gap_check/chi1/chi1_gaps_filled_ars1_without_nbases.bed -wa | uniq | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 3 -m
```

#### Gaps ARS1 filled in CHI1 compared to del, inv and dup calls

|Entry | Count|
|:-----|-----:|
|DEL   |     6|
|DUP   |     5|
|INV   |     8|

That's 18 gaps (sum - 1 for the same duplicate on cluster_19). Most of them are the same del/dup/inv calls listed above (only 4 new regions).

## Quick unplaced analysis

I just need to identify the reason why our degenerate contigs exist. I suspect that most are just very strong repeats.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/fosmid_check

```bash
# Degenerate count:
grep unplaced ../gap_check/papadum-v13.full.fa.fai | wc -l
29,315

# Checking how many have > 75% length devoted to repeats
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0.75){print "$up\n";}}' | wc -l
21,820

# How many have SOME satellite sequence
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed | grep Satellite |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0){print "$up\n";}}' | wc -l
25,023

# How many have centromeric/satellite sequence over 25% 
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed | grep Satellite |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0.25){print "$up\n";}}' | wc -l
24,500

# Satellite over 75% of the contig
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed | grep Satellite |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0.75){print "$up\n";}}' | wc -l
19,295

# Zero read alignments
samtools idxstats /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam | grep 'unplaced' | perl -lane 'if($F[2] == 0){print $F[0];}' | wc -l
123
```

OK, so most, right off the bat, are highly repetitive. Let's do a Venn analysis to confirm.

```bash
grep unplaced ../gap_check/papadum-v13.full.fa.fai | perl -lane 'print $F[0];' > ars1_degen_master.list
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0.75){print "$up\n";}}' > ars1_degen_rep_gt75.list
perl -e 'chomp(@ARGV); open(IN, "grep unplaced ../gap_check/papadum-v13.full.fa.fai |"); %lens; while(<IN>){chomp; @s = split(/\t/); $lens{$s[0]} = $s[1];} close IN; open(IN, "grep unplaced ../repeat_analysis/papadum-v13.full.repeatmask.out.bed | grep Satellite |"); %reps; while(<IN>){chomp; @s = split(/\t/); $reps{$s[0]} += $s[2] - $s[1];} close IN; foreach my $up (keys(%reps)){ my $ratio = $reps{$up} / $lens{$up}; if($ratio > 0){print "$up\n";}}' > ars1_degen_rep_satellite.list
samtools idxstats /mnt/nfs/nfs2/GoatData/Ilmn/papadum-v13/bwa-out/Goat-Ilmn-Freezev13.bam | grep 'unplaced' | perl -lane 'if($F[2] == 0){print $F[0];}' > ars1_degen_zero_reads.list

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl ars1_degen_master.list ars1_degen_rep_satellite.list ars1_degen_rep_gt75.list ars1_degen_zero_reads.list
File Number 1: ars1_degen_master.list
File Number 2: ars1_degen_rep_satellite.list
File Number 3: ars1_degen_rep_gt75.list
File Number 4: ars1_degen_zero_reads.list
Set     Count
1       3494	<- The only set that is not accounted for here
1;2     3886
1;2;3   21125
1;2;3;4 5
1;2;4   7
1;3     687
1;3;4   3
1;4     108

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1 ars1_degen_master.list ars1_degen_rep_satellite.list ars1_degen_rep_gt75.list ars1_degen_zero_reads.list > ars1_unaccounted_for.list
```

## GWAS analysis and plot generation

I need to select examples of remapped probes that show significant effects in Heather's GWAS. Let's start with the low-hanging fruit by taking the unmapped probes from the original CHIR_1.0 and checking their significance in the Peulh data.

In the remapped 56k chip, the first 1409 markers are unmapped in CHIR_1.0.

> pwd: /home/dbickhart/share/goat_assembly_paper/GWAS

```bash
head -n 1409 ../snp_probes/ADAPTmap_ARS1_updated_ignoreX.map | perl -lane 'if($F[0] == 0){next;}else{print "$F[0]\t$F[3]\t$F[1]";}' | perl -e 'chomp(@ARGV); %h; while(<STDIN>){chomp; @s = split(/\t/); $h{$s[2]} = "$s[0]:$s[1]";} open(IN, "< $ARGV[0]"); print "snp\tchr\tpos\tdistance\tpvalue\t-log10p\tFDR\tproportionVar\tcallRate\n"; while(<IN>){chomp; @s = split(/,/); if(exists($h{$s[0]}) && $s[5] > 2){print "$s[0]\t$s[1]\t$s[2]\t$s[3]\t$s[4]\t$s[5]\t$s[10]\t$s[11]\t$s[14]\n";}} close IN;' Peulh\ AGIN\ assembly.csv
snp     chr     pos     distance        pvalue  -log10p FDR     proportionVar   callRate
snp31007-scaffold3424-53240     3       103188580       0       0.00125982960141677     2.89968819152438        0.63533700850272        0.369811375711095       1
snp50485-scaffold7244-13356     6       3106462 0       0.00799190836500886     2.09734950436155        0.884076934167077       0.268327332225227       1
snp54629-scaffold834-105661     6       3237591 0       0.00564371471500465     2.24843494758081        0.806408447847567       0.28835576489127        1
snp40027-scaffold5107-301       7       107694942       0       0.000741254328141377    3.13003275779978        0.577717899776732       0.396680780134015       1
snp11230-scaffold1401-556817    12      87109949        0       2.03317637356759e-006   5.69182494563904        0.0522922797399718      0.632392457012523       1
snp13191-scaffold1505-153770    26      45556525        0       0.00916246880419353     2.03798749100726        0.916942087196324       0.260347300254181       1
snp4713-scaffold1153-307649     28      42140721        0       7.65446357828852e-006   5.11608523908096        0.0984344880008958      0.588609980477254       1
snp4712-scaffold1153-270557     28      42177859        0       7.65446357828852e-006   5.11608523908096        0.0787475904007166      0.588609980477254       1
snp24521-scaffold249-93325      29      38814227        0       0.00855510445876449     2.06777468330721        0.89809391480487        0.26435926350978        1
```

OK, we've got one marker to work with right now that is barely before the FDR cutoff: snp11230-scaffold1401-556817

Let's grep out all of the markers on chr12 from the dataset. It looks like this marker is the second to last on the chromosome and is flanked by a long distance from the end of the chr.

```R
data <- read.delim("Peulh AGIN assembly.csv", sep=",", header=TRUE)
chr12_markers <- data[data$Chromosome == 12, ]
chr12_markers$X.log10.P.Value. <- as.numeric(levels(chr12_markers$X.log10.P.Value.))[chr12_markers$X.log10.P.Value.]
plot(x = chr12_markers$Position, y =chr12_markers$X.log10.P.Value., col = ifelse(chr12_markers$Marker == "snp11230-scaffold1401-556817", "red", "black"))
```

Let's combine the data to draw out regions that have changed between the two assemblies.




