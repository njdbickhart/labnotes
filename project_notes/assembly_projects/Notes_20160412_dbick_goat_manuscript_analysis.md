# Goat assembly analysis
---
*4/12/2016*

These are my notes on the conclusion of the analysis of the goat reference genome. I wanted to segregate these from the main file, which was becoming quite lengthy!

## Table of Contents
* [Spearman rank correlation of Lachesis](#correlation)
* [Preparing data for NCBI](#ncbi)
* [Centromere repeat check](#centromere)

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
