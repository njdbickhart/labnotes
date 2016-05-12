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
		if($d->[3] != $mostchr || $d->[4] eq '.'){next;}
		push(@processed, [$d->[1], $d->[2] + $addl, $d->[4] + $addr]);
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

```

<a name="centromere"></a>
## Centromere repeat check

Serge identified several low complexity regions that may be centromeric. I'm interested to see if we've found it so I downloaded the consensus sequence from [an article](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-r10) that claims to have identified it (note: it's in the online supplementary notes). I'm going to pull the sequence from clusters 7 and 9 to see if it matches the centromere.

> Blade14: /mnt/iscsi/vnx_gliu_7/goat_assembly/centromere_check

```bash
# Just because it's easy, let's try bwa
bwa mem /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz centromere.fa > centromere.sam

samtools faidx /mnt/nfs/nfs2/GoatData/Goat-Genome-Assembly/Papadum-v13/papadum-v13.full.fa.gz cluster_7:112500444-112506199