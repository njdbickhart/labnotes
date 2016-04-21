# Goat assembly analysis
---
*4/12/2016*

These are my notes on the conclusion of the analysis of the goat reference genome. I wanted to segregate these from the main file, which was becoming quite lengthy!

## Table of Contents
* [Spearman rank correlation of Lachesis](#correlation)
* [Preparing data for NCBI](#ncbi)

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