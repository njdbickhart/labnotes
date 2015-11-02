# Resource testing alignment software
---
*11/2/2015*

These are my notes on generating runtime reports for Paul to compare with his variant calling software.

First, I need to modify the test fastq file that Paul generated. I'll do this with a one-shot script.

This is Paul's format:
> @chr 1  54897247         1           1 1   150 F

I'll preserve all of the numbers but remove the "F" and "R" designation of alignment for now.

> 3850: /work3/paul/UMD
> changeFqNames.pl

```perl
#!/usr/bin/perl

use strict;

my $usage = "perl $0 <input fq file 1> <input fq file 2> <output name file>  <output name base>\n";

chomp(@ARGV);
unless(scalar(@ARGV) == 4){
        print $usage;
        exit;
}

open(my $OUT, "> $ARGV[2]");
open(my $FQ1, "< $ARGV[0]");
open(my $FQ2, "< $ARGV[1]");

# Open revised fastq files
open(my $OQ1, "> $ARGV[3].1.fq");
open(my $OQ2, "> $ARGV[3].2.fq");
my $rcount = 1;
while(my $n1 = <$FQ1>){
        my $s1 = <$FQ1>;
        my $p1 = <$FQ1>;
        my $q1 = <$FQ1>;

        my $n2 = <$FQ2>;
        my $s2 = <$FQ2>;
        my $p2 = <$FQ2>;
        my $q2 = <$FQ2>;

        chomp($n1, $s1, $q1);
        chomp($n2, $s2, $q2);

        # Replace whitespace with underscores
        $n1 =~ s/\s+/_/g;
        $n2 =~ s/\s+/_/g;

 		# Take care of the case where the chromosome name was space delimited
        $n1 =~ s/\@//;
        $n2 =~ s/\@//;
        $n1 =~ s/chr_/chr/;
        $n2 =~ s/chr_/chr/;

        # Print read names to output file:
        print {$OUT} "$n1\t$n2\n";

        my $tile = int(rand(99999));
        my $id = int(rand(9999));
        my @rd1segs = split(/_/, $n1);
        my @rd2segs = split(/_/, $n2);
        $rd1segs[0] =~ s/\@//; $rd2segs[0] =~ s/\@//;
        #my $tn1 = "\@EASE:$rd1segs[0]$rd1segs[1]:$rd1segs[1]:$rd1segs[2]:$rd1segs[3]:$rd1segs[4]:$rcount 1:N:0:ATCACG";
        #my $tn2 = "\@EASE:$rd2segs[0]$rd2segs[1]:$rd2segs[1]:$rd2segs[2]:$rd2segs[3]:$rd2segs[4]:$rcount 2:N:0:ATCACG";

        my $tn1 = "\@EASE:$tile:$id:$rcount 1:N:0:ATCACG\n";
        my $tn2 = "\@EASE:$tile:$id:$rcount 2:N:0:ATCACG\n";

        $rcount++;
        print $OQ1 "$tn1\n$s1\n+\n$q1\n";
        print $OQ2 "$tn2\n$s2\n+\n$q2\n";
}

close ($OQ1, $OQ2, $FQ1, $FQ2, $OUT);
exit;
```

So, the script will pull the old names, remove the whitespace and replace it with underscores, print the full names to a separate tab file and then print new read headers that shouldn't be removed by BWA during alignment.

Hopefully it works! I'll have to run it on a directory that I have permissions though!

> 3850: /seq1/bickhart/paul_alignment

```bash
# running the script
perl ~/changeFqNames.pl /work3/paul/UMD/segments1.fq /work3/paul/UMD/segments2.fq original_segments_names.tab segments_revised
```

OK! Everything looks in order. I'm now going to test the files using my [bwa profile script](https://github.com/njdbickhart/perl_toolchain/blob/master/simulations/runBWAAlignmentResourceTest.pl). 

```bash
perl ~/perl_toolchain/simulations/runBWAAlignmentResourceTest.pl -f segments_revised.1.fq -r segments_revised.2.fq -l segments_test.log -g /seq1/reference/umd3_kary_unmask_ngap.fa -o segments_revised

