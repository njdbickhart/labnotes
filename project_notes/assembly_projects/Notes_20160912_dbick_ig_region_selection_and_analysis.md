# IG BAC selection
---
*9/12/2016*

These are my notes on selecting BACs to sequence from the CHORI 240 library. I'm going to use one-shot scripts to try to automate this process and generate BAC end alignments across the entirety of Tim's new reference cattle genome.

I am using the 31Mar2016 assembly.

## Aligning the data

This was fairly straightforward. 

> Blade14: /mnt/iscsi/vnx_gliu_7/john_assembled_contigs

```bash
bwa mem cattle_31Mar2016_bfJmS.fasta /mnt/nfs/nfs2/dbickhart/dominette_ilm/chori_bac_end_clones_formatted.fasta.gz > chori_bac_end_aligns.sam

perl -e '$c = 0; $a = 0; while(<>){chomp; @s = split(/\t/); if($s[0] =~ /^@/){next;} if($s[2] eq "*"){$c++;}else{$a++;}} print "$c\n$a\n";' < chori_bac_end_aligns.sam
1388	No aligns
306126	Alignments

gunzip -c /mnt/nfs/nfs2/dbickhart/dominette_ilm/chori_bac_end_clones_formatted.fasta.gz | grep '>' | wc -l
304710

# We've got some split alignments to deal with, but the alignment rate is quite good
# (304710 - 1388) / 304710 = 99.54%
```

## Processing the alignments into tab delimited formats

Now, I need to generate tables so that I can cross-reference Ig annotations with the BAC list.

```bash
perl process_chori_clone_aligns.pl chori_bac_end_aligns.sam chori_bac_end_aligns.tab

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f chori_bac_end_aligns.tab -c 1
Entry   Count
1       32067
2       135271
```

The data looks good. Here's the one-shot script that I used to generate the table:

```perl
#!/usr/bin/perl
# This is a one shot script to process the alignments for the CHORI BACs
# BAC end read names follow this format: CH240_9A5.2
# Output file format: clone read_count read_num chr pos read_num chr pos

use strict;

my $usage = "perl $0 <chori bac sam> <output tab file>\n";
chomp(@ARGV);

unless(scalar(@ARGV) == 2){
        print $usage;
        exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open input sam!\n$usage";
open(my $OUT, "> $ARGV[1]");

my %data; # {clone}->{id} = [chr, pos]
while(my $line = <$IN>){
        chomp $line;
        if($line =~ /^@/){next;}
        my @segs = split(/\t/, $line);

        # Remove sub alignments
        if($segs[1] & 2048 == 2048){
                next;
        }

        # Process read name and store in data hash
        my ($clone, $rn) = $segs[0] =~ /(.+)\.(\d*)/;
        # Need to deal with empty read number fields
        if($rn eq '' || length($rn) == 0){
                if(exists($data{$clone})){
                        $rn = 2;
                }else{
                        $rn = 1;
                }
        }
        $data{$clone}->{$rn} = [$segs[2], $segs[3]];
}
close $IN;

# Process the data
print {$OUT} "clone\tread_count\tread_num\tchr\tpos\tread_num\tchr\tpos\n";
foreach my $key (sort {$a cmp $b} keys(%data)){
        my $rcount = scalar(keys(%{$data{$key}}));
        print {$OUT} "$key\t$rcount";
        foreach my $read (sort {$a <=> $b} keys(%{$data{$key}})){
                print {$OUT} "\t$read\t" . $data{$key}->{$read}->[0] . "\t" . $data{$key}->{$read}->[1];
        }
        print {$OUT} "\n";
}
close $OUT;
exit;
```

## Organizing data and selecting clones

From a prior email, I know that John was interested in the following scaffolds from the assembly:

| Scaffold | Gene | Scaffold_length | BAC_hits |
| :--- | :--- | ---: | ---: |
ScbfJmS_217  |    IGHA | 67999220 | 4120
ScbfJmS_2469 |  IGHV | 114817 | 6 
ScbfJmS_1121 |  IGKC | 98976055 | 5897
ScbfJmS_1160 |  IGLC | 2395294 | 353
ScbfJmS_2336 |  IGLC | 18658 | 0
ScbfJmS_1659 |  IGLV | 42106859 | 2545
ScbfJmS_1914 |  IGLV | 11008 | 0
ScbfJmS_2054 |  IGLV | 14030 | 0
ScbfJmS_2297 |  IGLV | 33706 | 2
ScbfJmS_2373 |  IGLV | 34368 | 2
ScbfJmS_514  |  IGLV | 51142 | 4
ScbfJmS_318  |  LRC | 63079628 | 3752
ScbfJmS_446  |  LRC | 2620727 | 128

Some of these will be like trying to hit a needle with another needle! Still, let's see if we've gotten any hits in the BAC ends.

```bash
for i in ScbfJmS_217 ScbfJmS_2469 ScbfJmS_1121 ScbfJmS_1160 ScbfJmS_2336 ScbfJmS_1659 ScbfJmS_1914 ScbfJmS_2054 ScbfJmS_2297 ScbfJmS_2373 ScbfJmS_514 ScbfJmS_318 ScbfJmS_446; do echo $i; grep -P "$i\t" chori_bac_end_aligns.tab | wc -l; done
	ScbfJmS_217 4120
	ScbfJmS_2469 6
	ScbfJmS_1121 5897
	ScbfJmS_1160 353
	ScbfJmS_2336 0
	ScbfJmS_1659 2545
	ScbfJmS_1914 0
	ScbfJmS_2054 0
	ScbfJmS_2297 2
	ScbfJmS_2373 2
	ScbfJmS_514 4
	ScbfJmS_318 3752
	ScbfJmS_446 128

```

OK! We're somewhat in business! Some of these scaffolds are so small as to contain ONLY the BAC ends, but they're a start!

Interestingly, they also bridge the gap to other scaffolds, as is the case here:

```bash
grep 'ScbfJmS_2297' chori_bac_end_aligns.tab
CH240_196P13    2       1       ScbfJmS_2297    10740   2       ScbfJmS_1659    41941553
CH240_231B16    2       1       ScbfJmS_2297    5662    2       ScbfJmS_1659    42001414
```

**NOTE: The RPCI-42 clone that we need is: 567N23**

## SNP selection and curation

I need to test the feasibility of selecting SNPs from the new cattle assembly vs the individual contigs. Let's prepare the contigs for alignment and then align the entirety of John's holstein to them for subsequent variant calling.

Here's what I need to do:
* organize the contigs into a larger fasta file
* align dominette's data to them and pilon correct them
* align john's holstein data to the uncorrected contigs
* prepare a table with distinguishing SNPs for each contig

Let's start with the preparation:

> Blade14: /mnt/iscsi/vnx_gliu_7/immune_grant/assembly_wgs_align

```bash
# Getting John's holsteins into a format that lets me run my alignment program on them
ls /mnt/cifs/bickhart-qnap/john_holstein/*/Project_AIPL/*/*.fastq.gz | perl -e '%h; while(<>){chomp; @b = split(/\//); @s = split(/[\_\.]/, $b[-1]); $h{$s[0]}->{$b[5]}->{$s[4]}->{$s[2]}->{$s[3]} = join("/", @b);} foreach my $an (keys(%h)){foreach my $flow (keys(%{$h{$an}})){ foreach my $split (keys(%{$h{$an}->{$flow}})){ foreach my $lane (keys(%{$h{$an}->{$flow}->{$split}})){ print $h{$an}->{$flow}->{$split}->{$lane}->{"R1"} . "\t" . $h{$an}->{$flow}->{$split}->{$lane}->{"R2"} . "\t$flow\t$an\n";}}}}' > john_animals_list.spreadsheet.simple

# I used vim to edit this list to remove the single file entries

# preparing the contig fastas for pilon correction
cat /mnt/nfs/nfs2/dbickhart/transfer/tim_assemblies/LIB144*.fasta > tims_fastas_vector_trimmed.fa
# I used vim to change the fasta headers into the actual clone names
bwa index tims_fastas_vector_trimmed.fa; samtools faidx tims_fastas_vector_trimmed.fa

# I'm just going to use a subsection of the reads from the dominette lung tissue list here
bwa mem tims_fastas_vector_trimmed.fa /mnt/nfs/nfs2/SequenceData/Dominette/Dominette_NextSeq_data/Run_1/LIB18483_S1_L003_R1_001.fastq.gz /mnt/nfs/nfs2/SequenceData/Dominette/Dominette_NextSeq_data/Run_1/LIB18483_S1_L003_R2_001.fastq.gz | samtools view -bS -F 4 - | samtools sort -o tims_fastas_vector_trimmed.bam -T tims -

```

I need a figure to show a good comparison with the old assembly. Let's do an alignment of the LRC cluster to get the positions of these genes on UMD and the canu assembly.

```bash
# First extract the fastas
samtools faidx tims_fastas_vector_trimmed.fa RP42-168O11_LRC RP42-141D20_LRC RPCI42_65A17_LRC > tims_unpolished_lrcs.fa

# Canu alignment
bwa mem /mnt/nfs/nfs2/dbickhart/dominette_asm/canu/topolish.filledWithCanuAndPBJelly.fasta.gz tims_unpolished_lrcs.fa > tims_unpolished_lrcs.sam

perl -lane 'if($F[0] =~ /^@/){next;}else{$e = 0; while($F[5] =~ /(\d+)(\D{1})/g){if($2 eq "M" || $2 eq "D" || $2 eq "S" || $2 eq "X"){$e += $1;}} $e += $F[3]; $l = $e - $F[3]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$e\t$l";}' < tims_unpolished_lrcs.sam

# Highlighting a particular region now
samtools faidx /mnt/nfs/nfs2/dbickhart/dominette_asm/canu/topolish.filledWithCanuAndPBJelly.fasta 18:63325473-63465228 > canu_assembly_lrc_region.fa
samtools faidx tims_fastas_vector_trimmed.fa RP42-168O11_LRC > RP42-168O11_LRC.fa
bwa mem /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa RP42-168O11_LRC.fa > RP42-168O11_LRC.sam

perl -lane 'if($F[0] =~ /^@/){next;}else{$e = 0; while($F[5] =~ /(\d+)(\D{1})/g){if($2 eq "M" || $2 eq "D" || $2 eq "S" || $2 eq "X"){$e += $1;}} $e += $F[3]; $l = $e - $F[3]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$e\t$l";}' < RP42-168O11_LRC.sam

samtools faidx /mnt/iscsi/vnx_gliu_7/reference/umd3_kary_unmask_ngap.fa chr18:63342820-63479594 > umd3_lrc_region.fa

samtools index tims_fastas_vector_trimmed.bam
~/jdk1.8.0_05/bin/java -Xmx80g -jar ~/pilon-1.16.jar --genome tims_fastas_vector_trimmed.fa --frags tims_fastas_vector_trimmed.bam --diploid --nostrays --threads 10 --output tims_fastas_vector_pilon
```
