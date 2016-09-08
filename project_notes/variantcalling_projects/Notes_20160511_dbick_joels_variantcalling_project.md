# Identifying variants within Joel's bulls
---
*5/11/2016*

These are my notes on the alignment and variant calling performed on Joel's bulls as part of our BARD funded grant.

## Table of Contents
* [Organizing the data](#organizing)
* [Generating 1000 bulls SNP and INDEL annotations](#onethousand)
	* [Output tab file columns](#outheads)
* [Calling SNPs and INDELs in Canadian and US Holstein data](#finertune)
* [Plan B: generating variant site probabilities for overlap](#probability)
* [Reconciling bull lists](#reconcile)

<a name="organizing"></a>
## Organizing the data

Several of Joel's bulls have been sequenced already and I would like to perform the analysis right away to give him some data to work with. I will make a list of all previously sequenced animals that we have and then run my analysis pipeline when ready.

#### Here is the pipeline:
* SNP and INDEL calling with Samtools
* RD CNV calling with JaRMS
* PE+SR CNV calling with RAPTR-SV
* Variant annotation with my java program and SNPEff

#### Here are the animals and their attributed contributors
* Genome Canada
	* HOLCANM000005279989.bam      
	* HOLDEUM000000253642.bam      
	* HOLUSAM000002265005.bam      
	* HOLUSAM000017349617.bam
	* HOLCANM000006026421.bam      
	* HOLGBRM000000598172.bam      
	* HOLUSAM000002297473.bam      
	* HOLUSAM000123066734.bam
	* HOLCANM000100745543.bam      
	* HOLITAM006001001962.bam      
	* HOLUSAM000017129288.bam      
	* HOLUSAM000132973942.bam
* Our lab
	* HODEU000000253642
	* HOGBR000000598172
	* HOUSA000002040728
	* HOUSA000002290977
	* HOUSA000017349617
	* HOUSA000122358313
* 1000 bulls
	* HOUSA000129800008
	* HOITA006001001962
	* HOUSA000131823833
	* HOUSA000130588960
	* HOUSA000123066734
	* HOUSA000001697572
	* HOCAN000010705608
	* HOUSA000017129288
	* HOAUS000H00930377
	* HOUSA000002297473
	* HOCAN000006820564
	* HOUSA000060372887
	* HOUSA000002103297
	* HOCAN000006947936
	* HOCAN000006026421
	* HOCAN000005279989
	* HOUSA000002205082
	* HODEU000578448776
	* HOUSA000002147486
	* HOUSA000002183007
	* HOCAN000006961162

That makes **39** animals already available for some form of processing. I will explain the caveats about using the 1000 bulls data to Joel later. Here are a list of these animals on the 3850 for future reference.

> 3850: /seq1/bickhart/side_projects/joels_bulls

* **final_file_joel_bulls.list** -> Contains bam files that we currently maintain from Genome Canada and our project
* **1000_bulls_sequenced_joels_bulls.list** -> Contains the names (reformatted) of bulls that were in the 1000 bulls data files but we do not have the bams for.

I am missing 2 bulls from Joel's list. Here are my notes on how to reconcile the lists.

> 3850: /seq1/bickhart/side_projects/joels_bulls

```bash
perl -lane 'print $F[1];' < /work1/grw/grw1/Joel/With_Seq_41.txt > georges_list_41.txt
cat final_file_joel_bulls.list 1000_bulls_sequenced_joels_bulls.list > my_list_bulls.txt

# I used vim to format the my_list_bulls.txt file so that the names were uniform
# Still quite a few were missing
```

I need to go back to the source and use the last list of bulls to reconcile the differences. 

```bash
for i in ../../../1000_bulls_bams/HO*.bam; do name=`basename $i | cut -d'.' -f1`; echo $name; done > 100_bulls_holsteins.list

# ID bulls in 1000 bulls data
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt
mv group_1_2.txt 1000_bulls_joel_with_canada.list
mv group_2.txt bulls_minus_1000_data.list

# I realized that I can generate the association for all files at once
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: canadian_bulls.list
File Number 4: 100_bulls_holsteins.list
Set     Count
1       364
1;2     19
1;2;3   9
1;2;3;4 3
1;2;4   5
1;4     26
2       34
2;4     1
4       2

mv group_1_2.txt 1000_bulls_only_joel.list
cat group_1_2_4.txt group_1_2_3.txt group_1_2_3_4.txt group_2_4.txt > joel_already_controlled.list

# George has a list of bulls that have sequence, supposedly. Let's use that list instead
perl -lane 'print $F[1];' < With_Seq_41.txt > georges_list_sequenced.list
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt georges_list_sequenced.list canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: georges_list_sequenced.list
File Number 4: canadian_bulls.list
File Number 5: 100_bulls_holsteins.list
Set     Count
1       364
1;2;3   19
1;2;3;4 9
1;2;3;4;5       3
1;2;3;5 5
1;5     26
2       30
2;3     4
2;3;5   1
5       2

# OK, let's find out who is in the 2;3 group and where they're located
head group_2_3.txt
HODEU000000254210
HOUSA000002247437  <---- this was the ID that I found later through animal keys
HONLD000839380546
HODEU000000830287

for i in HODEU000000254210 HOUSA000002247437 HONLD000839380546 HODEU000000830287; do getids $i | perl -e '<>; <>; while(<>){chomp; @s = split(/\s+/); print "$s[1]\n";}' > $i.names; done
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl 1000_bulls_sequenced_reformatted.list HODEU000000254210.names
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: HODEU000000254210.names
Set     Count
1       425
1;2     1
2       5
# I did this for the rest, and only found discrepencies with HOUSA000002247437 and HONLD000839380546

# Testing George Liu's list for preferred ids
for i in `cat george_lius_list.txt`; do getids $i | perl -e '<>; <>; while(<>){chomp; @s = split(/\s+/); print "$s[1]\n";}'; done > george_lius_list_altnames.txt
grep HOUSA000002247437 george_lius_list_altnames.txt
grep HONLD000839380546 george_lius_list_altnames.txt

# Nothing
# HONLD000839380546 was in the 1000 bulls list, but I was unable to find it because of my permissive substitution script
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list Bulls_100_sons_1504_ID.txt georges_list_sequenced.list canadian_bulls.list 100_bulls_holsteins.list
File Number 1: 1000_bulls_sequenced_reformatted.list
File Number 2: Bulls_100_sons_1504_ID.txt
File Number 3: georges_list_sequenced.list
File Number 4: canadian_bulls.list
File Number 5: 100_bulls_holsteins.list
Set     Count
1       363
1;2;3   20
1;2;3;4 9
1;2;3;4;5       3
1;2;3;5 5
1;5     26
2       30
2;3     3
2;3;5   1
5       2

# So we remove HOUSA000002247437 but process the rest
grep -v 'HOUSA000002247437' group_2_3.txt > alt_1000_bulls_ids.txt
mv group_1_2_3.txt 1000_bulls_presumptive_list.list
cat group_2_3_5.txt group_1_2_3_5.txt group_1_2_3_4_5.txt group_1_2_3_4.txt > joels_bulls_we_already_have.list

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list HODEU000000254210.names
# This is HODEU000578194407
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o 1000_bulls_sequenced_reformatted.list HODEU000000830287.names
# This is HODEU002261530135
```

The files that contain the proper sequenced bulls are:
* /seq1/bickhart/side_projects/joels_bulls/1000_bulls_presumptive_list.list
* /seq1/bickhart/side_projects/joels_bulls/joels_bulls_we_already_have.list

<a name="onethousand"></a>
## Generating 1000 bulls SNP and INDEL annotations

I am going to create tabular excel output for Joel by splitting the 1000 bull genomes files into the tabular output that I've generated previously for Tad. Joel would like the entire chromosome from the 1kbulls VCF file. While this makes things easier for me, it may take longer to process. I'll get started immediately.

I just rewrote my [vcf subsectioning](https://github.com/njdbickhart/perl_toolchain/blob/master/vcf_utils/filterAndSubsectionVCFfile.pl) script to handle the processing of the vcf.

> 3850: /seq1/1kbulls_annotatedvcf

```bash
bgzip Chr5-Beagle-Run5.eff.vcf
bcftools index Chr5-Beagle-Run5.eff.vcf.gz

perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f Chr5-Beagle-Run5.eff.vcf.gz -o Chr5_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list

```

Here are the output data columns:

<a name="outheads"></a>
#### Output tab file columns

1. Chromosome
2. Base pair position
3. Reference allele
4. Alternate allele
5. Quality score (not present in 1k bulls data)
6. Type (can be either a "SNP" or an "INDEL")
7. Mutation (different classifications based on the mutation's presence in or outside of a gene)
8. Priority (can be "HIGH," "MODERATE," "LOW," or "MODIFIER" in decreasing order of interest)
9. Gene (gene name or ensembl accession)
10. AA or Amino Acid (only present in cases of nonsynonymous mutations)
11. Columns 11 through the last column are animal genotypes

I need to process the following, remaining chromosomes for Joel: 15, 3, 20, 6, 14, 16, 1, 18, 19, 7, 13, 28, 10.

I've created a short bash script to thread this to speed things up.

#### /seq1/1kbulls_annotatedvcf/process_script.sh
```bash
# $1 = chromosome

beaglevcf=${1}-Beagle-Run5.eff.vcf.gz
beagleuncomp=${1}-Beagle-Run5.eff.vcf
progout=${1}_joels_holstein_subsection.tab

# Uncompress and recompress files
gunzip $beaglevcf
bgzip $beagleuncomp
bcftools index $beaglevcf

perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f $beaglevcf -o $progout -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list
```

And here is the code that I'm using to farm out the jobs.

> 3850: /seq1/1kbulls_annotatedvcf/

```bash
for i in Chr15 Chr3 Chr20 Chr6 Chr14 Chr16 Chr1 Chr18 Chr7 Chr13 Chr28 Chr10; do echo $i; sh process_script.sh $i & done

# Chr19 was already bgzipped, so I will process that directly
perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f Chr19-Beagle-Run5.eff.vcf.gz -o Chr19_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_sequenced_joels_bulls_priorformat.list
```

#### I just redid the 1000 bulls list, so I need to reprocess the data so that Joel has the proper animals

```bash
for i in Chr10 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr28 Chr29 Chr6 Chr20 Chr1 Chr3; do echo $i; perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f ${i}-Beagle-Run5.eff.vcf.gz -o ${i}_joels_holstein_subsection.tab -a ../bickhart/side_projects/joels_bulls/1000_bulls_presumptive_list_reformatted.list; done

# All files were copied to Joel's seq1 directory
```

<a name="finertune"></a>
## Calling SNPs and INDELs in Canadian and US Holstein data

I am going to use Samtools to generate SNP and INDEl calls for the remaining data. I will need to discuss the format of CNV calls with him to see what type of input data he can run in his SAS data -- I suspect that problems will be encountered with overlapping variants (ie. CNVs that cover numerous SNPs).

> 3850 /seq1/

```bash
# I just need to identify the bulls that are part of the 100 bulls project that need to be recalled
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1 bickhart/side_projects/joels_bulls/joels_bulls_we_already_have.list bickhart/side_projects/joels_bulls/canadian_bulls.list
File Number 1: bickhart/side_projects/joels_bulls/joels_bulls_we_already_have.list
HOUSA000002290977
HOUSA000002040728
HOUSA000002147486
HOUSA000122358313
HOUSA000002103297
HOUSA000001697572

# The rest are the Canadian bulls
for i in ls /seq1/genome_canada/*.bam; do echo -n "$i "; done; echo
for i in Chr10 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr28 Chr29 Chr6 Chr20 Chr1 Chr3; do samtools mpileup -C50 -gf /seq1/reference/umd_3_1_reference_1000_bull_genomes.fa -uv -t DP -r $i /seq1/genome_canada/HOLCANM000005279989.bam /seq1/genome_canada/HOLCANM000006026421.bam /seq1/genome_canada/HOLCANM000100745543.bam /seq1/genome_canada/HOLDEUM000000253642.bam /seq1/genome_canada/HOLGBRM000000598172.bam /seq1/genome_canada/HOLITAM006001001962.bam /seq1/genome_canada/HOLUSAM000002265005.bam /seq1/genome_canada/HOLUSAM000002297473.bam /seq1/genome_canada/HOLUSAM000017129288.bam /seq1/genome_canada/HOLUSAM000017349617.bam /seq1/genome_canada/HOLUSAM000123066734.bam /seq1/genome_canada/HOLUSAM000132973942.bam /seq1/1000_bulls_bams/HOUSA000002290977.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002040728.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002147486.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000122358313.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002103297.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000001697572.reformatted.sorted.bam | bcftools call -vmO z -o /seq1/bickhart/side_projects/joels_bulls/resequenced_holstein_18_${i}.vcf.gz & done
```

Now, I need to do the fun and exciting part: combining vcf files. I've generated two scripts for this purpose. Here they are:

#### combine_tab_format_vcfs.pl

```perl
#!/usr/bin/perl
# This is a one-shot script designed to combine the genotypes from recalled snps and my
# tab format VCF condensed version

use strict;

my $usage = "perl $0 <input tab file> <input new vcf>\n";
chomp(@ARGV);
unless(scalar(@ARGV) == 2){
        print $usage;
        exit;
}

open(my $VCF, "< $ARGV[1]");
open(my $TAB, "< $ARGV[0]");

# Gather header information
my $tans; my $vans;
my $thead = <$TAB>;
chomp($thead);
my @theadsegs = split(/\t/, $thead);
$tans = scalar(@theadsegs) - 11;

my @vcfheadsegs;
while(my $line = <$VCF>){
        chomp $line;
        if($line =~ /^#CHROM/){
                @vcfheadsegs = split(/\t/, $line);
                $vans = scalar(@vcfheadsegs) - 9;
                last;
        }
}
my @combinedhead;
my $vcflist = scalar(@vcfheadsegs) -1;
push(@combinedhead, (@theadsegs, @vcfheadsegs[9 .. $vcflist]));
print join("\t", @combinedhead) . "\n";
my %vdata; # {chr}->{pos}->{ref} = [alt, qual, animal genotypes]
# Storing vcf file into memory
while(my $line = <$VCF>){
        chomp $line;
        my @vsegs = split(/\t/, $line);
        my $maxlines = scalar(@vsegs) - 1;
        $vdata{$vsegs[0]}->{$vsegs[1]}->{$vsegs[3]} = [$vsegs[4], $vsegs[5], @vsegs[9 .. $maxlines]];
}
close $VCF;

my $notintab = 0;

# Now to sequentially move through both files and process each line
while(my $tline = <$TAB>){
        chomp($tline);
        my @tsegs = split(/\t/, $tline);

        if(exists($vdata{$tsegs[0]}->{$tsegs[1]}->{$tsegs[2]})){
                # Same exact allele
                my $str = generateOutString(\@tsegs, $vdata{$tsegs[0]}->{$tsegs[1]}->{$tsegs[2]}, $tans, $vans, 0);
                push(@{$vdata{$tsegs[0]}->{$tsegs[1]}->{$tsegs[2]}}, -1);
                # the -1 is a signifier that we've used this variant call before
                print "$str\n";
        }else{
                $notintab++;
                my $str = generateOutString(\@tsegs, my $blank, $tans, $vans, 1);
                print "$str\n";
        }
}

close $TAB;
my $notmatch = 0;
foreach my $chr (keys(%vdata)){
        foreach my $pos (keys(%{$vdata{$chr}})){
                foreach my $ref (keys(%{$vdata{$chr}->{$pos}})){
                        my $vref = $vdata{$chr}->{$pos}->{$ref};
                        if($vref->[-1] eq "-1"){
                                next;
                        }else{
                                # We have a VCF entry that didn't match up with the orignal
                                $notmatch++;
                                my @vals;
                                push(@vals, ($chr, $pos, $ref, $vref->[0], $vref->[1], (length($ref) > 1 || length($vref->[0]) > 1)? "INDEL" : "SNP", ".", ".", ".", "."));
                                for(my $x = 0; $x < $tans; $x++){
                                        push(@vals, "./.");
                                }
                                for(my $x = 2; $x < scalar(@{$vref}); $x++){
                                        my $g = $vref->[$x];
                                        my @gsegs = split(":", $g);
                                        push(@vals, $gsegs[0]);
                                }
                                my $str = join("\t", @vals);
                                print "$str\n";
                        }
                }
        }
}

print STDERR "For: $ARGV[0] and $ARGV[1], $notmatch entries in the vcf did not match the tab file, and $notintab for viceversa\n";

exit;

sub generateOutString{
        my ($tref, $vref, $tans, $vans, $visempty) = @_;
        my $str;
        if($visempty){
                my @vals;
                push(@vals, @{$tref});
                for(my $x = 0; $x < $vans; $x++){
                        push(@vals, "./.");
                }
                $str = join("\t", @vals);
        }else{
                my @vals;
                if($vref->[0] ne $tref->[3]){
                        push(@vals, @{$tref});
                        for(my $x = 0; $x < $vans; $x++){
                             push(@vals, "./.");
                        }
                        $str = join("\t", @vals);
                        $str .= "\n";
                        @vals = ();

                        push(@vals, (@{$tref}[0 .. 2], $vref->[0], $vref->[1], (length($tref->[2]) > 1 || length($vref->[0]) > 1)? "INDEL" : "SNP", ".", ".", ".", "."));
                        for(my $x = 0; $x < $tans; $x++){
                                push(@vals, "./.");
                        }
                        for(my $x = 2; $x < scalar(@{$vref}); $x++){
                                my $g = $vref->[$x];
                                my @gsegs = split(":", $g);
                                push(@vals, $gsegs[0]);
                        }
                        $str .= join("\t", @vals);
                }else{
                        push(@vals, @{$tref});
                        for(my $x = 2; $x < scalar(@{$vref}); $x++){
                                my $g = $vref->[$x];
                                my @gsegs = split(":", $g);
                                push(@vals, $gsegs[0]);
                        }
                        $str = join("\t", @vals);
                }
       	}
        return $str;
}

```

And a sorting script:

#### sort_vcf_tab_format.pl

```perl
#!/usr/bin/perl
# This is another one shot script designed to sort my tab vcf output

use strict;
chomp(@ARGV);
my $usage = "perl $0 <input tab file>\n";

unless(scalar(@ARGV) == 1){
        print $usage;
        exit;
}

open(my $IN, "< $ARGV[0]");
my $head = <$IN>;
print $head;
my %data;
while(my $line = <$IN>){
        chomp $line;
        my @segs = split(/\t/, $line);
        push(@{$data{$segs[0]}->{$segs[1]}}, \@segs);
}
close $IN;

foreach my $chr (sort {$a cmp $b} keys(%data)){
        foreach my $pos (sort {$a <=> $b} keys(%{$data{$chr}})){
                foreach my $ref (@{$data{$chr}->{$pos}}){
                        print join("\t", @{$ref}) . "\n";
                }
        }
}
exit;
```


> 3850 /seq1/bickhart/side_projects/joel

```bash
for i in Chr10 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr1 Chr20 Chr28 Chr29 Chr3 Chr6; do perl combine_tab_format_vcfs.pl ../../../joel/${i}_joels_holstein_subsection.tab resequenced_holstein_18_${i}.vcf > joel_combined_holstein_${i}.tab; perl sort_vcf_tab_format.pl joel_combined_holstein_{$i}.tab > joel_combined_holstein_{$i}.sorted.tab; done
For: ../../../joel/Chr10_joels_holstein_subsection.tab and resequenced_holstein_18_Chr10.vcf, 133990 entries in the vcf did not match the tab file, and 1070521 for viceversa
For: ../../../joel/Chr13_joels_holstein_subsection.tab and resequenced_holstein_18_Chr13.vcf, 99519 entries in the vcf did not match the tab file, and 879606 for viceversa
For: ../../../joel/Chr14_joels_holstein_subsection.tab and resequenced_holstein_18_Chr14.vcf, 116600 entries in the vcf did not match the tab file, and 871460 for viceversa
For: ../../../joel/Chr15_joels_holstein_subsection.tab and resequenced_holstein_18_Chr15.vcf, 126369 entries in the vcf did not match the tab file, and 977319 for viceversa
For: ../../../joel/Chr16_joels_holstein_subsection.tab and resequenced_holstein_18_Chr16.vcf, 106784 entries in the vcf did not match the tab file, and 937587 for viceversa
For: ../../../joel/Chr17_joels_holstein_subsection.tab and resequenced_holstein_18_Chr17.vcf, 111746 entries in the vcf did not match the tab file, and 780052 for viceversa
For: ../../../joel/Chr18_joels_holstein_subsection.tab and resequenced_holstein_18_Chr18.vcf, 92108 entries in the vcf did not match the tab file, and 649237 for viceversa
For: ../../../joel/Chr19_joels_holstein_subsection.tab and resequenced_holstein_18_Chr19.vcf, 80263 entries in the vcf did not match the tab file, and 644141 for viceversa
For: ../../../joel/Chr1_joels_holstein_subsection.tab and resequenced_holstein_18_Chr1.vcf, 212633 entries in the vcf did not match the tab file, and 1674948 for viceversa
For: ../../../joel/Chr20_joels_holstein_subsection.tab and resequenced_holstein_18_Chr20.vcf, 95802 entries in the vcf did not match the tab file, and 779795 for viceversa
For: ../../../joel/Chr28_joels_holstein_subsection.tab and resequenced_holstein_18_Chr28.vcf, 66233 entries in the vcf did not match the tab file, and 531242 for viceversa
For: ../../../joel/Chr29_joels_holstein_subsection.tab and resequenced_holstein_18_Chr29.vcf, 81369 entries in the vcf did not match the tab file, and 611619 for viceversa
For: ../../../joel/Chr3_joels_holstein_subsection.tab and resequenced_holstein_18_Chr3.vcf, 151845 entries in the vcf did not match the tab file, and 1198002 for viceversa
For: ../../../joel/Chr6_joels_holstein_subsection.tab and resequenced_holstein_18_Chr6.vcf, 166544 entries in the vcf did not match the tab file, and 1191275 for viceversa
```

Hmm... that's far too much of a sacrifice here. Most of the variant sites aren't matching up.

Let's do a quick test where I extract the same animals and attempt to combine them in the same process.

```bash
for i in /seq1/genome_canada/*.bam; do name=`basename $i | cut -d'.' -f1`; echo "$name"; done > test_set_bulls.list
# I added the 100 bulls entries as well to this list

perl ~/perl_toolchain/vcf_utils/filterAndSubsectionVCFfile.pl -f /seq1/1kbulls_annotatedvcf/Chr28-Beagle-Run5.eff.vcf.gz -o test_Chr28_joels_holstein_subsection.tab -a test_set_bulls.list
perl combine_tab_format_vcfs.pl test_Chr28_joels_holstein_subsection.tab resequenced_holstein_18_Chr28.vcf > test_joel_combined_holstein_Chr28.tab
	For: test_Chr28_joels_holstein_subsection.tab and resequenced_holstein_18_Chr28.vcf, 66233 entries in the vcf did not match the tab file, and 531242 for viceversa

wc -l test_joel_combined_holstein_Chr28.tab
	839441 test_joel_combined_holstein_Chr28.tab
wc -l resequenced_holstein_18_Chr28.vcf
	311287 resequenced_holstein_18_Chr28.vcf
# Testing how many are not monomorphic reference
perl -lane 'if($F[0] eq "CHR"){next;}else{$c = 0; for($x = 10; $x < 22; $x++){if($F[$x] eq "0|0"){$c++;}} if($c >= 12){next;}else{print $_;}}' < test_Chr28_joels_holstein_subsection.tab | wc -l
	769,182

perl -lane 'if($F[0] eq "CHR"){print $_;}else{$c = 0; for($x = 10; $x < 22; $x++){if($F[$x] eq "0|0"){$c++;}} if($c >= 12){next;}else{print $_;}}' < test_Chr28_joels_holstein_subsection.tab > test_Chr28_joels_filtered_holstein_subsection.tab
perl combine_tab_format_vcfs.pl test_Chr28_joels_filtered_holstein_subsection.tab resequenced_holstein_18_Chr28.vcf > test_joel_filtered_holstein_Chr28.tab
	For: test_Chr28_joels_filtered_holstein_subsection.tab and resequenced_holstein_18_Chr28.vcf, 66326 entries in the vcf did not match the tab file, and 527654 for viceversa
```

OK, so let's calculate the venn stats for this. 66,233 unique for the sequence data, 531,242 unique for the 1000 bulls, and 308,199 shared sites. ~ 4,000 predicted monomorphic reference sites in the 1k bulls data are variant calls within the sequence data.

<a name="probability"></a>
## Plan B: generating variant site probabilities for overlap

OK, I need to revise my plan to account for the discrepancy of overlap here. I'm going to generate genotype probability scores and then use them to fill out the VCF tab file.

```bash
for i in Chr10 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr28 Chr29 Chr6 Chr20 Chr1 Chr3; do samtools mpileup -C50 -gf /seq1/reference/umd_3_1_reference_1000_bull_genomes.fa -uv -t DP -r $i /seq1/genome_canada/HOLCANM000005279989.bam /seq1/genome_canada/HOLCANM000006026421.bam /seq1/genome_canada/HOLCANM000100745543.bam /seq1/genome_canada/HOLDEUM000000253642.bam /seq1/genome_canada/HOLGBRM000000598172.bam /seq1/genome_canada/HOLITAM006001001962.bam /seq1/genome_canada/HOLUSAM000002265005.bam /seq1/genome_canada/HOLUSAM000002297473.bam /seq1/genome_canada/HOLUSAM000017129288.bam /seq1/genome_canada/HOLUSAM000017349617.bam /seq1/genome_canada/HOLUSAM000123066734.bam /seq1/genome_canada/HOLUSAM000132973942.bam /seq1/1000_bulls_bams/HOUSA000002290977.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002040728.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002147486.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000122358313.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002103297.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000001697572.reformatted.sorted.bam

samtools mpileup -C50 -gf /seq1/reference/umd_3_1_reference_1000_bull_genomes.fa -uv -t DP -r Chr5 /seq1/genome_canada/HOLCANM000005279989.bam /seq1/genome_canada/HOLCANM000006026421.bam /seq1/genome_canada/HOLCANM000100745543.bam /seq1/genome_canada/HOLDEUM000000253642.bam /seq1/genome_canada/HOLGBRM000000598172.bam /seq1/genome_canada/HOLITAM006001001962.bam /seq1/genome_canada/HOLUSAM000002265005.bam /seq1/genome_canada/HOLUSAM000002297473.bam /seq1/genome_canada/HOLUSAM000017129288.bam /seq1/genome_canada/HOLUSAM000017349617.bam /seq1/genome_canada/HOLUSAM000123066734.bam /seq1/genome_canada/HOLUSAM000132973942.bam /seq1/1000_bulls_bams/HOUSA000002290977.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002040728.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002147486.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000122358313.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000002103297.reformatted.sorted.bam /seq1/1000_bulls_bams/HOUSA000001697572.reformatted.sorted.bam | bcftools call -vmO z -o resequenced_holstein_18_Chr5.vcf.gz
```

<a name="reconcile"></a>
## Reconciling bull lists

There were more problems with the list of animals that have been sequenced. I'm going to try to resolve these by starting from scratch with George's list, the 1000 bulls list and Joel's new list.

First of all, let's stop with the comparison of international IDs as they are not unique identifiers. I will be jury rigging the 

> 3850: /seq1/bickhart/side_projects/joels_bulls/reformat_lists

```bash
# formatting Joel's email list
perl -lane 'print $F[1];' < raw_data > pauls_list_wgs_bulls.list
perl -lane 'open(IN, "echo $F[0] | ID2key |"); $l = <IN>; @s = split(/\s+/, $l); print "$s[1]"; close IN;' < pauls_list_wgs_bulls.list > pauls_list_wgs_bulls.ids

# copying datasets for venn comparisons
perl -lane 'open(IN, "echo $F[0] | ID2key |"); $l = <IN>; @s = split(/\s+/, $l); print "$s[1]"; close IN;' < ../1000_bulls_sequenced.list > 1000_list_wgs_bulls.ids
perl -lane 'open(IN, "echo $F[0] | ID2key |"); $l = <IN>; @s = split(/\s+/, $l); print "$s[1]"; close IN;' < ../Bulls_100_sons_1504_ID.list > georges_list_wgs_bulls.ids
perl -lane 'open(IN, "echo $F[0] | ID2key |"); $l = <IN>; @s = split(/\s+/, $l); print "$s[1]"; close IN;' < ../canadian_bulls.list > canadian_list_wgs_bulls.ids
perl -lane 'open(IN, "echo $F[0] | ID2key |"); $l = <IN>; @s = split(/\s+/, $l); print "$s[1]"; close IN;' < ../100_bulls_holsteins.list > liu_list_wgs_bulls.ids

perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl pauls_list_wgs_bulls.ids 1000_list_wgs_bulls.ids georges_list_wgs_bulls.ids canadian_list_wgs_bulls.ids liu_list_wgs_bulls.ids prior_list_wgs_bulls.ids
	File Number 1: pauls_list_wgs_bulls.ids
	File Number 2: 1000_list_wgs_bulls.ids
	File Number 3: georges_list_wgs_bulls.ids
	File Number 4: canadian_list_wgs_bulls.ids
	File Number 5: liu_list_wgs_bulls.ids
	File Number 6: prior_list_wgs_bulls.ids
	Set     Count
	1;2     1
	1;2;3;4;5;6     3
	1;2;3;4;6       9
	1;2;3;5;6       5
	1;2;3;6 20
	2       333
	2;3     1
	2;3;6   2
	2;5     26
	3       30
	3;5;6   1
	5       2

# Checking individual bulls now to see which ones are different among the lists
# Paul's list vs George's list
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1_2 pauls_list_wgs_bulls.ids 1000_list_wgs_bulls.ids georges_list_wgs_bulls.ids canadian_list_wgs_bulls.ids liu_list_wgs_bulls.ids
	File Number 1: pauls_list_wgs_bulls.ids
	File Number 2: 1000_list_wgs_bulls.ids
	51356830

echo "51356830" | key2ID2
getids HOUSA000137974489
	51356830 HOUSA000137974489 M    B   1  ** 2016-05-24 RONELEE TOYSTORY DOMAIN-ET

# George's list and the 1000 bulls, but not in my prior list of sequenced animals
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 2_3 pauls_l                          t_wgs_bulls.ids georges_list_wgs_bulls.ids canadian_list_wgs_bulls.ids liu_list_wgs_bulls.ids prior_list_wgs_bulls.ids
	File Number 2: 1000_list_wgs_bulls.ids
	File Number 3: georges_list_wgs_bulls.ids
	2539707

echo "2539707" | key2ID2
getids HOUSA000002247437
	2539707 HOUSA000002247437 M    B   1  ** 1999-08-25 528 ETAZON CELSIUS-ET


# Now converting to a tab delimited list of all bull sequenced (and by whom!)
perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o pauls_list_wgs_bulls.ids 1000_list_wgs_bulls.ids georges_list_wgs_bulls.ids canadian_list_wgs_bulls.ids liu_list_wgs_bulls.ids prior_list_wgs_bulls.ids

perl -lane '$type = "1000;Canada;Liu"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_1_2_3_4_5_6.txt > current_list_of_animals.tab
perl -lane '$type = "1000"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_1_2.txt >> current_list_of_animals.tab
perl -lane '$type = "1000;Canada"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_1_2_3_4_6.txt >> current_list_of_animals.tab
perl -lane '$type = "1000;Liu"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_1_2_3_5_6.txt >> current_list_of_animals.tab
perl -lane '$type = "1000"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_1_2_3_6.txt >> current_list_of_animals.tab
perl -lane '$type = "1000"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_2_3.txt >> current_list_of_animals.tab
perl -lane '$type = "1000"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_2_3_6.txt >> current_list_of_animals.tab
perl -lane '$type = "Liu"; open(IN, "echo $F[0] | key2ID2 |"); $l = <IN>; $l =~ s/^\s+//; @s = split(/\s+/, $l); print "$s[0]\t$s[1]\t$type"; close IN;' < group_3_5_6.txt >> current_list_of_animals.tab
```