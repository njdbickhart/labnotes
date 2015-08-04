# Placing SNPs using flanking sequence
---
*8/4/2015*

This is a protocol on how to generate genomic positions for SNP probes based on their flanking sequence.

The following required software must be on the user's path:
* BWA
* Perl
* Perl Modules:
	* Mouse
	* namespace::autoclean
* [My github repository for perl scripts](https://github.com/njdbickhart/perl_toolchain)

## Isolating flanking sequence.

This is highly dependent on the input file, but most file formats used by SNP chip designers flank the variant nucleotide with square brackets. Since the flanking sequence is normally 100bp+, I can use this in a short read alignment tool to generate coordinates fast.

George gave me a list of probes that need UMD3.1 reference genome assignment, so I'm going to pull their flanking sequence and generate a fasta file. 

Here is the file that contains the flanking probe sequence:
**/work1/grw/chips/GH2/GGP_HDv2_B_StrandReport_FDT_V1.csv**

> 3850: /home/bickhart

```bash
# First, I'm going to transfer the names of the probes to a text file
vim probenames.txt

# Now to read them in using Perl in order to extract the nucleotide sequence
perl -e 'chomp(@ARGV); my %h; open(IN, "< $ARGV[0]"); while(<IN>){chomp; $h{$_} = 1;}; close IN; open(IN, "< $ARGV[1]"); while(<IN>){chomp; @s = split(/,/); if(exists($h{$s[1]})){ @b = split(/[\[\]]/, $s[6]); $b[2] =~ tr/ACGT/TGCA/; $b[2] = reverse($b[2]); print ">$s[1].f\n$b[0]\n>$s[1].r\n$b[2]\n";}} close IN;' probenames.txt /work1/grw/chips/GH2/GGP_HDv2_B_StrandReport_FDT_V1.csv > probenames.fa

# Each flanking region is defined by ".f" and ".r" with ".f" sequence being the 5' flanking and ".r" as the 3' flanking
# ".r" sequence is reverse complemented so that the 3' end of both sequence regions should terminate with the variant site
```

## Alignment and probe positioning

Next, I'm going to align the data to the UMD3.1 reference genome and then use a script to parse the position

> 3850: /home/bickhart

```bash
# I'm going to use a previously formatted version of the UMD3.1 reference genome here
bwa mem /seq1/reference/umd3_kary_unmask_ngap.fa probenames.fa > probenames.sam

# Now that I have the sam file, I'm going to use a script I wrote to check the coordinate positions of the SNPs
perl perl_toolchain/snp_utilities/getSNPPositionFromSAM.pl probenames.sam probenames.tab

# let's check how it did:
perl -lane 'if($F[3] == 1){print $_; }' < probenames.tab | wc -l
	17
wc -l probenames.tab
	69 probenames.tab

# So 17 out of 68 had probe sequence that matched precisely
perl -lane 'if(abs($F[2] - $F[7]) == 1){print $_;}' < probenames.tab | wc -l
	38

# 38 of those differed by 1 base, so the forward direction base is likely correct.
```

All data was saved to **probenames.tab**