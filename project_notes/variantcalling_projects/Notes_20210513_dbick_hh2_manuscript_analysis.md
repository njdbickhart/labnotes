# HH2 manuscript analysis
---
*5/13/2021*

## Concordance analysis

| Locus | Coords | Description|
| :--- | :--- | :--- |
| IFT80 | chr1:107,087,310-107,210,812 | ARS-UCD1.2 coordinates of gene |
| IFT80 extend | chr1:107,056,435-107,241,687 | UCSC 50% coordinate expand |
| HH2 variant site | chr1:107172615| Actual variant site |

> Ceres: /lustre/project/bostauruscnv/bam

```bash
sbatch -N 1 -n 2 -p priority -q msn --mem=55000 --wrap="python3 ~/python_toolchain/sequenceData/vcfFileConcordanceAnalysis.py -v /lustre/project/bostauruscnv/bam/vcf/allseq1_annotated.vcf.gz -c 1 -s 107056435 -e 107241687 -l carriers.list -o ift80_concordance_analysis -m 0.93"

# Checking concordance at the exact site:
gunzip -c /lustre/project/bostauruscnv/bam/vcf/allseq1_annotated.vcf.gz | perl -e '@header; while(<>){chomp; if($_ =~ /##/){next;}elsif($_ =~ /#CHROM/){ @header = split(/\t/, $_);}else{@s = split(/\t/); if($s[0] eq "1" && $s[1] eq "107172615"){for($x = 9; $x < scalar(@s); $x++){@g = split(/[\/:;]/, $s[$x]); if(($g[0] ne "0" && $g[0] ne ".") || ($g[1] ne "0" && $g[1] ne ".")){print "$header[$x]\n";}}}}}' | head -n 15
HOLCANM000006193092
HOLUSAM000063685691
HOLUSAM000064178187
HOL840M003014467478
HOLUSAM000061365427
HOLUSAM000069983398
HOLUSAM000071813244
HOLUSAM000137033416

# there's one carrier that doesn't have a variant site at the location
perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl carriers.list divergents.list
File Number 1: carriers.list
File Number 2: divergents.list
Set     Count
1       1
1;2     8

perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -l 1 carriers.list divergents.list
File Number 1: carriers.list
HOLUSAM000061547476

gunzip -c /lustre/project/bostauruscnv/bam/vcf/allseq1_annotated.vcf.gz | grep '#CHROM' | perl -lane 'for($x = 0; $x < scalar(@F); $x++){print "$x\t$F[$x]";}' | grep 'HOLUSAM000061547476'
51      HOLUSAM000061547476

gunzip -c /lustre/project/bostauruscnv/bam/vcf/allseq1_annotated.vcf.gz | perl -e '@header; while(<>){chomp; if($_ =~ /##/){next;}elsif($_ =~ /#CHROM/){ @header = split(/\t/, $_);}else{@s = split(/\t/); if($s[0] eq "1" && $s[1] eq "107172615"){print "$s[51]\n";}}}' 

gunzip -c /lustre/project/bostauruscnv/bam/vcf/allseq1_annotated.vcf.gz | perl -e '@header; while(<>){chomp; if($_ =~ /##/){next;}elsif($_ =~ /#CHROM/){ @header = split(/\t/, $_);}else{@s = split(/\t/); if($s[0] eq "1" && $s[1] eq "107172615"){print "$s[51]\n";}}}'
0/0:14,0,120:25:21,4

# GT format: GT:PL:DP:AD
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

```