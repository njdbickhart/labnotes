# IGC variant selection and filtering
---
**9/5/2017**

These are my notes on the selection and filtering of sequence variants suitable for genotyping large populations of cattle. The overall goal is to select 50 SNP targets for downstream genotyping and then branch out for there.

## Table of Contents
* [Selecting regions for SNP analysis](#selecting)
* [Generating VCFs for genomic region subsets](#generating)
* [Including the last section of chr5 NKC](#last_sec)
* [Tabulating and interpreting Neogen results](#neogen)
* [Placing variant sites](#placing)

<a name="selecting"></a>
## Selecting regions for SNP analysis

One of the problems with this project is that I have two haplotypes that represent each respective IGC on the reference. I need to identify which regions of the reference correspond to the IGCs first. The good news is that I have generated liftover files to help with this! Here are the locations of the IGC regions on the UMD3.1 reference:

> Assembler 2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing

```bash
cat igc_regions.bed
	chr23   28300000        28750000        MHC
	chr18   63100000        63400000        LRC
	chr5    99500000        99850000        NKC
```

And now to liftover.

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver igc_regions.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/combined.umd3_to_ars-ucd.liftover.chain igc_regions.ars-ucd.v14.bed igc_regions.ars-ucd.v14.unmapped

wc -l igc_regions.ars-ucd.v14*
  0 igc_regions.ars-ucd.v14.bed
  6 igc_regions.ars-ucd.v14.unmapped
  6 total

# Allowing multiple output
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver -multiple igc_regions.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/combined.umd3_to_ars-ucd.liftover.chain igc_regions.ars-ucd.v14.bed igc_regions.ars-ucd.v14.unmapped

wc -l igc_regions.ars-ucd.v14*
 15 igc_regions.ars-ucd.v14.bed
  0 igc_regions.ars-ucd.v14.unmapped
 15 total
```

Well, wasn't going to be easy from the start! Looks like the NKC didn't map at all here, but it's probably an issue with the UMD3.1 reference. I may need to search for the haplotype myself through alignment of one of our contigs.

```bash
module load bwa
# Using a previously uploaded list of IGC haplotypes I have laying around
bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta /mnt/nfs/nfs2/dbickhart/igc/usda_cumulative_igc_haplotypes.fa > /mnt/nfs/nfs2/dbickhart/igc/usda_cumulative_igc_haplotypes.sam

perl ~/sperl/sequence_data_scripts/BriefSamOutFormat.pl -s /mnt/nfs/nfs2/dbickhart/igc/usda_cumulative_igc_haplotypes.sam | grep 'NKC'
LIB14370_NKC    16      5       98898487        99081393        219086  60
LIB14398_NKC    0       5       99382031        99508055        126020  60
LIB14414_NKC    16      5       99210422        99254870        107677  60
LIB14435_NKC    16      5       99365642        99504627        138986  60
LIB14363_NKC    0       5       98792567        98900149        205555  60
```

So let's manually summarize the data using the mapping coordinates generated:

#### ARS-UCDv14 base IGC coordinates:

| IGC | Chr | Start | End |
|:--- | :---| ---: | ---: |
|MHC | 23 | 28460024 | 28533695 |
|LRC | 18 | 63110866 | 63144659 |
|NKC | 5 | 98792567 | 99508055 |

So now we can generate a coordinate file for samtools mpileup so that we're only processing the regions we need to investigate.

```bash
# Manually edited region file
cat igc_ars_ucd_v14_regions.tab
	CH240_391K10_KIR        LRC1
	Domino_MHCclassI_gene2-5hap1_MHC        MHC1
	CH240_370M3_LILR_LRC    LRC2
	TPI4222_A14_MHCclassI_MHC       MHC2
	HF_LRC_hap1_KIR_LRC     LRC3
	LIB14427_MHC    MHC3
	LIB14413_LRC    LRC4
	23:28460024-28533695    MHCA
	18:63110866-63144659    LRCA
	5:98792567-99508055     NKCA
```

<a name="generating"></a>
## Generating VCFs for genomic region subsets

From the CDDR's natdb project, I have about 127 bams consisting of animals of diverse backgrounds in the Holstein breed. This will likely be the most diverse genetic background I can throw at this project, so this is our best shot at checking the frequency of our assembled haplotypes against the normal reference genome.

Let's start by preparing the BAMs for mpileup.

> Assembler 2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing

```bash
ls /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/aligns/*/*.merged.bam > igc_variant_list_bams.list

module load samtools

cat igc_variant_list_bams.list | xargs -I {} sbatch --ntasks-per-node=1 --mem=2000 --nodes=1 --wrap="samtools index {}"
```

Now that the BAMs are getting indexed, let's start writing the wrapper for mpileup. The script name is [samtoolsSelectiveMpileup.pl](https://github.com/njdbickhart/perl_toolchain/blob/master/sequence_data_scripts/samtoolsSelectiveMpileup.pl). 

And, finally, to run the script.

```bash
perl -lane 'system("sbatch /home/dbickhart/sperl/sequence_data_scripts/samtoolsSelectiveMpileup.pl -b igc_variant_list_bams.list -s $F[0] -n $F[1] -f /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta");' < igc_ars_ucd_v14_regions.tab
```

The script ran successfully, but now I need to run the variant calls through bcftools to finalize the calls.

```bash
module load samtools
module load bcftools

for i in *.bcf; do bcftools index $i; done
bcftools concat -o igc_concat.mpileup.bcf -O b LRC1.mpileup.bcf LRC2.mpileup.bcf LRC3.mpileup.bcf LRC4.mpileup.bcf LRCA.mpileup.bcf MHC1.mpileup.bcf MHC2.mpileup.bcf MHC3.mpileup.bcf MHCA.mpileup.bcf NKCA.mpileup.bcf
bcftools call -vmO z -o igc_125_animals.calls.vcf.gz igc_concat.mpileup.bcf

gunzip -c igc_125_animals.calls.vcf.gz | wc -l
bcftools stats -F /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta -s - igc_125_animals.calls.vcf.gz > igc_125_animals.calls.vcf.gz.stats

# Generating plots of variant call stats
mkdir plots
plot-vcfstats -p plots/ igc_125_animals.calls.vcf.gz.stats

# Processing each region separately
bcftools concat -o LRC.combined.bcf -O b LRC1.mpileup.bcf LRC2.mpileup.bcf LRC3.mpileup.bcf LRC4.mpileup.bcf LRCA.mpileup.bcf
bcftools concat -o MHC.combined.bcf -O b MHC1.mpileup.bcf MHC2.mpileup.bcf MHC3.mpileup.bcf MHCA.mpileup.bcf
```

Some interesting features that are derived from the stats of this dataset:
* Some animals have near 2.0 TS/TV ratios, suggesting that their alleles are likely well covered by the current assembly
* HWE is not observed in the distribution of alleles in the population, with a bimodal distribution of allele frequencies
* The majority of variants have a quality score of 998+, suggesting that larger sample sizes have given this dataset better quality

Generating a plot of variants per chromosome per animal

```R
data <- read.delim("igc_125_animals.calls.ids.calls.peranimal.tab", header=FALSE)
data.formated <- data[c("X18", "X23", "X5", "CH240_370M3_LILR_LRC", "CH240_391K10_KIR", "Domino_MHCclassI_gene2.5hap1_MHC", "HF_LRC_hap1_KIR_LRC", "LIB14413_LRC", "LIB14427_MHC", "TPI4222_A14_MHCclassI_MHC")]
rownames(data.formated) <- data$X

library(ggplot2)
library(reshape2)

data.formated.c <- data.frame(data.formated, count = c(1:125))
data.melt <- melt(data.formated.c, id.vars="count", measure.vars = colnames(data.formated))
ggplot(data.melt, aes(count, value, colour=variable)) + geom_point() + xlab(label="Holstein Bull Samples") + ylab(label="Variant counts per chromosome") + ggtitle(label="Variant counts per chromosome per animal")
dev.copy2pdf(file="per_sample_variant_counts.pdf", useDingbats=FALSE)

```

I need to assign variant sites to different haplotypes within each subgroup. I suspect that the best way to do this will be through hierarchical clustering, as I do not know the number of haplotypes within my samples (I don't trust the AIP haplotypes in these regions). I will attempt to generate clusters within each subgroup and then assign variant sites to each haplogroup. I will then filter variants based on aligned read MAPQ scores directly adjacent to the site, within 36bp of each other.

I think that I can use structure to assign clusters and then use that to try to backtrace marker assignments.

```bash
# first, to generate a suitable structure file for each vcf
module load plink/2.00alM-2017-05-22
# Running on the whole dataset to see how many variants are in LD
plink2 --vcf igc_125_animals.calls.vcf.gz --indep-pairwise 50 5 0.5 --allow-extra-chr

# OOPs! Each SNP needs an ID for this analysis!
bcftools annotate -o igc_125_animals.calls.ids.vcf.gz -O z -I 'ARS\_PIRBRIGHT\_%CHROM\_%POS' igc_125_animals.calls.vcf.gz

plink2 --vcf igc_125_animals.calls.ids.vcf.gz --indep-pairwise 50 5 0.5 --allow-extra-chr
wc -l plink2.prune.in plink2.prune.out
  20891 plink2.prune.in
  38937 plink2.prune.out
  59828 total

# That's a good first filter, let's now try to run some PLINK clustering to identify individual clusters
# Great, plink2 doesn't provide clustering algorithms!
# Thankfully, Plink1.90 does.

plink --vcf igc_125_animals.calls.ids.vcf.gz --extract plink2.prune.in --cluster --allow-extra-chr --threads 20 --double-id
# Only one "cluster" present... a little dissappointing. 
# Trying to calculate IBD
plink --vcf igc_125_animals.calls.ids.vcf.gz --extract plink2.prune.in --genome --allow-extra-chr --threads 20 --double-id

# Because of the small number of markers I have and the bias in marker selection (from the same polymorphic regions), the IBD estimates are pretty poor indicators of actual descent
# One last filter: HWE
plink --vcf igc_125_animals.calls.ids.vcf.gz --extract plink2.prune.in --hardy --allow-extra-chr --threads 20 --double-id

perl -lane 'if($F[8] < 0.05){print $_;}' < plink.hwe | wc -l
	5904 <- Let's test to see how many per "chromosome"

perl -lane 'if($F[8] < 0.05){print $F[0];}' < plink.hwe | perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -c 0 -f stdin -m
```
#### Genetics filtered variant list

|Entry                            | Count|
|:--------------------------------|-----:|
|18                               |    47|
|23                               |   220|
|CH240_370M3_LILR_LRC             |   182|
|CH240_391K10_KIR                 |   109|
|Domino_MHCclassI_gene2-5hap1_MHC |   760|
|HF_LRC_hap1_KIR_LRC              |   527|
|LIB14413_LRC                     |   186|
|LIB14427_MHC                     |   647|
|TPI4222_A14_MHCclassI_MHC        |  1328|

OK, now I need to filter these by virtue of their flanking sequence map quality scores. The scheme will be to test the upstream and downstream 36 bp to assess the average MAPQ of the aligned reads. If the MAPQ is low (< 12 =~ 5% chance mapping site is incorrect) then I'll drop the read. The pileup takes several seconds to load, so it is infeasible to subset each time for each marker position. I'm going to generate a large flatfile and use a second script to process it.

```bash
perl -lane 'print "$F[0]"; system("samtools mpileup -s -O -r $F[0] -b igc_variant_list_bams.list >> igc_variant_pileup_regions.tab");' < igc_ars_ucd_v14_regions.tab

# I need to remove SNP IDs that are not unique and HWEP violating SNPs first
perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f plink2.prune.in -c 0 | perl -lane 'if($F[1] > 1){print $F[0];}' > plink2.prune.in.dups
perl -lane 'if($F[0] eq "CHR"){next;} if($F[8] > 0.05){print $F[1];}' < plink.hwe > plink.hwe.remove

perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl plink2.prune.in plink2.prune.in.dups plink.hwe.remove
File Number 1: plink2.prune.in
File Number 2: plink2.prune.in.dups
File Number 3: plink.hwe.remove
Set     Count
1       5804
1;2     5
1;2;3   51
1;3     14975

perl ~/sperl/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o plink2.prune.in plink2.prune.in.dups plink.hwe.remove
perl ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v igc_125_animals.calls.ids.vcf.gz -l group_1.txt -m igc_variant_pileup_regions.tab -o igc_mapq_variant_sites.tab

perl -d ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v igc_125_animals.calls.ids.vcf.gz -l group_1.txt -m igc_variant_pileup_regions.tab -o igc_mapq_variant_sites.tab

# my script was too cumbersome and pushed the limits of Perls' processing limits. 
# I'm going to subset the data and work on it using other reliable tools
perl -ne '@F = split(/\t/); $e = $F[1] + 1; $sum = 0; $c = 0; for($x = 7; $x < scalar(@F);$x += 5){if($F[$x] ne "*"){@bsegs = split(/,/, $F[$x]); $c += scalar(@bsegs); foreach $j (@bsegs){$sum += $j;}}} $avg = ($c > 0)? $sum / $c : 0; $sum = 0; $c = 0; print "$F[0]\t$F[1]\t$e\t$avg\n";' < igc_variant_pileup_regions.tab > igc_variant_pileup_regions.scores.bed

# Now to select only the sites that we kept after HWE and LD pruning
perl -lane 'if($F[1] eq "SNP"){next;}elsif($F[8] < 0.05){print $F[1]}' < plink.hwe > igc_125_animals.calls.ids.calls.passfilter.snpids
perl -e 'chomp(@ARGV); open($IN, "< $ARGV[0]"); %snps; while(<$IN>){chomp; $snps{$_} = 1;} close $IN; open($IN, " gunzip -c $ARGV[1] |"); while(<$IN>){chomp; if($_ =~ /^#/){next;}else{@s = split(/\t/); if(exists($snps{$s[2]})){$e = $s[1] + 1; print "$s[0]\t$s[1]\t$e\t$s[2]\n";}}} close $IN;' igc_125_animals.calls.ids.calls.passfilter.snpids igc_125_animals.calls.ids.vcf.gz > igc_125_animals.calls.ids.calls.passfilter.bed

# I think that there are some duplicate IDs -- checking really quickly
cat igc_125_animals.calls.ids.calls.passfilter.bed | cut -f4 | sort | uniq | wc -l
5882
cat igc_125_animals.calls.ids.calls.passfilter.snpids | sort | uniq | wc -l
5882

# There are -- most likely due to multi-allelic sites. Let's find the duplicates and remove them
perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f igc_125_animals.calls.ids.calls.passfilter.bed -c 3 | perl -lane 'if($F[1] > 1){print $F[0];}' > test_filtered.dupids
perl -e '%exclude; chomp(@ARGV); open($IN, "< $ARGV[0]"); while(<$IN>){chomp; $exclude{$_} = 1;} close $IN; open($IN, "< $ARGV[1]"); while(<$IN>){chomp; @s = split(/\t/); unless(exists($exclude{$s[3]})){print "$_\n";}}' test_filtered.dupids igc_125_animals.calls.ids.calls.passfilter.bed > igc_125_animals.calls.ids.calls.passfilter.nodups.bed

# Now for the association of Mapq sites with the variants
perl ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v igc_125_animals.calls.ids.calls.passfilter.nodups.bed -l igc_125_animals.calls.ids.calls.passfilter.snpids -m igc_variant_pileup_regions.scores.bed -o igc_125_animals.calls.ids.calls.passfilter.nodups.scores.bed

# These are clearly higher quality markers in the best mapping locations
perl -lane 'if($F[4] > 80 && $F[5] > 80){print $_;}' < igc_125_animals.calls.ids.calls.passfilter.nodups.scores.bed | wc -l
1150

# I will now take each of these markers, extract flanking sequence and map them back to the reference genome to assess mapping site quality scores.
# If a 36 mer can't map reliably back to the reference, then we're in trouble!
perl -lane 'if($F[4] > 80 && $F[5] > 80){print $_;}' < igc_125_animals.calls.ids.calls.passfilter.nodups.scores.bed > igc_125_animals.calls.ids.calls.passfilter.gt80scores.bed

perl -lane '$s = $F[1] - 36; $e = $F[1] + 36; system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$s\-$F[1] >> igc_locs_mapping_test.fa"); system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$F[1]\-$e >> igc_locs_mapping_test.fa");' < igc_125_animals.calls.ids.calls.passfilter.gt80scores.bed

bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta igc_locs_mapping_test.fa > igc_locs_mapping_test.sam
perl grep_marker_locs.pl < igc_locs_mapping_test.noheader.sam | perl -lane 'if($F[2] > 20 || $F[3] > 20){print $_;}' > igc_locs_mapping_test.mapq.assoc.tab

# I only recovered sequence on a few chromosomes, but it is heavily filtered. 
perl ~/sperl/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f igc_locs_mapping_test.mapq.assoc.tab -c 0 -m
```

#### Final filtration of markers by chromosome
|Entry                     | Count|
|:-------------------------|-----:|
|18                        |     6|
|5                         |   186|
|CH240_370M3_LILR_LRC      |    10|
|CH240_391K10_KIR          |     1|
|HF_LRC_hap1_KIR_LRC       |    17|
|LIB14413_LRC              |     9|
|LIB14427_MHC              |    38|
|TPI4222_A14_MHCclassI_MHC |    74|

All three IGCs are represented, though some more so than others. Let's prepare the list, then compare them to the locations of the BovineHD markers to try to select equidistant markers.

Thankfully, my old script to recursively select markers has survived! But I need to adapt it on a chromosome by chromosome basis to select the different markers.

Now to generate the input data to the [script](https://github.com/njdbickhart/perl_toolchain/blob/master/snp_utilities/gmsSNPSelectionStrategyGreedy.pl):

```bash
perl generate_mapq_dataframe.pl igc_locs_mapping_test.mapq.assoc.tab igc_125_animals.calls.ids.vcf.gz > igc_locs_mapping_test.total.assoc.tab

# Now I want to take the BovineHD remappings and see if we have covered any of our target haplotypes
perl -lane 'open(IN, "gunzip -c ../cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/bovinehd_illumina.snplocs.ars-ucd14.bed.gz |"); while(<IN>){chomp; @s = split(/\t/); if($s[0] eq $F[0] && $s[1] < $F[2] && $s[1] >= $F[1]){print join("\t", @s);}} close IN;' < nkc_regions_ars14.bed > nkc_regions.bovineHDmarkers.bed

# There were NO mappings to the assembled region of chr23
# The chr18 coords can be refined to: 18:63110866-63132207, but the region is too difficult to map to
# The NKC regions are pretty well dispersed on chr5
perl -e '$min = 908089099; $max = 0; while(<>){chomp; @F = split(/\t/); if($F[0] eq "5"){if($min > $F[1]){$min = $F[1];} if($max < $F[2]){$max = $F[2];}}} print "$min\t$max\n";' < nkc_regions.bovineHDmarkers.bed
98797892        99505322

# Of the other chromosomes with VCF variants discovered, only CH240_370M3_LILR_LRC had Bovine HD markers (63 total).
# I'm just going to select 6 equally spaced markers without regard to the BovineHD markerset -- so long as they don't overlap there's good reason to think that they're different

perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 18 -i igc_locs_mapping_test.total.assoc.tab -s 63110866 -e 63144659 -d 500 -m 6 > LRCA_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 5 -i igc_locs_mapping_test.total.assoc.tab -s 98792567 -e 99508055 -d 500 -m 6 > NKCA_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 23 -i igc_locs_mapping_test.total.assoc.tab -s 28460024 -e 28533695 -d 500 -m 6 > MHCA_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c CH240_391K10_KIR -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 100779 -d 500 -m 6 > LRC1_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c Domino_MHCclassI_gene2-5hap1_MHC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 339066 -d 500 -m 6 > MHC1_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c CH240_370M3_LILR_LRC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 179560 -d 500 -m 6 > LRC2_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c TPI4222_A14_MHCclassI_MHC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 407922 -d 500 -m 6 > MHC2_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c HF_LRC_hap1_KIR_LRC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 370760 -d 500 -m 6 > LRC3_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c LIB14427_MHC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 204406 -d 500 -m 6 > MHC3_snp_selections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c LIB14413_LRC -i igc_locs_mapping_test.total.assoc.tab -s 1 -e 133818 -d 500 -m 6 > LRC4_snp_selections.tab

# Let's see how many markers we selected in total:
wc -l *_snp_selections.tab
   1 LRC1_snp_selections.tab
   5 LRC2_snp_selections.tab
   6 LRC3_snp_selections.tab
   5 LRC4_snp_selections.tab
   4 LRCA_snp_selections.tab
   0 MHC1_snp_selections.tab
   7 MHC2_snp_selections.tab
   7 MHC3_snp_selections.tab
   0 MHCA_snp_selections.tab
   6 NKCA_snp_selections.tab
  41 total

```


So, I've selected 41 variants for the first round. Let's pass the list to John to see if he agrees with the locations and/or wants to select different types of variants.

<a name="last_sec"></a>
## Including the last section of chr5 NKC

Doro and John rightly pointed out that we are missing markers from the tail end of the ARS-UCDv14 NKC region. The reason for this is that I incorrectly cropped out the tail end of this region from the 99.5 Mb to the 99.8 Mb region.

Let's call SNPs in this area and then filter them using the same criteria as before. I am hoping to select about 5 more markers with reasonable MAF and mappability scores. **Target region is: 5:99,508,055-99,800,000**

> Assembler2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing

```bash
# Generating the bcf file
sbatch ~/sperl/sequence_data_scripts/samtoolsSelectiveMpileup.pl -b igc_variant_list_bams.list -s 5:99508055-99800000 -n NKCA2 -f /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta

# I also need to process the LRC and MHC regions differently as well
# LRC: 18:62400000-63450000
# MHC: 23:28250235-28651950

sbatch ~/sperl/sequence_data_scripts/samtoolsSelectiveMpileup.pl -b igc_variant_list_bams.list -s 18:62400000-63450000 -n LRCA2 -f /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta
sbatch ~/sperl/sequence_data_scripts/samtoolsSelectiveMpileup.pl -b igc_variant_list_bams.list -s 23:28250235-28651950 -n MHCA2 -f /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta

bcftools index LRCA2.mpileup.bcf
bcftools index MHCA2.mpileup.bcf

# Skipping indels to avoid ID collisions downstream
bcftools call -vmO z --skip-variants indels -o igc_125_animals.LRCA2.vcf.gz LRCA2.mpileup.bcf
bcftools call -vmO z --skip-variants indels -o igc_125_animals.MHCA2.vcf.gz MHCA2.mpileup.bcf
bcftools call -vmO z --skip-variants indels -o igc_125_animals.NKC2.vcf.gz NKCA2.mpileup.bcf

# generating IDs for each variant site
for i in LRCA2 MHCA2; do echo $i; bcftools annotate -o igc_125_animals.${i}.ids.vcf.gz -O z -I 'ARS\_PIRBRIGHT\_%CHROM\_%POS' igc_125_animals.${i}.vcf.gz; done
bcftools annotate -o igc_125_animals.NKC2.ids.vcf.gz -O z -I 'ARS\_PIRBRIGHT\_%CHROM\_%POS' igc_125_animals.NKC2.vcf.gz

module load plink/2.00alM-2017-05-22
# saving the first pruned IDs
mv plink2.prune.in firsttry.plink2.prune.in
mv plink2.prune.out firsttry.plink2.prune.out
for i in LRCA2 MHCA2 NKC2; do echo $i; plink2 --vcf igc_125_animals.${i}.ids.vcf.gz --indep-pairwise 50 5 0.5 --allow-extra-chr; mv plink2.prune.in ${i}.prune.in; mv plink2.prune.out ${i}.prune.out; done

wc -l *.prune.*
  20891 firsttry.plink2.prune.in
  38937 firsttry.plink2.prune.out
   5753 LRCA2.prune.in
  12139 LRCA2.prune.out
   1023 MHCA2.prune.in
   2060 MHCA2.prune.out
     56 plink2.prune.in.dups

# removing markers that violate HWE
module unload plink/2.00alM-2017-05-22
module load plink/1.90b4.4-2017-05-21
for i in LRCA2 MHCA2 NKC2; do echo $i; plink --vcf igc_125_animals.${i}.ids.vcf.gz --extract ${i}.prune.in --hardy --allow-extra-chr --threads 20 --double-id; mv plink.hwe plink.${i}.hwe; done
for i in LRCA2 MHCA2 NKC2; do echo $i; perl -lane 'if($F[0] eq "CHR"){next;} if($F[8] > 0.05){print $F[1];}' < plink.${i}.hwe > plink.${i}.hwe.remove; done

wc -l *.remove
 15116 plink.hwe.remove
  4533 plink.LRCA2.hwe.remove  <- 1221 remaining
   287 plink.MHCA2.hwe.remove	<- 737 remaining

module load samtools
for i in "18:62400000-63450000" "23:28250235-28651950" "5:99508055-99800000"; do echo $i | perl -lane 'print "$F[0]"; system("samtools mpileup -s -O -r $F[0] -b igc_variant_list_bams.list >> igc_variant_pileup_regions.new.tab");'; done
for i in "5:99508055-99800000"; do echo $i | perl -lane 'print "$F[0]"; system("samtools mpileup -s -O -r $F[0] -b igc_variant_list_bams.list >> igc_variant_pileup_regions.NKC.tab");'; done

perl -ne '@F = split(/\t/); $e = $F[1] + 1; $sum = 0; $c = 0; for($x = 7; $x < scalar(@F);$x += 5){if($F[$x] ne "*"){@bsegs = split(/,/, $F[$x]); $c += scalar(@bsegs); foreach $j (@bsegs){$sum += $j;}}} $avg = ($c > 0)? $sum / $c : 0; $sum = 0; $c = 0; print "$F[0]\t$F[1]\t$e\t$avg\n";' < igc_variant_pileup_regions.new.tab > igc_variant_pileup_regions.new.scores.bed

for i in LRCA2 MHCA2; do echo $i; perl -lane 'if($F[0] eq "CHR"){next;} if($F[8] < 0.05){print $F[1];}' < plink.${i}.hwe >> new.region.hwe.keep; done
for i in NKC2; do echo $i; perl -lane 'if($F[0] eq "CHR"){next;} if($F[8] < 0.05){print $F[1];}' < plink.${i}.hwe >> ${i}.region.hwe.keep; done

for i in LRCA2 MHCA2; do echo $i; perl -e 'chomp(@ARGV); open($IN, "< $ARGV[0]"); %snps; while(<$IN>){chomp; $snps{$_} = 1;} close $IN; open($IN, " gunzip -c $ARGV[1] |"); while(<$IN>){chomp; if($_ =~ /^#/){next;}else{@s = split(/\t/); if(exists($snps{$s[2]})){$e = $s[1] + 1; print "$s[0]\t$s[1]\t$e\t$s[2]\n";}}} close $IN;' new.region.hwe.keep igc_125_animals.${i}.ids.vcf.gz > igc_125_animals.new.${i}.passfilter.bed; done
for i in NKC2; do echo $i; perl -e 'chomp(@ARGV); open($IN, "< $ARGV[0]"); %snps; while(<$IN>){chomp; $snps{$_} = 1;} close $IN; open($IN, " gunzip -c $ARGV[1] |"); while(<$IN>){chomp; if($_ =~ /^#/){next;}else{@s = split(/\t/); if(exists($snps{$s[2]})){$e = $s[1] + 1; print "$s[0]\t$s[1]\t$e\t$s[2]\n";}}} close $IN;' NKC2.region.hwe.keep igc_125_animals.${i}.ids.vcf.gz > igc_125_animals.new.${i}.passfilter.bed; done

# associating marker regions with mapq scores
for i in LRCA2 MHCA2; do echo $i;  perl ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v igc_125_animals.new.${i}.passfilter.bed -l new.region.hwe.keep -m igc_variant_pileup_regions.new.scores.bed -o igc_125_animals.new.calls.${i}.scores.bed; done
for i in NKC2; do echo $i; perl ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v igc_125_animals.new.${i}.passfilter.bed -l NKC2.region.hwe.keep -m igc_variant_pileup_regions.NKC.tab -o igc_125_animals.new.calls.${i}.scores.bed; done

# This produced nothing, unfortunately. Will need to change the score threshold for NKC
perl -lane 'if($F[4] > 80 && $F[5] > 80){print $_;}' < igc_125_animals.new.calls.NKC2.scores.bed > igc_125_animals.new.calls.NKC2.scores.gt80.bed

# produced only 36! Might as well work with the whole dataset
perl -lane 'if($F[4] > 18 && $F[5] > 18){print $_;}' < igc_125_animals.new.calls.NKC2.scores.bed > igc_125_animals.new.calls.NKC2.scores.gt18.bed

cat igc_125_animals.new.calls.LRCA2.scores.gt80.bed igc_125_animals.new.calls.MHCA2.scores.gt80.bed | perl -lane '$s = $F[1] - 36; $e = $F[1] + 36; system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$s\-$F[1] >> igc_locs_new_mapping_test.fa"); system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$F[1]\-$e >> igc_locs_new_mapping_test.fa");'

cat igc_125_animals.new.calls.NKC2.scores.bed |  perl -lane '$s = $F[1] - 36; $e = $F[1] + 36; system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$s\-$F[1] >> igc_locs_nkc_mapping_test.fa"); system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$F[1]\-$e >> igc_locs_nkc_mapping_test.fa");'

module load bwa; bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta igc_locs_new_mapping_test.fa > igc_locs_new_mapping_test.sam
module load bwa; bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta igc_locs_nkc_mapping_test.fa > igc_locs_nkc_mapping_test.sam
perl grep_marker_locs.pl < igc_locs_new_mapping_test.sam | perl -lane 'if($F[2] > 20 || $F[3] > 20){print $_;}' > igc_locs_new_mapping.mapq.assoc.tab
perl grep_marker_locs.pl < igc_locs_nkc_mapping_test.sam | perl -lane 'if($F[2] > 20 || $F[3] > 20){print $_;}' > igc_locs_nkc_mapping.mapq.assoc.tab

for i in LRCA2 MHCA2; do echo $i;  perl generate_mapq_dataframe.pl igc_locs_new_mapping.mapq.assoc.tab igc_125_animals.${i}.ids.vcf.gz >> igc_locs_new_mapping.total.assoc.tab; done
for i in NKC2; do echo $i;  perl generate_mapq_dataframe.pl igc_locs_nkc_mapping.mapq.assoc.tab igc_125_animals.${i}.ids.vcf.gz > igc_locs_nkc_mapping.total.assoc.tab; done

perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 18 -i igc_locs_new_mapping.total.assoc.tab -s 62400000 -e 63450000 -d 500 -m 6 > LRCA2.new.snpselections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 23 -i igc_locs_new_mapping.total.assoc.tab -s 28250235 -e 28651950 -d 500 -m 6 > MHCA2.new.snpselections.tab
perl ~/sperl/snp_utilities/gmsSNPSelectionStrategyGreedy.pl -c 5 -i igc_locs_nkc_mapping.total.assoc.tab -s 99508055 -e 99800000 -d 500 -m 6 > NKCA2.new.snpselections.tab
```

## Generating flanking sequence files

I need to convert Doro's marker picks into a format for Geneseek to use.

#### Doro's marker selections

> Assembler2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing

```bash
# first, purge the original vcf of all chr18, 5 and 23 references so that I can include the extended regions
gunzip -c igc_125_animals.calls.ids.vcf.gz | perl -lane 'if($F[0] =~ /^#/){print $_;}elsif($F[0] eq "5" || $F[0] eq "18" || $F[0] eq "23"){next;}else{print $_;}' > igc_125_animals.calls.ids.noautos.vcf
gzip igc_125_animals.calls.ids.noautos.vcf

# Now, merge the assoc.tab files together
cat igc_locs_mapping_test.total.assoc.tab igc_locs_new_mapping.total.assoc.tab igc_locs_nkc_mapping.total.assoc.tab > igc_locs_mapping_combined.total.assoc.tab

# Generating Doro's tab delimited list
#NOTE: many of the haplotype regions are missing... also I'm not sure what the A11 haplotype is in the MHC region
vim doro_marker_selections.tab
dos2unix doro_marker_selections.tab

# now to cross my fingers and hope that this works!
module load samtools
perl ~/sperl/sequence_data_scripts/createFlankingSequenceSNPList.pl -f /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta -v igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.NKC2.ids.vcf.gz -p doro_marker_selections.tab -m igc_locs_mapping_combined.total.assoc.tab -o doro_marker_selections.round1.geneseek.tab

# lots of bugs, but it did work. Need to debug though

# Doro sent me a new list of marker variants to run through:
vim doro_marker_selections_round2.tab
perl ~/sperl/sequence_data_scripts/createFlankingSequenceSNPList.pl -f /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta -v igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.NKC2.ids.vcf.gz -p doro_marker_selections_round2.tab -m igc_locs_mapping_combined.total.assoc.tab -o doro_marker_selections.round2.geneseek.tab

# I'm not getting any upstream/downstream GMS for the variants because they're not represented in the total.assoc.tab file.
# Trying to generate the data to associate with the markers
perl -lane '$e = $F[3] + 1; print "$F[2]\t$F[3]\t$e\t$F[0]";' < doro_marker_selections.round2.geneseek.tab > doro_marker_selections.round2.geneseek.bed
perl -lane 'print "$F[0]";' < doro_marker_selections.round2.geneseek.tab > doro_marker_selections.round2.geneseek.list
perl ~/sperl/sequence_data_scripts/assessSamtoolsMpileupMapQSNPMarkers.pl -v doro_marker_selections.round2.geneseek.bed -v doro_marker_selections.round2.geneseek.bed -l doro_marker_selections.round2.geneseek.list -m igc_variant_pileup_regions.new.scores.bed -o doro_marker_selections.round2.geneseek.scores.bed

# Damn, that was part of my initial filtering step, so it was unneeded
cat doro_marker_selections.round2.geneseek.scores.bed | perl -lane '$s = $F[1] - 36; $e = $F[1] + 36; system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$s\-$F[1] >> doro_locs_new_mapping_test.fa"); system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$F[1]\-$e >> doro_locs_new_mapping_test.fa");'
module load bwa; bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta doro_locs_new_mapping_test.fa > doro_locs_new_mapping_test.sam
perl grep_marker_locs.pl < doro_locs_new_mapping_test.sam > doro_locs_new_mapping_test.mapq.assoc.tab
perl generate_mapq_dataframe.pl doro_locs_new_mapping_test.mapq.assoc.tab igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.NKC2.ids.vcf.gz > doro_locs_new_mapping_test.total.assoc.tab

perl ~/sperl/sequence_data_scripts/createFlankingSequenceSNPList.pl -f /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta -v igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.NKC2.ids.vcf.gz -p doro_marker_selections_round2.tab -m doro_locs_new_mapping_test.total.assoc.tab -o doro_marker_selections.round2.geneseek.tab
```

Damn, I missed some markers from Doro's selections! Let me add them to the list and then generate the updated files.

```bash
# I found out why I missed those markers -- the VCF I used for the NKC region is from my old coordinates!
gunzip -c igc_125_animals.calls.ids.vcf.gz | perl -lane 'if($_ =~ /^#/){print $_;}elsif($F[0] eq "5"){print $_;}' > igc_125_animals.calls.ids.vcf.NKCA2.ids.vcf
# And the tail end calls from the other vcf...
gunzip -c igc_125_animals.NKC2.ids.vcf.gz | perl -lane 'if($F[0] eq "5" && $F[1] > 99508015){print $_;}' >> igc_125_animals.calls.ids.vcf.NKCA2.ids.vcf

vim doro_locs_new_mapping_test.fa
cat doro_updated_list.round3.tab | perl -lane '$s = $F[1] - 36; $e = $F[1] + 36; system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$s\-$F[1] >> doro_locs_new_mapping_test.fa"); system("samtools faidx /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0]:$F[1]\-$e >> doro_locs_new_mapping_test.fa");'
module load bwa; bwa mem /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta doro_locs_new_mapping_test.fa > doro_locs_new_mapping_test.sam
perl grep_marker_locs.pl < doro_locs_new_mapping_test.sam > doro_locs_new_mapping_test.mapq.assoc.tab
perl generate_mapq_dataframe.pl doro_locs_new_mapping_test.mapq.assoc.tab igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.calls.ids.vcf.NKCA2.ids.vcf.gz > doro_locs_new_mapping_test.total.assoc.tab

cat doro_updated_list.round3.tab doro_marker_selections_round2.tab > doro_marker_selections_round3.tab
perl ~/sperl/sequence_data_scripts/createFlankingSequenceSNPList.pl -f /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta -v igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.calls.ids.vcf.NKCA2.ids.vcf.gz -p doro_marker_selections_round2.tab -m doro_locs_new_mapping_test.total.assoc.tab -o doro_marker_selections.round3.geneseek.tab
```

#### My marker selections

Now I will format my marker selections using the same criteria to enter into the spreadsheet for geneseek.

```bash
# First, I need to generate a list of variant site and locations like Doro's list
for i in MHCA2.new.snpselections.tab NKCA2.new.snpselections.tab MHC3_snp_selections.tab MHC2_snp_selections.tab LRCA2.new.snpselections.tab LRC4_snp_selections.tab LRC3_snp_selections.tab LRC2_snp_selections.tab LRC1_snp_selections.tab; do perl -lane 'print $F[2];' < $i; done > derek_snp_selection.ids.list

for i in igc_125_animals.calls.ids.noautos.vcf.gz igc_125_animals.LRCA2.ids.vcf.gz igc_125_animals.MHCA2.ids.vcf.gz igc_125_animals.NKC2.ids.vcf.gz; do perl -e 'chomp @ARGV; %h; open(IN, "<$ARGV[0]"); while(<IN>){chomp; $h{$_} = 1;} close IN; open(IN, "gunzip -c $ARGV[1] |"); while(<IN>){chomp; if($_ =~ /^#/){next;} @s = split(/\t/); if(exists($h{$s[2]})){print "$s[0]\t$s[1]\t$s[3]\t$s[4]\t$s[2]\n";}} close IN;' derek_snp_selection.ids.list $i; done > derek_snp_selection.round1.tab

# Unfortunately, my previous decision NOT to remove INDELS off the bat rears its ugly head here
# I'll remove them manually (removes 4 markers)
vim derek_snp_selection.round1.tab

perl ~/sperl/sequence_data_scripts/createFlankingSequenceSNPList.pl -f /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/ARS-UCD1.0.14.clean.wIGCHaps.fasta -v igc_125_animals.calls.ids.noautos.vcf.gz,igc_125_animals.LRCA2.ids.vcf.gz,igc_125_animals.MHCA2.ids.vcf.gz,igc_125_animals.NKC2.ids.vcf.gz -p derek_snp_selection.round1.tab -m igc_locs_mapping_combined.total.assoc.tab -o derek_marker_selections.round1.geneseek.tab

wc -l derek_marker_selections.round1.geneseek.tab doro_marker_selections.round2.geneseek.tab
   41 derek_marker_selections.round1.geneseek.tab
  158 doro_marker_selections.round2.geneseek.tab
  199 total

# Some of Doro's selections have a high MAF and low mappability, but they're otherwise a good set
# We'll see how the group reacts before we filter
```

<a name="neogen"></a>
## Tabulating and interpreting Neogen results

We Just got the data back from Neogen. I want to run some preliminary stats on the calls and send them to the group.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results

```bash
# I wrote a script to condense the data into a format that will be useful for plotting
perl calculate_summary_stats.pl

perl -e '%c = ("18" => "LRC", "23" => "MHC", "5" => "NKC", "KIR" => "LRC"); $h = <>; chomp $h; @s = split(/\t/, $h); foreach $x (@s){ @bsegs = split(/_/, $x); $cat = ""; if(exists($c{$bsegs[0]})){$cat = $c{$bsegs[0]};}else{$cat = $bsegs[0];} print "$x\t$cat\n";}' < var_names.tab > var_names.cats.tab
```

> pwd: /home/dbickhart/share/grants/immune_gene_cluster_grant/variant_selection/neogen_first_round_results

```R
library(dplyr)
library(ggplot2)

samples <- read.delim("neogen_results_sampsummary.tab", header=TRUE)
samples.alt <- mutate(samples, CallRate = Genotyped / TotMarkers, HetRate = Hets / (Hom + Hets))
summary(samples.alt)
        Sample       TotMarkers      Genotyped        Missing      
 G0052_0001:   1   Min.   :49.00   Min.   :17.00   Min.   : 5.000  
 G0052_0002:   1   1st Qu.:67.00   1st Qu.:59.00   1st Qu.: 6.000  
 G0052_0031:   1   Median :67.00   Median :60.00   Median : 7.000  
 G0052_0033:   1   Mean   :66.95   Mean   :59.87   Mean   : 7.081  
 G0052_0034:   1   3rd Qu.:67.00   3rd Qu.:61.00   3rd Qu.: 8.000  
 G0052_0037:   1   Max.   :67.00   Max.   :62.00   Max.   :47.000  
 (Other)   :1791                                                   
      Hets           HomAlt          HomRef          NumAlt     
 Min.   : 0.00   Min.   : 5.00   Min.   : 7.00   Min.   :12.00  
 1st Qu.:10.00   1st Qu.:18.00   1st Qu.:26.00   1st Qu.:49.00  
 Median :12.00   Median :20.00   Median :28.00   Median :52.00  
 Mean   :11.66   Mean   :20.04   Mean   :28.18   Mean   :51.73  
 3rd Qu.:14.00   3rd Qu.:22.00   3rd Qu.:30.00   3rd Qu.:56.00  
 Max.   :22.00   Max.   :29.00   Max.   :39.00   Max.   :65.00  
                                                                
     NumRef         CallRate         HetRate      
 Min.   :18.00   Min.   :0.2537   Min.   :0.0000  
 1st Qu.:65.00   1st Qu.:0.8806   1st Qu.:0.1639  
 Median :68.00   Median :0.8955   Median :0.1967  
 Mean   :68.01   Mean   :0.8936   Mean   :0.1948  
 3rd Qu.:72.00   3rd Qu.:0.9104   3rd Qu.:0.2295  
 Max.   :88.00   Max.   :0.9254   Max.   :0.3684

cor(samples.alt[,c(2,3,4,5,6,7,8,9,10,11)])
            TotMarkers   Genotyped     Missing        Hets     HomAlt
TotMarkers  1.00000000  0.51358993 -0.40495366  0.12077102  0.2277562
Genotyped   0.51358993  1.00000000 -0.99251415  0.21364192  0.3524533
Missing    -0.40495366 -0.99251415  1.00000000 -0.21047028 -0.3431621
Hets        0.12077102  0.21364192 -0.21047028  1.00000000 -0.4422142
HomAlt      0.22775618  0.35245330 -0.34316213 -0.44221423  1.0000000
HomRef      0.23787454  0.56843306 -0.57187355 -0.25403239 -0.1572339
NumAlt      0.31825561  0.50679649 -0.49475145  0.04181606  0.8776333
NumRef      0.29626396  0.67361855 -0.67565002  0.20088900 -0.3618677
CallRate    0.51358993  1.00000000 -0.99251415  0.21364192  0.3524533
HetRate     0.04241602 -0.04856465  0.05778864  0.95852953 -0.5414084
               HomRef      NumAlt      NumRef    CallRate     HetRate
TotMarkers  0.2378745  0.31825561  0.29626396  0.51358993  0.04241602
Genotyped   0.5684331  0.50679649  0.67361855  1.00000000 -0.04856465
Missing    -0.5718735 -0.49475145 -0.67565002 -0.99251415  0.05778864
Hets       -0.2540324  0.04181606  0.20088900  0.21364192  0.95852953
HomAlt     -0.1572339  0.87763330 -0.36186774  0.35245330 -0.54140843
HomRef      1.0000000 -0.31091475  0.89644616  0.56843306 -0.40984963
NumAlt     -0.3109147  1.00000000 -0.29574732  0.50679649 -0.09084570
NumRef      0.8964462 -0.29574732  1.00000000  0.67361855  0.02407009
CallRate    0.5684331  0.50679649  0.67361855  1.00000000 -0.04856465
HetRate    -0.4098496 -0.09084570  0.02407009 -0.04856465  1.00000000


library(ggcorrplot)
ggcorrplot(cor(samples.alt[,c(2,3,4,5,6,7,8,9,10,11)]), hc.order=TRUE, type="upper", outline.col = "white")

# Now Let's load and intersect the sample information 
sampinfo <- read.delim("holstein_state_info.tab", header=TRUE)
samples.alt.full <- left_join(samples.alt, sampinfo)

# Now to generate a sample callrate plot
ggplot(samples.alt.full, aes(CallRate)) + stat_bin(binwidth=0.02, aes(fill=TB_pheno)) + stat_bin(binwidth=0.02, geom="text", aes(label=..count.., alpha=..count..), vjust=-1.5) + facet_grid(. ~ TB_pheno) + ggtitle(label="Custom Genotype Callrate Per TB_Phenotype")

# Now let's plot a scatterpot of the number of Alt alleles and number of Reference alleles per TB phenotype
allele <- ggplot(samples.alt.full, aes(x=NumRef, y=NumAlt, color=TB_pheno, alpha=CallRate)) + geom_count() + ggtitle(label="Reference and Alt Allele counts among samples")

# Final plot: a density of heterozygous calls 
ggplot(samples.alt.full, aes(x=Hets, fill=TB_pheno)) + geom_density(alpha=0.4) + xlab(label="Heterozygous Genotypes per sample") + ggtitle(label="Heterozygous alleles per TB Phenotype")
dev.copy2pdf(file="density_het_alleles_persample.pdf", useDingbats=FALSE)

## Now for the per marker site stats
markers <- read.delim("neogen_results_variantsummary.tab", header=TRUE)
markers.alt <- mutate(markers, CallRate = 1 - (Missing / (Genotyped + Missing)), p = NumRef / (NumRef + NumAlt), q = NumAlt / (NumRef + NumAlt))
markers.alt.full <- mutate(markers.alt, ExpAA = p ** 2 * (HomAlt + HomRef + Hets), ExpAa = 2 * p * q * (HomAlt + HomRef + Hets), Expaa = q ** 2 * (HomAlt + HomRef + Hets), chi = (((HomRef - ExpAA) ** 2) / ExpAA) + (((Hets - ExpAa) ** 2) / ExpAa) + (((HomAlt - Expaa)) ** 2 / Expaa))
markers.short <- select(markers.alt.full, -p, -q, -ExpAA, -ExpAa, -Expaa)

summary(markers.short)
         Marker     Genotyped       Missing            Hets       
 18_62404864: 1   Min.   :   0   Min.   :   0.0   Min.   :   0.0  
 18_62438492: 1   1st Qu.:1695   1st Qu.:  10.0   1st Qu.:   0.0  
 18_62460291: 1   Median :1777   Median :  19.0   Median : 158.0  
 18_62469715: 1   Mean   :1606   Mean   : 189.9   Mean   : 312.7  
 18_62559417: 1   3rd Qu.:1786   3rd Qu.: 101.5   3rd Qu.: 557.5  
 18_62575095: 1   Max.   :1791   Max.   :1797.0   Max.   :1790.0  
 (Other)    :61                                                   
     HomAlt           HomRef           NumAlt           NumRef    
 Min.   :   0.0   Min.   :   0.0   Min.   :   0.0   Min.   :   0  
 1st Qu.:   2.5   1st Qu.:  10.5   1st Qu.: 146.5   1st Qu.: 470  
 Median : 119.0   Median : 655.0   Median : 939.0   Median :2001  
 Mean   : 537.4   Mean   : 755.7   Mean   :1387.4   Mean   :1824  
 3rd Qu.:1001.0   3rd Qu.:1419.0   3rd Qu.:2542.0   3rd Qu.:3156  
 Max.   :1791.0   Max.   :1791.0   Max.   :3582.0   Max.   :3582  
                                                                  
    CallRate           chi           
 Min.   :0.0000   Min.   :   0.0001  
 1st Qu.:0.9435   1st Qu.:   0.4622  
 Median :0.9894   Median :   7.3879  
 Mean   :0.8943   Mean   : 299.4821  
 3rd Qu.:0.9944   3rd Qu.: 481.0365  
 Max.   :1.0000   Max.   :1790.0000  
                  NA's   :21

# First a simple correlation plot:
ggcorrplot(cor(markers.short[,c(2,3,4,5,6,7,8,9,10)], use ="complete.obs"), hc.order=TRUE, type="upper", outline.col = "white")

# Now to associate marker regions with the markers
markercats <- read.delim("var_names.cats.tab", header=FALSE)
colnames(markercats) <- c("Marker", "IGC")
markers.short.full <- left_join(markers.short, markercats)

# Now a callrate plot:
ggplot(markers.short.full, aes(CallRate)) + stat_bin(binwidth=0.1, aes(fill=IGC)) + stat_bin(binwidth=0.1, geom="text", aes(label=..count.., alpha=..count..), vjust=-1.5) + facet_grid(. ~ IGC) + ggtitle(label = "Custom Genotype Callrate per Marker Site")
dev.copy2pdf(file="genotype_callrate_permaker.pdf", useDingbats=FALSE)

# Now one more histogram 
markers.short.full <- mutate(markers.short.full, HetRate = Hets / (Hets + HomAlt + HomRef))
ggplot(markers.short.full, aes(x = HetRate, y = chi, alpha = CallRate, color = IGC)) + geom_point() + facet_grid( . ~ IGC) + geom_hline(yintercept=3.84) + xlab(label="Heterozygous rate (Hets / SampleNum)") + ylab(label="Chi squared HWE value (solid line is 95% confidence)") + ggtitle(label="Call rate distribution for markers and Hardy Weinberg Equillibrium test")
dev.copy2pdf(file="marker_status_hwe_permaker.pdf", useDingbats=FALSE)

# Final plot = let's separate markers into categories of useability
markers.use <- mutate(markers.short.full, Use = ifelse(HetRate > 0.025 & HetRate < 1.0 & CallRate > 0.85, TRUE, FALSE))
ggplot(markers.use, aes(Use, fill = IGC)) + stat_count() + stat_count(geom="text", aes(label=..count..), vjust=-1.5) + facet_grid(. ~ IGC) + xlab(label="Suitable marker site?") + ylab(label="Count of Useable Markers") + ggtitle(label="Marker sites likely to be suitable for analysis")
dev.copy2pdf(file="count_of_passing_markers.pdf", useDingbats=FALSE)

```

And finally, here is the analysis perl script used to partition the results:

#### calculate_summary_stats.pl

```perl
#!/usr/bin/perl
# This is a script designed to process the neogen results and calculate summary stats for plotting in R

use strict;
my $input = "neogen_results_5_18_2018.tab";
my $persamp = "neogen_results_sampsummary.tab";
my $pervar = "neogen_results_variantsummary.tab";

my %varRefAlleles = (
        "18_62404864" => "G",
        "18_62438492" => "T",
        "18_62460291" => "C",
        "18_62469715" => "C",
        "18_62559417" => "A",
        "18_62575095" => "A",
        "18_62644240" => "T",
        "18_62670367" => "G",
        "18_62723255" => "C",
        "18_62774779" => "C",
        "18_62812297" => "A",
        "18_62824331" => "G",
        "18_62914384" => "C",
        "18_62985809" => "C",
        "18_63089154" => "T",
        "18_63138493" => "C",
        "18_63222615" => "C",
        "18_63269795" => "C",
        "18_63308589" => "C",
        "18_63356424" => "C",
        "23_28481840" => "T",
        "23_28518797" => "A",
        "23_28535354" => "C",
        "23_28651088" => "A",
        "5_98985645" => "A",
        "5_99036250" => "C",
        "5_99073686" => "T",
        "5_99228484" => "C",
        "5_99266853" => "T",
        "5_99333595" => "A",
        "5_99392730" => "C",
        "5_99445508" => "A",
        "5_99560651" => "C",
        "5_99706728" => "A",
        "5_99763326" => "G",
        "5_99766613" => "C",
        "5_99780983" => "C",
        "KIR_41480" => "T",
        "LRC_106729" => "G",
        "LRC_118549" => "T",
        "LRC_174904" => "G",
        "LRC_61352" => "A",
        "LRC_65014" => "G",
        "LRC_72198" => "G",
        "LRC_81028" => "T",
        "LRC_82562" => "G",
        "LRC_98357" => "G",
        "LRC_98785" => "A",
        "MHC_103858" => "C",
        "MHC_115082" => "G",
        "MHC_121069" => "G",
        "MHC_143922" => "A",
        "MHC_154399" => "A",
        "MHC_208321" => "T",
        "MHC_260913" => "A",
        "MHC_286137" => "C",
        "MHC_317666" => "A",
        "MHC_324231" => "T",
        "MHC_356873" => "T",
        "MHC_359718" => "A",
        "MHC_380177" => "T",
        "MHC_400205" => "A",
        "MHC_43656" => "G",
        "MHC_59013" => "A",
        "MHC_63411" => "A",
        "MHC_86084" => "A",
        "MHC_9213" => "A");

my %sampStats; #{sample} = [variantsGtyped, missing, hetcalls, homcalls, homref, altcalls, refcalls]
my %varStats; #{variant} = [callcount, missingcount, hetcalls, homcalls, homref, altcalls, refcalls]

open(my $IN, "< $input") || die "Could not open input file: $input!\n";
my $head = <$IN>;
chomp $head;
my @hsegs = split(/\t/, $head);
my %varIdx;
for(my $x = 2; $x < scalar(@hsegs); $x++){
	$varIdx{$x} = $hsegs[$x];
	$varStats{$hsegs[$x]} = [0, 0, 0, 0, 0, 0, 0];
}
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my $sampGtype = 0; my $sampMiss = 0; my $sampHet = 0; my $sampHom = 0; my $sampHRef = 0; my $sampRef = 0; my $sampAlt = 0;
	for(my $x = 2; $x < scalar(@segs); $x++){
		my $var = $varIdx{$x};
		my $hallele = $varRefAlleles{$var};
		if(!defined($segs[$x]) || $segs[$x] eq ""){
			$sampMiss++;
			$varStats{$var}->[1] += 1;
		}else{
			$sampGtype++;
			$varStats{$var}->[0] += 1;
			if(length($segs[$x]) > 1){
				$sampHet++;
				$sampRef++; $sampAlt++;
				$varStats{$var}->[2] += 1;
				$varStats{$var}->[5] += 1;
				$varStats{$var}->[6] += 1;
			}elsif($segs[$x] eq $hallele){
				$sampHRef++;
				$sampRef += 2;
				$varStats{$var}->[4] += 1;
				$varStats{$var}->[6] += 2;
			}else{
				$sampHom++;
				$sampAlt += 2;
				$varStats{$var}->[3] += 1;
				$varStats{$var}->[5] += 2;
			}
		}
	}
	$sampStats{$segs[0]} = [$sampGtype, $sampMiss, $sampHet, $sampHom, $sampHRef, $sampAlt, $sampRef];
}
close $IN;

open(my $OUT, "> $persamp");
print {$OUT} "Sample\tTotMarkers\tGenotyped\tMissing\tHets\tHomAlt\tHomRef\tNumAlt\tNumRef\n";
foreach my $s (sort {$a cmp $b} keys(%sampStats)){
	my $entry = $sampStats{$s};
	my $tot = $entry->[0] + $entry->[1];
	print {$OUT} "$s\t$tot\t" . join("\t", @{$entry}) . "\n";
}
close $OUT;

open(my $OUT, "> $pervar");
#{variant} = [callcount, missingcount, hetcalls, homcalls]
print {$OUT} "Marker\tGenotyped\tMissing\tHets\tHomAlt\tHomRef\tNumAlt\tNumRef\n";
foreach my $v (sort{$a cmp $b} keys(%varStats)){
	print {$OUT} "$v\t" . join("\t", @{$varStats{$v}}) . "\n";
}
close $OUT;

exit;

```

## Flirting with graph format

These are my barebones notes on using fermi to assemble our contigs into graph format for plotting.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
perl ~/sperl/sequence_data_pipeline/splitFASTAintoOVLPChunks.pl -f total_assembled_bacs.fa -o total_assembled_bacs.ovlp.segs.fq

mv temp.fq total_assembled_bacs.ovlp.segs.fq

/mnt/nfs/nfs2/bickhart-users/binaries/fermikit/fermi.kit/fermi2.pl unitig -s 1m -t 16 -l 100 -p total_assembled_segs total_assembled_bacs.ovlp.segs.fq > total_assembled_segs.mak

make -f total_assembled_segs.mak

/mnt/nfs/nfs2/bickhart-users/binaries/mag2gfa/mag2gfa total_assembled_segs.mag.gz > total_assembled_segs.mag.gfa

# I ran bandage on my linux virtualbox. Looks OK, but is a mess! Let's try a targetted area first

samtools faidx total_assembled_bacs.fa
samtools faidx total_assembled_bacs.fa Domino_MHCclassI_gene2-5hap1_MHC TPI4222_A14_MHCclassI_MHC LIB14427_MHC > graph_mhc_regions.fa

samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta 23:28460024-28533695 >> graph_mhc_regions.fa

perl ~/sperl/sequence_data_pipeline/splitFASTAintoOVLPChunks.pl -f graph_mhc_regions.fa -o graph_mhc_regions.segs.fq
mv temp.fq graph_mhc_regions.segs.fq

/mnt/nfs/nfs2/bickhart-users/binaries/fermikit/fermi.kit/fermi2.pl unitig -s 1m -t 16 -l 100 -p graph_mhc_regions.segs.fermi graph_mhc_regions.segs.fq > graph_mhc_regions.segs.fermi.mak
make -f graph_mhc_regions.segs.fermi.mak
/mnt/nfs/nfs2/bickhart-users/binaries/mag2gfa/mag2gfa graph_mhc_regions.segs.fermi.mag.gz > graph_mhc_regions.segs.fermi.mag.gfa

```

<a name="placing"></a>
## Placing variant sites

I am going to give Kiranmayee mapping locations of the haplotypes and the SNP markers. The first thing that I realized is that I sent her coordinates from the v1.25 assembly, but the mapping coordinates I have are from the v14 assembly. Let's remap the BovineHD coordinates to that assembly first.

> Assembler2: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/

```bash
# Generating the tab file for comparisons
perl ~/sperl/assembly_scripts/alignAndOrderSnpProbes.pl -a ../../cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta -p /mnt/nfs/nfs2/dbickhart/dominette_asm/recombination/BovineHD_B1.probseq.umdcoords.fa -o ARS-UCD1.0.14.BovineHD

# Generating the bed file for marker site coordinates
perl -lane '$e = $F[3] + 1; print "$F[2]\t$F[3]\t$e\t$F[0]";' < ../doro_marker_selections.round3.geneseek.tab > marker_selections_arsv1.14.locs.bed

# Now to do the mappings. I know that the last haplotypes on the v14_igc assembly are the haplotypes that need to be aligned. Let's do the comparison of this assembly against the previous version that did not contain the haplotypes so as to avoid perfect identity overlaps
# Generating the list of haplotypes for alignment
samtools faidx ../../cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta CH240_391K10_KIR Domino_MHCclassI_gene2-5hap1_MHC CH240_370M3_LILR_LRC TPI4222_A14_MHCclassI_MHC HF_LRC_hap1_KIR_LRC LIB14427_MHC LIB14413_LRC > alt_haplotype_segments.fa

# Now, performing the alignment
/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -x asm5 -t 10 ../../cattle_asms/ars_ucd_114/ARS-UCD1.0.14.clean.fasta alt_haplotype_segments.fa > alt_haplotype_segments.paf

# Adding the cigar string to the paf file:
/mnt/nfs/nfs2/bickhart-users/binaries/minimap2/minimap2 -x asm5 -c -L -t 10 ../../cattle_asms/ars_ucd_114/ARS-UCD1.0.14.clean.fasta alt_haplotype_segments.fa > alt_haplotype_segments.cigar.paf
```