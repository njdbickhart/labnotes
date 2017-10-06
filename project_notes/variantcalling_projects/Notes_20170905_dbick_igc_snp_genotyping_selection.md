# IGC variant selection and filtering
---
**9/5/2017**

These are my notes on the selection and filtering of sequence variants suitable for genotyping large populations of cattle. The overall goal is to select 50 SNP targets for downstream genotyping and then branch out for there.

## Table of Contents
* [Selecting regions for SNP analysis](#selecting)
* [Generating VCFs for genomic region subsets](#generating)
* [Including the last section of chr5 NKC](#last_sec)

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

