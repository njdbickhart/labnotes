# IGC variant selection and filtering
---
**9/5/2017**

These are my notes on the selection and filtering of sequence variants suitable for genotyping large populations of cattle. The overall goal is to select 50 SNP targets for downstream genotyping and then branch out for there.

## Table of Contents


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

**TODO tomorrow: make a box plot of variants per region per animal**

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
|CHR                              |     1|
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
perl -lane '$e = $F[1] + 1; $sum = 0; $c = 0; for($x = 7; $x < scalar(@F);$x += 5){if($F[$x] ne "*"){@bsegs = split(/,/, $F[$x]); $c += scalar(@bsegs); foreach $j (@bsegs){$sum += $j;}}} $avg = ($c > 0)? $sum / $c : 0; $sum = 0; $c = 0; print "$F[0]\t$F[1]\t$e\t$avg";' < igc_variant_pileup_regions.tab > igc_variant_pileup_regions.scores.bed