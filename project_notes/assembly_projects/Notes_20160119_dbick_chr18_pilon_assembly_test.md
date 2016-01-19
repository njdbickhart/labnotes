# Creating accurate scaffolds from BAC-PACBIO sequencing
---
*1/19/2016*

These are my notes on taking existing BAC sequences from the chr18 bac clones, polishing them and assembling them into larger contiguous segments (or at least scaffolding them!).

The quivered BAC clone contigs are unrefined, and likely need pilon correction. I'm going to create a dummy fasta file (reference genome + extra segments to be polished) and Pilon-correct it. Then I'm going to use Velvet to try to do a quick graph-based assembly. If that fails, I'll go pure consensus-overlap (Celera?).

#### Here are the pipeline steps:
1. Identify regions on the original reference genome that are covered by the BAC sequences and remove them
2. Place BAC fasta entries at the end of the reference genome as extra "chrs"
3. Align reads from an individual to the segments
4. Run Pilon on the BAM file (or a subset of the BAM with only the clones we're interested in?)
5. Extract the corrected entries
6. Perform OLC or graph-based assembly

Let's start by getting things ready and proceeding in a step-wise fashion with the existing BAC sequences from John's data.

Files:
* LIB14363_unitig_3_oriented_vector_trimmed.fasta  
* LIB14414_unitig_273.fasta  
* LIB14435_unitig_94_vector_trim.fasta

> Blade14: /mnt/iscsi/vnx_gliu_7/john_assembled_contigs

```bash
# I'm going to test to see if I can identify the regions of the chromosome that need to be removed from a simple alignment
cat *.fasta > pilon_testrun/combined_chr18_fastas.fa
bwa mem ../reference/umd3_kary_unmask_ngap.fa pilon_testrun/combined_chr18_fastas.fa > pilon_testrun/combined_chr18_fastas.sam

perl -lane 'if($F[0] =~ /^@/ || $F[1] > 16){next;}else{print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t" . (length($F[9]) + $F[3]);}' < pilon_testrun/combined_chr18_fastas.sam
	16      chr18   57435857        60      57626394
	0       chr18   57565716        60      57716083
	16      chr18   57676772        60      57835648

# This is obviously not the full picture because the alignments aren't complete, but it's a good start!
# Creating the mask bed file so that we can exclude this region from alignment
echo -e "chr18\t57435857\t57835648" > pilon_testrun/chr18_mask_location.bed
~/bedtools-2.17.0/bin/maskFastaFromBed -fi ../reference/umd3_kary_unmask_ngap.fa -bed pilon_testrun/chr18_mask_location.bed -fo pilon_testrun/umd3_chr18_region_mask.fa

# Now to make and finalize the reference file
cat pilon_testrun/umd3_chr18_region_mask.fa pilon_testrun/combined_chr18_fastas.fa > pilon_testrun/umd3_chr18_masked_combined_unitigs.fa
```

> Blade14: /mnt/iscsi/vnx_gliu_7/john_assembled_contigs/pilon_testrun

```bash
bwa index umd3_chr18_masked_combined_unitigs.fa
samtools faidx umd3_chr18_masked_combined_unitigs.fa
perl ~/perl_toolchain/sequence_data_pipeline/runMergedBamPipeline.pl --fastqs ../Arlinda-Chief.10x.spreadsheet.tab --reference umd3_chr18_masked_combined_unitigs.fa --config test_pipeline.cnfg --threads 15 --output arlinda_chr18
```